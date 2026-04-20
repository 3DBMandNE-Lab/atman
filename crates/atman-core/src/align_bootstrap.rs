//! Subject-level bootstrap for cross-cohort archetype alignment.
//!
//! For each bootstrap iteration, resample subjects within every
//! cohort with replacement, re-run FastICA per cohort on the
//! resampled matrix, run the existing cosine-similarity archetype
//! alignment, and record — for each point-estimate archetype —
//! which cohorts the matched bootstrap archetype was recovered in.
//!
//! Output per point-estimate archetype: mean cohort count, fraction
//! of iterations it was recovered in every cohort (`prob_universal`),
//! fraction in at least two cohorts (`prob_multi`), and a percentile
//! CI over the cohort-count distribution.
//!
//! v1 scope (tracked in feature_requests.md Priority 3):
//! - Cosine similarity only (other metrics deferred).
//! - Matching is via representative-program best-cosine; each
//!   point-estimate archetype is represented by the loading vector
//!   from its first-cohort program, and every bootstrap archetype
//!   is represented by its first-cohort program. Match is the
//!   bootstrap archetype whose representative has maximum cosine
//!   similarity against the point-estimate representative, with a
//!   configurable floor (`match_tau`).
//! - No BCa CI; percentile only.
//! - No alignment-entropy metric.
//!
//! Determinism: all randomness comes from
//! [`crate::ica::Xoshiro256pp`] seeded per iteration via a
//! SplitMix64 derivation over `(seed, iter)`.

use crate::align::{build_archetypes, AlignMetric, AlignedProgram};
use crate::ica::{canonicalize_ica, fast_ica, CanonicalIca, Xoshiro256pp};
use crate::stats::cosine;

#[derive(Debug, Clone)]
pub struct CohortMatrix {
    pub label: String,
    /// Sample-major (rows = samples, columns = proteins).
    pub data: Vec<Vec<f64>>,
    /// Shared protein label ordering used across all cohorts for
    /// alignment. Every cohort's `data` columns must be indexed
    /// identically into this slice.
    pub protein_labels: Vec<String>,
}

#[derive(Debug, Clone, Copy)]
pub struct BootstrapParams {
    pub k: usize,
    pub n_boot: usize,
    pub seed: u64,
    pub top_n: usize,
    pub cosine_tau: f64,
    /// Floor on cosine similarity between the point-estimate
    /// representative and a bootstrap program for the bootstrap
    /// archetype containing that program to be considered a match.
    pub match_tau: f64,
    pub max_iter: usize,
    pub tol: f64,
    /// Bootstrap iterations abort with an error when any cohort has
    /// fewer than this many distinct subjects.
    pub min_subjects: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BootstrapRow {
    pub archetype_id: usize,
    pub observed_n_cohorts: usize,
    pub observed_cohorts: Vec<String>,
    pub bootstrap_mean_n_cohorts: f64,
    pub bootstrap_prob_universal: f64,
    pub bootstrap_prob_multi: f64,
    /// Percentile CI over the bootstrap `n_cohorts` distribution
    /// at 2.5% / 97.5% (the two-sided 95% band).
    pub ci_lower_n_cohorts: usize,
    pub ci_upper_n_cohorts: usize,
    /// Fraction of bootstrap iterations where **any** archetype was
    /// matched to this point-estimate archetype at all. When this is
    /// low the other statistics are noisy.
    pub bootstrap_match_rate: f64,
}

/// SplitMix64-derived sub-seed for iteration `iter` under top-level
/// `seed`. Same primitive as `ica_null::derive_sub_seed`.
fn derive_sub_seed(seed: u64, iter: usize) -> u64 {
    let mut z = seed.wrapping_add(0x9E3779B97F4A7C15_u64.wrapping_mul(iter as u64 + 1));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// Sample `n_samples` indices with replacement from `0..n` using
/// rejection-sampling on the RNG.
fn bootstrap_indices(rng: &mut Xoshiro256pp, n: usize, n_samples: usize) -> Vec<usize> {
    let mut out = Vec::with_capacity(n_samples);
    for _ in 0..n_samples {
        let bound = n as u64;
        loop {
            let v = rng.next_normal().to_bits();
            let limit = u64::MAX - u64::MAX % bound;
            if v < limit {
                out.push((v % bound) as usize);
                break;
            }
        }
    }
    out
}

/// Resample rows of `data` with replacement, returning a new
/// sample × protein matrix of the same shape.
pub fn resample_rows(data: &[Vec<f64>], rng: &mut Xoshiro256pp) -> Vec<Vec<f64>> {
    let n = data.len();
    if n == 0 {
        return Vec::new();
    }
    let idx = bootstrap_indices(rng, n, n);
    idx.iter().map(|&i| data[i].clone()).collect()
}

fn fit_and_canonicalize(
    data: &[Vec<f64>],
    k: usize,
    seed: u64,
    max_iter: usize,
    tol: f64,
) -> CanonicalIca {
    let result = fast_ica(data, k, seed, max_iter, tol);
    canonicalize_ica(&result)
}

/// Run the cross-cohort alignment on already-canonicalized programs,
/// returning `(programs, archetype_labels)` parallel to the input.
/// Uses cosine similarity with `tau` threshold; no reciprocal-best
/// filter (caller handles that upstream if desired).
fn align_canonicals(
    canonicals: &[(String, CanonicalIca)],
    protein_labels: &[String],
    top_n: usize,
    cosine_tau: f64,
) -> (Vec<AlignedProgram>, Vec<usize>) {
    let mut programs: Vec<AlignedProgram> = Vec::new();
    for (cohort, canon) in canonicals {
        for (c, loading) in canon.loadings.iter().enumerate() {
            programs.push(AlignedProgram {
                cohort: cohort.clone(),
                program: format!("program_{:02}", c + 1),
                category: None,
                values: loading.clone(),
            });
        }
    }
    let _ = protein_labels; // kept for a future labelled-output path
    let labels = build_archetypes(
        &programs,
        AlignMetric::Cosine,
        top_n,
        cosine_tau,
        false, // reciprocal_best
        false, // category_constraint
    );
    (programs, labels)
}

/// Group programs by their archetype label and return for each
/// multi-member archetype: (`archetype_id`, cohort set,
/// representative loading vector from the first program in the
/// archetype).
fn group_archetypes(
    programs: &[AlignedProgram],
    labels: &[usize],
) -> Vec<(usize, Vec<String>, Vec<f64>)> {
    use std::collections::BTreeMap;
    let mut by_label: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (idx, lab) in labels.iter().enumerate() {
        by_label.entry(*lab).or_default().push(idx);
    }
    let mut out = Vec::new();
    for (lab, members) in by_label {
        if members.len() < 2 {
            continue;
        }
        let cohorts: Vec<String> = {
            let mut s: Vec<String> =
                members.iter().map(|&i| programs[i].cohort.clone()).collect();
            s.sort();
            s.dedup();
            s
        };
        let representative = programs[members[0]].values.clone();
        out.push((lab, cohorts, representative));
    }
    out
}

pub fn align_bootstrap(
    cohorts: &[CohortMatrix],
    params: BootstrapParams,
) -> Result<Vec<BootstrapRow>, String> {
    if cohorts.len() < 2 {
        return Err("align bootstrap requires at least 2 cohorts".into());
    }
    if params.n_boot == 0 {
        return Err("--n-boot must be >= 1".into());
    }
    // Enforce consistent protein universe across cohorts.
    let universe = &cohorts[0].protein_labels;
    for c in cohorts.iter().skip(1) {
        if &c.protein_labels != universe {
            return Err(format!(
                "protein-label mismatch: cohort {:?} universe differs from first",
                c.label
            ));
        }
    }
    for c in cohorts {
        if c.data.len() < params.min_subjects {
            return Err(format!(
                "cohort {:?} has {} subjects; --min-subjects is {}",
                c.label,
                c.data.len(),
                params.min_subjects
            ));
        }
        if c.data.is_empty() || c.data[0].len() != universe.len() {
            return Err(format!(
                "cohort {:?} matrix shape disagrees with protein universe",
                c.label
            ));
        }
    }

    // Point-estimate fit on unresampled matrices.
    let pe_canons: Vec<(String, CanonicalIca)> = cohorts
        .iter()
        .map(|c| {
            let canon = fit_and_canonicalize(
                &c.data,
                params.k,
                params.seed,
                params.max_iter,
                params.tol,
            );
            (c.label.clone(), canon)
        })
        .collect();
    let (pe_programs, pe_labels) =
        align_canonicals(&pe_canons, universe, params.top_n, params.cosine_tau);
    let pe_archetypes = group_archetypes(&pe_programs, &pe_labels);
    if pe_archetypes.is_empty() {
        return Ok(Vec::new());
    }

    // Per-PE-archetype accumulators: (sum_n_cohorts, n_universal,
    // n_multi, n_cohorts_series).
    let n_total_cohorts = cohorts.len();
    let mut acc: Vec<BootstrapAcc> = pe_archetypes
        .iter()
        .map(|_| BootstrapAcc::new(params.n_boot))
        .collect();

    for iter in 0..params.n_boot {
        let sub_seed = derive_sub_seed(params.seed, iter);
        let mut rng = Xoshiro256pp::new(sub_seed);
        // Resample every cohort using the same rng stream so iterations
        // stay tied to a single sub-seed.
        let bs_canons: Vec<(String, CanonicalIca)> = cohorts
            .iter()
            .enumerate()
            .map(|(ci, c)| {
                let resampled = resample_rows(&c.data, &mut rng);
                let canon = fit_and_canonicalize(
                    &resampled,
                    params.k,
                    sub_seed.wrapping_add(0x1_0000_0000 + ci as u64),
                    params.max_iter,
                    params.tol,
                );
                (c.label.clone(), canon)
            })
            .collect();
        let (bs_programs, bs_labels) =
            align_canonicals(&bs_canons, universe, params.top_n, params.cosine_tau);
        let bs_archetypes = group_archetypes(&bs_programs, &bs_labels);

        // Match every PE archetype to at most one bootstrap archetype
        // by max cosine between representative loadings.
        for (ai, (_, _, pe_rep)) in pe_archetypes.iter().enumerate() {
            let mut best_sim = f64::NEG_INFINITY;
            let mut best_n_cohorts: Option<usize> = None;
            for (_, bs_cohorts, bs_rep) in &bs_archetypes {
                let sim = cosine(pe_rep, bs_rep).unwrap_or(0.0).abs();
                if sim > best_sim {
                    best_sim = sim;
                    best_n_cohorts = Some(bs_cohorts.len());
                }
            }
            if best_sim >= params.match_tau {
                if let Some(n) = best_n_cohorts {
                    acc[ai].observe(n, n_total_cohorts);
                }
            } else {
                acc[ai].observe_miss();
            }
        }
    }

    Ok(pe_archetypes
        .iter()
        .zip(acc.iter())
        .enumerate()
        .map(|(ai, ((_, pe_cohorts, _), a))| a.to_row(ai + 1, pe_cohorts.clone()))
        .collect())
}

struct BootstrapAcc {
    matches: Vec<usize>, // n_cohorts for each successful match
    iterations: usize,
    n_boot: usize,
}

impl BootstrapAcc {
    fn new(n_boot: usize) -> Self {
        Self {
            matches: Vec::with_capacity(n_boot),
            iterations: 0,
            n_boot,
        }
    }
    fn observe(&mut self, n_cohorts: usize, _total_cohorts: usize) {
        self.matches.push(n_cohorts);
        self.iterations += 1;
    }
    fn observe_miss(&mut self) {
        self.iterations += 1;
    }
    fn to_row(&self, archetype_id: usize, observed_cohorts: Vec<String>) -> BootstrapRow {
        let mut sorted = self.matches.clone();
        sorted.sort();
        let mean = if sorted.is_empty() {
            0.0
        } else {
            sorted.iter().sum::<usize>() as f64 / sorted.len() as f64
        };
        let match_rate = if self.iterations == 0 {
            0.0
        } else {
            sorted.len() as f64 / self.iterations as f64
        };
        let prob_universal = if self.n_boot == 0 {
            0.0
        } else {
            let total_cohorts = observed_cohorts.len().max(1);
            sorted.iter().filter(|&&n| n >= total_cohorts).count() as f64
                / self.n_boot as f64
        };
        let prob_multi = if self.n_boot == 0 {
            0.0
        } else {
            sorted.iter().filter(|&&n| n >= 2).count() as f64 / self.n_boot as f64
        };
        let (ci_lower, ci_upper) = if sorted.is_empty() {
            (0, 0)
        } else {
            let n = sorted.len();
            let lower_idx = ((0.025 * n as f64) as usize).min(n - 1);
            let upper_idx = ((0.975 * n as f64).ceil() as usize)
                .saturating_sub(1)
                .min(n - 1);
            (sorted[lower_idx], sorted[upper_idx])
        };
        BootstrapRow {
            archetype_id,
            observed_n_cohorts: observed_cohorts.len(),
            observed_cohorts,
            bootstrap_mean_n_cohorts: mean,
            bootstrap_prob_universal: prob_universal,
            bootstrap_prob_multi: prob_multi,
            ci_lower_n_cohorts: ci_lower,
            ci_upper_n_cohorts: ci_upper,
            bootstrap_match_rate: match_rate,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy_cohort(label: &str, n: usize, p: usize, seed: u64) -> CohortMatrix {
        let mut rng = Xoshiro256pp::new(seed);
        let data: Vec<Vec<f64>> = (0..n)
            .map(|_| (0..p).map(|_| rng.next_normal()).collect())
            .collect();
        CohortMatrix {
            label: label.into(),
            data,
            protein_labels: (0..p).map(|j| format!("P{j:03}")).collect(),
        }
    }

    #[test]
    fn resample_rows_length_matches_input() {
        let m = toy_cohort("A", 20, 10, 11);
        let mut rng = Xoshiro256pp::new(3);
        let r = resample_rows(&m.data, &mut rng);
        assert_eq!(r.len(), m.data.len());
        for row in &r {
            assert_eq!(row.len(), 10);
        }
    }

    #[test]
    fn align_bootstrap_refuses_single_cohort() {
        let a = toy_cohort("A", 20, 10, 1);
        let p = BootstrapParams {
            k: 2,
            n_boot: 1,
            seed: 0,
            top_n: 5,
            cosine_tau: 0.0,
            match_tau: 0.5,
            max_iter: 20,
            tol: 1e-2,
            min_subjects: 5,
        };
        let err = align_bootstrap(&[a], p).unwrap_err();
        assert!(err.contains("at least 2 cohorts"));
    }

    #[test]
    fn align_bootstrap_refuses_small_cohorts() {
        let a = toy_cohort("A", 3, 10, 1);
        let b = toy_cohort("B", 3, 10, 2);
        let p = BootstrapParams {
            k: 2,
            n_boot: 1,
            seed: 0,
            top_n: 5,
            cosine_tau: 0.0,
            match_tau: 0.5,
            max_iter: 20,
            tol: 1e-2,
            min_subjects: 5,
        };
        let err = align_bootstrap(&[a, b], p).unwrap_err();
        assert!(err.contains("min-subjects"), "unexpected: {err}");
    }

    #[test]
    fn align_bootstrap_refuses_mismatched_universes() {
        let mut a = toy_cohort("A", 20, 10, 1);
        let b = toy_cohort("B", 20, 10, 2);
        a.protein_labels[0] = "DIFFERENT".into();
        let p = BootstrapParams {
            k: 2,
            n_boot: 1,
            seed: 0,
            top_n: 5,
            cosine_tau: 0.0,
            match_tau: 0.5,
            max_iter: 20,
            tol: 1e-2,
            min_subjects: 5,
        };
        let err = align_bootstrap(&[a, b], p).unwrap_err();
        assert!(err.contains("mismatch"), "unexpected: {err}");
    }

    #[test]
    fn align_bootstrap_runs_and_returns_finite_stats_on_toy_cohorts() {
        let a = toy_cohort("A", 20, 15, 1);
        let b = toy_cohort("B", 20, 15, 2);
        let p = BootstrapParams {
            k: 2,
            n_boot: 5,
            seed: 20260418,
            top_n: 5,
            cosine_tau: 0.0,
            match_tau: 0.0,
            max_iter: 40,
            tol: 1e-2,
            min_subjects: 5,
        };
        let rows = align_bootstrap(&[a, b], p).unwrap();
        for r in &rows {
            assert!(r.bootstrap_mean_n_cohorts >= 0.0);
            assert!(r.bootstrap_prob_universal >= 0.0 && r.bootstrap_prob_universal <= 1.0);
            assert!(r.bootstrap_prob_multi >= 0.0 && r.bootstrap_prob_multi <= 1.0);
            assert!(r.bootstrap_match_rate >= 0.0 && r.bootstrap_match_rate <= 1.0);
        }
    }

    #[test]
    fn align_bootstrap_is_deterministic() {
        let a = toy_cohort("A", 20, 15, 1);
        let b = toy_cohort("B", 20, 15, 2);
        let p = BootstrapParams {
            k: 2,
            n_boot: 4,
            seed: 20260418,
            top_n: 5,
            cosine_tau: 0.0,
            match_tau: 0.0,
            max_iter: 40,
            tol: 1e-2,
            min_subjects: 5,
        };
        let x = align_bootstrap(&[a.clone(), b.clone()], p).unwrap();
        let y = align_bootstrap(&[a, b], p).unwrap();
        assert_eq!(x.len(), y.len());
        for (p, q) in x.iter().zip(y.iter()) {
            assert_eq!(p, q);
        }
    }
}
