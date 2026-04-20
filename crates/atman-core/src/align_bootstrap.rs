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
//! fraction in at least two cohorts (`prob_multi`), a percentile CI
//! over the cohort-count distribution, Shannon entropy over the
//! bootstrap n_cohorts histogram (`alignment_entropy`, higher = more
//! uncertain), and a BCa (bias-corrected accelerated) CI whose
//! acceleration is estimated by pooled subject-level jackknife.
//!
//! Matching is via representative-program best-cosine; each
//! point-estimate archetype is represented by the loading vector
//! from its first-cohort program, and every bootstrap/jackknife
//! archetype is represented likewise. Match is the archetype whose
//! representative has maximum absolute cosine similarity against the
//! point-estimate representative, with a configurable floor
//! (`match_tau`). Misses contribute zero-cohort observations to the
//! distribution (they are a valid outcome under resampling noise).
//!
//! Determinism: bootstrap randomness comes from
//! [`crate::ica::Xoshiro256pp`] seeded per iteration via a
//! SplitMix64 derivation over `(seed, iter)`. Jackknife is
//! deterministic — the single ICA seed per cohort is reused.

use crate::align::{build_archetypes, AlignMetric, AlignedProgram};
use crate::ica::{canonicalize_ica, fast_ica, CanonicalIca, Xoshiro256pp};
use crate::stats::cosine;
use statrs::distribution::{ContinuousCDF, Normal};

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
    /// Shannon entropy (in bits) of the empirical distribution of
    /// `n_cohorts` across bootstrap iterations, treating each missed
    /// match as `n_cohorts = 0`. Low entropy ⇒ the archetype's
    /// cohort coverage is stable under subject resampling.
    pub alignment_entropy: f64,
    /// BCa (bias-corrected accelerated) bootstrap CI lower bound at
    /// the two-sided 95% level. Acceleration is estimated by
    /// pooled-subject jackknife. When the BCa denominator is
    /// non-positive (non-monotone tail) the CI falls back to the
    /// percentile CI and `bca_fallback_to_percentile` is true.
    pub bca_lower_n_cohorts: f64,
    pub bca_upper_n_cohorts: f64,
    pub bca_fallback_to_percentile: bool,
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

/// Match a point-estimate archetype to the jackknife/bootstrap
/// archetype whose representative has maximum absolute cosine
/// similarity. Returns `None` when the best similarity is below
/// `match_tau` or the resampled pipeline produced no archetypes.
fn match_pe_archetype(
    pe_rep: &[f64],
    bs_archetypes: &[(usize, Vec<String>, Vec<f64>)],
    match_tau: f64,
) -> Option<usize> {
    let mut best_sim = f64::NEG_INFINITY;
    let mut best_n: Option<usize> = None;
    for (_, bs_cohorts, bs_rep) in bs_archetypes {
        let sim = cosine(pe_rep, bs_rep).unwrap_or(0.0).abs();
        if sim > best_sim {
            best_sim = sim;
            best_n = Some(bs_cohorts.len());
        }
    }
    if best_sim >= match_tau { best_n } else { None }
}

/// Shannon entropy (in bits) of the empirical histogram of a
/// non-negative integer-valued series, treating each distinct value
/// as a category. Returns 0 on an empty series.
fn shannon_entropy_bits(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    use std::collections::BTreeMap;
    let mut counts: BTreeMap<i64, usize> = BTreeMap::new();
    for &v in values {
        *counts.entry(v.round() as i64).or_insert(0) += 1;
    }
    let n = values.len() as f64;
    let mut h = 0.0;
    for &c in counts.values() {
        if c == 0 {
            continue;
        }
        let p = c as f64 / n;
        h -= p * p.log2();
    }
    h
}

/// BCa (bias-corrected accelerated) CI at confidence level
/// `1 - alpha`. `bootstrap` is the B-vector of resampled θ values;
/// `jackknife` is the leave-one-subject-out vector; `theta_hat` is
/// the point estimate. Returns `(lower, upper, fallback_to_percentile)`.
fn bca_ci(
    bootstrap: &[f64],
    jackknife: &[f64],
    theta_hat: f64,
    alpha: f64,
) -> (f64, f64, bool) {
    let b = bootstrap.len();
    if b == 0 {
        return (theta_hat, theta_hat, true);
    }
    let mut sorted = bootstrap.to_vec();
    sorted.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let normal = match Normal::new(0.0, 1.0) {
        Ok(n) => n,
        Err(_) => return (sorted[0], sorted[b - 1], true),
    };
    // Percentile fallback values, computed once.
    let pct_lo = sorted[((alpha / 2.0) * b as f64) as usize];
    let pct_hi = sorted[(((1.0 - alpha / 2.0) * b as f64) as usize).min(b - 1)];

    // z0: bias correction from the fraction of bootstrap < θ̂.
    let below = bootstrap.iter().filter(|&&v| v < theta_hat).count();
    let prop = (below as f64 / b as f64).clamp(
        1.0 / (b as f64 + 1.0),
        1.0 - 1.0 / (b as f64 + 1.0),
    );
    let z0 = normal.inverse_cdf(prop);

    // Acceleration via jackknife. Needs ≥ 2 jackknife values and
    // positive spread; otherwise fall back to bias-corrected-only
    // by setting a = 0.
    let a = if jackknife.len() >= 2 {
        let jbar = jackknife.iter().sum::<f64>() / jackknife.len() as f64;
        let num: f64 = jackknife.iter().map(|&t| (jbar - t).powi(3)).sum();
        let den_sq: f64 = jackknife.iter().map(|&t| (jbar - t).powi(2)).sum();
        if den_sq > 0.0 {
            num / (6.0 * den_sq.powf(1.5))
        } else {
            0.0
        }
    } else {
        0.0
    };

    let z_lo = normal.inverse_cdf(alpha / 2.0);
    let z_hi = normal.inverse_cdf(1.0 - alpha / 2.0);
    let denom_lo = 1.0 - a * (z0 + z_lo);
    let denom_hi = 1.0 - a * (z0 + z_hi);
    if !denom_lo.is_finite() || !denom_hi.is_finite() || denom_lo <= 0.0 || denom_hi <= 0.0 {
        return (pct_lo, pct_hi, true);
    }
    let alpha1 = normal.cdf(z0 + (z0 + z_lo) / denom_lo);
    let alpha2 = normal.cdf(z0 + (z0 + z_hi) / denom_hi);
    if !alpha1.is_finite() || !alpha2.is_finite() {
        return (pct_lo, pct_hi, true);
    }
    let lo_idx = ((alpha1 * b as f64) as usize).min(b - 1);
    let hi_idx = ((alpha2 * b as f64) as usize).min(b - 1);
    (sorted[lo_idx], sorted[hi_idx], false)
}

/// Run the point-estimate ICA + alignment pipeline on a set of
/// cohort matrices. Returns, for each input PE archetype, the
/// matched jackknife archetype's cohort count (or 0 on miss).
fn point_estimate_match_counts(
    cohorts: &[CohortMatrix],
    params: &BootstrapParams,
    universe: &[String],
    pe_archetypes: &[(usize, Vec<String>, Vec<f64>)],
) -> Vec<f64> {
    let canons: Vec<(String, CanonicalIca)> = cohorts
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
    let (progs, labels) =
        align_canonicals(&canons, universe, params.top_n, params.cosine_tau);
    let archetypes = group_archetypes(&progs, &labels);
    pe_archetypes
        .iter()
        .map(|(_, _, pe_rep)| {
            match_pe_archetype(pe_rep, &archetypes, params.match_tau)
                .map(|n| n as f64)
                .unwrap_or(0.0)
        })
        .collect()
}

/// Subject-level jackknife over the pooled cohort space. For each
/// subject across all cohorts, drop them and re-run the point-
/// estimate pipeline, recording the matched-archetype cohort count
/// per PE archetype. Returns a `pe_archetypes.len() × N_total` matrix
/// (outer = PE archetype, inner = jackknife replicate).
fn jackknife_n_cohorts(
    cohorts: &[CohortMatrix],
    params: &BootstrapParams,
    universe: &[String],
    pe_archetypes: &[(usize, Vec<String>, Vec<f64>)],
) -> Vec<Vec<f64>> {
    let mut jack: Vec<Vec<f64>> = pe_archetypes.iter().map(|_| Vec::new()).collect();
    for (ci, cohort) in cohorts.iter().enumerate() {
        // Skip cohorts below min_subjects + 1 (leave-one-out would
        // take the cohort below the published minimum). Contributing
        // zero jackknife replicates from such a cohort is honest —
        // its effective acceleration contribution was going to be
        // unreliable anyway.
        if cohort.data.len() <= params.min_subjects {
            continue;
        }
        for drop_idx in 0..cohort.data.len() {
            let reduced: Vec<CohortMatrix> = cohorts
                .iter()
                .enumerate()
                .map(|(cj, c)| {
                    let data = if cj == ci {
                        c.data
                            .iter()
                            .enumerate()
                            .filter(|&(k, _)| k != drop_idx)
                            .map(|(_, r)| r.clone())
                            .collect()
                    } else {
                        c.data.clone()
                    };
                    CohortMatrix {
                        label: c.label.clone(),
                        data,
                        protein_labels: c.protein_labels.clone(),
                    }
                })
                .collect();
            let counts = point_estimate_match_counts(&reduced, params, universe, pe_archetypes);
            for (ai, n) in counts.iter().enumerate() {
                jack[ai].push(*n);
            }
        }
    }
    jack
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
            match match_pe_archetype(pe_rep, &bs_archetypes, params.match_tau) {
                Some(n) => acc[ai].observe(n, n_total_cohorts),
                None => acc[ai].observe_miss(),
            }
        }
    }

    // Subject-level jackknife for BCa acceleration.
    let jack = jackknife_n_cohorts(cohorts, &params, universe, &pe_archetypes);

    Ok(pe_archetypes
        .iter()
        .zip(acc.iter())
        .zip(jack.iter())
        .enumerate()
        .map(|(ai, (((_, pe_cohorts, _), a), jack_row))| {
            a.to_row(ai + 1, pe_cohorts.clone(), jack_row)
        })
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
    fn to_row(
        &self,
        archetype_id: usize,
        observed_cohorts: Vec<String>,
        jackknife: &[f64],
    ) -> BootstrapRow {
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

        // Full n_boot-length bootstrap series for entropy + BCa:
        // misses count as n_cohorts = 0 (they are a valid outcome
        // under resampling noise, not a missing observation).
        let mut full: Vec<f64> = self.matches.iter().map(|&n| n as f64).collect();
        full.resize(self.n_boot, 0.0);
        let entropy = shannon_entropy_bits(&full);
        let theta_hat = observed_cohorts.len() as f64;
        let (bca_lo, bca_hi, bca_fallback) = bca_ci(&full, jackknife, theta_hat, 0.05);

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
            alignment_entropy: entropy,
            bca_lower_n_cohorts: bca_lo,
            bca_upper_n_cohorts: bca_hi,
            bca_fallback_to_percentile: bca_fallback,
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
            assert!(r.alignment_entropy.is_finite() && r.alignment_entropy >= 0.0);
            assert!(r.bca_lower_n_cohorts.is_finite());
            assert!(r.bca_upper_n_cohorts.is_finite());
            assert!(r.bca_lower_n_cohorts <= r.bca_upper_n_cohorts);
        }
    }

    #[test]
    fn shannon_entropy_bits_uniform_equals_log2_k() {
        let v: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0];
        // Uniform over 4 categories ⇒ H = log2(4) = 2 bits.
        assert!((shannon_entropy_bits(&v) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn shannon_entropy_bits_degenerate_is_zero() {
        let v: Vec<f64> = vec![1.0; 8];
        assert!(shannon_entropy_bits(&v).abs() < 1e-12);
    }

    #[test]
    fn bca_ci_symmetric_data_approximates_percentile() {
        // Symmetric distribution ⇒ z0 ≈ 0; skew-free jackknife ⇒ a ≈ 0.
        // BCa then reduces to the percentile CI.
        let bootstrap: Vec<f64> = (0..1001).map(|i| i as f64 / 1000.0).collect();
        let jackknife: Vec<f64> = (0..100).map(|i| 0.5 + (i as f64 - 49.5) * 0.001).collect();
        let (lo, hi, fallback) = bca_ci(&bootstrap, &jackknife, 0.5, 0.05);
        assert!(!fallback);
        assert!((lo - 0.025).abs() < 0.02, "lo={lo}");
        assert!((hi - 0.975).abs() < 0.02, "hi={hi}");
    }

    #[test]
    fn bca_ci_zero_spread_jackknife_reduces_to_bias_corrected() {
        // All jackknife values equal ⇒ den_sq = 0 ⇒ a = 0. BCa still
        // runs with bias correction only; no fallback.
        let bootstrap: Vec<f64> = (0..1001).map(|i| i as f64 / 1000.0).collect();
        let jackknife: Vec<f64> = vec![0.5; 50];
        let (lo, hi, fallback) = bca_ci(&bootstrap, &jackknife, 0.5, 0.05);
        assert!(!fallback);
        assert!(lo < hi);
    }

    #[test]
    fn bca_ci_empty_bootstrap_falls_back() {
        let (lo, hi, fallback) = bca_ci(&[], &[1.0, 2.0], 1.5, 0.05);
        assert!(fallback);
        assert_eq!(lo, 1.5);
        assert_eq!(hi, 1.5);
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
