//! Subject-level bootstrap for cross-cohort archetype alignment.
//!
//! For each bootstrap iteration, resample subjects within every
//! cohort with replacement, re-run the configured [`Decomposition`]
//! (FastICA or single-seed NMF) per cohort on the resampled matrix,
//! run the existing cosine-similarity archetype alignment, and
//! record — for each point-estimate archetype — which cohorts the
//! matched bootstrap archetype was recovered in.
//!
//! Output per point-estimate archetype, over TWO DIFFERENT
//! POPULATIONS. Nothing in the column names says which, so it is
//! stated here and on every field of [`BootstrapRow`].
//!
//! Over MATCHED REPLICATES ONLY, meaning the iterations where this
//! archetype was matched at all: the mean cohort count
//! (`bootstrap_mean_n_cohorts`) and the percentile CI
//! (`ci_lower_n_cohorts`, `ci_upper_n_cohorts`).
//!
//! Over ALL `n_boot` ITERATIONS, with every miss entered as zero
//! cohorts: the fraction recovered in every cohort
//! (`prob_universal`), the fraction in at least two cohorts
//! (`prob_multi`), the Shannon entropy of the n_cohorts histogram
//! (`alignment_entropy`, higher = more uncertain), and a BCa
//! (bias-corrected accelerated) CI whose acceleration is estimated by
//! pooled subject-level jackknife.
//!
//! `bootstrap_match_rate` is the fraction that matched at all, and it
//! is the bridge between the two. A percentile CI read without it
//! looks far stronger than the archetype is: one matched twice in 200
//! resamples, spanning six cohorts both times, reports [6, 6]. Two
//! archetypes differing four-fold in reproducibility can carry the
//! identical interval.
//!
//! Matching is via representative-program best-cosine; each
//! point-estimate archetype is represented by the loading vector
//! from its first-cohort program, and every bootstrap/jackknife
//! archetype is represented likewise. Match is the archetype whose
//! representative has maximum absolute cosine similarity against the
//! point-estimate representative, with a configurable floor
//! (`match_tau`). A miss is a valid outcome under resampling noise
//! rather than a missing observation, so it enters `prob_universal`,
//! `prob_multi`, `alignment_entropy` and the BCa CI as a zero-cohort
//! observation. It is ABSENT from the mean and the percentile CI,
//! which is what makes those two conditional on recovery.
//!
//! Determinism: bootstrap randomness comes from
//! [`crate::ica::Xoshiro256pp`] seeded per iteration via a
//! SplitMix64 derivation over `(seed, iter)`. Jackknife is
//! deterministic — the single decomposition seed per cohort is
//! reused.

use crate::align::{build_archetypes, similarity, AlignMetric, AlignedProgram};
use crate::ica::{canonicalize_ica, fast_ica, CanonicalIca, Xoshiro256pp};
use crate::nmf::{apply_transform, nmf, BetaLoss, Init, NmfConfig, Transform as NmfTransform};
use crate::rng::derive_sub_seed;
use rayon::prelude::*;
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

/// Per-resample decomposition method for [`align_bootstrap`]. Chosen once
/// per invocation and applied identically to the point estimate, every
/// bootstrap resample, and every jackknife replicate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Decomposition {
    /// FastICA (the original, and still default, method). `max_iter` /
    /// `tol` are FastICA's fixed-point iteration knobs.
    Ica { max_iter: usize, tol: f64 },
    /// Single-seed multiplicative-updates NMF. `transform` is
    /// re-applied to a fresh copy of the per-resample matrix immediately
    /// before decomposition — never cached across resamples, the same
    /// recompute-per-call rule `align project`'s `shift-min` transform
    /// uses. With `Transform::None`, a negative entry is rejected with a
    /// clear error rather than tripping `atman_core::nmf::nmf`'s internal
    /// non-negativity assertion.
    Nmf {
        beta_loss: BetaLoss,
        init: Init,
        max_iter: usize,
        tol: f64,
        transform: NmfTransform,
    },
}

#[derive(Debug, Clone, Copy)]
pub struct BootstrapParams {
    pub k: usize,
    pub n_boot: usize,
    pub seed: u64,
    pub top_n: usize,
    pub cosine_tau: f64,
    /// Floor on similarity between the point-estimate representative
    /// and a bootstrap/jackknife program for the archetype containing
    /// that program to be considered a match. Units follow `metric`
    /// (cosine ∈ [0,1], Jaccard ∈ [0,1], |Spearman| ∈ [0,1]).
    pub match_tau: f64,
    pub decomposition: Decomposition,
    /// Bootstrap iterations abort with an error when any cohort has
    /// fewer than this many distinct subjects.
    pub min_subjects: usize,
    /// Similarity metric used both to align bootstrap programs into
    /// archetypes and to match bootstrap/jackknife archetypes back to
    /// the point-estimate archetypes. Cosine (default) is fastest and
    /// sign-invariant; Jaccard is set-based over the top-`top_n`
    /// loadings by `|value|`; Spearman uses the absolute rank
    /// correlation.
    pub metric: AlignMetric,
    /// Two-sided CI alpha. Both the percentile CI on
    /// `n_cohorts` and the BCa CI use `alpha/2` and `1 - alpha/2`
    /// quantiles. `0.05` ⇒ 95% CI, and it must lie in `(0, 1)`.
    ///
    /// The two intervals are taken over different populations. The
    /// percentile CI uses matched replicates only. The BCa CI uses all
    /// `n_boot` iterations, with each miss entered as zero. See
    /// [`BootstrapRow::ci_lower_n_cohorts`].
    pub ci_alpha: f64,
    /// Worker threads for the bootstrap and jackknife loops; `0` means
    /// one per available core. Output does not depend on this value:
    /// every iteration draws from its own `derive_sub_seed(seed, iter)`
    /// stream and results are folded in iteration order.
    pub threads: usize,
}

/// One point-estimate archetype's bootstrap summary.
///
/// The fields split over two populations and each one says which.
/// Read `bootstrap_match_rate` before any conditional field.
#[derive(Debug, Clone, PartialEq)]
pub struct BootstrapRow {
    pub archetype_id: usize,
    pub observed_n_cohorts: usize,
    pub observed_cohorts: Vec<String>,
    /// Mean cohort count OVER MATCHED REPLICATES ONLY. This is the
    /// archetype's mean span given that it was recovered, and it says
    /// nothing about how often it was recovered.
    pub bootstrap_mean_n_cohorts: f64,
    /// Fraction of ALL `n_boot` iterations where the matched
    /// archetype spanned every cohort. Misses count against it.
    pub bootstrap_prob_universal: f64,
    /// Fraction of ALL `n_boot` iterations where the matched
    /// archetype spanned at least two cohorts. Misses count against
    /// it.
    pub bootstrap_prob_multi: f64,
    /// Percentile CI at the `alpha/2` and `1 - alpha/2` quantiles,
    /// where `alpha` comes from `BootstrapParams::ci_alpha`,
    /// OVER MATCHED REPLICATES ONLY.
    ///
    /// The interval is CONDITIONAL ON RECOVERY. It gives the span of
    /// the archetype in the resamples that recovered it, not the
    /// probability that a resample recovers it. An archetype matched
    /// twice in 200 resamples, spanning six cohorts both times,
    /// reports `[6, 6]` — indistinguishable here from one matched in
    /// every resample.
    ///
    /// A criterion of the form `ci_lower_n_cohorts >= k` is therefore
    /// close to untestable, because a rarely matched archetype passes
    /// it. Gate on `bootstrap_match_rate` or
    /// `bootstrap_prob_universal` instead, which are unconditional.
    pub ci_lower_n_cohorts: usize,
    /// Upper percentile bound. Also OVER MATCHED REPLICATES ONLY —
    /// see [`BootstrapRow::ci_lower_n_cohorts`].
    pub ci_upper_n_cohorts: usize,
    /// Fraction of bootstrap iterations where **any** archetype was
    /// matched to this point-estimate archetype at all.
    ///
    /// This is the denominator the conditional fields are missing.
    /// When it is low, `bootstrap_mean_n_cohorts` and the percentile
    /// CI rest on few replicates and overstate the archetype.
    /// `align bootstrap` warns on stderr when any archetype falls
    /// below 0.50.
    pub bootstrap_match_rate: f64,
    /// Shannon entropy (in bits) of the empirical distribution of
    /// `n_cohorts` over ALL `n_boot` ITERATIONS, treating each missed
    /// match as `n_cohorts = 0`.
    ///
    /// NOT MONOTONE IN REPRODUCIBILITY. When an archetype spans the
    /// same cohorts whenever it is matched, the histogram is two-point
    /// and this reduces to the binary entropy of
    /// `bootstrap_match_rate`. It peaks at a match rate of 0.5 and
    /// falls to zero at BOTH ends. A barely recovered archetype
    /// therefore scores LOW — measured at `n_boot = 200`, a match rate
    /// of 0.025 gives 0.169 against 0.732 at 0.205, so the less
    /// reproducible archetype looks the more stable one.
    ///
    /// Low entropy means the cohort count is repeatable. It means the
    /// archetype is stable only when read together with a high
    /// `bootstrap_match_rate`. Do not rank on it alone.
    pub alignment_entropy: f64,
    /// BCa (bias-corrected accelerated) bootstrap CI lower bound at
    /// the level set by `BootstrapParams::ci_alpha`, over ALL `n_boot`
    /// ITERATIONS with each miss entered as zero. Unlike the
    /// percentile CI it is NOT conditional on recovery, so across the
    /// middle of the range a low match rate pulls its lower bound to
    /// zero.
    ///
    /// DEGENERATE AT BOTH TAILS. `n_cohorts` is a heavily tied
    /// two-point vector, so `z0` is driven by the miss fraction alone
    /// and beyond roughly 0.97 in either direction both endpoints run
    /// off the same end. Measured at `n_boot = 200`, `ci_alpha = 0.05`:
    /// the lower bound reads 6.0 at a match rate of 0.025, and the
    /// interval collapses to [0, 0] at 0.99. `bca_fallback_to_percentile`
    /// does not fire, because it tests the acceleration denominator and
    /// not this. See
    /// `tests::bca_endpoints_are_degenerate_at_extreme_match_rates`.
    /// Acceleration is
    /// estimated by pooled-subject jackknife. When the BCa
    /// denominator is non-positive (non-monotone tail) the CI falls
    /// back to the percentile CI and `bca_fallback_to_percentile` is
    /// true.
    pub bca_lower_n_cohorts: f64,
    pub bca_upper_n_cohorts: f64,
    pub bca_fallback_to_percentile: bool,
}

/// Sample `n_samples` indices with replacement from `0..n` using the
/// reviewed unbiased bounded draw on the raw u64 stream.
fn bootstrap_indices(rng: &mut Xoshiro256pp, n: usize, n_samples: usize) -> Vec<usize> {
    (0..n_samples).map(|_| rng.bounded(n)).collect()
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

/// Reject a decomposition's loadings if any entry is NaN or ±∞, rather
/// than letting it flow silently into cosine similarity, archetype
/// grouping, and the bootstrap/BCa accumulators downstream. Degenerate
/// resamples (duplicate rows, zero-variance columns) can in principle
/// drive an iterative decomposition to diverge; this is the shared
/// loud-failure gate for that, applied identically to both
/// [`Decomposition`] variants. `context` names the cohort and call site
/// (point estimate / bootstrap iteration / jackknife replicate) so the
/// error is debuggable.
fn check_finite_loadings(loadings: Vec<Vec<f64>>, context: &str) -> Result<Vec<Vec<f64>>, String> {
    if loadings.iter().flatten().any(|v| !v.is_finite()) {
        return Err(format!(
            "align bootstrap: decomposition produced non-finite loadings for {context}; \
             this usually means a degenerate resample (duplicate rows / zero-variance \
             column) drove the fit to diverge. Consider a smaller --k, a different --seed, \
             or (for --decomposition nmf) --transform exp2-clip to bound the input scale."
        ));
    }
    Ok(loadings)
}

/// Fit the configured [`Decomposition`] on `data` (n × p, rows =
/// samples, columns = the shared protein universe) and return its
/// `k × p` feature loadings.
///
/// This is the single dispatch point used for the point estimate,
/// every bootstrap resample, and every jackknife replicate — the
/// same per-call behavior at all three sites. For [`Decomposition::Nmf`],
/// `transform` is applied to a fresh copy of `data` right before
/// decomposition (see [`Decomposition::Nmf`] for why it is never
/// cached across resamples). Both branches route their result through
/// [`check_finite_loadings`] before returning. `context` names the
/// cohort and call site, used only in the (rare) error messages.
fn fit_loadings(
    data: &[Vec<f64>],
    k: usize,
    seed: u64,
    decomposition: &Decomposition,
    context: &str,
) -> Result<Vec<Vec<f64>>, String> {
    match decomposition {
        Decomposition::Ica { max_iter, tol } => check_finite_loadings(
            fit_and_canonicalize(data, k, seed, *max_iter, *tol).loadings,
            context,
        ),
        Decomposition::Nmf {
            beta_loss,
            init,
            max_iter,
            tol,
            transform,
        } => {
            let mut transformed = data.to_vec();
            apply_transform(&mut transformed, transform);
            if matches!(transform, NmfTransform::None) {
                for row in &transformed {
                    for &v in row {
                        if v < 0.0 {
                            return Err(format!(
                                "align bootstrap: NMF (--transform none) requires \
                                 non-negative input, but found value {v} for {context}. \
                                 Apply a non-negativity-preserving transform (--transform \
                                 exp2-clip or --transform shift-min) before running \
                                 align bootstrap --decomposition nmf."
                            ));
                        }
                    }
                }
            }
            let cfg = NmfConfig {
                k,
                beta_loss: *beta_loss,
                init: *init,
                max_iter: *max_iter,
                tol: *tol,
                seed,
            };
            check_finite_loadings(nmf(&transformed, &cfg).h, context)
        }
    }
}

/// Run the cross-cohort alignment on per-cohort loadings, returning
/// `(programs, archetype_labels)` parallel to the input. Uses `metric`
/// similarity with `tau` threshold; no reciprocal-best filter (caller
/// handles that upstream if desired).
fn align_program_loadings(
    per_cohort: &[(String, Vec<Vec<f64>>)],
    protein_labels: &[String],
    top_n: usize,
    cosine_tau: f64,
    metric: AlignMetric,
) -> (Vec<AlignedProgram>, Vec<usize>) {
    let mut programs: Vec<AlignedProgram> = Vec::new();
    for (cohort, loadings) in per_cohort {
        for (c, loading) in loadings.iter().enumerate() {
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
        &programs, metric, top_n, cosine_tau, false, // reciprocal_best
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
            let mut s: Vec<String> = members
                .iter()
                .map(|&i| programs[i].cohort.clone())
                .collect();
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
/// archetype whose representative has maximum `metric` similarity
/// (`|cosine|` / Jaccard-top-N / `|Spearman|`, all in `[0, 1]`).
/// Returns `None` when the best similarity is below `match_tau` or
/// the resampled pipeline produced no archetypes.
fn match_pe_archetype(
    pe_rep: &[f64],
    bs_archetypes: &[(usize, Vec<String>, Vec<f64>)],
    match_tau: f64,
    metric: AlignMetric,
    top_n: usize,
) -> Option<usize> {
    let mut best_sim = f64::NEG_INFINITY;
    let mut best_n: Option<usize> = None;
    for (_, bs_cohorts, bs_rep) in bs_archetypes {
        let sim = similarity(pe_rep, bs_rep, metric, top_n);
        if sim > best_sim {
            best_sim = sim;
            best_n = Some(bs_cohorts.len());
        }
    }
    if best_sim >= match_tau {
        best_n
    } else {
        None
    }
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
fn bca_ci(bootstrap: &[f64], jackknife: &[f64], theta_hat: f64, alpha: f64) -> (f64, f64, bool) {
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
    let prop =
        (below as f64 / b as f64).clamp(1.0 / (b as f64 + 1.0), 1.0 - 1.0 / (b as f64 + 1.0));
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
    context: &str,
) -> Result<Vec<f64>, String> {
    let loadings: Vec<(String, Vec<Vec<f64>>)> = cohorts
        .iter()
        .map(|c| {
            fit_loadings(
                &c.data,
                params.k,
                params.seed,
                &params.decomposition,
                &format!("cohort {:?} ({context})", c.label),
            )
            .map(|l| (c.label.clone(), l))
        })
        .collect::<Result<Vec<_>, String>>()?;
    let (progs, labels) = align_program_loadings(
        &loadings,
        universe,
        params.top_n,
        params.cosine_tau,
        params.metric,
    );
    let archetypes = group_archetypes(&progs, &labels);
    Ok(pe_archetypes
        .iter()
        .map(|(_, _, pe_rep)| {
            match_pe_archetype(
                pe_rep,
                &archetypes,
                params.match_tau,
                params.metric,
                params.top_n,
            )
            .map(|n| n as f64)
            .unwrap_or(0.0)
        })
        .collect())
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
    pool: &rayon::ThreadPool,
) -> Result<Vec<Vec<f64>>, String> {
    // Replicates in the serial order (cohort-major, dropped subject
    // minor). Cohorts at or below min_subjects are skipped: leave-one-out
    // would take them below the published minimum, and contributing
    // zero replicates from such a cohort is honest — its acceleration
    // contribution was going to be unreliable anyway.
    let replicates: Vec<(usize, usize)> = cohorts
        .iter()
        .enumerate()
        .filter(|(_, c)| c.data.len() > params.min_subjects)
        .flat_map(|(ci, c)| (0..c.data.len()).map(move |drop_idx| (ci, drop_idx)))
        .collect();
    let per_replicate: Vec<Result<Vec<f64>, String>> = pool.install(|| {
        replicates
            .par_iter()
            .map(|&(ci, drop_idx)| {
                jackknife_replicate(cohorts, params, universe, pe_archetypes, ci, drop_idx)
            })
            .collect()
    });
    let mut jack: Vec<Vec<f64>> = pe_archetypes.iter().map(|_| Vec::new()).collect();
    for outcome in per_replicate {
        for (ai, n) in outcome?.iter().enumerate() {
            jack[ai].push(*n);
        }
    }
    Ok(jack)
}

/// One leave-one-subject-out replicate: drop subject `drop_idx` from
/// cohort `ci`, refit every cohort, and return the per-PE-archetype
/// matched cohort counts.
fn jackknife_replicate(
    cohorts: &[CohortMatrix],
    params: &BootstrapParams,
    universe: &[String],
    pe_archetypes: &[(usize, Vec<String>, Vec<f64>)],
    ci: usize,
    drop_idx: usize,
) -> Result<Vec<f64>, String> {
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
    point_estimate_match_counts(
        &reduced,
        params,
        universe,
        pe_archetypes,
        &format!(
            "jackknife drop cohort {:?} subject #{drop_idx}",
            cohorts[ci].label
        ),
    )
}

/// One bootstrap iteration: resample every cohort from the iteration's
/// own sub-seed stream, refit, align, group, and match each point-
/// estimate archetype. Returns, per PE archetype, `Some(n_cohorts)` of
/// the matched bootstrap archetype or `None` on a miss.
fn bootstrap_iteration(
    cohorts: &[CohortMatrix],
    params: &BootstrapParams,
    universe: &[String],
    pe_archetypes: &[(usize, Vec<String>, Vec<f64>)],
    iter: usize,
) -> Result<Vec<Option<usize>>, String> {
    let sub_seed = derive_sub_seed(params.seed, iter);
    let mut rng = Xoshiro256pp::new(sub_seed);
    // Resample every cohort using the same rng stream so iterations
    // stay tied to a single sub-seed.
    let bs_loadings: Vec<(String, Vec<Vec<f64>>)> = cohorts
        .iter()
        .enumerate()
        .map(|(ci, c)| {
            let resampled = resample_rows(&c.data, &mut rng);
            let cohort_seed = sub_seed.wrapping_add(0x1_0000_0000 + ci as u64);
            fit_loadings(
                &resampled,
                params.k,
                cohort_seed,
                &params.decomposition,
                &format!("cohort {:?} (bootstrap iter {iter})", c.label),
            )
            .map(|l| (c.label.clone(), l))
        })
        .collect::<Result<Vec<_>, String>>()?;
    let (bs_programs, bs_labels) = align_program_loadings(
        &bs_loadings,
        universe,
        params.top_n,
        params.cosine_tau,
        params.metric,
    );
    let bs_archetypes = group_archetypes(&bs_programs, &bs_labels);

    // Match every PE archetype to at most one bootstrap archetype
    // by max cosine between representative loadings.
    Ok(pe_archetypes
        .iter()
        .map(|(_, _, pe_rep)| {
            match_pe_archetype(
                pe_rep,
                &bs_archetypes,
                params.match_tau,
                params.metric,
                params.top_n,
            )
        })
        .collect())
}

/// How selective the two tau gates are on the point-estimate program
/// vectors they operate over.
///
/// `align bootstrap` applies `metric` twice: `cosine_tau` groups
/// bootstrap programs into archetypes, and `match_tau` matches
/// bootstrap and jackknife archetypes back to the point estimate. Both
/// run on the same loading vectors, so a gate that admits nearly
/// everything is deciding recurrence without thresholding — which is
/// the bootstrap's entire job. Non-negative loadings hit this because
/// cosine cannot go below zero on them.
#[derive(Debug, Clone, Copy)]
pub struct TauDiagnostics {
    pub n_pairs: usize,
    pub n_above_cosine_tau: usize,
    pub n_above_match_tau: usize,
    pub median_similarity: f64,
}

/// [`align_bootstrap`] plus the selectivity of its two gates.
pub fn align_bootstrap_with_diagnostics(
    cohorts: &[CohortMatrix],
    params: BootstrapParams,
) -> Result<(Vec<BootstrapRow>, Option<TauDiagnostics>), String> {
    let diagnostics = std::cell::RefCell::new(None);
    let rows = align_bootstrap_inner(cohorts, params, Some(&diagnostics))?;
    Ok((rows, diagnostics.into_inner()))
}

pub fn align_bootstrap(
    cohorts: &[CohortMatrix],
    params: BootstrapParams,
) -> Result<Vec<BootstrapRow>, String> {
    align_bootstrap_inner(cohorts, params, None)
}

fn align_bootstrap_inner(
    cohorts: &[CohortMatrix],
    params: BootstrapParams,
    diagnostics: Option<&std::cell::RefCell<Option<TauDiagnostics>>>,
) -> Result<Vec<BootstrapRow>, String> {
    if cohorts.len() < 2 {
        return Err("align bootstrap requires at least 2 cohorts".into());
    }
    if params.n_boot == 0 {
        return Err("--n-boot must be >= 1".into());
    }
    if !(params.ci_alpha.is_finite() && params.ci_alpha > 0.0 && params.ci_alpha < 1.0) {
        return Err(format!(
            "--ci-alpha must lie in (0, 1); got {}",
            params.ci_alpha
        ));
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
    let pe_loadings: Vec<(String, Vec<Vec<f64>>)> = cohorts
        .iter()
        .map(|c| {
            fit_loadings(
                &c.data,
                params.k,
                params.seed,
                &params.decomposition,
                &format!("cohort {:?} (point estimate)", c.label),
            )
            .map(|l| (c.label.clone(), l))
        })
        .collect::<Result<Vec<_>, String>>()?;
    let (pe_programs, pe_labels) = align_program_loadings(
        &pe_loadings,
        universe,
        params.top_n,
        params.cosine_tau,
        params.metric,
    );
    // Measure what the two gates actually admit on these vectors,
    // before using either of them to decide anything.
    if let Some(sink) = diagnostics {
        let mut sims: Vec<f64> = Vec::new();
        for i in 0..pe_programs.len() {
            for j in (i + 1)..pe_programs.len() {
                if pe_programs[i].cohort == pe_programs[j].cohort {
                    continue; // within-cohort pairs are never gated
                }
                let v = crate::align::similarity(
                    &pe_programs[i].values,
                    &pe_programs[j].values,
                    params.metric,
                    params.top_n,
                );
                if v.is_finite() {
                    sims.push(v);
                }
            }
        }
        if !sims.is_empty() {
            let mut sorted = sims.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            *sink.borrow_mut() = Some(TauDiagnostics {
                n_pairs: sims.len(),
                n_above_cosine_tau: sims.iter().filter(|v| **v >= params.cosine_tau).count(),
                n_above_match_tau: sims.iter().filter(|v| **v >= params.match_tau).count(),
                median_similarity: sorted[sorted.len() / 2],
            });
        }
    }
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

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(params.threads)
        .build()
        .map_err(|e| {
            format!(
                "align bootstrap: building thread pool ({} threads): {e}",
                params.threads
            )
        })?;

    // Iterations are independent given their sub-seeds, so they run as a
    // parallel map; the fold below walks the results in iteration order,
    // which keeps every accumulator (and the first reported error)
    // identical to the serial loop for any thread count.
    let per_iter: Vec<Result<Vec<Option<usize>>, String>> = pool.install(|| {
        (0..params.n_boot)
            .into_par_iter()
            .map(|iter| bootstrap_iteration(cohorts, &params, universe, &pe_archetypes, iter))
            .collect()
    });
    for outcome in per_iter {
        for (ai, matched) in outcome?.into_iter().enumerate() {
            match matched {
                Some(n) => acc[ai].observe(n, n_total_cohorts),
                None => acc[ai].observe_miss(),
            }
        }
    }

    // Subject-level jackknife for BCa acceleration.
    let jack = jackknife_n_cohorts(cohorts, &params, universe, &pe_archetypes, &pool)?;

    Ok(pe_archetypes
        .iter()
        .zip(acc.iter())
        .zip(jack.iter())
        .enumerate()
        .map(|(ai, (((_, pe_cohorts, _), a), jack_row))| {
            a.to_row(ai + 1, pe_cohorts.clone(), jack_row, params.ci_alpha)
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
        ci_alpha: f64,
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
            sorted.iter().filter(|&&n| n >= total_cohorts).count() as f64 / self.n_boot as f64
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
            let lower_idx = (((ci_alpha / 2.0) * n as f64) as usize).min(n - 1);
            let upper_idx = (((1.0 - ci_alpha / 2.0) * n as f64).ceil() as usize)
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
        let (bca_lo, bca_hi, bca_fallback) = bca_ci(&full, jackknife, theta_hat, ci_alpha);

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

    /// Builds an accumulator matched in `m` of 200 resamples, always
    /// spanning all six cohorts when matched. This is the shape the
    /// GBM manuscript session hit, and it isolates the match rate as
    /// the only thing that varies.
    fn acc_matched(m: usize) -> BootstrapRow {
        let mut acc = BootstrapAcc::new(200);
        for _ in 0..m {
            acc.observe(6, 6);
        }
        for _ in 0..(200 - m) {
            acc.observe_miss();
        }
        let cohorts: Vec<String> = (1..=6).map(|i| format!("c{i}")).collect();
        acc.to_row(1, cohorts, &[], 0.05)
    }

    /// Pins which population each summary field is over.
    ///
    /// The percentile CI and the mean use matched replicates only,
    /// while `prob_universal`, `prob_multi` and the BCa bound use all
    /// `n_boot`. Reported by the GBM manuscript session: an NMF
    /// archetype at match rate 0.205 and an ICA archetype at 0.805
    /// carried identical [6, 6] intervals — a four-fold difference in
    /// reproducibility that the interval could not show.
    #[test]
    fn percentile_ci_conditions_on_recovery_but_prob_universal_does_not() {
        let nmf = acc_matched(41); // 0.205
        let ica = acc_matched(161); // 0.805

        // Conditional on recovery: identical, and blind to the gap.
        assert_eq!(
            (nmf.ci_lower_n_cohorts, nmf.ci_upper_n_cohorts),
            (ica.ci_lower_n_cohorts, ica.ci_upper_n_cohorts),
            "the percentile CI cannot separate a 0.205 from a 0.805 \
             match rate, which is why it must not be gated on"
        );
        assert_eq!((nmf.ci_lower_n_cohorts, nmf.ci_upper_n_cohorts), (6, 6));
        assert!((nmf.bootstrap_mean_n_cohorts - 6.0).abs() < 1e-12);
        assert!((ica.bootstrap_mean_n_cohorts - 6.0).abs() < 1e-12);

        // Unconditional: these carry the four-fold difference.
        assert!((nmf.bootstrap_match_rate - 0.205).abs() < 1e-12);
        assert!((ica.bootstrap_match_rate - 0.805).abs() < 1e-12);
        assert!((nmf.bootstrap_prob_universal - 0.205).abs() < 1e-12);
        assert!((ica.bootstrap_prob_universal - 0.805).abs() < 1e-12);
    }

    /// An archetype matched in every resample carries the identical
    /// percentile interval, which is the point: the conditional
    /// columns cannot tell it from a rarely matched one.
    #[test]
    fn an_always_matched_archetype_carries_the_same_percentile_ci() {
        let rare = acc_matched(5); // 0.025
        let always = acc_matched(200); // 1.000

        assert_eq!(
            (rare.ci_lower_n_cohorts, rare.ci_upper_n_cohorts),
            (always.ci_lower_n_cohorts, always.ci_upper_n_cohorts)
        );
        assert!((rare.bootstrap_match_rate - 0.025).abs() < 1e-12);
        assert!((always.bootstrap_match_rate - 1.0).abs() < 1e-12);
    }

    /// `alignment_entropy` is NOT monotone in reproducibility.
    ///
    /// When an archetype spans every cohort whenever it is matched,
    /// the n_cohorts histogram is two-point over {0, 6} and the
    /// entropy is the binary entropy of the match rate. It peaks at
    /// 0.5 and falls to zero at BOTH ends, so a barely recovered
    /// archetype scores lower — looks more stable — than a moderately
    /// recovered one. Never rank on entropy alone.
    #[test]
    fn alignment_entropy_is_not_monotone_in_reproducibility() {
        let barely = acc_matched(5); // 0.025
        let reported = acc_matched(41); // 0.205
        let even = acc_matched(100); // 0.500
        let nearly_always = acc_matched(199); // 0.995

        assert!(even.alignment_entropy > reported.alignment_entropy);
        assert!(
            barely.alignment_entropy < reported.alignment_entropy,
            "an archetype recovered in 2.5% of resamples reports LOWER \
             entropy ({:.3}) than one recovered in 20.5% ({:.3})",
            barely.alignment_entropy,
            reported.alignment_entropy
        );
        assert!(nearly_always.alignment_entropy < even.alignment_entropy);
        // The two ends are indistinguishable by entropy alone.
        assert!(
            (barely.alignment_entropy - nearly_always.alignment_entropy).abs() < 0.13,
            "0.025 and 0.995 match rates give near-identical entropy: \
             {:.3} vs {:.3}",
            barely.alignment_entropy,
            nearly_always.alignment_entropy
        );
    }

    /// Pins a KNOWN DEFECT rather than a desired behaviour.
    ///
    /// The BCa endpoints are degenerate at both tails of the match
    /// rate. `n_cohorts` here is a heavily tied two-point vector, so
    /// the bias correction `z0` is driven by the miss fraction alone.
    /// Beyond roughly 0.97 in either direction both endpoints run off
    /// the same end of the sorted vector. Measured at `n_boot = 200`,
    /// `ci_alpha = 0.05`: the lower bound reads 6.0 at a match rate of
    /// 0.025, and the interval collapses to [0, 0] at 0.99 — an
    /// archetype recovered almost every time.
    ///
    /// `bca_fallback_to_percentile` does not fire, because it tests
    /// the acceleration denominator and not this. Recorded so a change
    /// to the BCa endpoints is a deliberate act with a visible
    /// diff, not a silent one.
    #[test]
    fn bca_endpoints_are_degenerate_at_extreme_match_rates() {
        let barely = acc_matched(5); // 0.025
        assert!((barely.bca_lower_n_cohorts - 6.0).abs() < 1e-9);
        assert!(!barely.bca_fallback_to_percentile);

        let nearly_always = acc_matched(198); // 0.990
        assert!((nearly_always.bca_lower_n_cohorts - 0.0).abs() < 1e-9);
        assert!((nearly_always.bca_upper_n_cohorts - 0.0).abs() < 1e-9);
        assert!(!nearly_always.bca_fallback_to_percentile);

        // In between it behaves: misses drag the lower bound to zero.
        for m in [10usize, 41, 100, 161, 190] {
            let r = acc_matched(m);
            assert!(
                (r.bca_lower_n_cohorts - 0.0).abs() < 1e-9,
                "match rate {} gave BCa lower {}",
                r.bootstrap_match_rate,
                r.bca_lower_n_cohorts
            );
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
            decomposition: Decomposition::Ica {
                max_iter: 20,
                tol: 1e-2,
            },
            min_subjects: 5,
            metric: AlignMetric::Cosine,
            ci_alpha: 0.05,
            threads: 1,
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
            decomposition: Decomposition::Ica {
                max_iter: 20,
                tol: 1e-2,
            },
            min_subjects: 5,
            metric: AlignMetric::Cosine,
            ci_alpha: 0.05,
            threads: 1,
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
            decomposition: Decomposition::Ica {
                max_iter: 20,
                tol: 1e-2,
            },
            min_subjects: 5,
            metric: AlignMetric::Cosine,
            ci_alpha: 0.05,
            threads: 1,
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
            decomposition: Decomposition::Ica {
                max_iter: 40,
                tol: 1e-2,
            },
            min_subjects: 5,
            metric: AlignMetric::Cosine,
            ci_alpha: 0.05,
            threads: 1,
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

    /// The bootstrap and jackknife loops are embarrassingly parallel
    /// (one sub-seed per iteration). Whatever the thread count, the
    /// fold must happen in iteration order so the rows are identical.
    #[test]
    fn align_bootstrap_rows_are_identical_across_thread_counts() {
        let cohorts = [
            toy_cohort("A", 24, 12, 1),
            toy_cohort("B", 22, 12, 2),
            toy_cohort("C", 26, 12, 3),
        ];
        let mk = |threads: usize| BootstrapParams {
            k: 2,
            n_boot: 12,
            seed: 20260906,
            top_n: 5,
            cosine_tau: 0.0,
            match_tau: 0.0,
            decomposition: Decomposition::Ica {
                max_iter: 40,
                tol: 1e-2,
            },
            min_subjects: 5,
            metric: AlignMetric::Cosine,
            ci_alpha: 0.05,
            threads,
        };
        let serial = align_bootstrap(&cohorts, mk(1)).unwrap();
        assert!(!serial.is_empty(), "toy cohorts should yield archetypes");
        let three = align_bootstrap(&cohorts, mk(3)).unwrap();
        let all_cores = align_bootstrap(&cohorts, mk(0)).unwrap();
        assert_eq!(serial, three, "threads=3 must reproduce threads=1 exactly");
        assert_eq!(
            serial, all_cores,
            "threads=0 (all cores) must reproduce threads=1 exactly"
        );
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
    fn align_bootstrap_accepts_every_metric_without_changing_row_count() {
        let a = toy_cohort("A", 20, 15, 1);
        let b = toy_cohort("B", 20, 15, 2);
        let base = BootstrapParams {
            k: 2,
            n_boot: 3,
            seed: 20260418,
            top_n: 5,
            cosine_tau: 0.0,
            match_tau: 0.0,
            decomposition: Decomposition::Ica {
                max_iter: 40,
                tol: 1e-2,
            },
            min_subjects: 5,
            metric: AlignMetric::Cosine,
            ci_alpha: 0.05,
            threads: 1,
        };
        let cos_rows = align_bootstrap(&[a.clone(), b.clone()], base).unwrap();
        let mut jac = base;
        jac.metric = AlignMetric::Jaccard;
        let jac_rows = align_bootstrap(&[a.clone(), b.clone()], jac).unwrap();
        let mut spr = base;
        spr.metric = AlignMetric::Spearman;
        let spr_rows = align_bootstrap(&[a, b], spr).unwrap();
        // Different metric ⇒ different archetype alignment ⇒ row
        // count can differ across metrics, but each must produce
        // valid, finite-stat rows.
        for rows in [&cos_rows, &jac_rows, &spr_rows] {
            for r in rows.iter() {
                assert!(r.alignment_entropy.is_finite());
                assert!(r.bootstrap_prob_multi >= 0.0 && r.bootstrap_prob_multi <= 1.0);
            }
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
            decomposition: Decomposition::Ica {
                max_iter: 40,
                tol: 1e-2,
            },
            min_subjects: 5,
            metric: AlignMetric::Cosine,
            ci_alpha: 0.05,
            threads: 1,
        };
        let x = align_bootstrap(&[a.clone(), b.clone()], p).unwrap();
        let y = align_bootstrap(&[a, b], p).unwrap();
        assert_eq!(x.len(), y.len());
        for (p, q) in x.iter().zip(y.iter()) {
            assert_eq!(p, q);
        }
    }

    #[test]
    fn check_finite_loadings_rejects_nan() {
        let bad = vec![vec![1.0, f64::NAN], vec![2.0, 3.0]];
        let err = check_finite_loadings(bad, "cohort \"A\" (test)").unwrap_err();
        assert!(err.contains("non-finite"), "unexpected: {err}");
        assert!(err.contains("cohort \"A\" (test)"), "unexpected: {err}");
    }

    #[test]
    fn check_finite_loadings_rejects_infinite() {
        let bad = vec![vec![f64::INFINITY, 0.0]];
        let err = check_finite_loadings(bad, "cohort \"B\" (test)").unwrap_err();
        assert!(err.contains("non-finite"), "unexpected: {err}");
    }

    #[test]
    fn check_finite_loadings_rejects_negative_infinity() {
        let bad = vec![vec![0.0, f64::NEG_INFINITY]];
        let err = check_finite_loadings(bad, "cohort \"C\" (test)").unwrap_err();
        assert!(err.contains("non-finite"), "unexpected: {err}");
    }

    #[test]
    fn check_finite_loadings_accepts_finite_matrix() {
        let good = vec![vec![1.0, 2.0], vec![-3.5, 0.0]];
        let out = check_finite_loadings(good.clone(), "cohort \"D\" (test)").unwrap();
        assert_eq!(out, good);
    }
}
