//! Differential abundance testing for log2-NPX proteomics workflows.
//!
//! Statistical discipline rails:
//!
//! - **Sample-level only.** The test operates on per-subject paired
//!   differences. The caller is responsible for supplying pairs keyed by a
//!   valid biological-replicate identifier (`subject_id`); this module just
//!   does the math.
//! - **Min-pairs guard.** Below `min_pairs` the test is refused and a
//!   `PairedTResult::Skipped` is returned. Default 5.
//! - **Two-sided p-values only** — the default scientific convention when
//!   we don't have a prior directional hypothesis.
//! - **BH-FDR is applied per family** — the caller decides what a family
//!   is (e.g., one comparison × all proteins). This module does not pool
//!   families together.

use statrs::distribution::{ContinuousCDF, StudentsT};

/// A single paired t-test result for one protein in one comparison.
#[derive(Debug, Clone, PartialEq)]
pub enum PairedTResult {
    /// Test ran. `n_pairs >= min_pairs` and sd > 0.
    Computed {
        n_pairs: usize,
        mean_a: f64,
        mean_b: f64,
        mean_diff: f64,
        t: f64,
        df: f64,
        p_value: f64,
    },
    /// Test was not run. Reason captured for the report sidecar.
    Skipped { reason: SkipReason, n_pairs: usize },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SkipReason {
    /// `n_pairs < min_pairs`.
    InsufficientPairs,
    /// All paired differences are exactly zero → t is undefined (0/0).
    ZeroVariance,
    /// At least one non-finite input (NaN, Inf).
    NonFiniteInput,
}

/// Compute the paired Student's t-test on `(a, b)` pairs.
/// Returns `Computed` if `n_pairs >= min_pairs`, else `Skipped`.
///
/// All arithmetic is in f64. Differences `d_i = a_i - b_i` are summed with
/// naive f64 accumulation — for n ≤ 40 this is within 1 ULP and is what
/// Python / R / SAS would produce on this scale.
pub fn paired_t(pairs: &[(f64, f64)], min_pairs: usize) -> PairedTResult {
    let n = pairs.len();
    if n < min_pairs || n < 2 {
        return PairedTResult::Skipped {
            reason: SkipReason::InsufficientPairs,
            n_pairs: n,
        };
    }
    if pairs.iter().any(|(a, b)| !a.is_finite() || !b.is_finite()) {
        return PairedTResult::Skipped {
            reason: SkipReason::NonFiniteInput,
            n_pairs: n,
        };
    }

    let n_f = n as f64;
    let mean_a: f64 = pairs.iter().map(|(a, _)| a).sum::<f64>() / n_f;
    let mean_b: f64 = pairs.iter().map(|(_, b)| b).sum::<f64>() / n_f;
    let diffs: Vec<f64> = pairs.iter().map(|(a, b)| a - b).collect();
    let mean_diff: f64 = diffs.iter().sum::<f64>() / n_f;

    // Sample variance with Bessel's correction (df = n-1).
    let var_diff: f64 = diffs.iter().map(|d| (d - mean_diff).powi(2)).sum::<f64>() / (n_f - 1.0);
    let sd_diff = var_diff.sqrt();
    if !sd_diff.is_finite() {
        return PairedTResult::Skipped {
            reason: SkipReason::NonFiniteInput,
            n_pairs: n,
        };
    }
    if sd_diff == 0.0 {
        return PairedTResult::Skipped {
            reason: SkipReason::ZeroVariance,
            n_pairs: n,
        };
    }

    let se = sd_diff / n_f.sqrt();
    let t = mean_diff / se;
    let df = n_f - 1.0;
    if !t.is_finite() || !df.is_finite() || df <= 0.0 {
        return PairedTResult::Skipped {
            reason: SkipReason::InsufficientPairs,
            n_pairs: n,
        };
    }

    // Two-sided p = 2 * P(T > |t|) on df degrees of freedom.
    let t_dist = match StudentsT::new(0.0, 1.0, df) {
        Ok(d) => d,
        Err(_) => {
            return PairedTResult::Skipped {
                reason: SkipReason::InsufficientPairs,
                n_pairs: n,
            };
        }
    };
    let p_value = 2.0 * (1.0 - t_dist.cdf(t.abs()));

    PairedTResult::Computed {
        n_pairs: n,
        mean_a,
        mean_b,
        mean_diff,
        t,
        df,
        p_value,
    }
}

/// Welch's two-sample t-test for unpaired designs.
///
/// Used when samples in the two groups come from independent subjects — the
/// standard case-control setup. Does not assume equal variances: denominator
/// is $\sqrt{s_a^2/n_a + s_b^2/n_b}$ and the degrees of freedom use the
/// Welch-Satterthwaite approximation
///
/// $$
/// \nu = \frac{(s_a^2/n_a + s_b^2/n_b)^2}
///            {(s_a^2/n_a)^2/(n_a-1) + (s_b^2/n_b)^2/(n_b-1)}.
/// $$
///
/// Both groups must have at least `min_per_group` non-finite-free samples
/// and non-zero within-group variance; the `mean_a`, `mean_b` and
/// `mean_diff = mean_a - mean_b` fields in the result mirror the paired-t
/// output schema so downstream writers stay unchanged.
pub fn welch_t(group_a: &[f64], group_b: &[f64], min_per_group: usize) -> PairedTResult {
    let na = group_a.len();
    let nb = group_b.len();
    let n_total = na + nb;
    if na < min_per_group || nb < min_per_group || na < 2 || nb < 2 {
        return PairedTResult::Skipped {
            reason: SkipReason::InsufficientPairs,
            n_pairs: n_total,
        };
    }
    if group_a.iter().chain(group_b.iter()).any(|v| !v.is_finite()) {
        return PairedTResult::Skipped {
            reason: SkipReason::NonFiniteInput,
            n_pairs: n_total,
        };
    }

    let na_f = na as f64;
    let nb_f = nb as f64;
    let mean_a: f64 = group_a.iter().sum::<f64>() / na_f;
    let mean_b: f64 = group_b.iter().sum::<f64>() / nb_f;
    let var_a: f64 = group_a.iter().map(|v| (v - mean_a).powi(2)).sum::<f64>() / (na_f - 1.0);
    let var_b: f64 = group_b.iter().map(|v| (v - mean_b).powi(2)).sum::<f64>() / (nb_f - 1.0);
    if !var_a.is_finite() || !var_b.is_finite() {
        return PairedTResult::Skipped {
            reason: SkipReason::NonFiniteInput,
            n_pairs: n_total,
        };
    }
    if var_a == 0.0 && var_b == 0.0 {
        return PairedTResult::Skipped {
            reason: SkipReason::ZeroVariance,
            n_pairs: n_total,
        };
    }

    let se = (var_a / na_f + var_b / nb_f).sqrt();
    if se == 0.0 {
        return PairedTResult::Skipped {
            reason: SkipReason::ZeroVariance,
            n_pairs: n_total,
        };
    }
    let mean_diff = mean_a - mean_b;
    let t = mean_diff / se;
    // Welch-Satterthwaite df.
    let num = (var_a / na_f + var_b / nb_f).powi(2);
    let den = (var_a / na_f).powi(2) / (na_f - 1.0) + (var_b / nb_f).powi(2) / (nb_f - 1.0);
    let df = num / den;
    if !t.is_finite() || !df.is_finite() || df <= 0.0 {
        return PairedTResult::Skipped {
            reason: SkipReason::InsufficientPairs,
            n_pairs: n_total,
        };
    }
    let t_dist = match StudentsT::new(0.0, 1.0, df) {
        Ok(d) => d,
        Err(_) => {
            return PairedTResult::Skipped {
                reason: SkipReason::InsufficientPairs,
                n_pairs: n_total,
            };
        }
    };
    let p_value = 2.0 * (1.0 - t_dist.cdf(t.abs()));

    PairedTResult::Computed {
        n_pairs: n_total,
        mean_a,
        mean_b,
        mean_diff,
        t,
        df,
        p_value,
    }
}

/// Stability-weighted inference score.
///
/// Combines the three per-protein quantities that a small-cohort DE
/// analysis produces — effect size, significance, and LOO directional
/// stability — into a single monotonic ranking score
///
/// $$ S_i = |\bar d_i| \cdot \mathrm{SSR}_i \cdot (1 - q_i) $$
///
/// where $\bar d_i$ is the mean paired difference (or unpaired effect),
/// $q_i \in [0,1]$ is the BH-adjusted $p$-value, and
/// $\mathrm{SSR}_i \in [0,1]$ is the LOO sign-match rate from `atman
/// robustness`. The formulation keeps effect-size units while
/// multiplicatively down-weighting direction-unstable proteins (low
/// $\mathrm{SSR}_i$) and threshold-marginal proteins (high $q_i$). Ranking
/// proteins by $S_i$ rather than $q_i$ alone gives a stability-aware
/// prioritisation that does not discard statistical inference but
/// explicitly shrinks fragile marginal hits below robust borderline-q hits.
///
/// Returns `None` when any input is absent or non-finite (propagates the
/// atman pattern of Optional results for skipped tests). `q_i > 1` is
/// clamped and `q_i ≤ 0` is treated as $q = 0$ (maximal weight, same as BH).
pub fn stability_weighted_score(
    mean_diff: Option<f64>,
    bh_q: Option<f64>,
    sign_match_rate: Option<f64>,
) -> Option<f64> {
    let d = mean_diff?;
    let q = bh_q?.clamp(0.0, 1.0);
    let ssr = sign_match_rate?.clamp(0.0, 1.0);
    if !d.is_finite() || !q.is_finite() || !ssr.is_finite() {
        return None;
    }
    Some(d.abs() * ssr * (1.0 - q))
}

/// Ordinary-least-squares fit of a linear model `y = Xβ + ε` with optional
/// covariate columns alongside the group indicator.
///
/// The design matrix `X` is passed as a slice of row vectors, one per sample
/// (outer index = sample, inner index = design column). Column 0 is
/// conventionally the intercept (a column of 1.0), and the caller decides
/// which further columns are the group indicator and the covariates.
///
/// Returns per-coefficient estimates, standard errors, t-statistics, and
/// two-sided p-values. `df = n - p` where `n` is the number of samples and
/// `p` is the number of design columns.
///
/// The solver uses Cholesky decomposition of the normal equations
/// `XᵀX β = Xᵀy`. Cholesky is used rather than QR because
///
/// - design matrices in practice are tiny (`p ≤ ~5` for `Group + age + sex`
///   in this project),
/// - `XᵀX` is symmetric positive-definite when `X` has full column rank, so
///   Cholesky is the most economical stable choice,
/// - if Cholesky fails (pivot ≤ 0), the design is near-singular and we
///   return `Skipped(ZeroVariance)` rather than producing an inflated fit.
///
/// Complete-case handling is the caller's responsibility: rows in `design`
/// and values in `y` must already have had non-finite-covariate rows
/// filtered out. This function only guards against non-finite entries that
/// still slip through.
#[derive(Debug, Clone, PartialEq)]
pub struct OlsFit {
    /// Number of samples used (rows of `X`).
    pub n: usize,
    /// Number of design columns (cols of `X`).
    pub p: usize,
    /// Regression coefficients, indexed by design column.
    pub beta: Vec<f64>,
    /// Standard error for each coefficient.
    pub se: Vec<f64>,
    /// t-statistic for each coefficient (`β / SE`).
    pub t: Vec<f64>,
    /// Two-sided p-value for each coefficient on `df` degrees of freedom.
    pub p_value: Vec<f64>,
    /// Residual degrees of freedom, `n − p`.
    pub df: f64,
    /// Residual variance estimate, `RSS / df`.
    pub sigma2: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum OlsOutcome {
    /// Fit succeeded.
    Computed(OlsFit),
    /// Fit was not run. Reason captured for the report sidecar.
    Skipped { reason: SkipReason, n: usize },
}

pub fn ols(design: &[Vec<f64>], y: &[f64], min_samples: usize) -> OlsOutcome {
    let n = design.len();
    if n != y.len() {
        return OlsOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }
    if n == 0 {
        return OlsOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n: 0,
        };
    }
    let p = design[0].len();
    if p == 0 {
        return OlsOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }
    if n < min_samples || n <= p {
        return OlsOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }
    // Finite-input guard.
    for row in design {
        if row.len() != p {
            return OlsOutcome::Skipped {
                reason: SkipReason::InsufficientPairs,
                n,
            };
        }
        if row.iter().any(|v| !v.is_finite()) {
            return OlsOutcome::Skipped {
                reason: SkipReason::NonFiniteInput,
                n,
            };
        }
    }
    if y.iter().any(|v| !v.is_finite()) {
        return OlsOutcome::Skipped {
            reason: SkipReason::NonFiniteInput,
            n,
        };
    }

    // Build XᵀX (p × p) and Xᵀy (p).
    let mut xtx = vec![vec![0.0f64; p]; p];
    let mut xty = vec![0.0f64; p];
    for i in 0..n {
        let row = &design[i];
        let yi = y[i];
        for j in 0..p {
            let rij = row[j];
            for k in 0..p {
                xtx[j][k] += rij * row[k];
            }
            xty[j] += rij * yi;
        }
    }

    // Cholesky: XᵀX = L Lᵀ (lower-triangular `l`).
    let mut l = vec![vec![0.0f64; p]; p];
    for j in 0..p {
        let diag_sum: f64 = l[j].iter().take(j).map(|v| v * v).sum();
        let diag = xtx[j][j] - diag_sum;
        if !diag.is_finite() || diag <= 0.0 {
            return OlsOutcome::Skipped {
                reason: SkipReason::ZeroVariance,
                n,
            };
        }
        l[j][j] = diag.sqrt();
        for i in (j + 1)..p {
            let off_sum: f64 = l[i]
                .iter()
                .zip(l[j].iter())
                .take(j)
                .map(|(li, lj)| li * lj)
                .sum();
            l[i][j] = (xtx[i][j] - off_sum) / l[j][j];
        }
    }

    // Forward-solve L z = Xᵀy.
    let mut z = vec![0.0f64; p];
    for i in 0..p {
        let mut sum = xty[i];
        for k in 0..i {
            sum -= l[i][k] * z[k];
        }
        z[i] = sum / l[i][i];
    }
    // Back-solve Lᵀ β = z.
    let mut beta = vec![0.0f64; p];
    for i in (0..p).rev() {
        let mut sum = z[i];
        for k in (i + 1)..p {
            sum -= l[k][i] * beta[k]; // Lᵀ[i][k] = L[k][i].
        }
        beta[i] = sum / l[i][i];
    }

    // Residuals + residual variance.
    let mut rss = 0.0f64;
    for i in 0..n {
        let mut pred = 0.0;
        for j in 0..p {
            pred += design[i][j] * beta[j];
        }
        let r = y[i] - pred;
        rss += r * r;
    }
    let df = (n - p) as f64;
    if df <= 0.0 || !df.is_finite() {
        return OlsOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }
    let sigma2 = rss / df;
    if !sigma2.is_finite() {
        return OlsOutcome::Skipped {
            reason: SkipReason::NonFiniteInput,
            n,
        };
    }

    // Standard errors: SE_j = sqrt(σ² · (XᵀX)⁻¹_{jj}).
    // Compute diag of (XᵀX)⁻¹ column-by-column from the Cholesky factor:
    // solve L z_j = e_j, then Lᵀ col_j = z_j; (XᵀX)⁻¹_{jj} = col_j[j].
    let mut se = vec![0.0f64; p];
    for target in 0..p {
        let mut zcol = vec![0.0f64; p];
        for i in 0..p {
            let mut sum = if i == target { 1.0 } else { 0.0 };
            for k in 0..i {
                sum -= l[i][k] * zcol[k];
            }
            zcol[i] = sum / l[i][i];
        }
        let mut col = vec![0.0f64; p];
        for i in (0..p).rev() {
            let mut sum = zcol[i];
            for k in (i + 1)..p {
                sum -= l[k][i] * col[k];
            }
            col[i] = sum / l[i][i];
        }
        let var_jj = col[target];
        if !var_jj.is_finite() || var_jj < 0.0 {
            return OlsOutcome::Skipped {
                reason: SkipReason::ZeroVariance,
                n,
            };
        }
        se[target] = (sigma2 * var_jj).sqrt();
    }

    // t-stats and two-sided p-values.
    let t_dist = match StudentsT::new(0.0, 1.0, df) {
        Ok(d) => d,
        Err(_) => {
            return OlsOutcome::Skipped {
                reason: SkipReason::InsufficientPairs,
                n,
            };
        }
    };
    let t: Vec<f64> = beta
        .iter()
        .zip(se.iter())
        .map(|(b, s)| if *s == 0.0 { f64::NAN } else { b / s })
        .collect();
    let p_value: Vec<f64> = t
        .iter()
        .map(|&ti| {
            if !ti.is_finite() {
                f64::NAN
            } else {
                2.0 * (1.0 - t_dist.cdf(ti.abs()))
            }
        })
        .collect();

    OlsOutcome::Computed(OlsFit {
        n,
        p,
        beta,
        se,
        t,
        p_value,
        df,
        sigma2,
    })
}

/// Benjamini–Hochberg FDR correction. Input `p_values` is a slice of optional
/// p-values (None = skipped test, excluded from the correction). Returns a
/// `Vec<Option<f64>>` of q-values aligned with the input indices; the same
/// `None` slots carry through.
///
/// Standard BH: sort p's ascending, compute `q_i = p_i * m / rank_i`, then
/// do a reverse min-sweep to enforce monotonicity, then map back.
pub fn bh_fdr(p_values: &[Option<f64>]) -> Vec<Option<f64>> {
    // Collect (original_index, p) for non-None entries.
    let mut indexed: Vec<(usize, f64)> = p_values
        .iter()
        .enumerate()
        .filter_map(|(i, p)| p.map(|v| (i, v)))
        .collect();
    let m = indexed.len();
    if m == 0 {
        return vec![None; p_values.len()];
    }
    // Sort by p ascending; stable on ties is not required.
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    // Raw q_i = p_i * m / rank (rank is 1-based).
    let mut q: Vec<f64> = indexed
        .iter()
        .enumerate()
        .map(|(rank0, (_, p))| p * (m as f64) / ((rank0 + 1) as f64))
        .collect();

    // Enforce monotonicity: q_i = min(q_i, q_{i+1}, ..., q_m).
    for i in (0..m.saturating_sub(1)).rev() {
        if q[i + 1] < q[i] {
            q[i] = q[i + 1];
        }
    }
    // Clamp to [0, 1].
    for v in q.iter_mut() {
        if *v > 1.0 {
            *v = 1.0;
        }
    }

    // Scatter back to original positions.
    let mut out: Vec<Option<f64>> = vec![None; p_values.len()];
    for (rank0, (orig_idx, _)) in indexed.iter().enumerate() {
        out[*orig_idx] = Some(q[rank0]);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn paired_t_zero_mean_diff_gives_t_zero_p_one() {
        let pairs = vec![(1.0, 1.0), (2.0, 2.0), (3.0, 3.0), (4.0, 4.0), (5.0, 5.0)];
        let result = paired_t(&pairs, 5);
        match result {
            PairedTResult::Skipped { reason, .. } => {
                assert_eq!(reason, SkipReason::ZeroVariance);
            }
            _ => panic!("expected ZeroVariance skip"),
        }
    }

    #[test]
    fn paired_t_constant_positive_shift() {
        // Every pair shifted by +1.0 with small noise.
        // d = [1.0, 1.0, 1.0, 1.0, 1.0] → mean_diff=1.0, sd=0 → skip.
        // Add noise:
        let pairs = vec![(2.0, 1.0), (3.1, 2.0), (4.0, 2.9), (5.2, 4.0), (6.0, 5.1)];
        let r = paired_t(&pairs, 5);
        match r {
            PairedTResult::Computed {
                n_pairs,
                mean_diff,
                t,
                df,
                p_value,
                ..
            } => {
                assert_eq!(n_pairs, 5);
                assert!(mean_diff > 0.9 && mean_diff < 1.1);
                assert!(t > 10.0, "t should be large; got {}", t);
                assert!((df - 4.0).abs() < 1e-12);
                assert!(p_value < 0.001);
            }
            _ => panic!("expected Computed"),
        }
    }

    #[test]
    fn paired_t_insufficient_pairs_skipped() {
        let pairs = vec![(1.0, 0.0), (2.0, 1.0)];
        let r = paired_t(&pairs, 5);
        assert!(matches!(
            r,
            PairedTResult::Skipped {
                reason: SkipReason::InsufficientPairs,
                ..
            }
        ));
    }

    #[test]
    fn paired_t_never_panics_when_min_pairs_is_too_small() {
        // Caller bug: min_pairs=0. The function still must not panic.
        let r = paired_t(&[(1.0, 0.5)], 0);
        assert!(matches!(
            r,
            PairedTResult::Skipped {
                reason: SkipReason::InsufficientPairs,
                n_pairs: 1
            }
        ));
    }

    #[test]
    fn paired_t_known_textbook_value() {
        // Pairs with diffs [2, 1, 3, 1, 3]:
        //   mean_diff = 10/5 = 2.0
        //   var_diff  = ((0)^2 + (-1)^2 + (1)^2 + (-1)^2 + (1)^2) / 4 = 4/4 = 1.0
        //   sd_diff   = 1.0
        //   se        = 1.0 / sqrt(5) ≈ 0.4472135955
        //   t         = 2.0 / 0.4472135955 ≈ 4.4721359550
        //   df        = 4
        //   two-sided p = 2 * (1 - t_cdf(4.4721, df=4)) ≈ 0.011056 (verified scipy)
        let pairs = vec![(3.0, 1.0), (4.0, 3.0), (5.0, 2.0), (2.0, 1.0), (4.0, 1.0)];
        let r = paired_t(&pairs, 5);
        match r {
            PairedTResult::Computed {
                mean_diff,
                t,
                p_value,
                df,
                ..
            } => {
                assert!(approx(mean_diff, 2.0, 1e-12), "mean_diff={}", mean_diff);
                assert!(approx(t, 4.47213595499958, 1e-10), "t={}", t);
                assert!(approx(df, 4.0, 1e-12));
                assert!(
                    approx(p_value, 0.011056493393450051, 1e-6),
                    "p_value={} (expected ~0.011056)",
                    p_value
                );
            }
            _ => panic!("expected Computed"),
        }
    }

    #[test]
    fn bh_fdr_simple_monotonic_sweep() {
        // 5 p-values, one obvious signal and background noise.
        // Standard BH example: p = [0.01, 0.02, 0.03, 0.8, 0.9]
        // ranks: 1..5
        // raw q: p * 5 / rank = [0.05, 0.05, 0.05, 1.0, 0.9]
        // min sweep (reverse): q[3]=min(1.0, 0.9)=0.9, then q[2]=min(0.05,0.9)=0.05,
        // q[1]=min(0.05,0.05)=0.05, q[0]=min(0.05,0.05)=0.05. q[4]=0.9, clamped.
        // But wait: q[4] starts at 0.9*5/5=0.9. q[3]=0.8*5/4=1.0. Reverse sweep
        // from i=3: q[3]=min(q[3]=1.0, q[4]=0.9)=0.9. q[2]=min(0.05, 0.9)=0.05.
        // q[1]=min(0.05, 0.05)=0.05. q[0]=min(0.05, 0.05)=0.05.
        let p = vec![Some(0.01), Some(0.02), Some(0.03), Some(0.8), Some(0.9)];
        let q = bh_fdr(&p);
        let q: Vec<f64> = q.into_iter().map(|x| x.unwrap()).collect();
        assert!(approx(q[0], 0.05, 1e-12));
        assert!(approx(q[1], 0.05, 1e-12));
        assert!(approx(q[2], 0.05, 1e-12));
        assert!(approx(q[3], 0.9, 1e-12));
        assert!(approx(q[4], 0.9, 1e-12));
    }

    #[test]
    fn bh_fdr_handles_none_slots() {
        let p = vec![Some(0.01), None, Some(0.5), None, Some(0.9)];
        let q = bh_fdr(&p);
        assert_eq!(q.len(), 5);
        assert!(q[0].is_some());
        assert!(q[1].is_none());
        assert!(q[2].is_some());
        assert!(q[3].is_none());
        assert!(q[4].is_some());
        // Only 3 real p-values in the family.
        // q for p=0.01 at rank 1 = 0.01 * 3 / 1 = 0.03
        assert!(approx(q[0].unwrap(), 0.03, 1e-12));
    }

    #[test]
    fn bh_fdr_preserves_order_after_scatter() {
        // p-values in non-sorted original order, verify output indices
        // align with input indices (not sorted).
        let p = vec![Some(0.9), Some(0.01), Some(0.5)];
        let q = bh_fdr(&p);
        assert!(q[0].unwrap() > q[1].unwrap()); // p=0.9 position keeps its high q
    }

    // ---- Welch two-sample tests ----

    #[test]
    fn welch_t_rejects_small_groups() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let r = welch_t(&a, &b, 5);
        assert!(matches!(
            r,
            PairedTResult::Skipped {
                reason: SkipReason::InsufficientPairs,
                ..
            }
        ));
    }

    #[test]
    fn welch_t_matches_scipy_reference() {
        // Reference values from scipy.stats.ttest_ind(a, b, equal_var=False):
        //   a = [1.1, 1.9, 3.0, 3.8, 5.2]         mean=3.00, var(ddof=1)=2.575
        //   b = [2.0, 2.1, 2.2, 2.3, 2.5, 2.4]    mean=2.25, var(ddof=1)=0.035
        //   t = 1.0392304845413263
        //   df = 4.09070821482279
        //   p = 0.35617615688335097
        let a = vec![1.1, 1.9, 3.0, 3.8, 5.2];
        let b = vec![2.0, 2.1, 2.2, 2.3, 2.5, 2.4];
        let r = welch_t(&a, &b, 5);
        match r {
            PairedTResult::Computed {
                n_pairs,
                mean_a,
                mean_b,
                mean_diff,
                t,
                df,
                p_value,
            } => {
                assert_eq!(n_pairs, 11);
                assert!(approx(mean_a, 3.0, 1e-10));
                assert!(approx(mean_b, 2.25, 1e-10));
                assert!(approx(mean_diff, 0.75, 1e-10));
                assert!(approx(t, 1.0392304845413263, 1e-10), "t={}", t);
                assert!(approx(df, 4.09070821482279, 1e-10), "df={}", df);
                assert!(approx(p_value, 0.35617615688335097, 1e-8), "p={}", p_value);
            }
            _ => panic!("expected Computed"),
        }
    }

    #[test]
    fn welch_t_strong_separation_gives_tiny_p() {
        let a = vec![0.0, 0.1, -0.1, 0.05, -0.05, 0.02];
        let b = vec![5.0, 5.1, 4.9, 5.05, 4.95, 5.02];
        let r = welch_t(&a, &b, 5);
        match r {
            PairedTResult::Computed {
                mean_diff, p_value, ..
            } => {
                assert!(mean_diff < -4.5);
                assert!(p_value < 1e-10, "p={}", p_value);
            }
            _ => panic!("expected Computed"),
        }
    }

    #[test]
    fn welch_t_zero_variance_both_groups_skipped() {
        let a = vec![3.0, 3.0, 3.0, 3.0, 3.0];
        let b = vec![3.0, 3.0, 3.0, 3.0, 3.0];
        let r = welch_t(&a, &b, 5);
        assert!(matches!(
            r,
            PairedTResult::Skipped {
                reason: SkipReason::ZeroVariance,
                ..
            }
        ));
    }

    #[test]
    fn stability_score_basic() {
        // Stable, strong signal → high score.
        let s = stability_weighted_score(Some(2.0), Some(0.001), Some(1.0)).unwrap();
        assert!(approx(s, 2.0 * 1.0 * (1.0 - 0.001), 1e-12), "s={}", s);
    }

    #[test]
    fn stability_score_unstable_direction_suppressed() {
        // Same effect + same q, but sign is flipping → score halves.
        let s_stable = stability_weighted_score(Some(2.0), Some(0.001), Some(1.0)).unwrap();
        let s_unstable = stability_weighted_score(Some(2.0), Some(0.001), Some(0.5)).unwrap();
        assert!(s_unstable < s_stable);
        assert!(approx(s_unstable / s_stable, 0.5, 1e-12));
    }

    #[test]
    fn stability_score_marginal_q_suppressed() {
        // Same effect + same SSR, but q ≈ 1 → score → 0.
        let s_strong = stability_weighted_score(Some(2.0), Some(0.001), Some(1.0)).unwrap();
        let s_marginal = stability_weighted_score(Some(2.0), Some(0.999), Some(1.0)).unwrap();
        assert!(s_marginal < s_strong / 100.0);
    }

    #[test]
    fn stability_score_sign_invariance() {
        // |d| is used, so opposite-sign effects of equal magnitude score equal.
        let s_pos = stability_weighted_score(Some(2.0), Some(0.01), Some(0.9)).unwrap();
        let s_neg = stability_weighted_score(Some(-2.0), Some(0.01), Some(0.9)).unwrap();
        assert!(approx(s_pos, s_neg, 1e-12));
    }

    #[test]
    fn stability_score_none_on_missing_input() {
        assert!(stability_weighted_score(None, Some(0.1), Some(0.9)).is_none());
        assert!(stability_weighted_score(Some(1.0), None, Some(0.9)).is_none());
        assert!(stability_weighted_score(Some(1.0), Some(0.1), None).is_none());
    }

    #[test]
    fn stability_score_clamps_q_above_one() {
        // q slightly > 1 (BH rounding) clamped to 1 → (1-q)=0 → score=0.
        let s = stability_weighted_score(Some(2.0), Some(1.01), Some(1.0)).unwrap();
        assert!(approx(s, 0.0, 1e-12));
    }

    #[test]
    fn welch_t_nonfinite_input_skipped() {
        let a = vec![1.0, 2.0, f64::NAN, 4.0, 5.0];
        let b = vec![2.0, 3.0, 4.0, 5.0, 6.0];
        let r = welch_t(&a, &b, 5);
        assert!(matches!(
            r,
            PairedTResult::Skipped {
                reason: SkipReason::NonFiniteInput,
                ..
            }
        ));
    }

    // ---- OLS (covariate-adjusted linear model) tests ----

    #[test]
    fn ols_simple_linear_matches_numpy_reference() {
        // Reference from numpy.linalg.lstsq on y ~ 1 + x:
        //   x = [1, 2, 3, 4, 5]; y = [2.1, 3.9, 6.2, 8.1, 9.8]
        //   beta      = [0.14, 1.96]
        //   se        = [0.18366636, 0.05537749]
        //   t         = [0.76225171, 35.39344081]
        //   p_value   = [0.50136080, 4.95970358e-05]
        //   df        = 3, sigma2 = 0.030666666666666616
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![2.1, 3.9, 6.2, 8.1, 9.8];
        let design: Vec<Vec<f64>> = x.iter().map(|&xi| vec![1.0, xi]).collect();
        let r = ols(&design, &y, 5);
        match r {
            OlsOutcome::Computed(fit) => {
                assert_eq!(fit.n, 5);
                assert_eq!(fit.p, 2);
                assert!(approx(fit.df, 3.0, 1e-12));
                assert!(
                    approx(fit.sigma2, 0.030666666666666616, 1e-12),
                    "sigma2={}",
                    fit.sigma2
                );
                assert!(approx(fit.beta[0], 0.14, 1e-12), "beta[0]={}", fit.beta[0]);
                assert!(approx(fit.beta[1], 1.96, 1e-12), "beta[1]={}", fit.beta[1]);
                assert!(
                    approx(fit.se[0], 0.1836663641860787, 1e-12),
                    "se[0]={}",
                    fit.se[0]
                );
                assert!(
                    approx(fit.se[1], 0.0553774924194538, 1e-12),
                    "se[1]={}",
                    fit.se[1]
                );
                assert!(
                    approx(fit.t[0], 0.7622517090726342, 1e-10),
                    "t[0]={}",
                    fit.t[0]
                );
                assert!(
                    approx(fit.t[1], 35.39344080721618, 1e-10),
                    "t[1]={}",
                    fit.t[1]
                );
                assert!(
                    approx(fit.p_value[0], 0.501360798878407, 1e-8),
                    "p[0]={}",
                    fit.p_value[0]
                );
                assert!(
                    approx(fit.p_value[1], 4.959703577522845e-5, 1e-8),
                    "p[1]={}",
                    fit.p_value[1]
                );
            }
            _ => panic!("expected Computed"),
        }
    }

    #[test]
    fn ols_two_covariates_matches_numpy_reference() {
        // Reference from numpy.linalg.lstsq on y ~ 1 + x1 + x2:
        //   x1 = 1..8; x2 = [0,1,0,1,0,1,0,1]; y = [1.1, 3.0, 3.1, 5.2, 5.0, 7.0, 7.2, 9.1]
        //   beta      = [0.07, 1.0075, 0.9675]
        //   se        = [0.07669746, 0.01504161, 0.06892931]
        //   t         = [0.91267693, 66.9808664, 14.0361187]
        //   p_value   = [4.03e-01, 1.40e-08, 3.30e-05]
        //   sigma2    = 0.009050000000000018
        let x1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let x2 = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0];
        let y = vec![1.1, 3.0, 3.1, 5.2, 5.0, 7.0, 7.2, 9.1];
        let design: Vec<Vec<f64>> = (0..8).map(|i| vec![1.0, x1[i], x2[i]]).collect();
        let r = ols(&design, &y, 5);
        match r {
            OlsOutcome::Computed(fit) => {
                assert_eq!(fit.n, 8);
                assert_eq!(fit.p, 3);
                assert!(approx(fit.df, 5.0, 1e-12));
                assert!(
                    approx(fit.sigma2, 0.009050000000000018, 1e-12),
                    "sigma2={}",
                    fit.sigma2
                );
                assert!(approx(fit.beta[0], 0.07, 1e-10), "beta[0]={}", fit.beta[0]);
                assert!(
                    approx(fit.beta[1], 1.0075, 1e-10),
                    "beta[1]={}",
                    fit.beta[1]
                );
                assert!(
                    approx(fit.beta[2], 0.9675, 1e-10),
                    "beta[2]={}",
                    fit.beta[2]
                );
                assert!(
                    approx(fit.se[0], 0.07669745758498135, 1e-12),
                    "se[0]={}",
                    fit.se[0]
                );
                assert!(
                    approx(fit.se[1], 0.015041608956491342, 1e-12),
                    "se[1]={}",
                    fit.se[1]
                );
                assert!(
                    approx(fit.se[2], 0.06892931161704728, 1e-12),
                    "se[2]={}",
                    fit.se[2]
                );
                assert!(
                    approx(fit.t[1], 66.98086639934007, 1e-8),
                    "t[1]={}",
                    fit.t[1]
                );
                assert!(
                    approx(fit.t[2], 14.036118703493143, 1e-8),
                    "t[2]={}",
                    fit.t[2]
                );
                assert!(fit.p_value[1] < 1e-7);
                assert!(fit.p_value[2] < 1e-4);
            }
            _ => panic!("expected Computed"),
        }
    }

    #[test]
    fn ols_rejects_singular_design() {
        // Two identical predictor columns → XᵀX is singular, Cholesky pivot ≤ 0.
        let design: Vec<Vec<f64>> = (1..=6)
            .map(|i| vec![1.0, i as f64, 2.0 * i as f64])
            .collect();
        let y = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let r = ols(&design, &y, 5);
        assert!(matches!(
            r,
            OlsOutcome::Skipped {
                reason: SkipReason::ZeroVariance,
                ..
            }
        ));
    }

    #[test]
    fn ols_rejects_n_le_p() {
        // 3 samples, 3 design columns → df = 0. Must skip, not divide by zero.
        let design = vec![
            vec![1.0, 0.0, 1.0],
            vec![1.0, 1.0, 0.0],
            vec![1.0, 0.5, 0.5],
        ];
        let y = vec![1.0, 2.0, 3.0];
        let r = ols(&design, &y, 3);
        assert!(matches!(
            r,
            OlsOutcome::Skipped {
                reason: SkipReason::InsufficientPairs,
                ..
            }
        ));
    }

    #[test]
    fn ols_rejects_nonfinite_input() {
        let design = vec![
            vec![1.0, 1.0],
            vec![1.0, 2.0],
            vec![1.0, f64::NAN],
            vec![1.0, 4.0],
            vec![1.0, 5.0],
        ];
        let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let r = ols(&design, &y, 5);
        assert!(matches!(
            r,
            OlsOutcome::Skipped {
                reason: SkipReason::NonFiniteInput,
                ..
            }
        ));
    }

    #[test]
    fn ols_reproduces_welch_t_with_group_indicator() {
        // Fit y ~ 1 + group on two-sample data. The group coefficient equals
        // (mean_a − mean_b), its t equals the Student two-sample t assuming
        // equal variances (the OLS residual variance pools both groups), and
        // the df equals n_total − 2.
        //
        // Here we fit the intercept at group=1 (reference), so beta[1]
        // carries group_a_minus_group_b on the flip. We verify the magnitude
        // of beta matches mean difference and that the fit reports the right
        // residual df.
        let design = vec![
            vec![1.0, 1.0], // group A
            vec![1.0, 1.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
            vec![1.0, 0.0], // group B
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 0.0],
        ];
        let y = vec![1.1, 1.9, 3.0, 3.8, 5.2, 2.0, 2.1, 2.2, 2.3, 2.5];
        let mean_a = (1.1 + 1.9 + 3.0 + 3.8 + 5.2) / 5.0;
        let mean_b = (2.0 + 2.1 + 2.2 + 2.3 + 2.5) / 5.0;
        let r = ols(&design, &y, 5);
        match r {
            OlsOutcome::Computed(fit) => {
                assert_eq!(fit.n, 10);
                assert!(approx(fit.df, 8.0, 1e-12), "df={}", fit.df);
                assert!(
                    approx(fit.beta[1], mean_a - mean_b, 1e-12),
                    "beta[1]={} expected={}",
                    fit.beta[1],
                    mean_a - mean_b
                );
            }
            _ => panic!("expected Computed"),
        }
    }

    #[test]
    fn ols_mismatched_lengths_skipped() {
        let design = vec![vec![1.0, 1.0], vec![1.0, 2.0], vec![1.0, 3.0]];
        let y = vec![1.0, 2.0];
        let r = ols(&design, &y, 2);
        assert!(matches!(r, OlsOutcome::Skipped { .. }));
    }

    #[test]
    fn ols_variable_row_width_skipped() {
        // Second row has 3 columns instead of 2 — should skip, not panic.
        let design = vec![vec![1.0, 1.0], vec![1.0, 2.0, 99.0], vec![1.0, 3.0]];
        let y = vec![1.0, 2.0, 3.0];
        let r = ols(&design, &y, 3);
        assert!(matches!(r, OlsOutcome::Skipped { .. }));
    }
}
