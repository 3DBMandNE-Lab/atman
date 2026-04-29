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

use statrs::distribution::{ContinuousCDF, FisherSnedecor, StudentsT};

/// Numerically stable two-sided Student's t p-value.
///
/// Use the distribution survival function directly rather than
/// `2 * (1 - cdf(|t|))`; the subtraction path rounds to zero for large but
/// still representable tails such as `t = 39, df = 326`.
pub fn two_sided_t_p_value(t: f64, df: f64) -> Option<f64> {
    if !t.is_finite() || !df.is_finite() || df <= 0.0 {
        return None;
    }
    let dist = StudentsT::new(0.0, 1.0, df).ok()?;
    Some((2.0 * dist.sf(t.abs())).clamp(0.0, 1.0))
}

/// Numerically stable right-tail F-test p-value.
pub fn right_tail_f_p_value(f: f64, df_num: f64, df_den: f64) -> Option<f64> {
    if !f.is_finite()
        || !df_num.is_finite()
        || !df_den.is_finite()
        || df_num <= 0.0
        || df_den <= 0.0
    {
        return None;
    }
    let dist = FisherSnedecor::new(df_num, df_den).ok()?;
    Some(dist.sf(f).clamp(0.0, 1.0))
}

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
    let p_value = (2.0 * t_dist.sf(t.abs())).clamp(0.0, 1.0);

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
    let p_value = (2.0 * t_dist.sf(t.abs())).clamp(0.0, 1.0);

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
    /// Residual variance estimate, `RSS / df` (σ²_res).
    pub sigma2: f64,
    /// For mixed-effects fits only: ratio of random-intercept variance
    /// to residual variance (`σ²_u / σ²_res`). `None` for plain OLS.
    /// Consumers that need the absolute random-intercept variance
    /// compute `variance_ratio * sigma2`.
    pub variance_ratio: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum OlsOutcome {
    /// Fit succeeded.
    Computed(OlsFit),
    /// Fit was not run. Reason captured for the report sidecar.
    Skipped { reason: SkipReason, n: usize },
}

/// Lower-triangular Cholesky factor `L` of a symmetric positive-definite
/// matrix `a`, such that `a = L · Lᵀ`. Returns `None` when `a` is not SPD
/// (a non-positive pivot is encountered).
///
/// `a` is expected to be `n × n` row-major. Only the lower triangle is
/// read; the upper triangle is ignored.
pub(crate) fn cholesky_lower(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    if n == 0 {
        return Some(Vec::new());
    }
    if a.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut l = vec![vec![0.0_f64; n]; n];
    for j in 0..n {
        let diag_sum: f64 = l[j].iter().take(j).map(|v| v * v).sum();
        let diag = a[j][j] - diag_sum;
        if !diag.is_finite() || diag <= 0.0 {
            return None;
        }
        l[j][j] = diag.sqrt();
        for i in (j + 1)..n {
            let off_sum: f64 = l[i]
                .iter()
                .zip(l[j].iter())
                .take(j)
                .map(|(li, lj)| li * lj)
                .sum();
            l[i][j] = (a[i][j] - off_sum) / l[j][j];
        }
    }
    Some(l)
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
    let Some(l) = cholesky_lower(&xtx) else {
        return OlsOutcome::Skipped {
            reason: SkipReason::ZeroVariance,
            n,
        };
    };

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
                (2.0 * t_dist.sf(ti.abs())).clamp(0.0, 1.0)
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
        variance_ratio: None,
    })
}

/// Random-intercept linear mixed model with one grouping factor.
///
/// This fits `y = Xβ + u_group + ε`, with `u ~ N(0, σ_u²)` and
/// `ε ~ N(0, σ_e²)`, by profiling the REML objective over
/// `λ = σ_u² / σ_e²`. For a fixed `λ`, the marginal covariance is
/// `V = σ_e² (I + λZZᵀ)`, and `β` is estimated by GLS. The reported
/// coefficient standard errors use the profiled REML scale estimate and
/// `(XᵀV⁻¹X)⁻¹`; p-values use a conservative residual df of `n - p`.
///
/// Initial scope is intentionally narrow: random intercept only, one grouping
/// factor, and no random slopes.
pub fn mixed_random_intercept(
    design: &[Vec<f64>],
    y: &[f64],
    groups: &[String],
    min_samples: usize,
) -> OlsOutcome {
    let n = design.len();
    if n != y.len() || n != groups.len() || n == 0 {
        return OlsOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }
    let p = design[0].len();
    if p == 0 || n < min_samples || n <= p {
        return OlsOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }
    for row in design {
        if row.len() != p || row.iter().any(|v| !v.is_finite()) {
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

    let group_index = group_indices(groups);
    if group_index.len() < 2 || group_index.iter().all(|idx| idx.len() < 2) {
        return OlsOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }

    let candidates: [f64; 8] = [-18.0, -12.0, -8.0, -4.0, 0.0, 4.0, 8.0, 12.0];
    let mut best_t = candidates[0];
    let mut best_obj = f64::INFINITY;
    for t in candidates {
        let obj = reml_objective(t.exp(), design, y, &group_index);
        if obj < best_obj {
            best_obj = obj;
            best_t = t;
        }
    }
    let mut left = best_t - 4.0;
    let mut right = best_t + 4.0;
    let gr = (5.0_f64.sqrt() - 1.0) / 2.0;
    let mut c = right - gr * (right - left);
    let mut d = left + gr * (right - left);
    let mut fc = reml_objective(c.exp(), design, y, &group_index);
    let mut fd = reml_objective(d.exp(), design, y, &group_index);
    for _ in 0..60 {
        if fc < fd {
            right = d;
            d = c;
            fd = fc;
            c = right - gr * (right - left);
            fc = reml_objective(c.exp(), design, y, &group_index);
        } else {
            left = c;
            c = d;
            fc = fd;
            d = left + gr * (right - left);
            fd = reml_objective(d.exp(), design, y, &group_index);
        }
    }
    let lambda = ((left + right) / 2.0).exp();
    let Some(gls) = gls_fit(lambda, design, y, &group_index) else {
        return OlsOutcome::Skipped {
            reason: SkipReason::ZeroVariance,
            n,
        };
    };
    let df = (n - p) as f64;
    if !(df.is_finite() && df > 0.0) {
        return OlsOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }
    let sigma2 = gls.rss / df;
    if !(sigma2.is_finite() && sigma2 >= 0.0) {
        return OlsOutcome::Skipped {
            reason: SkipReason::NonFiniteInput,
            n,
        };
    }
    let se: Vec<f64> = gls
        .xtvix_inv_diag
        .iter()
        .map(|v| (sigma2 * v).sqrt())
        .collect();
    let t: Vec<f64> = gls
        .beta
        .iter()
        .zip(se.iter())
        .map(|(b, s)| if *s == 0.0 { f64::NAN } else { b / s })
        .collect();
    let t_dist = match StudentsT::new(0.0, 1.0, df) {
        Ok(d) => d,
        Err(_) => {
            return OlsOutcome::Skipped {
                reason: SkipReason::InsufficientPairs,
                n,
            };
        }
    };
    let p_value = t
        .iter()
        .map(|&ti| {
            if ti.is_finite() {
                (2.0 * t_dist.sf(ti.abs())).clamp(0.0, 1.0)
            } else {
                f64::NAN
            }
        })
        .collect();

    OlsOutcome::Computed(OlsFit {
        n,
        p,
        beta: gls.beta,
        se,
        t,
        p_value,
        df,
        sigma2,
        variance_ratio: Some(lambda),
    })
}

struct GlsFit {
    beta: Vec<f64>,
    xtvix_inv_diag: Vec<f64>,
    rss: f64,
    logdet_xtvix: f64,
}

fn group_indices(groups: &[String]) -> Vec<Vec<usize>> {
    let mut map: std::collections::BTreeMap<&str, Vec<usize>> = std::collections::BTreeMap::new();
    for (idx, group) in groups.iter().enumerate() {
        map.entry(group.as_str()).or_default().push(idx);
    }
    map.into_values().collect()
}

fn reml_objective(lambda: f64, design: &[Vec<f64>], y: &[f64], group_index: &[Vec<usize>]) -> f64 {
    let Some(fit) = gls_fit(lambda, design, y, group_index) else {
        return f64::INFINITY;
    };
    let n = design.len();
    let p = design[0].len();
    if n <= p || fit.rss <= 0.0 || !fit.rss.is_finite() {
        return f64::INFINITY;
    }
    let logdet_v0: f64 = group_index
        .iter()
        .map(|idx| (1.0 + lambda * idx.len() as f64).ln())
        .sum();
    let df = (n - p) as f64;
    df * (fit.rss / df).ln() + logdet_v0 + fit.logdet_xtvix
}

fn gls_fit(
    lambda: f64,
    design: &[Vec<f64>],
    y: &[f64],
    group_index: &[Vec<usize>],
) -> Option<GlsFit> {
    let p = design[0].len();
    let mut wx = vec![vec![0.0; p]; design.len()];
    let mut wy = vec![0.0; y.len()];
    apply_compound_symmetry_inverse(lambda, design, y, group_index, &mut wx, &mut wy);

    let mut xtvix = vec![vec![0.0; p]; p];
    let mut xtviy = vec![0.0; p];
    for i in 0..design.len() {
        for j in 0..p {
            xtviy[j] += design[i][j] * wy[i];
            for k in 0..p {
                xtvix[j][k] += design[i][j] * wx[i][k];
            }
        }
    }
    let chol = cholesky(&xtvix)?;
    let beta = solve_cholesky(&chol, &xtviy);
    let inv_diag = inverse_diag_from_cholesky(&chol);
    let logdet_xtvix = 2.0
        * chol
            .iter()
            .enumerate()
            .map(|(i, row)| row[i].ln())
            .sum::<f64>();

    let residuals: Vec<f64> = design
        .iter()
        .zip(y.iter())
        .map(|(row, yi)| {
            let pred: f64 = row.iter().zip(beta.iter()).map(|(x, b)| x * b).sum();
            yi - pred
        })
        .collect();
    let mut wres = vec![0.0; residuals.len()];
    apply_compound_symmetry_inverse_to_vector(lambda, &residuals, group_index, &mut wres);
    let rss = residuals
        .iter()
        .zip(wres.iter())
        .map(|(r, wr)| r * wr)
        .sum();
    Some(GlsFit {
        beta,
        xtvix_inv_diag: inv_diag,
        rss,
        logdet_xtvix,
    })
}

fn apply_compound_symmetry_inverse(
    lambda: f64,
    design: &[Vec<f64>],
    y: &[f64],
    group_index: &[Vec<usize>],
    wx: &mut [Vec<f64>],
    wy: &mut [f64],
) {
    for idx in group_index {
        let m = idx.len() as f64;
        let factor = lambda / (1.0 + lambda * m);
        let mut sum_y = 0.0;
        let mut sum_x = vec![0.0; design[0].len()];
        for &i in idx {
            sum_y += y[i];
            for (col, value) in design[i].iter().enumerate() {
                sum_x[col] += value;
            }
        }
        for &i in idx {
            wy[i] = y[i] - factor * sum_y;
            for (col, value) in design[i].iter().enumerate() {
                wx[i][col] = value - factor * sum_x[col];
            }
        }
    }
}

fn apply_compound_symmetry_inverse_to_vector(
    lambda: f64,
    values: &[f64],
    group_index: &[Vec<usize>],
    out: &mut [f64],
) {
    for idx in group_index {
        let m = idx.len() as f64;
        let factor = lambda / (1.0 + lambda * m);
        let sum: f64 = idx.iter().map(|&i| values[i]).sum();
        for &i in idx {
            out[i] = values[i] - factor * sum;
        }
    }
}

fn cholesky(matrix: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let p = matrix.len();
    let mut l = vec![vec![0.0; p]; p];
    for j in 0..p {
        let diag_sum: f64 = l[j].iter().take(j).map(|v| v * v).sum();
        let diag = matrix[j][j] - diag_sum;
        if !diag.is_finite() || diag <= 0.0 {
            return None;
        }
        l[j][j] = diag.sqrt();
        for i in (j + 1)..p {
            let off_sum: f64 = l[i]
                .iter()
                .zip(l[j].iter())
                .take(j)
                .map(|(li, lj)| li * lj)
                .sum();
            l[i][j] = (matrix[i][j] - off_sum) / l[j][j];
        }
    }
    Some(l)
}

fn solve_cholesky(l: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let p = rhs.len();
    let mut z = vec![0.0; p];
    for i in 0..p {
        let mut sum = rhs[i];
        for (k, zk) in z.iter().enumerate().take(i) {
            sum -= l[i][k] * zk;
        }
        z[i] = sum / l[i][i];
    }
    let mut out = vec![0.0; p];
    for i in (0..p).rev() {
        let mut sum = z[i];
        for k in (i + 1)..p {
            sum -= l[k][i] * out[k];
        }
        out[i] = sum / l[i][i];
    }
    out
}

fn inverse_diag_from_cholesky(l: &[Vec<f64>]) -> Vec<f64> {
    let p = l.len();
    let mut out = vec![0.0; p];
    for target in 0..p {
        let mut rhs = vec![0.0; p];
        rhs[target] = 1.0;
        let col = solve_cholesky(l, &rhs);
        out[target] = col[target];
    }
    out
}

/// One fitted linear contrast on an OLS design.
#[derive(Debug, Clone, PartialEq)]
pub struct ContrastResult {
    /// `cᵀβ`.
    pub estimate: f64,
    /// `sqrt(cᵀ (XᵀX)⁻¹ c · σ²)`.
    pub se: f64,
    /// `estimate / se`.
    pub t: f64,
    /// Two-sided Student-t p on `df` degrees of freedom.
    pub p_value: f64,
}

/// One-contrast inference on the already-fit OLS coefficients.
///
/// `design` is the n × p sample-major design matrix used for the
/// original fit. `beta` is the fitted coefficients (length p).
/// `contrast_weights` is a length-p weight vector `c`; the test
/// evaluates `H₀: cᵀβ = 0`. `sigma2` is the residual variance
/// estimate from the fit (`RSS / df`). `df` is the residual
/// degrees of freedom (`n − p`).
///
/// Numerically, we rebuild `XᵀX`, Cholesky-factor it, solve
/// `XᵀX z = c`, and read off `cᵀ(XᵀX)⁻¹c = cᵀz`. Returns `None`
/// when the design is singular, sizes are inconsistent, or
/// `df <= 0`.
pub fn contrast_inference(
    design: &[Vec<f64>],
    beta: &[f64],
    contrast_weights: &[f64],
    sigma2: f64,
    df: f64,
) -> Option<ContrastResult> {
    if design.is_empty() || beta.is_empty() || contrast_weights.is_empty() {
        return None;
    }
    let p = beta.len();
    if design[0].len() != p || contrast_weights.len() != p {
        return None;
    }
    if !df.is_finite() || df <= 0.0 || !sigma2.is_finite() || sigma2 < 0.0 {
        return None;
    }
    // Rebuild XᵀX (p × p).
    let mut xtx = vec![vec![0.0_f64; p]; p];
    for row in design {
        if row.len() != p || row.iter().any(|v| !v.is_finite()) {
            return None;
        }
        for i in 0..p {
            for j in 0..p {
                xtx[i][j] += row[i] * row[j];
            }
        }
    }
    let l = cholesky_lower(&xtx)?;
    // Solve XᵀX z = c via forward/back substitution on L.
    let mut z_forward = vec![0.0_f64; p];
    for i in 0..p {
        let mut sum = contrast_weights[i];
        for k in 0..i {
            sum -= l[i][k] * z_forward[k];
        }
        z_forward[i] = sum / l[i][i];
    }
    let mut z = vec![0.0_f64; p];
    for i in (0..p).rev() {
        let mut sum = z_forward[i];
        for k in (i + 1)..p {
            sum -= l[k][i] * z[k];
        }
        z[i] = sum / l[i][i];
    }
    // Variance scale `cᵀ (XᵀX)⁻¹ c = cᵀz`.
    let mut var_scale = 0.0_f64;
    for i in 0..p {
        var_scale += contrast_weights[i] * z[i];
    }
    if !var_scale.is_finite() || var_scale < 0.0 {
        return None;
    }
    let se = (var_scale * sigma2).sqrt();
    let mut estimate = 0.0_f64;
    for i in 0..p {
        estimate += contrast_weights[i] * beta[i];
    }
    let t = if se == 0.0 { f64::NAN } else { estimate / se };
    let p_value = match StudentsT::new(0.0, 1.0, df) {
        Ok(dist) if t.is_finite() => (2.0 * dist.sf(t.abs())).clamp(0.0, 1.0),
        _ => f64::NAN,
    };
    Some(ContrastResult {
        estimate,
        se,
        t,
        p_value,
    })
}

/// Omnibus F-test of the joint hypothesis `β_j = 0` for every
/// `j` in `factor_columns`. Uses the standard
/// `F = (β_S' · [A⁻¹]_SS · β_S) / (k · σ²)` form where `A = XᵀX`,
/// `S` is the column subset, and `k = |S|`. `df_num = k`,
/// `df_den = df`.
///
/// Returns `None` when the design is singular or dimensions are
/// inconsistent.
#[derive(Debug, Clone, PartialEq)]
pub struct OmnibusF {
    pub f_statistic: f64,
    pub df_num: usize,
    pub df_den: f64,
    pub p_value: f64,
}

pub fn omnibus_f_test(
    design: &[Vec<f64>],
    beta: &[f64],
    factor_columns: &[usize],
    sigma2: f64,
    df: f64,
) -> Option<OmnibusF> {
    if design.is_empty() || beta.is_empty() || factor_columns.is_empty() {
        return None;
    }
    let p = beta.len();
    if design[0].len() != p || !sigma2.is_finite() || sigma2 <= 0.0 || !(df.is_finite() && df > 0.0)
    {
        return None;
    }
    for &c in factor_columns {
        if c >= p {
            return None;
        }
    }

    // XᵀX and its full inverse (only the factor-column submatrix is
    // needed, but p is small in practice).
    let mut xtx = vec![vec![0.0_f64; p]; p];
    for row in design {
        if row.len() != p || row.iter().any(|v| !v.is_finite()) {
            return None;
        }
        for i in 0..p {
            for j in 0..p {
                xtx[i][j] += row[i] * row[j];
            }
        }
    }
    let l = cholesky_lower(&xtx)?;
    // Build (XᵀX)⁻¹ column-by-column from L.
    let mut inv = vec![vec![0.0_f64; p]; p];
    (0..p).for_each(|target| {
        let mut z = vec![0.0_f64; p];
        for i in 0..p {
            let mut sum = if i == target { 1.0 } else { 0.0 };
            for k in 0..i {
                sum -= l[i][k] * z[k];
            }
            z[i] = sum / l[i][i];
        }
        let mut col = vec![0.0_f64; p];
        for i in (0..p).rev() {
            let mut sum = z[i];
            for k in (i + 1)..p {
                sum -= l[k][i] * col[k];
            }
            col[i] = sum / l[i][i];
        }
        for (i, row) in inv.iter_mut().enumerate() {
            row[target] = col[i];
        }
    });
    // Extract principal submatrix at factor_columns.
    let k = factor_columns.len();
    let mut a_ss = vec![vec![0.0_f64; k]; k];
    for (a, &ia) in factor_columns.iter().enumerate() {
        for (b, &ib) in factor_columns.iter().enumerate() {
            a_ss[a][b] = inv[ia][ib];
        }
    }
    // Invert the k × k submatrix via its own Cholesky.
    let ls = cholesky_lower(&a_ss)?;
    // Build β_S and solve A_ss · x = β_S.
    let beta_s: Vec<f64> = factor_columns.iter().map(|&i| beta[i]).collect();
    let mut z_fw = vec![0.0_f64; k];
    for i in 0..k {
        let mut sum = beta_s[i];
        for kk in 0..i {
            sum -= ls[i][kk] * z_fw[kk];
        }
        z_fw[i] = sum / ls[i][i];
    }
    let mut x = vec![0.0_f64; k];
    for i in (0..k).rev() {
        let mut sum = z_fw[i];
        for kk in (i + 1)..k {
            sum -= ls[kk][i] * x[kk];
        }
        x[i] = sum / ls[i][i];
    }
    let quad: f64 = beta_s.iter().zip(x.iter()).map(|(b, xx)| b * xx).sum();
    let f = quad / (k as f64 * sigma2);
    if !f.is_finite() || f < 0.0 {
        return None;
    }
    let p_value = match FisherSnedecor::new(k as f64, df) {
        Ok(dist) => dist.sf(f).clamp(0.0, 1.0),
        Err(_) => f64::NAN,
    };
    Some(OmnibusF {
        f_statistic: f,
        df_num: k,
        df_den: df,
        p_value,
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
    fn contrast_matches_single_coefficient_se_on_standard_design() {
        // A 3-column intercept + two indicator design with 9 samples.
        // Contrast c = (0, 1, 0) picks out β_1 directly. The resulting
        // se/t/p should match what ols() already returns for β_1.
        let design: Vec<Vec<f64>> = (0..9)
            .map(|i| {
                let cond = i / 3; // 0, 1, 2 each ×3
                vec![
                    1.0,
                    if cond == 1 { 1.0 } else { 0.0 },
                    if cond == 2 { 1.0 } else { 0.0 },
                ]
            })
            .collect();
        let y: Vec<f64> = (0..9)
            .map(|i| {
                let cond = i / 3;
                cond as f64 * 2.0 + (i as f64 * 0.01).sin()
            })
            .collect();
        let fit = match ols(&design, &y, 2) {
            OlsOutcome::Computed(f) => f,
            _ => panic!("ols should fit"),
        };
        let contrast = vec![0.0, 1.0, 0.0];
        let r = contrast_inference(&design, &fit.beta, &contrast, fit.sigma2, fit.df).unwrap();
        // Contrast estimate should match β_1 and SE should match se[1].
        assert!(approx(r.estimate, fit.beta[1], 1e-10));
        assert!(approx(r.se, fit.se[1], 1e-10));
        assert!(approx(r.t, fit.t[1], 1e-10));
        assert!(approx(r.p_value, fit.p_value[1], 1e-10));
    }

    #[test]
    fn contrast_difference_of_two_levels_matches_hand_computation() {
        // Same 3-column design. Contrast c = (0, 1, -1) is the
        // difference between level 1 and level 2. Recovery: the
        // estimate equals β_1 − β_2.
        let design: Vec<Vec<f64>> = (0..12)
            .map(|i| {
                let cond = i / 4;
                vec![
                    1.0,
                    if cond == 1 { 1.0 } else { 0.0 },
                    if cond == 2 { 1.0 } else { 0.0 },
                ]
            })
            .collect();
        let y: Vec<f64> = (0..12)
            .map(|i| {
                let cond = i / 4;
                cond as f64 * 1.5 + (i as f64 * 0.1).cos() * 0.05
            })
            .collect();
        let fit = match ols(&design, &y, 2) {
            OlsOutcome::Computed(f) => f,
            _ => panic!("ols should fit"),
        };
        let contrast = vec![0.0, 1.0, -1.0];
        let r = contrast_inference(&design, &fit.beta, &contrast, fit.sigma2, fit.df).unwrap();
        let expected_est = fit.beta[1] - fit.beta[2];
        assert!(
            approx(r.estimate, expected_est, 1e-10),
            "estimate {} vs expected {}",
            r.estimate,
            expected_est
        );
        // The SE of (β_1 − β_2) must be positive and finite.
        assert!(r.se.is_finite() && r.se > 0.0);
    }

    #[test]
    fn omnibus_f_under_h0_has_expected_distribution_p_range() {
        // y is random noise with no dependency on the factor; with
        // n=30 and 3 levels, p should be large (bounded well away
        // from 0). Using a deterministic pseudorandom sequence.
        let design: Vec<Vec<f64>> = (0..30)
            .map(|i| {
                let cond = i / 10;
                vec![
                    1.0,
                    if cond == 1 { 1.0 } else { 0.0 },
                    if cond == 2 { 1.0 } else { 0.0 },
                ]
            })
            .collect();
        let y: Vec<f64> = (0..30)
            .map(|i| {
                // Fresh deterministic noise — no cond dependence.
                let u = ((i as f64 * 0.97).sin() * 7.0).cos();
                u + (i as f64).ln_1p() * 0.15
            })
            .collect();
        let fit = match ols(&design, &y, 2) {
            OlsOutcome::Computed(f) => f,
            _ => panic!(),
        };
        let omni = omnibus_f_test(&design, &fit.beta, &[1, 2], fit.sigma2, fit.df).unwrap();
        assert_eq!(omni.df_num, 2);
        assert_eq!(omni.df_den, 27.0); // n - p = 30 - 3
        assert!(
            omni.p_value > 0.05,
            "expected null F p>0.05, got {}",
            omni.p_value
        );
    }

    #[test]
    fn omnibus_f_large_on_planted_factor_effect() {
        // Strong planted effect on the factor columns.
        let design: Vec<Vec<f64>> = (0..12)
            .map(|i| {
                let cond = i / 4;
                vec![
                    1.0,
                    if cond == 1 { 1.0 } else { 0.0 },
                    if cond == 2 { 1.0 } else { 0.0 },
                ]
            })
            .collect();
        let y: Vec<f64> = (0..12)
            .map(|i| {
                let cond = i / 4;
                cond as f64 * 3.0 + (i as f64 * 0.1).cos() * 0.05
            })
            .collect();
        let fit = match ols(&design, &y, 2) {
            OlsOutcome::Computed(f) => f,
            _ => panic!(),
        };
        let omni = omnibus_f_test(&design, &fit.beta, &[1, 2], fit.sigma2, fit.df).unwrap();
        assert!(
            omni.f_statistic > 10.0,
            "expected F>10, got {}",
            omni.f_statistic
        );
        assert!(
            omni.p_value < 0.001,
            "expected p<0.001, got {}",
            omni.p_value
        );
    }

    #[test]
    fn two_sided_t_p_value_uses_stable_tail_for_extreme_t() {
        let p = two_sided_t_p_value(39.4186141848897, 326.12788865688213).expect("finite p-value");
        assert!(p.is_finite() && p > 0.0, "p={p}");
        assert!(
            (p.log10() + 125.366440245).abs() < 1e-6,
            "log10(p)={}",
            p.log10()
        );
    }

    #[test]
    fn right_tail_f_p_value_uses_stable_tail_for_extreme_f() {
        let p = right_tail_f_p_value(500.0, 2.0, 300.0).expect("finite p-value");
        assert!(p.is_finite() && p > 0.0, "p={p}");
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

    #[test]
    fn cholesky_lower_factors_known_spd_matrix() {
        // A = L Lᵀ with L = [[2, 0], [6, 1]] ⇒ A = [[4, 12], [12, 37]]
        let a = vec![vec![4.0, 12.0], vec![12.0, 37.0]];
        let l = super::cholesky_lower(&a).expect("SPD");
        assert!((l[0][0] - 2.0).abs() < 1e-12);
        assert!((l[1][0] - 6.0).abs() < 1e-12);
        assert!((l[1][1] - 1.0).abs() < 1e-12);
        assert!(l[0][1].abs() < 1e-12); // upper triangle must be zero
    }

    #[test]
    fn cholesky_lower_returns_none_for_non_spd() {
        let a = vec![vec![1.0, 2.0], vec![2.0, 1.0]]; // indefinite
        assert!(super::cholesky_lower(&a).is_none());
    }

    #[test]
    fn cholesky_lower_handles_empty_matrix() {
        let a: Vec<Vec<f64>> = Vec::new();
        let l = super::cholesky_lower(&a).expect("empty is OK");
        assert!(l.is_empty());
    }
}
