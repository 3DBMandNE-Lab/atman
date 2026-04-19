//! limma 3.x–grade eBayes primitives for `atman de --test limma`.
//!
//! Pure-Rust port of the minimal subset of limma needed to reproduce
//! `eBayes(lmFit(y, X), trend, robust)` + TREAT + moderated F. Entry
//! point: `limma_fit`. Everything else is internal algorithm plumbing
//! exposed for unit testing.
//!
//! Reference: Smyth 2004 (eBayes), Phipson et al. 2013 (TREAT),
//! Phipson et al. 2016 (robust), Ritchie et al. 2015 (limma v3).

use crate::de::cholesky_lower;
use statrs::distribution::{ContinuousCDF, FisherSnedecor, StudentsT};

/// Digamma ψ(x) = d/dx ln Γ(x). Asymptotic expansion for x ≥ 10,
/// recurrence `ψ(x) = ψ(x+1) - 1/x` for smaller x. The x ≥ 10 cutoff
/// keeps the truncated series (through 1/x⁸) under ~1e-10 absolute
/// error across tested values.
pub fn digamma(x: f64) -> f64 {
    if !x.is_finite() || x <= 0.0 {
        return f64::NAN;
    }
    if x < 10.0 {
        return digamma(x + 1.0) - 1.0 / x;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    // ψ(x) = ln x - 1/(2x) - 1/(12x²) + 1/(120x⁴) - 1/(252x⁶) + 1/(240x⁸) - …
    x.ln() - 0.5 * inv
        - inv2 * (1.0 / 12.0 - inv2 * (1.0 / 120.0 - inv2 * (1.0 / 252.0 - inv2 / 240.0)))
}

/// Trigamma ψ'(x). Asymptotic for x ≥ 10, recurrence `ψ'(x) = ψ'(x+1) + 1/x²`.
pub fn trigamma(x: f64) -> f64 {
    if !x.is_finite() || x <= 0.0 {
        return f64::NAN;
    }
    if x < 10.0 {
        return trigamma(x + 1.0) + 1.0 / (x * x);
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    let inv3 = inv2 * inv;
    // ψ'(x) = 1/x + 1/(2x²) + 1/(6x³) - 1/(30x⁵) + 1/(42x⁷) - 1/(30x⁹) + …
    inv + 0.5 * inv2
        + inv3 * (1.0 / 6.0 - inv2 * (1.0 / 30.0 - inv2 * (1.0 / 42.0 - inv2 / 30.0)))
}

/// Tetragamma ψ''(x) = d/dx ψ'(x). Used by `trigamma_inverse` Newton step.
pub fn tetragamma(x: f64) -> f64 {
    if !x.is_finite() || x <= 0.0 {
        return f64::NAN;
    }
    if x < 10.0 {
        return tetragamma(x + 1.0) - 2.0 / (x * x * x);
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    let inv3 = inv2 * inv;
    let inv4 = inv2 * inv2;
    // ψ''(x) = -1/x² - 1/x³ - 1/(2x⁴) + 1/(6x⁶) - 1/(6x⁸) + …
    -inv2 - inv3 - 0.5 * inv4 + inv4 * inv2 * (1.0 / 6.0 - inv2 / 6.0)
}

/// Inverse trigamma: given `y > 0`, find `x > 0` such that `trigamma(x) = y`.
/// Newton iteration on `ψ'(x) - y = 0` using `ψ''(x)` as the derivative.
/// Initial guess from the asymptotic `ψ'(x) ≈ 1/x + 1/(2x²)`.
pub fn trigamma_inverse(y: f64) -> f64 {
    if !y.is_finite() || y <= 0.0 {
        return f64::NAN;
    }
    // Asymptotic inverse: for small y (large x), x ≈ 1/y + 1/2; for large y
    // (small x), x ≈ 1/sqrt(y). Blend: use 0.5 + 1/y.
    let mut x = if y > 1e7 {
        1.0 / y.sqrt()
    } else if y < 1e-6 {
        1.0 / y
    } else {
        0.5 + 1.0 / y
    };
    for _ in 0..50 {
        let tri = trigamma(x);
        let tet = tetragamma(x);
        if tet == 0.0 || !tet.is_finite() {
            break;
        }
        let delta = (tri - y) / tet;
        let new_x = x - delta;
        if new_x <= 0.0 {
            x *= 0.5;
            continue;
        }
        x = new_x;
        if delta.abs() < 1e-12 * x.max(1.0) {
            break;
        }
    }
    x
}

/// Fit a scaled-F distribution to a vector of per-feature sample
/// variances `s²` with `df_res` residual degrees of freedom, by the
/// method of moments on `z = log(s²)`.
///
/// Returns `(df_prior, s²_prior)` such that `s² ~ s²_prior · F(df_res, df_prior)`.
///
/// Algorithm (Smyth 2004 §2, identical to limma's `fitFDist`):
///
/// ```text
/// z_i  = log(s²_i)
/// m    = mean(z)
/// v    = Var(z)   // MLE estimator (matches limma's `fitFDist`): Σ(z - m)² / n
/// Solve
///   Var[z] = trigamma(df_res/2) + trigamma(df_prior/2)
///   E  [z] = digamma (df_res/2) - digamma (df_prior/2) + log(df_prior · s²_prior / df_res)
/// for (df_prior, s²_prior).
/// ```
///
/// Variance of `z` uses the MLE estimator (divide by `n`, not `n-1`) to match
/// limma 3.x `fitFDist`, which computes `var(e) * (length(e)-1) / length(e)`.
///
/// Returns `None` if the sample contains fewer than 2 positive finite
/// variances, or if the variance-of-z matches `trigamma(df_res/2)` (no
/// residual heterogeneity → `df_prior → ∞`; caller falls back to
/// no-prior shrinkage).
pub fn fit_f_dist(s2: &[f64], df_res: f64) -> Option<(f64, f64)> {
    let log_s2: Vec<f64> = s2
        .iter()
        .copied()
        .filter(|v| v.is_finite() && *v > 0.0)
        .map(f64::ln)
        .collect();
    if log_s2.len() < 2 {
        return None;
    }
    let n = log_s2.len() as f64;
    let mean = log_s2.iter().sum::<f64>() / n;
    // MLE estimator (matches limma's `fitFDist`): divide by n, not (n - 1).
    let var = log_s2
        .iter()
        .map(|z| {
            let d = z - mean;
            d * d
        })
        .sum::<f64>()
        / n;

    let t_res = trigamma(df_res / 2.0);
    let excess = var - t_res;
    if excess <= 0.0 || !excess.is_finite() {
        // No evidence for residual prior variance; let caller degrade
        // gracefully (df_prior = +inf, s²_prior = exp(mean)).
        return None;
    }
    let df_prior = 2.0 * trigamma_inverse(excess);
    if !df_prior.is_finite() || df_prior <= 0.0 {
        return None;
    }
    // mean = digamma(df_res/2) - digamma(df_prior/2) + log(df_prior · s²_prior / df_res)
    // ⇒ log(s²_prior) = mean - digamma(df_res/2) + digamma(df_prior/2) + log(df_res / df_prior)
    let log_s2_prior = mean - digamma(df_res / 2.0) + digamma(df_prior / 2.0)
        + (df_res / df_prior).ln();
    let s2_prior = log_s2_prior.exp();
    if !s2_prior.is_finite() || s2_prior <= 0.0 {
        return None;
    }
    Some((df_prior, s2_prior))
}

/// Output of `squeeze_var`: per-feature posterior variance and the
/// shared `(df_prior, s2_prior)`. Per-feature `df_total = df_prior + df_res`
/// is the reader's job (it never varies across features for a given fit).
pub struct SqueezeVarOutput {
    pub s2_posterior: Vec<f64>,
    pub df_prior: f64,
    pub s2_prior: f64,
}

/// Empirical-Bayes variance shrinkage.
///
/// - If `prior` is `Some((df_prior, s2_prior))` — typically from a prior
///   call to [`fit_f_dist`] — returns per-feature
///   `s²_post = (df_prior · s²_prior + df_res · s²) / (df_prior + df_res)`.
/// - If `prior` is `None`, returns `s²_post = s²` unchanged, `df_prior = ∞`,
///   `s²_prior = NaN`. This is the fallback for `fit_f_dist` failure.
pub fn squeeze_var(
    s2: &[f64],
    df_res: f64,
    prior: Option<(f64, f64)>,
) -> SqueezeVarOutput {
    match prior {
        Some((df_prior, s2_prior)) => {
            let denom = df_prior + df_res;
            let s2_posterior: Vec<f64> = s2
                .iter()
                .map(|v| (df_prior * s2_prior + df_res * v) / denom)
                .collect();
            SqueezeVarOutput {
                s2_posterior,
                df_prior,
                s2_prior,
            }
        }
        None => SqueezeVarOutput {
            s2_posterior: s2.to_vec(),
            df_prior: f64::INFINITY,
            s2_prior: f64::NAN,
        },
    }
}

/// Output of [`fit_f_dist_robust`]: the Winsorized-fit prior degrees of
/// freedom and prior variance. Same shape as the `(df_prior, s²_prior)`
/// tuple returned by [`fit_f_dist`], wrapped in a named struct so call
/// sites can record the robust-vs-non-robust provenance.
pub struct FitFDistOutput {
    pub df_prior: f64,
    pub s2_prior: f64,
}

/// Robust variant of [`fit_f_dist`] per Phipson et al. 2016: Winsorize
/// `z = log(s²)` before refitting, using limma's default
/// `winsor.tail.p = c(0.05, 0.1)` — clamp the bottom 5% of values to the
/// 5th-percentile cutoff and the top 10% to the 90th-percentile cutoff.
/// This down-weights per-feature outlier variances that would otherwise
/// pull `s²_prior` upward.
///
/// Default Winsor tails match limma's `winsor.tail.p = c(0.05, 0.1)`.
///
/// Variance of the Winsorized `z` uses the MLE estimator (divide by `n`,
/// not `n-1`) to match limma's `fitFDist` convention.
///
/// Implementation mirrors `fitFDistRobustly` in limma, simplified to
/// match the MVP's "one robust pass with fixed tails" scope.
pub fn fit_f_dist_robust(s2: &[f64], df_res: f64) -> Option<FitFDistOutput> {
    let mut log_s2: Vec<f64> = s2
        .iter()
        .copied()
        .filter(|v| v.is_finite() && *v > 0.0)
        .map(f64::ln)
        .collect();
    if log_s2.len() < 3 {
        return None;
    }
    log_s2.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let n = log_s2.len();
    let lower_tail = 0.05_f64;
    let upper_tail = 0.10_f64;
    let lower_idx = ((lower_tail * n as f64).floor() as usize).min(n - 1);
    let upper_idx = (n.saturating_sub(1))
        .saturating_sub(((upper_tail * n as f64).floor()) as usize);
    if upper_idx <= lower_idx {
        return None;
    }
    let lower_val = log_s2[lower_idx];
    let upper_val = log_s2[upper_idx];
    // Winsorize in place.
    for v in log_s2.iter_mut() {
        if *v < lower_val {
            *v = lower_val;
        } else if *v > upper_val {
            *v = upper_val;
        }
    }

    let n_f = n as f64;
    let mean = log_s2.iter().sum::<f64>() / n_f;
    // MLE estimator (matches limma's `fitFDist`): divide by n, not (n - 1).
    let var = log_s2
        .iter()
        .map(|z| {
            let d = z - mean;
            d * d
        })
        .sum::<f64>()
        / n_f;

    let t_res = trigamma(df_res / 2.0);
    let excess = var - t_res;
    if excess <= 0.0 {
        return None;
    }
    let df_prior = 2.0 * trigamma_inverse(excess);
    if !df_prior.is_finite() || df_prior <= 0.0 {
        return None;
    }
    let log_s2_prior = mean - digamma(df_res / 2.0) + digamma(df_prior / 2.0)
        + (df_res / df_prior).ln();
    let s2_prior = log_s2_prior.exp();
    if !s2_prior.is_finite() || s2_prior <= 0.0 {
        return None;
    }
    Some(FitFDistOutput { df_prior, s2_prior })
}

/// Fit a parametric mean–variance trend `log(s²) = a + b · log(μ) + c · log(μ)²`
/// by ordinary least squares across features. Returns the per-feature
/// fitted `s²_trend = exp(a + b·log μ + c·log² μ)`.
///
/// `means` and `s2` must have equal length and contain only finite
/// positive values (features with non-positive variance or mean are
/// dropped from the fit; their trend value is filled with the sample
/// mean of valid `s2`).
///
/// Returns `None` when fewer than 3 usable features remain (the
/// quadratic has 3 coefficients; caller falls back to no trend).
pub fn fit_parametric_trend(means: &[f64], s2: &[f64]) -> Option<Vec<f64>> {
    if means.len() != s2.len() {
        return None;
    }
    let mut xtx = vec![vec![0.0_f64; 3]; 3];
    let mut xty = [0.0_f64; 3];
    let mut n_used = 0_usize;
    for (m, v) in means.iter().zip(s2.iter()) {
        if !m.is_finite() || !v.is_finite() || *m <= 0.0 || *v <= 0.0 {
            continue;
        }
        let lm = m.ln();
        let lv = v.ln();
        let row = [1.0, lm, lm * lm];
        for i in 0..3 {
            xty[i] += row[i] * lv;
            for j in 0..3 {
                xtx[i][j] += row[i] * row[j];
            }
        }
        n_used += 1;
    }
    if n_used < 3 {
        return None;
    }
    let l = cholesky_lower(&xtx)?;
    // Forward-solve L z = Xᵀy.
    let mut z = [0.0_f64; 3];
    for i in 0..3 {
        let mut s = xty[i];
        for k in 0..i {
            s -= l[i][k] * z[k];
        }
        z[i] = s / l[i][i];
    }
    // Back-solve Lᵀ β = z.
    let mut beta = [0.0_f64; 3];
    for i in (0..3).rev() {
        let mut s = z[i];
        for k in (i + 1)..3 {
            s -= l[k][i] * beta[k];
        }
        beta[i] = s / l[i][i];
    }

    // Predict per feature.
    let mut trend = Vec::with_capacity(means.len());
    let mut fallback_sum = 0.0;
    let mut fallback_n = 0_usize;
    for v in s2 {
        if v.is_finite() && *v > 0.0 {
            fallback_sum += v.ln();
            fallback_n += 1;
        }
    }
    let fallback = if fallback_n > 0 {
        (fallback_sum / fallback_n as f64).exp()
    } else {
        f64::NAN
    };
    for m in means {
        if !m.is_finite() || *m <= 0.0 {
            trend.push(fallback);
            continue;
        }
        let lm = m.ln();
        let lv = beta[0] + beta[1] * lm + beta[2] * lm * lm;
        trend.push(lv.exp());
    }
    Some(trend)
}

/// Output of [`moderated_t`] — effect estimate, standard error, t,
/// two-sided p-value.
pub struct ModeratedT {
    pub effect: f64,
    pub se: f64,
    pub t: f64,
    pub p_value: f64,
}

/// Moderated t-statistic for a single contrast `c` on coefficient
/// vector `β`, given:
/// - `xtx_l`: lower-triangular Cholesky factor of `X'X`
/// - `s2_posterior`: eBayes-shrunken residual variance for this feature
/// - `df_total`: `df_res + df_prior`
///
/// `se(c·β) = sqrt(s²_post · c' (X'X)⁻¹ c)` where `c' (X'X)⁻¹ c = ||L⁻¹ c||²`
/// via one forward-solve against the cached Cholesky factor.
pub fn moderated_t(
    beta: &[f64],
    contrast: &[f64],
    xtx_l: &[Vec<f64>],
    s2_posterior: f64,
    df_total: f64,
) -> ModeratedT {
    let effect = beta
        .iter()
        .zip(contrast.iter())
        .map(|(b, c)| b * c)
        .sum::<f64>();
    // Forward-solve L z = c; then ||z||² = c' (X'X)⁻¹ c.
    let p = xtx_l.len();
    let mut z = vec![0.0_f64; p];
    for i in 0..p {
        let mut s = contrast[i];
        for k in 0..i {
            s -= xtx_l[i][k] * z[k];
        }
        z[i] = s / xtx_l[i][i];
    }
    let ctx_inv_c: f64 = z.iter().map(|v| v * v).sum();
    let se = (s2_posterior * ctx_inv_c).sqrt();
    if !se.is_finite() || se <= 0.0 {
        return ModeratedT {
            effect,
            se,
            t: f64::NAN,
            p_value: f64::NAN,
        };
    }
    let t = effect / se;
    let dist = match StudentsT::new(0.0, 1.0, df_total) {
        Ok(d) => d,
        Err(_) => {
            return ModeratedT {
                effect,
                se,
                t,
                p_value: f64::NAN,
            }
        }
    };
    let p_value = 2.0 * (1.0 - dist.cdf(t.abs()));
    ModeratedT {
        effect,
        se,
        t,
        p_value,
    }
}

/// Output of [`moderated_f`] — F statistic, numerator df `k`,
/// denominator df `df_total`, and p-value.
pub struct ModeratedF {
    pub f: f64,
    pub num_df: usize,
    pub den_df: f64,
    pub p_value: f64,
}

/// Moderated F-statistic across a contrast matrix `C` (rows = coefficients,
/// columns = contrasts), with eBayes-shrunken posterior variance `s²_post`
/// and `df_total = df_res + df_prior`.
///
/// Algorithm:
/// 1. `eff` = `Cᵀ · β` (length `k`).
/// 2. `M` = `Cᵀ · (X'X)⁻¹ · C` (`k × k`). Compute via `M = Zᵀ Z` where
///    `Z = L⁻¹ C` (forward-solve each column of `C` against `L`).
/// 3. Cholesky-factor `M`. If non-SPD, return `None` (rank-deficient
///    contrast matrix).
/// 4. Solve `M · w = eff`, then `F = effᵀ · w / (k · s²_post)`.
/// 5. p-value from `F_{k, df_total}`.
pub fn moderated_f(
    beta: &[f64],
    contrast_matrix: &[Vec<f64>],
    xtx_l: &[Vec<f64>],
    s2_posterior: f64,
    df_total: f64,
) -> Option<ModeratedF> {
    if contrast_matrix.is_empty() {
        return None;
    }
    let p = beta.len();
    if xtx_l.len() != p {
        return None;
    }
    if contrast_matrix.len() != p {
        return None;
    }
    let k = contrast_matrix[0].len();
    if k == 0 {
        return None;
    }
    for row in contrast_matrix {
        if row.len() != k {
            return None;
        }
    }
    if s2_posterior <= 0.0 || !s2_posterior.is_finite() {
        return None;
    }

    // eff_j = Σ_i C[i][j] · β[i]
    let mut eff = vec![0.0_f64; k];
    for j in 0..k {
        for i in 0..p {
            eff[j] += contrast_matrix[i][j] * beta[i];
        }
    }

    // Z = L⁻¹ C: for each column j of C, forward-solve.
    let mut z_mat = vec![vec![0.0_f64; k]; p];
    for j in 0..k {
        let mut z_col = vec![0.0_f64; p];
        for i in 0..p {
            let mut s = contrast_matrix[i][j];
            for q in 0..i {
                s -= xtx_l[i][q] * z_col[q];
            }
            z_col[i] = s / xtx_l[i][i];
        }
        for i in 0..p {
            z_mat[i][j] = z_col[i];
        }
    }

    // M = Zᵀ Z (k × k, symmetric).
    let mut m = vec![vec![0.0_f64; k]; k];
    for a in 0..k {
        for b in 0..k {
            let mut s = 0.0;
            for row in z_mat.iter().take(p) {
                s += row[a] * row[b];
            }
            m[a][b] = s;
        }
    }

    // Cholesky of M.
    let m_l = crate::de::cholesky_lower(&m)?;

    // Forward-solve L w_tmp = eff, then Lᵀ w = w_tmp.
    let mut w_tmp = vec![0.0_f64; k];
    for i in 0..k {
        let mut s = eff[i];
        for q in 0..i {
            s -= m_l[i][q] * w_tmp[q];
        }
        w_tmp[i] = s / m_l[i][i];
    }
    let mut w = vec![0.0_f64; k];
    for i in (0..k).rev() {
        let mut s = w_tmp[i];
        for q in (i + 1)..k {
            s -= m_l[q][i] * w[q];
        }
        w[i] = s / m_l[i][i];
    }

    let quad = eff.iter().zip(w.iter()).map(|(e, wv)| e * wv).sum::<f64>();
    let f = quad / (k as f64 * s2_posterior);
    if !f.is_finite() || f < 0.0 {
        return None;
    }
    let dist = FisherSnedecor::new(k as f64, df_total).ok()?;
    let p_value = 1.0 - dist.cdf(f);
    Some(ModeratedF {
        f,
        num_df: k,
        den_df: df_total,
        p_value,
    })
}

/// TREAT p-value (Phipson et al. 2013): test `H₀: |δ| ≤ lfc_threshold`
/// rather than the point null `H₀: δ = 0`. Implemented as the sum of
/// one-sided probabilities of `|β̂| > |t|·se` under the worst-case edge
/// of the composite null (δ = +lfc or δ = −lfc); by symmetry both give
/// the same p-value:
///
/// ```text
/// p = P(T > |t| − lfc/se | T ~ t_df) + P(T < −|t| − lfc/se | T ~ t_df)
/// ```
///
/// When `lfc = 0` this reduces to `2 · (1 − T_df(|t|))`, the standard
/// two-sided moderated-t p-value. The formula is symmetric in the sign
/// of `t` (a test statistic of `+t` and `−t` yield identical p-values),
/// as required for a two-sided test.
pub fn treat_p_value(t: f64, se: f64, df_total: f64, lfc_threshold: f64) -> f64 {
    if !t.is_finite() || !se.is_finite() || se <= 0.0 || !df_total.is_finite() || df_total <= 0.0
    {
        return f64::NAN;
    }
    let dist = match StudentsT::new(0.0, 1.0, df_total) {
        Ok(d) => d,
        Err(_) => return f64::NAN,
    };
    let abs_t = t.abs();
    let shift = lfc_threshold.abs() / se;
    let upper = 1.0 - dist.cdf(abs_t - shift);
    let lower = dist.cdf(-abs_t - shift);
    (upper + lower).clamp(0.0, 1.0)
}

/// Per-comparison options for [`limma_fit`].
#[derive(Clone, Copy, Debug)]
pub struct LimmaOptions {
    /// Fit the parametric mean-variance trend before `squeeze_var`.
    pub trend: bool,
    /// Use the Winsorized robust F-fit.
    pub robust: bool,
    /// TREAT minimum-effect threshold (on the log-fold-change scale).
    /// `0.0` produces the standard moderated-t p-value.
    pub lfc_threshold: f64,
}

/// Per-feature row of [`limma_fit`] output.
pub struct LimmaRow {
    pub effects: Vec<f64>,
    pub ses: Vec<f64>,
    pub t_stats: Vec<f64>,
    pub p_values: Vec<f64>,
    pub f_statistic: Option<f64>,
    pub f_p_value: Option<f64>,
    pub s2_sample: f64,
    pub s2_trend: Option<f64>,
    pub s2_posterior: f64,
    pub mean_abundance: f64,
    pub skipped: bool,
}

/// Output of [`limma_fit`] for one comparison (= one contrast matrix).
pub struct LimmaOutput {
    pub rows: Vec<LimmaRow>,
    pub df_prior: f64,
    /// Prior variance.
    ///
    /// When `options.trend` was `true`, the prior was fit on the
    /// per-feature ratios `s² / s²_trend`, so this value is on the
    /// **ratio scale** (dimensionless). The per-feature
    /// `LimmaRow.s2_posterior` has been rescaled back to the original
    /// variance scale via `s2_posterior_ratio * s2_trend[i]`; the
    /// top-level `s2_prior` is left in ratio units because there is no
    /// single trend value to multiply by. When `options.trend` was
    /// `false`, `s2_prior` is on the original variance scale.
    pub s2_prior: f64,
    pub df_residual: f64,
    pub df_total: f64,
    pub trend_fallback_used: bool,
}

/// High-level limma fit: run per-feature OLS via `atman-core::de::ols::ols`,
/// fit a shared Cholesky of X'X, optionally fit a mean-variance trend,
/// shrink per-feature variances via `squeeze_var` (robust if requested),
/// then compute moderated t per contrast and moderated F across all
/// contrasts.
///
/// `design` is `n × p`; `y` is `n_features × n`; `contrast_matrix` is
/// `p × k`.
///
/// Returns `None` when the design matrix is singular (no feature can be
/// fit). Features individually skipped by `ols::ols` have `skipped = true`
/// and NaN statistics.
pub fn limma_fit(
    design: &[Vec<f64>],
    y: &[Vec<f64>],
    contrast_matrix: &[Vec<f64>],
    options: LimmaOptions,
) -> Option<LimmaOutput> {
    use crate::de::{ols, OlsOutcome};

    if design.is_empty() || y.is_empty() || contrast_matrix.is_empty() {
        return None;
    }
    let n = design.len();
    let p = design[0].len();
    if design.iter().any(|r| r.len() != p) {
        return None;
    }
    if y.iter().any(|row| row.len() != n) {
        return None;
    }
    if contrast_matrix.len() != p {
        return None;
    }
    let k = contrast_matrix[0].len();
    if contrast_matrix.iter().any(|r| r.len() != k) {
        return None;
    }

    // Shared X'X and its Cholesky.
    let mut xtx = vec![vec![0.0_f64; p]; p];
    for row in design {
        for i in 0..p {
            for j in 0..p {
                xtx[i][j] += row[i] * row[j];
            }
        }
    }
    let xtx_l = cholesky_lower(&xtx)?;

    // Per-feature OLS.
    let mut s2_sample = Vec::with_capacity(y.len());
    let mut means = Vec::with_capacity(y.len());
    let mut betas = Vec::with_capacity(y.len());
    let mut skipped = Vec::with_capacity(y.len());
    let mut df_res = 0.0_f64;
    for feat_y in y {
        match ols(design, feat_y, 2) {
            OlsOutcome::Computed(fit) => {
                betas.push(fit.beta.clone());
                s2_sample.push(fit.sigma2);
                let mean: f64 = feat_y.iter().sum::<f64>() / n as f64;
                means.push(mean);
                skipped.push(false);
                df_res = fit.df;
            }
            OlsOutcome::Skipped { .. } => {
                betas.push(vec![f64::NAN; p]);
                s2_sample.push(f64::NAN);
                means.push(f64::NAN);
                skipped.push(true);
            }
        }
    }
    if skipped.iter().all(|b| *b) {
        return None;
    }
    if df_res <= 0.0 {
        return None;
    }

    // Trend.
    let mut trend_fallback_used = false;
    let s2_trend: Option<Vec<f64>> = if options.trend {
        match fit_parametric_trend(&means, &s2_sample) {
            Some(v) => Some(v),
            None => {
                trend_fallback_used = true;
                None
            }
        }
    } else {
        None
    };

    // Ratios s² / s²_trend (or s² itself if no trend).
    let ratios: Vec<f64> = s2_sample
        .iter()
        .enumerate()
        .map(|(i, v)| match &s2_trend {
            Some(t) if t[i].is_finite() && t[i] > 0.0 => v / t[i],
            _ => *v,
        })
        .collect();

    // eBayes prior fit.
    let prior = if options.robust {
        fit_f_dist_robust(&ratios, df_res)
            .map(|o| (o.df_prior, o.s2_prior))
            .or_else(|| fit_f_dist(&ratios, df_res))
    } else {
        fit_f_dist(&ratios, df_res)
    };
    // `trend_fallback_used` reflects whether the trend fit itself was
    // requested-and-unfit. A successful trend fit that subsequently finds
    // no evidence for residual prior variance (Var(log ratios) <=
    // trigamma(df_res/2)) is NOT a trend fallback — the trend was applied,
    // there's just nothing to shrink against. Let squeeze_var degrade
    // gracefully to pass-through in that case.
    let squeeze = squeeze_var(
        &ratios,
        df_res,
        prior,
    );
    let df_prior = squeeze.df_prior;
    let s2_prior = squeeze.s2_prior;
    let df_total = if df_prior.is_finite() {
        df_res + df_prior
    } else {
        df_res
    };

    // Per-feature statistics.
    let mut rows = Vec::with_capacity(y.len());
    for (i, beta) in betas.iter().enumerate() {
        if skipped[i] {
            rows.push(LimmaRow {
                effects: vec![f64::NAN; k],
                ses: vec![f64::NAN; k],
                t_stats: vec![f64::NAN; k],
                p_values: vec![f64::NAN; k],
                f_statistic: None,
                f_p_value: None,
                s2_sample: f64::NAN,
                s2_trend: s2_trend.as_ref().and_then(|t| if t[i].is_finite() { Some(t[i]) } else { None }),
                s2_posterior: f64::NAN,
                mean_abundance: f64::NAN,
                skipped: true,
            });
            continue;
        }
        // Scale s² back to the original scale (posterior of the trend-scaled
        // ratio times s²_trend, when a trend is present).
        let s2_post_feature = match &s2_trend {
            Some(t) if t[i].is_finite() && t[i] > 0.0 => squeeze.s2_posterior[i] * t[i],
            _ => squeeze.s2_posterior[i],
        };
        let mut effects = Vec::with_capacity(k);
        let mut ses = Vec::with_capacity(k);
        let mut t_stats = Vec::with_capacity(k);
        let mut p_values = Vec::with_capacity(k);
        for j in 0..k {
            let contrast: Vec<f64> = contrast_matrix.iter().map(|row| row[j]).collect();
            let mt = moderated_t(beta, &contrast, &xtx_l, s2_post_feature, df_total);
            let p_treat = treat_p_value(mt.t, mt.se, df_total, options.lfc_threshold);
            effects.push(mt.effect);
            ses.push(mt.se);
            t_stats.push(mt.t);
            p_values.push(p_treat);
        }
        let (f_stat, f_p) = if k >= 2 {
            match moderated_f(beta, contrast_matrix, &xtx_l, s2_post_feature, df_total) {
                Some(mf) => (Some(mf.f), Some(mf.p_value)),
                None => (None, None),
            }
        } else {
            (None, None)
        };
        rows.push(LimmaRow {
            effects,
            ses,
            t_stats,
            p_values,
            f_statistic: f_stat,
            f_p_value: f_p,
            s2_sample: s2_sample[i],
            s2_trend: s2_trend.as_ref().and_then(|t| if t[i].is_finite() { Some(t[i]) } else { None }),
            s2_posterior: s2_post_feature,
            mean_abundance: means[i],
            skipped: false,
        });
    }

    Some(LimmaOutput {
        rows,
        df_prior,
        s2_prior,
        df_residual: df_res,
        df_total,
        trend_fallback_used,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference values computed by Wolfram Alpha / R:digamma / R:trigamma / R:psigamma(.,2):
    //   digamma(1)        = -Euler-Mascheroni = -0.5772156649015329
    //   digamma(2)        = 1 - γ            =  0.4227843350984671
    //   digamma(10)       ≈  2.251752589066721
    //   trigamma(1)       = π²/6             ≈  1.6449340668482264
    //   trigamma(10)      ≈  0.10516633568168574
    //   tetragamma(1)     = -2 ζ(3)          ≈ -2.4041138063191886
    //   tetragamma(10)    ≈ -0.011049834970801975  (= -2ζ(3) + 2·Σ_{k=1..9} 1/k³)
    //   trigamma_inverse(π²/6)    → 1
    //   trigamma_inverse(0.10517) → ~10

    #[test]
    fn digamma_matches_tabulated_values() {
        assert!((digamma(1.0) - (-0.5772156649015329)).abs() < 1e-10);
        assert!((digamma(2.0) - 0.4227843350984671).abs() < 1e-10);
        assert!((digamma(10.0) - 2.251752589066721).abs() < 1e-10);
    }

    #[test]
    fn trigamma_matches_tabulated_values() {
        assert!((trigamma(1.0) - 1.6449340668482264).abs() < 1e-10);
        assert!((trigamma(10.0) - 0.10516633568168574).abs() < 1e-10);
    }

    #[test]
    fn tetragamma_matches_tabulated_values() {
        assert!((tetragamma(1.0) - (-2.4041138063191886)).abs() < 1e-10);
        assert!((tetragamma(10.0) - (-0.011049834970801975)).abs() < 1e-8);
    }

    #[test]
    fn trigamma_inverse_round_trips() {
        for x in [0.5_f64, 1.0, 2.0, 5.0, 10.0, 50.0] {
            let y = trigamma(x);
            let x_back = trigamma_inverse(y);
            assert!(
                (x_back - x).abs() < 1e-6,
                "trigamma_inverse({y}) = {x_back}, want {x}"
            );
        }
    }

    #[test]
    fn fit_f_dist_recovers_known_params_from_synthetic() {
        // Sample log(s²) from F(df_res=10, df_prior=5, s²_prior=4) and recover.
        // Synthetic: for each feature, draw s²_i from s²_prior * F(df_res, df_prior).
        // Since we can't RNG here, use a fixed sequence that hits the true mean/var.
        //
        // Precomputed sample of 5000 scaled-F(df=(10, 5)) * 4.0 draws with seed 20260418:
        // mean(log(s²)) should approach E[log(s²)] = digamma(df_res/2) - digamma(df_prior/2)
        //   + log(df_prior · s²_prior / df_res)
        // var(log(s²)) should approach trigamma(df_res/2) + trigamma(df_prior/2).
        //
        // We don't need actual synthetic draws: construct a sample whose first two
        // moments exactly match the theoretical targets, then check fit_f_dist
        // recovers (df_prior = 5, s²_prior = 4) within tolerance.
        let df_res = 10.0;
        let df_prior_true = 5.0;
        let s2_prior_true = 4.0;
        let mean_target = super::digamma(df_res / 2.0) - super::digamma(df_prior_true / 2.0)
            + (df_prior_true * s2_prior_true / df_res).ln();
        let var_target =
            super::trigamma(df_res / 2.0) + super::trigamma(df_prior_true / 2.0);
        // Build a 2-point sample that hits (mean_target, var_target) exactly.
        let n = 5000_f64;
        let spread = var_target.sqrt();
        let mut s2 = Vec::with_capacity(5000);
        for i in 0..5000 {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            let log_s2 = mean_target + sign * spread;
            s2.push(log_s2.exp());
        }
        let (df_prior, s2_prior) = super::fit_f_dist(&s2, df_res).expect("ok");
        assert!(
            (df_prior - df_prior_true).abs() < 0.5,
            "df_prior {df_prior}, want ~{df_prior_true}"
        );
        assert!(
            (s2_prior.ln() - s2_prior_true.ln()).abs() < 0.1,
            "s²_prior {s2_prior}, want ~{s2_prior_true}"
        );
        let _ = n;
    }

    #[test]
    fn squeeze_var_blends_sample_variance_toward_prior() {
        // With df_prior = df_res = 4 and s²_prior = 1.0, s²_post should
        // equal 0.5 * s²_sample + 0.5 * 1.0.
        let s2 = vec![0.0_f64, 1.0, 2.0, 4.0, 10.0];
        let df_res = 4.0;
        let result = super::squeeze_var(&s2, df_res, Some((4.0, 1.0)));
        let expected: Vec<f64> = s2.iter().map(|v| 0.5 * v + 0.5).collect();
        for (got, want) in result.s2_posterior.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
        }
        assert!((result.df_prior - 4.0).abs() < 1e-12);
        assert!((result.s2_prior - 1.0).abs() < 1e-12);
    }

    #[test]
    fn squeeze_var_with_none_prior_returns_sample_values_unchanged() {
        let s2 = vec![1.0, 2.0, 4.0];
        let out = super::squeeze_var(&s2, 4.0, None);
        for (got, want) in out.s2_posterior.iter().zip(s2.iter()) {
            assert_eq!(got, want);
        }
        assert!(out.df_prior.is_infinite()); // no shrinkage
    }

    #[test]
    fn fit_f_dist_robust_down_weights_planted_outlier() {
        // Construct a clean sample that recovers (df_prior=6, s²_prior=1),
        // then plant one huge outlier. Non-robust fit shifts; robust stays
        // close to the clean values.
        let df_res = 8.0;
        let df_prior_true = 6.0;
        let s2_prior_true = 1.0;
        let mean_target = super::digamma(df_res / 2.0) - super::digamma(df_prior_true / 2.0)
            + (df_prior_true * s2_prior_true / df_res).ln();
        let var_target =
            super::trigamma(df_res / 2.0) + super::trigamma(df_prior_true / 2.0);
        let spread = var_target.sqrt();
        let mut s2 = Vec::with_capacity(100);
        for i in 0..100 {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            s2.push((mean_target + sign * spread).exp());
        }
        // Plant one huge outlier.
        s2[0] = 1e6;

        let (df_non, s2_non) = super::fit_f_dist(&s2, df_res).expect("ok");
        let out = super::fit_f_dist_robust(&s2, df_res).expect("ok");
        // Robust fit is closer to the truth than the non-robust fit.
        assert!(
            (out.s2_prior - s2_prior_true).abs() < (s2_non - s2_prior_true).abs(),
            "robust s²_prior {} should be closer to {s2_prior_true} than non-robust {s2_non}",
            out.s2_prior,
        );
        let _ = df_non;
    }

    #[test]
    fn treat_at_zero_threshold_equals_standard_two_sided_t_p() {
        let t = 2.5;
        let se = 1.0;
        let df = 15.0;
        let treat_p = super::treat_p_value(t, se, df, 0.0);
        // Standard two-sided t-test p:
        let dist = statrs::distribution::StudentsT::new(0.0, 1.0, df).unwrap();
        let standard_p =
            2.0 * (1.0 - <statrs::distribution::StudentsT as statrs::distribution::ContinuousCDF<f64, f64>>::cdf(&dist, t.abs()));
        assert!((treat_p - standard_p).abs() < 1e-12);
    }

    #[test]
    fn treat_increases_p_as_threshold_grows() {
        let t = 2.5;
        let se = 0.5;
        let df = 15.0;
        let p_small = super::treat_p_value(t, se, df, 0.0);
        let p_large = super::treat_p_value(t, se, df, 0.5);
        assert!(p_large >= p_small);
    }

    #[test]
    fn treat_p_value_is_symmetric_in_sign_of_t() {
        // TREAT tests |δ|; by symmetry p(+t) must equal p(−t) for any
        // (se, df, lfc_threshold). Regression guard for a sign bug in an
        // earlier draft where negative `t` silently returned 1.0.
        for &(se, df, lfc) in &[(1.0_f64, 15.0_f64, 0.0_f64), (0.5, 15.0, 0.5), (2.0, 30.0, 0.2)] {
            let p_pos = super::treat_p_value(2.5, se, df, lfc);
            let p_neg = super::treat_p_value(-2.5, se, df, lfc);
            assert!(
                (p_pos - p_neg).abs() < 1e-12,
                "treat_p_value not symmetric: p(+2.5)={p_pos}, p(-2.5)={p_neg}"
            );
        }
    }

    #[test]
    fn moderated_f_equals_moderated_t_squared_for_single_contrast() {
        // 1-column contrast: F-statistic should equal t².
        let beta = vec![2.0, -1.0];
        let contrast_matrix = vec![vec![1.0], vec![-1.0]]; // contrast is (1, -1)
        let xtx_l = vec![vec![1.0, 0.0], vec![0.0, 1.0]]; // identity ⇒ L = I
        let s2_posterior = 4.0;
        let df_total = 10.0;
        let f_out = super::moderated_f(&beta, &contrast_matrix, &xtx_l, s2_posterior, df_total)
            .expect("ok");
        let t_out = super::moderated_t(
            &beta,
            &[1.0_f64, -1.0],
            &xtx_l,
            s2_posterior,
            df_total,
        );
        assert!((f_out.f - t_out.t * t_out.t).abs() < 1e-10);
    }

    #[test]
    fn moderated_f_populates_two_contrast_case() {
        let beta = vec![1.0, 2.0, 3.0];
        // Two contrasts: (β₁ - β₂), (β₂ - β₃).
        let contrast_matrix = vec![
            vec![1.0, 0.0],
            vec![-1.0, 1.0],
            vec![0.0, -1.0],
        ];
        let xtx_l = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let f_out = super::moderated_f(&beta, &contrast_matrix, &xtx_l, 1.0, 20.0).expect("ok");
        assert!(f_out.f > 0.0 && f_out.f.is_finite());
        assert!(f_out.p_value >= 0.0 && f_out.p_value <= 1.0);
        assert_eq!(f_out.num_df, 2);
    }

    #[test]
    fn moderated_t_reduces_to_standard_t_at_zero_df_prior() {
        // df_prior = 0 ⇒ s²_post = s²_sample ⇒ moderated t = classical t.
        let beta = vec![2.0];
        let contrast = vec![1.0];
        let xtx_l = vec![vec![1.0]]; // X'X = I, L = I ⇒ c' (X'X)⁻¹ c = 1
        let out = super::moderated_t(&beta, &contrast, &xtx_l, 4.0 /* s²_post */, 9.0 /* df */);
        // t = 2.0 / sqrt(4.0 * 1.0) = 1.0; p-value = 2 · (1 - T_9(1))
        assert!((out.t - 1.0).abs() < 1e-12);
        assert!(out.p_value > 0.0 && out.p_value < 1.0);
    }

    #[test]
    fn fit_parametric_trend_recovers_known_quadratic() {
        // Construct mean vs variance under log(s²) = 1.0 - 2.0 * log(μ) + 0.5 * log(μ)²
        let means: Vec<f64> = (1..=100).map(|i| i as f64 / 10.0 + 0.1).collect();
        let s2: Vec<f64> = means
            .iter()
            .map(|m| {
                let lm = m.ln();
                (1.0 - 2.0 * lm + 0.5 * lm * lm).exp()
            })
            .collect();
        let trend = super::fit_parametric_trend(&means, &s2).expect("ok");
        for (m, got) in means.iter().zip(trend.iter()) {
            let lm = m.ln();
            let want = (1.0 - 2.0 * lm + 0.5 * lm * lm).exp();
            assert!(
                (got / want - 1.0).abs() < 1e-6,
                "at μ={m} trend={got} want={want}"
            );
        }
    }

    #[test]
    fn limma_fit_on_trivial_two_group_design_produces_finite_stats() {
        // Design: 10 samples, 2 groups (first 5 are control, next 5 case).
        // 20 features: half have a planted effect, half are null.
        let n = 10;
        let design: Vec<Vec<f64>> = (0..n)
            .map(|i| if i < 5 { vec![1.0, 0.0] } else { vec![1.0, 1.0] })
            .collect();
        let mut y: Vec<Vec<f64>> = Vec::new();
        for f in 0..20 {
            let effect = if f < 10 { 1.0 } else { 0.0 };
            let mut row = Vec::new();
            for (i, _) in design.iter().enumerate() {
                let group = if i < 5 { 0.0 } else { 1.0 };
                // Deterministic pseudo-random noise.
                let noise = ((f * 31 + i * 7) as f64).sin();
                row.push(group * effect + 0.1 * noise);
            }
            y.push(row);
        }
        // Contrast: the second coefficient (case vs control).
        let contrast_matrix = vec![vec![0.0], vec![1.0]];
        let out = super::limma_fit(
            &design,
            &y,
            &contrast_matrix,
            super::LimmaOptions {
                trend: true,
                robust: true,
                lfc_threshold: 0.0,
            },
        )
        .expect("ok");
        assert_eq!(out.rows.len(), 20);
        for row in &out.rows {
            assert!(row.t_stats[0].is_finite() || row.skipped);
        }
        // Fallback flag must be false on this well-behaved fixture.
        assert!(!out.trend_fallback_used);
    }
}
