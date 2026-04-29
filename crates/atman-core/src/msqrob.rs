//! Peptide-level ridge-regularized linear mixed model for `atman de --test msqrob`.
//!
//! One model per protein. Rows are peptide-level measurements stacked over
//! samples and peptides belonging to that protein. The model is
//!
//! ```text
//! y_ij = X_i β + u_j + ε_ij
//!   u_j ~ N(0, σ²_peptide)
//!   ε_ij ~ N(0, σ²_res)
//! ```
//!
//! with an optional ridge penalty `λ · ‖β_nonintercept‖²` on the fixed
//! effects. The random intercept `u_j` is indexed by peptide. `β` holds
//! the condition effect (and any formula covariates). The penalty is the
//! msqrob2 "ridge on the condition effect" formulation generalized to all
//! non-intercept fixed-effect columns.
//!
//! Fitting:
//!
//! 1. REML profile: 1D golden-section search over `log(τ)` where
//!    `τ = σ²_peptide / σ²_res`, minimizing
//!    `df · log(σ²̂) + Σ_g log(1 + τ · m_g) + log|XᵀV⁻¹X + λP|`.
//! 2. At the optimum, compute penalized GLS β, standard errors from
//!    `σ²̂ · [(XᵀV⁻¹X + λP)⁻¹]_jj`, t, two-sided p-values on `n − p` df.
//!
//! The compound-symmetry inverse `V⁻¹` is applied analytically, so no
//! matrix larger than `p × p` is ever formed. All arithmetic is in f64
//! and deterministic — no randomness in the fit itself.

use crate::de::{cholesky_lower, SkipReason};
use statrs::distribution::{ContinuousCDF, StudentsT};

/// One msqrob result for one protein in one comparison.
#[derive(Debug, Clone, PartialEq)]
pub enum MsqrobOutcome {
    Computed(MsqrobFit),
    Skipped { reason: SkipReason, n: usize },
}

#[derive(Debug, Clone, PartialEq)]
pub struct MsqrobFit {
    /// Number of peptide-level observations used (rows of X).
    pub n: usize,
    /// Number of design columns (cols of X).
    pub p: usize,
    /// Distinct peptides observed for this protein.
    pub n_peptides: usize,
    /// Regression coefficients, indexed by design column.
    pub beta: Vec<f64>,
    /// Standard error for each coefficient from the penalized posterior.
    pub se: Vec<f64>,
    /// t-statistic = β / SE.
    pub t: Vec<f64>,
    /// Two-sided p-value on `n − p` df.
    pub p_value: Vec<f64>,
    /// Residual degrees of freedom, `n − p`.
    pub df: f64,
    /// Residual variance estimate `σ²̂ = RSS_weighted / df`.
    pub sigma2: f64,
    /// Fitted variance ratio `τ = σ²_peptide / σ²_res`.
    pub peptide_variance_ratio: f64,
    /// Diagonal of `(XᵀV⁻¹X + λP)⁻¹`. Exposed so the caller can
    /// recompute SE under an empirical-Bayes shrunk variance without
    /// refitting.
    pub inv_diag_penalized: Vec<f64>,
}

/// Apply limma-style empirical-Bayes variance shrinkage across a set
/// of msqrob fits. Feeds per-protein `sigma2` values into
/// `limma::fit_f_dist` and updates each fit's `sigma2`, `se`, `t`,
/// `df`, and `p_value` with the shrunk scale and augmented degrees
/// of freedom. This is the cross-protein variance-stabilization step
/// that msqrob2 applies via `squeezeVarRob` and brings our per-protein
/// inference into line with the reference implementation when the
/// number of fitted proteins is large enough for the F-distribution
/// prior to be well-determined.
///
/// Returns `Some((df_prior, s2_prior))` when a finite prior was fit
/// (≥ 2 positive finite variances and detectable excess over the
/// residual `trigamma(df_res/2)`), else `None` and leaves fits
/// untouched.
pub fn squeeze_variance(fits: &mut [MsqrobFit]) -> Option<(f64, f64)> {
    if fits.len() < 2 {
        return None;
    }
    let s2: Vec<f64> = fits.iter().map(|f| f.sigma2).collect();
    // Use the median residual df across fits as the reference df for
    // the F-dist fit. limma uses the common residual df when all
    // features share one; msqrob2 per-feature residual df's can vary,
    // so the median is a conservative representative.
    let mut dfs: Vec<f64> = fits.iter().map(|f| f.df).collect();
    dfs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let df_ref = dfs[dfs.len() / 2];
    if !(df_ref.is_finite() && df_ref > 0.0) {
        return None;
    }
    let (df_prior, s2_prior) = crate::limma::fit_f_dist(&s2, df_ref)?;
    if !(df_prior.is_finite() && df_prior > 0.0 && s2_prior.is_finite() && s2_prior >= 0.0) {
        return None;
    }
    for fit in fits.iter_mut() {
        let df_total = fit.df + df_prior;
        let s2_post = (fit.df * fit.sigma2 + df_prior * s2_prior) / df_total;
        if !(s2_post.is_finite() && s2_post >= 0.0 && df_total.is_finite() && df_total > 0.0) {
            continue;
        }
        let t_dist = match StudentsT::new(0.0, 1.0, df_total) {
            Ok(d) => d,
            Err(_) => continue,
        };
        for j in 0..fit.beta.len() {
            let inv = fit.inv_diag_penalized[j];
            if !(inv.is_finite() && inv >= 0.0) {
                continue;
            }
            let se = (s2_post * inv).sqrt();
            fit.se[j] = se;
            let t = if se == 0.0 {
                f64::NAN
            } else {
                fit.beta[j] / se
            };
            fit.t[j] = t;
            fit.p_value[j] = if t.is_finite() {
                (2.0 * t_dist.sf(t.abs())).clamp(0.0, 1.0)
            } else {
                f64::NAN
            };
        }
        fit.sigma2 = s2_post;
        fit.df = df_total;
    }
    Some((df_prior, s2_prior))
}

/// Fit the msqrob2-class model for one protein.
///
/// `design`: one row per peptide-level observation, columns indexed by
///   `[intercept, condition, covariate_1, ..., covariate_k]` (caller's
///   responsibility to include an intercept column at index 0).
/// `y`: log2 abundance per observation, same length as `design`.
/// `peptide_indices`: per-observation peptide id (dense, 0..n_peptides).
/// `ridge_lambda`: L2 penalty on non-intercept fixed-effect coefficients.
///   `0.0` disables the penalty and recovers the vanilla random-intercept
///   LMM. Negative values are refused.
/// `min_peptides`: refuse to fit if fewer distinct peptides are observed.
/// `min_samples`: refuse to fit if fewer total observations are present.
/// `robust`: when `true`, run Huber IRWLS on the GLS-whitened residuals
///   after the initial Gaussian fit. Matches `msqrob2::msqrob(robust =
///   TRUE)` which uses `MASS::rlm(psi = psi.huber)` with threshold
///   `k = 1.345` on the whitened residuals and reports the MAD-based
///   robust scale as the residual variance. The random-intercept
///   variance ratio τ is held at the REML-profile optimum from the
///   Gaussian fit — msqrob2's `lme4` + `rlm` pipeline freezes it the
///   same way.
pub fn fit_msqrob(
    design: &[Vec<f64>],
    y: &[f64],
    peptide_indices: &[usize],
    ridge_lambda: f64,
    min_peptides: usize,
    min_samples: usize,
    robust: bool,
) -> MsqrobOutcome {
    let n = design.len();
    if n != y.len() || n != peptide_indices.len() || n == 0 {
        return MsqrobOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }
    let p = design[0].len();
    if p == 0 || n < min_samples || n <= p {
        return MsqrobOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }
    for row in design {
        if row.len() != p || row.iter().any(|v| !v.is_finite()) {
            return MsqrobOutcome::Skipped {
                reason: SkipReason::NonFiniteInput,
                n,
            };
        }
    }
    if y.iter().any(|v| !v.is_finite()) {
        return MsqrobOutcome::Skipped {
            reason: SkipReason::NonFiniteInput,
            n,
        };
    }
    if !ridge_lambda.is_finite() || ridge_lambda < 0.0 {
        return MsqrobOutcome::Skipped {
            reason: SkipReason::NonFiniteInput,
            n,
        };
    }

    let peptide_index = peptide_indices_to_groups(peptide_indices);
    let n_peptides = peptide_index.len();
    if n_peptides < min_peptides {
        return MsqrobOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }

    // REML 1D search over log(τ). Initial grid → golden-section refine.
    let candidates: [f64; 8] = [-18.0, -12.0, -8.0, -4.0, 0.0, 4.0, 8.0, 12.0];
    let mut best_t = candidates[0];
    let mut best_obj = f64::INFINITY;
    for t in candidates {
        let obj = reml_objective(t.exp(), design, y, &peptide_index, ridge_lambda);
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
    let mut fc = reml_objective(c.exp(), design, y, &peptide_index, ridge_lambda);
    let mut fd = reml_objective(d.exp(), design, y, &peptide_index, ridge_lambda);
    for _ in 0..60 {
        if fc < fd {
            right = d;
            d = c;
            fd = fc;
            c = right - gr * (right - left);
            fc = reml_objective(c.exp(), design, y, &peptide_index, ridge_lambda);
        } else {
            left = c;
            c = d;
            fc = fd;
            d = left + gr * (right - left);
            fd = reml_objective(d.exp(), design, y, &peptide_index, ridge_lambda);
        }
    }
    let tau = ((left + right) / 2.0).exp();

    let Some(mut fit) = penalized_gls_fit(tau, design, y, &peptide_index, ridge_lambda) else {
        return MsqrobOutcome::Skipped {
            reason: SkipReason::ZeroVariance,
            n,
        };
    };
    let df = (n - p) as f64;
    if !(df.is_finite() && df > 0.0) {
        return MsqrobOutcome::Skipped {
            reason: SkipReason::InsufficientPairs,
            n,
        };
    }
    let mut sigma2 = fit.rss_weighted / df;
    if !(sigma2.is_finite() && sigma2 >= 0.0) {
        return MsqrobOutcome::Skipped {
            reason: SkipReason::NonFiniteInput,
            n,
        };
    }
    // Huber IRWLS on the GLS-whitened residuals when robust = true.
    // Iterates β via weighted normal equations with Huber ψ threshold
    // k = 1.345 and MAD-based scale σ̂ = MAD × 1.4826. The random-
    // intercept variance ratio τ is frozen at the Gaussian REML
    // optimum — msqrob2 does the same via `lme4` upstream of `rlm`.
    if robust {
        let mut beta = fit.beta.clone();
        let mut wx = vec![vec![0.0_f64; p]; n];
        let mut wy = vec![0.0_f64; n];
        apply_compound_symmetry_inverse(tau, design, y, &peptide_index, &mut wx, &mut wy);
        const K: f64 = 1.345;
        const MAX_ITER: usize = 50;
        const TOL: f64 = 1e-6;
        let mut robust_sigma = f64::NAN;
        for _ in 0..MAX_ITER {
            let residuals: Vec<f64> = wy
                .iter()
                .zip(wx.iter())
                .map(|(yi, xr)| yi - xr.iter().zip(beta.iter()).map(|(a, b)| a * b).sum::<f64>())
                .collect();
            let mut abs_r: Vec<f64> = residuals.iter().map(|r| r.abs()).collect();
            abs_r.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let mad = if abs_r.is_empty() {
                0.0
            } else if abs_r.len().is_multiple_of(2) {
                0.5 * (abs_r[abs_r.len() / 2 - 1] + abs_r[abs_r.len() / 2])
            } else {
                abs_r[abs_r.len() / 2]
            };
            let sigma = mad * 1.4826;
            if !(sigma > 0.0 && sigma.is_finite()) {
                break;
            }
            robust_sigma = sigma;
            let weights: Vec<f64> = residuals
                .iter()
                .map(|r| {
                    let u = r.abs() / sigma;
                    if u > K {
                        K / u
                    } else {
                        1.0
                    }
                })
                .collect();
            let mut xtwx = vec![vec![0.0_f64; p]; p];
            let mut xtwy = vec![0.0_f64; p];
            for i in 0..n {
                let wi = weights[i];
                for j in 0..p {
                    xtwy[j] += wi * wx[i][j] * wy[i];
                    for kk in 0..p {
                        xtwx[j][kk] += wi * wx[i][j] * wx[i][kk];
                    }
                }
            }
            for (j, row) in xtwx.iter_mut().enumerate().skip(1) {
                row[j] += ridge_lambda;
            }
            let Some(l) = cholesky_lower(&xtwx) else {
                break;
            };
            let beta_new = solve_cholesky(&l, &xtwy);
            let delta = beta
                .iter()
                .zip(beta_new.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            beta = beta_new;
            fit.inv_diag_penalized = inverse_diag_from_cholesky(&l);
            if delta < TOL {
                break;
            }
        }
        fit.beta = beta;
        if robust_sigma.is_finite() && robust_sigma > 0.0 {
            sigma2 = robust_sigma * robust_sigma;
        }
    }
    let se: Vec<f64> = fit
        .inv_diag_penalized
        .iter()
        .map(|v| (sigma2 * v).sqrt())
        .collect();
    let t: Vec<f64> = fit
        .beta
        .iter()
        .zip(se.iter())
        .map(|(b, s)| if *s == 0.0 { f64::NAN } else { b / s })
        .collect();
    let t_dist = match StudentsT::new(0.0, 1.0, df) {
        Ok(d) => d,
        Err(_) => {
            return MsqrobOutcome::Skipped {
                reason: SkipReason::InsufficientPairs,
                n,
            };
        }
    };
    let p_value: Vec<f64> = t
        .iter()
        .map(|&ti| {
            if ti.is_finite() {
                (2.0 * t_dist.sf(ti.abs())).clamp(0.0, 1.0)
            } else {
                f64::NAN
            }
        })
        .collect();

    MsqrobOutcome::Computed(MsqrobFit {
        n,
        p,
        n_peptides,
        beta: fit.beta,
        se,
        t,
        p_value,
        df,
        sigma2,
        peptide_variance_ratio: tau,
        inv_diag_penalized: fit.inv_diag_penalized,
    })
}

struct PenalizedGls {
    beta: Vec<f64>,
    /// Diagonal of `(XᵀV⁻¹X + λP)⁻¹`, used for SE.
    inv_diag_penalized: Vec<f64>,
    /// Residual quadratic form `(y − Xβ)ᵀ V⁻¹ (y − Xβ)`.
    rss_weighted: f64,
    /// log|XᵀV⁻¹X + λP|, used by REML.
    logdet_penalized: f64,
}

/// Build the peptide groups. Each inner `Vec<usize>` is the observation
/// indices belonging to one peptide. Groups with zero observations are
/// dropped.
fn peptide_indices_to_groups(peptide_indices: &[usize]) -> Vec<Vec<usize>> {
    if peptide_indices.is_empty() {
        return Vec::new();
    }
    let max_id = *peptide_indices.iter().max().unwrap();
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); max_id + 1];
    for (obs, &p) in peptide_indices.iter().enumerate() {
        buckets[p].push(obs);
    }
    buckets.into_iter().filter(|b| !b.is_empty()).collect()
}

/// Apply `V⁻¹` to `design` and `y` in-place into `wx`, `wy` where
/// `V = I + τ · ZZᵀ` (so `V⁻¹` is compound-symmetric within peptide).
fn apply_compound_symmetry_inverse(
    tau: f64,
    design: &[Vec<f64>],
    y: &[f64],
    peptide_index: &[Vec<usize>],
    wx: &mut [Vec<f64>],
    wy: &mut [f64],
) {
    let p = design[0].len();
    for idx in peptide_index {
        let m = idx.len() as f64;
        let factor = tau / (1.0 + tau * m);
        let mut sum_y = 0.0;
        let mut sum_x = vec![0.0; p];
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

/// Apply `V⁻¹` to a residual vector.
fn apply_compound_symmetry_inverse_vec(
    tau: f64,
    values: &[f64],
    peptide_index: &[Vec<usize>],
    out: &mut [f64],
) {
    for idx in peptide_index {
        let m = idx.len() as f64;
        let factor = tau / (1.0 + tau * m);
        let sum: f64 = idx.iter().map(|&i| values[i]).sum();
        for &i in idx {
            out[i] = values[i] - factor * sum;
        }
    }
}

/// Forward then back substitution for a Cholesky factor `L` where
/// `L Lᵀ = A`. Returns `A⁻¹ · rhs`.
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

/// Diagonal of the inverse matrix given its Cholesky factor.
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

/// Penalized GLS fit at a fixed variance ratio `τ` and ridge `λ`.
fn penalized_gls_fit(
    tau: f64,
    design: &[Vec<f64>],
    y: &[f64],
    peptide_index: &[Vec<usize>],
    ridge_lambda: f64,
) -> Option<PenalizedGls> {
    let n = design.len();
    let p = design[0].len();
    let mut wx = vec![vec![0.0; p]; n];
    let mut wy = vec![0.0; n];
    apply_compound_symmetry_inverse(tau, design, y, peptide_index, &mut wx, &mut wy);

    // Build XᵀV⁻¹X and XᵀV⁻¹y from (X, wx, wy).
    let mut xtvix = vec![vec![0.0; p]; p];
    let mut xtviy = vec![0.0; p];
    for i in 0..n {
        for j in 0..p {
            xtviy[j] += design[i][j] * wy[i];
            for k in 0..p {
                xtvix[j][k] += design[i][j] * wx[i][k];
            }
        }
    }
    // Add ridge to non-intercept diagonal. Convention: column 0 is the
    // intercept and is never penalized.
    for (j, row) in xtvix.iter_mut().enumerate().skip(1) {
        row[j] += ridge_lambda;
    }

    let l = cholesky_lower(&xtvix)?;
    let beta = solve_cholesky(&l, &xtviy);
    let inv_diag = inverse_diag_from_cholesky(&l);
    let logdet = 2.0
        * l.iter()
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
    let mut wres = vec![0.0; n];
    apply_compound_symmetry_inverse_vec(tau, &residuals, peptide_index, &mut wres);
    let rss_weighted = residuals
        .iter()
        .zip(wres.iter())
        .map(|(r, wr)| r * wr)
        .sum();

    Some(PenalizedGls {
        beta,
        inv_diag_penalized: inv_diag,
        rss_weighted,
        logdet_penalized: logdet,
    })
}

/// REML-style profile objective in `τ`. Lower is better.
fn reml_objective(
    tau: f64,
    design: &[Vec<f64>],
    y: &[f64],
    peptide_index: &[Vec<usize>],
    ridge_lambda: f64,
) -> f64 {
    let Some(fit) = penalized_gls_fit(tau, design, y, peptide_index, ridge_lambda) else {
        return f64::INFINITY;
    };
    let n = design.len();
    let p = design[0].len();
    if n <= p || fit.rss_weighted <= 0.0 || !fit.rss_weighted.is_finite() {
        return f64::INFINITY;
    }
    let df = (n - p) as f64;
    let sigma2 = fit.rss_weighted / df;
    let logdet_v0: f64 = peptide_index
        .iter()
        .map(|idx| (1.0 + tau * idx.len() as f64).ln())
        .sum();
    df * sigma2.ln() + logdet_v0 + fit.logdet_penalized
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    /// When all peptides agree (τ → 0, no peptide effect), msqrob with
    /// λ=0 collapses to plain OLS. Check β_condition recovers the
    /// injected effect on a small designed fixture.
    #[test]
    fn recovers_injected_condition_effect_no_peptide_noise() {
        // 4 peptides × 6 samples = 24 obs. Group A (samples 0..3) has
        // abundance 10 per peptide, group B (samples 3..6) has 12. No
        // peptide-specific shift, no noise.
        let mut design = Vec::new();
        let mut y = Vec::new();
        let mut peptide = Vec::new();
        for pep in 0..4usize {
            for sample in 0..6usize {
                let cond = if sample < 3 { 0.0 } else { 1.0 };
                design.push(vec![1.0, cond]);
                y.push(if sample < 3 { 10.0 } else { 12.0 });
                peptide.push(pep);
            }
        }
        let out = fit_msqrob(&design, &y, &peptide, 0.0, 2, 4, false);
        match out {
            MsqrobOutcome::Computed(fit) => {
                assert_eq!(fit.n, 24);
                assert_eq!(fit.n_peptides, 4);
                // Intercept = group-A mean = 10; condition = +2.
                assert!(approx(fit.beta[0], 10.0, 1e-8), "beta[0] = {}", fit.beta[0]);
                assert!(approx(fit.beta[1], 2.0, 1e-8), "beta[1] = {}", fit.beta[1]);
            }
            _ => panic!("expected Computed"),
        }
    }

    /// REML profile should recover a known peptide variance ratio.
    /// Peptide-specific intercepts spread ±0.5; residual noise is 0.
    /// With zero residual variance we expect τ → very large (peptide
    /// variance dominates). Check at least that τ > 1 (peptide effect
    /// detected as larger than residual).
    #[test]
    fn reml_detects_large_peptide_variance_ratio() {
        let peptide_offsets = [-0.5, -0.2, 0.1, 0.6];
        let mut design = Vec::new();
        let mut y = Vec::new();
        let mut peptide = Vec::new();
        for (pep, offset) in peptide_offsets.iter().enumerate() {
            for sample in 0..6usize {
                let cond = if sample < 3 { 0.0 } else { 1.0 };
                design.push(vec![1.0, cond]);
                // Baseline 10 + peptide offset + 1.5 * condition
                y.push(10.0 + offset + if sample < 3 { 0.0 } else { 1.5 });
                peptide.push(pep);
            }
        }
        let out = fit_msqrob(&design, &y, &peptide, 0.0, 2, 4, false);
        match out {
            MsqrobOutcome::Computed(fit) => {
                assert!(
                    fit.peptide_variance_ratio > 1.0,
                    "expected τ > 1, got {}",
                    fit.peptide_variance_ratio
                );
                // Condition effect still recovered cleanly.
                assert!(approx(fit.beta[1], 1.5, 1e-4), "beta[1] = {}", fit.beta[1]);
            }
            _ => panic!("expected Computed"),
        }
    }

    /// Ridge shrinks β toward 0 when λ > 0. Check that the condition
    /// effect is attenuated relative to λ=0 on identical data.
    #[test]
    fn ridge_shrinks_condition_effect_toward_zero() {
        let mut design = Vec::new();
        let mut y = Vec::new();
        let mut peptide = Vec::new();
        for pep in 0..3usize {
            for sample in 0..4usize {
                let cond = if sample < 2 { 0.0 } else { 1.0 };
                design.push(vec![1.0, cond]);
                y.push(if sample < 2 { 5.0 } else { 8.0 });
                peptide.push(pep);
            }
        }
        let out0 = fit_msqrob(&design, &y, &peptide, 0.0, 2, 4, false);
        let out_lambda = fit_msqrob(&design, &y, &peptide, 5.0, 2, 4, false);
        match (out0, out_lambda) {
            (MsqrobOutcome::Computed(f0), MsqrobOutcome::Computed(fl)) => {
                assert!(
                    fl.beta[1].abs() < f0.beta[1].abs(),
                    "ridge should shrink: λ=0 β={}, λ=5 β={}",
                    f0.beta[1],
                    fl.beta[1]
                );
                assert!(fl.beta[1] > 0.0, "direction preserved");
            }
            _ => panic!("expected two Computed outcomes"),
        }
    }

    /// Intercept is never penalized. Check that fitting the same data
    /// with varying λ does not change β[0] meaningfully more than β[1].
    #[test]
    fn intercept_is_not_penalized() {
        let mut design = Vec::new();
        let mut y = Vec::new();
        let mut peptide = Vec::new();
        for pep in 0..3usize {
            for sample in 0..6usize {
                let cond = if sample < 3 { 0.0 } else { 1.0 };
                design.push(vec![1.0, cond]);
                y.push(if sample < 3 { 10.0 } else { 11.5 });
                peptide.push(pep);
            }
        }
        let f0 = match fit_msqrob(&design, &y, &peptide, 0.0, 2, 4, false) {
            MsqrobOutcome::Computed(f) => f,
            _ => panic!(),
        };
        let fl = match fit_msqrob(&design, &y, &peptide, 10.0, 2, 4, false) {
            MsqrobOutcome::Computed(f) => f,
            _ => panic!(),
        };
        let d_intercept = (f0.beta[0] - fl.beta[0]).abs();
        let d_condition = (f0.beta[1] - fl.beta[1]).abs();
        assert!(
            d_condition > d_intercept,
            "condition should move more than intercept under ridge: \
             Δintercept={}, Δcondition={}",
            d_intercept,
            d_condition
        );
    }

    /// Fewer peptides than `min_peptides` → skipped.
    #[test]
    fn refuses_below_min_peptides() {
        let design = vec![
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ];
        let y = vec![1.0, 2.0, 1.1, 2.1];
        let peptide = vec![0, 0, 0, 0];
        let out = fit_msqrob(&design, &y, &peptide, 0.0, 2, 4, false);
        assert!(matches!(
            out,
            MsqrobOutcome::Skipped {
                reason: SkipReason::InsufficientPairs,
                ..
            }
        ));
    }

    /// Non-finite input → Skipped(NonFiniteInput).
    #[test]
    fn refuses_nonfinite_y() {
        let design = vec![vec![1.0, 0.0]; 8];
        let mut y = vec![1.0; 8];
        y[3] = f64::NAN;
        let peptide = vec![0, 0, 1, 1, 2, 2, 3, 3];
        let out = fit_msqrob(&design, &y, &peptide, 0.0, 2, 4, false);
        assert!(matches!(
            out,
            MsqrobOutcome::Skipped {
                reason: SkipReason::NonFiniteInput,
                ..
            }
        ));
    }

    /// Determinism: two runs on identical input produce byte-equal
    /// numerical output.
    #[test]
    fn fit_is_deterministic() {
        let mut design = Vec::new();
        let mut y = Vec::new();
        let mut peptide = Vec::new();
        for pep in 0..4usize {
            for sample in 0..6usize {
                let cond = if sample < 3 { 0.0 } else { 1.0 };
                design.push(vec![1.0, cond]);
                y.push(
                    10.0 + (pep as f64) * 0.3
                        + if sample < 3 { 0.0 } else { 1.2 }
                        + ((sample + pep) as f64).sin() * 0.05,
                );
                peptide.push(pep);
            }
        }
        let a = fit_msqrob(&design, &y, &peptide, 0.5, 2, 4, false);
        let b = fit_msqrob(&design, &y, &peptide, 0.5, 2, 4, false);
        assert_eq!(a, b);
    }

    /// Huber IRWLS (`robust = true`) downweights outliers: inject a
    /// single peptide-level outlier into a clean fixture and verify
    /// the condition coefficient moves toward the no-outlier fit.
    #[test]
    fn robust_huber_downweights_outlier_peptide() {
        let peptide_offsets = [0.0_f64, 0.0, 0.0, 0.0];
        let true_effect = 2.0_f64;
        let mut design = Vec::new();
        let mut y_clean = Vec::new();
        let mut peptide = Vec::new();
        for (pep, off) in peptide_offsets.iter().enumerate() {
            for sample in 0..8usize {
                let cond = if sample < 4 { 0.0 } else { 1.0 };
                design.push(vec![1.0, cond]);
                y_clean.push(
                    10.0 + off + true_effect * cond + ((sample + pep) as f64 * 0.13).sin() * 0.05,
                );
                peptide.push(pep);
            }
        }
        let mut y_outlier = y_clean.clone();
        // Single heavy outlier on one peptide in group A.
        y_outlier[0] += 12.0;
        let clean_robust = match fit_msqrob(&design, &y_clean, &peptide, 0.0, 2, 4, true) {
            MsqrobOutcome::Computed(f) => f,
            _ => panic!("clean-data fit failed"),
        };
        let outlier_gauss = match fit_msqrob(&design, &y_outlier, &peptide, 0.0, 2, 4, false) {
            MsqrobOutcome::Computed(f) => f,
            _ => panic!("gaussian-outlier fit failed"),
        };
        let outlier_robust = match fit_msqrob(&design, &y_outlier, &peptide, 0.0, 2, 4, true) {
            MsqrobOutcome::Computed(f) => f,
            _ => panic!("robust-outlier fit failed"),
        };
        // Robust fit should pull the condition coefficient back
        // toward the clean truth compared to the Gaussian fit on
        // the same outlier-contaminated data.
        let gauss_err = (outlier_gauss.beta[1] - clean_robust.beta[1]).abs();
        let robust_err = (outlier_robust.beta[1] - clean_robust.beta[1]).abs();
        assert!(
            robust_err < gauss_err,
            "robust beta[1]={} should be closer to clean beta[1]={} than gaussian beta[1]={} \
             (robust_err={}, gauss_err={})",
            outlier_robust.beta[1],
            clean_robust.beta[1],
            outlier_gauss.beta[1],
            robust_err,
            gauss_err
        );
    }

    /// Negative ridge λ is refused.
    #[test]
    fn refuses_negative_lambda() {
        let design = vec![vec![1.0, 0.0]; 8];
        let y = vec![1.0; 8];
        let peptide = vec![0, 0, 1, 1, 2, 2, 3, 3];
        let out = fit_msqrob(&design, &y, &peptide, -1.0, 2, 4, false);
        assert!(matches!(
            out,
            MsqrobOutcome::Skipped {
                reason: SkipReason::NonFiniteInput,
                ..
            }
        ));
    }
}
