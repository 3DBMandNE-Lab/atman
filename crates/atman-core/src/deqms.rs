//! DEqMS peptide-count-weighted variance prior for `atman de
//! --test limma --peptide-metadata ...`.
//!
//! DEqMS (Zhu et al. 2020, MCP) extends limma's empirical-Bayes
//! variance shrinkage by conditioning the prior on per-protein
//! peptide count: proteins with more peptides get a tighter prior
//! because more peptides make the residual-variance estimate less
//! noisy. The implementation mirrors `limma::spectraCounteBayes`:
//!
//! 1. Tricube-kernel moving average of per-feature `log(σ²_residual)`
//!    against `log(peptide_count + 1)` → per-feature trend variance
//!    `s²_trend[i]`.
//! 2. Trend-scaled ratios `r[i] = σ²_residual[i] / s²_trend[i]`.
//! 3. Fit a scaled-F prior to the ratios (re-using
//!    `limma::fit_f_dist`) to obtain `(df_prior, s²_prior_ratio)`.
//! 4. Per-feature posterior on the ratio scale, multiplied back by
//!    `s²_trend[i]` for the original variance scale.
//!
//! The smoother follows `statmod::tricubeMovingAverage`: rank-order
//! observations on the covariate, form a window of `ceil(span · n)`
//! neighbours on each side, and take a tricube-weighted mean.
//!
//! Entry point: [`deqms_shrink`] accepts per-feature residual
//! variances and peptide counts and returns per-feature shrunk
//! variances plus the scalar `(df_prior, s²_prior_ratio)` prior.
//! Skipped features (NaN σ²) are passed through unchanged.
//!
//! Reference: Zhu Y et al., "DEqMS: A Method for Accurate Variance
//! Estimation in Differential Protein Expression Analysis," Molecular
//! & Cellular Proteomics (2020).

use crate::limma::fit_f_dist;

/// Per-feature DEqMS output aligned with the input index.
#[derive(Debug, Clone, PartialEq)]
pub struct DeqmsShrinkage {
    /// Tricube-smoothed trend variance per feature.
    pub s2_trend: Vec<f64>,
    /// Posterior variance per feature on the original variance scale.
    pub s2_posterior: Vec<f64>,
    /// Scalar prior degrees of freedom from the F-dist fit on the
    /// trend-scaled ratios. `+Inf` means no shrinkage was applied
    /// (the ratios showed no residual heterogeneity, or the smoother
    /// degenerated).
    pub df_prior: f64,
    /// Scalar prior variance on the **ratio scale** (dimensionless).
    /// The per-feature posterior on that scale is implicit; callers
    /// use `s2_posterior` above for the original scale.
    pub s2_prior_ratio: f64,
    /// True when the smoother could not produce a usable trend and
    /// the function fell through to a flat prior (equivalent to
    /// `fit_f_dist` on the raw variances).
    pub trend_fallback_used: bool,
}

/// Run the full DEqMS shrinkage step.
///
/// `s2_sample`: per-feature residual variance σ²̂. NaN values mark
///   skipped features and are propagated as NaN in the output.
/// `peptide_counts`: per-feature integer peptide counts, aligned
///   with `s2_sample`.
/// `df_residual`: residual degrees of freedom shared across features
///   (limma's per-feature OLS uses `n − p`, the same for every
///   feature because the design is shared).
/// `span`: tricube smoother bandwidth, expressed as a fraction of
///   the observation count. `0.5` matches `limma::squeezeVar`'s
///   default.
pub fn deqms_shrink(
    s2_sample: &[f64],
    peptide_counts: &[u32],
    df_residual: f64,
    span: f64,
) -> DeqmsShrinkage {
    assert_eq!(s2_sample.len(), peptide_counts.len());
    let n = s2_sample.len();
    let mut trend = vec![f64::NAN; n];
    let mut posterior = vec![f64::NAN; n];

    // Log1p covariate: `log(count + 1)` over finite, positive-variance
    // features. `fit_deqms_trend` fits the smoother on those only
    // and scatters the result back to the full-length `trend`.
    let usable: Vec<usize> = s2_sample
        .iter()
        .enumerate()
        .filter(|(_, v)| v.is_finite() && **v > 0.0)
        .map(|(i, _)| i)
        .collect();
    if usable.len() < 2 {
        return DeqmsShrinkage {
            s2_trend: trend,
            s2_posterior: posterior,
            df_prior: f64::INFINITY,
            s2_prior_ratio: f64::NAN,
            trend_fallback_used: true,
        };
    }
    let xs: Vec<f64> = usable
        .iter()
        .map(|&i| ((peptide_counts[i] as f64) + 1.0).ln())
        .collect();
    let ys_log: Vec<f64> = usable.iter().map(|&i| s2_sample[i].ln()).collect();
    let fitted_log = tricube_moving_average(&xs, &ys_log, span);
    let mut trend_fallback_used = false;
    if fitted_log.iter().any(|v| !v.is_finite()) {
        trend_fallback_used = true;
    }
    for (local, &global) in usable.iter().enumerate() {
        let v = fitted_log[local];
        if v.is_finite() {
            trend[global] = v.exp();
        }
    }

    // Trend-scaled ratios for the F-dist fit.
    let ratios: Vec<f64> = (0..n)
        .filter_map(|i| {
            if !trend[i].is_finite() || trend[i] <= 0.0 {
                return None;
            }
            if !s2_sample[i].is_finite() {
                return None;
            }
            Some(s2_sample[i] / trend[i])
        })
        .collect();

    let prior = fit_f_dist(&ratios, df_residual);
    let (df_prior, s2_prior_ratio) = match prior {
        Some((d, s)) if d.is_finite() && s.is_finite() && s > 0.0 => (d, s),
        _ => {
            trend_fallback_used = true;
            // Fall back to no shrinkage: posterior = sample when the
            // prior cannot be fit.
            for (i, &v) in s2_sample.iter().enumerate() {
                posterior[i] = v;
            }
            return DeqmsShrinkage {
                s2_trend: trend,
                s2_posterior: posterior,
                df_prior: f64::INFINITY,
                s2_prior_ratio: f64::NAN,
                trend_fallback_used,
            };
        }
    };

    // Per-feature posterior on the ratio scale, scaled back to the
    // original variance scale.
    for i in 0..n {
        if !trend[i].is_finite() || !s2_sample[i].is_finite() {
            continue;
        }
        let ratio_i = s2_sample[i] / trend[i];
        let post_ratio =
            (df_residual * ratio_i + df_prior * s2_prior_ratio) / (df_residual + df_prior);
        posterior[i] = post_ratio * trend[i];
    }

    DeqmsShrinkage {
        s2_trend: trend,
        s2_posterior: posterior,
        df_prior,
        s2_prior_ratio,
        trend_fallback_used,
    }
}

/// Tricube-kernel moving-average smoother, matching
/// `statmod::tricubeMovingAverage`.
///
/// `x` and `y` must have the same length; both are unsorted on input
/// and the output is returned in the original order.
///
/// For each point `i`, the smoother finds the rank-ordered window of
/// `w = ceil(span · n)` nearest neighbours on each side, computes
/// tricube weights `w_j = (1 − (|r_j − r_i| / max_rank_dist)³)³` on
/// them, and returns the weighted mean of the corresponding `y`
/// values. `span = 0.5` matches the `limma::squeezeVar` default.
///
/// Ties in `x` are handled by stable sort on the original index.
pub fn tricube_moving_average(x: &[f64], y: &[f64], span: f64) -> Vec<f64> {
    assert_eq!(x.len(), y.len());
    let n = x.len();
    let mut out = vec![f64::NAN; n];
    if n == 0 {
        return out;
    }
    if n == 1 {
        out[0] = y[0];
        return out;
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        x[a].partial_cmp(&x[b])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.cmp(&b))
    });

    // Window half-width in rank units. `ceil(span · n)` matches
    // statmod's default for `span = 0.5`; minimum 1 to keep the
    // smoother well-defined on tiny samples.
    let half = ((span * n as f64).ceil() as usize).max(1);

    for rank_i in 0..n {
        let lo = rank_i.saturating_sub(half);
        let hi = (rank_i + half).min(n - 1);
        let max_dist = (rank_i - lo).max(hi - rank_i).max(1) as f64;
        let mut num = 0.0_f64;
        let mut den = 0.0_f64;
        for rank_j in lo..=hi {
            let d = ((rank_j as isize - rank_i as isize).unsigned_abs()) as f64;
            let u = (d / max_dist).min(1.0);
            let w = (1.0 - u * u * u).max(0.0).powi(3);
            let j_orig = order[rank_j];
            num += w * y[j_orig];
            den += w;
        }
        let fitted = if den > 0.0 { num / den } else { f64::NAN };
        out[order[rank_i]] = fitted;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tricube smoother on constant y reproduces the constant
    /// regardless of span.
    #[test]
    fn tricube_constant_y_is_constant_output() {
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![7.5; 6];
        for span in [0.1_f64, 0.3, 0.5, 0.9] {
            let fitted = tricube_moving_average(&x, &y, span);
            for v in &fitted {
                assert!((v - 7.5).abs() < 1e-12, "span {span} drifted to {v}");
            }
        }
    }

    /// Tricube smoother on a linearly-ordered input preserves monotone
    /// behaviour: fitted[i+1] ≥ fitted[i] when y is monotone in x.
    #[test]
    fn tricube_preserves_monotonicity_on_increasing_input() {
        let x: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let y: Vec<f64> = (0..20).map(|i| (i as f64).ln_1p()).collect();
        let fitted = tricube_moving_average(&x, &y, 0.3);
        for i in 1..fitted.len() {
            assert!(
                fitted[i] >= fitted[i - 1] - 1e-9,
                "monotonicity broken at {i}: {fitted:?}"
            );
        }
    }

    /// Tricube output matches hand-computed values on a tiny
    /// deterministic grid. This pins the kernel normalization.
    #[test]
    fn tricube_exact_center_point_matches_hand_value() {
        // 5 points, span = 1.0 → every point uses all 5 as neighbours.
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let fitted = tricube_moving_average(&x, &y, 1.0);
        // At the center (rank_i=2), distances are {0,1,1,2,2}, max=2,
        // u = {0, 0.5, 0.5, 1, 1}, weights = {(1-0)^3=1,
        // (1-0.125)^3, (1-0.125)^3, 0, 0}.
        let w0 = 1.0_f64;
        let w1 = (1.0 - 0.125_f64).powi(3);
        let num = w0 * 3.0 + w1 * 2.0 + w1 * 4.0;
        let den = w0 + 2.0 * w1;
        let expected = num / den;
        assert!(
            (fitted[2] - expected).abs() < 1e-12,
            "center fitted {} != expected {}",
            fitted[2],
            expected
        );
    }

    /// High-count proteins should get a lower `s2_trend` than
    /// low-count proteins when the input variances track the
    /// peptide-count trend that DEqMS is designed to capture.
    #[test]
    fn deqms_assigns_lower_trend_to_higher_peptide_count() {
        let counts: Vec<u32> = (1..=10).collect();
        // Variance decreasing with peptide count, with small noise.
        let s2: Vec<f64> = counts
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let base = 2.0 / (*c as f64);
                let jitter = if i % 2 == 0 { 0.05 } else { -0.03 };
                base + jitter
            })
            .collect();
        let out = deqms_shrink(&s2, &counts, 8.0, 0.5);
        // Trend at count=1 should exceed trend at count=10.
        assert!(
            out.s2_trend[0] > out.s2_trend[9],
            "expected decreasing trend: s2_trend[0]={}, s2_trend[9]={}",
            out.s2_trend[0],
            out.s2_trend[9]
        );
        // Posterior should shrink noisy low-count features toward the
        // trend: the peak sample value (index 0, highest variance)
        // must be pulled down by the prior.
        assert!(out.s2_posterior[0] <= s2[0] + 1e-12);
    }

    /// NaN residual variances (skipped features) pass through as NaN
    /// without breaking the fit on the remaining features.
    #[test]
    fn deqms_propagates_nan_skips_without_corrupting_fit() {
        let mut s2: Vec<f64> = (1..=10).map(|c| 1.0 / c as f64).collect();
        s2[3] = f64::NAN;
        let counts: Vec<u32> = (1..=10).collect();
        let out = deqms_shrink(&s2, &counts, 6.0, 0.5);
        assert!(out.s2_posterior[3].is_nan());
        for (i, v) in out.s2_posterior.iter().enumerate() {
            if i == 3 {
                continue;
            }
            assert!(v.is_finite() && *v > 0.0, "post[{i}] = {v}");
        }
    }

    /// Fewer than two usable features → no shrinkage, `df_prior = +Inf`.
    #[test]
    fn deqms_refuses_below_two_usable_features() {
        let s2 = vec![1.0];
        let counts = vec![2u32];
        let out = deqms_shrink(&s2, &counts, 4.0, 0.5);
        assert!(out.trend_fallback_used);
        assert!(out.df_prior.is_infinite());
    }

    /// Determinism: identical inputs always produce bit-equal output
    /// arrays on the finite channels. `assert_eq` can't compare the
    /// NaN-valued `s2_prior_ratio` when the smoother degenerates,
    /// so we assert the per-feature vectors and the flags directly.
    #[test]
    fn deqms_is_deterministic() {
        let counts: Vec<u32> = (1..=12).collect();
        let s2: Vec<f64> = counts.iter().map(|c| 0.5 + 1.0 / *c as f64).collect();
        let a = deqms_shrink(&s2, &counts, 7.0, 0.5);
        let b = deqms_shrink(&s2, &counts, 7.0, 0.5);
        assert_eq!(a.s2_trend, b.s2_trend);
        assert_eq!(a.s2_posterior, b.s2_posterior);
        assert_eq!(
            a.trend_fallback_used, b.trend_fallback_used,
            "fallback flags should match"
        );
        // df_prior is a finite float or +Inf; bit-equality via `==`
        // works for both.
        assert_eq!(a.df_prior, b.df_prior);
    }
}
