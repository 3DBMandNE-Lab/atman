//! Archetype variance decomposition for `atman decompose variance`.
//!
//! For each archetype activation vector (length n_samples), fit a
//! linear mixed model (random intercept for one grouping factor +
//! optional fixed effects) and partition the observed variance into
//! per-fixed-factor contributions, random-intercept variance, and
//! residual variance.
//!
//! Reuses [`crate::de::mixed_random_intercept`] (or [`crate::de::ols`]
//! when no grouping factor is supplied) as the numerical engine; this
//! module adds only the partitioning math on top.
//!
//! **Type III sums of squares.** For each fixed factor `f` whose
//! columns in the design matrix span `c..d`, the variance
//! attributable to `f` after controlling for every other term is
//!
//! ```text
//! SS_III(f) = βᵀ_f · [(XᵀX)⁻¹_{ff}]⁻¹ · β_f
//! F_f       = SS_III(f) / (k · σ²_residual)
//! ```
//!
//! where `k = d − c` is the factor's degrees of freedom. This is
//! the Wald form of the classical Type III F-test: the variance
//! reduction you'd observe by dropping the factor from the full
//! model, partialling out every other covariate. It matches R's
//! `car::Anova(lm, type = 3)` on balanced OLS designs.
//!
//! Random-intercept variance is `σ²_u = variance_ratio · σ²_res`.
//! Residual variance is `σ²_res`. The intraclass correlation for
//! the random factor is `ICC = σ²_u / (σ²_u + σ²_res)`.
//!
//! In the mixed-effects case (`group_labels = Some(...)`), the
//! Type III F above uses the *OLS* Hessian `XᵀX` rather than the
//! GLS Hessian `Xᵀ V⁻¹ X`. Proper mixed-model Type III requires
//! the full GLS normal-equations matrix, which the current
//! `mixed_random_intercept` entry point does not expose; for mixed
//! fits the F and p-value are reported as NaN and the per-coefficient
//! max-|t| / min-p summaries remain the only factor-level inference.
//! Users who need mixed Type III should run `atman de --test mixed`
//! and apply the contrast machinery there.

use crate::de::{mixed_random_intercept, ols, omnibus_f_test, OlsOutcome};

#[derive(Debug, Clone, PartialEq)]
pub struct FixedFactor {
    /// Label (e.g. `"cohort"`, `"condition"`) used as an output
    /// column prefix.
    pub name: String,
    /// Column range in the design matrix (half-open: `[start, end)`).
    pub columns: std::ops::Range<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct VarianceRow {
    pub archetype_id: String,
    pub total_var: f64,
    pub var_residual: f64,
    pub var_random: Option<f64>,
    pub icc_random: Option<f64>,
    /// One entry per fixed factor in the input order. `var` is the
    /// Type-I projection variance; `max_abs_t` and `min_p` are the
    /// per-coefficient Wald stats across the factor's columns, kept
    /// so callers have a quick per-factor summary without the full
    /// per-coefficient table.
    pub per_factor: Vec<FactorRow>,
    /// `true` when the underlying fit skipped this archetype (e.g.
    /// singular design, too few samples, non-finite input). The
    /// other fields are still populated with `NaN` / `0.0` sentinels
    /// when skipped, and [`VarianceRow::skip_reason`] is `Some`.
    pub skipped: bool,
    pub skip_reason: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FactorRow {
    pub name: String,
    /// Type III sum of squares (`βᵀ_f · [(XᵀX)⁻¹_{ff}]⁻¹ · β_f`).
    /// `NaN` for mixed fits (see module docstring).
    pub ss_type3: f64,
    /// Factor-level omnibus F-statistic under the Type III test.
    /// `NaN` for mixed fits.
    pub f_statistic: f64,
    pub df_num: usize,
    pub df_den: f64,
    /// Two-sided F-distribution p-value of the Type III omnibus.
    pub p_value: f64,
    /// Max `|t|` across the factor's per-coefficient Wald stats.
    /// Kept for quick per-coefficient inspection.
    pub max_abs_t: f64,
    pub min_p: f64,
    pub n_coefficients: usize,
}

/// Decompose the variance of every row of `activations` (one row per
/// archetype) under the linear mixed model `y ~ X β + u_g + ε`.
///
/// `design` rows align with the columns of `activations`. If
/// `group_labels` is `Some`, fit a random-intercept mixed model on
/// the supplied group; otherwise fit plain OLS and the returned
/// `var_random` / `icc_random` are `None`.
pub fn decompose_archetype_variance(
    archetype_ids: &[String],
    activations: &[Vec<f64>],
    design: &[Vec<f64>],
    fixed_factors: &[FixedFactor],
    group_labels: Option<&[String]>,
    min_samples: usize,
) -> Result<Vec<VarianceRow>, String> {
    if archetype_ids.len() != activations.len() {
        return Err(format!(
            "archetype_ids ({}) vs activations rows ({}) mismatch",
            archetype_ids.len(),
            activations.len()
        ));
    }
    if activations.is_empty() {
        return Ok(Vec::new());
    }
    let n_samples = design.len();
    if n_samples == 0 {
        return Err("empty design matrix".into());
    }
    if activations.iter().any(|row| row.len() != n_samples) {
        return Err("activation row length must match design sample count".into());
    }
    if let Some(g) = group_labels {
        if g.len() != n_samples {
            return Err(format!(
                "group_labels ({}) vs design rows ({}) mismatch",
                g.len(),
                n_samples
            ));
        }
    }
    for f in fixed_factors {
        if f.columns.start >= f.columns.end || f.columns.end > design[0].len() {
            return Err(format!(
                "fixed factor {:?} columns {}..{} out of range (p={})",
                f.name,
                f.columns.start,
                f.columns.end,
                design[0].len()
            ));
        }
    }

    let mut out = Vec::with_capacity(activations.len());
    for (i, y) in activations.iter().enumerate() {
        let archetype_id = archetype_ids[i].clone();
        let total_var = sample_variance(y);
        let outcome = match group_labels {
            Some(g) => mixed_random_intercept(design, y, g, min_samples),
            None => ols(design, y, min_samples),
        };
        let fit = match outcome {
            OlsOutcome::Computed(f) => f,
            OlsOutcome::Skipped { reason, .. } => {
                out.push(VarianceRow {
                    archetype_id,
                    total_var,
                    var_residual: f64::NAN,
                    var_random: None,
                    icc_random: None,
                    per_factor: fixed_factors
                        .iter()
                        .map(|f| FactorRow {
                            name: f.name.clone(),
                            ss_type3: f64::NAN,
                            f_statistic: f64::NAN,
                            df_num: f.columns.end - f.columns.start,
                            df_den: f64::NAN,
                            p_value: f64::NAN,
                            max_abs_t: f64::NAN,
                            min_p: f64::NAN,
                            n_coefficients: f.columns.end - f.columns.start,
                        })
                        .collect(),
                    skipped: true,
                    skip_reason: Some(format!("{:?}", reason)),
                });
                continue;
            }
        };
        let var_residual = fit.sigma2;
        let var_random = fit.variance_ratio.map(|tau| tau * fit.sigma2);
        let icc_random = var_random.map(|vr| {
            let total = vr + var_residual;
            if total > 0.0 {
                vr / total
            } else {
                f64::NAN
            }
        });
        let is_mixed = group_labels.is_some();
        let per_factor: Vec<FactorRow> = fixed_factors
            .iter()
            .map(|f| {
                let (max_abs_t, min_p) =
                    summarize_factor_coefficients(&fit.t, &fit.p_value, f.columns.clone());
                let cols: Vec<usize> = f.columns.clone().collect();
                // Type III Wald F via the OLS Hessian. For mixed
                // fits this is NOT correct (would need the GLS
                // Hessian), so we emit NaN and keep the per-coef
                // summary as the only factor-level signal.
                let (ss_type3, f_stat, df_num, df_den, p_value) = if is_mixed {
                    (f64::NAN, f64::NAN, cols.len(), fit.df, f64::NAN)
                } else {
                    match omnibus_f_test(design, &fit.beta, &cols, fit.sigma2, fit.df) {
                        Some(om) => {
                            let ss = om.f_statistic * (om.df_num as f64) * fit.sigma2;
                            (ss, om.f_statistic, om.df_num, om.df_den, om.p_value)
                        }
                        None => (f64::NAN, f64::NAN, cols.len(), fit.df, f64::NAN),
                    }
                };
                FactorRow {
                    name: f.name.clone(),
                    ss_type3,
                    f_statistic: f_stat,
                    df_num,
                    df_den,
                    p_value,
                    max_abs_t,
                    min_p,
                    n_coefficients: f.columns.end - f.columns.start,
                }
            })
            .collect();
        out.push(VarianceRow {
            archetype_id,
            total_var,
            var_residual,
            var_random,
            icc_random,
            per_factor,
            skipped: false,
            skip_reason: None,
        });
    }
    Ok(out)
}

fn sample_variance(values: &[f64]) -> f64 {
    let n = values.len();
    if n < 2 {
        return 0.0;
    }
    let mean: f64 = values.iter().sum::<f64>() / n as f64;
    values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64
}

fn summarize_factor_coefficients(
    t: &[f64],
    p_value: &[f64],
    columns: std::ops::Range<usize>,
) -> (f64, f64) {
    let mut max_abs_t = f64::NEG_INFINITY;
    let mut min_p = f64::INFINITY;
    for c in columns {
        let tc = t[c];
        let pc = p_value[c];
        if tc.is_finite() && tc.abs() > max_abs_t {
            max_abs_t = tc.abs();
        }
        if pc.is_finite() && pc < min_p {
            min_p = pc;
        }
    }
    let max_abs_t = if max_abs_t.is_finite() {
        max_abs_t
    } else {
        f64::NAN
    };
    let min_p = if min_p.is_finite() { min_p } else { f64::NAN };
    (max_abs_t, min_p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_variance_is_n_minus_one_denominator() {
        // var of [1,2,3,4,5] with ddof=1: ((1-3)^2+(2-3)^2+...+(5-3)^2)/4 = 10/4 = 2.5
        let v = sample_variance(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        assert!((v - 2.5).abs() < 1e-12, "got {v}");
    }

    #[test]
    fn type3_ss_matches_wald_form_on_single_column_factor() {
        // For a one-column factor, the Type III F equals (β / se)²
        // and SS_III = F · σ² — both derivable from the OLS fit
        // the test itself builds. Ensures the omnibus_f_test path
        // agrees with the t² identity.
        let design: Vec<Vec<f64>> = (0..10)
            .map(|i| vec![1.0, if i < 5 { 0.0 } else { 1.0 }])
            .collect();
        let y: Vec<f64> = (0..10)
            .map(|i| if i < 5 { 0.0 } else { 2.5 } + (i as f64 * 0.1).sin() * 0.2)
            .collect();
        let factors = vec![FixedFactor {
            name: "condition".into(),
            columns: 1..2,
        }];
        let ids = vec!["archetype_01".into()];
        let rows = decompose_archetype_variance(&ids, &[y], &design, &factors, None, 2).unwrap();
        let r = &rows[0].per_factor[0];
        assert_eq!(r.df_num, 1);
        // F = t²: the omnibus on a single column must reproduce the
        // squared Wald t-statistic to machine precision.
        let implied_t2 = r.max_abs_t.powi(2);
        assert!(
            (r.f_statistic - implied_t2).abs() < 1e-8,
            "F {} vs t² {}",
            r.f_statistic,
            implied_t2
        );
        // Single-column F-p equals two-sided t-p.
        assert!((r.p_value - r.min_p).abs() < 1e-8);
        // SS_III = F · σ²
        let expected_ss = r.f_statistic * rows[0].var_residual;
        assert!(
            (r.ss_type3 - expected_ss).abs() < 1e-8,
            "SS_III {} vs F·σ² {}",
            r.ss_type3,
            expected_ss
        );
    }

    #[test]
    fn decompose_recovers_zero_random_variance_on_pure_fixed_effect_data() {
        // n=10, condition column = [0,0,0,0,0,1,1,1,1,1]. activations
        // = 1.0 * condition + tiny noise. No group structure.
        let design: Vec<Vec<f64>> = (0..10)
            .map(|i| vec![1.0, if i < 5 { 0.0 } else { 1.0 }])
            .collect();
        let y: Vec<f64> = (0..10)
            .map(|i| if i < 5 { 0.0 } else { 1.0 } + (i as f64) * 1e-5)
            .collect();
        let factors = vec![
            FixedFactor {
                name: "intercept".into(),
                columns: 0..1,
            },
            FixedFactor {
                name: "condition".into(),
                columns: 1..2,
            },
        ];
        let ids = vec!["archetype_01".into()];
        let rows = decompose_archetype_variance(&ids, &[y], &design, &factors, None, 2).unwrap();
        assert_eq!(rows.len(), 1);
        let r = &rows[0];
        assert!(!r.skipped);
        assert!(r.var_random.is_none());
        assert!(r.icc_random.is_none());
        let cond_factor = r.per_factor.iter().find(|f| f.name == "condition").unwrap();
        // On 10 near-noise-free samples of y = 1[i≥5] + tiny, the
        // condition coefficient is near 1 so Type III F is huge
        // and the Type III SS dominates the residual variance.
        assert!(
            cond_factor.f_statistic > 100.0,
            "condition F should be large on near-clean signal; got {}",
            cond_factor.f_statistic
        );
        assert!(
            cond_factor.p_value < 1e-5,
            "condition p should clear 1e-5; got {}",
            cond_factor.p_value
        );
        assert!(cond_factor.ss_type3 > 0.1);
        assert!(
            r.var_residual < 1e-5,
            "residual should be near zero on clean data: {}",
            r.var_residual
        );
    }

    #[test]
    fn decompose_detects_nonzero_random_variance_with_grouping() {
        // Two groups of 6 each, activations = group_intercept + noise.
        let mut design = Vec::new();
        let mut y = Vec::new();
        let mut groups = Vec::new();
        let group_a_offset = 1.0;
        let group_b_offset = -1.0;
        for i in 0..12 {
            let in_a = i < 6;
            design.push(vec![1.0]);
            y.push(if in_a { group_a_offset } else { group_b_offset } + (i as f64) * 1e-4);
            groups.push(if in_a { "A" } else { "B" }.to_string());
        }
        let factors = vec![FixedFactor {
            name: "intercept".into(),
            columns: 0..1,
        }];
        let ids = vec!["archetype_01".into()];
        let rows =
            decompose_archetype_variance(&ids, &[y.clone()], &design, &factors, Some(&groups), 2)
                .unwrap();
        assert_eq!(rows.len(), 1);
        let r = &rows[0];
        assert!(!r.skipped, "skip reason: {:?}", r.skip_reason);
        let var_r = r.var_random.expect("mixed fit should populate var_random");
        let icc = r.icc_random.expect("icc");
        assert!(var_r > 0.1, "random var should be substantial: {var_r}");
        assert!(
            icc > 0.5,
            "ICC should exceed 0.5 on group-dominated data: {icc}"
        );
    }

    #[test]
    fn decompose_propagates_skip_on_singular_design() {
        // p > n: design has 3 columns, only 2 samples.
        let design = vec![vec![1.0, 0.0, 0.0], vec![1.0, 1.0, 1.0]];
        let y = vec![vec![1.0, 2.0]];
        let factors = vec![FixedFactor {
            name: "x".into(),
            columns: 0..3,
        }];
        let ids = vec!["archetype_01".into()];
        let rows = decompose_archetype_variance(&ids, &y, &design, &factors, None, 2).unwrap();
        assert!(rows[0].skipped);
        assert!(rows[0].skip_reason.is_some());
    }

    #[test]
    fn decompose_is_deterministic() {
        let design: Vec<Vec<f64>> = (0..8)
            .map(|i| vec![1.0, if i < 4 { 0.0 } else { 1.0 }])
            .collect();
        let y: Vec<Vec<f64>> = (0..3)
            .map(|k| (0..8).map(|i| (i + k) as f64 * 0.3).collect())
            .collect();
        let factors = vec![
            FixedFactor {
                name: "intercept".into(),
                columns: 0..1,
            },
            FixedFactor {
                name: "condition".into(),
                columns: 1..2,
            },
        ];
        let ids: Vec<String> = (1..=3).map(|i| format!("archetype_{i:02}")).collect();
        let a = decompose_archetype_variance(&ids, &y, &design, &factors, None, 2).unwrap();
        let b = decompose_archetype_variance(&ids, &y, &design, &factors, None, 2).unwrap();
        assert_eq!(a, b);
    }
}
