//! Cross-method ensemble DE aggregation. Pure functions; no IO.
//!
//! Input: per-method per-protein `(mean_diff, p_value, bh_q)` triples
//! keyed by `(comparison_label, assay_id)`. Output: per-protein
//! ensemble summary with Stouffer-combined p, majority sign, and a
//! VALIDATED / PROVISIONAL / INSUFFICIENT grade.
//!
//! Conventions:
//! - `mean_diff` uses atman's `mean_a − mean_b` convention for every
//!   method; the aggregation trusts that contract.
//! - `None` p-values / `None` mean-diff are treated as
//!   method-not-applicable for that protein and excluded from counts.
//! - Significance uses per-method BH-q against a caller-supplied
//!   threshold (default 0.05).

use statrs::distribution::{ContinuousCDF, Normal};

#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleInput {
    pub method: String,
    pub mean_diff: Option<f64>,
    pub p_value: Option<f64>,
    pub bh_q: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnsembleGrade {
    Validated,
    Provisional,
    Insufficient,
}

impl EnsembleGrade {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Validated => "VALIDATED",
            Self::Provisional => "PROVISIONAL",
            Self::Insufficient => "INSUFFICIENT",
        }
    }
}

/// Thresholds for [`assign_grade`] and [`aggregate_per_protein`].
///
/// `q_threshold` is used two ways: once inside
/// [`aggregate_per_protein`] to count per-method significant
/// findings (informational — `n_significant`), and once inside
/// [`assign_grade`] as the `ensemble_q` bar for VALIDATED /
/// PROVISIONAL assignment.
#[derive(Debug, Clone, Copy)]
pub struct GradeThresholds {
    /// Per-method BH-q threshold (counts `n_significant`) AND the
    /// ensemble-BH q-threshold for VALIDATED / PROVISIONAL.
    pub q_threshold: f64,
    /// Sign-consistency fraction required for VALIDATED.
    pub validated_sign_fraction: f64,
    /// Sign-consistency fraction required for PROVISIONAL. Below
    /// this, INSUFFICIENT.
    pub provisional_sign_fraction: f64,
}

impl Default for GradeThresholds {
    fn default() -> Self {
        Self {
            q_threshold: 0.05,
            validated_sign_fraction: 1.00,
            provisional_sign_fraction: 0.50,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleRow {
    pub n_applied: usize,
    pub n_significant: usize,
    pub n_sign_consistent: usize,
    /// +1.0, -1.0, or 0.0 when no applicable rows.
    pub majority_sign: f64,
    pub ensemble_p: Option<f64>,
    pub methods_applied: Vec<String>,
}

/// Assign a VALIDATED / PROVISIONAL / INSUFFICIENT grade to one
/// protein given its ensemble BH-q (computed across proteins within
/// a comparison) and its sign-consistency fraction.
///
/// Grade criteria:
/// - **INSUFFICIENT** when no method fit the protein
///   (`n_applied == 0`) or when `ensemble_q` is missing or
///   `>= thresholds.q_threshold`.
/// - **VALIDATED** when `ensemble_q < q_threshold` AND the
///   sign-consistency fraction reaches
///   `validated_sign_fraction` (default 1.00).
/// - **PROVISIONAL** when `ensemble_q < q_threshold` AND the
///   sign-consistency fraction reaches
///   `provisional_sign_fraction` (default 0.50).
/// - Otherwise INSUFFICIENT.
pub fn assign_grade(
    ensemble_q: Option<f64>,
    n_applied: usize,
    n_sign_consistent: usize,
    thresholds: GradeThresholds,
) -> EnsembleGrade {
    if n_applied == 0 {
        return EnsembleGrade::Insufficient;
    }
    let Some(q) = ensemble_q else {
        return EnsembleGrade::Insufficient;
    };
    if !q.is_finite() || q >= thresholds.q_threshold {
        return EnsembleGrade::Insufficient;
    }
    let sign_frac = n_sign_consistent as f64 / n_applied as f64;
    if sign_frac >= thresholds.validated_sign_fraction {
        EnsembleGrade::Validated
    } else if sign_frac >= thresholds.provisional_sign_fraction {
        EnsembleGrade::Provisional
    } else {
        EnsembleGrade::Insufficient
    }
}

/// Stouffer's Z-score combination. Returns `None` when no finite
/// p-values are present. Equal weights; exact-zero p clamped to 1e-300
/// and exact-one p clamped to `1 - 1e-16` to keep `Φ⁻¹` finite.
pub fn combine_stouffer(p_values: &[f64]) -> Option<f64> {
    let normal = Normal::new(0.0, 1.0).expect("N(0,1) construction is infallible");
    let mut z_sum = 0.0;
    let mut k = 0usize;
    for &p in p_values {
        if !p.is_finite() {
            continue;
        }
        let p_clamped = p.clamp(1e-300, 1.0 - 1e-16);
        z_sum += normal.inverse_cdf(1.0 - p_clamped);
        k += 1;
    }
    if k == 0 {
        return None;
    }
    let z_combined = z_sum / (k as f64).sqrt();
    Some(1.0 - normal.cdf(z_combined))
}

/// Aggregate per-method rows for one (comparison, protein) into an
/// [`EnsembleRow`]. Does NOT assign a grade — the caller applies
/// [`assign_grade`] after BH-FDR on the ensemble p-values within the
/// comparison has produced an `ensemble_q`.
///
/// `n_significant` is the count of methods whose per-method BH-q was
/// below `thresholds.q_threshold`. It is informational only, not a
/// grade determinant in the current design.
pub fn aggregate_per_protein(inputs: &[EnsembleInput], thresholds: GradeThresholds) -> EnsembleRow {
    let mut applied_p: Vec<f64> = Vec::new();
    let mut applied_signs: Vec<f64> = Vec::new();
    let mut applied_methods: Vec<String> = Vec::new();
    let mut n_significant = 0usize;
    for row in inputs {
        let (Some(p), Some(m)) = (row.p_value, row.mean_diff) else {
            continue;
        };
        if !p.is_finite() || !m.is_finite() {
            continue;
        }
        applied_p.push(p);
        applied_signs.push(m.signum());
        applied_methods.push(row.method.clone());
        if let Some(q) = row.bh_q {
            if q.is_finite() && q < thresholds.q_threshold {
                n_significant += 1;
            }
        }
    }
    let n_applied = applied_p.len();
    if n_applied == 0 {
        return EnsembleRow {
            n_applied: 0,
            n_significant: 0,
            n_sign_consistent: 0,
            majority_sign: 0.0,
            ensemble_p: None,
            methods_applied: Vec::new(),
        };
    }
    let pos = applied_signs.iter().filter(|s| **s > 0.0).count();
    let neg = applied_signs.iter().filter(|s| **s < 0.0).count();
    let (majority_sign, n_sign_consistent) = if pos >= neg { (1.0, pos) } else { (-1.0, neg) };
    let ensemble_p = combine_stouffer(&applied_p);
    EnsembleRow {
        n_applied,
        n_significant,
        n_sign_consistent,
        majority_sign,
        ensemble_p,
        methods_applied: applied_methods,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk(method: &str, mean: f64, p: f64, q: f64) -> EnsembleInput {
        EnsembleInput {
            method: method.into(),
            mean_diff: Some(mean),
            p_value: Some(p),
            bh_q: Some(q),
        }
    }

    #[test]
    fn stouffer_two_equal_p_values_exact() {
        // Stouffer of (p, p) = Φ(−√2 · Φ⁻¹(p))
        let p = 0.1_f64;
        let combined = combine_stouffer(&[p, p]).unwrap();
        let normal = Normal::new(0.0, 1.0).unwrap();
        let expected = 1.0 - normal.cdf((2.0_f64).sqrt() * normal.inverse_cdf(1.0 - p));
        assert!(
            (combined - expected).abs() < 1e-9,
            "got {combined}, want {expected}"
        );
    }

    #[test]
    fn stouffer_returns_none_on_empty_and_all_nan() {
        assert!(combine_stouffer(&[]).is_none());
        assert!(combine_stouffer(&[f64::NAN, f64::NAN]).is_none());
    }

    #[test]
    fn aggregate_collects_counts_and_majority_sign() {
        let rows = vec![
            mk("welch-t", 1.5, 1e-5, 1e-4),
            mk("ols", 1.4, 5e-6, 8e-5),
            mk("limma", 1.6, 1e-6, 3e-5),
            mk("msqrob", 1.3, 2e-5, 5e-4),
        ];
        let out = aggregate_per_protein(&rows, GradeThresholds::default());
        assert_eq!(out.n_applied, 4);
        assert_eq!(out.n_significant, 4);
        assert_eq!(out.n_sign_consistent, 4);
        assert!(out.majority_sign > 0.0);
        // Grade is not computed by aggregate_per_protein anymore; caller
        // combines this with ensemble_q via `assign_grade`.
    }

    #[test]
    fn aggregate_handles_mixed_signs() {
        let rows = vec![
            mk("welch-t", 1.5, 1e-5, 1e-4),
            mk("limma", -1.4, 1e-5, 1e-4),
            mk("msqrob", 1.2, 1e-5, 1e-4),
        ];
        let out = aggregate_per_protein(&rows, GradeThresholds::default());
        // Majority is positive (2 vs 1), sign_consistent count = 2.
        assert_eq!(out.n_sign_consistent, 2);
        assert_eq!(out.majority_sign, 1.0);
    }

    #[test]
    fn skipped_methods_are_excluded_from_counts() {
        let rows = vec![
            mk("welch-t", 1.5, 1e-5, 1e-4),
            EnsembleInput {
                method: "msqrob".into(),
                mean_diff: None,
                p_value: None,
                bh_q: None,
            },
        ];
        let out = aggregate_per_protein(&rows, GradeThresholds::default());
        assert_eq!(out.n_applied, 1);
        assert_eq!(out.methods_applied, vec!["welch-t"]);
    }

    #[test]
    fn empty_input_yields_zero_counts_and_none_ensemble_p() {
        let out = aggregate_per_protein(&[], GradeThresholds::default());
        assert_eq!(out.n_applied, 0);
        assert_eq!(out.majority_sign, 0.0);
        assert!(out.ensemble_p.is_none());
    }

    #[test]
    fn assign_grade_validated_requires_significant_ensemble_q_and_full_sign() {
        let t = GradeThresholds::default();
        assert_eq!(assign_grade(Some(1e-4), 4, 4, t), EnsembleGrade::Validated);
    }

    #[test]
    fn assign_grade_provisional_when_sign_partial() {
        let t = GradeThresholds::default();
        // 2/3 sign-consistent → below validated (1.0) threshold but
        // at/above provisional (0.5). ensemble_q significant.
        assert_eq!(
            assign_grade(Some(1e-4), 3, 2, t),
            EnsembleGrade::Provisional
        );
    }

    #[test]
    fn assign_grade_insufficient_when_ensemble_q_fails() {
        let t = GradeThresholds::default();
        // ensemble_q above threshold regardless of sign consistency.
        assert_eq!(
            assign_grade(Some(0.5), 4, 4, t),
            EnsembleGrade::Insufficient
        );
        // Missing ensemble_q.
        assert_eq!(assign_grade(None, 4, 4, t), EnsembleGrade::Insufficient);
        // No applicable methods.
        assert_eq!(
            assign_grade(Some(1e-4), 0, 0, t),
            EnsembleGrade::Insufficient
        );
    }

    #[test]
    fn assign_grade_insufficient_when_majority_sign_fails_provisional_threshold() {
        let t = GradeThresholds::default();
        // 2/5 sign-consistent → below provisional (0.5).
        assert_eq!(
            assign_grade(Some(1e-4), 5, 2, t),
            EnsembleGrade::Insufficient
        );
    }
}
