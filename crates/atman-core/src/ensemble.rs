//! Cross-method ensemble DE aggregation. Pure functions; no IO.
//!
//! Input: per-method per-protein `(mean_diff, p_value, bh_q)` triples
//! keyed by `(comparison_label, assay_id)`. Output: per-protein
//! ensemble summary with a method-agreement / sign-consistency grade,
//! plus a Stouffer-combined p reported as a heuristic only.
//!
//! Conventions:
//! - `mean_diff` uses atman's `mean_a − mean_b` convention for every
//!   method; the aggregation trusts that contract.
//! - `None` p-values / `None` mean-diff are treated as
//!   method-not-applicable for that protein and excluded from counts.
//! - Significance uses per-method BH-q against a caller-supplied
//!   threshold (default 0.05).
//!
//! ## Why the grade is NOT driven by a combined p-value
//!
//! Every ensemble method (`paired-t`, `welch-t`, `ols`, `mixed`,
//! `limma`, `msqrob`) is fit on the SAME abundance matrix, so their
//! per-method p-values are strongly positively correlated. Stouffer's
//! method (see [`combine_stouffer`]) assumes the inputs are
//! INDEPENDENT; under positive correlation it is anti-conservative —
//! the combined p (and the BH-q derived from it) is smaller than the
//! evidence warrants. Feeding that into a VALIDATED grade would
//! over-state confidence.
//!
//! The grade therefore uses only quantities that do not assume
//! independence: how many methods *individually* clear per-method
//! BH-q significance (`n_significant`) and how consistent their
//! effect sign is (`n_sign_consistent` / `n_applied`). `ensemble_p`
//! and the `ensemble_q` derived from it are retained as a
//! convenience ranking heuristic ONLY — they are NOT a calibrated
//! p-value, NOT used by [`assign_grade`], and must not be read as
//! evidence of statistical significance.

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
/// `q_threshold` is the per-method BH-q bar applied inside
/// [`aggregate_per_protein`] to count how many methods individually
/// call the protein significant (`n_significant`). That count, with
/// sign consistency, drives the grade — the combined `ensemble_q` is
/// no longer a grade determinant (see module docs).
#[derive(Debug, Clone, Copy)]
pub struct GradeThresholds {
    /// Per-method BH-q threshold used to count `n_significant`.
    pub q_threshold: f64,
    /// Both the per-method significant fraction AND the sign-consistency
    /// fraction must reach this for VALIDATED.
    pub validated_sign_fraction: f64,
    /// Both the per-method significant fraction AND the sign-consistency
    /// fraction must reach this for PROVISIONAL. Below this,
    /// INSUFFICIENT.
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
/// protein from how many methods *independently* agree, NOT from a
/// combined p-value.
///
/// The grade reads two independence-free quantities:
/// - the per-method significant fraction `n_significant / n_applied`
///   — how many of the fitted methods individually clear per-method
///   BH-q significance; and
/// - the sign-consistency fraction `n_sign_consistent / n_applied`
///   — how many fitted methods share the majority effect direction.
///
/// This deliberately avoids the Stouffer-combined `ensemble_q`,
/// which is anti-conservative here because all methods share the
/// same abundance matrix (see module docs).
///
/// Grade criteria:
/// - **INSUFFICIENT** when no method fit the protein
///   (`n_applied == 0`).
/// - **VALIDATED** when BOTH the significant fraction AND the
///   sign-consistency fraction reach `validated_sign_fraction`
///   (default 1.00 — every method significant and agreeing on sign).
/// - **PROVISIONAL** when BOTH fractions reach
///   `provisional_sign_fraction` (default 0.50 — a majority of
///   methods significant and agreeing on sign).
/// - Otherwise INSUFFICIENT.
pub fn assign_grade(
    n_significant: usize,
    n_applied: usize,
    n_sign_consistent: usize,
    thresholds: GradeThresholds,
) -> EnsembleGrade {
    if n_applied == 0 {
        return EnsembleGrade::Insufficient;
    }
    let sig_frac = n_significant as f64 / n_applied as f64;
    let sign_frac = n_sign_consistent as f64 / n_applied as f64;
    if sig_frac >= thresholds.validated_sign_fraction
        && sign_frac >= thresholds.validated_sign_fraction
    {
        EnsembleGrade::Validated
    } else if sig_frac >= thresholds.provisional_sign_fraction
        && sign_frac >= thresholds.provisional_sign_fraction
    {
        EnsembleGrade::Provisional
    } else {
        EnsembleGrade::Insufficient
    }
}

/// Stouffer's Z-score combination. Returns `None` when no finite
/// p-values are present. Equal weights; exact-zero p clamped to 1e-300
/// and exact-one p clamped to `1 - 1e-16` to keep `Φ⁻¹` finite.
///
/// HEURISTIC ONLY. Stouffer's method assumes the combined p-values
/// are independent. The ensemble methods are all fit on the same
/// abundance matrix, so their p-values are positively correlated and
/// this combination is anti-conservative. The result (and the
/// `ensemble_q` BH-adjusted from it) is exposed purely as a ranking
/// convenience — it is NOT a calibrated p-value and is NOT used to
/// assign grades (see [`assign_grade`] and module docs).
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
    Some(normal.sf(z_combined))
}

/// Aggregate per-method rows for one (comparison, protein) into an
/// [`EnsembleRow`]. Does NOT assign a grade — the caller applies
/// [`assign_grade`] after BH-FDR on the ensemble p-values within the
/// comparison has produced an `ensemble_q`.
///
/// `n_significant` is the count of methods whose per-method BH-q was
/// below `thresholds.q_threshold`. It (with sign consistency) is the
/// grade determinant — see [`assign_grade`]. The Stouffer
/// `ensemble_p` produced here is a heuristic only (see
/// [`combine_stouffer`]).
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
    fn assign_grade_validated_requires_all_methods_significant_and_full_sign() {
        let t = GradeThresholds::default();
        // assign_grade(n_significant, n_applied, n_sign_consistent, t).
        // 4/4 significant AND 4/4 sign-consistent → VALIDATED.
        assert_eq!(assign_grade(4, 4, 4, t), EnsembleGrade::Validated);
    }

    #[test]
    fn assign_grade_provisional_when_significant_and_sign_are_majority() {
        let t = GradeThresholds::default();
        // 2/3 significant AND 2/3 sign-consistent → below validated
        // (1.0) but at/above provisional (0.5).
        assert_eq!(assign_grade(2, 3, 2, t), EnsembleGrade::Provisional);
    }

    #[test]
    fn assign_grade_does_not_validate_on_sign_alone() {
        let t = GradeThresholds::default();
        // All methods agree on sign but only 1/4 is significant →
        // significant fraction 0.25 < provisional 0.5 → INSUFFICIENT.
        // (Under the old combined-p design this could have validated.)
        assert_eq!(assign_grade(1, 4, 4, t), EnsembleGrade::Insufficient);
    }

    #[test]
    fn assign_grade_insufficient_when_no_methods_applied() {
        let t = GradeThresholds::default();
        assert_eq!(assign_grade(0, 0, 0, t), EnsembleGrade::Insufficient);
    }

    #[test]
    fn assign_grade_insufficient_when_sign_fails_provisional_threshold() {
        let t = GradeThresholds::default();
        // 5/5 significant but only 2/5 sign-consistent → below
        // provisional (0.5) on sign → INSUFFICIENT.
        assert_eq!(assign_grade(5, 5, 2, t), EnsembleGrade::Insufficient);
    }
}
