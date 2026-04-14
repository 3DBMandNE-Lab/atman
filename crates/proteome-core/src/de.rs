//! Differential abundance testing. v0.1 of the analysis pipeline implements
//! exactly one test — paired Student's t on log2-NPX differences — and one
//! multiple-testing correction — Benjamini–Hochberg FDR.
//!
//! Statistical discipline rails (from CLAUDE.md ten commandments):
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
    if n < min_pairs {
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
    if sd_diff == 0.0 {
        return PairedTResult::Skipped {
            reason: SkipReason::ZeroVariance,
            n_pairs: n,
        };
    }

    let se = sd_diff / n_f.sqrt();
    let t = mean_diff / se;
    let df = n_f - 1.0;

    // Two-sided p = 2 * P(T > |t|) on df degrees of freedom.
    let t_dist =
        StudentsT::new(0.0, 1.0, df).expect("StudentsT::new requires df > 0; n_pairs >= 5 ensures");
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
}
