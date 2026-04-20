//! Recovery scoring for `atman bench decompose`.
//!
//! Given a ground-truth loading matrix `planted[k × p]` and a
//! recovered loading matrix `recovered[k' × p]` (not necessarily the
//! same `k'`, and always sign- and permutation-ambiguous under ICA),
//! match the rows pairwise via reciprocal-best absolute cosine and
//! compute per-archetype recovery metrics:
//!
//! - `archetype_correlation`: `|pearson(planted_row, matched_row)|`.
//! - `recovery_jaccard`: Jaccard overlap between the two rows'
//!   top-N by `|loading|` sets. Default `N = 50`, caller-tunable.
//!
//! Unmatched planted archetypes (no reciprocal-best partner) return
//! zeros on every metric. This is the world's-smallest benchmark
//! harness — scope is deliberately narrow so the whole surface
//! lives in ≤ 300 LoC of dependency-free math.

use crate::stats;

#[derive(Debug, Clone)]
pub struct RecoveryScore {
    pub planted_archetype: usize,
    pub matched_recovered: Option<usize>,
    pub archetype_correlation: f64,
    pub recovery_jaccard: f64,
    pub matched_cosine: f64,
}

/// Score every row of `planted` against its reciprocal-best match in
/// `recovered`. Both matrices are row-major: `values[k][p]`.
pub fn score_against_planted(
    planted: &[Vec<f64>],
    recovered: &[Vec<f64>],
    top_n: usize,
) -> Vec<RecoveryScore> {
    let k_planted = planted.len();
    let k_recovered = recovered.len();
    if k_planted == 0 {
        return Vec::new();
    }
    if k_recovered == 0 {
        return (0..k_planted)
            .map(|i| RecoveryScore {
                planted_archetype: i,
                matched_recovered: None,
                archetype_correlation: 0.0,
                recovery_jaccard: 0.0,
                matched_cosine: 0.0,
            })
            .collect();
    }
    // Absolute cosine between every (planted, recovered) pair.
    let mut sim = vec![vec![0.0_f64; k_recovered]; k_planted];
    for i in 0..k_planted {
        for j in 0..k_recovered {
            let c = stats::cosine(&planted[i], &recovered[j]).unwrap_or(0.0).abs();
            sim[i][j] = c;
        }
    }
    // Reciprocal-best: for each planted i, best j maximizes sim[i][.];
    // valid if that j's own argmax over i also returns i.
    let best_for_planted: Vec<Option<usize>> =
        (0..k_planted).map(|i| argmax(&sim[i])).collect();
    let best_for_recovered: Vec<Option<usize>> = (0..k_recovered)
        .map(|j| argmax_col(&sim, j, k_planted))
        .collect();

    let mut out = Vec::with_capacity(k_planted);
    for i in 0..k_planted {
        let matched = match best_for_planted[i] {
            Some(j) if best_for_recovered.get(j).copied().flatten() == Some(i) => Some(j),
            _ => None,
        };
        match matched {
            Some(j) => {
                let cos = sim[i][j];
                let corr = stats::pearson(&planted[i], &recovered[j])
                    .map(|r| r.abs())
                    .unwrap_or(0.0);
                let jacc = top_n_abs_jaccard(&planted[i], &recovered[j], top_n);
                out.push(RecoveryScore {
                    planted_archetype: i,
                    matched_recovered: Some(j),
                    archetype_correlation: corr,
                    recovery_jaccard: jacc,
                    matched_cosine: cos,
                });
            }
            None => out.push(RecoveryScore {
                planted_archetype: i,
                matched_recovered: None,
                archetype_correlation: 0.0,
                recovery_jaccard: 0.0,
                matched_cosine: 0.0,
            }),
        }
    }
    out
}

fn argmax(values: &[f64]) -> Option<usize> {
    values
        .iter()
        .enumerate()
        .fold(None, |acc, (i, &v)| match acc {
            None => Some((i, v)),
            Some((_, best)) if v > best => Some((i, v)),
            other => other,
        })
        .map(|(i, _)| i)
}

fn argmax_col(sim: &[Vec<f64>], col: usize, n_rows: usize) -> Option<usize> {
    sim.iter()
        .take(n_rows)
        .enumerate()
        .fold((None, f64::NEG_INFINITY), |(best, best_val), (i, row)| {
            if row[col] > best_val {
                (Some(i), row[col])
            } else {
                (best, best_val)
            }
        })
        .0
}

fn top_n_abs_jaccard(a: &[f64], b: &[f64], top_n: usize) -> f64 {
    if top_n == 0 || a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let set_a = top_n_abs_set(a, top_n);
    let set_b = top_n_abs_set(b, top_n);
    let inter = set_a.intersection(&set_b).count();
    let union = set_a.union(&set_b).count();
    if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    }
}

fn top_n_abs_set(values: &[f64], top_n: usize) -> std::collections::BTreeSet<usize> {
    use std::cmp::Ordering;
    let mut ranked: Vec<(usize, f64)> = values
        .iter()
        .enumerate()
        .map(|(i, v)| (i, v.abs()))
        .collect();
    ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
    ranked
        .into_iter()
        .take(top_n)
        .map(|(i, _)| i)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_recovery_scores_all_ones() {
        let planted = vec![vec![1.0, 0.9, 0.1, 0.0], vec![0.0, 0.1, 0.9, 1.0]];
        let recovered = planted.clone();
        let scores = score_against_planted(&planted, &recovered, 2);
        assert_eq!(scores.len(), 2);
        for s in &scores {
            assert!((s.archetype_correlation - 1.0).abs() < 1e-9);
            assert!((s.recovery_jaccard - 1.0).abs() < 1e-9);
            assert!(s.matched_recovered.is_some());
        }
    }

    #[test]
    fn sign_flipped_recovery_scores_one_under_abs_cosine() {
        let planted = vec![vec![1.0, 0.9, 0.1, 0.0]];
        let recovered = vec![vec![-1.0, -0.9, -0.1, 0.0]];
        let s = &score_against_planted(&planted, &recovered, 2)[0];
        assert!((s.archetype_correlation - 1.0).abs() < 1e-9);
        assert!((s.recovery_jaccard - 1.0).abs() < 1e-9);
    }

    #[test]
    fn unrelated_recovery_scores_near_zero_jaccard() {
        let planted = vec![vec![1.0, 0.9, 0.0, 0.0]];
        let recovered = vec![vec![0.0, 0.0, 1.0, 0.9]];
        let s = &score_against_planted(&planted, &recovered, 2)[0];
        // They reciprocal-best match (only one candidate each), but
        // the top-2 sets are disjoint so jaccard = 0.
        assert!(s.matched_recovered.is_some());
        assert!(s.recovery_jaccard < 1e-9);
    }

    #[test]
    fn empty_recovered_yields_zero_scores() {
        let planted = vec![vec![1.0, 0.9, 0.1, 0.0]];
        let recovered: Vec<Vec<f64>> = Vec::new();
        let scores = score_against_planted(&planted, &recovered, 2);
        assert_eq!(scores.len(), 1);
        assert!(scores[0].matched_recovered.is_none());
        assert_eq!(scores[0].archetype_correlation, 0.0);
    }

    #[test]
    fn reciprocal_best_rejects_one_to_many_match() {
        // Two planted archetypes, one recovered. Both planted prefer
        // the single recovered, but only one is the recovered's top
        // pick back — the other is unmatched.
        let planted = vec![
            vec![1.0, 0.5, 0.0], // prefers recovered[0] (cosine ~1)
            vec![0.6, 0.5, 0.0], // also prefers recovered[0] (smaller cosine)
        ];
        let recovered = vec![vec![1.0, 0.5, 0.0]];
        let scores = score_against_planted(&planted, &recovered, 2);
        let matched: Vec<_> = scores
            .iter()
            .filter(|s| s.matched_recovered.is_some())
            .collect();
        assert_eq!(matched.len(), 1, "exactly one reciprocal-best match");
    }
}
