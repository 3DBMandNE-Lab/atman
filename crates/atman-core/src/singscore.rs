//! Per-sample signature scoring (singscore; Foroutan et al. 2018).
//!
//! For each sample, proteins are ranked ascending by abundance (ties broken
//! by average rank, matching R's `rank(..., ties.method = "average")`).
//! Each gene-set signature is then scored as the centered, scale-normalized
//! mean rank of its members within the per-sample ranking:
//!
//! ```text
//!     score = (2 * mean_rank - n - 1) / (2 * (n - m))
//! ```
//!
//! where `n` is the count of proteins observed in the sample (non-missing)
//! and `m` is the count of signature genes observed in that sample. The
//! score lies in `[-0.5, 0.5]`: `+0.5` when all signature genes rank at the
//! top of the sample, `-0.5` when all rank at the bottom, `0` when the
//! signature ranks are centered. This matches `singscore::simpleScore` with
//! `centerScore = TRUE` (its default).
//!
//! Missing values are handled per-sample: a protein with missing abundance
//! in sample `s` is excluded from `s`'s ranking and from any signature it
//! belongs to in `s`. Signatures whose `m < min_set_size` in a given sample
//! are emitted with `score = NaN` so the sample row is preserved (downstream
//! consumers can decide to filter or impute).
//!
//! The algorithm is fully deterministic — no randomness, no permutations.
//! All sample ordering, signature ordering, and tie-breaking are stable.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

/// Per-sample, per-signature score row.
#[derive(Debug, Clone)]
pub struct SingscoreRow {
    pub sample_id: String,
    pub set_name: String,
    pub score: f64,
    pub set_size_declared: usize,
    pub set_size_observed: usize,
    pub n_observed_in_sample: usize,
}

/// Score every (sample, signature) pair.
///
/// `abundance` maps `sample_id -> gene_symbol -> abundance` (non-finite or
/// missing values are excluded by the caller; only finite values should
/// appear in the inner map). `gene_sets` is the standard `set_name -> gene_symbol set`
/// representation. `min_set_size` is the minimum number of signature genes
/// that must be observed in a given sample for that signature's score to
/// be computed; below that, the row is emitted with `score = NaN`.
///
/// Output is sorted by `(set_name, sample_id)` for stable downstream joins.
pub fn singscore(
    abundance: &BTreeMap<String, BTreeMap<String, f64>>,
    gene_sets: &BTreeMap<String, BTreeSet<String>>,
    min_set_size: usize,
) -> Vec<SingscoreRow> {
    let mut out = Vec::with_capacity(abundance.len() * gene_sets.len());
    for (sample_id, gene_to_value) in abundance {
        let n = gene_to_value.len();
        if n == 0 {
            for (set_name, set_genes) in gene_sets {
                out.push(SingscoreRow {
                    sample_id: sample_id.clone(),
                    set_name: set_name.clone(),
                    score: f64::NAN,
                    set_size_declared: set_genes.len(),
                    set_size_observed: 0,
                    n_observed_in_sample: 0,
                });
            }
            continue;
        }
        let ranks = average_ranks(gene_to_value);
        for (set_name, set_genes) in gene_sets {
            let observed_ranks: Vec<f64> = set_genes
                .iter()
                .filter_map(|g| ranks.get(g.as_str()).copied())
                .collect();
            let m = observed_ranks.len();
            let score = if m < min_set_size || m >= n {
                f64::NAN
            } else {
                let mean_rank: f64 = observed_ranks.iter().sum::<f64>() / m as f64;
                centered_singscore(mean_rank, n, m)
            };
            out.push(SingscoreRow {
                sample_id: sample_id.clone(),
                set_name: set_name.clone(),
                score,
                set_size_declared: set_genes.len(),
                set_size_observed: m,
                n_observed_in_sample: n,
            });
        }
    }
    out.sort_by(|a, b| {
        a.set_name
            .cmp(&b.set_name)
            .then_with(|| a.sample_id.cmp(&b.sample_id))
    });
    out
}

/// Centered, scale-normalized singscore: matches `singscore::simpleScore`
/// with `centerScore = TRUE`.
///
/// Range: `[-0.5, 0.5]`.
fn centered_singscore(mean_rank: f64, n: usize, m: usize) -> f64 {
    debug_assert!(
        n > m,
        "n must exceed m for centered singscore to be defined"
    );
    let denom = 2.0 * (n - m) as f64;
    (2.0 * mean_rank - (n + 1) as f64) / denom
}

/// Average-rank assignment matching R's `rank(x, ties.method = "average")`.
/// Ties are assigned the mean of the rank positions they would have
/// occupied in a strict ordering. Returns `gene_symbol -> rank` (1-based).
fn average_ranks(values: &BTreeMap<String, f64>) -> BTreeMap<&str, f64> {
    let mut entries: Vec<(&str, f64)> = values.iter().map(|(g, v)| (g.as_str(), *v)).collect();
    entries.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));
    let mut out = BTreeMap::new();
    let mut i = 0usize;
    while i < entries.len() {
        let mut j = i + 1;
        while j < entries.len() && entries[j].1 == entries[i].1 {
            j += 1;
        }
        // Tied block: positions [i+1, j] in 1-based ranks; average is (i+1 + j) / 2.
        let avg_rank = ((i + 1 + j) as f64) / 2.0;
        for entry in &entries[i..j] {
            out.insert(entry.0, avg_rank);
        }
        i = j;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_with(values: &[(&str, f64)]) -> BTreeMap<String, f64> {
        values.iter().map(|(g, v)| (g.to_string(), *v)).collect()
    }

    fn one_set(name: &str, genes: &[&str]) -> BTreeMap<String, BTreeSet<String>> {
        let mut m = BTreeMap::new();
        m.insert(
            name.to_string(),
            genes.iter().map(|g| g.to_string()).collect(),
        );
        m
    }

    #[test]
    fn signature_at_top_scores_positive_half() {
        // 10 proteins; signature is the top 3 (highest values).
        let mut abundance = BTreeMap::new();
        let mut s1 = BTreeMap::new();
        for i in 0..10 {
            s1.insert(format!("G{:02}", i), i as f64);
        }
        abundance.insert("S1".to_string(), s1);
        let sets = one_set("top3", &["G07", "G08", "G09"]);
        let out = singscore(&abundance, &sets, 1);
        assert_eq!(out.len(), 1);
        // Mean rank of top 3 = 9. n=10, m=3.
        // score = (2*9 - 11) / (2*7) = 7/14 = 0.5
        assert!((out[0].score - 0.5).abs() < 1e-12, "got {}", out[0].score);
        assert_eq!(out[0].set_size_observed, 3);
    }

    #[test]
    fn signature_at_bottom_scores_negative_half() {
        let mut abundance = BTreeMap::new();
        let mut s1 = BTreeMap::new();
        for i in 0..10 {
            s1.insert(format!("G{:02}", i), i as f64);
        }
        abundance.insert("S1".to_string(), s1);
        let sets = one_set("bot3", &["G00", "G01", "G02"]);
        let out = singscore(&abundance, &sets, 1);
        // Mean rank = 2. score = (2*2 - 11) / 14 = -7/14 = -0.5
        assert!((out[0].score + 0.5).abs() < 1e-12, "got {}", out[0].score);
    }

    #[test]
    fn signature_in_middle_scores_zero() {
        let mut abundance = BTreeMap::new();
        let mut s1 = BTreeMap::new();
        for i in 0..11 {
            s1.insert(format!("G{:02}", i), i as f64);
        }
        abundance.insert("S1".to_string(), s1);
        // Signature = single middle gene (rank 6 of 11). score = (12-12)/(2*10) = 0.
        let sets = one_set("mid", &["G05"]);
        let out = singscore(&abundance, &sets, 1);
        assert!(out[0].score.abs() < 1e-12, "got {}", out[0].score);
    }

    #[test]
    fn average_ranks_handles_ties_correctly() {
        // values: G1=1, G2=2, G3=2, G4=2, G5=3
        // strict sorted positions: G1(1), G2/G3/G4(2,3,4), G5(5)
        // average ranks: G1=1, G2=G3=G4=3, G5=5
        let m = sample_with(&[
            ("G1", 1.0),
            ("G2", 2.0),
            ("G3", 2.0),
            ("G4", 2.0),
            ("G5", 3.0),
        ]);
        let r = average_ranks(&m);
        assert!((r["G1"] - 1.0).abs() < 1e-12);
        assert!((r["G2"] - 3.0).abs() < 1e-12);
        assert!((r["G3"] - 3.0).abs() < 1e-12);
        assert!((r["G4"] - 3.0).abs() < 1e-12);
        assert!((r["G5"] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn missing_genes_in_sample_are_excluded() {
        // S1 has G0..G9; S2 has G0..G4 only. Signature = {G5, G6, G7}.
        // S1 mean rank = mean(6,7,8) = 7, m=3, n=10. score = (14-11)/14 = 0.214...
        // S2 has m=0 in sample → score = NaN with min_set_size=1.
        let mut abundance = BTreeMap::new();
        let mut s1 = BTreeMap::new();
        for i in 0..10 {
            s1.insert(format!("G{i}"), i as f64);
        }
        let mut s2 = BTreeMap::new();
        for i in 0..5 {
            s2.insert(format!("G{i}"), i as f64);
        }
        abundance.insert("S1".to_string(), s1);
        abundance.insert("S2".to_string(), s2);
        let sets = one_set("upper", &["G5", "G6", "G7"]);
        let out = singscore(&abundance, &sets, 1);
        assert_eq!(out.len(), 2);
        let s1_row = out.iter().find(|r| r.sample_id == "S1").unwrap();
        let s2_row = out.iter().find(|r| r.sample_id == "S2").unwrap();
        // S1: mean_rank = 7, score = (14 - 11) / 14 = 3/14
        assert!((s1_row.score - 3.0 / 14.0).abs() < 1e-12);
        assert_eq!(s1_row.set_size_observed, 3);
        assert!(s2_row.score.is_nan());
        assert_eq!(s2_row.set_size_observed, 0);
    }

    #[test]
    fn min_set_size_filter_emits_nan_with_observed_count() {
        let mut abundance = BTreeMap::new();
        let mut s1 = BTreeMap::new();
        for i in 0..10 {
            s1.insert(format!("G{i}"), i as f64);
        }
        abundance.insert("S1".to_string(), s1);
        // Set has 5 genes but only 1 observed in S1.
        let sets = one_set("partial", &["G3", "X1", "X2", "X3", "X4"]);
        let out = singscore(&abundance, &sets, 3);
        assert_eq!(out.len(), 1);
        assert!(out[0].score.is_nan());
        assert_eq!(out[0].set_size_observed, 1);
        assert_eq!(out[0].set_size_declared, 5);
    }

    #[test]
    fn output_is_sorted_by_set_then_sample() {
        let mut abundance = BTreeMap::new();
        for sid in ["S2", "S1", "S3"] {
            let mut s = BTreeMap::new();
            for i in 0..6 {
                s.insert(format!("G{i}"), i as f64);
            }
            abundance.insert(sid.to_string(), s);
        }
        let mut sets = BTreeMap::new();
        for sn in ["beta", "alpha"] {
            sets.insert(sn.to_string(), {
                let mut s = BTreeSet::new();
                s.insert("G0".to_string());
                s.insert("G1".to_string());
                s
            });
        }
        let out = singscore(&abundance, &sets, 1);
        assert_eq!(out.len(), 6);
        let order: Vec<(&str, &str)> = out
            .iter()
            .map(|r| (r.set_name.as_str(), r.sample_id.as_str()))
            .collect();
        assert_eq!(
            order,
            vec![
                ("alpha", "S1"),
                ("alpha", "S2"),
                ("alpha", "S3"),
                ("beta", "S1"),
                ("beta", "S2"),
                ("beta", "S3"),
            ]
        );
    }
}
