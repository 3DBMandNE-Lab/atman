//! Differential coexpression across multiple cohorts.
//!
//! Given per-cohort feature × subject matrices on the same feature
//! universe, computes one signed correlation matrix per cohort and then
//! summarizes how each protein-pair edge varies across cohorts in one of
//! three modes:
//!
//! - **edge-pairwise** (DGCA-style): per-edge × per-cohort-pair Fisher-z
//!   test of `H_0: corr_A = corr_B`. Output cardinality is
//!   `O(n_edges × n_cohort_pairs)`. Good for asking "where do exactly
//!   these two cohorts diverge?", less good for visualizing.
//! - **edge-summary**: per-edge cross-cohort summary (mean, sd, min/max
//!   absolute correlation, sign-flip count, conservation and divergence
//!   scores). Output cardinality is `O(n_edges)`. Best fit for
//!   "conserved vs. cohort-specific" pan-cancer maps.
//! - **module**: per-module within-module connectivity per cohort, with
//!   cross-cohort summary statistics. Requires user-supplied modules
//!   (e.g. MSigDB Hallmarks). Output cardinality is `O(n_modules)`.
//!
//! All correlations are signed Pearson (default) or signed Spearman.
//! Sample sizes per cohort are tracked separately because the Fisher-z
//! test and confidence intervals depend on `n - 3`.

use std::collections::{BTreeMap, BTreeSet};

use crate::stats::{mean, pearson, spearman};

/// Signed correlation metric for differential coexpression. Distinct
/// from [`crate::network::SimilarityMetric`], which always takes
/// absolute values. Differential coexpression specifically needs the
/// sign so that sign-flip edges across cohorts can be detected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SignedMetric {
    Pearson,
    Spearman,
}

impl SignedMetric {
    pub fn as_str(self) -> &'static str {
        match self {
            SignedMetric::Pearson => "pearson",
            SignedMetric::Spearman => "spearman",
        }
    }
}

/// One cohort's feature × subject matrix. `data[i][s]` is feature `i`'s
/// value for subject `s`. NaNs / non-finites are dropped pairwise per
/// edge (complete-case per pair).
pub struct CohortData<'a> {
    pub label: String,
    pub feature_labels: &'a [String],
    pub data: &'a [Vec<f64>],
}

/// Per-cohort signed correlation matrix on the shared feature universe.
/// `values[i][j]` is the signed correlation between feature `i` and
/// feature `j` in this cohort. `n_pairs[i][j]` is the number of subjects
/// with finite values for both features (= `n_eff` for the Fisher-z
/// formula). Diagonal is `1.0` / `n_subjects`.
pub struct CohortCorrelations {
    pub label: String,
    pub features: Vec<String>,
    pub values: Vec<Vec<f64>>,
    pub n_pairs: Vec<Vec<usize>>,
}

/// Compute signed pairwise correlations for a cohort.
///
/// Returns `None` if fewer than 2 features or fewer than `min_overlap`
/// per-pair finite samples for any pair (we keep the partial matrix and
/// flag the sparse pairs via `n_pairs`).
pub fn signed_pairwise_correlations(
    cohort: &CohortData<'_>,
    metric: SignedMetric,
    min_overlap: usize,
) -> Option<CohortCorrelations> {
    if cohort.feature_labels.len() != cohort.data.len() || cohort.data.len() < 2 {
        return None;
    }
    let n_features = cohort.data.len();
    let mut values = vec![vec![0.0_f64; n_features]; n_features];
    let mut n_pairs = vec![vec![0usize; n_features]; n_features];
    for i in 0..n_features {
        values[i][i] = 1.0;
        n_pairs[i][i] = cohort.data[i].iter().filter(|v| v.is_finite()).count();
    }
    for i in 0..n_features {
        for j in (i + 1)..n_features {
            let pairs: Vec<(f64, f64)> = cohort.data[i]
                .iter()
                .zip(cohort.data[j].iter())
                .filter(|(a, b)| a.is_finite() && b.is_finite())
                .map(|(a, b)| (*a, *b))
                .collect();
            n_pairs[i][j] = pairs.len();
            n_pairs[j][i] = pairs.len();
            if pairs.len() < min_overlap {
                continue;
            }
            let xs: Vec<f64> = pairs.iter().map(|p| p.0).collect();
            let ys: Vec<f64> = pairs.iter().map(|p| p.1).collect();
            let r = match metric {
                SignedMetric::Pearson => pearson(&xs, &ys),
                SignedMetric::Spearman => spearman(&xs, &ys),
            };
            if let Some(r) = r {
                values[i][j] = r;
                values[j][i] = r;
            }
        }
    }
    Some(CohortCorrelations {
        label: cohort.label.clone(),
        features: cohort.feature_labels.to_vec(),
        values,
        n_pairs,
    })
}

/// Edge-level pairwise differential row (mode `edge-pairwise`).
#[derive(Debug, Clone)]
pub struct EdgePairwiseRow {
    pub feature_a: String,
    pub feature_b: String,
    pub cohort_a: String,
    pub cohort_b: String,
    pub n_a: usize,
    pub n_b: usize,
    pub corr_a: f64,
    pub corr_b: f64,
    pub z_diff: f64,
    pub p_value: f64,
}

/// Edge-level cross-cohort summary row (mode `edge-summary`).
#[derive(Debug, Clone)]
pub struct EdgeSummaryRow {
    pub feature_a: String,
    pub feature_b: String,
    pub n_cohorts: usize,
    pub per_cohort_corrs: Vec<(String, f64)>,
    pub mean_corr: f64,
    pub sd_corr: f64,
    pub min_abs_corr: f64,
    pub max_abs_corr: f64,
    pub range_corr: f64,
    pub n_sign_flips: usize,
    pub conservation_score: f64,
    pub divergence_score: f64,
}

/// Per-module rewiring row (mode `module`).
#[derive(Debug, Clone)]
pub struct ModuleRewiringRow {
    pub set_name: String,
    pub set_size_declared: usize,
    pub set_size_observed: usize,
    pub per_cohort_connectivity: Vec<(String, f64)>,
    pub mean_connectivity: f64,
    pub sd_connectivity: f64,
    pub range_connectivity: f64,
    pub rewiring_score: f64,
}

/// Pairwise Fisher-z differential coexpression test across every
/// (edge, cohort-pair) combination.
///
/// For each protein pair `(i, j)` and each ordered pair of distinct
/// cohorts `(A, B)` with `A < B` lexicographically:
/// `Z = (atanh(r_A) - atanh(r_B)) / sqrt(1/(n_A - 3) + 1/(n_B - 3))`,
/// `p = 2 * (1 - Φ(|Z|))`. Edges with `n < 4` in either cohort are
/// skipped; edges with `|r| ≥ 1.0 - eps` are clamped to `±(1 - 1e-9)`
/// before the atanh to keep the test finite.
pub fn edge_pairwise_differential(
    cohorts: &[&CohortCorrelations],
) -> Vec<EdgePairwiseRow> {
    if cohorts.len() < 2 {
        return Vec::new();
    }
    let features = &cohorts[0].features;
    let n_features = features.len();
    debug_assert!(cohorts.iter().all(|c| c.features == *features));
    let mut out = Vec::new();
    for i in 0..n_features {
        for j in (i + 1)..n_features {
            for a in 0..cohorts.len() {
                for b in (a + 1)..cohorts.len() {
                    let n_a = cohorts[a].n_pairs[i][j];
                    let n_b = cohorts[b].n_pairs[i][j];
                    if n_a < 4 || n_b < 4 {
                        continue;
                    }
                    let r_a = clamp_corr(cohorts[a].values[i][j]);
                    let r_b = clamp_corr(cohorts[b].values[i][j]);
                    let z_a = r_a.atanh();
                    let z_b = r_b.atanh();
                    let se = (1.0 / ((n_a - 3) as f64) + 1.0 / ((n_b - 3) as f64)).sqrt();
                    let z_diff = (z_a - z_b) / se;
                    let p_value = two_sided_normal_p(z_diff.abs());
                    out.push(EdgePairwiseRow {
                        feature_a: features[i].clone(),
                        feature_b: features[j].clone(),
                        cohort_a: cohorts[a].label.clone(),
                        cohort_b: cohorts[b].label.clone(),
                        n_a,
                        n_b,
                        corr_a: cohorts[a].values[i][j],
                        corr_b: cohorts[b].values[i][j],
                        z_diff,
                        p_value,
                    });
                }
            }
        }
    }
    out
}

/// Per-edge summary across cohorts.
///
/// For each protein pair, computes mean / sd / range / sign-flip count
/// of the per-cohort correlation vector, plus two scores:
///
/// - `conservation_score = min(|r_k|)` over cohorts. High if every
///   cohort has a meaningful (signed) coexpression at this edge.
/// - `divergence_score = sd(r) / (|mean(r)| + 1e-3)`. A coefficient of
///   variation; high if the per-cohort correlations are dispersed
///   relative to their mean.
///
/// Edges where any cohort has fewer than `min_overlap_per_cohort`
/// observations are skipped.
pub fn edge_summary_differential(
    cohorts: &[&CohortCorrelations],
    min_overlap_per_cohort: usize,
) -> Vec<EdgeSummaryRow> {
    if cohorts.len() < 2 {
        return Vec::new();
    }
    let features = &cohorts[0].features;
    let n_features = features.len();
    debug_assert!(cohorts.iter().all(|c| c.features == *features));
    let mut out = Vec::new();
    for i in 0..n_features {
        for j in (i + 1)..n_features {
            if cohorts
                .iter()
                .any(|c| c.n_pairs[i][j] < min_overlap_per_cohort)
            {
                continue;
            }
            let per_cohort_corrs: Vec<(String, f64)> = cohorts
                .iter()
                .map(|c| (c.label.clone(), c.values[i][j]))
                .collect();
            let r_vec: Vec<f64> = per_cohort_corrs.iter().map(|(_, r)| *r).collect();
            let m = mean(&r_vec);
            let sd = sample_sd(&r_vec, m);
            let abs_corrs: Vec<f64> = r_vec.iter().map(|r| r.abs()).collect();
            let min_abs = abs_corrs
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min);
            let max_abs = abs_corrs
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            let range = max_abs - min_abs;
            let n_sign_flips = sign_flip_count(&r_vec);
            let conservation = min_abs;
            let divergence = sd / (m.abs() + 1e-3);
            out.push(EdgeSummaryRow {
                feature_a: features[i].clone(),
                feature_b: features[j].clone(),
                n_cohorts: cohorts.len(),
                per_cohort_corrs,
                mean_corr: m,
                sd_corr: sd,
                min_abs_corr: min_abs,
                max_abs_corr: max_abs,
                range_corr: range,
                n_sign_flips,
                conservation_score: conservation,
                divergence_score: divergence,
            });
        }
    }
    out
}

/// Per-module within-module connectivity across cohorts.
///
/// For each module `M` (gene set), within-module connectivity in cohort
/// `k` is the mean of `|r_k(i, j)|` over all unordered pairs `(i, j)`
/// with both `i, j` in `M ∩ features`. The rewiring score is the sd of
/// the per-cohort connectivity vector — high when a module's internal
/// covariance is reorganized across tumor types.
pub fn module_rewiring(
    cohorts: &[&CohortCorrelations],
    modules: &BTreeMap<String, BTreeSet<String>>,
) -> Vec<ModuleRewiringRow> {
    if cohorts.is_empty() || modules.is_empty() {
        return Vec::new();
    }
    let features = &cohorts[0].features;
    let feature_index: BTreeMap<&str, usize> = features
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i))
        .collect();
    let mut out = Vec::new();
    for (set_name, set_genes) in modules {
        let positions: Vec<usize> = set_genes
            .iter()
            .filter_map(|g| feature_index.get(g.as_str()).copied())
            .collect();
        let m = positions.len();
        if m < 2 {
            continue;
        }
        let mut per_cohort = Vec::with_capacity(cohorts.len());
        for cohort in cohorts {
            let mut acc = 0.0_f64;
            let mut n = 0usize;
            for a in 0..positions.len() {
                for b in (a + 1)..positions.len() {
                    let pi = positions[a];
                    let pj = positions[b];
                    if cohort.n_pairs[pi][pj] >= 2 {
                        acc += cohort.values[pi][pj].abs();
                        n += 1;
                    }
                }
            }
            let conn = if n > 0 { acc / n as f64 } else { 0.0 };
            per_cohort.push((cohort.label.clone(), conn));
        }
        let conn_vec: Vec<f64> = per_cohort.iter().map(|(_, c)| *c).collect();
        let mean_conn = mean(&conn_vec);
        let sd_conn = sample_sd(&conn_vec, mean_conn);
        let max_conn = conn_vec.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let min_conn = conn_vec.iter().copied().fold(f64::INFINITY, f64::min);
        out.push(ModuleRewiringRow {
            set_name: set_name.clone(),
            set_size_declared: set_genes.len(),
            set_size_observed: m,
            per_cohort_connectivity: per_cohort,
            mean_connectivity: mean_conn,
            sd_connectivity: sd_conn,
            range_connectivity: max_conn - min_conn,
            rewiring_score: sd_conn,
        });
    }
    out
}

fn clamp_corr(r: f64) -> f64 {
    let eps = 1e-9;
    r.clamp(-1.0 + eps, 1.0 - eps)
}

fn sign_flip_count(values: &[f64]) -> usize {
    values
        .windows(2)
        .filter(|pair| {
            let a = pair[0];
            let b = pair[1];
            (a > 0.0 && b < 0.0) || (a < 0.0 && b > 0.0)
        })
        .count()
}

fn sample_sd(values: &[f64], mu: f64) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let var = values.iter().map(|v| (v - mu).powi(2)).sum::<f64>() / (values.len() - 1) as f64;
    var.sqrt()
}

/// Two-sided normal-tail p-value: `2 * (1 - Φ(z))` for `z >= 0`. Uses
/// the standard `erfc` formulation via `statrs` to avoid pulling a new
/// dep — but `statrs` is already in the workspace, so we use its
/// `Normal::cdf` directly.
fn two_sided_normal_p(z: f64) -> f64 {
    use statrs::distribution::{ContinuousCDF, Normal};
    if !z.is_finite() {
        return f64::NAN;
    }
    let n = Normal::new(0.0, 1.0).expect("standard normal");
    let p = 2.0 * (1.0 - n.cdf(z.abs()));
    p.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cohort(label: &str, features: &[String], cols: Vec<Vec<f64>>) -> CohortCorrelations {
        let cd = CohortData {
            label: label.to_string(),
            feature_labels: features,
            data: Box::leak(Box::new(cols)),
        };
        signed_pairwise_correlations(&cd, SignedMetric::Pearson, 2).unwrap()
    }

    #[test]
    fn signed_correlations_preserve_sign() {
        let features: Vec<String> = vec!["A".into(), "B".into()];
        let data = vec![
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
            vec![5.0, 4.0, 3.0, 2.0, 1.0], // perfectly anti-correlated
        ];
        let cd = CohortData {
            label: "c1".into(),
            feature_labels: &features,
            data: &data,
        };
        let cc = signed_pairwise_correlations(&cd, SignedMetric::Pearson, 2).unwrap();
        assert!((cc.values[0][1] - (-1.0)).abs() < 1e-12);
        assert!((cc.values[1][0] - (-1.0)).abs() < 1e-12);
        assert_eq!(cc.n_pairs[0][1], 5);
    }

    #[test]
    fn edge_summary_flags_sign_flip_across_cohorts() {
        let features = vec!["A".into(), "B".into(), "C".into()];
        let pos = cohort(
            "pos",
            &features,
            vec![
                vec![1.0, 2.0, 3.0, 4.0, 5.0],
                vec![1.5, 2.1, 3.4, 3.9, 5.2],
                vec![5.0, 4.0, 3.0, 2.0, 1.0],
            ],
        );
        let neg = cohort(
            "neg",
            &features,
            vec![
                vec![1.0, 2.0, 3.0, 4.0, 5.0],
                vec![5.2, 3.9, 3.4, 2.1, 1.5],
                vec![5.0, 4.0, 3.0, 2.0, 1.0],
            ],
        );
        let summary = edge_summary_differential(&[&pos, &neg], 2);
        let ab = summary
            .iter()
            .find(|r| r.feature_a == "A" && r.feature_b == "B")
            .unwrap();
        assert_eq!(ab.n_sign_flips, 1, "A-B should sign-flip across cohorts");
        let ac = summary
            .iter()
            .find(|r| r.feature_a == "A" && r.feature_b == "C")
            .unwrap();
        assert_eq!(ac.n_sign_flips, 0, "A-C is anti-correlated in both cohorts");
    }

    #[test]
    fn edge_summary_conservation_high_when_all_cohorts_agree() {
        let features = vec!["A".into(), "B".into()];
        let c1 = cohort(
            "c1",
            &features,
            vec![
                vec![1.0, 2.0, 3.0, 4.0, 5.0],
                vec![1.1, 2.1, 3.0, 4.1, 4.9],
            ],
        );
        let c2 = cohort(
            "c2",
            &features,
            vec![
                vec![1.0, 2.0, 3.0, 4.0, 5.0],
                vec![1.0, 2.2, 3.1, 3.9, 5.1],
            ],
        );
        let s = edge_summary_differential(&[&c1, &c2], 2);
        let ab = &s[0];
        assert!(
            ab.conservation_score > 0.95,
            "conservation should be high when both cohorts strongly agree, got {}",
            ab.conservation_score
        );
        assert_eq!(ab.n_sign_flips, 0);
    }

    #[test]
    fn edge_pairwise_fisher_z_detects_strong_difference() {
        let n = 50usize;
        let xs: Vec<f64> = (0..n).map(|i| i as f64).collect();
        // Cohort A: y perfectly correlates with x (r ≈ +1).
        let ys_pos: Vec<f64> = xs.iter().map(|&x| x + 0.05).collect();
        // Cohort B: y perfectly anti-correlates (r ≈ -1).
        let ys_neg: Vec<f64> = xs.iter().map(|&x| -x + 0.05).collect();
        let features = vec!["X".into(), "Y".into()];
        let a = cohort("A", &features, vec![xs.clone(), ys_pos]);
        let b = cohort("B", &features, vec![xs.clone(), ys_neg]);
        let rows = edge_pairwise_differential(&[&a, &b]);
        assert_eq!(rows.len(), 1);
        let r = &rows[0];
        assert!(r.corr_a > 0.99, "cohort A should be ~+1, got {}", r.corr_a);
        assert!(r.corr_b < -0.99, "cohort B should be ~-1, got {}", r.corr_b);
        assert!(
            r.p_value < 1e-10,
            "Fisher-z should be highly significant at n=50 with r flipping sign, got {}",
            r.p_value
        );
    }

    #[test]
    fn module_rewiring_high_when_within_module_corrs_diverge() {
        let features = vec!["A".into(), "B".into(), "C".into()];
        // Cohort 1: all three strongly correlated with each other.
        let c1 = cohort(
            "c1",
            &features,
            vec![
                vec![1.0, 2.0, 3.0, 4.0, 5.0],
                vec![1.0, 2.1, 3.0, 4.1, 5.0],
                vec![1.0, 2.0, 3.1, 3.9, 5.1],
            ],
        );
        // Cohort 2: A independent of B and C.
        let c2 = cohort(
            "c2",
            &features,
            vec![
                vec![5.0, 1.0, 4.0, 2.0, 3.0],
                vec![1.0, 2.1, 3.0, 4.1, 5.0],
                vec![1.0, 2.0, 3.1, 3.9, 5.1],
            ],
        );
        let mut modules = BTreeMap::new();
        let mut s = BTreeSet::new();
        s.insert("A".to_string());
        s.insert("B".to_string());
        s.insert("C".to_string());
        modules.insert("triplet".to_string(), s);
        let rows = module_rewiring(&[&c1, &c2], &modules);
        assert_eq!(rows.len(), 1);
        let r = &rows[0];
        assert_eq!(r.set_size_observed, 3);
        let conn1 = r.per_cohort_connectivity[0].1;
        let conn2 = r.per_cohort_connectivity[1].1;
        assert!(
            conn1 > 0.9,
            "cohort 1 connectivity should be high (all three correlated), got {conn1}"
        );
        assert!(
            conn2 < conn1,
            "cohort 2 connectivity should be lower (A randomized), got {conn2}"
        );
        assert!(r.rewiring_score > 0.0);
    }

    #[test]
    fn min_overlap_filter_skips_low_n_edges() {
        let features = vec!["A".into(), "B".into()];
        let c1 = cohort(
            "c1",
            &features,
            vec![
                vec![1.0, 2.0, 3.0, 4.0, 5.0],
                vec![f64::NAN, f64::NAN, 3.0, 4.0, f64::NAN],
            ],
        );
        let c2 = cohort(
            "c2",
            &features,
            vec![
                vec![1.0, 2.0, 3.0, 4.0, 5.0],
                vec![1.1, 2.1, 3.0, 4.1, 5.0],
            ],
        );
        // c1's A-B edge has n=2; require min 3 → edge dropped.
        let s = edge_summary_differential(&[&c1, &c2], 3);
        assert!(s.is_empty(), "edge should be filtered for low overlap");
    }
}
