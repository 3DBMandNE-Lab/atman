//! Cross-cohort program alignment primitives.
//!
//! Given per-cohort program loading vectors (over a shared gene-symbol universe
//! by union), compute pairwise-cohort program-to-program similarity under one
//! of three metrics (top-N Jaccard, cosine, Spearman), apply optional
//! reciprocal-best and category-agreement filters, and assemble archetypes
//! as connected components of the filtered bipartite graph.

use std::collections::{BTreeMap, BTreeSet};

use crate::stats;
pub use crate::stats::jaccard_top_n;

/// Similarity metrics.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum AlignMetric {
    Jaccard,
    Cosine,
    /// Cosine after mean-centring each loading vector.
    ///
    /// Plain cosine has a positivity floor on non-negative loadings: two
    /// NMF programs cannot have a negative cosine, and in practice they
    /// sit high simply from sharing that baseline. Measured across six
    /// CPTAC cohorts on a shared 4,375-gene universe, the median
    /// cross-cohort cosine was +0.813 for NMF against −0.004 for ICA, so
    /// a `--tau 0.30` gate admitted 99.9% of NMF pairs and 7.7% of ICA
    /// pairs. The same threshold was selective for one method and inert
    /// for the other, which makes a method comparison at a shared tau
    /// meaningless.
    ///
    /// Centring removes the floor rather than compensating for it: the
    /// centred NMF median is −0.002 with 9.4% above 0.30, against ICA's
    /// 7.7%. Raising tau for NMF instead would keep the floor and demand
    /// extreme similarity on top of it.
    CosineCentered,
    Spearman,
}

impl AlignMetric {
    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "jaccard" => Some(AlignMetric::Jaccard),
            "cosine" => Some(AlignMetric::Cosine),
            "cosine-centered" | "cosine-centred" => Some(AlignMetric::CosineCentered),
            "spearman" => Some(AlignMetric::Spearman),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            AlignMetric::Jaccard => "jaccard",
            AlignMetric::Cosine => "cosine",
            AlignMetric::CosineCentered => "cosine-centered",
            AlignMetric::Spearman => "spearman",
        }
    }
}

/// One program in the shared-label universe. `values` aligns to a shared label
/// ordering used across all cohorts. Missing labels contribute 0.
#[derive(Debug, Clone)]
pub struct AlignedProgram {
    pub cohort: String,
    pub program: String,
    pub category: Option<String>,
    pub values: Vec<f64>,
}

/// Compute pairwise-cohort similarity between two program vectors. Cosine and
/// Spearman return the *absolute* correlation so sign-ambiguous ICA programs
/// match regardless of signing. When the underlying metric is undefined
/// (e.g. zero-norm vector for cosine, constant input for Spearman) the
/// returned value is `NaN` rather than a silent `0.0` — downstream
/// alignment filters drop non-finite similarities, so a failed compute is
/// distinguishable from a true zero similarity at the value level.
pub fn similarity(a: &[f64], b: &[f64], metric: AlignMetric, top_n: usize) -> f64 {
    match metric {
        AlignMetric::Jaccard => stats::jaccard_top_n(a, b, top_n),
        AlignMetric::Cosine => stats::cosine(a, b).map(f64::abs).unwrap_or(f64::NAN),
        AlignMetric::CosineCentered => {
            let centre = |v: &[f64]| -> Vec<f64> {
                if v.is_empty() {
                    return Vec::new();
                }
                let mean = v.iter().sum::<f64>() / v.len() as f64;
                v.iter().map(|x| x - mean).collect()
            };
            stats::cosine(&centre(a), &centre(b))
                .map(f64::abs)
                .unwrap_or(f64::NAN)
        }
        AlignMetric::Spearman => stats::spearman(a, b).map(f64::abs).unwrap_or(f64::NAN),
    }
}

/// Pairwise edge accepted by the alignment filters.
#[derive(Debug, Clone, PartialEq)]
pub struct AlignEdge {
    pub program_a_index: usize,
    pub program_b_index: usize,
    pub similarity: f64,
}

/// For two cohorts (slices of programs), compute best-matches under `metric`
/// with optional reciprocal-best and category-agreement filtering.
pub fn pairwise_edges(
    a_programs: &[&AlignedProgram],
    b_programs: &[&AlignedProgram],
    metric: AlignMetric,
    top_n: usize,
    tau: f64,
    reciprocal_best: bool,
    category_constraint: bool,
) -> Vec<AlignEdge> {
    if a_programs.is_empty() || b_programs.is_empty() {
        return Vec::new();
    }
    let mut sim = vec![vec![0.0_f64; b_programs.len()]; a_programs.len()];
    for (i, a) in a_programs.iter().enumerate() {
        for (j, b) in b_programs.iter().enumerate() {
            if category_constraint && !categories_agree(&a.category, &b.category) {
                sim[i][j] = f64::NEG_INFINITY;
                continue;
            }
            sim[i][j] = similarity(&a.values, &b.values, metric, top_n);
        }
    }
    let best_a: Vec<Option<usize>> = (0..a_programs.len()).map(|i| argmax(&sim[i])).collect();
    let best_b: Vec<Option<usize>> = (0..b_programs.len()).map(|j| argmax_col(&sim, j)).collect();
    let mut edges = Vec::new();
    for i in 0..a_programs.len() {
        if let Some(j) = best_a[i] {
            if sim[i][j] < tau || !sim[i][j].is_finite() {
                continue;
            }
            if reciprocal_best && best_b[j] != Some(i) {
                continue;
            }
            edges.push(AlignEdge {
                program_a_index: i,
                program_b_index: j,
                similarity: sim[i][j],
            });
        }
    }
    edges
}

fn argmax(row: &[f64]) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;
    for (j, &v) in row.iter().enumerate() {
        if !v.is_finite() {
            continue;
        }
        match best {
            None => best = Some((j, v)),
            Some((_, cur)) if v > cur => best = Some((j, v)),
            _ => {}
        }
    }
    best.map(|(j, _)| j)
}

fn argmax_col(sim: &[Vec<f64>], col: usize) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;
    for (i, row) in sim.iter().enumerate() {
        let v = row[col];
        if !v.is_finite() {
            continue;
        }
        match best {
            None => best = Some((i, v)),
            Some((_, cur)) if v > cur => best = Some((i, v)),
            _ => {}
        }
    }
    best.map(|(i, _)| i)
}

fn categories_agree(a: &Option<String>, b: &Option<String>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => x == y,
        _ => false,
    }
}

/// Union-find archetype assembly across all programs. Returns archetype id per program.
pub fn archetypes_union_find(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for (a, b) in edges {
        let ra = find(&mut parent, *a);
        let rb = find(&mut parent, *b);
        if ra != rb {
            parent[ra] = rb;
        }
    }
    (0..n).map(|i| find(&mut parent, i)).collect()
}

/// How selective a `tau` actually is on this data, over every
/// cross-cohort program pair.
///
/// A threshold that admits nearly every candidate is not thresholding;
/// the structure is then decided entirely by reciprocal-best matching
/// and the gate is a no-op. That is invisible in the output, and it
/// makes a sweep over `tau` look robust when it is merely saturated.
///
/// Returns `(n_pairs, n_above_tau, median_similarity)`.
pub fn tau_selectivity(
    programs: &[AlignedProgram],
    metric: AlignMetric,
    top_n: usize,
    tau: f64,
) -> (usize, usize, f64) {
    let mut by_cohort: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (idx, p) in programs.iter().enumerate() {
        by_cohort.entry(p.cohort.clone()).or_default().push(idx);
    }
    let cohorts: Vec<&String> = by_cohort.keys().collect();
    let mut sims: Vec<f64> = Vec::new();
    for i in 0..cohorts.len() {
        for j in (i + 1)..cohorts.len() {
            for &ia in &by_cohort[cohorts[i]] {
                for &ib in &by_cohort[cohorts[j]] {
                    let v = similarity(&programs[ia].values, &programs[ib].values, metric, top_n);
                    if v.is_finite() {
                        sims.push(v);
                    }
                }
            }
        }
    }
    let n_above = sims.iter().filter(|v| **v >= tau).count();
    let median = if sims.is_empty() {
        f64::NAN
    } else {
        let mut sorted = sims.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        sorted[sorted.len() / 2]
    };
    (sims.len(), n_above, median)
}

/// Given a set of programs across cohorts already aligned to a shared label
/// universe, produce archetype ids (connected components of the accepted-edge
/// graph over all cohort pairs).
pub fn build_archetypes(
    programs: &[AlignedProgram],
    metric: AlignMetric,
    top_n: usize,
    tau: f64,
    reciprocal_best: bool,
    category_constraint: bool,
) -> Vec<usize> {
    let mut by_cohort: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (idx, p) in programs.iter().enumerate() {
        by_cohort.entry(p.cohort.clone()).or_default().push(idx);
    }
    let cohorts: Vec<&String> = by_cohort.keys().collect();
    let mut edges_global: Vec<(usize, usize)> = Vec::new();
    for i in 0..cohorts.len() {
        for j in (i + 1)..cohorts.len() {
            let indices_a = &by_cohort[cohorts[i]];
            let indices_b = &by_cohort[cohorts[j]];
            let a_refs: Vec<&AlignedProgram> = indices_a.iter().map(|k| &programs[*k]).collect();
            let b_refs: Vec<&AlignedProgram> = indices_b.iter().map(|k| &programs[*k]).collect();
            let edges = pairwise_edges(
                &a_refs,
                &b_refs,
                metric,
                top_n,
                tau,
                reciprocal_best,
                category_constraint,
            );
            for e in edges {
                edges_global.push((indices_a[e.program_a_index], indices_b[e.program_b_index]));
            }
        }
    }
    archetypes_union_find(programs.len(), &edges_global)
}

/// Count archetype statistics: total archetypes, multi-cohort, universal, plus
/// per-archetype category recovery rate when categories are supplied.
#[derive(Debug, Clone)]
pub struct ArchetypeSummary {
    pub n_archetypes: usize,
    pub n_multi_cohort: usize,
    pub n_universal: usize,
    pub category_recovery: f64,
}

pub fn summarize_archetypes(programs: &[AlignedProgram], labels: &[usize]) -> ArchetypeSummary {
    let cohorts: BTreeSet<&str> = programs.iter().map(|p| p.cohort.as_str()).collect();
    let total_cohorts = cohorts.len();
    let mut by_arch: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (idx, lab) in labels.iter().enumerate() {
        by_arch.entry(*lab).or_default().push(idx);
    }
    // Count only multi-member archetypes (singletons excluded from "alignment" counts).
    let mut n_archetypes = 0usize;
    let mut n_multi_cohort = 0usize;
    let mut n_universal = 0usize;
    let mut category_agree_count = 0usize;
    let mut category_total = 0usize;
    for (_lab, members) in by_arch.iter() {
        if members.len() < 2 {
            continue;
        }
        n_archetypes += 1;
        let cohorts: BTreeSet<&str> = members
            .iter()
            .map(|i| programs[*i].cohort.as_str())
            .collect();
        if cohorts.len() >= 2 {
            n_multi_cohort += 1;
        }
        if cohorts.len() == total_cohorts {
            n_universal += 1;
        }
        let cats: Vec<&str> = members
            .iter()
            .filter_map(|i| programs[*i].category.as_deref())
            .collect();
        if cats.len() >= 2 {
            category_total += 1;
            let first = cats[0];
            if cats.iter().all(|c| *c == first) {
                category_agree_count += 1;
            }
        }
    }
    let category_recovery = if category_total == 0 {
        f64::NAN
    } else {
        category_agree_count as f64 / category_total as f64
    };
    ArchetypeSummary {
        n_archetypes,
        n_multi_cohort,
        n_universal,
        category_recovery,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn similarity_delegates_to_metric() {
        // Perfect monotone pair -> |spearman| = 1.
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![10.0, 20.0, 30.0, 40.0];
        assert!((similarity(&a, &b, AlignMetric::Spearman, 4) - 1.0).abs() < 1e-12);
        // Anti-aligned cosine still reports |cosine| = 1 for archetype matching.
        let c = vec![1.0, 2.0, 3.0];
        let d = vec![-1.0, -2.0, -3.0];
        assert!((similarity(&c, &d, AlignMetric::Cosine, 3) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn reciprocal_best_filters_out_asymmetric_matches() {
        let labels = 5;
        let mk = |cohort: &str, name: &str, values: Vec<f64>| AlignedProgram {
            cohort: cohort.to_string(),
            program: name.to_string(),
            category: None,
            values,
        };
        let a1 = mk("A", "p1", vec![1.0, 0.0, 0.0, 0.0, 0.0]);
        let a2 = mk("A", "p2", vec![0.0, 1.0, 0.0, 0.0, 0.0]);
        let b1 = mk("B", "q1", vec![1.0, 0.9, 0.0, 0.0, 0.0]);
        let a_refs = vec![&a1, &a2];
        let b_refs = vec![&b1];
        // Without reciprocal-best, both p1 and p2 match q1 (their best).
        let loose = pairwise_edges(
            &a_refs,
            &b_refs,
            AlignMetric::Cosine,
            labels,
            0.0,
            false,
            false,
        );
        assert_eq!(loose.len(), 2);
        // With reciprocal-best, only p1 (q1's best) is kept.
        let strict = pairwise_edges(
            &a_refs,
            &b_refs,
            AlignMetric::Cosine,
            labels,
            0.0,
            true,
            false,
        );
        assert_eq!(strict.len(), 1);
        assert_eq!(strict[0].program_a_index, 0);
    }

    #[test]
    fn build_archetypes_multi_cohort() {
        let programs = vec![
            AlignedProgram {
                cohort: "A".to_string(),
                program: "pA".to_string(),
                category: Some("astrocyte".to_string()),
                values: vec![1.0, 0.0, 0.0, 0.0],
            },
            AlignedProgram {
                cohort: "B".to_string(),
                program: "pB".to_string(),
                category: Some("astrocyte".to_string()),
                values: vec![0.9, 0.1, 0.0, 0.0],
            },
            AlignedProgram {
                cohort: "C".to_string(),
                program: "pC".to_string(),
                category: Some("astrocyte".to_string()),
                values: vec![0.95, 0.0, 0.0, 0.05],
            },
            AlignedProgram {
                cohort: "A".to_string(),
                program: "qA".to_string(),
                category: Some("neuronal".to_string()),
                values: vec![0.0, 0.0, 1.0, 0.0],
            },
        ];
        let labels = build_archetypes(&programs, AlignMetric::Cosine, 2, 0.5, false, false);
        assert_eq!(labels[0], labels[1]);
        assert_eq!(labels[1], labels[2]);
        assert_ne!(labels[0], labels[3]);
        let summary = summarize_archetypes(&programs, &labels);
        assert_eq!(summary.n_multi_cohort, 1);
        assert_eq!(summary.n_universal, 1);
        assert!((summary.category_recovery - 1.0).abs() < 1e-12);
    }
}

#[cfg(test)]
mod centered_cosine_tests {
    use super::*;

    fn prog(cohort: &str, name: &str, values: Vec<f64>) -> AlignedProgram {
        AlignedProgram {
            cohort: cohort.to_string(),
            program: name.to_string(),
            values,
            category: None,
        }
    }

    /// Plain cosine cannot go below zero on non-negative vectors, so two
    /// unrelated non-negative programs still score high and a shared tau
    /// stops thresholding. Centring removes the floor.
    #[test]
    fn centered_cosine_removes_the_non_negative_positivity_floor() {
        // Two non-negative vectors with disjoint "signal" on a shared
        // positive baseline: unrelated, but both strictly positive.
        let a: Vec<f64> = (0..200)
            .map(|i| 1.0 + if i < 20 { 3.0 } else { 0.0 })
            .collect();
        let b: Vec<f64> = (0..200)
            .map(|i| 1.0 + if (100..120).contains(&i) { 3.0 } else { 0.0 })
            .collect();

        let raw = similarity(&a, &b, AlignMetric::Cosine, 20);
        let centered = similarity(&a, &b, AlignMetric::CosineCentered, 20);
        assert!(
            raw > 0.5,
            "plain cosine should sit high on the shared baseline, got {raw}"
        );
        assert!(
            centered < raw - 0.3,
            "centring should remove the baseline: raw {raw}, centered {centered}"
        );
    }

    /// Centring must not flatten genuine agreement.
    #[test]
    fn centered_cosine_keeps_real_similarity_high() {
        let a: Vec<f64> = (0..200)
            .map(|i| 1.0 + if i < 20 { 3.0 } else { 0.0 })
            .collect();
        let b: Vec<f64> = (0..200)
            .map(|i| 1.0 + if i < 20 { 2.7 } else { 0.05 })
            .collect();
        let centered = similarity(&a, &b, AlignMetric::CosineCentered, 20);
        assert!(
            centered > 0.9,
            "programs sharing their signal should stay similar after centring, got {centered}"
        );
    }

    /// The saturation measurement reports what a tau actually admits.
    #[test]
    fn tau_selectivity_detects_an_inert_threshold() {
        let base = |offset: usize| -> Vec<f64> {
            (0..200)
                .map(|i| {
                    1.0 + if (offset..offset + 20).contains(&i) {
                        3.0
                    } else {
                        0.0
                    }
                })
                .collect()
        };
        let programs = vec![
            prog("c1", "p1", base(0)),
            prog("c1", "p2", base(40)),
            prog("c2", "p1", base(80)),
            prog("c2", "p2", base(120)),
        ];
        let (n_pairs, n_above_raw, median_raw) =
            tau_selectivity(&programs, AlignMetric::Cosine, 20, 0.30);
        assert_eq!(n_pairs, 4, "two cohorts of two programs give four pairs");
        assert_eq!(
            n_above_raw, 4,
            "plain cosine admits every pair here, which is the inert case: median {median_raw}"
        );
        let (_, n_above_centered, _) =
            tau_selectivity(&programs, AlignMetric::CosineCentered, 20, 0.30);
        assert!(
            n_above_centered < n_above_raw,
            "centring should make the same tau selective again"
        );
    }
}
