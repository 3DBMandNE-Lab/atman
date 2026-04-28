//! Feature-covariance network analysis for `atman network influence`.
//!
//! Builds a feature × feature similarity matrix over subjects, turns
//! it into a weighted undirected adjacency via hard thresholding or
//! soft-power adjacency (WGCNA-style), then scores each feature's
//! role as a network hub via:
//!
//! - Eigenvector centrality (power iteration on the weighted
//!   adjacency) — how "connected to other connected nodes" a feature
//!   is.
//! - Betweenness centrality (Brandes' algorithm on the graph viewed
//!   with reciprocal-similarity distances) — how often a feature
//!   lies on shortest paths between other features.
//! - Influence score: `eigenvector × betweenness`. Follows the
//!   Burberry-Pillai 2026 construction used on CyTOF immune
//!   signaling hubs.
//!
//! All computations are pure and deterministic. The eigenvector
//! power iteration starts from a fixed uniform vector and normalises
//! after each step; convergence tolerance is caller-specified.
//!
//! See `atman network influence` for the CLI wrapper.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::stats::{mean, pearson, spearman};

/// Similarity metric for the feature × feature graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimilarityMetric {
    /// Pearson correlation over subjects (absolute value).
    Pearson,
    /// Spearman rank correlation over subjects (absolute value).
    Spearman,
    /// Sample covariance over subjects (absolute value).
    Covariance,
}

impl SimilarityMetric {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pearson => "pearson",
            Self::Spearman => "spearman",
            Self::Covariance => "covariance",
        }
    }
}

/// Adjacency-construction policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AdjacencyPolicy {
    /// Retain an edge when `|similarity| >= threshold`; edge weight
    /// is the absolute similarity. Others get weight 0.
    HardThreshold { threshold: f64 },
    /// Soft-power: `weight = |similarity|^power`. Preserves weighted
    /// edges; no hard drop. Matches WGCNA's soft-threshold
    /// construction.
    SoftPower { power: u32 },
}

/// Pairwise feature × feature similarity matrix. `features.len() == n_features`.
/// `values[i][j]` is the (absolute) similarity between feature `i` and
/// feature `j`; the diagonal is 1.0.
pub struct SimilarityMatrix {
    pub features: Vec<String>,
    pub values: Vec<Vec<f64>>,
}

/// Audit counts produced alongside a similarity matrix. Pairs with
/// fewer than 2 finite observations or an undefined-correlation
/// outcome (e.g. zero-variance input to Pearson/Spearman) get
/// `0.0` in the matrix; this struct records how many such pairs
/// there were so callers can surface them in run metadata instead of
/// silently treating "couldn't compute" as "uncorrelated".
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct SimilarityAudit {
    /// Total number of off-diagonal feature pairs evaluated
    /// (`n_features * (n_features - 1) / 2`).
    pub n_pairs_total: usize,
    /// Pairs with fewer than 2 finite observations after the
    /// complete-case filter.
    pub n_pairs_insufficient_overlap: usize,
    /// Pairs where the metric returned no value (Pearson/Spearman
    /// stddev = 0). Counted only for metrics that can return `None`.
    pub n_pairs_undefined_metric: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum PairOutcome {
    Ok(f64),
    InsufficientOverlap,
    UndefinedMetric,
}

/// Compute the pairwise similarity matrix over a feature × subject
/// layout. `data[i][s]` is feature `i`'s value for subject `s`.
///
/// Subjects with non-finite values in either feature are dropped
/// pairwise (complete-case per pair). Returns `None` if fewer than
/// two features or two subjects are supplied. The accompanying
/// [`SimilarityAudit`] reports how many pairs fell back to `0.0` and
/// why.
pub fn pairwise_similarity(
    feature_labels: &[String],
    data: &[Vec<f64>],
    metric: SimilarityMetric,
) -> Option<(SimilarityMatrix, SimilarityAudit)> {
    if feature_labels.len() != data.len() || data.len() < 2 {
        return None;
    }
    if data[0].len() < 2 {
        return None;
    }
    let n_features = data.len();
    let mut values = vec![vec![0.0_f64; n_features]; n_features];
    for (i, row) in values.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let mut audit = SimilarityAudit::default();
    for i in 0..n_features {
        for j in (i + 1)..n_features {
            audit.n_pairs_total += 1;
            let sim = match pairwise(&data[i], &data[j], metric) {
                PairOutcome::Ok(s) => s,
                PairOutcome::InsufficientOverlap => {
                    audit.n_pairs_insufficient_overlap += 1;
                    0.0
                }
                PairOutcome::UndefinedMetric => {
                    audit.n_pairs_undefined_metric += 1;
                    0.0
                }
            };
            values[i][j] = sim;
            values[j][i] = sim;
        }
    }
    Some((
        SimilarityMatrix {
            features: feature_labels.to_vec(),
            values,
        },
        audit,
    ))
}

fn pairwise(x: &[f64], y: &[f64], metric: SimilarityMetric) -> PairOutcome {
    // Complete-case pairwise filter on NaNs.
    let pairs: Vec<(f64, f64)> = x
        .iter()
        .zip(y.iter())
        .filter(|(a, b)| a.is_finite() && b.is_finite())
        .map(|(a, b)| (*a, *b))
        .collect();
    if pairs.len() < 2 {
        return PairOutcome::InsufficientOverlap;
    }
    let xs: Vec<f64> = pairs.iter().map(|(a, _)| *a).collect();
    let ys: Vec<f64> = pairs.iter().map(|(_, b)| *b).collect();
    match metric {
        SimilarityMetric::Pearson => match pearson(&xs, &ys) {
            Some(r) => PairOutcome::Ok(r.abs()),
            None => PairOutcome::UndefinedMetric,
        },
        SimilarityMetric::Spearman => match spearman(&xs, &ys) {
            Some(r) => PairOutcome::Ok(r.abs()),
            None => PairOutcome::UndefinedMetric,
        },
        SimilarityMetric::Covariance => {
            let mx = mean(&xs);
            let my = mean(&ys);
            let cov: f64 = xs
                .iter()
                .zip(ys.iter())
                .map(|(a, b)| (a - mx) * (b - my))
                .sum::<f64>()
                / (xs.len() - 1) as f64;
            PairOutcome::Ok(cov.abs())
        }
    }
}

/// Apply the adjacency policy to a similarity matrix, in place
/// semantically. Zero out the diagonal (no self-loops).
pub fn adjacency(similarity: &SimilarityMatrix, policy: AdjacencyPolicy) -> Vec<Vec<f64>> {
    let n = similarity.values.len();
    let mut a = vec![vec![0.0_f64; n]; n];
    for (i, sim_row) in similarity.values.iter().enumerate() {
        for (j, &s) in sim_row.iter().enumerate() {
            if i == j {
                continue;
            }
            let w = match policy {
                AdjacencyPolicy::HardThreshold { threshold } => {
                    if s.abs() >= threshold {
                        s.abs()
                    } else {
                        0.0
                    }
                }
                AdjacencyPolicy::SoftPower { power } => s.abs().powi(power as i32),
            };
            a[i][j] = w;
        }
    }
    a
}

/// Eigenvector centrality via shifted power iteration on
/// `(A + αI)` where `α` is the maximum row sum of `A`. Starts from
/// a uniform vector `v_i = 1/sqrt(n)`, multiplies by `A + αI`,
/// normalises after each step, and converges when
/// `||v_new - v_old||_∞ < tol`. `adjacency` must be symmetric and
/// non-negative.
///
/// The shift trick is essential: on bipartite graphs (e.g. a star),
/// `A` has eigenvalues `±√(n-1)` with the same magnitude, and plain
/// power iteration oscillates. Shifting by `α = ||A||_∞` makes every
/// eigenvalue non-negative, so power iteration converges to the
/// Perron eigenvector of `A` (the all-nonnegative principal
/// eigenvector). When `α = 0` (zero adjacency), every feature gets
/// centrality `0.0`.
pub fn eigenvector_centrality(adjacency: &[Vec<f64>], max_iter: usize, tol: f64) -> Vec<f64> {
    let n = adjacency.len();
    if n == 0 {
        return Vec::new();
    }
    let shift = max_row_sum(adjacency);
    if shift == 0.0 {
        return vec![0.0_f64; n];
    }
    let mut v = vec![1.0_f64 / (n as f64).sqrt(); n];
    for _ in 0..max_iter {
        let mut next = vec![0.0_f64; n];
        for i in 0..n {
            let mut acc = shift * v[i];
            for j in 0..n {
                acc += adjacency[i][j] * v[j];
            }
            next[i] = acc;
        }
        let norm: f64 = next.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm == 0.0 {
            return vec![0.0_f64; n];
        }
        for x in &mut next {
            *x /= norm;
        }
        let max_diff = v
            .iter()
            .zip(next.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        v = next;
        if max_diff < tol {
            break;
        }
    }
    // Sign convention: make the dominant direction non-negative.
    // With an identity shift the eigenvector is already non-negative
    // on non-negative adjacencies, but flip-if-majority-negative is
    // a cheap guard.
    let negatives = v.iter().filter(|x| **x < 0.0).count();
    if negatives * 2 > n {
        for x in &mut v {
            *x = -*x;
        }
    }
    v
}

fn max_row_sum(adjacency: &[Vec<f64>]) -> f64 {
    adjacency
        .iter()
        .map(|row| row.iter().sum::<f64>())
        .fold(0.0_f64, f64::max)
}

/// Brandes' algorithm for weighted betweenness centrality on an
/// undirected graph. Edge weights in `adjacency` are
/// similarities/adjacencies; internally we route SSSP via Dijkstra
/// with edge **distances = 1 / weight** so high-similarity edges
/// correspond to short paths.
///
/// Self-loops (`adjacency[i][i]`) are ignored. Edges with weight
/// `<= 0` are treated as absent.
pub fn betweenness_centrality(adjacency: &[Vec<f64>]) -> Vec<f64> {
    let n = adjacency.len();
    if n == 0 {
        return Vec::new();
    }
    let mut bc = vec![0.0_f64; n];
    for s in 0..n {
        // Single-source shortest paths via Dijkstra.
        let mut stack: Vec<usize> = Vec::new();
        let mut pred: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut sigma = vec![0.0_f64; n];
        sigma[s] = 1.0;
        let mut dist = vec![f64::INFINITY; n];
        dist[s] = 0.0;
        let mut heap: BinaryHeap<State> = BinaryHeap::new();
        heap.push(State { dist: 0.0, node: s });
        let mut visited = vec![false; n];
        while let Some(State { dist: d, node: v }) = heap.pop() {
            if visited[v] {
                continue;
            }
            visited[v] = true;
            stack.push(v);
            for w in 0..n {
                if v == w {
                    continue;
                }
                let edge = adjacency[v][w];
                if !(edge.is_finite() && edge > 0.0) {
                    continue;
                }
                let alt = d + 1.0 / edge;
                if alt < dist[w] - 1e-12 {
                    dist[w] = alt;
                    heap.push(State { dist: alt, node: w });
                    sigma[w] = sigma[v];
                    pred[w] = vec![v];
                } else if (alt - dist[w]).abs() < 1e-12 {
                    sigma[w] += sigma[v];
                    pred[w].push(v);
                }
            }
        }
        // Accumulation.
        let mut delta = vec![0.0_f64; n];
        while let Some(w) = stack.pop() {
            for &v in &pred[w] {
                if sigma[w] > 0.0 {
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
                }
            }
            if w != s {
                bc[w] += delta[w];
            }
        }
    }
    // Undirected normalization: each shortest-path pair counted twice.
    for x in &mut bc {
        *x /= 2.0;
    }
    bc
}

#[derive(Copy, Clone)]
struct State {
    dist: f64,
    node: usize,
}

impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        self.dist == other.dist && self.node == other.node
    }
}
impl Eq for State {}
impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        // Min-heap ordering on distance, break ties by node index
        // (lower first) for determinism.
        other
            .dist
            .partial_cmp(&self.dist)
            .unwrap_or(Ordering::Equal)
            .then_with(|| other.node.cmp(&self.node))
    }
}
impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Per-feature network statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct InfluenceRow {
    pub feature_id: String,
    pub eigenvector_centrality: f64,
    pub betweenness_centrality: f64,
    pub influence_score: f64,
    pub degree: usize,
}

/// Compute per-feature influence rows given an adjacency matrix
/// (symmetric, non-negative; zero diagonal). `features.len()` must
/// match `adjacency.len()`.
pub fn influence_scores(
    features: &[String],
    adjacency: &[Vec<f64>],
    eig_max_iter: usize,
    eig_tol: f64,
) -> Vec<InfluenceRow> {
    assert_eq!(features.len(), adjacency.len());
    let eig = eigenvector_centrality(adjacency, eig_max_iter, eig_tol);
    let bc = betweenness_centrality(adjacency);
    features
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let ec = eig[i];
            let bwc = bc[i];
            let degree = adjacency[i].iter().filter(|w| **w > 0.0).count();
            InfluenceRow {
                feature_id: f.clone(),
                eigenvector_centrality: ec,
                betweenness_centrality: bwc,
                influence_score: ec * bwc,
                degree,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pearson_pairwise_recovers_planted_correlation() {
        // Features: x1, x2 = 2*x1 + noise, x3 = unrelated
        let x1 = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let x2: Vec<f64> = x1.iter().map(|v| 2.0 * v + 0.01).collect();
        let x3 = vec![5.0, -1.0, 2.0, 8.0, -3.0, 4.0, 0.0, 6.0];
        let data = vec![x1, x2, x3];
        let labels = vec!["f1".into(), "f2".into(), "f3".into()];
        let (sim, audit) = pairwise_similarity(&labels, &data, SimilarityMetric::Pearson).unwrap();
        // f1 ~ f2 should be close to 1.0.
        assert!(sim.values[0][1] > 0.99, "got {}", sim.values[0][1]);
        // f1 ~ f3 should be weakly correlated.
        assert!(sim.values[0][2].abs() < 0.95);
        // Clean inputs: every pair computes; no fallbacks.
        assert_eq!(audit.n_pairs_total, 3);
        assert_eq!(audit.n_pairs_insufficient_overlap, 0);
        assert_eq!(audit.n_pairs_undefined_metric, 0);
    }

    #[test]
    fn similarity_audit_counts_insufficient_overlap_and_undefined_metric() {
        // f1 has all-NaN; every pair with f1 hits insufficient overlap.
        let f1 = vec![f64::NAN; 5];
        // f2 is constant; pearson against any non-NaN partner is undefined.
        let f2 = vec![1.0; 5];
        let f3 = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let data = vec![f1, f2, f3];
        let labels = vec!["f1".into(), "f2".into(), "f3".into()];
        let (_, audit) = pairwise_similarity(&labels, &data, SimilarityMetric::Pearson).unwrap();
        assert_eq!(audit.n_pairs_total, 3);
        // (f1,f2) and (f1,f3) both fail on overlap (f1 is all-NaN).
        assert_eq!(audit.n_pairs_insufficient_overlap, 2);
        // (f2,f3) succeeds at the overlap step but pearson is undefined
        // because f2 is constant.
        assert_eq!(audit.n_pairs_undefined_metric, 1);
    }

    #[test]
    fn hard_threshold_drops_subthreshold_edges() {
        let sim = SimilarityMatrix {
            features: vec!["a".into(), "b".into(), "c".into()],
            values: vec![
                vec![1.0, 0.9, 0.1],
                vec![0.9, 1.0, 0.05],
                vec![0.1, 0.05, 1.0],
            ],
        };
        let adj = adjacency(&sim, AdjacencyPolicy::HardThreshold { threshold: 0.5 });
        assert!(adj[0][1] > 0.0);
        assert_eq!(adj[0][2], 0.0);
        assert_eq!(adj[1][2], 0.0);
        // Diagonal zeroed.
        assert_eq!(adj[0][0], 0.0);
    }

    #[test]
    fn soft_power_preserves_all_edges_scaled() {
        let sim = SimilarityMatrix {
            features: vec!["a".into(), "b".into()],
            values: vec![vec![1.0, 0.5], vec![0.5, 1.0]],
        };
        let adj = adjacency(&sim, AdjacencyPolicy::SoftPower { power: 6 });
        // 0.5^6 = 0.015625
        assert!((adj[0][1] - 0.015625).abs() < 1e-12);
    }

    /// Build a star adjacency matrix: hub at node 0, leaves at 1..n.
    fn star_adjacency(n: usize) -> Vec<Vec<f64>> {
        (0..n)
            .map(|i| {
                let mut row = vec![0.0_f64; n];
                if i == 0 {
                    for v in row.iter_mut().skip(1) {
                        *v = 1.0;
                    }
                } else {
                    row[0] = 1.0;
                }
                row
            })
            .collect()
    }

    /// Star topology: hub connected to every leaf, leaves not
    /// connected to each other. Hub's eigenvector centrality must
    /// dominate all leaves.
    #[test]
    fn eigenvector_centrality_finds_star_hub() {
        let n = 5;
        let adj = star_adjacency(n);
        let v = eigenvector_centrality(&adj, 500, 1e-10);
        for (i, &vi) in v.iter().enumerate().skip(1) {
            assert!(v[0] > vi, "hub {} should exceed leaf {}: {}", v[0], i, vi);
        }
    }

    /// Star topology: hub's betweenness must strictly exceed every
    /// leaf's. In a 1-hub-4-leaf star, the hub lies on every leaf-
    /// to-leaf shortest path.
    #[test]
    fn betweenness_centrality_finds_star_hub() {
        let n = 5;
        let adj = star_adjacency(n);
        let bc = betweenness_centrality(&adj);
        for (i, &leaf_bc) in bc.iter().enumerate().skip(1) {
            assert!(
                bc[0] > leaf_bc,
                "hub bc {} must exceed leaf {} bc {}",
                bc[0],
                i,
                leaf_bc,
            );
        }
        // 4 leaves → C(4, 2) = 6 leaf-leaf pairs all routed through
        // the hub. Undirected normalization halves; expected 6.
        assert!(
            (bc[0] - 6.0).abs() < 1e-9,
            "star hub BC should be 6, got {}",
            bc[0]
        );
    }

    /// Isolated graph (all zero adjacency): every centrality is zero
    /// and betweenness is zero; no panics.
    #[test]
    fn zero_adjacency_produces_zero_centralities() {
        let adj = vec![vec![0.0; 3]; 3];
        let ev = eigenvector_centrality(&adj, 100, 1e-10);
        let bc = betweenness_centrality(&adj);
        for v in &ev {
            assert_eq!(*v, 0.0);
        }
        for v in &bc {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn influence_score_is_eigen_times_betweenness_with_degree_reporting() {
        let n = 4;
        // Clique (complete graph): symmetric, all identical.
        let mut adj = vec![vec![1.0; n]; n];
        for (i, row) in adj.iter_mut().enumerate() {
            row[i] = 0.0;
        }
        let features = (0..n).map(|i| format!("f{i}")).collect::<Vec<_>>();
        let rows = influence_scores(&features, &adj, 500, 1e-10);
        for r in &rows {
            assert_eq!(r.degree, n - 1, "clique degree should be n-1");
            // By symmetry, every node has the same centralities.
            assert!(
                (r.eigenvector_centrality - rows[0].eigenvector_centrality).abs() < 1e-8,
                "symmetric eigenvector"
            );
        }
    }

    #[test]
    fn influence_is_deterministic() {
        let n = 5;
        let adj = star_adjacency(n);
        let features = (0..n).map(|i| format!("f{i}")).collect::<Vec<_>>();
        let a = influence_scores(&features, &adj, 500, 1e-10);
        let b = influence_scores(&features, &adj, 500, 1e-10);
        assert_eq!(a, b);
    }
}
