//! Data-driven module discovery on a subject × protein matrix.
//!
//! Soft-thresholded signed-|r| adjacency + topological overlap matrix
//! (TOM) + UPGMA hierarchical clustering on `1 − TOM`, with a fixed
//! height cut and a minimum module-size filter. Matches the core of
//! WGCNA's workflow without importing the R package.
//!
//! v1 methods:
//! - `wgcna-soft`: soft-power `|r|^β` adjacency; auto-β picks the
//!   smallest integer β ∈ {1..20} whose scale-free topology `R² ≥ 0.8`
//!   and slope < 0.
//! - `hard-threshold`: binary adjacency at `|r| ≥ threshold`; connected
//!   components become modules.
//!
//! Consensus clustering (resample subjects, aggregate co-assignments)
//! and dynamic tree-cut are documented follow-ons. v1 uses a fixed
//! `--cut-height` on the UPGMA dendrogram plus a `--min-module-size`
//! floor; features not assigned to any surviving module land in the
//! `grey` catch-all, matching WGCNA convention.

use crate::stats;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Similarity {
    Pearson,
    Spearman,
}

#[derive(Debug, Clone)]
pub struct ModuleAssignment {
    pub feature: String,
    /// `"grey"` for unassigned features; otherwise `"M01"`, `"M02"`, …
    pub module: String,
}

#[derive(Debug, Clone)]
pub struct ModuleReportRow {
    pub module: String,
    pub size: usize,
    pub mean_within_abs_correlation: f64,
    pub hub_feature: String,
    /// Fraction of within-module variance explained by the module's
    /// first principal component (the WGCNA "module eigengene"
    /// variance-explained metric).
    pub eigenprotein_pc1_variance_explained: f64,
}

#[derive(Debug, Clone)]
pub struct SoftPowerSweepRow {
    pub beta: usize,
    pub scale_free_r_squared: f64,
    pub slope: f64,
    pub mean_k: f64,
}

#[derive(Debug, Clone)]
pub struct DiscoveryResult {
    pub modules: Vec<ModuleAssignment>,
    pub report: Vec<ModuleReportRow>,
    pub soft_power_sweep: Vec<SoftPowerSweepRow>,
    /// Chosen β (0 for `hard-threshold`).
    pub soft_power_chosen: usize,
    /// Full outcome of the automatic sweep, when one ran. `None` for an
    /// explicit `--soft-power` and for `hard-threshold`.
    pub soft_power_selection: Option<SoftPowerSelection>,
    /// Audit counts for the |similarity| matrix construction.
    pub similarity_audit: SimilarityAudit,
}

/// Audit counts produced by [`pairwise_abs_similarity`]. Pairs whose
/// metric is undefined (e.g. zero-variance Pearson/Spearman input)
/// fall back to `0.0` in the matrix; this struct records how many.
/// The caller is expected to surface non-zero values in the run
/// sidecar so a silent zero is never confused with a true uncorrelated
/// pair.
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct SimilarityAudit {
    /// Total off-diagonal feature pairs evaluated
    /// (`p * (p - 1) / 2`).
    pub n_pairs_total: usize,
    /// Pairs where Pearson/Spearman returned `None`.
    pub n_pairs_undefined_metric: usize,
}

/// Pairwise |similarity| over subjects. `data` is subject-major:
/// `data[i][j]` is subject i's value at feature j. Returns a
/// `p × p` matrix with 1.0 on the diagonal alongside an audit of
/// fallback occurrences.
pub fn pairwise_abs_similarity(
    data: &[Vec<f64>],
    sim: Similarity,
) -> (Vec<Vec<f64>>, SimilarityAudit) {
    if data.is_empty() {
        return (Vec::new(), SimilarityAudit::default());
    }
    let p = data[0].len();
    // Transpose to column-major: `col[j][i]` = subject i at feature j.
    let mut col: Vec<Vec<f64>> = vec![Vec::with_capacity(data.len()); p];
    for row in data {
        for (j, v) in row.iter().enumerate() {
            col[j].push(*v);
        }
    }
    let mut out = vec![vec![0.0_f64; p]; p];
    let mut audit = SimilarityAudit::default();
    for i in 0..p {
        out[i][i] = 1.0;
        for j in (i + 1)..p {
            audit.n_pairs_total += 1;
            let r = match sim {
                Similarity::Pearson => stats::pearson(&col[i], &col[j]),
                Similarity::Spearman => stats::spearman(&col[i], &col[j]),
            };
            let v = match r {
                Some(value) => value.abs(),
                None => {
                    audit.n_pairs_undefined_metric += 1;
                    0.0
                }
            };
            out[i][j] = v;
            out[j][i] = v;
        }
    }
    (out, audit)
}

/// Soft-threshold: `a_ij = |r_ij|^β` element-wise. Diagonal stays at 1.
pub fn soft_adjacency(abs_sim: &[Vec<f64>], beta: usize) -> Vec<Vec<f64>> {
    let p = abs_sim.len();
    let mut out = vec![vec![0.0_f64; p]; p];
    let b = beta as i32;
    for i in 0..p {
        for j in 0..p {
            if i == j {
                out[i][j] = 1.0;
            } else {
                out[i][j] = abs_sim[i][j].powi(b);
            }
        }
    }
    out
}

/// Topological overlap matrix. For signed-weighted adjacency `a`:
///
/// ```text
/// TOM_ij = (Σ_u a_iu · a_uj + a_ij) / (min(k_i, k_j) + 1 − a_ij)
/// ```
///
/// where `k_i = Σ_{u ≠ i} a_iu`. Diagonal = 1.
pub fn compute_tom(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let p = a.len();
    // k_i = row sum excluding diagonal.
    let k: Vec<f64> = (0..p)
        .map(|i| (0..p).map(|u| if u == i { 0.0 } else { a[i][u] }).sum())
        .collect();
    let mut out = vec![vec![0.0_f64; p]; p];
    for i in 0..p {
        out[i][i] = 1.0;
        for j in (i + 1)..p {
            // Σ_{u ≠ i, j} a_iu · a_uj.
            let l_ij: f64 = a[i]
                .iter()
                .zip(a.iter())
                .enumerate()
                .filter(|&(u, _)| u != i && u != j)
                .map(|(_, (&aiu, row_u))| aiu * row_u[j])
                .sum();
            let num = l_ij + a[i][j];
            let den = k[i].min(k[j]) + 1.0 - a[i][j];
            let tom = if den > 0.0 { num / den } else { 0.0 };
            out[i][j] = tom;
            out[j][i] = tom;
        }
    }
    out
}

/// Scale-free topology R² + slope for a given adjacency. Histograms
/// the connectivity `k_i = Σ_{u≠i} a_iu` into `n_bins` log-bins,
/// fits `log10 p(k) = α + slope · log10 k`, returns `(R², slope, mean_k)`.
pub fn scale_free_r_squared(a: &[Vec<f64>], n_bins: usize) -> (f64, f64, f64) {
    let p = a.len();
    if p == 0 {
        return (0.0, 0.0, 0.0);
    }
    let ks: Vec<f64> = (0..p)
        .map(|i| (0..p).map(|u| if u == i { 0.0 } else { a[i][u] }).sum())
        .collect();
    let mean_k = ks.iter().sum::<f64>() / ks.len() as f64;
    let max_k = ks.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if max_k <= 0.0 || n_bins < 2 {
        return (0.0, 0.0, mean_k);
    }
    // Uniform log-bins on (0, max_k].
    let lo = (ks
        .iter()
        .cloned()
        .filter(|v| *v > 0.0)
        .fold(f64::INFINITY, f64::min))
    .max(1e-6);
    let hi = max_k;
    let log_lo = lo.ln();
    let log_hi = hi.ln();
    let width = (log_hi - log_lo) / n_bins as f64;
    if !(width.is_finite() && width > 0.0) {
        return (0.0, 0.0, mean_k);
    }
    let mut counts = vec![0usize; n_bins];
    let bin_centers: Vec<f64> = (0..n_bins)
        .map(|b| (log_lo + (b as f64 + 0.5) * width).exp())
        .collect();
    for &k in &ks {
        if k <= 0.0 {
            continue;
        }
        let idx = (((k.ln() - log_lo) / width) as i64)
            .max(0)
            .min(n_bins as i64 - 1) as usize;
        counts[idx] += 1;
    }
    // Fit log10 p vs log10 k on non-empty bins.
    let xs: Vec<f64> = counts
        .iter()
        .zip(bin_centers.iter())
        .filter(|(&c, _)| c > 0)
        .map(|(_, &k)| k.log10())
        .collect();
    let ys: Vec<f64> = counts
        .iter()
        .zip(bin_centers.iter())
        .filter(|(&c, _)| c > 0)
        .map(|(&c, _)| (c as f64 / p as f64).log10())
        .collect();
    if xs.len() < 3 {
        return (0.0, 0.0, mean_k);
    }
    // Simple OLS slope + R².
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for (x, y) in xs.iter().zip(ys.iter()) {
        sxx += (x - mx) * (x - mx);
        sxy += (x - mx) * (y - my);
        syy += (y - my) * (y - my);
    }
    if sxx <= 0.0 {
        return (0.0, 0.0, mean_k);
    }
    let slope = sxy / sxx;
    let r2 = if syy > 0.0 {
        (sxy * sxy) / (sxx * syy)
    } else {
        0.0
    };
    (r2, slope, mean_k)
}

/// Outcome of an automatic soft-power sweep.
///
/// `criterion_met` is the field that matters: β alone cannot tell a
/// reader whether the scale-free criterion was satisfied or whether the
/// sweep fell back, and those are very different runs.
#[derive(Debug, Clone)]
pub struct SoftPowerSelection {
    /// β to use for the adjacency.
    pub beta: usize,
    /// Whether some β in range actually met `R² ≥ r2_target` with
    /// negative slope. When false, `beta` came from the fallback rule.
    pub criterion_met: bool,
    /// Best `R²` seen anywhere in the sweep, and the β that produced it.
    /// Reported whether or not the criterion was met, so a reader can
    /// see how far short the data fell.
    pub best_r_squared: f64,
    pub best_r_squared_beta: usize,
    /// Names the fallback rule when `criterion_met` is false.
    pub fallback_rule: Option<&'static str>,
}

/// WGCNA's documented default soft power for when no β satisfies the
/// scale-free criterion, as a function of sample count.
///
/// From the WGCNA FAQ's table for **unsigned** networks, which is what
/// `|r|^β` builds: fewer than 20 samples → 9, 20–30 → 8, 30–40 → 7,
/// above 40 → 6.
///
/// The point is that the fallback must be a *typical* power, not an
/// extreme one. Falling back to the β with the largest `R²` sounds
/// reasonable and is a trap: when the scale-free fit degrades
/// monotonically — which is what a bulk proteome tends to do — the
/// largest `R²` sits at β = 1, which is no soft-thresholding at all.
/// That keeps every weak correlation in the adjacency and collapses the
/// whole feature set into the grey catch-all, i.e. it inverts the
/// purpose of the sweep at exactly the moment the sweep failed.
pub fn wgcna_default_soft_power(n_samples: usize) -> usize {
    match n_samples {
        0..=19 => 9,
        20..=30 => 8,
        31..=40 => 7,
        _ => 6,
    }
}

/// Sweep β ∈ 1..=max_beta and pick the smallest whose scale-free
/// R² ≥ r2_target AND slope < 0.
///
/// When no β qualifies, falls back to [`wgcna_default_soft_power`]
/// (clamped to the swept range) and says so in the returned
/// [`SoftPowerSelection`].
pub fn auto_soft_power(
    abs_sim: &[Vec<f64>],
    max_beta: usize,
    r2_target: f64,
    n_bins: usize,
    n_samples: usize,
) -> (SoftPowerSelection, Vec<SoftPowerSweepRow>) {
    let mut rows = Vec::with_capacity(max_beta);
    let mut best: Option<(usize, f64)> = None;
    let mut chosen: Option<usize> = None;
    for beta in 1..=max_beta {
        let adj = soft_adjacency(abs_sim, beta);
        let (r2, slope, mean_k) = scale_free_r_squared(&adj, n_bins);
        rows.push(SoftPowerSweepRow {
            beta,
            scale_free_r_squared: r2,
            slope,
            mean_k,
        });
        if chosen.is_none() && r2 >= r2_target && slope < 0.0 {
            chosen = Some(beta);
        }
        if best.map(|(_, r)| r2 > r).unwrap_or(true) {
            best = Some((beta, r2));
        }
    }
    let (best_beta, best_r2) = best.unwrap_or((1, f64::NAN));
    let selection = match chosen {
        Some(beta) => SoftPowerSelection {
            beta,
            criterion_met: true,
            best_r_squared: best_r2,
            best_r_squared_beta: best_beta,
            fallback_rule: None,
        },
        None => {
            let fallback = wgcna_default_soft_power(n_samples).clamp(1, max_beta.max(1));
            SoftPowerSelection {
                beta: fallback,
                criterion_met: false,
                best_r_squared: best_r2,
                best_r_squared_beta: best_beta,
                fallback_rule: Some("wgcna-default-by-sample-count"),
            }
        }
    };
    (selection, rows)
}

/// Average-linkage (UPGMA) hierarchical clustering on a symmetric
/// `p × p` dissimilarity matrix. Returns a vector of merge events
/// `(left_cluster, right_cluster, height, cluster_size)`; cluster
/// indices `< p` are leaves, indices `≥ p` are internal merges
/// (id = p + event_index).
pub fn upgma(dissim: &[Vec<f64>]) -> Vec<(usize, usize, f64, usize)> {
    let p = dissim.len();
    if p == 0 {
        return Vec::new();
    }
    // Active cluster ids and sizes.
    let mut active: Vec<usize> = (0..p).collect();
    let mut sizes: Vec<usize> = vec![1; p];
    let mut next_id = p;
    // Distance matrix indexed by active position.
    let mut d: Vec<Vec<f64>> = dissim.to_vec();
    let mut merges = Vec::with_capacity(p.saturating_sub(1));
    while active.len() > 1 {
        // Find the smallest off-diagonal distance.
        let (mut a, mut b, mut best) = (0usize, 1usize, f64::INFINITY);
        for (i, row) in d.iter().enumerate() {
            for (j, &dij) in row.iter().enumerate().skip(i + 1) {
                if dij < best {
                    best = dij;
                    a = i;
                    b = j;
                }
            }
        }
        let ca = active[a];
        let cb = active[b];
        let sa = sizes[a];
        let sb = sizes[b];
        let merged_size = sa + sb;
        merges.push((ca, cb, best, merged_size));
        // Compute new row: UPGMA weighted average.
        let mut new_row = vec![0.0_f64; active.len()];
        for k in 0..active.len() {
            if k == a || k == b {
                continue;
            }
            new_row[k] = (sa as f64 * d[a][k] + sb as f64 * d[b][k]) / merged_size as f64;
        }
        // Replace a's row/col with the merged row; remove b.
        for k in 0..active.len() {
            d[a][k] = new_row[k];
            d[k][a] = new_row[k];
        }
        d[a][a] = 0.0;
        // Drop b.
        d.remove(b);
        for row in &mut d {
            row.remove(b);
        }
        active[a] = next_id;
        sizes[a] = merged_size;
        active.remove(b);
        sizes.remove(b);
        next_id += 1;
    }
    merges
}

/// Cut an UPGMA tree at a fixed `height`: any merge with height ≤
/// `height` joins its two sub-clusters. Returns a feature → cluster-id
/// mapping for each of the `p` leaves.
pub fn cut_tree_by_height(
    merges: &[(usize, usize, f64, usize)],
    p: usize,
    height: f64,
) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..p + merges.len()).collect();
    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] == x {
            return x;
        }
        let r = find(parent, parent[x]);
        parent[x] = r;
        r
    }
    for (event_idx, (a, b, h, _)) in merges.iter().enumerate() {
        if *h <= height {
            let ra = find(&mut parent, *a);
            let rb = find(&mut parent, *b);
            let merged_id = p + event_idx;
            parent[ra] = merged_id;
            parent[rb] = merged_id;
        }
    }
    (0..p).map(|i| find(&mut parent, i)).collect()
}

/// Run the full discovery pipeline on a subject × protein matrix.
pub fn discover(
    data: &[Vec<f64>],
    feature_labels: &[String],
    sim: Similarity,
    method: DiscoveryMethod,
    cut_height: f64,
    min_module_size: usize,
) -> Result<DiscoveryResult, String> {
    if data.is_empty() || feature_labels.is_empty() {
        return Err("empty input matrix".into());
    }
    let p = feature_labels.len();
    if data.iter().any(|row| row.len() != p) {
        return Err("non-rectangular input matrix".into());
    }
    let (abs_sim, similarity_audit) = pairwise_abs_similarity(data, sim);

    let (adjacency, sweep, chosen_beta, soft_power_selection) = match method {
        DiscoveryMethod::WgcnaSoft {
            beta: None,
            r2_target,
            n_bins,
            max_beta,
        } => {
            let (selection, rows) =
                auto_soft_power(&abs_sim, max_beta, r2_target, n_bins, data.len());
            let picked = selection.beta;
            (
                soft_adjacency(&abs_sim, picked),
                rows,
                picked,
                Some(selection),
            )
        }
        DiscoveryMethod::WgcnaSoft { beta: Some(b), .. } => {
            (soft_adjacency(&abs_sim, b), Vec::new(), b, None)
        }
        DiscoveryMethod::HardThreshold { threshold } => {
            let mut adj = vec![vec![0.0_f64; p]; p];
            for i in 0..p {
                adj[i][i] = 1.0;
                for j in (i + 1)..p {
                    let v = if abs_sim[i][j] >= threshold { 1.0 } else { 0.0 };
                    adj[i][j] = v;
                    adj[j][i] = v;
                }
            }
            (adj, Vec::new(), 0, None)
        }
    };
    let tom = compute_tom(&adjacency);
    let mut dissim = vec![vec![0.0_f64; p]; p];
    for i in 0..p {
        for j in 0..p {
            dissim[i][j] = if i == j {
                0.0
            } else {
                (1.0 - tom[i][j]).max(0.0)
            };
        }
    }
    let merges = upgma(&dissim);
    let assignments = cut_tree_by_height(&merges, p, cut_height);

    // Rename cluster ids to compact M01/M02/…; filter by min_module_size
    // (below-threshold members go to "grey").
    use std::collections::BTreeMap;
    let mut cluster_members: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, c) in assignments.iter().enumerate() {
        cluster_members.entry(*c).or_default().push(i);
    }
    let mut module_for: Vec<String> = vec!["grey".into(); p];
    let mut next_module = 1usize;
    let mut kept_modules: Vec<(String, Vec<usize>)> = Vec::new();
    for (_c, members) in cluster_members {
        if members.len() >= min_module_size {
            let name = format!("M{next_module:02}");
            for &i in &members {
                module_for[i] = name.clone();
            }
            kept_modules.push((name, members));
            next_module += 1;
        }
    }
    let grey_members: Vec<usize> = (0..p).filter(|i| module_for[*i] == "grey").collect();
    if !grey_members.is_empty() {
        kept_modules.push(("grey".into(), grey_members));
    }
    let modules: Vec<ModuleAssignment> = feature_labels
        .iter()
        .enumerate()
        .map(|(i, label)| ModuleAssignment {
            feature: label.clone(),
            module: module_for[i].clone(),
        })
        .collect();

    // Per-module report.
    let mut report = Vec::with_capacity(kept_modules.len());
    for (name, members) in &kept_modules {
        let size = members.len();
        let mean_abs_r = {
            if size < 2 {
                0.0
            } else {
                let mut s = 0.0;
                let mut c = 0;
                for i in 0..members.len() {
                    for j in (i + 1)..members.len() {
                        s += abs_sim[members[i]][members[j]];
                        c += 1;
                    }
                }
                if c > 0 {
                    s / c as f64
                } else {
                    0.0
                }
            }
        };
        // Hub = member with highest sum of adjacency to other module members.
        let hub = members
            .iter()
            .map(|&i| {
                let k: f64 = members
                    .iter()
                    .filter(|&&j| j != i)
                    .map(|&j| adjacency[i][j])
                    .sum();
                (i, k)
            })
            .fold((members[0], f64::NEG_INFINITY), |acc, x| {
                if x.1 > acc.1 {
                    x
                } else {
                    acc
                }
            })
            .0;
        let pc1_var_explained = eigenprotein_pc1_variance_explained(data, members);
        report.push(ModuleReportRow {
            module: name.clone(),
            size,
            mean_within_abs_correlation: mean_abs_r,
            hub_feature: feature_labels[hub].clone(),
            eigenprotein_pc1_variance_explained: pc1_var_explained,
        });
    }

    Ok(DiscoveryResult {
        modules,
        report,
        soft_power_sweep: sweep,
        soft_power_chosen: chosen_beta,
        soft_power_selection,
        similarity_audit,
    })
}

/// PC1 variance fraction via power iteration on the within-module
/// sample × feature submatrix's covariance. Cheap and deterministic.
fn eigenprotein_pc1_variance_explained(data: &[Vec<f64>], members: &[usize]) -> f64 {
    if members.len() < 2 || data.is_empty() {
        return 0.0;
    }
    let n = data.len();
    let m = members.len();
    // Center columns.
    let mut x = vec![vec![0.0_f64; m]; n];
    for (j, &mi) in members.iter().enumerate() {
        let sum: f64 = data.iter().map(|row| row[mi]).sum();
        let mean = sum / n as f64;
        for (xi_row, data_row) in x.iter_mut().zip(data.iter()) {
            xi_row[j] = data_row[mi] - mean;
        }
    }
    // Total variance = trace(XᵀX) / (n − 1).
    let total: f64 = x
        .iter()
        .map(|row| row.iter().map(|v| v * v).sum::<f64>())
        .sum();
    if total <= 0.0 {
        return 0.0;
    }
    // Power iteration on XᵀX (m × m).
    let mut v = vec![1.0_f64 / (m as f64).sqrt(); m];
    for _ in 0..80 {
        // u = X v
        let u: Vec<f64> = x
            .iter()
            .map(|row| row.iter().zip(v.iter()).map(|(a, b)| a * b).sum())
            .collect();
        // v_new = Xᵀ u
        let mut v_new = vec![0.0_f64; m];
        for (j, vj) in v_new.iter_mut().enumerate() {
            *vj = x.iter().zip(u.iter()).map(|(row, &ui)| row[j] * ui).sum();
        }
        let norm = v_new.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm == 0.0 {
            return 0.0;
        }
        for vj in &mut v_new {
            *vj /= norm;
        }
        let delta: f64 = v.iter().zip(v_new.iter()).map(|(a, b)| (a - b).abs()).sum();
        v = v_new;
        if delta < 1e-10 {
            break;
        }
    }
    // Eigenvalue = v' XᵀX v = ||X v||².
    let u: Vec<f64> = x
        .iter()
        .map(|row| row.iter().zip(v.iter()).map(|(a, b)| a * b).sum())
        .collect();
    let pc1_var = u.iter().map(|v| v * v).sum::<f64>();
    if total > 0.0 {
        (pc1_var / total).clamp(0.0, 1.0)
    } else {
        0.0
    }
}

#[derive(Debug, Clone, Copy)]
pub enum DiscoveryMethod {
    WgcnaSoft {
        /// Explicit β; `None` ⇒ auto-sweep.
        beta: Option<usize>,
        r2_target: f64,
        n_bins: usize,
        max_beta: usize,
    },
    HardThreshold {
        threshold: f64,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_block_fixture(seed: u64) -> (Vec<Vec<f64>>, Vec<String>) {
        // 40 subjects × (20 block-1 + 20 block-2 + 40 background) features.
        // Two planted blocks of 20 correlated features each, plus 40
        // uncorrelated noise features. Block scores driven by a shared
        // per-subject latent so within-block |r| is high; background
        // features are independent N(0, 1) per subject/feature.
        use std::num::Wrapping;
        let n = 40usize;
        let block = 20usize;
        let noise = 40usize;
        let p = 2 * block + noise;
        let mut state = Wrapping(seed);
        let mut next = || {
            state = state * Wrapping(6364136223846793005_u64) + Wrapping(1442695040888963407_u64);
            let u = ((state.0 >> 33) as f64) / (u32::MAX as f64);
            // Crude Box-Muller on the fly using a second draw.
            state = state * Wrapping(6364136223846793005_u64) + Wrapping(1442695040888963407_u64);
            let v = ((state.0 >> 33) as f64) / (u32::MAX as f64);
            let u = u.max(1e-12);
            (-2.0_f64 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
        };
        let mut data = vec![vec![0.0_f64; p]; n];
        for row in data.iter_mut() {
            let z1 = next();
            let z2 = next();
            for slot in row.iter_mut().take(block) {
                *slot = 2.0 * z1 + 0.3 * next();
            }
            for slot in row.iter_mut().skip(block).take(block) {
                *slot = 2.0 * z2 + 0.3 * next();
            }
            for slot in row.iter_mut().skip(2 * block).take(noise) {
                *slot = next();
            }
        }
        let labels: Vec<String> = (0..p).map(|j| format!("F{j:03}")).collect();
        (data, labels)
    }

    #[test]
    fn upgma_recovers_two_planted_blocks() {
        let (data, labels) = two_block_fixture(123);
        let out = discover(
            &data,
            &labels,
            Similarity::Pearson,
            DiscoveryMethod::WgcnaSoft {
                beta: None,
                r2_target: 0.8,
                n_bins: 10,
                max_beta: 20,
            },
            0.5,
            10,
        )
        .unwrap();
        // Two non-grey modules expected.
        let non_grey: Vec<&ModuleReportRow> =
            out.report.iter().filter(|r| r.module != "grey").collect();
        assert!(
            non_grey.len() >= 2,
            "expected ≥2 non-grey modules, got {}: {:?}",
            non_grey.len(),
            out.report
        );
        // Each should have ≥ 15 members — not perfect recovery (noise
        // can steal members) but clean separation.
        let big_modules = non_grey.iter().filter(|r| r.size >= 15).count();
        assert!(
            big_modules >= 2,
            "expected ≥2 modules of size ≥ 15, got {}: {:?}",
            big_modules,
            non_grey
                .iter()
                .map(|r| (&r.module, r.size))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn eigenprotein_variance_is_between_zero_and_one() {
        let (data, labels) = two_block_fixture(456);
        let out = discover(
            &data,
            &labels,
            Similarity::Pearson,
            DiscoveryMethod::WgcnaSoft {
                beta: Some(6),
                r2_target: 0.8,
                n_bins: 10,
                max_beta: 20,
            },
            0.5,
            10,
        )
        .unwrap();
        for r in &out.report {
            assert!(
                (0.0..=1.0 + 1e-9).contains(&r.eigenprotein_pc1_variance_explained),
                "{} eigenprotein variance out of range: {}",
                r.module,
                r.eigenprotein_pc1_variance_explained
            );
            if r.module != "grey" && r.size >= 5 {
                // Non-grey block should have high PC1 fraction.
                assert!(
                    r.eigenprotein_pc1_variance_explained > 0.5,
                    "{} should have high PC1 variance explained; got {}",
                    r.module,
                    r.eigenprotein_pc1_variance_explained
                );
            }
        }
    }

    #[test]
    fn empty_input_is_rejected() {
        let empty: Vec<Vec<f64>> = Vec::new();
        let labels: Vec<String> = Vec::new();
        let err = discover(
            &empty,
            &labels,
            Similarity::Pearson,
            DiscoveryMethod::HardThreshold { threshold: 0.5 },
            0.5,
            3,
        )
        .unwrap_err();
        assert!(err.contains("empty"));
    }

    #[test]
    fn hard_threshold_path_produces_nonempty_result() {
        let (data, labels) = two_block_fixture(789);
        let out = discover(
            &data,
            &labels,
            Similarity::Pearson,
            DiscoveryMethod::HardThreshold { threshold: 0.3 },
            0.5,
            10,
        )
        .unwrap();
        assert!(!out.modules.is_empty());
    }

    /// A matrix whose scale-free fit only gets worse with β — which is
    /// what a dense bulk-proteome correlation structure tends to do —
    /// must NOT fall back to β = 1. β = 1 is no soft-thresholding at
    /// all: it keeps every weak correlation, which collapses everything
    /// into the grey catch-all and yields zero real modules.
    #[test]
    fn soft_power_fallback_is_not_the_minimum_when_the_criterion_fails() {
        // Dense, nearly-uniform similarity: no power produces a
        // scale-free degree distribution.
        let p = 40;
        let abs_sim: Vec<Vec<f64>> = (0..p)
            .map(|i| {
                (0..p)
                    .map(|j| {
                        if i == j {
                            1.0
                        } else {
                            let (a, b) = if i < j { (i, j) } else { (j, i) };
                            0.55 + 0.02 * (((a * 31 + b * 17) % 7) as f64 / 7.0)
                        }
                    })
                    .collect()
            })
            .collect();
        // r2_target = 0.99 is unreachable here, forcing the fallback.
        let (sel, rows) = auto_soft_power(&abs_sim, 20, 0.99, 10, 93);
        assert_eq!(rows.len(), 20);
        assert!(
            !sel.criterion_met,
            "fixture must not satisfy the criterion; best R² was {}",
            sel.best_r_squared
        );
        assert_ne!(
            sel.beta, 1,
            "falling back to β = 1 disables soft-thresholding entirely"
        );
        assert_eq!(
            sel.beta,
            wgcna_default_soft_power(93),
            "fallback should be WGCNA's default for this sample count"
        );
        assert_eq!(sel.fallback_rule, Some("wgcna-default-by-sample-count"));
    }

    /// The fallback follows WGCNA's documented unsigned-network table.
    #[test]
    fn wgcna_default_soft_power_matches_the_published_table() {
        assert_eq!(wgcna_default_soft_power(12), 9);
        assert_eq!(wgcna_default_soft_power(25), 8);
        assert_eq!(wgcna_default_soft_power(35), 7);
        assert_eq!(wgcna_default_soft_power(93), 6);
    }

    /// A `criterion_met` selection is reported as such, and reports the
    /// best R² regardless. Uses the same two-block fixture as the
    /// range test below, which has genuine block structure.
    #[test]
    fn soft_power_selection_reports_whether_the_criterion_was_met() {
        let (data, labels) = two_block_fixture(321);
        let out = discover(
            &data,
            &labels,
            Similarity::Pearson,
            DiscoveryMethod::WgcnaSoft {
                beta: None,
                r2_target: 0.8,
                n_bins: 10,
                max_beta: 20,
            },
            0.5,
            10,
        )
        .unwrap();
        let sel = out
            .soft_power_selection
            .expect("an auto sweep must report its selection");
        assert_eq!(sel.beta, out.soft_power_chosen);
        assert!(
            sel.best_r_squared.is_finite(),
            "best R² must be reported whether or not the criterion was met"
        );
        // Whichever way this fixture falls, the two fields must agree:
        // a met criterion has no fallback rule, and vice versa.
        assert_eq!(sel.criterion_met, sel.fallback_rule.is_none());
    }

    /// An explicit `--soft-power` is not a selection, so no selection is
    /// reported and no fallback can be silently attributed to it.
    #[test]
    fn explicit_soft_power_reports_no_selection() {
        let (data, labels) = two_block_fixture(321);
        let out = discover(
            &data,
            &labels,
            Similarity::Pearson,
            DiscoveryMethod::WgcnaSoft {
                beta: Some(6),
                r2_target: 0.8,
                n_bins: 10,
                max_beta: 20,
            },
            0.5,
            10,
        )
        .unwrap();
        assert_eq!(out.soft_power_chosen, 6);
        assert!(out.soft_power_selection.is_none());
    }

    #[test]
    fn soft_power_auto_picks_beta_between_1_and_max_beta() {
        let (data, labels) = two_block_fixture(321);
        let out = discover(
            &data,
            &labels,
            Similarity::Pearson,
            DiscoveryMethod::WgcnaSoft {
                beta: None,
                r2_target: 0.8,
                n_bins: 10,
                max_beta: 20,
            },
            0.5,
            10,
        )
        .unwrap();
        assert!(out.soft_power_chosen >= 1 && out.soft_power_chosen <= 20);
        assert_eq!(out.soft_power_sweep.len(), 20);
    }
}
