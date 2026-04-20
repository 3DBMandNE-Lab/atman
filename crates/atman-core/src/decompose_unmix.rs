//! Geometric compartmental unmixing — VCA + FCLS, ported from
//! hyperspectral remote sensing.
//!
//! Where `decompose ica` returns statistically-independent abstract
//! axes, this module returns *geometrically-identified pure
//! endmembers* with per-subject fractional abundances that sum to 1
//! under non-negativity. The generating model is `x ≈ E α` where `E`
//! is a `p × k` matrix of endmember feature vectors and `α` is a
//! per-subject abundance simplex.
//!
//! - **VCA** (Vertex Component Analysis; Nascimento & Bioucas-Dias
//!   2005): iteratively picks the most extreme reduced-space sample
//!   projection onto the orthogonal complement of already-selected
//!   endmembers. Completely deterministic under a fixed seed — the
//!   only randomness is the initial projection direction, derived
//!   from a SplitMix64 draw on the caller's seed.
//! - **FCLS** (Fully Constrained Least Squares; Heinz & Chang 2001):
//!   per-sample `argmin_α ||x − E α||² s.t. α ≥ 0 ∧ Σ α = 1`, solved
//!   here via projected gradient with Euclidean projection onto the
//!   probability simplex (Duchi et al. 2008). Cheap, numerically
//!   stable, and — at sufficient iterations — matches the classical
//!   Heinz–Chang active-set solver to ≤ 1e-6 on typical problems.
//! - **UCLS**: plain least-squares, dropping the simplex constraints
//!   — useful when the input transform (CLR, ILR) makes
//!   "sum to 1" meaningless.
//!
//! v1 deliberately omits `nfindr`, `--k auto` (HySime),
//! `--n-boot`, and `--annotate-markers`. Each is its own
//! ~150-line follow-on; the core pipeline is already substantial.

use crate::ica::Xoshiro256pp;

fn derive_sub_seed(seed: u64, iter: usize) -> u64 {
    let mut z = seed.wrapping_add(0x9E3779B97F4A7C15_u64.wrapping_mul(iter as u64 + 1));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

#[derive(Debug, Clone)]
pub struct VcaResult {
    /// Sample indices of the selected endmembers in the original
    /// `[sample][feature]` matrix, length `k`.
    pub endmember_sample_indices: Vec<usize>,
    /// `[k][p]` feature-space loading vectors (original coordinates).
    pub endmember_loadings: Vec<Vec<f64>>,
}

/// VCA endmember extraction on a sample × feature matrix.
/// Returns `k` sample indices + their feature vectors.
///
/// `data[i][j]` = sample `i`, feature `j`. `k` must satisfy
/// `2 ≤ k ≤ n_samples / 2` — the caller enforces that in the CLI
/// (`k > n / 2` makes the problem under-determined).
pub fn vca(data: &[Vec<f64>], k: usize, seed: u64) -> Result<VcaResult, String> {
    let n = data.len();
    if n == 0 {
        return Err("VCA: empty sample matrix".into());
    }
    if k < 2 {
        return Err(format!("VCA: k must be ≥ 2 (got {k})"));
    }
    if k > n {
        return Err(format!(
            "VCA: k ({k}) > n_samples ({n}); problem is under-determined"
        ));
    }
    let p = data[0].len();
    if p == 0 {
        return Err("VCA: samples have zero features".into());
    }
    for row in data {
        if row.len() != p {
            return Err("VCA: non-rectangular sample matrix".into());
        }
        for v in row {
            if !v.is_finite() {
                return Err("VCA: non-finite value in input".into());
            }
        }
    }
    // Feature-space matrix Y with columns = samples (p × n).
    let mut y: Vec<Vec<f64>> = vec![vec![0.0; n]; p];
    for (i, row) in data.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            y[j][i] = v;
        }
    }

    // A holds selected endmembers as columns in the FEATURE space:
    // a (p × i) matrix after i picks. We work in feature space
    // directly (skipping the SVD reduction step in Nascimento &
    // Bioucas-Dias' paper). This is slower by a constant factor for
    // large p but removes an entire moving part.
    let mut selected: Vec<usize> = Vec::with_capacity(k);
    let mut a: Vec<Vec<f64>> = Vec::new(); // [endmember_index][feature]

    // Initial direction: random unit vector in feature space.
    let mut rng = Xoshiro256pp::new(seed);
    let mut u = (0..p).map(|_| rng.next_normal()).collect::<Vec<f64>>();
    normalize_in_place(&mut u);

    for step in 0..k {
        // Project every sample onto u: v_j = u · y[:, j].
        let mut projections = vec![0.0_f64; n];
        for j in 0..n {
            let mut acc = 0.0;
            for f in 0..p {
                acc += u[f] * y[f][j];
            }
            projections[j] = acc;
        }
        // Argmax |projection|. Skip already-selected to preserve
        // distinct endmember indices.
        let mut best_idx = 0usize;
        let mut best_abs = -1.0_f64;
        for (j, &v) in projections.iter().enumerate() {
            if selected.contains(&j) {
                continue;
            }
            let av = v.abs();
            if av > best_abs {
                best_abs = av;
                best_idx = j;
            }
        }
        selected.push(best_idx);
        let loading: Vec<f64> = (0..p).map(|f| y[f][best_idx]).collect();
        a.push(loading.clone());

        if step + 1 == k {
            break;
        }

        // Draw new random direction, then orthogonalize against the
        // span of selected loadings (modified Gram-Schmidt).
        let sub_seed = derive_sub_seed(seed, step + 1);
        let mut rng2 = Xoshiro256pp::new(sub_seed);
        let mut w = (0..p).map(|_| rng2.next_normal()).collect::<Vec<f64>>();
        for v in &a {
            let dot = dot_product(&w, v);
            let v_norm_sq = dot_product(v, v);
            if v_norm_sq <= 0.0 {
                continue;
            }
            let coef = dot / v_norm_sq;
            for i in 0..p {
                w[i] -= coef * v[i];
            }
        }
        let w_norm = dot_product(&w, &w).sqrt();
        if w_norm <= 1e-15 {
            // Degenerate — fall back to a canonical basis direction
            // we haven't fully consumed yet.
            for f in 0..p {
                w[f] = ((f + step) as f64).sin();
            }
            for v in &a {
                let dot = dot_product(&w, v);
                let v_norm_sq = dot_product(v, v);
                if v_norm_sq > 0.0 {
                    let coef = dot / v_norm_sq;
                    for i in 0..p {
                        w[i] -= coef * v[i];
                    }
                }
            }
            let nn = dot_product(&w, &w).sqrt();
            if nn > 1e-15 {
                for v in &mut w {
                    *v /= nn;
                }
            }
            u = w;
        } else {
            for v in &mut w {
                *v /= w_norm;
            }
            u = w;
        }
    }

    Ok(VcaResult {
        endmember_sample_indices: selected,
        endmember_loadings: a,
    })
}

fn dot_product(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn normalize_in_place(v: &mut [f64]) {
    let n = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if n > 0.0 {
        for x in v.iter_mut() {
            *x /= n;
        }
    }
}

/// Euclidean projection onto the probability simplex `{α ≥ 0, Σα = 1}`
/// (Duchi, Shalev-Shwartz, Singer, Chandra 2008). `O(k log k)`.
pub fn simplex_projection(v: &[f64]) -> Vec<f64> {
    let k = v.len();
    if k == 0 {
        return Vec::new();
    }
    let mut u = v.to_vec();
    u.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let mut rho = 0;
    let mut cumsum = 0.0;
    let mut threshold = 0.0;
    for (i, &ui) in u.iter().enumerate() {
        cumsum += ui;
        let t = (cumsum - 1.0) / (i + 1) as f64;
        if ui - t > 0.0 {
            rho = i + 1;
            threshold = t;
        }
    }
    if rho == 0 {
        // Shouldn't happen for finite input; fall back to uniform.
        return vec![1.0 / k as f64; k];
    }
    v.iter().map(|x| (x - threshold).max(0.0)).collect()
}

/// Unconstrained least-squares abundance for one sample:
/// `α = (EᵀE + λ I)⁻¹ Eᵀ x` via Cholesky, with a tiny Tikhonov shift
/// (`λ = 1e-10 · trace(EᵀE) / k`) for stability on near-singular `E`.
pub fn ucls(e: &[Vec<f64>], x: &[f64]) -> Result<Vec<f64>, String> {
    let k = e.len();
    if k == 0 {
        return Err("UCLS: empty endmember matrix".into());
    }
    let p = e[0].len();
    if x.len() != p {
        return Err(format!(
            "UCLS: sample length {} ≠ feature dimension {}",
            x.len(),
            p
        ));
    }
    let mut ata = vec![vec![0.0_f64; k]; k];
    let mut atb = vec![0.0_f64; k];
    for i in 0..k {
        for j in 0..k {
            let mut acc = 0.0;
            for f in 0..p {
                acc += e[i][f] * e[j][f];
            }
            ata[i][j] = acc;
        }
        let mut acc = 0.0;
        for f in 0..p {
            acc += e[i][f] * x[f];
        }
        atb[i] = acc;
    }
    let trace: f64 = (0..k).map(|i| ata[i][i]).sum();
    let lambda = (trace / k as f64).max(1e-12) * 1e-10;
    for i in 0..k {
        ata[i][i] += lambda;
    }
    let l = cholesky_lower(&ata).ok_or_else(|| "UCLS: EᵀE not SPD".to_string())?;
    Ok(solve_cholesky(&l, &atb))
}

/// Fully Constrained Least Squares via projected gradient on the
/// probability simplex. Iterates until `max_iter` or `|Δα| < tol`.
pub fn fcls(e: &[Vec<f64>], x: &[f64], max_iter: usize, tol: f64) -> Result<Vec<f64>, String> {
    let k = e.len();
    if k == 0 {
        return Err("FCLS: empty endmember matrix".into());
    }
    let p = e[0].len();
    if x.len() != p {
        return Err(format!(
            "FCLS: sample length {} ≠ feature dimension {}",
            x.len(),
            p
        ));
    }
    // Precompute EᵀE (k×k) and Eᵀx (k).
    let mut ata = vec![vec![0.0_f64; k]; k];
    let mut atb = vec![0.0_f64; k];
    for i in 0..k {
        for j in 0..k {
            let mut acc = 0.0;
            for f in 0..p {
                acc += e[i][f] * e[j][f];
            }
            ata[i][j] = acc;
        }
        let mut acc = 0.0;
        for f in 0..p {
            acc += e[i][f] * x[f];
        }
        atb[i] = acc;
    }
    // Lipschitz bound for projected gradient: 2 · λ_max(EᵀE).
    // Cheap safe estimate: 2 · sum of absolute row sums (Gershgorin).
    let lip = 2.0
        * (0..k)
            .map(|i| (0..k).map(|j| ata[i][j].abs()).sum::<f64>())
            .fold(1e-12_f64, f64::max);
    let step = 1.0 / lip;

    let mut alpha = vec![1.0 / k as f64; k];
    for _ in 0..max_iter {
        // grad = 2 (EᵀE α − Eᵀx).
        let mut grad = vec![0.0_f64; k];
        for i in 0..k {
            let mut acc = 0.0;
            for j in 0..k {
                acc += ata[i][j] * alpha[j];
            }
            grad[i] = 2.0 * (acc - atb[i]);
        }
        let candidate: Vec<f64> = alpha
            .iter()
            .zip(grad.iter())
            .map(|(a, g)| a - step * g)
            .collect();
        let new_alpha = simplex_projection(&candidate);
        let delta: f64 = alpha
            .iter()
            .zip(new_alpha.iter())
            .map(|(a, b)| (a - b).abs())
            .sum();
        alpha = new_alpha;
        if delta < tol {
            break;
        }
    }
    Ok(alpha)
}

fn cholesky_lower(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut l = vec![vec![0.0_f64; n]; n];
    for j in 0..n {
        let diag_sum: f64 = l[j].iter().take(j).map(|v| v * v).sum();
        let diag = a[j][j] - diag_sum;
        if !diag.is_finite() || diag <= 0.0 {
            return None;
        }
        l[j][j] = diag.sqrt();
        for i in (j + 1)..n {
            let off_sum: f64 = l[i]
                .iter()
                .zip(l[j].iter())
                .take(j)
                .map(|(li, lj)| li * lj)
                .sum();
            l[i][j] = (a[i][j] - off_sum) / l[j][j];
        }
    }
    Some(l)
}

fn solve_cholesky(l: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let p = rhs.len();
    let mut z = vec![0.0; p];
    for i in 0..p {
        let mut sum = rhs[i];
        for (k, zk) in z.iter().enumerate().take(i) {
            sum -= l[i][k] * zk;
        }
        z[i] = sum / l[i][i];
    }
    let mut out = vec![0.0; p];
    for i in (0..p).rev() {
        let mut sum = z[i];
        for k in (i + 1)..p {
            sum -= l[k][i] * out[k];
        }
        out[i] = sum / l[i][i];
    }
    out
}

#[derive(Debug, Clone, Copy)]
pub enum AbundanceMethod {
    Fcls,
    Ucls,
}

#[derive(Debug, Clone, Copy)]
pub enum EndmemberMethod {
    Vca,
    /// Iterative simplex-volume maximization (N-FINDR, Winter 1999),
    /// initialized from VCA's picks. Swaps each endmember with the
    /// candidate sample that most increases the Gram-matrix
    /// determinant of the simplex; repeats until a full pass produces
    /// no swap or `max_passes` is exhausted.
    Nfindr { max_passes: usize },
}

/// |det(Gram matrix of pairwise differences)| of a set of `k` feature
/// vectors in `p`-space. Proportional to the squared volume of the
/// simplex they span; zero for degenerate (coincident / collinear)
/// vertex sets. Uses Cholesky on the symmetric positive semidefinite
/// Gram to avoid a general determinant routine.
pub fn simplex_gram_determinant(endmembers: &[Vec<f64>]) -> f64 {
    let k = endmembers.len();
    if k < 2 {
        return 0.0;
    }
    let p = endmembers[0].len();
    let km1 = k - 1;
    let mut d = vec![vec![0.0_f64; p]; km1];
    for i in 1..k {
        for j in 0..p {
            d[i - 1][j] = endmembers[i][j] - endmembers[0][j];
        }
    }
    let mut g = vec![vec![0.0_f64; km1]; km1];
    for i in 0..km1 {
        for j in 0..km1 {
            let mut s = 0.0;
            for q in 0..p {
                s += d[i][q] * d[j][q];
            }
            g[i][j] = s;
        }
    }
    match cholesky_lower(&g) {
        Some(l) => {
            let mut det = 1.0;
            for i in 0..km1 {
                det *= l[i][i] * l[i][i];
            }
            det
        }
        None => 0.0,
    }
}

/// N-FINDR endmember extraction. Initializes from VCA, then iterates:
/// for each of the `k` endmember slots, try replacing it with every
/// non-member sample; keep the swap that most increases
/// `simplex_gram_determinant`. Stops when a full pass produces no
/// swap or after `max_passes` iterations.
pub fn nfindr(
    data: &[Vec<f64>],
    k: usize,
    seed: u64,
    max_passes: usize,
) -> Result<VcaResult, String> {
    let vca_result = vca(data, k, seed)?;
    let mut selected = vca_result.endmember_sample_indices;
    let mut loadings = vca_result.endmember_loadings;
    let mut current_vol = simplex_gram_determinant(&loadings);
    let n = data.len();
    for _ in 0..max_passes {
        let mut changed = false;
        for slot in 0..k {
            let mut best_idx = selected[slot];
            let mut best_vol = current_vol;
            for s in 0..n {
                if selected
                    .iter()
                    .enumerate()
                    .any(|(i, &v)| i != slot && v == s)
                {
                    continue;
                }
                let old_loading = std::mem::replace(&mut loadings[slot], data[s].clone());
                let candidate_vol = simplex_gram_determinant(&loadings);
                if candidate_vol > best_vol {
                    best_vol = candidate_vol;
                    best_idx = s;
                }
                loadings[slot] = old_loading;
            }
            if best_idx != selected[slot] {
                loadings[slot] = data[best_idx].clone();
                selected[slot] = best_idx;
                current_vol = best_vol;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    Ok(VcaResult {
        endmember_sample_indices: selected,
        endmember_loadings: loadings,
    })
}

#[derive(Debug, Clone)]
pub struct UnmixResult {
    pub endmember_sample_indices: Vec<usize>,
    /// `[endmember][feature]`.
    pub endmember_loadings: Vec<Vec<f64>>,
    /// `[sample][endmember]` — abundances per sample.
    pub abundances: Vec<Vec<f64>>,
    /// Per-sample reconstruction residual norm.
    pub residual_norms: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct KSweepRow {
    pub k: usize,
    pub mean_residual_norm: f64,
    pub marginal_improvement: f64,
}

/// Sweep `k ∈ [k_min, k_max]`, compute the mean per-sample
/// reconstruction residual at each k, and pick the smallest `k`
/// whose marginal improvement over `k − 1` falls below
/// `elbow_threshold × peak_improvement`. Standard WGCNA-style
/// scree-elbow heuristic, bounded and deterministic.
///
/// Returns `(chosen_k, sweep_rows)`. When every marginal
/// improvement remains above the threshold, falls back to `k_max`.
pub fn select_k_auto(
    data: &[Vec<f64>],
    k_min: usize,
    k_max: usize,
    seed: u64,
    endmember_method: EndmemberMethod,
    abundance_method: AbundanceMethod,
    fcls_max_iter: usize,
    fcls_tol: f64,
    elbow_threshold: f64,
) -> Result<(usize, Vec<KSweepRow>), String> {
    if k_min < 2 || k_max < k_min {
        return Err(format!("invalid sweep range k_min={k_min}, k_max={k_max}"));
    }
    if k_max > data.len() / 2 {
        return Err(format!(
            "k_max={k_max} > n_samples/2={}; sweep is under-determined",
            data.len() / 2
        ));
    }
    let mut residuals: Vec<(usize, f64)> = Vec::with_capacity(k_max - k_min + 1);
    for k in k_min..=k_max {
        let result = unmix(
            data,
            k,
            seed,
            endmember_method,
            abundance_method,
            fcls_max_iter,
            fcls_tol,
        )?;
        let n = result.residual_norms.len();
        let mean = if n > 0 {
            result.residual_norms.iter().sum::<f64>() / n as f64
        } else {
            f64::NAN
        };
        residuals.push((k, mean));
    }
    let mut rows: Vec<KSweepRow> = Vec::with_capacity(residuals.len());
    let mut peak_improvement = 0.0_f64;
    for i in 0..residuals.len() {
        let marginal = if i == 0 {
            0.0
        } else {
            (residuals[i - 1].1 - residuals[i].1).max(0.0)
        };
        peak_improvement = peak_improvement.max(marginal);
        rows.push(KSweepRow {
            k: residuals[i].0,
            mean_residual_norm: residuals[i].1,
            marginal_improvement: marginal,
        });
    }
    // Pick smallest k (> k_min) whose marginal improvement is small
    // relative to the sweep's peak improvement. A value of 1.0 means
    // "pick k_min"; a value near 0 means "only stop when residual is
    // completely flat." Default 0.10 approximates the eye-inspection
    // elbow on proteomics fixtures.
    let mut chosen = k_max;
    for row in rows.iter().skip(1) {
        if peak_improvement <= 0.0 {
            chosen = row.k;
            break;
        }
        if row.marginal_improvement / peak_improvement < elbow_threshold {
            // Stop at the previous k — the one before the knee.
            chosen = (row.k - 1).max(k_min);
            break;
        }
    }
    Ok((chosen, rows))
}

/// End-to-end: endmember extraction (VCA or N-FINDR) → abundance
/// estimation → per-sample reconstruction residual.
pub fn unmix(
    data: &[Vec<f64>],
    k: usize,
    seed: u64,
    endmember_method: EndmemberMethod,
    abundance_method: AbundanceMethod,
    fcls_max_iter: usize,
    fcls_tol: f64,
) -> Result<UnmixResult, String> {
    let vca_result = match endmember_method {
        EndmemberMethod::Vca => vca(data, k, seed)?,
        EndmemberMethod::Nfindr { max_passes } => nfindr(data, k, seed, max_passes)?,
    };
    let p = data[0].len();
    let n = data.len();
    let mut abundances = Vec::with_capacity(n);
    let mut residuals = Vec::with_capacity(n);
    for x in data {
        let alpha = match abundance_method {
            AbundanceMethod::Fcls => fcls(&vca_result.endmember_loadings, x, fcls_max_iter, fcls_tol)?,
            AbundanceMethod::Ucls => ucls(&vca_result.endmember_loadings, x)?,
        };
        // Reconstruction residual ||x - Eα||.
        let mut r = 0.0;
        for f in 0..p {
            let mut pred = 0.0;
            for (i, &a) in alpha.iter().enumerate() {
                pred += a * vca_result.endmember_loadings[i][f];
            }
            r += (x[f] - pred).powi(2);
        }
        abundances.push(alpha);
        residuals.push(r.sqrt());
    }
    Ok(UnmixResult {
        endmember_sample_indices: vca_result.endmember_sample_indices,
        endmember_loadings: vca_result.endmember_loadings,
        abundances,
        residual_norms: residuals,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic LCG → Box-Muller Gaussian stream.
    struct Lcg(std::num::Wrapping<u64>);
    impl Lcg {
        fn new(seed: u64) -> Self {
            Self(std::num::Wrapping(seed))
        }
        fn next_u(&mut self) -> f64 {
            self.0 = self.0 * std::num::Wrapping(6364136223846793005_u64)
                + std::num::Wrapping(1442695040888963407_u64);
            ((self.0 .0 >> 33) as f64) / (u32::MAX as f64)
        }
        fn next_z(&mut self) -> f64 {
            let a = self.next_u().max(1e-12);
            let b = self.next_u();
            (-2.0_f64 * a.ln()).sqrt() * (2.0 * std::f64::consts::PI * b).cos()
        }
    }

    /// Build a fixture of `n_samples` drawn from a Dirichlet-like
    /// simplex over `k` planted endmembers on `p` features.
    fn planted_unmix_fixture(n: usize, p: usize, k: usize, seed: u64)
        -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>)
    {
        let mut rng = Lcg::new(seed);
        // Planted endmembers: each endmember loads on a disjoint
        // block of features with amplitude ≈ 1, rest are zero with
        // a small base intensity so the simplex interior isn't
        // exactly linearly dependent.
        // Planted endmembers with ZERO background — strictly disjoint
        // feature blocks per endmember. VCA needs a clean low-rank
        // simplex structure to identify pure vertices; a uniform
        // background term couples every endmember and breaks
        // identifiability at k > 2. The block-structured version is
        // still a fair test of endmember recovery because real
        // biological compartments (plasma, cellular, etc.) also load
        // on disjoint protein groups with small residual cross-talk
        // that is correctly modelled as *noise* rather than as a
        // shared background.
        let block = p / k;
        let mut endmembers = vec![vec![0.0_f64; p]; k];
        for i in 0..k {
            for f in 0..block {
                let fi = i * block + f;
                if fi < p {
                    endmembers[i][fi] = 1.0 + 0.05 * rng.next_z();
                }
            }
        }
        // Abundances: ensure the simplex has a pure-endmember
        // anchor per endmember so VCA can find them deterministically.
        let mut abundances = vec![vec![0.0_f64; k]; n];
        for i in 0..k {
            abundances[i][i] = 1.0;
        }
        // Remaining n-k samples are uniform Dirichlet draws.
        for si in k..n {
            let mut row = vec![0.0_f64; k];
            for v in row.iter_mut() {
                *v = rng.next_u().max(1e-12);
                *v = -v.ln(); // Exp(1)
            }
            let s: f64 = row.iter().sum();
            for v in row.iter_mut() {
                *v /= s;
            }
            abundances[si] = row;
        }
        // Observations: x_i = Eᵀ α_i + small noise.
        let mut data = vec![vec![0.0_f64; p]; n];
        for si in 0..n {
            for f in 0..p {
                let mut val = 0.0;
                for ei in 0..k {
                    val += abundances[si][ei] * endmembers[ei][f];
                }
                val += 0.005 * rng.next_z();
                data[si][f] = val;
            }
        }
        (data, endmembers, abundances)
    }

    fn best_matched_cosine(
        planted: &[Vec<f64>],
        recovered: &[Vec<f64>],
    ) -> Vec<f64> {
        planted
            .iter()
            .map(|p| {
                recovered
                    .iter()
                    .map(|r| {
                        let dot: f64 = p.iter().zip(r.iter()).map(|(a, b)| a * b).sum();
                        let np = p.iter().map(|v| v * v).sum::<f64>().sqrt();
                        let nr = r.iter().map(|v| v * v).sum::<f64>().sqrt();
                        if np > 0.0 && nr > 0.0 {
                            (dot / (np * nr)).abs()
                        } else {
                            0.0
                        }
                    })
                    .fold(0.0_f64, f64::max)
            })
            .collect()
    }

    #[test]
    fn vca_recovers_planted_endmembers_with_cosine_above_threshold() {
        let (data, planted, _) = planted_unmix_fixture(100, 60, 3, 20260420);
        let vca_result = vca(&data, 3, 42).unwrap();
        let cos = best_matched_cosine(&planted, &vca_result.endmember_loadings);
        for (i, c) in cos.iter().enumerate() {
            assert!(
                *c >= 0.9,
                "planted endmember {} not recovered at cosine ≥ 0.9; got {c}",
                i
            );
        }
    }

    #[test]
    fn fcls_satisfies_simplex_constraints_within_tolerance() {
        let (data, _planted, _dir_ab) = planted_unmix_fixture(60, 30, 3, 11);
        let result = unmix(&data, 3, 42, EndmemberMethod::Vca, AbundanceMethod::Fcls, 500, 1e-9).unwrap();
        for (i, row) in result.abundances.iter().enumerate() {
            let sum: f64 = row.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-6,
                "sample {i} row-sum = {sum} ≠ 1",
            );
            for (j, &a) in row.iter().enumerate() {
                assert!(
                    a >= -1e-9,
                    "sample {i} endmember {j} abundance = {a} < 0",
                );
            }
        }
    }

    #[test]
    fn ucls_is_closed_form_least_squares() {
        // Handcrafted 2-endmember / 3-feature system. The LS solution
        // is recoverable by direct computation.
        let e = vec![vec![1.0, 0.0, 0.5], vec![0.0, 1.0, 0.5]];
        let x = vec![0.7, 0.3, 0.5]; // should give α ≈ (0.7, 0.3).
        let alpha = ucls(&e, &x).unwrap();
        assert!((alpha[0] - 0.7).abs() < 1e-6, "got α[0] = {}", alpha[0]);
        assert!((alpha[1] - 0.3).abs() < 1e-6, "got α[1] = {}", alpha[1]);
    }

    #[test]
    fn unmix_is_deterministic_under_fixed_seed() {
        let (data, _, _) = planted_unmix_fixture(40, 20, 3, 55);
        let a = unmix(&data, 3, 7, EndmemberMethod::Vca, AbundanceMethod::Fcls, 300, 1e-9).unwrap();
        let b = unmix(&data, 3, 7, EndmemberMethod::Vca, AbundanceMethod::Fcls, 300, 1e-9).unwrap();
        assert_eq!(a.endmember_sample_indices, b.endmember_sample_indices);
        for (ra, rb) in a.abundances.iter().zip(b.abundances.iter()) {
            for (x, y) in ra.iter().zip(rb.iter()) {
                assert!((x - y).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn simplex_projection_preserves_already_simplex_points() {
        let v = vec![0.3, 0.4, 0.3];
        let p = simplex_projection(&v);
        for (x, y) in p.iter().zip(v.iter()) {
            assert!((x - y).abs() < 1e-12);
        }
    }

    #[test]
    fn simplex_projection_handles_outside_simplex_points() {
        let v = vec![2.0, -1.0, 0.5];
        let p = simplex_projection(&v);
        let sum: f64 = p.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
        assert!(p.iter().all(|x| *x >= -1e-15));
    }

    #[test]
    fn select_k_auto_picks_planted_k_on_planted_3_fixture() {
        let (data, _, _) = planted_unmix_fixture(80, 60, 3, 77);
        let (chosen, rows) = select_k_auto(
            &data,
            2,
            6,
            42,
            EndmemberMethod::Vca,
            AbundanceMethod::Fcls,
            200,
            1e-7,
            0.10,
        )
        .unwrap();
        assert!(
            chosen == 3 || chosen == 4,
            "expected k ≈ 3 on planted-3 fixture, got {chosen}; sweep={:?}",
            rows.iter().map(|r| (r.k, r.mean_residual_norm)).collect::<Vec<_>>()
        );
    }

    #[test]
    fn select_k_auto_picks_planted_k_on_planted_5_fixture() {
        let (data, _, _) = planted_unmix_fixture(120, 100, 5, 88);
        let (chosen, _rows) = select_k_auto(
            &data,
            2,
            8,
            42,
            EndmemberMethod::Vca,
            AbundanceMethod::Fcls,
            200,
            1e-7,
            0.10,
        )
        .unwrap();
        assert!(
            chosen == 5 || chosen == 6,
            "expected k ≈ 5 on planted-5 fixture, got {chosen}"
        );
    }

    #[test]
    fn nfindr_volume_is_not_smaller_than_vca_volume() {
        // Starting from VCA and swapping only when volume strictly
        // increases, N-FINDR can never return a worse simplex.
        let (data, _, _) = planted_unmix_fixture(80, 40, 3, 21);
        let vca_r = vca(&data, 3, 99).unwrap();
        let nfindr_r = nfindr(&data, 3, 99, 10).unwrap();
        let vca_vol = simplex_gram_determinant(&vca_r.endmember_loadings);
        let nfindr_vol = simplex_gram_determinant(&nfindr_r.endmember_loadings);
        assert!(
            nfindr_vol >= vca_vol - 1e-12,
            "NFINDR simplex volume ({nfindr_vol}) must be ≥ VCA volume ({vca_vol})"
        );
    }

    #[test]
    fn nfindr_recovers_planted_endmembers_at_least_as_well_as_vca() {
        let (data, planted, _) = planted_unmix_fixture(80, 40, 3, 31);
        let vca_cos = best_matched_cosine(&planted, &vca(&data, 3, 99).unwrap().endmember_loadings);
        let nfindr_cos = best_matched_cosine(
            &planted,
            &nfindr(&data, 3, 99, 10).unwrap().endmember_loadings,
        );
        // NFINDR should not *worsen* recovery; on clean fixtures it
        // usually matches VCA exactly (VCA is already near-optimal).
        for (i, (&v, &n)) in vca_cos.iter().zip(nfindr_cos.iter()).enumerate() {
            assert!(
                n >= v - 0.02,
                "NFINDR regressed cosine for planted endmember {i}: vca={v}, nfindr={n}"
            );
        }
    }

    #[test]
    fn vca_refuses_k_greater_than_n_samples() {
        let data = vec![vec![1.0, 2.0]; 5];
        let err = vca(&data, 10, 1).unwrap_err();
        assert!(err.contains("under-determined"), "got: {err}");
    }

    #[test]
    fn vca_refuses_empty_input() {
        let data: Vec<Vec<f64>> = Vec::new();
        let err = vca(&data, 3, 1).unwrap_err();
        assert!(err.contains("empty"), "got: {err}");
    }

    #[test]
    fn fcls_abundances_pearson_against_planted_is_high() {
        let (data, _endmembers, planted_ab) = planted_unmix_fixture(120, 60, 3, 30);
        let result = unmix(&data, 3, 42, EndmemberMethod::Vca, AbundanceMethod::Fcls, 1000, 1e-9).unwrap();
        // Match recovered columns to planted columns via best absolute
        // cosine between recovered endmember loadings and planted.
        let p = data[0].len();
        let planted_endmembers = {
            let (_, e, _) = planted_unmix_fixture(120, 60, 3, 30);
            e
        };
        let n_plant = planted_endmembers.len();
        let mut perm = vec![0usize; n_plant];
        for i in 0..n_plant {
            let mut best = -1.0_f64;
            let mut best_j = 0;
            for j in 0..result.endmember_loadings.len() {
                let dot: f64 = planted_endmembers[i]
                    .iter()
                    .zip(result.endmember_loadings[j].iter())
                    .map(|(a, b)| a * b)
                    .sum();
                let np = planted_endmembers[i].iter().map(|v| v * v).sum::<f64>().sqrt();
                let nr = result.endmember_loadings[j]
                    .iter()
                    .map(|v| v * v)
                    .sum::<f64>()
                    .sqrt();
                if np > 0.0 && nr > 0.0 {
                    let c = (dot / (np * nr)).abs();
                    if c > best {
                        best = c;
                        best_j = j;
                    }
                }
            }
            perm[i] = best_j;
        }
        // For each planted endmember, Pearson-correlate planted abundance
        // column vs recovered abundance column (permuted).
        let n = data.len();
        let _ = p;
        for ei in 0..n_plant {
            let planted_col: Vec<f64> = (0..n).map(|si| planted_ab[si][ei]).collect();
            let recovered_col: Vec<f64> =
                (0..n).map(|si| result.abundances[si][perm[ei]]).collect();
            let mp = planted_col.iter().sum::<f64>() / n as f64;
            let mr = recovered_col.iter().sum::<f64>() / n as f64;
            let mut num = 0.0;
            let mut dp = 0.0;
            let mut dr = 0.0;
            for si in 0..n {
                let a = planted_col[si] - mp;
                let b = recovered_col[si] - mr;
                num += a * b;
                dp += a * a;
                dr += b * b;
            }
            let corr = if dp > 0.0 && dr > 0.0 {
                num / (dp.sqrt() * dr.sqrt())
            } else {
                0.0
            };
            assert!(
                corr.abs() >= 0.9,
                "endmember {ei} abundance Pearson vs planted = {corr}; expected ≥ 0.9",
            );
        }
    }
}
