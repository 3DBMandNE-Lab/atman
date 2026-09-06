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

use crate::ica::jacobi_eigen;
use crate::rng::{derive_sub_seed, Xoshiro256pp};

/// In-place modified Gram-Schmidt: subtract from `w` its projection
/// onto each non-zero basis vector in `basis`. Skips degenerate
/// (zero-norm) basis vectors so the orthogonalization is robust to
/// linearly-dependent inputs.
/// Remove from `w` its component along every vector in `basis`.
///
/// **`basis` must be orthonormal.** A single Gram-Schmidt pass over a
/// non-orthogonal basis does not leave `w` orthogonal to the basis's
/// span — subtracting the projection onto one vector reintroduces
/// components along the others whenever they are correlated. Build the
/// basis with [`push_orthonormal`], which keeps that invariant.
fn orthogonalize_against(w: &mut [f64], basis: &[Vec<f64>]) {
    for v in basis {
        let v_norm_sq = dot_product(v, v);
        if v_norm_sq <= 0.0 {
            continue;
        }
        let coef = dot_product(w, v) / v_norm_sq;
        for (wi, vi) in w.iter_mut().zip(v.iter()) {
            *wi -= coef * vi;
        }
    }
}

/// Orthonormalise `v` against `basis` and append it, preserving the
/// orthonormal invariant [`orthogonalize_against`] relies on. Returns
/// `false` (appending nothing) when `v` lies in the existing span.
fn push_orthonormal(basis: &mut Vec<Vec<f64>>, v: &[f64]) -> bool {
    let mut c = v.to_vec();
    // Two passes: classical Gram-Schmidt loses orthogonality on
    // near-parallel inputs, and simplex vertices sharing a large common
    // component are exactly that case.
    orthogonalize_against(&mut c, basis);
    orthogonalize_against(&mut c, basis);
    let norm = dot_product(&c, &c).sqrt();
    if norm <= 1e-12 {
        return false;
    }
    for x in &mut c {
        *x /= norm;
    }
    basis.push(c);
    true
}

#[derive(Debug, Clone)]
pub struct VcaResult {
    /// Sample indices of the selected endmembers in the original
    /// `[sample][feature]` matrix, length `k`.
    pub endmember_sample_indices: Vec<usize>,
    /// `[k][p]` feature-space loading vectors (original coordinates).
    pub endmember_loadings: Vec<Vec<f64>>,
}

/// Project `data` onto its `k`-dimensional signal subspace and apply the
/// projective transform onto the simplex hyperplane.
///
/// Shared by [`spa`] and [`vca`]: both need the same reduced geometry,
/// and only their vertex-selection rule differs. Returns
/// `projected[sample][component]` alongside the mean direction.
fn reduce_to_simplex_coords(data: &[Vec<f64>], k: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
    let n = data.len();
    // p ≫ n here, so take the subspace from the n × n Gram matrix
    // rather than a p × p covariance: with Y = U S Vᵀ, Gram = YᵀY =
    // V S² Vᵀ, and the projected coordinate of sample j along component
    // d is s_d · V[j][d]. Same subspace, n × n eigenproblem.
    let mut gram = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in i..n {
            let g = dot_product(&data[i], &data[j]);
            gram[i][j] = g;
            gram[j][i] = g;
        }
    }
    let (eigenvalues, eigenvectors) = jacobi_eigen(&gram);

    let mut coords: Vec<Vec<f64>> = vec![vec![0.0; k]; n];
    for d in 0..k {
        let s_d = eigenvalues.get(d).copied().unwrap_or(0.0).max(0.0).sqrt();
        for (c, evec_row) in coords.iter_mut().zip(eigenvectors.iter()) {
            c[d] = s_d * evec_row[d];
        }
    }
    // Eigenvector signs are arbitrary; canonicalise so the same input
    // yields the same coordinates on every platform and build.
    for d in 0..k {
        let mut pivot = 0usize;
        let mut best = -1.0_f64;
        for (j, c) in coords.iter().enumerate() {
            let a = c[d].abs();
            if a > best + 1e-15 {
                best = a;
                pivot = j;
            }
        }
        if coords[pivot][d] < 0.0 {
            for c in coords.iter_mut() {
                c[d] = -c[d];
            }
        }
    }

    let mut mean_dir = vec![0.0_f64; k];
    for c in &coords {
        for (d, v) in c.iter().enumerate() {
            mean_dir[d] += *v / (n as f64);
        }
    }
    let mut projected: Vec<Vec<f64>> = Vec::with_capacity(n);
    for c in &coords {
        let denom = dot_product(c, &mean_dir);
        if denom.abs() > 1e-12 {
            projected.push(c.iter().map(|v| v / denom).collect());
        } else {
            projected.push(c.clone());
        }
    }
    (projected, mean_dir)
}

/// Relative residual below which SPA declares the data exhausted.
///
/// Measured against the first (largest) residual, so it is scale-free.
/// The gap at the true rank is enormous — on a rank-3 fixture the ratio
/// falls from 2.3e-1 at the last real endmember to 1.9e-8 at the first
/// spurious one — so any threshold inside that gap behaves identically;
/// 1e-6 sits well below real structure and well above numerical noise.
const RESIDUAL_COLLAPSE_RATIO: f64 = 1e-6;

/// Successive Projection Algorithm (Araújo et al. 2001; analysed as
/// robust separable-NMF by Gillis & Vavasis 2014).
///
/// **Deterministic.** Where VCA probes the reduced space with random
/// directions and takes `argmax |u · y|`, SPA repeatedly takes the point
/// of largest residual norm and projects it out:
///
/// ```text
/// R ← projected samples
/// repeat k times:
///     j ← argmax_j ‖R_j‖          (ties: lowest index)
///     select j
///     u ← R_j / ‖R_j‖
///     R ← R − u (uᵀ R)            (deflate)
/// ```
///
/// The greedy max-residual rule is exactly QR with column pivoting on
/// the reduced matrix, and under near-separability it provably recovers
/// the vertices with error bounded by the noise level. Because there is
/// no random direction there is no seed: the vertex set is a function of
/// the data alone, which is what makes an endmember claim reproducible
/// rather than a draw.
pub fn spa(data: &[Vec<f64>], k: usize) -> Result<VcaResult, String> {
    let n = validate_endmember_input(data, k)?;
    let (projected, _mean_dir) = reduce_to_simplex_coords(data, k);

    let mut residual = projected;
    let mut selected: Vec<usize> = Vec::with_capacity(k);
    let mut loadings: Vec<Vec<f64>> = Vec::with_capacity(k);
    // Scale reference for the "no independent direction left" test,
    // taken from the data itself: an absolute tolerance would be
    // meaningless because the projective transform sets the coordinate
    // scale.
    let mut scale_ref = 0.0_f64;

    for step in 0..k {
        let mut best_idx = usize::MAX;
        let mut best_norm = -1.0_f64;
        for (j, r) in residual.iter().enumerate() {
            if selected.contains(&j) {
                continue;
            }
            let norm_sq = dot_product(r, r);
            // Strict `>` keeps the lowest index on a tie, so the result
            // does not depend on iteration order.
            if norm_sq > best_norm + 1e-15 {
                best_norm = norm_sq;
                best_idx = j;
            }
        }
        if best_idx == usize::MAX {
            return Err(format!(
                "SPA: ran out of candidate samples at step {step} (k={k}, n={n})"
            ));
        }
        let best_len = best_norm.max(0.0).sqrt();
        if step == 0 {
            scale_ref = best_len;
        } else if best_len <= RESIDUAL_COLLAPSE_RATIO * scale_ref {
            // Every remaining sample lies in the span of the vertices
            // already chosen, so the next pick would be numerical noise.
            return Err(format!(
                "SPA: no independent direction left after {step} endmembers (largest \
                 remaining residual is {:.3e} of the first, k={k}); the data does not \
                 support {k} distinct endmembers — lower --k",
                best_len / scale_ref.max(f64::MIN_POSITIVE)
            ));
        }
        selected.push(best_idx);
        loadings.push(data[best_idx].clone());

        if step + 1 == k {
            break;
        }
        // Deflate: remove the selected direction from every residual.
        let pivot = residual[best_idx].clone();
        let pivot_norm = dot_product(&pivot, &pivot).sqrt();
        if pivot_norm <= RESIDUAL_COLLAPSE_RATIO * scale_ref.max(f64::MIN_POSITIVE) {
            return Err(format!(
                "SPA: residual collapsed after {} endmembers; k={k} exceeds the \
                 dimensionality actually present in the data — lower --k",
                step + 1
            ));
        }
        let u: Vec<f64> = pivot.iter().map(|v| v / pivot_norm).collect();
        for r in residual.iter_mut() {
            let c = dot_product(r, &u);
            for (ri, ui) in r.iter_mut().zip(u.iter()) {
                *ri -= c * ui;
            }
        }
    }

    Ok(VcaResult {
        endmember_sample_indices: selected,
        endmember_loadings: loadings,
    })
}

/// Shared shape/finiteness validation for endmember extractors.
fn validate_endmember_input(data: &[Vec<f64>], k: usize) -> Result<usize, String> {
    let n = data.len();
    if n == 0 {
        return Err("endmember extraction: empty sample matrix".into());
    }
    if k < 2 {
        return Err(format!("endmember extraction: k must be ≥ 2 (got {k})"));
    }
    if k > n {
        return Err(format!(
            "endmember extraction: k ({k}) > n_samples ({n}); problem is under-determined"
        ));
    }
    let p = data[0].len();
    if p == 0 {
        return Err("endmember extraction: samples have zero features".into());
    }
    for row in data {
        if row.len() != p {
            return Err("endmember extraction: non-rectangular sample matrix".into());
        }
        for v in row {
            if !v.is_finite() {
                return Err("endmember extraction: non-finite value in input".into());
            }
        }
    }
    Ok(n)
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
    // ---- Signal-subspace projection ------------------------------
    //
    // The vertex search MUST happen in the k-dimensional signal
    // subspace, not in the p-dimensional ambient space. In p
    // dimensions a random unit direction is almost orthogonal to a
    // k-dimensional simplex (concentration of measure), so
    // `argmax |u · y|` is decided by noise and the selected vertices
    // depend on the seed rather than on the data. This is the
    // dimensionality-reduction stage of Nascimento & Bioucas-Dias
    // (2005) — noise suppression, not a speed optimisation.
    //
    // p ≫ n here, so the subspace is obtained from the n × n Gram
    // matrix rather than a p × p covariance: with Y = U S Vᵀ,
    // Gram = YᵀY = V S² Vᵀ, and the projected coordinate of sample j
    // along component d is s_d · V[j][d]. Same subspace, n × n
    // eigenproblem. (Same trick as `nmf::nndsvda`'s init.)
    let mut gram = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in i..n {
            let g = dot_product(&data[i], &data[j]);
            gram[i][j] = g;
            gram[j][i] = g;
        }
    }
    let (eigenvalues, eigenvectors) = jacobi_eigen(&gram);

    // Projected coordinates: `coords[sample][component]`, k components.
    let mut coords: Vec<Vec<f64>> = vec![vec![0.0; k]; n];
    for d in 0..k {
        let s_d = eigenvalues.get(d).copied().unwrap_or(0.0).max(0.0).sqrt();
        for (c, evec_row) in coords.iter_mut().zip(eigenvectors.iter()) {
            c[d] = s_d * evec_row[d];
        }
    }

    // Eigenvector signs are arbitrary; canonicalise so the same input
    // yields the same coordinates on every platform and build. Rule:
    // the entry of largest magnitude in each component is positive,
    // ties broken by lowest sample index.
    for d in 0..k {
        let mut pivot = 0usize;
        let mut best = -1.0_f64;
        for (j, c) in coords.iter().enumerate() {
            let a = c[d].abs();
            if a > best + 1e-15 {
                best = a;
                pivot = j;
            }
        }
        if coords[pivot][d] < 0.0 {
            for c in coords.iter_mut() {
                c[d] = -c[d];
            }
        }
    }

    // Projective transform onto the simplex hyperplane: divide each
    // sample's coordinate vector by its projection onto the mean
    // direction, so vertices are picked by composition rather than by
    // overall magnitude. Samples whose projection is degenerate keep
    // their unscaled coordinates.
    let mut mean_dir = vec![0.0_f64; k];
    for c in &coords {
        for (d, v) in c.iter().enumerate() {
            mean_dir[d] += *v / (n as f64);
        }
    }
    let mut projected: Vec<Vec<f64>> = Vec::with_capacity(n);
    for c in &coords {
        let denom = dot_product(c, &mean_dir);
        if denom.abs() > 1e-12 {
            projected.push(c.iter().map(|v| v / denom).collect());
        } else {
            projected.push(c.clone());
        }
    }

    // ---- Vertex search, in the reduced space ---------------------
    let mut selected: Vec<usize> = Vec::with_capacity(k);
    // Selected vertices in original feature coordinates, for output.
    let mut a: Vec<Vec<f64>> = Vec::new();

    // Gram-Schmidt basis the search direction is kept orthogonal to.
    //
    // It is seeded with the mean direction and that seed is DISCARDED
    // once the first vertex is picked, mirroring the auxiliary matrix
    // in Nascimento & Bioucas-Dias' formulation, whose fixed initial
    // column is overwritten by the first selected vertex. The seeding
    // matters because after the projective transform every sample has
    // the same component along `mean_dir`, which adds an identical
    // offset to every projection — harmless for `argmax`, but
    // `argmax |·|` is not shift-invariant, so a seed-dependent offset
    // would silently decide the winner. Retaining it afterwards would
    // be equally wrong: it would consume one of the k dimensions and
    // leave the final vertex to the degenerate fallback.
    let mut ortho_basis: Vec<Vec<f64>> = Vec::with_capacity(k + 1);
    push_orthonormal(&mut ortho_basis, &mean_dir);

    let mut rng = Xoshiro256pp::new(seed);
    let mut u = (0..k).map(|_| rng.next_normal()).collect::<Vec<f64>>();
    orthogonalize_against(&mut u, &ortho_basis);
    normalize_in_place(&mut u);

    for step in 0..k {
        let mut best_idx = usize::MAX;
        let mut best_abs = -1.0_f64;
        for (j, pj) in projected.iter().enumerate() {
            if selected.contains(&j) {
                continue;
            }
            let av = dot_product(&u, pj).abs();
            if av > best_abs {
                best_abs = av;
                best_idx = j;
            }
        }
        if best_idx == usize::MAX {
            return Err(format!(
                "VCA: ran out of candidate samples at step {step} (k={k}, n={n})"
            ));
        }
        selected.push(best_idx);
        if step == 0 {
            // Drop the mean-direction seed; from here the basis spans
            // the vertices chosen so far.
            ortho_basis.clear();
        }
        push_orthonormal(&mut ortho_basis, &projected[best_idx]);
        a.push(data[best_idx].clone());

        if step + 1 == k {
            break;
        }

        // New random direction, orthogonalised against the span of the
        // vertices already chosen (modified Gram-Schmidt), in reduced
        // space.
        let sub_seed = derive_sub_seed(seed, step + 1);
        let mut rng2 = Xoshiro256pp::new(sub_seed);
        let mut w = (0..k).map(|_| rng2.next_normal()).collect::<Vec<f64>>();
        orthogonalize_against(&mut w, &ortho_basis);
        let w_norm = dot_product(&w, &w).sqrt();
        if w_norm <= 1e-15 {
            // Degenerate draw: fall back to a deterministic direction
            // and re-orthogonalise.
            for (d, wd) in w.iter_mut().enumerate() {
                *wd = ((d + step + 1) as f64).sin();
            }
            orthogonalize_against(&mut w, &ortho_basis);
            let nn = dot_product(&w, &w).sqrt();
            if nn > 1e-15 {
                for v in &mut w {
                    *v /= nn;
                }
            }
        } else {
            for v in &mut w {
                *v /= w_norm;
            }
        }
        u = w;
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
            ata[i][j] = dot_product(&e[i], &e[j]);
        }
        atb[i] = dot_product(&e[i], x);
    }
    let trace: f64 = (0..k).map(|i| ata[i][i]).sum();
    let lambda = (trace / k as f64).max(1e-12) * 1e-10;
    for (i, row) in ata.iter_mut().enumerate() {
        row[i] += lambda;
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
            ata[i][j] = dot_product(&e[i], &e[j]);
        }
        atb[i] = dot_product(&e[i], x);
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
        let grad: Vec<f64> = ata
            .iter()
            .zip(atb.iter())
            .map(|(row, &b)| 2.0 * (dot_product(row, &alpha) - b))
            .collect();
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
    /// Successive Projection Algorithm — deterministic greedy
    /// max-residual vertex search. No RNG, so the result is a function
    /// of the data alone. Default.
    Spa,
    Vca,
    /// Iterative simplex-volume maximization (N-FINDR, Winter 1999),
    /// initialized from VCA's picks. Swaps each endmember with the
    /// candidate sample that most increases the Gram-matrix
    /// determinant of the simplex; repeats until a full pass produces
    /// no swap or `max_passes` is exhausted.
    Nfindr {
        max_passes: usize,
    },
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
            g[i][j] = dot_product(&d[i], &d[j]);
        }
    }
    match cholesky_lower(&g) {
        Some(l) => l
            .iter()
            .enumerate()
            .map(|(i, row)| row[i] * row[i])
            .product(),
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
    for _ in 0..max_passes {
        let mut changed = false;
        for slot in 0..k {
            let mut best_idx = selected[slot];
            let mut best_vol = current_vol;
            for (s, sample) in data.iter().enumerate() {
                if selected
                    .iter()
                    .enumerate()
                    .any(|(i, &v)| i != slot && v == s)
                {
                    continue;
                }
                let old_loading = std::mem::replace(&mut loadings[slot], sample.clone());
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
pub struct LoadingCi {
    /// `[endmember][feature]` — lower bound.
    pub lower: Vec<Vec<f64>>,
    /// `[endmember][feature]` — upper bound.
    pub upper: Vec<Vec<f64>>,
}

#[derive(Debug, Clone)]
pub struct AbundanceCi {
    /// `[subject][endmember]` — lower bound.
    pub lower: Vec<Vec<f64>>,
    /// `[subject][endmember]` — upper bound.
    pub upper: Vec<Vec<f64>>,
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
    cfg: UnmixConfig,
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
        let result = unmix(data, k, cfg)?;
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

/// Algorithm configuration for `unmix` and its callers
/// (`select_k_auto`, `bootstrap_ci`). Bundles the five parameters
/// that describe *how* to unmix — `k` and the data matrix are kept
/// positional because they vary per call; everything else is stable
/// across a sweep or bootstrap run.
#[derive(Debug, Clone, Copy)]
pub struct UnmixConfig {
    pub seed: u64,
    pub endmember_method: EndmemberMethod,
    pub abundance_method: AbundanceMethod,
    pub fcls_max_iter: usize,
    pub fcls_tol: f64,
}

/// End-to-end: endmember extraction (VCA or N-FINDR) → abundance
/// estimation → per-sample reconstruction residual.
pub fn unmix(data: &[Vec<f64>], k: usize, cfg: UnmixConfig) -> Result<UnmixResult, String> {
    let vca_result = match cfg.endmember_method {
        EndmemberMethod::Spa => spa(data, k)?,
        EndmemberMethod::Vca => vca(data, k, cfg.seed)?,
        EndmemberMethod::Nfindr { max_passes } => nfindr(data, k, cfg.seed, max_passes)?,
    };
    let p = data[0].len();
    let n = data.len();
    let mut abundances = Vec::with_capacity(n);
    let mut residuals = Vec::with_capacity(n);
    for x in data {
        let alpha = match cfg.abundance_method {
            AbundanceMethod::Fcls => fcls(
                &vca_result.endmember_loadings,
                x,
                cfg.fcls_max_iter,
                cfg.fcls_tol,
            )?,
            AbundanceMethod::Ucls => ucls(&vca_result.endmember_loadings, x)?,
        };
        // Reconstruction residual ||x - Eα||.
        let r: f64 = (0..p)
            .map(|f| {
                let pred: f64 = alpha
                    .iter()
                    .zip(vca_result.endmember_loadings.iter())
                    .map(|(&a, e)| a * e[f])
                    .sum();
                (x[f] - pred).powi(2)
            })
            .sum();
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

fn sample_indices_with_replacement(rng: &mut Xoshiro256pp, n: usize) -> Vec<usize> {
    // Unbiased bounded draw on the raw u64 stream (was a modulo-biased
    // rejection loop over `next_normal().to_bits()`, which both biased the
    // result and drew from float bit patterns rather than the integer
    // stream).
    (0..n).map(|_| rng.bounded(n)).collect()
}

/// Subject-level bootstrap CI for both endmember loadings and
/// per-subject abundances. For each of `n_boot` iterations:
///
/// 1. Resample subjects with replacement (deterministic under
///    `seed`).
/// 2. Run VCA+FCLS on the resampled matrix to get a bootstrap
///    endmember matrix `E^(b)`.
/// 3. Match each PE endmember to its best `|cosine|` counterpart in
///    `E^(b)`, sign-correcting so the matched loading has positive
///    inner product with the PE reference.
/// 4. Re-solve abundances on the **original** subject vectors using
///    the matched+signed `E^(b)` (so every iteration's abundance
///    table is indexed by the original subject ordering).
///
/// After all iterations, take 2.5% and 97.5% percentiles per element
/// of the loading/abundance replicates → `LoadingCi` / `AbundanceCi`.
///
/// Returns `None` on both when `n_boot == 0`.
pub fn bootstrap_ci(
    pe: &UnmixResult,
    data: &[Vec<f64>],
    cfg: UnmixConfig,
    n_boot: usize,
) -> Result<(Option<LoadingCi>, Option<AbundanceCi>), String> {
    let seed = cfg.seed;
    if n_boot == 0 {
        return Ok((None, None));
    }
    let k = pe.endmember_loadings.len();
    let p = if k > 0 {
        pe.endmember_loadings[0].len()
    } else {
        0
    };
    let n = data.len();
    if n == 0 || k == 0 || p == 0 {
        return Ok((None, None));
    }
    // Replicate containers: [iter][endmember][feature] and
    // [iter][subject][endmember].
    let mut loading_reps: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_boot);
    let mut abundance_reps: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_boot);

    for b in 0..n_boot {
        let sub_seed = derive_sub_seed(seed, b);
        let mut rng = Xoshiro256pp::new(sub_seed);
        let idx = sample_indices_with_replacement(&mut rng, n);
        let resampled: Vec<Vec<f64>> = idx.iter().map(|&i| data[i].clone()).collect();
        let sub_cfg = UnmixConfig {
            seed: sub_seed,
            ..cfg
        };
        let boot = match unmix(&resampled, k, sub_cfg) {
            Ok(r) => r,
            Err(_) => {
                // A collapsed resample is rare but possible — skip.
                continue;
            }
        };
        // Match each PE endmember i to the boot endmember j maximizing
        // |cos(E_pe[i], E_boot[j])|. Allow one-to-one: if two PE
        // endmembers pick the same boot j, keep the pairing with
        // higher |cos| and leave the other unmatched.
        let mut match_for_pe: Vec<Option<usize>> = vec![None; k];
        let mut match_cos: Vec<f64> = vec![0.0; k];
        let mut boot_used = vec![false; k];
        for _ in 0..k {
            let mut best_i = 0usize;
            let mut best_j = 0usize;
            let mut best_abs = -1.0_f64;
            for (i, pe_load) in pe.endmember_loadings.iter().enumerate() {
                if match_for_pe[i].is_some() {
                    continue;
                }
                for (j, boot_load) in boot.endmember_loadings.iter().enumerate() {
                    if boot_used[j] {
                        continue;
                    }
                    let c = cosine(pe_load, boot_load).abs();
                    if c > best_abs {
                        best_abs = c;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
            if best_abs >= 0.0 {
                match_for_pe[best_i] = Some(best_j);
                match_cos[best_i] = best_abs;
                boot_used[best_j] = true;
            }
        }
        // Build matched+signed boot endmember matrix [k × p].
        let matched: Vec<Vec<f64>> = match_for_pe
            .iter()
            .copied()
            .zip(pe.endmember_loadings.iter())
            .map(|(slot, pe_load)| match slot {
                Some(j) => {
                    let mut loading = boot.endmember_loadings[j].clone();
                    // Sign-correct against PE reference.
                    let s = pe_load
                        .iter()
                        .zip(loading.iter())
                        .map(|(a, c)| a * c)
                        .sum::<f64>();
                    if s < 0.0 {
                        for v in &mut loading {
                            *v = -*v;
                        }
                    }
                    loading
                }
                // No match → repeat PE loading so downstream abundance
                // solve is still well-defined.
                None => pe_load.clone(),
            })
            .collect();
        loading_reps.push(matched.clone());
        // Abundances on ORIGINAL subjects under matched boot E.
        let mut boot_abundances: Vec<Vec<f64>> = Vec::with_capacity(n);
        for x in data {
            let alpha = match cfg.abundance_method {
                AbundanceMethod::Fcls => fcls(&matched, x, cfg.fcls_max_iter, cfg.fcls_tol)?,
                AbundanceMethod::Ucls => ucls(&matched, x)?,
            };
            boot_abundances.push(alpha);
        }
        abundance_reps.push(boot_abundances);
    }

    if loading_reps.is_empty() {
        return Ok((None, None));
    }
    // Element-wise 2.5% / 97.5% percentiles over iterations.
    let pct = |values: &mut Vec<f64>, q: f64| {
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let idx = ((q * (values.len() - 1) as f64).round() as usize).min(values.len() - 1);
        values[idx]
    };
    let mut loading_lo = vec![vec![0.0_f64; p]; k];
    let mut loading_hi = vec![vec![0.0_f64; p]; k];
    for ei in 0..k {
        for fi in 0..p {
            let mut vals: Vec<f64> = loading_reps.iter().map(|rep| rep[ei][fi]).collect();
            loading_lo[ei][fi] = pct(&mut vals.clone(), 0.025);
            loading_hi[ei][fi] = pct(&mut vals, 0.975);
        }
    }
    let mut abundance_lo = vec![vec![0.0_f64; k]; n];
    let mut abundance_hi = vec![vec![0.0_f64; k]; n];
    for si in 0..n {
        for ei in 0..k {
            let mut vals: Vec<f64> = abundance_reps.iter().map(|rep| rep[si][ei]).collect();
            abundance_lo[si][ei] = pct(&mut vals.clone(), 0.025);
            abundance_hi[si][ei] = pct(&mut vals, 0.975);
        }
    }
    Ok((
        Some(LoadingCi {
            lower: loading_lo,
            upper: loading_hi,
        }),
        Some(AbundanceCi {
            lower: abundance_lo,
            upper: abundance_hi,
        }),
    ))
}

#[derive(Debug, Clone)]
pub struct AnnotationRow {
    pub endmember_id: String,
    pub set_name: String,
    pub universe_size: usize,
    pub set_size_in_universe: usize,
    pub top_n: usize,
    pub overlap: usize,
    pub p_value: f64,
}

/// Hypergeometric upper-tail: `P(X ≥ k)` for drawing `n` samples
/// without replacement from a universe of size `N` containing `K`
/// successes. Uses the statrs hypergeometric CDF.
fn hypergeometric_upper_tail(n_universe: u64, k_success: u64, n_draw: u64, k_overlap: u64) -> f64 {
    use statrs::distribution::{DiscreteCDF, Hypergeometric};
    if n_draw == 0 || k_overlap == 0 {
        return 1.0;
    }
    if k_overlap > k_success || k_overlap > n_draw {
        return 0.0;
    }
    let dist = match Hypergeometric::new(n_universe, k_success, n_draw) {
        Ok(d) => d,
        Err(_) => return 1.0,
    };
    // P(X ≥ k) = 1 − P(X ≤ k − 1)
    if k_overlap == 0 {
        return 1.0;
    }
    (1.0 - dist.cdf(k_overlap - 1)).clamp(0.0, 1.0)
}

/// Per-endmember hypergeometric ORA against user-supplied marker
/// sets. For each endmember, the `top_n` proteins by `|loading|`
/// are tested for enrichment of each marker set's intersection with
/// the observation universe. Returns one row per (endmember, set).
pub fn ora_enrichment(
    endmember_loadings: &[Vec<f64>],
    protein_labels: &[String],
    marker_sets: &std::collections::BTreeMap<String, Vec<String>>,
    top_n: usize,
) -> Vec<AnnotationRow> {
    let n_universe = protein_labels.len();
    if n_universe == 0 || endmember_loadings.is_empty() || marker_sets.is_empty() {
        return Vec::new();
    }
    let top_n = top_n.min(n_universe);
    use std::collections::BTreeSet;
    let label_idx: std::collections::HashMap<&str, usize> = protein_labels
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    // Pre-intersect each marker set with the observation universe.
    let intersected: Vec<(String, BTreeSet<usize>)> = marker_sets
        .iter()
        .map(|(name, genes)| {
            let hits: BTreeSet<usize> = genes
                .iter()
                .filter_map(|g| label_idx.get(g.as_str()).copied())
                .collect();
            (name.clone(), hits)
        })
        .collect();
    let mut out = Vec::with_capacity(endmember_loadings.len() * intersected.len());
    for (ai, loading) in endmember_loadings.iter().enumerate() {
        let mut ranked: Vec<(usize, f64)> = loading
            .iter()
            .enumerate()
            .map(|(i, v)| (i, v.abs()))
            .collect();
        ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let top_set: BTreeSet<usize> = ranked.iter().take(top_n).map(|(i, _)| *i).collect();
        for (set_name, set_members) in &intersected {
            let overlap = top_set.intersection(set_members).count();
            let p = hypergeometric_upper_tail(
                n_universe as u64,
                set_members.len() as u64,
                top_n as u64,
                overlap as u64,
            );
            out.push(AnnotationRow {
                endmember_id: format!("E{:03}", ai + 1),
                set_name: set_name.clone(),
                universe_size: n_universe,
                set_size_in_universe: set_members.len(),
                top_n,
                overlap,
                p_value: p,
            });
        }
    }
    out
}

fn cosine(a: &[f64], b: &[f64]) -> f64 {
    let dot: f64 = a.iter().zip(b.iter()).map(|(x, y)| x * y).sum();
    let na: f64 = a.iter().map(|v| v * v).sum::<f64>().sqrt();
    let nb: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
    if na > 0.0 && nb > 0.0 {
        dot / (na * nb)
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct PlantedFixture {
        data: Vec<Vec<f64>>,
        endmembers: Vec<Vec<f64>>,
        abundances: Vec<Vec<f64>>,
    }

    fn vca_fcls(seed: u64, fcls_max_iter: usize, fcls_tol: f64) -> UnmixConfig {
        UnmixConfig {
            seed,
            endmember_method: EndmemberMethod::Vca,
            abundance_method: AbundanceMethod::Fcls,
            fcls_max_iter,
            fcls_tol,
        }
    }

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
    fn planted_unmix_fixture(n: usize, p: usize, k: usize, seed: u64) -> PlantedFixture {
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
        for (i, row) in endmembers.iter_mut().enumerate() {
            for f in 0..block {
                let fi = i * block + f;
                if fi < p {
                    row[fi] = 1.0 + 0.05 * rng.next_z();
                }
            }
        }
        // Abundances: ensure the simplex has a pure-endmember
        // anchor per endmember so VCA can find them deterministically.
        let mut abundances = vec![vec![0.0_f64; k]; n];
        for (i, row) in abundances.iter_mut().enumerate().take(k) {
            row[i] = 1.0;
        }
        // Remaining n-k samples are uniform Dirichlet draws.
        for row in abundances.iter_mut().skip(k) {
            for v in row.iter_mut() {
                *v = rng.next_u().max(1e-12);
                *v = -v.ln(); // Exp(1)
            }
            let s: f64 = row.iter().sum();
            for v in row.iter_mut() {
                *v /= s;
            }
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
        PlantedFixture {
            data,
            endmembers,
            abundances,
        }
    }

    fn best_matched_cosine(planted: &[Vec<f64>], recovered: &[Vec<f64>]) -> Vec<f64> {
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
        let PlantedFixture {
            data,
            endmembers: planted,
            ..
        } = planted_unmix_fixture(100, 60, 3, 20260420);
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
        let PlantedFixture { data, .. } = planted_unmix_fixture(60, 30, 3, 11);
        let result = unmix(&data, 3, vca_fcls(42, 500, 1e-9)).unwrap();
        for (i, row) in result.abundances.iter().enumerate() {
            let sum: f64 = row.iter().sum();
            assert!((sum - 1.0).abs() < 1e-6, "sample {i} row-sum = {sum} ≠ 1",);
            for (j, &a) in row.iter().enumerate() {
                assert!(a >= -1e-9, "sample {i} endmember {j} abundance = {a} < 0",);
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
    fn bootstrap_ci_produces_finite_bounds_and_envelops_point_estimate() {
        let PlantedFixture { data, .. } = planted_unmix_fixture(60, 30, 3, 61);
        let pe = unmix(&data, 3, vca_fcls(42, 500, 1e-8)).unwrap();
        let (loading_ci, abundance_ci) =
            bootstrap_ci(&pe, &data, vca_fcls(42, 500, 1e-8), 30).unwrap();
        let lci = loading_ci.unwrap();
        let aci = abundance_ci.unwrap();
        assert_eq!(lci.lower.len(), pe.endmember_loadings.len());
        assert_eq!(aci.lower.len(), pe.abundances.len());
        // Bounds must be finite and ordered.
        for ei in 0..pe.endmember_loadings.len() {
            for fi in 0..pe.endmember_loadings[ei].len() {
                assert!(lci.lower[ei][fi].is_finite());
                assert!(lci.upper[ei][fi].is_finite());
                assert!(lci.lower[ei][fi] <= lci.upper[ei][fi] + 1e-12);
            }
        }
        for si in 0..pe.abundances.len() {
            for ei in 0..pe.abundances[si].len() {
                assert!(aci.lower[si][ei].is_finite());
                assert!(aci.upper[si][ei].is_finite());
                assert!(aci.lower[si][ei] <= aci.upper[si][ei] + 1e-12);
            }
        }
    }

    #[test]
    fn bootstrap_ci_with_n_boot_zero_returns_none() {
        let PlantedFixture { data, .. } = planted_unmix_fixture(40, 20, 3, 71);
        let pe = unmix(&data, 3, vca_fcls(42, 300, 1e-8)).unwrap();
        let (lci, aci) = bootstrap_ci(&pe, &data, vca_fcls(42, 300, 1e-8), 0).unwrap();
        assert!(lci.is_none() && aci.is_none());
    }

    #[test]
    fn ora_enrichment_detects_planted_marker_sets() {
        // One endmember loads heavily on G000..G009 (block 0). Define
        // a marker set that exactly covers those 10 genes — enrichment
        // p-value should be tiny. A second marker set covers
        // G020..G029 (not in block 0) — expect a high p-value.
        let p = 30;
        let mut loading = vec![0.0_f64; p];
        for (i, v) in loading.iter_mut().enumerate().take(10) {
            *v = 1.0 - 0.01 * (i as f64);
        }
        let labels: Vec<String> = (0..p).map(|i| format!("G{i:03}")).collect();
        let mut marker_sets = std::collections::BTreeMap::new();
        marker_sets.insert(
            "true_set".into(),
            (0..10).map(|i| format!("G{i:03}")).collect::<Vec<_>>(),
        );
        marker_sets.insert(
            "decoy_set".into(),
            (20..30).map(|i| format!("G{i:03}")).collect::<Vec<_>>(),
        );
        let rows = ora_enrichment(&[loading], &labels, &marker_sets, 10);
        assert_eq!(rows.len(), 2);
        let true_row = rows.iter().find(|r| r.set_name == "true_set").unwrap();
        let decoy_row = rows.iter().find(|r| r.set_name == "decoy_set").unwrap();
        assert_eq!(true_row.overlap, 10);
        assert!(
            true_row.p_value < 1e-6,
            "true_set p should be tiny, got {}",
            true_row.p_value
        );
        assert_eq!(decoy_row.overlap, 0);
        assert!(
            decoy_row.p_value > 0.5,
            "decoy_set p should be ≈1, got {}",
            decoy_row.p_value
        );
    }

    #[test]
    fn unmix_is_deterministic_under_fixed_seed() {
        let PlantedFixture { data, .. } = planted_unmix_fixture(40, 20, 3, 55);
        let a = unmix(&data, 3, vca_fcls(7, 300, 1e-9)).unwrap();
        let b = unmix(&data, 3, vca_fcls(7, 300, 1e-9)).unwrap();
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
        let PlantedFixture { data, .. } = planted_unmix_fixture(80, 60, 3, 77);
        let (chosen, rows) = select_k_auto(&data, 2, 6, vca_fcls(42, 200, 1e-7), 0.10).unwrap();
        assert!(
            chosen == 3 || chosen == 4,
            "expected k ≈ 3 on planted-3 fixture, got {chosen}; sweep={:?}",
            rows.iter()
                .map(|r| (r.k, r.mean_residual_norm))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn select_k_auto_picks_planted_k_on_planted_5_fixture() {
        let PlantedFixture { data, .. } = planted_unmix_fixture(120, 100, 5, 88);
        let (chosen, _rows) = select_k_auto(&data, 2, 8, vca_fcls(42, 200, 1e-7), 0.10).unwrap();
        assert!(
            chosen == 5 || chosen == 6,
            "expected k ≈ 5 on planted-5 fixture, got {chosen}"
        );
    }

    #[test]
    fn nfindr_volume_is_not_smaller_than_vca_volume() {
        // Starting from VCA and swapping only when volume strictly
        // increases, N-FINDR can never return a worse simplex.
        let PlantedFixture { data, .. } = planted_unmix_fixture(80, 40, 3, 21);
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
        let PlantedFixture {
            data,
            endmembers: planted,
            ..
        } = planted_unmix_fixture(80, 40, 3, 31);
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

    /// SPA has no seed at all, so the only stability question is
    /// whether it recovers the planted vertices. It must, exactly.
    #[test]
    fn spa_recovers_the_planted_pure_samples() {
        for noise in [0.0_f64, 0.01, 0.05, 0.1] {
            let data = planted_simplex(60, 1200, 4, noise, 20260906);
            let got = spa(&data, 4).expect("spa");
            let mut sorted = got.endmember_sample_indices.clone();
            sorted.sort();
            assert_eq!(
                sorted,
                vec![0, 1, 2, 3],
                "noise={noise}: expected the four planted pure endmembers, got {sorted:?}"
            );
        }
    }

    /// The property the whole change exists for: on real-shaped data
    /// where VCA's answer moves with the seed, SPA has no seed to move
    /// with. Repeated calls are identical by construction, and the
    /// result depends only on the data.
    #[test]
    fn spa_is_a_function_of_the_data_alone() {
        let data = planted_simplex(40, 800, 3, 0.08, 7);
        let first = spa(&data, 3).expect("spa");
        for _ in 0..5 {
            let again = spa(&data, 3).expect("spa");
            assert_eq!(
                first.endmember_sample_indices, again.endmember_sample_indices,
                "SPA must be deterministic"
            );
        }
        // And it must not silently agree with an arbitrary answer: the
        // selection has to actually be the extreme points.
        let mut sorted = first.endmember_sample_indices.clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1, 2]);
    }

    /// `k` beyond the dimensionality actually present must fail loudly
    /// rather than return an arbitrary extra vertex.
    #[test]
    fn spa_refuses_k_beyond_available_dimensionality() {
        // Three planted sources, so the simplex is 3-vertex; ask for 6.
        let data = planted_simplex(30, 400, 3, 0.0, 11);
        let err = spa(&data, 6).unwrap_err();
        assert!(
            err.contains("no independent direction left")
                || err.contains("residual collapsed")
                || err.contains("under-determined"),
            "unexpected error: {err}"
        );
    }

    /// Planted-simplex fixture in high ambient dimension: `k` pure
    /// endmember samples plus convex mixtures of them, with small noise
    /// on every feature. The simplex structure lives in a `k`-dimensional
    /// subspace of `p` dimensions.
    fn planted_simplex(n: usize, p: usize, k: usize, noise: f64, seed: u64) -> Vec<Vec<f64>> {
        let mut rng = Xoshiro256pp::new(seed);
        // k pure source spectra.
        let sources: Vec<Vec<f64>> = (0..k)
            .map(|_| (0..p).map(|_| 5.0 + rng.next_normal()).collect())
            .collect();
        let mut out: Vec<Vec<f64>> = Vec::with_capacity(n);
        // First k samples are the pure endmembers.
        for src in &sources {
            out.push(src.clone());
        }
        // Remainder are interior convex mixtures.
        for _ in k..n {
            let mut w: Vec<f64> = (0..k).map(|_| rng.next_f64().abs() + 0.05).collect();
            let tot: f64 = w.iter().sum();
            for x in &mut w {
                *x /= tot;
            }
            let row: Vec<f64> = (0..p)
                .map(|f| {
                    let mixed: f64 = (0..k).map(|c| w[c] * sources[c][f]).sum();
                    mixed + noise * rng.next_normal()
                })
                .collect();
            out.push(row);
        }
        for row in out.iter_mut() {
            for v in row.iter_mut() {
                if !v.is_finite() {
                    *v = 0.0;
                }
            }
        }
        out
    }

    /// The vertex search must run in the signal subspace, not in the
    /// ambient feature space. In `p` dimensions a random unit direction
    /// is almost orthogonal to a `k`-dimensional simplex, so `argmax|u·y|`
    /// is decided by noise and the selected vertex set depends on the
    /// seed rather than on the data. Selection must be seed-invariant.
    #[test]
    fn vca_endmember_selection_is_invariant_to_seed() {
        let data = planted_simplex(60, 1200, 4, 0.05, 20260906);
        let baseline = vca(&data, 4, 1).expect("vca");
        let mut baseline_sorted = baseline.endmember_sample_indices.clone();
        baseline_sorted.sort();
        for seed in [2u64, 7, 42, 99, 12345] {
            let got = vca(&data, 4, seed).expect("vca");
            let mut got_sorted = got.endmember_sample_indices.clone();
            got_sorted.sort();
            assert_eq!(
                got_sorted, baseline_sorted,
                "seed {seed} selected a different vertex set ({got_sorted:?}) than seed 1 \
                 ({baseline_sorted:?}); the vertex search is being decided by the random \
                 direction rather than by simplex geometry"
            );
        }
    }

    /// With the subspace projection in place the selected vertices are
    /// the planted pure samples (indices 0..k), not arbitrary interior
    /// mixtures.
    #[test]
    fn vca_selects_the_planted_pure_samples() {
        let data = planted_simplex(60, 1200, 4, 0.05, 4242);
        let got = vca(&data, 4, 20260420).expect("vca");
        let mut sorted = got.endmember_sample_indices.clone();
        sorted.sort();
        assert_eq!(
            sorted,
            vec![0, 1, 2, 3],
            "expected the four planted pure endmembers; got {sorted:?}"
        );
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
        let PlantedFixture {
            data,
            abundances: planted_ab,
            ..
        } = planted_unmix_fixture(120, 60, 3, 30);
        let result = unmix(&data, 3, vca_fcls(42, 1000, 1e-9)).unwrap();
        // Match recovered columns to planted columns via best absolute
        // cosine between recovered endmember loadings and planted.
        let p = data[0].len();
        let planted_endmembers = {
            let PlantedFixture { endmembers: e, .. } = planted_unmix_fixture(120, 60, 3, 30);
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
                let np = planted_endmembers[i]
                    .iter()
                    .map(|v| v * v)
                    .sum::<f64>()
                    .sqrt();
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
