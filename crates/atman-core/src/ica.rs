//! FastICA with symmetric decorrelation and multi-seed stability primitives.
//!
//! The algorithm follows the standard FastICA formulation with the
//! `log(cosh)` contrast function: center, PCA-whiten to `k` components,
//! iterate `W <- (W W^T)^{-1/2} W` after fixed-point updates
//! `w_i <- E[x g(w_i^T x)] - E[g'(w_i^T x)] w_i`.
//!
//! Everything is pure Rust and deterministic: the built-in Xoshiro256++
//! PRNG means each `(k, seed)` pair yields a byte-exact decomposition.

use std::cmp::Ordering;

pub use crate::stats::jaccard_top_n;

// The reviewed Xoshiro256++ now lives in `crate::rng`; re-export it here so
// the historical `crate::ica::Xoshiro256pp` path keeps working for the core
// decomposition family (ica_null, align_bootstrap, nmf, decompose_unmix,
// gsea) without churning their imports.
pub use crate::rng::Xoshiro256pp;

/// Eigendecomposition of a symmetric `n x n` matrix via cyclic Jacobi rotations.
/// Returns `(values, vectors)` with columns of `vectors` as eigenvectors,
/// sorted in descending eigenvalue order.
pub fn jacobi_eigen(matrix: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = matrix.len();
    assert!(n > 0 && matrix.iter().all(|row| row.len() == n));
    let mut a = matrix.to_vec();
    let mut v: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();

    let max_sweeps = 200;
    for _sweep in 0..max_sweeps {
        let off: f64 = a
            .iter()
            .enumerate()
            .flat_map(|(p, row)| row.iter().skip(p + 1).map(|v| v.abs()))
            .sum();
        if off <= 1e-14 * (n as f64) {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p][q];
                if apq.abs() < 1e-18 {
                    continue;
                }
                let app = a[p][p];
                let aqq = a[q][q];
                let theta = (aqq - app) / (2.0 * apq);
                let t = if theta.abs() > 1e15 {
                    1.0 / (2.0 * theta)
                } else {
                    let sign = if theta >= 0.0 { 1.0 } else { -1.0 };
                    sign / (theta.abs() + (1.0 + theta * theta).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                a[p][p] = app - t * apq;
                a[q][q] = aqq + t * apq;
                a[p][q] = 0.0;
                a[q][p] = 0.0;
                for r in 0..n {
                    if r != p && r != q {
                        let arp = a[r][p];
                        let arq = a[r][q];
                        a[r][p] = c * arp - s * arq;
                        a[p][r] = a[r][p];
                        a[r][q] = s * arp + c * arq;
                        a[q][r] = a[r][q];
                    }
                    let vrp = v[r][p];
                    let vrq = v[r][q];
                    v[r][p] = c * vrp - s * vrq;
                    v[r][q] = s * vrp + c * vrq;
                }
            }
        }
    }

    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&i, &j| a[j][j].partial_cmp(&a[i][i]).unwrap_or(Ordering::Equal));
    let values: Vec<f64> = indices.iter().map(|&i| a[i][i]).collect();
    let vectors: Vec<Vec<f64>> = (0..n)
        .map(|r| indices.iter().map(|&i| v[r][i]).collect())
        .collect();
    (values, vectors)
}

/// PCA whitening result: `whitened = X_centered * whitening^T` has unit-variance
/// columns (approx). `unwhitening` maps sources back to the centered feature space.
pub struct Whitening {
    /// Column means of the input matrix.
    pub mean: Vec<f64>,
    /// `k x p` whitening matrix.
    pub whitening: Vec<Vec<f64>>,
    /// `p x k` pseudo-inverse of `whitening` for mapping sources back to features.
    pub unwhitening: Vec<Vec<f64>>,
    /// Top `k` eigenvalues of the sample covariance.
    pub eigenvalues: Vec<f64>,
    /// Full eigenvalue spectrum (sorted descending). Useful for K selection.
    pub full_spectrum: Vec<f64>,
    /// `n x k` whitened data.
    pub whitened: Vec<Vec<f64>>,
}

/// Standard PCA-based whitening for `n x p` matrix `x` to `k` components.
/// Reliability-weighted whitening: the same decomposition as
/// [`pca_whiten`], with each cell's contribution to the mean and the
/// covariance scaled by `cell_weights[i][j]`.
///
/// This is where a per-cell detection model belongs. Weighting a
/// FastICA contrast function per cell has no clean interpretation,
/// because the contrast acts on the projection `w·x` and a single cell
/// has no separate identity there — which is why the previous attempt
/// collapsed per-cell weights to a per-sample scalar and downweighted
/// whole low-abundance samples, losing heavy-tailed source recovery.
/// The covariance is different: it is a sum over cells, so a cell that
/// is probably a non-detection can contribute less to the second-order
/// structure without any sample being discarded.
///
/// Weighted moments:
/// `mu_a = Σ_i w_ia x_ia / Σ_i w_ia`, and
/// `C_ab = Σ_i w_ia w_ib (x_ia − mu_a)(x_ib − mu_b) / Σ_i w_ia w_ib`.
///
/// Uniform weights reproduce [`pca_whiten`] up to floating point.
pub fn pca_whiten_weighted(x: &[Vec<f64>], k: usize, cell_weights: &[Vec<f64>]) -> Whitening {
    let n = x.len();
    assert!(n > 1 && k > 0);
    let p = x[0].len();
    assert!(x.iter().all(|r| r.len() == p));
    assert!(
        cell_weights.len() == n && cell_weights.iter().all(|r| r.len() == p),
        "cell_weights must be n x p"
    );

    // Weighted column means.
    let mut mean = vec![0.0_f64; p];
    for j in 0..p {
        let mut num = 0.0;
        let mut den = 0.0;
        for i in 0..n {
            let w = cell_weights[i][j].max(0.0);
            num += w * x[i][j];
            den += w;
        }
        mean[j] = if den > 0.0 { num / den } else { 0.0 };
    }
    let xc: Vec<Vec<f64>> = x
        .iter()
        .map(|row| row.iter().enumerate().map(|(j, v)| v - mean[j]).collect())
        .collect();

    // Weighted covariance in feature space. Always p x p here: the
    // sample-space Gram shortcut does not carry per-cell weights,
    // because a weighted inner product between two SAMPLES is not the
    // same object as the weighted covariance between two FEATURES.
    let mut cov = vec![vec![0.0_f64; p]; p];
    for a in 0..p {
        for b in a..p {
            let mut num = 0.0;
            let mut sum_w = 0.0;
            let mut sum_w2 = 0.0;
            for i in 0..n {
                let w = cell_weights[i][a].max(0.0) * cell_weights[i][b].max(0.0);
                num += w * xc[i][a] * xc[i][b];
                sum_w += w;
                sum_w2 += w * w;
            }
            // Unbiased denominator for reliability weights:
            // `Σw − Σw²/Σw`, which is exactly `n − 1` when every weight
            // is 1. Using `Σw` instead would make uniform weights
            // disagree with `pca_whiten` by `n/(n−1)`, and any
            // comparison between the two paths would then be measuring
            // that instead of the weighting.
            let den = if sum_w > 0.0 {
                sum_w - sum_w2 / sum_w
            } else {
                0.0
            };
            let v = if den > 0.0 { num / den } else { 0.0 };
            cov[a][b] = v;
            cov[b][a] = v;
        }
    }
    let (vals, vecs) = jacobi_eigen(&cov);
    let mut whitening_pxk = vec![vec![0.0_f64; k]; p];
    for (i, row) in whitening_pxk.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = vecs[i][j];
        }
    }
    finish_whitening(xc, mean, vals, whitening_pxk, k, p)
}

pub fn pca_whiten(x: &[Vec<f64>], k: usize) -> Whitening {
    let n = x.len();
    assert!(n > 1 && k > 0);
    let p = x[0].len();
    assert!(x.iter().all(|r| r.len() == p));

    let mut mean = vec![0.0_f64; p];
    for row in x {
        for (j, v) in row.iter().enumerate() {
            mean[j] += *v;
        }
    }
    for m in mean.iter_mut() {
        *m /= n as f64;
    }

    let xc: Vec<Vec<f64>> = x
        .iter()
        .map(|row| row.iter().enumerate().map(|(j, v)| v - mean[j]).collect())
        .collect();

    // Work on the smaller of the two Gram matrices: X^T X (p x p) or X X^T (n x n).
    let (eigs, whitening_pxk) = if p <= n {
        let mut cov = vec![vec![0.0_f64; p]; p];
        for row in &xc {
            for i in 0..p {
                let xi = row[i];
                for j in 0..p {
                    cov[i][j] += xi * row[j];
                }
            }
        }
        let scale = (n as f64 - 1.0).max(1.0);
        for row in cov.iter_mut() {
            for v in row.iter_mut() {
                *v /= scale;
            }
        }
        let (vals, vecs) = jacobi_eigen(&cov);
        // `vecs` is p x p; take first k columns.
        let mut whitening_pxk = vec![vec![0.0_f64; k]; p];
        for i in 0..p {
            for j in 0..k {
                whitening_pxk[i][j] = vecs[i][j];
            }
        }
        (vals, whitening_pxk)
    } else {
        // Gram matrix X X^T (n x n); convert eigenvectors back to feature space.
        let mut gram = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            for j in i..n {
                let acc: f64 = xc[i].iter().zip(xc[j].iter()).map(|(a, b)| a * b).sum();
                gram[i][j] = acc;
                gram[j][i] = acc;
            }
        }
        let scale = (n as f64 - 1.0).max(1.0);
        for row in gram.iter_mut() {
            for v in row.iter_mut() {
                *v /= scale;
            }
        }
        let (vals, vecs_n) = jacobi_eigen(&gram);
        // Feature-space eigenvector V = X^T U / sqrt((n-1) * lambda).
        let mut vecs_p = vec![vec![0.0_f64; k]; p];
        for j in 0..k {
            let lam = vals[j].max(1e-18);
            let denom = (scale * lam).sqrt();
            for i in 0..p {
                let mut acc = 0.0_f64;
                for t in 0..n {
                    acc += xc[t][i] * vecs_n[t][j];
                }
                vecs_p[i][j] = acc / denom;
            }
        }
        (vals, vecs_p)
    };

    finish_whitening(xc, mean, eigs, whitening_pxk, k, p)
}

/// Shared tail of the whitening constructions: build the whitening and
/// unwhitening matrices from the eigenpairs and project the centred
/// data. Split out so the plain and reliability-weighted paths cannot
/// drift apart.
fn finish_whitening(
    xc: Vec<Vec<f64>>,
    mean: Vec<f64>,
    eigs: Vec<f64>,
    whitening_pxk: Vec<Vec<f64>>,
    k: usize,
    p: usize,
) -> Whitening {
    let n = xc.len();
    let top_eigs: Vec<f64> = eigs.iter().take(k).copied().collect();
    // whitening = diag(1/sqrt(lambda)) * V^T  (k x p)
    let mut whitening = vec![vec![0.0_f64; p]; k];
    for j in 0..k {
        let scale = 1.0 / top_eigs[j].max(1e-18).sqrt();
        for i in 0..p {
            whitening[j][i] = whitening_pxk[i][j] * scale;
        }
    }
    // unwhitening = V * diag(sqrt(lambda))  (p x k)
    let mut unwhitening = vec![vec![0.0_f64; k]; p];
    for j in 0..k {
        let scale = top_eigs[j].max(1e-18).sqrt();
        for i in 0..p {
            unwhitening[i][j] = whitening_pxk[i][j] * scale;
        }
    }
    let mut whitened = vec![vec![0.0_f64; k]; n];
    for (i, row) in xc.iter().enumerate() {
        for j in 0..k {
            let mut acc = 0.0_f64;
            for t in 0..p {
                acc += whitening[j][t] * row[t];
            }
            whitened[i][j] = acc;
        }
    }
    Whitening {
        mean,
        whitening,
        unwhitening,
        eigenvalues: top_eigs,
        full_spectrum: eigs,
        whitened,
    }
}

/// Output of one FastICA run.
pub struct IcaResult {
    /// `k x p` unmixing matrix; `sources = (x - mean) * unmixing^T`.
    pub unmixing: Vec<Vec<f64>>,
    /// `p x k` mixing matrix; `x - mean ≈ sources * mixing^T`.
    pub mixing: Vec<Vec<f64>>,
    /// `n x k` source activations per sample.
    pub sources: Vec<Vec<f64>>,
    /// Column means used for centering.
    pub mean: Vec<f64>,
    /// Iterations until convergence (or `max_iter`).
    pub n_iterations: usize,
    /// Final convergence metric (max absolute deviation from identity of W*W_old^T).
    pub final_tol: f64,
    /// Whitening eigenvalues (top `k`).
    pub eigenvalues: Vec<f64>,
    /// Full eigenvalue spectrum.
    pub full_spectrum: Vec<f64>,
}

/// Run FastICA on `x` (n x p) for `k` components with seed, max iterations, and tolerance.
pub fn fast_ica(x: &[Vec<f64>], k: usize, seed: u64, max_iter: usize, tol: f64) -> IcaResult {
    fast_ica_weighted(x, k, seed, max_iter, tol, &[])
}

/// Run weighted FastICA on `x` (n x p) for `k` components.
///
/// `cell_weights` is a flat row-major slice of per-cell weights with shape n × p
/// (i.e., `cell_weights[i * p + j]` is the weight for sample `i`, feature `j`).
/// An empty slice (`&[]`) is equivalent to all-ones weights, which collapses to
/// standard FastICA (the MAR-collapse parity contract).
///
/// The weights enter the fixed-point updates via a weighted expectation:
///   `w_c ← (Σ_i w_i * g(w_c^T x̃_i) * x̃_i) / Σ_i w_i_avg
///           - (Σ_i w_i_avg * g'(w_c^T x̃_i)) / Σ_i w_i_avg * w_c`
/// where `w_i_avg` is the mean of the per-cell weights for sample `i` (averaged
/// across features), so the sphering objective uses the same per-sample weight
/// for all whitened dimensions.  When all `cell_weights` equal 1 the update
/// reduces identically to the unweighted `fast_ica` above.
///
/// The whitening step uses the *unweighted* covariance.  Weighted whitening
/// (i.e., incorporating per-cell missingness into the sphering step) is a
/// deferred feature — under the mild missingness regime of Phase G the gain
/// is negligible relative to the weighted fixed-point iteration.
/// FastICA with reliability-weighted whitening.
///
/// `cell_weights` is `n x p` and scales each cell's contribution to the
/// whitening moments. The fixed-point iteration itself runs unweighted
/// on the whitened data, which is the point: a per-cell weight has a
/// clean meaning in a covariance and none in a contrast function.
///
/// Pass an empty slice to fall back to [`fast_ica`].
pub fn fast_ica_whitening_weighted(
    x: &[Vec<f64>],
    k: usize,
    seed: u64,
    max_iter: usize,
    tol: f64,
    cell_weights: &[Vec<f64>],
) -> IcaResult {
    if cell_weights.is_empty() {
        return fast_ica(x, k, seed, max_iter, tol);
    }
    let n = x.len();
    fast_ica_from_whitening_weighted(
        pca_whiten_weighted(x, k, cell_weights),
        k,
        seed,
        max_iter,
        tol,
        &vec![1.0_f64; n],
    )
}

pub fn fast_ica_weighted(
    x: &[Vec<f64>],
    k: usize,
    seed: u64,
    max_iter: usize,
    tol: f64,
    cell_weights: &[f64],
) -> IcaResult {
    let n = x.len();
    let p = if n > 0 { x[0].len() } else { 0 };

    // Validate weight dimensions when non-empty.
    let use_weights = !cell_weights.is_empty();
    if use_weights {
        assert_eq!(
            cell_weights.len(),
            n * p,
            "cell_weights must have length n*p = {} but got {}",
            n * p,
            cell_weights.len()
        );
    }

    // Compute per-sample scalar weight: mean of per-cell weights for that sample.
    // When weights are uniform (or not provided) every sample_weight is 1.0,
    // so the weighted mean equals the unweighted mean and the update is identical.
    let sample_weights: Vec<f64> = if use_weights {
        (0..n)
            .map(|i| {
                let start = i * p;
                let sum: f64 = cell_weights[start..start + p].iter().sum();
                sum / p as f64
            })
            .collect()
    } else {
        vec![1.0_f64; n]
    };

    let w = pca_whiten(x, k);
    fast_ica_from_whitening_weighted(w, k, seed, max_iter, tol, &sample_weights)
}

/// FastICA fixed point on an already-whitened matrix.
///
/// Split out so the plain path, the sample-weighted path and the
/// reliability-weighted-whitening path share one iteration and cannot
/// drift. `sample_weights` weights the CONTRAST; per-cell reliability
/// belongs in the whitening that produced `w`, not here.
fn fast_ica_from_whitening_weighted(
    w: Whitening,
    k: usize,
    seed: u64,
    max_iter: usize,
    tol: f64,
    sample_weights: &[f64],
) -> IcaResult {
    let n = w.whitened.len();
    let p = w.mean.len();
    let total_weight: f64 = sample_weights.iter().sum();
    let mut rng = Xoshiro256pp::new(seed);
    let mut weights = vec![vec![0.0_f64; k]; k];
    for row in weights.iter_mut() {
        for v in row.iter_mut() {
            *v = rng.next_normal();
        }
    }
    symmetric_decorrelate(&mut weights);

    let mut final_tol = f64::INFINITY;
    let mut iter_used = 0;
    for iteration in 0..max_iter {
        iter_used = iteration + 1;
        let mut next = vec![vec![0.0_f64; k]; k];
        // For each component, compute weighted E[x g(w^T x)] - weighted E[g'(w^T x)] w.
        for c in 0..k {
            let w_c: Vec<f64> = weights[c].clone();
            let mut sum_gxx = vec![0.0_f64; k];
            let mut sum_gp = 0.0_f64;
            for (i, row) in w.whitened.iter().enumerate() {
                let sw = sample_weights[i];
                let mut wx = 0.0_f64;
                for (j, xj) in row.iter().enumerate() {
                    wx += w_c[j] * xj;
                }
                let gwx = wx.tanh();
                let gpwx = 1.0 - gwx * gwx;
                sum_gp += sw * gpwx;
                for (j, xj) in row.iter().enumerate() {
                    sum_gxx[j] += sw * gwx * xj;
                }
            }
            let mean_gp = sum_gp / total_weight;
            for j in 0..k {
                next[c][j] = sum_gxx[j] / total_weight - mean_gp * w_c[j];
            }
        }
        symmetric_decorrelate(&mut next);
        let dev = max_identity_deviation(&next, &weights);
        weights = next;
        final_tol = dev;
        if dev < tol {
            break;
        }
    }

    // sources = whitened * W^T (n x k).
    let mut sources = vec![vec![0.0_f64; k]; n];
    for (i, row) in w.whitened.iter().enumerate() {
        for c in 0..k {
            let mut acc = 0.0_f64;
            for j in 0..k {
                acc += row[j] * weights[c][j];
            }
            sources[i][c] = acc;
        }
    }
    // unmixing (k x p) = W * whitening
    let unmixing: Vec<Vec<f64>> = weights
        .iter()
        .map(|w_row| {
            (0..p)
                .map(|j| {
                    w_row
                        .iter()
                        .zip(w.whitening.iter())
                        .map(|(&wt, white_row)| wt * white_row[j])
                        .sum()
                })
                .collect()
        })
        .collect();
    // mixing (p x k) = unwhitening * W^T  (since unwhitening = V * D^{1/2} and W is orthonormal)
    let mixing: Vec<Vec<f64>> = w
        .unwhitening
        .iter()
        .map(|uw_row| {
            weights
                .iter()
                .map(|w_row| uw_row.iter().zip(w_row.iter()).map(|(&u, &v)| u * v).sum())
                .collect()
        })
        .collect();

    IcaResult {
        unmixing,
        mixing,
        sources,
        mean: w.mean,
        n_iterations: iter_used,
        final_tol,
        eigenvalues: w.eigenvalues,
        full_spectrum: w.full_spectrum,
    }
}

/// Number of components required to reach at least `target` cumulative variance,
/// clamped to `[k_min, k_max]`.
pub fn select_k_cumulative_variance(
    spectrum: &[f64],
    target: f64,
    k_min: usize,
    k_max: usize,
) -> usize {
    let total: f64 = spectrum.iter().map(|v| v.max(0.0)).sum();
    if total <= 0.0 {
        return k_min.max(1);
    }
    let mut cum = 0.0_f64;
    let mut chosen = k_min;
    for (idx, lam) in spectrum.iter().enumerate() {
        cum += lam.max(0.0) / total;
        if cum >= target {
            chosen = idx + 1;
            break;
        }
        chosen = idx + 1;
    }
    chosen.clamp(k_min.max(1), k_max.max(k_min.max(1)))
}

/// Symmetric decorrelation: W <- (W W^T)^{-1/2} W, in place.
fn symmetric_decorrelate(w: &mut Vec<Vec<f64>>) {
    let k = w.len();
    let mut wwt = vec![vec![0.0_f64; k]; k];
    for i in 0..k {
        for j in i..k {
            let acc = row_dot(&w[i], &w[j]);
            wwt[i][j] = acc;
            wwt[j][i] = acc;
        }
    }
    let (vals, vecs) = jacobi_eigen(&wwt);
    // (WW^T)^{-1/2} = V diag(1/sqrt(lambda)) V^T
    let inv_sqrt: Vec<Vec<f64>> = vecs
        .iter()
        .map(|vi_row| {
            vecs.iter()
                .map(|vj_row| {
                    vi_row
                        .iter()
                        .zip(vj_row.iter())
                        .zip(vals.iter())
                        .map(|((vi, vj), &lam)| vi * vj / lam.max(1e-18).sqrt())
                        .sum()
                })
                .collect()
        })
        .collect();
    let p = w[0].len();
    let next: Vec<Vec<f64>> = inv_sqrt
        .iter()
        .map(|inv_row| {
            (0..p)
                .map(|j| {
                    inv_row
                        .iter()
                        .zip(w.iter())
                        .map(|(&inv, w_row)| inv * w_row[j])
                        .sum()
                })
                .collect()
        })
        .collect();
    *w = next;
}

fn row_dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn max_identity_deviation(w_new: &[Vec<f64>], w_old: &[Vec<f64>]) -> f64 {
    let mut max_dev = 0.0_f64;
    for (i, new_row) in w_new.iter().enumerate() {
        for (j, old_row) in w_old.iter().enumerate() {
            let acc = row_dot(new_row, old_row);
            let expected = if i == j { 1.0 } else { 0.0 };
            let dev = (acc.abs() - expected).abs();
            if dev > max_dev {
                max_dev = dev;
            }
        }
    }
    max_dev
}

/// Canonicalized ICA result: sign-flipped so each program's
/// max-|loading| entry is positive, and programs sorted by descending
/// max |loading|. Makes multi-seed stability comparisons meaningful
/// under ICA's inherent sign and permutation ambiguity.
pub struct CanonicalIca {
    /// `k` rows, each of length `p` (protein loadings per program).
    pub loadings: Vec<Vec<f64>>,
    /// `k` rows, each of length `n` (per-sample activations).
    pub activations: Vec<Vec<f64>>,
}

/// Canonicalize an [`IcaResult`]: sign-flip each program so its
/// max |loading| entry is positive, then sort programs by descending
/// max |loading|. The result is stable across runs that produce the
/// same components modulo ICA's sign / permutation gauge.
pub fn canonicalize_ica(result: &IcaResult) -> CanonicalIca {
    let k = result.mixing[0].len();
    let p = result.mixing.len();

    let mut signs = vec![1.0_f64; k];
    let mut loadings: Vec<Vec<f64>> = vec![vec![0.0; p]; k];
    for (c, sign) in signs.iter_mut().enumerate() {
        let mut max_abs = 0.0_f64;
        let mut argmax = 0usize;
        for (j, row) in result.mixing.iter().enumerate() {
            if row[c].abs() > max_abs {
                max_abs = row[c].abs();
                argmax = j;
            }
        }
        if result.mixing[argmax][c] < 0.0 {
            *sign = -1.0;
        }
        for (j, slot) in loadings[c].iter_mut().enumerate() {
            *slot = *sign * result.mixing[j][c];
        }
    }
    let activations: Vec<Vec<f64>> = signs
        .iter()
        .enumerate()
        .map(|(c, &sign)| {
            result
                .sources
                .iter()
                .map(|src_row| sign * src_row[c])
                .collect()
        })
        .collect();
    let mut order: Vec<usize> = (0..k).collect();
    order.sort_by(|&a, &b| {
        let ma = loadings[a].iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
        let mb = loadings[b].iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
        mb.partial_cmp(&ma).unwrap_or(std::cmp::Ordering::Equal)
    });
    let ordered_loadings: Vec<Vec<f64>> = order.iter().map(|&i| loadings[i].clone()).collect();
    let ordered_activations: Vec<Vec<f64>> =
        order.iter().map(|&i| activations[i].clone()).collect();
    CanonicalIca {
        loadings: ordered_loadings,
        activations: ordered_activations,
    }
}

/// Per-reference-program mean of the best Jaccard top-N overlap
/// achieved by any alternative-seed program. Mirrors the stability
/// metric written to `stability.tsv` by `atman decompose ica` and
/// used as the "observed stability" in `atman decompose null`.
pub fn compute_stability_scores(
    reference: &CanonicalIca,
    alternatives: &[CanonicalIca],
    top_n: usize,
) -> Vec<f64> {
    if alternatives.is_empty() {
        return vec![f64::NAN; reference.loadings.len()];
    }
    reference
        .loadings
        .iter()
        .map(|ref_row| {
            let sum_best: f64 = alternatives
                .iter()
                .map(|alt| {
                    alt.loadings
                        .iter()
                        .map(|alt_row| jaccard_top_n(ref_row, alt_row, top_n))
                        .fold(0.0_f64, f64::max)
                })
                .sum();
            sum_best / alternatives.len() as f64
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xoshiro_is_deterministic() {
        let mut a = Xoshiro256pp::new(20260418);
        let mut b = Xoshiro256pp::new(20260418);
        for _ in 0..128 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn jacobi_recovers_eigenvalues_of_diagonal_matrix() {
        let m = vec![
            vec![3.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 2.0],
        ];
        let (vals, _) = jacobi_eigen(&m);
        assert!((vals[0] - 3.0).abs() < 1e-10);
        assert!((vals[1] - 2.0).abs() < 1e-10);
        assert!((vals[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn fast_ica_recovers_two_sources() {
        // Two independent non-gaussian sources mixed linearly; FastICA should recover.
        let n = 400;
        let mut rng = Xoshiro256pp::new(7);
        let mut x = vec![vec![0.0_f64; 2]; n];
        for row in x.iter_mut() {
            // Uniform sources (sub-gaussian, non-normal).
            let s1 = (rng.next_f64() - 0.5) * 2.0;
            let s2 = (rng.next_f64() - 0.5) * 2.0;
            // Mixing: [[2, 1], [1, 3]]
            row[0] = 2.0 * s1 + 1.0 * s2;
            row[1] = 1.0 * s1 + 3.0 * s2;
        }
        let r = fast_ica(&x, 2, 42, 300, 1e-6);
        assert_eq!(r.sources.len(), n);
        assert_eq!(r.sources[0].len(), 2);
        // Recovery check: each true source should correlate strongly with one recovered source.
        let mut s1 = vec![0.0; n];
        let mut s2 = vec![0.0; n];
        let mut rng2 = Xoshiro256pp::new(7);
        for i in 0..n {
            s1[i] = (rng2.next_f64() - 0.5) * 2.0;
            s2[i] = (rng2.next_f64() - 0.5) * 2.0;
        }
        let recov0: Vec<f64> = r.sources.iter().map(|row| row[0]).collect();
        let recov1: Vec<f64> = r.sources.iter().map(|row| row[1]).collect();
        let c00 = abs_corr(&s1, &recov0);
        let c01 = abs_corr(&s1, &recov1);
        let c10 = abs_corr(&s2, &recov0);
        let c11 = abs_corr(&s2, &recov1);
        let best_1 = c00.max(c01);
        let best_2 = c10.max(c11);
        assert!(best_1 > 0.9, "source 1 best recovery = {best_1}");
        assert!(best_2 > 0.9, "source 2 best recovery = {best_2}");
    }

    #[test]
    fn fast_ica_is_deterministic_per_seed() {
        let n = 120;
        let mut rng = Xoshiro256pp::new(1);
        let mut x = vec![vec![0.0_f64; 3]; n];
        for row in x.iter_mut() {
            for v in row.iter_mut() {
                *v = rng.next_normal();
            }
        }
        let r1 = fast_ica(&x, 2, 17, 200, 1e-6);
        let r2 = fast_ica(&x, 2, 17, 200, 1e-6);
        for i in 0..r1.sources.len() {
            for c in 0..2 {
                assert_eq!(r1.sources[i][c].to_bits(), r2.sources[i][c].to_bits());
            }
        }
    }

    #[test]
    fn select_k_respects_clamp() {
        let spectrum = vec![5.0, 3.0, 1.5, 0.4, 0.1];
        // total=10, want >=0.8 -> need k=2 (5+3=8/10=0.80), but k_min=3 clamps up.
        assert_eq!(select_k_cumulative_variance(&spectrum, 0.8, 3, 10), 3);
        // k_max clamp.
        assert_eq!(select_k_cumulative_variance(&spectrum, 0.99, 1, 3), 3);
    }

    fn abs_corr(a: &[f64], b: &[f64]) -> f64 {
        let n = a.len() as f64;
        let ma = a.iter().sum::<f64>() / n;
        let mb = b.iter().sum::<f64>() / n;
        let mut sxx = 0.0;
        let mut syy = 0.0;
        let mut sxy = 0.0;
        for (x, y) in a.iter().zip(b.iter()) {
            let dx = x - ma;
            let dy = y - mb;
            sxx += dx * dx;
            syy += dy * dy;
            sxy += dx * dy;
        }
        (sxy / (sxx.sqrt() * syy.sqrt())).abs()
    }
}

#[cfg(test)]
mod weighted_whitening_tests {
    use super::*;

    /// Uniform weights must reproduce the unweighted whitening. If they
    /// do not, the weighted path is computing a different object and any
    /// comparison against plain ICA is confounded by that rather than by
    /// the weighting.
    #[test]
    fn uniform_weights_reproduce_plain_whitening() {
        let mut rng = Xoshiro256pp::new(11);
        let n = 40;
        let p = 12;
        let x: Vec<Vec<f64>> = (0..n)
            .map(|_| (0..p).map(|_| rng.next_normal()).collect())
            .collect();
        let ones = vec![vec![1.0_f64; p]; n];
        let plain = pca_whiten(&x, 3);
        let weighted = pca_whiten_weighted(&x, 3, &ones);
        for (a, b) in plain.mean.iter().zip(weighted.mean.iter()) {
            assert!((a - b).abs() < 1e-9, "means differ: {a} vs {b}");
        }
        for (ra, rb) in plain.whitened.iter().zip(weighted.whitened.iter()) {
            for (a, b) in ra.iter().zip(rb.iter()) {
                assert!(
                    (a.abs() - b.abs()).abs() < 1e-6,
                    "whitened coordinates differ: {a} vs {b}"
                );
            }
        }
    }

    /// A zero-weighted cell must not influence the moments. This is the
    /// property that makes the weighting mean anything: an undetected
    /// cell whose imputed value is arbitrary should not steer the
    /// covariance.
    #[test]
    fn zero_weighted_cells_do_not_move_the_moments() {
        let mut rng = Xoshiro256pp::new(7);
        let n = 30;
        let p = 8;
        let mut x: Vec<Vec<f64>> = (0..n)
            .map(|_| (0..p).map(|_| rng.next_normal()).collect())
            .collect();
        let mut weights = vec![vec![1.0_f64; p]; n];
        // Corrupt some cells and zero their weight.
        for i in 0..n {
            if i % 5 == 0 {
                x[i][3] = 1e6;
                weights[i][3] = 0.0;
            }
        }
        let w = pca_whiten_weighted(&x, 2, &weights);
        assert!(
            w.mean[3].abs() < 5.0,
            "a zero-weighted outlier still moved the mean: {}",
            w.mean[3]
        );
        assert!(
            w.eigenvalues.iter().all(|v| v.is_finite() && *v < 1e6),
            "a zero-weighted outlier still dominates the spectrum: {:?}",
            w.eigenvalues
        );
    }
}
