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

/// Deterministic Xoshiro256++ PRNG seeded via SplitMix64.
pub struct Xoshiro256pp {
    state: [u64; 4],
}

impl Xoshiro256pp {
    pub fn new(seed: u64) -> Self {
        let mut sm = seed.wrapping_add(0x9E3779B97F4A7C15);
        let mut next_seed = || {
            sm = sm.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = sm;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        };
        let state = [next_seed(), next_seed(), next_seed(), next_seed()];
        Self { state }
    }

    fn next_u64(&mut self) -> u64 {
        let result = self.state[0]
            .wrapping_add(self.state[3])
            .rotate_left(23)
            .wrapping_add(self.state[0]);
        let t = self.state[1] << 17;
        self.state[2] ^= self.state[0];
        self.state[3] ^= self.state[1];
        self.state[1] ^= self.state[2];
        self.state[0] ^= self.state[3];
        self.state[2] ^= t;
        self.state[3] = self.state[3].rotate_left(45);
        result
    }

    /// Uniform `[0, 1)` double.
    fn next_f64(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) * (1.0_f64 / ((1u64 << 53) as f64))
    }

    /// Standard normal draw via Box-Muller. Paired draws; we discard one.
    pub fn next_normal(&mut self) -> f64 {
        let mut u1 = self.next_f64();
        while u1 <= f64::MIN_POSITIVE {
            u1 = self.next_f64();
        }
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

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
        let mut off = 0.0_f64;
        for p in 0..n {
            for q in (p + 1)..n {
                off += a[p][q].abs();
            }
        }
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
                let mut acc = 0.0_f64;
                for t in 0..p {
                    acc += xc[i][t] * xc[j][t];
                }
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
    let w = pca_whiten(x, k);
    let n = w.whitened.len();

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
        // For each component, compute E[x g(w^T x)] - E[g'(w^T x)] w.
        for c in 0..k {
            let w_c: Vec<f64> = weights[c].clone();
            let mut sum_gxx = vec![0.0_f64; k];
            let mut sum_gp = 0.0_f64;
            for row in &w.whitened {
                let mut wx = 0.0_f64;
                for (j, xj) in row.iter().enumerate() {
                    wx += w_c[j] * xj;
                }
                let gwx = wx.tanh();
                let gpwx = 1.0 - gwx * gwx;
                sum_gp += gpwx;
                for (j, xj) in row.iter().enumerate() {
                    sum_gxx[j] += gwx * xj;
                }
            }
            let mean_gp = sum_gp / n as f64;
            for j in 0..k {
                next[c][j] = sum_gxx[j] / n as f64 - mean_gp * w_c[j];
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
    let p = w.whitening[0].len();
    let mut unmixing = vec![vec![0.0_f64; p]; k];
    for c in 0..k {
        for j in 0..p {
            let mut acc = 0.0_f64;
            for t in 0..k {
                acc += weights[c][t] * w.whitening[t][j];
            }
            unmixing[c][j] = acc;
        }
    }
    // mixing (p x k) = unwhitening * W^T  (since unwhitening = V * D^{1/2} and W is orthonormal)
    let mut mixing = vec![vec![0.0_f64; k]; p];
    for j in 0..p {
        for c in 0..k {
            let mut acc = 0.0_f64;
            for t in 0..k {
                acc += w.unwhitening[j][t] * weights[c][t];
            }
            mixing[j][c] = acc;
        }
    }

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

/// Jaccard similarity between the top-`n` |loading| index sets.
pub fn jaccard_top_n(a: &[f64], b: &[f64], top_n: usize) -> f64 {
    let set_a = top_abs_indices(a, top_n);
    let set_b = top_abs_indices(b, top_n);
    if set_a.is_empty() && set_b.is_empty() {
        return 0.0;
    }
    let inter = set_a.iter().filter(|i| set_b.contains(i)).count();
    let union = set_a.len() + set_b.len() - inter;
    if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    }
}

fn top_abs_indices(values: &[f64], top_n: usize) -> Vec<usize> {
    let mut indexed: Vec<(usize, f64)> = values
        .iter()
        .enumerate()
        .map(|(i, v)| (i, v.abs()))
        .collect();
    indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
    let keep = indexed.len().min(top_n);
    let mut out: Vec<usize> = indexed.into_iter().take(keep).map(|(i, _)| i).collect();
    out.sort_unstable();
    out
}

/// Symmetric decorrelation: W <- (W W^T)^{-1/2} W, in place.
fn symmetric_decorrelate(w: &mut Vec<Vec<f64>>) {
    let k = w.len();
    let mut wwt = vec![vec![0.0_f64; k]; k];
    for i in 0..k {
        for j in i..k {
            let mut acc = 0.0_f64;
            for t in 0..w[0].len() {
                acc += w[i][t] * w[j][t];
            }
            wwt[i][j] = acc;
            wwt[j][i] = acc;
        }
    }
    let (vals, vecs) = jacobi_eigen(&wwt);
    // (WW^T)^{-1/2} = V diag(1/sqrt(lambda)) V^T
    let mut inv_sqrt = vec![vec![0.0_f64; k]; k];
    for i in 0..k {
        for j in 0..k {
            let mut acc = 0.0_f64;
            for t in 0..k {
                let lam = vals[t].max(1e-18);
                acc += vecs[i][t] * (1.0 / lam.sqrt()) * vecs[j][t];
            }
            inv_sqrt[i][j] = acc;
        }
    }
    let p = w[0].len();
    let mut next = vec![vec![0.0_f64; p]; k];
    for i in 0..k {
        for j in 0..p {
            let mut acc = 0.0_f64;
            for t in 0..k {
                acc += inv_sqrt[i][t] * w[t][j];
            }
            next[i][j] = acc;
        }
    }
    *w = next;
}

fn max_identity_deviation(w_new: &[Vec<f64>], w_old: &[Vec<f64>]) -> f64 {
    let k = w_new.len();
    let mut max_dev = 0.0_f64;
    for i in 0..k {
        for j in 0..k {
            let mut acc = 0.0_f64;
            for t in 0..w_new[0].len() {
                acc += w_new[i][t] * w_old[j][t];
            }
            let expected = if i == j { 1.0 } else { 0.0 };
            let dev = (acc.abs() - expected).abs();
            if dev > max_dev {
                max_dev = dev;
            }
        }
    }
    max_dev
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
    fn jaccard_top_n_handles_basic_cases() {
        let a = vec![0.9, 0.1, -0.8, 0.2];
        let b = vec![0.85, 0.15, 0.0, -0.78];
        // top-2 of a: {0, 2}; top-2 of b: {0, 3}. intersection=1, union=3.
        let j = jaccard_top_n(&a, &b, 2);
        assert!((j - 1.0 / 3.0).abs() < 1e-12);
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
