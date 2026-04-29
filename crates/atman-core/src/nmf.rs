//! Multiplicative-updates NMF (Brunet 2004, Lee & Seung 1999/2001).
//!
//! Frobenius variant: minimises ||X - WH||_F^2 via the standard
//! multiplicative update rules. Initialization options are NNDSVDa
//! (Boutsidis & Gallopoulos 2008, deterministic) and Random
//! (Xoshiro256++ seeded, reproducible given `seed`).
//!
//! Data layout matches the ICA convention in `ica.rs`: all matrices
//! are `Vec<Vec<f64>>` (outer = rows). The public API converts from/to
//! that representation; the inner math uses index loops rather than
//! ndarray to avoid introducing new workspace dependencies.

use crate::ica::{jacobi_eigen, Xoshiro256pp};

// ── public types ────────────────────────────────────────────────────────────

/// Loss variant for NMF.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BetaLoss {
    Frobenius,
    /// Reserved for Task C; not yet implemented.
    KullbackLeibler,
}

/// Initialization strategy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Init {
    /// NNDSVDa: SVD-based init with zeros replaced by data mean.
    /// (Boutsidis & Gallopoulos 2008, deterministic given the input.)
    Nndsvda,
    /// Random non-negative init from `seed` (uniform [0, 1] scaled by
    /// `sqrt(mean(X) / k)`). Matches the Xoshiro256++ seeding in `ica.rs`.
    Random,
}

/// Configuration for [`nmf`].
#[derive(Debug, Clone, Copy)]
pub struct NmfConfig {
    pub k: usize,
    pub beta_loss: BetaLoss,
    pub init: Init,
    pub max_iter: usize,
    pub tol: f64,
    pub seed: u64,
}

/// Output of one NMF run.
#[derive(Debug)]
pub struct NmfResult {
    /// n × k sample activations (W).
    pub w: Vec<Vec<f64>>,
    /// k × p feature loadings (H).
    pub h: Vec<Vec<f64>>,
    /// Iterations actually run.
    pub n_iter: usize,
    /// Final reconstruction loss (Frobenius norm of residual).
    pub final_error: f64,
    /// Whether the algorithm converged before `max_iter`.
    pub converged: bool,
}

// ── entry point ─────────────────────────────────────────────────────────────

/// Multiplicative-updates NMF (Brunet 2004; Lee & Seung 1999/2001).
///
/// `x` is an n × p matrix (rows = samples, columns = features).
/// `cfg.k` is the number of components.
///
/// # Panics
/// Panics if `x` is empty, rows have unequal lengths, any entry is
/// negative, or `k` exceeds `min(n, p)`.
pub fn nmf(x: &[Vec<f64>], cfg: &NmfConfig) -> NmfResult {
    let n = x.len();
    assert!(n > 0, "nmf: x must be non-empty");
    let p = x[0].len();
    assert!(p > 0, "nmf: x must have at least one feature");
    assert!(
        x.iter().all(|r| r.len() == p),
        "nmf: all rows must have equal length"
    );
    assert!(
        x.iter().all(|r| r.iter().all(|&v| v >= 0.0)),
        "nmf: all entries of x must be non-negative"
    );
    assert!(cfg.k > 0 && cfg.k <= n.min(p), "nmf: k out of range");

    match cfg.beta_loss {
        BetaLoss::KullbackLeibler => {
            unimplemented!("KL-divergence NMF is reserved for Task C")
        }
        BetaLoss::Frobenius => {}
    }

    let (mut w, mut h) = match cfg.init {
        Init::Nndsvda => init_nndsvda(x, n, p, cfg.k),
        Init::Random => init_random(x, n, p, cfg.k, cfg.seed),
    };

    frobenius_mu(x, &mut w, &mut h, n, p, cfg.k, cfg.max_iter, cfg.tol)
}

// ── Frobenius multiplicative updates ────────────────────────────────────────

const EPSILON: f64 = 1e-9;

/// Run Lee & Seung multiplicative updates for the Frobenius loss.
/// Modifies `w` and `h` in place; returns the final [`NmfResult`].
fn frobenius_mu(
    x: &[Vec<f64>],
    w: &mut Vec<Vec<f64>>,
    h: &mut Vec<Vec<f64>>,
    n: usize,
    p: usize,
    k: usize,
    max_iter: usize,
    tol: f64,
) -> NmfResult {
    let mut prev_error = frobenius_error(x, w, h, n, p, k);
    let mut n_iter = 0usize;
    let mut converged = false;

    for iter in 0..max_iter {
        n_iter = iter + 1;

        // ── update H ──────────────────────────────────────────────────────
        // H_new[a, j] = H[a, j] * (W^T X)[a, j] / (W^T W H)[a, j]
        //
        // WtX[a, j]  = sum_i  W[i, a] * X[i, j]
        // WtWH[a, j] = sum_i (W^T W)[a, i] * H[i, j]
        //            = sum_i (sum_l W[l, a] * W[l, i]) * H[i, j]

        // precompute WtX (k × p)
        let wtx = mat_mul_at_b(w, x, k, n, p);
        // precompute WtW (k × k)
        let wtw = mat_mul_at_a(w, k, n);
        // WtWH (k × p) = WtW @ H
        let wtwh = mat_mul(k, k, &wtw, k, p, h);

        for a in 0..k {
            for j in 0..p {
                let num = wtx[a][j];
                let den = wtwh[a][j] + EPSILON;
                h[a][j] *= num / den;
                // clip negatives that floating-point errors might introduce
                if h[a][j] < 0.0 {
                    h[a][j] = 0.0;
                }
            }
        }

        // ── update W ──────────────────────────────────────────────────────
        // W_new[i, a] = W[i, a] * (X H^T)[i, a] / (W H H^T)[i, a]

        // XHt (n × k) = X @ H^T
        let xht = mat_mul_a_bt(x, h, n, p, k);
        // HHt (k × k) = H @ H^T
        let hht = mat_mul_a_at(h, k, p);
        // WHHt (n × k) = W @ HHt
        let whht = mat_mul(n, k, w, k, k, &hht);

        for i in 0..n {
            for a in 0..k {
                let num = xht[i][a];
                let den = whht[i][a] + EPSILON;
                w[i][a] *= num / den;
                if w[i][a] < 0.0 {
                    w[i][a] = 0.0;
                }
            }
        }

        // ── convergence check ─────────────────────────────────────────────
        let error = frobenius_error(x, w, h, n, p, k);
        if (prev_error - error).abs() < tol {
            converged = true;
            prev_error = error;
            break;
        }
        prev_error = error;
    }

    NmfResult {
        w: w.clone(),
        h: h.clone(),
        n_iter,
        final_error: prev_error,
        converged,
    }
}

/// ||X - WH||_F
fn frobenius_error(
    x: &[Vec<f64>],
    w: &[Vec<f64>],
    h: &[Vec<f64>],
    n: usize,
    p: usize,
    k: usize,
) -> f64 {
    let mut acc = 0.0_f64;
    for i in 0..n {
        for j in 0..p {
            let wh_ij: f64 = (0..k).map(|a| w[i][a] * h[a][j]).sum();
            let r = x[i][j] - wh_ij;
            acc += r * r;
        }
    }
    acc.sqrt()
}

// ── initializations ─────────────────────────────────────────────────────────

/// NNDSVDa initialisation (Boutsidis & Gallopoulos 2008).
///
/// For each of the top-k SVD triplets (u_r, s_r, v_r):
///   - split u_r into u+ = max(u, 0) and u- = max(-u, 0)
///   - split v_r into v+ = max(v, 0) and v- = max(-v, 0)
///   - choose the sign that gives the larger Frobenius contribution
///     (||u+|| * ||v+|| vs ||u-|| * ||v-||)
///   - set W[:, r] and H[r, :] accordingly, scaled by sqrt(s_r)
/// Remaining zeros (from the 0-norms) are replaced by mean(X) ("a" variant).
fn init_nndsvda(x: &[Vec<f64>], n: usize, p: usize, k: usize) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    // We compute the truncated SVD of X using the Gram-matrix trick from `pca_whiten`.
    // Singular values are sqrt(eigenvalues of X X^T or X^T X, whichever is smaller).
    // We need the top-k left (U, n×k) and right (V, p×k) singular vectors.

    // Flatten x for easier access
    let (u_nk, sigma_k, vt_kp) = truncated_svd(x, n, p, k);

    let mean_x: f64 = {
        let total: f64 = x.iter().flat_map(|r| r.iter()).sum();
        let count = (n * p) as f64;
        if count > 0.0 { total / count } else { 0.0 }
    };

    let mut w = vec![vec![0.0_f64; k]; n]; // n × k
    let mut h = vec![vec![0.0_f64; p]; k]; // k × p

    for r in 0..k {
        let s = sigma_k[r].max(0.0).sqrt(); // scale = sqrt(sigma_r)

        // u_r is the r-th column of U (length n)
        let u_r: Vec<f64> = (0..n).map(|i| u_nk[i][r]).collect();
        // v_r is the r-th row of Vt (length p)
        let v_r: Vec<f64> = vt_kp[r].clone();

        let (up, un) = pos_neg_split(&u_r);
        let (vp, vn) = pos_neg_split(&v_r);

        let norm_up = frobenius_norm_vec(&up);
        let norm_un = frobenius_norm_vec(&un);
        let norm_vp = frobenius_norm_vec(&vp);
        let norm_vn = frobenius_norm_vec(&vn);

        let pos_contrib = norm_up * norm_vp;
        let neg_contrib = norm_un * norm_vn;

        let (w_r, h_r): (Vec<f64>, Vec<f64>) = if pos_contrib >= neg_contrib {
            let scale = s * (pos_contrib).sqrt().max(EPSILON).recip();
            let w_r = up.iter().map(|v| v * s / (norm_up.max(EPSILON))).collect();
            let h_r = vp.iter().map(|v| v * s / (norm_vp.max(EPSILON))).collect();
            let _ = scale;
            (w_r, h_r)
        } else {
            let w_r = un.iter().map(|v| v * s / (norm_un.max(EPSILON))).collect();
            let h_r = vn.iter().map(|v| v * s / (norm_vn.max(EPSILON))).collect();
            (w_r, h_r)
        };

        for i in 0..n {
            w[i][r] = w_r[i];
        }
        for j in 0..p {
            h[r][j] = h_r[j];
        }
    }

    // "a" variant: replace zeros with mean(X)
    for i in 0..n {
        for a in 0..k {
            if w[i][a] == 0.0 {
                w[i][a] = mean_x;
            }
        }
    }
    for a in 0..k {
        for j in 0..p {
            if h[a][j] == 0.0 {
                h[a][j] = mean_x;
            }
        }
    }

    (w, h)
}

/// Random init: uniform [0, 1) scaled by `sqrt(mean(X) / k)`.
fn init_random(
    x: &[Vec<f64>],
    n: usize,
    p: usize,
    k: usize,
    seed: u64,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let mean_x: f64 = {
        let total: f64 = x.iter().flat_map(|r| r.iter()).sum();
        let count = (n * p) as f64;
        if count > 0.0 { total / count } else { 1.0 }
    };
    let scale = (mean_x / k as f64).sqrt().max(EPSILON);
    let mut rng = Xoshiro256pp::new(seed);

    // next_normal gives Box-Muller standard normal; we want uniform [0,1).
    // We derive uniform from the raw u64 the same way Xoshiro256pp::next_f64 does
    // internally (mantissa extraction), accessed through the public next_normal
    // isn't ideal; instead we draw |normal| which is half-normal ≥ 0 and widely
    // used as a non-negative random init.  This is reproducible given `seed`.
    let w: Vec<Vec<f64>> = (0..n)
        .map(|_| (0..k).map(|_| rng.next_normal().abs() * scale).collect())
        .collect();
    let h: Vec<Vec<f64>> = (0..k)
        .map(|_| (0..p).map(|_| rng.next_normal().abs() * scale).collect())
        .collect();
    (w, h)
}

// ── truncated SVD via Gram matrix + Jacobi ───────────────────────────────────

/// Compute the top-k truncated SVD of x (n × p) via the smaller Gram matrix.
/// Returns (U_nk, sigma_k, Vt_kp) where:
///   U_nk[i][r]  = i-th element of r-th left  singular vector
///   sigma_k[r]  = r-th singular value (not squared)
///   Vt_kp[r][j] = j-th element of r-th right singular vector
fn truncated_svd(
    x: &[Vec<f64>],
    n: usize,
    p: usize,
    k: usize,
) -> (Vec<Vec<f64>>, Vec<f64>, Vec<Vec<f64>>) {
    if n <= p {
        // Gram matrix X X^T (n × n), eigendecompose to get U and sigma.
        let gram = gram_xxt(x, n);
        let (eigs, evecs_n) = jacobi_eigen(&gram); // sorted descending

        // sigma_r = sqrt(max(lambda_r, 0))
        let sigma: Vec<f64> = eigs.iter().take(k).map(|&e| e.max(0.0).sqrt()).collect();

        // U (n × k): columns are the top-k eigenvectors of X X^T
        let u_nk: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..k).map(|r| evecs_n[i][r]).collect())
            .collect();

        // V (p × k): V_r = X^T U_r / sigma_r
        let mut vt_kp = vec![vec![0.0_f64; p]; k];
        for r in 0..k {
            let s = sigma[r].max(EPSILON);
            for j in 0..p {
                let mut acc = 0.0_f64;
                for i in 0..n {
                    acc += x[i][j] * evecs_n[i][r];
                }
                vt_kp[r][j] = acc / s;
            }
        }

        (u_nk, sigma, vt_kp)
    } else {
        // Gram matrix X^T X (p × p), eigendecompose to get V and sigma.
        let gram = gram_xtx(x, n, p);
        let (eigs, evecs_p) = jacobi_eigen(&gram);

        let sigma: Vec<f64> = eigs.iter().take(k).map(|&e| e.max(0.0).sqrt()).collect();

        // V (p × k)
        let vt_kp: Vec<Vec<f64>> = (0..k)
            .map(|r| (0..p).map(|j| evecs_p[j][r]).collect())
            .collect();

        // U (n × k): U_r = X V_r / sigma_r
        let mut u_nk = vec![vec![0.0_f64; k]; n];
        for r in 0..k {
            let s = sigma[r].max(EPSILON);
            for i in 0..n {
                let mut acc = 0.0_f64;
                for j in 0..p {
                    acc += x[i][j] * evecs_p[j][r];
                }
                u_nk[i][r] = acc / s;
            }
        }

        (u_nk, sigma, vt_kp)
    }
}

fn gram_xxt(x: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    let mut g = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in i..n {
            let acc: f64 = x[i].iter().zip(x[j].iter()).map(|(a, b)| a * b).sum();
            g[i][j] = acc;
            g[j][i] = acc;
        }
    }
    g
}

fn gram_xtx(x: &[Vec<f64>], n: usize, p: usize) -> Vec<Vec<f64>> {
    let mut g = vec![vec![0.0_f64; p]; p];
    for i in 0..p {
        for j in i..p {
            let mut acc = 0.0_f64;
            for r in 0..n {
                acc += x[r][i] * x[r][j];
            }
            g[i][j] = acc;
            g[j][i] = acc;
        }
    }
    g
}

// ── small matrix helpers ─────────────────────────────────────────────────────

/// (rows_a × cols_b) = (rows_a × shared) @ (shared × cols_b)
fn mat_mul(rows_a: usize, shared: usize, a: &[Vec<f64>], _shared2: usize, cols_b: usize, b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut out = vec![vec![0.0_f64; cols_b]; rows_a];
    for i in 0..rows_a {
        for j in 0..cols_b {
            let mut acc = 0.0_f64;
            for l in 0..shared {
                acc += a[i][l] * b[l][j];
            }
            out[i][j] = acc;
        }
    }
    out
}

/// A^T @ B: (cols_a × cols_b) from (rows × cols_a) and (rows × cols_b)
fn mat_mul_at_b(a: &[Vec<f64>], b: &[Vec<f64>], cols_a: usize, rows: usize, cols_b: usize) -> Vec<Vec<f64>> {
    let mut out = vec![vec![0.0_f64; cols_b]; cols_a];
    for r in 0..rows {
        for i in 0..cols_a {
            let ai = a[r][i];
            for j in 0..cols_b {
                out[i][j] += ai * b[r][j];
            }
        }
    }
    out
}

/// A^T @ A: (cols × cols) from (rows × cols)
fn mat_mul_at_a(a: &[Vec<f64>], cols: usize, rows: usize) -> Vec<Vec<f64>> {
    let mut out = vec![vec![0.0_f64; cols]; cols];
    for r in 0..rows {
        for i in 0..cols {
            let ai = a[r][i];
            for j in i..cols {
                let v = ai * a[r][j];
                out[i][j] += v;
                if i != j {
                    out[j][i] += v;
                }
            }
        }
    }
    out
}

/// A @ B^T: (rows_a × rows_b) from (rows_a × shared) and (rows_b × shared)
fn mat_mul_a_bt(a: &[Vec<f64>], b: &[Vec<f64>], rows_a: usize, _shared: usize, rows_b: usize) -> Vec<Vec<f64>> {
    let shared = a[0].len();
    let mut out = vec![vec![0.0_f64; rows_b]; rows_a];
    for i in 0..rows_a {
        for j in 0..rows_b {
            let mut acc = 0.0_f64;
            for l in 0..shared {
                acc += a[i][l] * b[j][l];
            }
            out[i][j] = acc;
        }
    }
    out
}

/// A @ A^T: (rows × rows) from (rows × cols)
fn mat_mul_a_at(a: &[Vec<f64>], rows: usize, cols: usize) -> Vec<Vec<f64>> {
    let mut out = vec![vec![0.0_f64; rows]; rows];
    for i in 0..rows {
        for j in i..rows {
            let acc: f64 = (0..cols).map(|l| a[i][l] * a[j][l]).sum();
            out[i][j] = acc;
            out[j][i] = acc;
        }
    }
    out
}

// ── misc helpers ─────────────────────────────────────────────────────────────

fn pos_neg_split(v: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let pos: Vec<f64> = v.iter().map(|&x| x.max(0.0)).collect();
    let neg: Vec<f64> = v.iter().map(|&x| (-x).max(0.0)).collect();
    (pos, neg)
}

fn frobenius_norm_vec(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

// ── tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Convenience: build a Vec<Vec<f64>> from a 2-D slice literal.
    fn mat(rows: &[&[f64]]) -> Vec<Vec<f64>> {
        rows.iter().map(|r| r.to_vec()).collect()
    }

    #[test]
    fn frobenius_recovers_synthetic_rank_1() {
        // Rank-1 synthetic: x = w * h with w = [1, 2, 3]^T, h = [1.5, 2.5]
        let x = mat(&[&[1.5, 2.5], &[3.0, 5.0], &[4.5, 7.5]]);
        let cfg = NmfConfig {
            k: 1,
            beta_loss: BetaLoss::Frobenius,
            init: Init::Nndsvda,
            max_iter: 500,
            tol: 1e-10,
            seed: 42,
        };
        let r = nmf(&x, &cfg);
        // Compute reconstruction error manually.
        let n = x.len();
        let p = x[0].len();
        let err: f64 = {
            let mut acc = 0.0_f64;
            for i in 0..n {
                for j in 0..p {
                    let wh: f64 = (0..cfg.k).map(|a| r.w[i][a] * r.h[a][j]).sum();
                    let d = x[i][j] - wh;
                    acc += d * d;
                }
            }
            acc.sqrt()
        };
        assert!(err < 1e-3, "rank-1 reconstruction err = {err:.6}");
        assert!(r.converged, "expected convergence");
    }

    #[test]
    fn frobenius_recovers_synthetic_rank_2_random_init() {
        // 6×4 matrix from rank-2 ground truth.
        // W_true = [[1,0],[0,1],[1,1],[2,0],[0,2],[1,2]], H_true = [[1,0,1,2],[0,1,2,1]]
        let w_true: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 1.0],
            vec![2.0, 0.0],
            vec![0.0, 2.0],
            vec![1.0, 2.0],
        ];
        let h_true: Vec<Vec<f64>> = vec![vec![1.0, 0.0, 1.0, 2.0], vec![0.0, 1.0, 2.0, 1.0]];
        // x = w_true @ h_true
        let x: Vec<Vec<f64>> = (0..6)
            .map(|i| (0..4).map(|j| w_true[i][0] * h_true[0][j] + w_true[i][1] * h_true[1][j]).collect())
            .collect();

        let cfg = NmfConfig {
            k: 2,
            beta_loss: BetaLoss::Frobenius,
            init: Init::Random,
            max_iter: 2000,
            tol: 1e-10,
            seed: 42,
        };
        let r = nmf(&x, &cfg);
        let n = x.len();
        let p = x[0].len();
        let err: f64 = {
            let mut acc = 0.0_f64;
            for i in 0..n {
                for j in 0..p {
                    let wh: f64 = (0..cfg.k).map(|a| r.w[i][a] * r.h[a][j]).sum();
                    let d = x[i][j] - wh;
                    acc += d * d;
                }
            }
            acc.sqrt()
        };
        assert!(err < 1e-2, "rank-2 reconstruction err = {err:.6}");
    }

    #[test]
    fn frobenius_returns_non_negative() {
        let x = mat(&[&[1.0, 2.0], &[2.0, 4.0]]);
        let cfg = NmfConfig {
            k: 1,
            beta_loss: BetaLoss::Frobenius,
            init: Init::Nndsvda,
            max_iter: 200,
            tol: 1e-9,
            seed: 42,
        };
        let r = nmf(&x, &cfg);
        assert!(
            r.w.iter().all(|row| row.iter().all(|&v| v >= 0.0)),
            "W contains negatives"
        );
        assert!(
            r.h.iter().all(|row| row.iter().all(|&v| v >= 0.0)),
            "H contains negatives"
        );
    }

    #[test]
    fn frobenius_deterministic_under_nndsvda() {
        // Same input + NNDSVDa init + same params → byte-identical W, H, n_iter.
        let x = mat(&[
            &[1.0, 2.0, 3.0],
            &[4.0, 5.0, 6.0],
            &[7.0, 8.0, 9.0],
        ]);
        let cfg = NmfConfig {
            k: 2,
            beta_loss: BetaLoss::Frobenius,
            init: Init::Nndsvda,
            max_iter: 100,
            tol: 1e-9,
            seed: 42,
        };
        let r1 = nmf(&x, &cfg);
        let r2 = nmf(&x, &cfg);
        assert_eq!(r1.n_iter, r2.n_iter);
        for i in 0..3 {
            for a in 0..2 {
                assert_eq!(
                    r1.w[i][a].to_bits(),
                    r2.w[i][a].to_bits(),
                    "W[{i}][{a}] differs"
                );
            }
        }
        for a in 0..2 {
            for j in 0..3 {
                assert_eq!(
                    r1.h[a][j].to_bits(),
                    r2.h[a][j].to_bits(),
                    "H[{a}][{j}] differs"
                );
            }
        }
    }

    #[test]
    fn frobenius_deterministic_under_random_init() {
        let x = mat(&[
            &[1.0, 2.0, 3.0],
            &[4.0, 5.0, 6.0],
            &[7.0, 8.0, 9.0],
        ]);
        let cfg = NmfConfig {
            k: 2,
            beta_loss: BetaLoss::Frobenius,
            init: Init::Random,
            max_iter: 50,
            tol: 1e-9,
            seed: 1234,
        };
        let r1 = nmf(&x, &cfg);
        let r2 = nmf(&x, &cfg);
        assert_eq!(r1.n_iter, r2.n_iter);
        for i in 0..3 {
            for a in 0..2 {
                assert_eq!(r1.w[i][a].to_bits(), r2.w[i][a].to_bits());
            }
        }
    }
}
