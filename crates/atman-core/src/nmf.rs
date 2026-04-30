//! Multiplicative-updates NMF (Brunet 2004, Lee & Seung 1999/2001).
//!
//! Frobenius variant: minimises ||X - WH||_F^2 via the standard
//! multiplicative update rules.  KL variant: minimises the generalised
//! KL divergence D_KL(X || WH) = sum X*log(X/WH) - X + WH via the Lee
//! & Seung 2001 (NIPS) multiplicative update rules.
//!
//! Initialization options are NNDSVDa (Boutsidis & Gallopoulos 2008,
//! deterministic) and Random (Xoshiro256++ seeded, reproducible given
//! `seed`).
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
    /// Lee & Seung 2001 (NIPS) KL-divergence multiplicative updates.
    /// Minimises D_KL(X || WH) = sum_ij [ X_ij * log(X_ij / (WH)_ij) - X_ij + (WH)_ij ].
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

    let (mut w, mut h) = match cfg.init {
        Init::Nndsvda => init_nndsvda(x, n, p, cfg.k),
        Init::Random => init_random(x, n, p, cfg.k, cfg.seed),
    };

    match cfg.beta_loss {
        BetaLoss::Frobenius => frobenius_mu(x, &mut w, &mut h, n, p, cfg.k, cfg.max_iter, cfg.tol),
        BetaLoss::KullbackLeibler => kl_mu(x, &mut w, &mut h, n, p, cfg.k, cfg.max_iter, cfg.tol),
    }
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

// ── KL-divergence multiplicative updates ────────────────────────────────────

/// Run Lee & Seung 2001 (NIPS) multiplicative updates for the KL-divergence loss.
///
/// Loss: D_KL(X || WH) = sum_ij [ X_ij * log(X_ij / (WH)_ij) - X_ij + (WH)_ij ]
///
/// Update rules (element-wise):
///   H_kj ← H_kj * (sum_i  W_ik * X_ij / (WH)_ij) / (sum_i W_ik)
///   W_ik ← W_ik * (sum_j  H_kj * X_ij / (WH)_ij) / (sum_j H_kj)
fn kl_mu(
    x: &[Vec<f64>],
    w: &mut Vec<Vec<f64>>,
    h: &mut Vec<Vec<f64>>,
    n: usize,
    p: usize,
    k: usize,
    max_iter: usize,
    tol: f64,
) -> NmfResult {
    let mut prev_loss = kl_loss(x, w, h, n, p, k);
    let mut n_iter = 0usize;
    let mut converged = false;

    for iter in 0..max_iter {
        n_iter = iter + 1;

        // ── update H ──────────────────────────────────────────────────────
        // H_kj ← H_kj * (sum_i W_ik * X_ij / (WH)_ij) / (sum_i W_ik)
        //
        // Precompute WH (n × p).
        let wh = mat_mul(n, k, w, k, p, h);

        // Numerator for each (k, j): sum_i W[i][k] * X[i][j] / WH[i][j]
        // Denominator for each k:    sum_i W[i][k]
        let mut h_num = vec![vec![0.0_f64; p]; k];
        let mut w_col_sum = vec![0.0_f64; k];

        for i in 0..n {
            for kk in 0..k {
                let w_ik = w[i][kk];
                w_col_sum[kk] += w_ik;
                for j in 0..p {
                    let wh_ij = wh[i][j] + EPSILON;
                    h_num[kk][j] += w_ik * x[i][j] / wh_ij;
                }
            }
        }

        for kk in 0..k {
            let den = w_col_sum[kk] + EPSILON;
            for j in 0..p {
                h[kk][j] *= h_num[kk][j] / den;
                if h[kk][j] < 0.0 {
                    h[kk][j] = 0.0;
                }
            }
        }

        // ── update W ──────────────────────────────────────────────────────
        // W_ik ← W_ik * (sum_j H_kj * X_ij / (WH)_ij) / (sum_j H_kj)
        //
        // Recompute WH after H update.
        let wh = mat_mul(n, k, w, k, p, h);

        // Numerator for each (i, k): sum_j H[k][j] * X[i][j] / WH[i][j]
        // Denominator for each k:    sum_j H[k][j]
        let mut w_num = vec![vec![0.0_f64; k]; n];
        let mut h_row_sum = vec![0.0_f64; k];

        for kk in 0..k {
            for j in 0..p {
                h_row_sum[kk] += h[kk][j];
            }
        }

        for i in 0..n {
            for j in 0..p {
                let wh_ij = wh[i][j] + EPSILON;
                let x_ij = x[i][j];
                for kk in 0..k {
                    w_num[i][kk] += h[kk][j] * x_ij / wh_ij;
                }
            }
        }

        for i in 0..n {
            for kk in 0..k {
                let den = h_row_sum[kk] + EPSILON;
                w[i][kk] *= w_num[i][kk] / den;
                if w[i][kk] < 0.0 {
                    w[i][kk] = 0.0;
                }
            }
        }

        // ── convergence check ─────────────────────────────────────────────
        let loss = kl_loss(x, w, h, n, p, k);
        if (prev_loss - loss).abs() < tol {
            converged = true;
            prev_loss = loss;
            break;
        }
        prev_loss = loss;
    }

    NmfResult {
        w: w.clone(),
        h: h.clone(),
        n_iter,
        final_error: prev_loss,
        converged,
    }
}

/// Generalised KL divergence D_KL(X || WH).
///
/// D = sum_ij [ X_ij * log(X_ij / (WH)_ij) - X_ij + (WH)_ij ]
/// Convention: 0 * log(0) = 0 (handled via the EPSILON floor on WH).
fn kl_loss(
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
            let wh_ij: f64 = (0..k).map(|kk| w[i][kk] * h[kk][j]).sum::<f64>().max(EPSILON);
            let x_ij = x[i][j];
            if x_ij > 0.0 {
                acc += x_ij * (x_ij / wh_ij).ln() - x_ij + wh_ij;
            } else {
                // x_ij == 0: contribution is 0*log(0) - 0 + WH = WH
                acc += wh_ij;
            }
        }
    }
    acc
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

// ── multi-seed stability ─────────────────────────────────────────────────────

/// Stability metric variant used by [`multi_seed_nmf`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StabilityMetric {
    /// Jaccard overlap of the top-N absolute-loading indices.
    JaccardTopN,
}

/// Configuration for [`multi_seed_nmf`].
#[derive(Debug, Clone, Copy)]
pub struct MultiSeedConfig {
    /// Number of independent random seeds to run.
    pub n_seeds: usize,
    /// Base seed; seed `i` is `seed_base + i as u64`.
    pub seed_base: u64,
    /// Stability metric used for program clustering across seeds.
    pub stability_metric: StabilityMetric,
    /// Top-N loadings used by the Jaccard metric.
    pub stability_top_n: usize,
}

/// One row in the multi-seed stability output.
#[derive(Debug, Clone)]
pub struct StableProgramRow {
    /// Label, e.g. `"program_01"`.
    pub program_id: String,
    /// Fraction of seeds in which this program was recovered (0.0–1.0).
    pub stable_seed_fraction: f64,
    /// Mean loading vector averaged over all seeds in which the program appeared.
    /// Length = number of features (`p`).
    pub mean_loading: Vec<f64>,
}

/// Run NMF `ms_cfg.n_seeds` times with seeds `ms_cfg.seed_base + i`,
/// cluster programs across seeds by Jaccard top-N overlap on H rows
/// (reference = seed 0), and return per-program stability rows.
///
/// NMF loadings are non-negative by construction so there is no sign
/// ambiguity — clustering directly compares H rows without sign flipping.
///
/// Each result uses `nmf_cfg.init` and `nmf_cfg.beta_loss`; `nmf_cfg.seed`
/// is overridden per iteration to `ms_cfg.seed_base + i`.
pub fn multi_seed_nmf(
    x: &[Vec<f64>],
    nmf_cfg: &NmfConfig,
    ms_cfg: &MultiSeedConfig,
) -> Vec<StableProgramRow> {
    assert!(ms_cfg.n_seeds >= 1, "n_seeds must be at least 1");

    let k = nmf_cfg.k;
    let p = if x.is_empty() { 0 } else { x[0].len() };

    // ── run all seeds ────────────────────────────────────────────────────────
    let mut all_h: Vec<Vec<Vec<f64>>> = Vec::with_capacity(ms_cfg.n_seeds);
    for i in 0..ms_cfg.n_seeds {
        let seed = ms_cfg.seed_base.wrapping_add(i as u64);
        let cfg_i = NmfConfig { seed, ..*nmf_cfg };
        let result = nmf(x, &cfg_i);
        // result.h is k × p; store as Vec<Vec<f64>> where each inner Vec is one program row.
        all_h.push(result.h);
    }

    // ── build per-reference-program clusters ─────────────────────────────────
    // Reference = seed 0.  For each alternative seed, greedily match its k
    // programs to reference programs (best Jaccard, each program matched once).
    let reference = &all_h[0];

    // For each reference program `a`, collect loadings from each seed that
    // matched it (plus the reference itself).
    let mut seed_loadings_per_ref: Vec<Vec<Vec<f64>>> = (0..k)
        .map(|a| vec![reference[a].clone()])
        .collect();

    for seed_idx in 1..ms_cfg.n_seeds {
        let alt = &all_h[seed_idx];

        // Greedy matching: for each reference program (in order), find the
        // best-matching alternative program not yet assigned.
        let mut used = vec![false; k];
        for ref_a in 0..k {
            let ref_row = &reference[ref_a];
            let mut best_j = 0.0_f64;
            let mut best_alt = usize::MAX;
            for alt_a in 0..k {
                if used[alt_a] {
                    continue;
                }
                let j = jaccard_top_n_vec(ref_row, &alt[alt_a], ms_cfg.stability_top_n);
                if j > best_j {
                    best_j = j;
                    best_alt = alt_a;
                }
            }
            if best_alt != usize::MAX {
                used[best_alt] = true;
                seed_loadings_per_ref[ref_a].push(alt[best_alt].clone());
            }
        }
    }

    // ── compute stability fraction and mean loadings ─────────────────────────
    let n_seeds_f = ms_cfg.n_seeds as f64;
    let mut rows: Vec<StableProgramRow> = Vec::with_capacity(k);

    for (a, matched_loadings) in seed_loadings_per_ref.into_iter().enumerate() {
        let n_present = matched_loadings.len(); // always >= 1 (includes seed 0)
        let stable_seed_fraction = n_present as f64 / n_seeds_f;

        // Mean loading = element-wise average over matched loadings (no sign flip needed).
        let mut mean_loading = vec![0.0_f64; p];
        for row in &matched_loadings {
            for (j, &v) in row.iter().enumerate() {
                mean_loading[j] += v;
            }
        }
        let n_f = n_present as f64;
        for v in mean_loading.iter_mut() {
            *v /= n_f;
        }

        let prog_idx = a + 1;
        let program_id = format!("program_{:02}", prog_idx);

        rows.push(StableProgramRow {
            program_id,
            stable_seed_fraction,
            mean_loading,
        });
    }

    rows
}

/// Jaccard overlap of top-`top_n` absolute-value index sets for two `Vec<f64>` rows.
/// Duplicates the `jaccard_top_n` helper from `stats.rs` to avoid a crate-level
/// dependency loop; kept private to this module.
fn jaccard_top_n_vec(a: &[f64], b: &[f64], top_n: usize) -> f64 {
    let set_a = top_abs_indices_vec(a, top_n);
    let set_b = top_abs_indices_vec(b, top_n);
    if set_a.is_empty() && set_b.is_empty() {
        return 0.0;
    }
    let inter = set_a.iter().filter(|i| set_b.contains(i)).count();
    let union_sz = set_a.len() + set_b.len() - inter;
    if union_sz == 0 {
        0.0
    } else {
        inter as f64 / union_sz as f64
    }
}

/// Top-`n` indices by absolute value, sorted ascending (for set membership).
fn top_abs_indices_vec(v: &[f64], n: usize) -> Vec<usize> {
    let n = n.min(v.len());
    if n == 0 {
        return Vec::new();
    }
    let mut indexed: Vec<(usize, f64)> = v.iter().enumerate().map(|(i, &x)| (i, x.abs())).collect();
    indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let mut out: Vec<usize> = indexed.iter().take(n).map(|&(i, _)| i).collect();
    out.sort_unstable();
    out
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

    // ── multi-seed stability tests ───────────────────────────────────────────

    #[test]
    fn multi_seed_nmf_recovers_k_stable_programs_on_synthetic_rank_3() {
        // Build a rank-3 synthetic matrix: X = W_true @ H_true
        // W_true: 20 × 3, H_true: 3 × 10
        let n = 20;
        let p = 10;
        let k_true = 3;

        // Construct distinct non-negative basis rows for H.
        let h_true: Vec<Vec<f64>> = vec![
            // program 1: loads mainly on features 0–2
            vec![5.0, 4.0, 3.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            // program 2: loads mainly on features 3–5
            vec![0.1, 0.1, 0.1, 5.0, 4.0, 3.0, 0.1, 0.1, 0.1, 0.1],
            // program 3: loads mainly on features 7–9
            vec![0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 5.0, 4.0, 3.0],
        ];

        // W_true: cyclic pattern so each program gets reasonable activation.
        let mut w_true = vec![vec![0.0_f64; k_true]; n];
        for i in 0..n {
            w_true[i][i % k_true] = 2.0 + (i as f64 * 0.1);
            w_true[i][(i + 1) % k_true] = 0.5;
        }

        // Build X = W_true @ H_true.
        let x: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..p)
                    .map(|j| {
                        (0..k_true)
                            .map(|a| w_true[i][a] * h_true[a][j])
                            .sum::<f64>()
                    })
                    .collect()
            })
            .collect();

        let nmf_cfg = NmfConfig {
            k: k_true,
            beta_loss: BetaLoss::Frobenius,
            init: Init::Random,
            max_iter: 600,
            tol: 1e-8,
            seed: 0, // will be overridden per-seed
        };
        let ms_cfg = MultiSeedConfig {
            n_seeds: 5,
            seed_base: 100,
            stability_metric: StabilityMetric::JaccardTopN,
            stability_top_n: 3,
        };

        let rows = multi_seed_nmf(&x, &nmf_cfg, &ms_cfg);
        assert_eq!(rows.len(), k_true, "should return k rows");

        let n_stable = rows
            .iter()
            .filter(|r| r.stable_seed_fraction >= 0.6)
            .count();
        assert!(
            n_stable >= k_true,
            "expected all {} programs stable (>=0.6), got {}: {:?}",
            k_true,
            n_stable,
            rows.iter()
                .map(|r| r.stable_seed_fraction)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn multi_seed_nmf_is_deterministic_given_seed_base() {
        let x = mat(&[
            &[3.0, 2.0, 0.1, 0.1],
            &[2.5, 3.0, 0.1, 0.1],
            &[0.1, 0.1, 4.0, 3.0],
            &[0.1, 0.1, 3.0, 4.5],
            &[1.5, 1.0, 2.0, 2.0],
            &[1.0, 1.5, 2.5, 2.0],
        ]);
        let nmf_cfg = NmfConfig {
            k: 2,
            beta_loss: BetaLoss::Frobenius,
            init: Init::Random,
            max_iter: 200,
            tol: 1e-8,
            seed: 0,
        };
        let ms_cfg = MultiSeedConfig {
            n_seeds: 3,
            seed_base: 42,
            stability_metric: StabilityMetric::JaccardTopN,
            stability_top_n: 2,
        };

        let r1 = multi_seed_nmf(&x, &nmf_cfg, &ms_cfg);
        let r2 = multi_seed_nmf(&x, &nmf_cfg, &ms_cfg);

        assert_eq!(r1.len(), r2.len());
        for (a, b) in r1.iter().zip(r2.iter()) {
            assert_eq!(a.program_id, b.program_id);
            assert_eq!(
                a.stable_seed_fraction.to_bits(),
                b.stable_seed_fraction.to_bits()
            );
            assert_eq!(a.mean_loading.len(), b.mean_loading.len());
            for (va, vb) in a.mean_loading.iter().zip(b.mean_loading.iter()) {
                assert_eq!(
                    va.to_bits(),
                    vb.to_bits(),
                    "mean_loading differs: {} vs {}",
                    va,
                    vb
                );
            }
        }
    }

    // ── KL-divergence tests ──────────────────────────────────────────────────

    #[test]
    fn kl_recovers_synthetic_rank_1() {
        // Rank-1 synthetic: x = w * h with w = [1, 2, 3]^T, h = [1.5, 2.5].
        // KL MU is also exact on rank-1 nonneg data.
        let x = mat(&[&[1.5, 2.5], &[3.0, 5.0], &[4.5, 7.5]]);
        let cfg = NmfConfig {
            k: 1,
            beta_loss: BetaLoss::KullbackLeibler,
            init: Init::Nndsvda,
            max_iter: 500,
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
        assert!(err < 1e-3, "KL rank-1 reconstruction err = {err:.6}");
        assert!(r.converged, "expected convergence");
    }

    #[test]
    fn kl_returns_non_negative() {
        let x = mat(&[&[1.0, 2.0], &[2.0, 4.0]]);
        let cfg = NmfConfig {
            k: 1,
            beta_loss: BetaLoss::KullbackLeibler,
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
    fn kl_deterministic_under_nndsvda() {
        // Same input + NNDSVDa init + same params → byte-identical W, H, n_iter.
        let x = mat(&[
            &[1.0, 2.0, 3.0],
            &[4.0, 5.0, 6.0],
            &[7.0, 8.0, 9.0],
        ]);
        let cfg = NmfConfig {
            k: 2,
            beta_loss: BetaLoss::KullbackLeibler,
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
}
