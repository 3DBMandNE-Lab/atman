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

/// Default clamp radius `c` for [`Transform::Exp2Clip`] (`--transform-clamp`
/// on `decompose nmf` / `align project`). Shared so both CLI surfaces default
/// to the identical value.
pub const DEFAULT_EXP2_CLIP_CLAMP: f64 = 6.0;

/// Pre-decomposition input transform for `decompose nmf` / `align project`.
///
/// NMF requires non-negative input; Atman's canonical abundance is often a
/// signed log2 ratio (e.g. CPTAC's log2-ratio-to-reference scale), so these
/// transforms restore non-negativity ahead of the decomposition. Pure,
/// deterministic, RNG-free elementwise operations — see [`apply_transform`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Transform {
    /// No-op. Caller is responsible for verifying non-negativity.
    None,
    /// `y = 2^clamp(x, -clamp, +clamp)`. Restores a non-negative ratio scale
    /// from log2-ratio input while winsorizing rare extreme tails.
    Exp2Clip { clamp: f64 },
    /// `y = x - min(X)` over the whole matrix. Sensitivity alternative to
    /// `Exp2Clip`. The shift is always recomputed on whatever matrix is
    /// passed in — callers projecting a new cohort must NOT reuse a
    /// training-time shift; the resulting [`TransformRecord::shift`] is the
    /// value actually applied to *this* call's matrix.
    ShiftMin,
}

/// Resolved parameters of an [`apply_transform`] call, for sidecar
/// provenance. Exactly one of `clamp` / `shift` is `Some`, matching the
/// `Transform` variant applied (`None` for both when `Transform::None`).
#[derive(Debug, Clone, PartialEq)]
pub struct TransformRecord {
    /// `"none"`, `"exp2-clip"`, or `"shift-min"`.
    pub name: String,
    /// The clamp radius used, when `Transform::Exp2Clip`.
    pub clamp: Option<f64>,
    /// The global matrix minimum subtracted, when `Transform::ShiftMin`.
    pub shift: Option<f64>,
}

/// Applies `t` to `x` in place, elementwise, and returns a record of the
/// resolved parameters for sidecar provenance. Deterministic, no RNG.
///
/// - `Transform::None`: no-op; `x` is left untouched.
/// - `Transform::Exp2Clip { clamp }`: assumes `clamp > 0` (validated by
///   callers); each entry is clamped to `[-clamp, +clamp]` then exponentiated
///   with base 2.
/// - `Transform::ShiftMin`: computes `m = min(x)` over every entry of the
///   whole matrix, then subtracts `m` from every entry (so the transformed
///   matrix's minimum is exactly `0.0`).
///
/// # Panics
/// `Transform::ShiftMin` panics if `x` is empty (no elements to take a
/// minimum over) — the same precondition `nmf` itself already enforces on
/// its input.
pub fn apply_transform(x: &mut [Vec<f64>], t: &Transform) -> TransformRecord {
    match t {
        Transform::None => TransformRecord {
            name: "none".to_string(),
            clamp: None,
            shift: None,
        },
        Transform::Exp2Clip { clamp } => {
            for row in x.iter_mut() {
                for v in row.iter_mut() {
                    let clamped = v.max(-*clamp).min(*clamp);
                    *v = clamped.exp2();
                }
            }
            TransformRecord {
                name: "exp2-clip".to_string(),
                clamp: Some(*clamp),
                shift: None,
            }
        }
        Transform::ShiftMin => {
            let m = x
                .iter()
                .flat_map(|row| row.iter().copied())
                .fold(f64::INFINITY, f64::min);
            assert!(
                m.is_finite(),
                "apply_transform(ShiftMin): matrix has no finite entries"
            );
            for row in x.iter_mut() {
                for v in row.iter_mut() {
                    *v -= m;
                }
            }
            TransformRecord {
                name: "shift-min".to_string(),
                clamp: None,
                shift: Some(m),
            }
        }
    }
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
// Numeric kernel: the matrix dims (n, p, k) and solver params (max_iter, tol)
// are intrinsic to the update and clearer as flat args than bundled in a struct.
#[allow(clippy::too_many_arguments)]
fn frobenius_mu(
    x: &[Vec<f64>],
    w: &mut [Vec<f64>],
    h: &mut [Vec<f64>],
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
        w: w.to_vec(),
        h: h.to_vec(),
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
    // One reconstructed row at a time, accumulated over `a` in the OUTER loop.
    //
    // The obvious formulation computes `(0..k).map(|a| w[i][a] * h[a][j]).sum()` per (i, j),
    // which walks `h[a][j]` down a column. With `Vec<Vec<f64>>` every step of that walk is a
    // pointer chase into a different heap allocation, so each output element costs k cache
    // misses. This is the hottest loop in the Frobenius path — n·p·k multiply-adds, paid on
    // EVERY iteration for the convergence check — so its layout dominates the runtime.
    //
    // Accumulating into a row buffer makes both `h_a[j]` and `wh_row[j]` contiguous and lets
    // the inner loop vectorise. The additions over `a` still run in ascending order and the
    // residual accumulation still walks (i, j) row-major, so the sequence of floating-point
    // operations is unchanged and the result is bit-identical.
    let mut acc = 0.0_f64;
    let mut wh_row = vec![0.0_f64; p];
    for i in 0..n {
        wh_row.iter_mut().for_each(|v| *v = 0.0);
        let w_i = &w[i];
        for (a, h_a) in h.iter().enumerate().take(k) {
            let w_ia = w_i[a];
            for (dst, h_aj) in wh_row.iter_mut().zip(h_a.iter()).take(p) {
                *dst += w_ia * h_aj;
            }

        }
        let x_i = &x[i];
        for (j, wh) in wh_row.iter().enumerate().take(p) {
            let r = x_i[j] - wh;
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
// Numeric kernel: the matrix dims (n, p, k) and solver params (max_iter, tol)
// are intrinsic to the update and clearer as flat args than bundled in a struct.
#[allow(clippy::too_many_arguments)]
fn kl_mu(
    x: &[Vec<f64>],
    w: &mut [Vec<f64>],
    h: &mut [Vec<f64>],
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
            for &v in &h[kk] {
                h_row_sum[kk] += v;
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
        w: w.to_vec(),
        h: h.to_vec(),
        n_iter,
        final_error: prev_loss,
        converged,
    }
}

/// Generalised KL divergence D_KL(X || WH).
///
/// D = sum_ij [ X_ij * log(X_ij / (WH)_ij) - X_ij + (WH)_ij ]
/// Convention: 0 * log(0) = 0 (handled via the EPSILON floor on WH).
fn kl_loss(x: &[Vec<f64>], w: &[Vec<f64>], h: &[Vec<f64>], n: usize, p: usize, k: usize) -> f64 {
    let mut acc = 0.0_f64;
    for i in 0..n {
        for j in 0..p {
            let wh_ij: f64 = (0..k)
                .map(|kk| w[i][kk] * h[kk][j])
                .sum::<f64>()
                .max(EPSILON);
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
///
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
        if count > 0.0 {
            total / count
        } else {
            0.0
        }
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
        h[r][..p].copy_from_slice(&h_r[..p]);
    }

    // "a" variant: replace zeros with mean(X)
    for v in w.iter_mut().flat_map(|row| row.iter_mut()) {
        if *v == 0.0 {
            *v = mean_x;
        }
    }
    for v in h.iter_mut().flat_map(|row| row.iter_mut()) {
        if *v == 0.0 {
            *v = mean_x;
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
        if count > 0.0 {
            total / count
        } else {
            1.0
        }
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
            // `r` indexes rows of x at two distinct columns (i and j) per step.
            #[allow(clippy::needless_range_loop)]
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
fn mat_mul(
    rows_a: usize,
    shared: usize,
    a: &[Vec<f64>],
    _shared2: usize,
    cols_b: usize,
    b: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    // i-l-j, not i-j-l. The inner loop must walk `j`, because `b[l][j]` and `out[i][j]` are
    // then both contiguous; with `j` in the middle the inner loop walks `b[l][j]` down a
    // column, which on `Vec<Vec<f64>>` is a pointer chase into a different allocation per
    // element. `mat_mul_at_b` below already has this order, which is why it is not a
    // bottleneck and this one was.
    //
    // The additions into `out[i][j]` still run over `l` in ascending order, exactly as the
    // scalar accumulator did, so the sequence of floating-point operations is unchanged.
    let mut out = vec![vec![0.0_f64; cols_b]; rows_a];
    for i in 0..rows_a {
        let a_i = &a[i];
        let out_i = &mut out[i];
        for (l, b_l) in b.iter().enumerate().take(shared) {
            let a_il = a_i[l];
            for (dst, b_lj) in out_i.iter_mut().zip(b_l.iter()).take(cols_b) {
                *dst += a_il * b_lj;
            }
        }
    }
    out
}

/// A^T @ B: (cols_a × cols_b) from (rows × cols_a) and (rows × cols_b)
fn mat_mul_at_b(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    cols_a: usize,
    rows: usize,
    cols_b: usize,
) -> Vec<Vec<f64>> {
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
    // `r` indexes rows of a; the inner loops read a[r] at two columns (i and j).
    #[allow(clippy::needless_range_loop)]
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
fn mat_mul_a_bt(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    rows_a: usize,
    _shared: usize,
    rows_b: usize,
) -> Vec<Vec<f64>> {
    // Transpose B once, then accumulate along `l` into rows_b independent accumulators.
    //
    // The natural form is a dot product per (i, j): `acc += a[i][l] * b[j][l]` over `l`.
    // Both operands are contiguous, so locality is fine — but it is a floating-point
    // REDUCTION, and the compiler may not vectorise it because doing so would reassociate
    // the sum. What is left is a serial multiply-add chain whose loop-carried dependency is
    // several cycles per element, so with `shared` in the thousands this loop is
    // latency-bound rather than throughput-bound and becomes the hottest call in the
    // multiplicative-update path once the others are fixed.
    //
    // Accumulating into `out[i][0..rows_b]` instead gives rows_b independent dependency
    // chains, which hides that latency and lets the inner loop vectorise. Each output still
    // accumulates over `l` in ascending order, exactly as the scalar accumulator did, so the
    // sequence of floating-point operations per output element is unchanged.
    //
    // The transpose costs rows_b × shared writes once per call, against rows_a × rows_b ×
    // shared multiply-adds for the product itself — negligible at any shape where this
    // function is hot.
    let shared = a[0].len();
    let mut bt = vec![vec![0.0_f64; rows_b]; shared];
    for (j, b_j) in b.iter().enumerate().take(rows_b) {
        for (l, b_jl) in b_j.iter().enumerate().take(shared) {
            bt[l][j] = *b_jl;
        }
    }
    let mut out = vec![vec![0.0_f64; rows_b]; rows_a];
    for i in 0..rows_a {
        let a_i = &a[i];
        let out_i = &mut out[i];
        for (l, bt_l) in bt.iter().enumerate().take(shared) {
            let a_il = a_i[l];
            for (dst, bt_lj) in out_i.iter_mut().zip(bt_l.iter()).take(rows_b) {
                *dst += a_il * bt_lj;
            }
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
    let mut seed_loadings_per_ref: Vec<Vec<Vec<f64>>> =
        (0..k).map(|a| vec![reference[a].clone()]).collect();

    for alt in all_h.iter().skip(1) {
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

// ── k-selection ─────────────────────────────────────────────────────────────

/// Method for automatic k-selection in NMF.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KSelection {
    /// Use `nmf_cfg.k` as-is; no sweep is performed.
    Fixed,
    /// Select k by cophenetic-correlation knee (Brunet 2004).
    CopheneticKnee,
    /// Select k by reconstruction-error (RSS) elbow.
    RssKnee,
}

/// One row in the k-sweep output.
#[derive(Debug, Clone)]
pub struct KSweepRow {
    pub k: usize,
    /// Mean cophenetic correlation across seeds (consensus-matrix based).
    pub cophenetic: f64,
    /// Mean ||X - WH||_F^2 across seeds.
    pub mean_rss: f64,
    /// Mean KL loss across seeds (None when beta_loss == Frobenius).
    pub mean_kl: Option<f64>,
}

/// Result of automatic k-selection.
#[derive(Debug, Clone)]
pub struct KSelectionResult {
    pub selected_k: usize,
    pub per_k: Vec<KSweepRow>,
}

/// Run automatic k-selection for NMF.
///
/// For `KSelection::Fixed`: trivial — returns a single-row sweep at `nmf_cfg.k`.
/// For `KSelection::CopheneticKnee`: sweeps k ∈ [k_min, k_max], builds
///   the per-k consensus matrix (fraction of seeds where two samples co-cluster),
///   and selects the k with the highest cophenetic correlation (with tie-breaking
///   toward the left/smaller k).
/// For `KSelection::RssKnee`: sweeps k ∈ [k_min, k_max] and selects the k at
///   the maximum distance from the line connecting the first and last RSS points
///   (standard elbow / max-perpendicular-distance method).
pub fn select_k(
    x: &[Vec<f64>],
    nmf_cfg: &NmfConfig,
    ms_cfg: &MultiSeedConfig,
    k_min: usize,
    k_max: usize,
    method: KSelection,
) -> KSelectionResult {
    if method == KSelection::Fixed {
        // Trivial: run a single NMF at the configured k.
        let result = nmf(x, nmf_cfg);
        let rss = result.final_error * result.final_error;
        let mean_kl = if nmf_cfg.beta_loss == BetaLoss::KullbackLeibler {
            Some(result.final_error)
        } else {
            None
        };
        return KSelectionResult {
            selected_k: nmf_cfg.k,
            per_k: vec![KSweepRow {
                k: nmf_cfg.k,
                cophenetic: f64::NAN,
                mean_rss: rss,
                mean_kl,
            }],
        };
    }

    assert!(k_min >= 2, "k_min must be >= 2");
    assert!(k_max >= k_min, "k_max must be >= k_min");

    let n = x.len();
    let ks: Vec<usize> = (k_min..=k_max).collect();
    let mut rows: Vec<KSweepRow> = Vec::with_capacity(ks.len());

    for &k in &ks {
        let cfg_k = NmfConfig { k, ..*nmf_cfg };

        // Run n_seeds NMF runs, collect W matrices and errors.
        let mut rss_sum = 0.0_f64;
        let mut kl_sum = 0.0_f64;
        // co_occurrence[i][j] = number of seeds where sample i and j are assigned
        // to the same component (argmax of W row).
        let mut co_count = vec![vec![0u32; n]; n];

        for seed_i in 0..ms_cfg.n_seeds {
            let seed = ms_cfg.seed_base.wrapping_add(seed_i as u64);
            let cfg_s = NmfConfig { seed, ..cfg_k };
            let res = nmf(x, &cfg_s);

            // Compute per-seed RSS (Frobenius squared).
            let n_p = if x.is_empty() { 0 } else { x[0].len() };
            let mut rss = 0.0_f64;
            // i and j index the X / W / H matrices in lockstep (x[i][j], w[i][a], h[a][j]).
            #[allow(clippy::needless_range_loop)]
            for i in 0..n {
                for j in 0..n_p {
                    let wh: f64 = (0..k).map(|a| res.w[i][a] * res.h[a][j]).sum();
                    let d = x[i][j] - wh;
                    rss += d * d;
                }
            }
            rss_sum += rss;

            // KL: only meaningful if BetaLoss::KullbackLeibler.
            if nmf_cfg.beta_loss == BetaLoss::KullbackLeibler {
                kl_sum += res.final_error;
            }

            // Assign each sample to the component with the highest W value.
            let assignments: Vec<usize> = res
                .w
                .iter()
                .map(|row| {
                    row.iter()
                        .enumerate()
                        .max_by(|(_, a), (_, b)| {
                            a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
                        })
                        .map(|(idx, _)| idx)
                        .unwrap_or(0)
                })
                .collect();

            // Update co-occurrence.
            for i in 0..n {
                for j in i..n {
                    if assignments[i] == assignments[j] {
                        co_count[i][j] += 1;
                        co_count[j][i] += 1;
                    }
                }
            }
        }

        let n_seeds_f = ms_cfg.n_seeds as f64;
        let mean_rss = rss_sum / n_seeds_f;
        let mean_kl = if nmf_cfg.beta_loss == BetaLoss::KullbackLeibler {
            Some(kl_sum / n_seeds_f)
        } else {
            None
        };

        // Build consensus matrix C[i][j] = co_count[i][j] / n_seeds.
        // Cophenetic correlation: Spearman r between upper-triangle of (1 - C)
        // and upper-triangle of cophenetic distances from average-linkage
        // hierarchical clustering of (1 - C).
        let cophenetic = if n < 2 {
            1.0
        } else {
            let consensus_dist: Vec<Vec<f64>> = (0..n)
                .map(|i| {
                    (0..n)
                        .map(|j| 1.0 - co_count[i][j] as f64 / n_seeds_f)
                        .collect()
                })
                .collect();
            let cophenetic_dist = average_linkage_cophenetic(&consensus_dist, n);
            spearman_upper_tri(&consensus_dist, &cophenetic_dist, n)
        };

        rows.push(KSweepRow {
            k,
            cophenetic,
            mean_rss,
            mean_kl,
        });
    }

    let selected_k = match method {
        KSelection::Fixed => unreachable!(),
        KSelection::CopheneticKnee => {
            // The cophenetic curve typically decreases monotonically with k
            // (adding components splits clusters that were already coherent).
            // We want the "elbow": the k where the curve first begins to
            // plateau, i.e., the point of maximum perpendicular distance from
            // the line connecting the first and last cophenetic values.
            // This is the same max-perpendicular-distance elbow as the RSS
            // method, but applied to the cophenetic curve (descending).
            let n_pts = rows.len();
            if n_pts <= 1 {
                rows.first().map(|r| r.k).unwrap_or(k_min)
            } else {
                let x0 = 0.0_f64;
                let y0 = rows[0].cophenetic;
                let x1 = (n_pts - 1) as f64;
                let y1 = rows[n_pts - 1].cophenetic;

                let a = y1 - y0;
                let b = x0 - x1;
                let c = x1 * y0 - x0 * y1;
                let denom = (a * a + b * b).sqrt().max(1e-30);

                let mut best_k = rows[0].k;
                let mut best_dist = f64::NEG_INFINITY;
                for (idx, row) in rows.iter().enumerate() {
                    if row.cophenetic.is_nan() {
                        continue;
                    }
                    let xi = idx as f64;
                    let yi = row.cophenetic;
                    let dist = (a * xi + b * yi + c).abs() / denom;
                    if dist > best_dist {
                        best_dist = dist;
                        best_k = row.k;
                    }
                }
                best_k
            }
        }
        KSelection::RssKnee => {
            // Standard elbow: max perpendicular distance from the line
            // connecting the (k_min, rss_first) and (k_max, rss_last) points.
            let n_pts = rows.len();
            if n_pts <= 1 {
                rows.first().map(|r| r.k).unwrap_or(k_min)
            } else {
                let x0 = 0.0_f64;
                let y0 = rows[0].mean_rss;
                let x1 = (n_pts - 1) as f64;
                let y1 = rows[n_pts - 1].mean_rss;

                // Line from (x0, y0) to (x1, y1): ax + by + c = 0
                // a = y1 - y0, b = x0 - x1, c = x1*y0 - x0*y1
                let a = y1 - y0;
                let b = x0 - x1;
                let c = x1 * y0 - x0 * y1;
                let denom = (a * a + b * b).sqrt().max(1e-30);

                let mut best_k = rows[0].k;
                let mut best_dist = f64::NEG_INFINITY;
                for (idx, row) in rows.iter().enumerate() {
                    let xi = idx as f64;
                    let yi = row.mean_rss;
                    let dist = (a * xi + b * yi + c).abs() / denom;
                    if dist > best_dist {
                        best_dist = dist;
                        best_k = row.k;
                    }
                }
                best_k
            }
        }
    };

    KSelectionResult {
        selected_k,
        per_k: rows,
    }
}

/// Compute cophenetic distances from average-linkage hierarchical clustering
/// of a pre-computed distance matrix `d` (n × n, symmetric, zero diagonal).
///
/// Returns an n × n cophenetic distance matrix where `coph[i][j]` is the
/// height at which samples i and j first merge.
///
/// This is a minimal O(n^3) implementation (sufficient for n ≤ 200).
fn average_linkage_cophenetic(d: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    // We use a standard O(n^3) agglomerative algorithm with active-cluster tracking.
    // Each cluster stores the set of original sample indices it contains.
    // The distance between two clusters is the average of all pairwise distances
    // between their members (true average linkage).

    // coph[i][j] = merge height of the pair (i, j); initialized to 0.
    let mut coph = vec![vec![0.0_f64; n]; n];

    // Cluster membership: cluster k -> Vec of original indices.
    // Start with n singleton clusters.
    let mut clusters: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    // Active flag per cluster slot.
    let mut active: Vec<bool> = vec![true; n];

    // Current distances between active clusters (we update in place).
    let mut dist = d.to_vec();

    for _step in 0..(n - 1) {
        // Find the two active clusters with minimum distance.
        let mut best = f64::INFINITY;
        let mut best_i = 0;
        let mut best_j = 0;
        for i in 0..dist.len() {
            if !active[i] {
                continue;
            }
            for j in (i + 1)..dist.len() {
                if !active[j] {
                    continue;
                }
                if dist[i][j] < best {
                    best = dist[i][j];
                    best_i = i;
                    best_j = j;
                }
            }
        }

        let merge_height = best;

        // Record cophenetic height for all pairs (a, b) where a ∈ cluster_i, b ∈ cluster_j.
        for &a in &clusters[best_i] {
            for &b in &clusters[best_j] {
                coph[a][b] = merge_height;
                coph[b][a] = merge_height;
            }
        }

        // Merge cluster best_j into best_i.
        // New distance from merged cluster to any other active cluster c:
        //   avg-linkage: (|ci| * d(ci, c) + |cj| * d(cj, c)) / (|ci| + |cj|)
        let size_i = clusters[best_i].len() as f64;
        let size_j = clusters[best_j].len() as f64;
        let total = size_i + size_j;

        let n_slots = dist.len();
        for c in 0..n_slots {
            if !active[c] || c == best_i || c == best_j {
                continue;
            }
            let new_d = (size_i * dist[best_i][c] + size_j * dist[best_j][c]) / total;
            dist[best_i][c] = new_d;
            dist[c][best_i] = new_d;
        }

        // Move best_j members into best_i.
        let members_j = clusters[best_j].clone();
        clusters[best_i].extend(members_j);
        active[best_j] = false;
    }

    coph
}

/// Spearman correlation between the upper triangles of two n×n symmetric matrices.
fn spearman_upper_tri(a: &[Vec<f64>], b: &[Vec<f64>], n: usize) -> f64 {
    let m = n * (n - 1) / 2;
    if m == 0 {
        return 1.0;
    }
    let mut va = Vec::with_capacity(m);
    let mut vb = Vec::with_capacity(m);
    for i in 0..n {
        for j in (i + 1)..n {
            va.push(a[i][j]);
            vb.push(b[i][j]);
        }
    }
    pearson_on_ranks(&va, &vb)
}

/// Convert values to ranks (average ranks for ties), then compute Pearson r.
fn pearson_on_ranks(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    let ra = rank_vector(a);
    let rb = rank_vector(b);
    pearson_r(&ra, &rb)
}

/// Assign average ranks to a slice of f64 values.
fn rank_vector(v: &[f64]) -> Vec<f64> {
    let n = v.len();
    // Sort by value, keeping original indices.
    let mut indexed: Vec<(usize, f64)> = v.iter().enumerate().map(|(i, &x)| (i, x)).collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let mut ranks = vec![0.0_f64; n];
    let mut i = 0;
    while i < n {
        // Find the run of equal values.
        let mut j = i;
        while j < n && (indexed[j].1 - indexed[i].1).abs() < f64::EPSILON * 1e3 {
            j += 1;
        }
        // Average rank for this group (1-based).
        let avg_rank = (i + j + 1) as f64 / 2.0;
        for k in i..j {
            ranks[indexed[k].0] = avg_rank;
        }
        i = j;
    }
    ranks
}

/// Pearson r between two equal-length slices.
fn pearson_r(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len();
    if n < 2 {
        return 0.0;
    }
    let mean_a = a.iter().sum::<f64>() / n as f64;
    let mean_b = b.iter().sum::<f64>() / n as f64;
    let mut num = 0.0_f64;
    let mut sa = 0.0_f64;
    let mut sb = 0.0_f64;
    for i in 0..n {
        let da = a[i] - mean_a;
        let db = b[i] - mean_b;
        num += da * db;
        sa += da * da;
        sb += db * db;
    }
    let denom = (sa * sb).sqrt();
    if denom < 1e-30 {
        0.0
    } else {
        num / denom
    }
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
            // i and j index the X / W / H matrices in lockstep.
            #[allow(clippy::needless_range_loop)]
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
            .map(|i| {
                (0..4)
                    .map(|j| w_true[i][0] * h_true[0][j] + w_true[i][1] * h_true[1][j])
                    .collect()
            })
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
            // i and j index the X / W / H matrices in lockstep.
            #[allow(clippy::needless_range_loop)]
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
        let x = mat(&[&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0], &[7.0, 8.0, 9.0]]);
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
        let x = mat(&[&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0], &[7.0, 8.0, 9.0]]);
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
            // i and j index the X / W / H matrices in lockstep.
            #[allow(clippy::needless_range_loop)]
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
        let x = mat(&[&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0], &[7.0, 8.0, 9.0]]);
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

    // ── k-selection tests ────────────────────────────────────────────────────

    /// Build a rank-3 synthetic matrix: 20 samples × 10 features.
    fn make_rank3_synthetic() -> Vec<Vec<f64>> {
        let n = 20usize;
        let p = 10usize;
        let k_true = 3usize;

        let h_true: Vec<Vec<f64>> = vec![
            vec![5.0, 4.0, 3.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            vec![0.1, 0.1, 0.1, 5.0, 4.0, 3.0, 0.1, 0.1, 0.1, 0.1],
            vec![0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 5.0, 4.0, 3.0],
        ];
        let mut w_true = vec![vec![0.0_f64; k_true]; n];
        for i in 0..n {
            w_true[i][i % k_true] = 2.0 + (i as f64 * 0.1);
            w_true[i][(i + 1) % k_true] = 0.5;
        }
        (0..n)
            .map(|i| {
                (0..p)
                    .map(|j| {
                        (0..k_true)
                            .map(|a| w_true[i][a] * h_true[a][j])
                            .sum::<f64>()
                    })
                    .collect()
            })
            .collect()
    }

    #[test]
    fn select_k_cophenetic_picks_planted_k_on_synthetic_rank_3() {
        let x = make_rank3_synthetic();
        let nmf_cfg = NmfConfig {
            k: 3, // will be overridden per-k
            beta_loss: BetaLoss::Frobenius,
            init: Init::Random,
            max_iter: 600,
            tol: 1e-8,
            seed: 0,
        };
        let ms_cfg = MultiSeedConfig {
            n_seeds: 10,
            seed_base: 200,
            stability_metric: StabilityMetric::JaccardTopN,
            stability_top_n: 3,
        };

        let result = select_k(&x, &nmf_cfg, &ms_cfg, 2, 6, KSelection::CopheneticKnee);

        assert_eq!(
            result.per_k.len(),
            5,
            "per_k should have k_max - k_min + 1 = 5 rows"
        );
        assert_eq!(
            result.selected_k, 3,
            "cophenetic knee should select k=3 on rank-3 synthetic; got k={}, cophenetic values: {:?}",
            result.selected_k,
            result.per_k.iter().map(|r| (r.k, r.cophenetic)).collect::<Vec<_>>()
        );
    }

    #[test]
    fn select_k_rss_picks_planted_k_on_synthetic_rank_3() {
        let x = make_rank3_synthetic();
        let nmf_cfg = NmfConfig {
            k: 3,
            beta_loss: BetaLoss::Frobenius,
            init: Init::Random,
            max_iter: 600,
            tol: 1e-8,
            seed: 0,
        };
        let ms_cfg = MultiSeedConfig {
            n_seeds: 10,
            seed_base: 300,
            stability_metric: StabilityMetric::JaccardTopN,
            stability_top_n: 3,
        };

        let result = select_k(&x, &nmf_cfg, &ms_cfg, 2, 6, KSelection::RssKnee);

        // Knee detection has off-by-one ambiguity; accept k ∈ {3, 4}.
        assert!(
            result.selected_k == 3 || result.selected_k == 4,
            "RSS knee should select k ∈ {{3, 4}} on rank-3 synthetic; got k={}, rss values: {:?}",
            result.selected_k,
            result
                .per_k
                .iter()
                .map(|r| (r.k, r.mean_rss))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn select_k_fixed_passes_through() {
        let x = mat(&[
            &[3.0, 2.0, 0.1, 0.1],
            &[2.5, 3.0, 0.1, 0.1],
            &[0.1, 0.1, 4.0, 3.0],
            &[0.1, 0.1, 3.0, 4.5],
            &[1.5, 1.0, 2.0, 2.0],
        ]);
        let nmf_cfg = NmfConfig {
            k: 2,
            beta_loss: BetaLoss::Frobenius,
            init: Init::Random,
            max_iter: 200,
            tol: 1e-8,
            seed: 42,
        };
        let ms_cfg = MultiSeedConfig {
            n_seeds: 3,
            seed_base: 42,
            stability_metric: StabilityMetric::JaccardTopN,
            stability_top_n: 2,
        };

        let result = select_k(&x, &nmf_cfg, &ms_cfg, 2, 5, KSelection::Fixed);

        assert_eq!(
            result.selected_k, nmf_cfg.k,
            "Fixed should pass through cfg.k"
        );
        assert_eq!(
            result.per_k.len(),
            1,
            "Fixed should return exactly 1 per_k row"
        );
        assert_eq!(result.per_k[0].k, nmf_cfg.k);
    }

    // ── input transform tests (`--transform exp2-clip` / `shift-min`) ───────

    #[test]
    fn exp2_clip_transform_clamps_then_exponentiates() {
        let mut x = mat(&[&[-8.0, 0.0], &[1.0, 30.0]]);
        let record = apply_transform(&mut x, &Transform::Exp2Clip { clamp: 6.0 });

        assert_eq!(record.name, "exp2-clip");
        assert_eq!(record.clamp, Some(6.0));
        assert_eq!(record.shift, None);

        let expected = mat(&[&[2f64.powf(-6.0), 1.0], &[2.0, 2f64.powf(6.0)]]);
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    (x[i][j] - expected[i][j]).abs() < 1e-12,
                    "x[{i}][{j}] = {}, expected {}",
                    x[i][j],
                    expected[i][j]
                );
            }
        }
    }

    #[test]
    fn shift_min_transform_subtracts_global_min_and_records_it() {
        let mut x = mat(&[&[-2.0, 1.0]]);
        let record = apply_transform(&mut x, &Transform::ShiftMin);

        assert_eq!(record.name, "shift-min");
        assert_eq!(record.clamp, None);
        assert_eq!(record.shift, Some(-2.0));
        assert_eq!(x, mat(&[&[0.0, 3.0]]));
    }

    #[test]
    fn none_transform_is_a_no_op() {
        let mut x = mat(&[&[-1.0, 2.0], &[3.0, -4.0]]);
        let original = x.clone();
        let record = apply_transform(&mut x, &Transform::None);

        assert_eq!(record.name, "none");
        assert_eq!(record.clamp, None);
        assert_eq!(record.shift, None);
        assert_eq!(x, original, "Transform::None must not mutate its input");
    }
}
