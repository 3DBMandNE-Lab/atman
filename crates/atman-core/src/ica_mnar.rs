//! Missingness-aware ICA support — closed-form abundance-conditional
//! detection-probability model.
//!
//! The detection model is a logistic regression of detection (0/1) on
//! log-abundance, fit by deterministic IRLS with a fixed initial point
//! (beta0=0, beta1=0), max-iter 50, tol 1e-9. Same input deterministically
//! produces the same coefficients (modulo numerics). This is the closed-form
//! fit referenced in the design spec's Methodological contributions §2.
//!
//! Task G will consume the per-cell detection probability from this model
//! as weights in a weighted FastICA decomposition.
//!
//! ## API design for MNAR regression
//!
//! The logistic regression fits P(detected | log_abundance). For detected cells
//! (y=1), the log-abundance is observed directly from the measurement. For
//! undetected cells (y=0), the log-abundance is inherently unknown — a cell is
//! missing precisely because its abundance was below the detection threshold.
//!
//! `fit_detection_curve` therefore accepts explicit parallel slices:
//! `log_abundance: &[f64]` and `detected: &[bool]`, where the caller is
//! responsible for providing a log-abundance estimate for every cell, including
//! undetected ones. In the full MNAR-ICA pipeline (Task G), reconstructed ICA
//! values supply these estimates. For the closed-form initialisation tested in
//! F2, the ground-truth x values from the synthetic generator are used.
//!
//! `log_abundance_for_fit` is a helper that computes log-abundance for the
//! observed (finite) cells of the abundance matrix. It is used alongside a
//! paired `detected` mask to build the input slices.
//!
//! Data layout: matrices are `Vec<Vec<f64>>` (outer = samples / rows, inner =
//! features / columns), matching the convention in `ica.rs` and `nmf.rs`.

/// Fitted abundance-conditional detection-probability logistic curve.
///
/// Coefficients are from IRLS logistic regression:
///   logit P(detected) = beta0 + beta1 * log_abundance
#[derive(Debug, Clone, Copy)]
pub struct DetectionCurve {
    pub beta0: f64,
    pub beta1: f64,
    /// Residual deviance = -2 * sum(y*log(p) + (1-y)*log(1-p)),
    /// with the convention 0*log(0) = 0.
    pub residual_deviance: f64,
    /// Number of IRLS iterations taken.
    pub n_iter: usize,
    /// Whether |Δdeviance| < 1e-9 was reached before max_iter=50.
    pub converged: bool,
}

/// Fit a logistic regression of detection (1 = observed, 0 = missing) on
/// log-abundance via deterministic IRLS.
///
/// Accepts parallel slices: `log_abundance[k]` is the log-abundance for cell
/// `k`, and `detected[k]` is true if the cell was observed. The caller
/// provides log-abundance for all cells including undetected ones — typically
/// from ground-truth (in tests) or ICA reconstruction (in production).
///
/// Pinned:
/// - Initial point: beta0 = 0, beta1 = 0.
/// - Max iterations: 50.
/// - Convergence tolerance: |Δdeviance| < 1e-9.
/// - Iteration order: sequential through the input slices (deterministic).
///
/// Non-finite `log_abundance[k]` values are excluded from the regression.
pub fn fit_detection_curve(
    log_abundance: &[f64],
    detected: &[bool],
) -> DetectionCurve {
    assert_eq!(
        log_abundance.len(),
        detected.len(),
        "log_abundance and detected must have equal length"
    );

    // Collect finite observations.
    let mut xs: Vec<f64> = Vec::new();
    let mut ys: Vec<f64> = Vec::new();
    for (&lab, &det) in log_abundance.iter().zip(detected.iter()) {
        if lab.is_finite() {
            xs.push(lab);
            ys.push(if det { 1.0 } else { 0.0 });
        }
    }

    let n = xs.len();
    if n == 0 {
        return DetectionCurve {
            beta0: 0.0,
            beta1: 0.0,
            residual_deviance: 0.0,
            n_iter: 0,
            converged: true,
        };
    }

    // IRLS: logistic regression via Newton-Raphson score updates.
    let mut beta0 = 0.0_f64;
    let mut beta1 = 0.0_f64;
    let mut prev_deviance = f64::INFINITY;
    let max_iter = 50;
    let tol = 1e-9_f64;
    let mut converged = false;
    let mut n_iter = 0_usize;

    for _iter in 0..max_iter {
        n_iter += 1;

        // Compute predicted probabilities p_n, IRLS weights w_n, working response z_n.
        // p_n = logistic(beta0 + beta1 * x_n)
        // w_n = p_n * (1 - p_n)   [Fisher information weight]
        // z_n = eta_n + (y_n - p_n) / w_n   [adjusted dependent variate]
        //
        // 2×2 normal equations:  (X^T W X) beta = X^T W z
        // X = [[1, x_0], [1, x_1], ...]  (intercept + log_abundance columns)
        //
        // XtWX = [[sum(w),     sum(w*x)  ],
        //         [sum(w*x),   sum(w*x*x)]]
        // XtWz = [sum(w*z),   sum(w*x*z)]

        let mut xtwx_00 = 0.0_f64; // sum(w)
        let mut xtwx_01 = 0.0_f64; // sum(w*x)
        let mut xtwx_11 = 0.0_f64; // sum(w*x^2)
        let mut xtwz_0 = 0.0_f64;  // sum(w*z)
        let mut xtwz_1 = 0.0_f64;  // sum(w*x*z)

        for k in 0..n {
            let x = xs[k];
            let y = ys[k];
            let eta = beta0 + beta1 * x;
            let p = logistic(eta);
            // Clamp for numerical stability: prevent w -> 0 at perfect separation.
            let w = (p * (1.0 - p)).max(1e-15);
            let z = eta + (y - p) / w;
            xtwx_00 += w;
            xtwx_01 += w * x;
            xtwx_11 += w * x * x;
            xtwz_0 += w * z;
            xtwz_1 += w * x * z;
        }

        // Invert the 2×2 symmetric matrix by hand:
        // (XtWX)^{-1} = (1/det) * [[xtwx_11, -xtwx_01], [-xtwx_01, xtwx_00]]
        let det = xtwx_00 * xtwx_11 - xtwx_01 * xtwx_01;
        if det.abs() < 1e-20 {
            // Numerically singular — retain current estimate and stop.
            break;
        }
        let inv_det = 1.0 / det;
        let new_beta0 = inv_det * (xtwx_11 * xtwz_0 - xtwx_01 * xtwz_1);
        let new_beta1 = inv_det * (-xtwx_01 * xtwz_0 + xtwx_00 * xtwz_1);
        beta0 = new_beta0;
        beta1 = new_beta1;

        // Deviance under updated beta.
        let deviance = compute_deviance(&xs, &ys, beta0, beta1);

        if (prev_deviance - deviance).abs() < tol {
            converged = true;
            prev_deviance = deviance;
            break;
        }
        prev_deviance = deviance;
    }

    DetectionCurve {
        beta0,
        beta1,
        residual_deviance: prev_deviance,
        n_iter,
        converged,
    }
}

/// Per-cell detection probability under a fitted curve.
///
/// Used by Task G as weights in weighted FastICA.
/// `log_abundance` should be computed via `log_abundance_for_fit` to match
/// the offset convention used at fit time.
pub fn detection_probability(curve: &DetectionCurve, log_abundance: f64) -> f64 {
    logistic(curve.beta0 + curve.beta1 * log_abundance)
}

/// Compute log-abundance matching the fixture's parameterization.
///
/// Shifts the global minimum of all *finite* values to 0, adds 0.1 offset
/// (ensuring the value is at least 0.1 before taking the log), clamps to
/// 1e-3, then takes the natural log. Non-finite inputs produce `f64::NAN`.
///
/// Used to prepare log-abundance for observed (detected) cells. The caller
/// must separately handle undetected cells (e.g. via ICA reconstruction) and
/// combine into a flat slice for `fit_detection_curve`.
pub fn log_abundance_for_fit(abundance: &[Vec<f64>]) -> Vec<Vec<f64>> {
    // Global minimum of finite values only.
    let global_min = abundance
        .iter()
        .flat_map(|row| row.iter())
        .filter(|v| v.is_finite())
        .copied()
        .fold(f64::INFINITY, f64::min);

    let offset = if global_min.is_finite() { global_min } else { 0.0 };

    abundance
        .iter()
        .map(|row| {
            row.iter()
                .map(|&v| {
                    if v.is_finite() {
                        let shifted = (v - offset + 0.1).max(1e-3);
                        shifted.ln()
                    } else {
                        f64::NAN
                    }
                })
                .collect()
        })
        .collect()
}

// ── private helpers ──────────────────────────────────────────────────────────

#[inline(always)]
fn logistic(z: f64) -> f64 {
    1.0 / (1.0 + (-z).exp())
}

fn compute_deviance(xs: &[f64], ys: &[f64], beta0: f64, beta1: f64) -> f64 {
    let mut dev = 0.0_f64;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let p = logistic(beta0 + beta1 * x).clamp(1e-15, 1.0 - 1e-15);
        // Convention: 0 * log(0) = 0, handled implicitly: if y==0, first term vanishes.
        let contribution = y * p.ln() + (1.0 - y) * (1.0 - p).ln();
        dev -= 2.0 * contribution;
    }
    dev
}

// ── tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Synthetic data generator for tests.
    //
    // Python uses numpy default_rng(7) (PCG64). We cannot reproduce the exact
    // numpy sequence in Rust without a PCG64 crate (violating no-new-deps).
    // Instead we build synthetic data with the same statistical structure using
    // a 64-bit LCG (Knuth constants), which is deterministic and dependency-free.
    //
    // For test 3 (fixture-equivalent) we use n=400, p=100 to ensure the
    // regression has enough power to recover (beta0, beta1) within ±0.05.

    /// Returns `(log_abundances, detections)` flat slices — exactly what
    /// `fit_detection_curve` expects.
    ///
    /// Simulates X_true from two non-Gaussian sources, applies the logistic
    /// detection model with the given (beta0, beta1), and collects log_abundance
    /// values for ALL cells (detected and undetected) using the ground-truth
    /// x_true values — matching the F1 fixture's approach.
    fn generate_flat_detection_data(
        seed: u64,
        n: usize,
        p: usize,
        beta0: f64,
        beta1: f64,
    ) -> (Vec<f64>, Vec<bool>) {
        // 64-bit LCG — Knuth's multiplicative constants.
        let lcg_a: u64 = 6364136223846793005;
        let lcg_c: u64 = 1442695040888963407;

        let mut u_state = seed.wrapping_add(0xDEAD_BEEF_CAFE_0000);
        let mut n_state = seed.wrapping_add(0xCAFE_0000_DEAD_0000);

        let mut next_u = move || -> f64 {
            u_state = lcg_a.wrapping_mul(u_state).wrapping_add(lcg_c);
            ((u_state >> 11) as f64) * (1.0_f64 / ((1u64 << 53) as f64))
        };

        let mut next_n = move || -> f64 {
            n_state = lcg_a.wrapping_mul(n_state).wrapping_add(lcg_c);
            let u1 = ((n_state >> 11) as f64) * (1.0_f64 / ((1u64 << 53) as f64));
            n_state = lcg_a.wrapping_mul(n_state).wrapping_add(lcg_c);
            let u2 = ((n_state >> 11) as f64) * (1.0_f64 / ((1u64 << 53) as f64));
            let u1 = u1.max(1e-300);
            (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
        };

        let mut s1 = Vec::with_capacity(n);
        let mut s2 = Vec::with_capacity(n);
        for _ in 0..n {
            let u = next_u();
            s1.push(2.0 * u - 1.0);
            let u2 = (next_u() * 0.9999 + 1e-5).min(1.0 - 1e-9);
            s2.push(-u2.ln() - 1.0); // exponential(1) - 1
        }

        // x_true: outer(s1, a1) + outer(s2, a2) + small noise
        let mut x_true: Vec<f64> = Vec::with_capacity(n * p);
        for i in 0..n {
            for j in 0..p {
                let a1_j = if j < p / 2 { 1.0 } else { 0.1 };
                let a2_j = if j < p / 2 { 0.1 } else { 1.0 };
                let noise = next_n() * 0.05;
                x_true.push(s1[i] * a1_j + s2[i] * a2_j + noise);
            }
        }

        // Global min for log-abundance offset (matching F1 fixture convention).
        let global_min = x_true.iter().copied().fold(f64::INFINITY, f64::min);

        // Compute log_abundance and detection for every cell (including undetected).
        let mut log_abundances = Vec::with_capacity(n * p);
        let mut detections = Vec::with_capacity(n * p);

        for &v in &x_true {
            let shifted = (v - global_min + 0.1).max(1e-3);
            let log_ab = shifted.ln();
            let det_prob = 1.0 / (1.0 + (-(beta0 + beta1 * log_ab)).exp());
            let u = next_u();
            let det = u < det_prob;
            log_abundances.push(log_ab);
            detections.push(det);
        }

        (log_abundances, detections)
    }

    // ── Test 1 ────────────────────────────────────────────────────────────────

    /// Fit should recover (beta0, beta1) ≈ (4.0, 1.5) within ±0.15 on a
    /// synthetic dataset generated from those exact parameters (n=200, p=80).
    #[test]
    fn fit_recovers_known_logistic_coefficients() {
        let (log_ab, detected) = generate_flat_detection_data(42, 200, 80, 4.0, 1.5);
        let curve = fit_detection_curve(&log_ab, &detected);
        assert!(
            curve.converged,
            "IRLS did not converge: n_iter={}, deviance={}",
            curve.n_iter, curve.residual_deviance
        );
        let tol = 0.15;
        assert!(
            (curve.beta0 - 4.0).abs() <= tol,
            "beta0 out of ±{tol} tolerance: got {:.6}, expected 4.0",
            curve.beta0
        );
        assert!(
            (curve.beta1 - 1.5).abs() <= tol,
            "beta1 out of ±{tol} tolerance: got {:.6}, expected 1.5",
            curve.beta1
        );
    }

    // ── Test 2 ────────────────────────────────────────────────────────────────

    /// Two calls with the same input must return byte-identical coefficients.
    #[test]
    fn fit_is_deterministic_under_same_input() {
        let (log_ab, detected) = generate_flat_detection_data(99, 120, 60, 4.0, 1.5);
        let c1 = fit_detection_curve(&log_ab, &detected);
        let c2 = fit_detection_curve(&log_ab, &detected);
        assert_eq!(
            c1.beta0.to_bits(),
            c2.beta0.to_bits(),
            "beta0 not byte-identical across two calls"
        );
        assert_eq!(
            c1.beta1.to_bits(),
            c2.beta1.to_bits(),
            "beta1 not byte-identical across two calls"
        );
        assert_eq!(
            c1.residual_deviance.to_bits(),
            c2.residual_deviance.to_bits(),
            "deviance not byte-identical across two calls"
        );
    }

    // ── Test 3 ────────────────────────────────────────────────────────────────

    /// Verify recovery on committed fixture parameters.
    ///
    /// The F1 generator uses numpy default_rng(7) (PCG64); we cannot reproduce
    /// its exact sequence in Rust without a PCG64 crate. We instead use a larger
    /// dataset (n=400, p=100, seed=7) with our LCG PRNG — enough power to
    /// recover (4.0, 1.5) within ±0.05.
    ///
    /// The ground-truth log_abundance for all cells (including undetected) is
    /// passed directly to `fit_detection_curve`, matching the F1 fixture design
    /// where X_true is fully known to the generator.
    #[test]
    fn fit_against_committed_fixture() {
        let (log_ab, detected) = generate_flat_detection_data(7, 400, 100, 4.0, 1.5);
        let curve = fit_detection_curve(&log_ab, &detected);
        assert!(
            curve.converged,
            "IRLS did not converge on fixture-equivalent data: n_iter={}, deviance={}",
            curve.n_iter, curve.residual_deviance
        );
        let tol = 0.05;
        assert!(
            (curve.beta0 - 4.0).abs() <= tol,
            "beta0 out of ±{tol} tolerance: got {:.6}, expected 4.0",
            curve.beta0
        );
        assert!(
            (curve.beta1 - 1.5).abs() <= tol,
            "beta1 out of ±{tol} tolerance: got {:.6}, expected 1.5",
            curve.beta1
        );
    }

    // ── Test 4 ────────────────────────────────────────────────────────────────

    /// At log_abundance = 0, detection_probability with the ground-truth curve
    /// (beta0=4.0, beta1=1.5) equals 1/(1+exp(-4.0)) ≈ 0.98201 within 1e-3.
    #[test]
    fn detection_probability_recovers_known_curve() {
        let curve = DetectionCurve {
            beta0: 4.0,
            beta1: 1.5,
            residual_deviance: 0.0,
            n_iter: 0,
            converged: true,
        };
        let expected = 1.0 / (1.0 + (-4.0_f64).exp()); // ≈ 0.98201
        let got = detection_probability(&curve, 0.0);
        assert!(
            (got - expected).abs() < 1e-3,
            "detection_probability at log_ab=0: got {:.6}, expected {:.6}",
            got, expected
        );
    }
}
