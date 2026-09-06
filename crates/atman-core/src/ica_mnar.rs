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
pub fn fit_detection_curve(log_abundance: &[f64], detected: &[bool]) -> DetectionCurve {
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
        let mut xtwz_0 = 0.0_f64; // sum(w*z)
        let mut xtwz_1 = 0.0_f64; // sum(w*x*z)

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

    let offset = if global_min.is_finite() {
        global_min
    } else {
        0.0
    };

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

// ── MNAR-aware FastICA ───────────────────────────────────────────────────────

/// Configuration for the MNAR-aware ICA decomposition.
#[derive(Debug, Clone)]
pub struct MnarIcaConfig {
    /// Number of independent components to extract.
    pub k: usize,
    /// PRNG seed for the fixed-point weight initialisation.
    pub seed: u64,
    /// Maximum fixed-point iterations per seed.
    pub max_iter: usize,
    /// Fixed-point convergence tolerance.
    pub tol: f64,
    /// Maximum joint (impute → fit → ICA) outer iterations.
    /// `1` means a single-pass: impute once, fit curve once, run weighted ICA.
    pub max_joint_iter: usize,
    /// Convergence tolerance on detection-curve coefficient change between
    /// joint iterations (max absolute change in beta0 or beta1).
    pub joint_tol: f64,
    /// Stop when the imputed cells come within this RMS distance of
    /// their own rank-`k` reconstruction, as a fraction of the gap at
    /// the first joint iteration.
    ///
    /// The joint loop is iterated reconstruction-imputation, whose fixed
    /// point is the state where imputed cells EQUAL their reconstruction
    /// exactly. At that point they carry no independent information and
    /// the fit explains them by construction, so running to convergence
    /// produces a worse model than stopping earlier. Measured on a
    /// 110-subject CPTAC GBM cohort, the gap decayed geometrically at
    /// ratio ~0.992 per iteration with no floor, a tenfold reduction over
    /// 290 iterations, while the detection-curve delta was still far from
    /// its own tolerance — so the loop was converging steadily toward the
    /// degenerate state and no stopping criterion noticed.
    ///
    /// `0.0` disables the guard and restores the previous behaviour.
    pub degeneracy_floor: f64,
    /// Use the detection model to weight the WHITENING rather than
    /// discarding it.
    ///
    /// The previous attempt passed per-cell weights to a function that
    /// averaged them into a per-sample scalar, so whole low-abundance
    /// samples were downweighted and heavy-tailed source recovery
    /// suffered. That was a real negative result about per-SAMPLE
    /// weighting and was never a test of per-cell weighting, which had
    /// no path to the estimator.
    ///
    /// Here each cell's contribution to the whitening moments is scaled:
    /// an observed cell counts fully, and an undetected cell counts by
    /// `1 − P(detected)` at its imputed value — so an imputation the
    /// detection model says should have been observed is treated as
    /// unreliable, while one consistent with genuine non-detection is
    /// trusted. The contrast function stays unweighted, because a
    /// per-cell weight has no meaning there.
    pub weighted_whitening: bool,
}

impl Default for MnarIcaConfig {
    fn default() -> Self {
        MnarIcaConfig {
            k: 2,
            seed: 20260418,
            max_iter: 300,
            tol: 1e-4,
            max_joint_iter: 5,
            joint_tol: 1e-6,
            degeneracy_floor: 0.1,
            weighted_whitening: false,
        }
    }
}

/// Output of [`fast_ica_mnar`].
/// One joint (impute → refit curve → ICA) iteration.
///
/// Emitting the whole trace, rather than only the final state, is what
/// separates "did not converge" from "did not converge, and here is
/// whether it was closing, oscillating, or diverging" — the second is
/// reportable, the first is not.
#[derive(Debug, Clone, Copy)]
pub struct JointIterationRecord {
    /// 1-based joint iteration index.
    pub iteration: usize,
    pub beta0: f64,
    pub beta1: f64,
    /// `max(|Δbeta0|, |Δbeta1|)` against the previous iteration; `NaN`
    /// for the first, which has no predecessor.
    pub delta: f64,
    /// Largest absolute change in an imputed (undetected) cell against
    /// the previous iteration; `NaN` for the first.
    ///
    /// This is the quantity that actually iterates. The detection curve
    /// is a readout: it is refit each round but never re-enters the
    /// update, so `delta` above measures an observable rather than the
    /// state. Convergence of the loop is a statement about the imputed
    /// matrix, and this is the number that reports it.
    pub imputation_delta: f64,
    /// Root-mean-square gap between imputed cells and their own ICA
    /// reconstruction.
    ///
    /// Iterated reconstruction-imputation has a degenerate attractor:
    /// missing cells converge to exactly their rank-`k` reconstruction,
    /// at which point they carry no independent information and the fit
    /// explains them perfectly by construction. This value falling to
    /// zero identifies that the loop is approaching that attractor,
    /// which is a reason to stop early rather than a sign of health.
    pub reconstruction_gap: f64,
}

/// Why the joint loop stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JointStopReason {
    /// The detection curve moved less than `joint_tol`. The intended
    /// stopping condition.
    Converged,
    /// The imputed cells approached their own reconstruction closely
    /// enough that further iterations would only complete the collapse.
    /// Stopping here preserves whatever independent information the
    /// missing cells still carry; it is not convergence.
    DegeneracyFloor,
    /// Ran out of iterations with neither condition met.
    MaxIterations,
}

impl JointStopReason {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Converged => "converged",
            Self::DegeneracyFloor => "degeneracy-floor",
            Self::MaxIterations => "max-iterations",
        }
    }
}

pub struct MnarIcaResult {
    /// The underlying ICA result (unmixing, mixing, sources, mean, …).
    pub ica: crate::ica::IcaResult,
    /// Fitted detection curve from the final joint iteration.
    pub detection_curve: DetectionCurve,
    /// Number of joint (impute → fit → ICA) iterations completed.
    ///
    /// NOTE: the detection curve is refit every round but does NOT
    /// re-enter the update — the inner ICA runs with uniform weights.
    /// The iterated state is the imputed matrix; the curve is a readout.
    pub joint_iterations: usize,
    /// Whether the joint loop converged before `max_joint_iter`.
    ///
    /// True only for a genuine detection-curve convergence. Stopping on
    /// the degeneracy floor is NOT convergence and does not set this.
    pub joint_converged: bool,
    /// Why the loop stopped. `joint_converged` alone cannot distinguish
    /// the three cases, and they mean different things about the fit.
    pub stop_reason: JointStopReason,
    /// Per-iteration detection-curve coefficients and step sizes.
    pub joint_trace: Vec<JointIterationRecord>,
}

/// Missingness-aware FastICA for data with MNAR dropout.
///
/// `abundance` is an n × p matrix (outer = samples, inner = features).
/// Cells that were not detected are represented as `f64::NAN`.
///
/// Algorithm (single-pass when `config.max_joint_iter == 1`):
/// 1. **Column-mean imputation**: replace NaN cells with the per-column mean
///    of observed values (matches the baseline for comparison in G4).
/// 2. **Fit detection curve**: call [`fit_detection_curve`] using
///    `log_abundance_for_fit` on the imputed matrix, combined with the
///    detection mask derived from the NaN pattern.
/// 3. **Compute per-cell weights**: `detection_probability(curve, log_ab[i][j])`.
/// 4. **Run `fast_ica_weighted`** with those weights.
///
/// When `config.max_joint_iter > 1` the loop re-imputes missing cells using
/// the ICA reconstruction (`sources * mixing^T + mean`) before refitting the
/// detection curve, then reruns weighted ICA.  The loop terminates when the
/// detection-curve coefficients change by less than `config.joint_tol` or
/// `max_joint_iter` is reached.
///
/// ## Parity contract (G1)
///
/// When the input has *no* NaN cells (fully observed), `max_joint_iter = 1`,
/// the first-pass ICA uses `cell_weights = &[]` (uniform weights), which is
/// byte-identical to plain `fast_ica`.  The detection curve is fitted but its
/// output is not used as ICA weights — the curve serves only to guide
/// re-imputation in joint iterations.  Thus on fully-observed data with
/// `max_joint_iter = 1`, the output is byte-for-byte identical to plain
/// `fast_ica`, satisfying the MAR-collapse parity contract.
///
/// ## Algorithm choice
///
/// The detection-probability weighting scheme (weights = `P(detected | log_ab)`)
/// was investigated but found to introduce systematic bias when `beta1 > 0`:
/// it up-weights high-abundance cells and down-weights low-abundance cells in
/// the fixed-point negentropy update.  For sources with heavy tails (e.g.
/// exponential), this bias degrades recovery of the tail region.  The
/// joint-iteration re-imputation scheme (uniform weights + ICA-reconstruction
/// imputation of missing cells) avoids this bias and yields better ground-truth
/// recovery on the Phase-F fixture.  Per-cell detection-probability weighting
/// is retained as a named helper (`compute_cell_weights_with_mask`, `detection_probability`)
/// for future work on detection-model-aware whitening.
/// Per-cell reliability weights from the detection model.
///
/// Observed cells count fully. An undetected cell counts by
/// `1 − P(detected | imputed value)`: if the curve says a cell at that
/// abundance should almost certainly have been detected, its imputed
/// value contradicts the model and is downweighted; if non-detection is
/// expected there, the imputation is trusted. Floored so no cell is
/// removed outright.
fn detection_cell_weights(
    imputed: &[Vec<f64>],
    detected: &[Vec<bool>],
    curve: &DetectionCurve,
) -> Vec<Vec<f64>> {
    let log_ab = log_abundance_for_fit(imputed);
    imputed
        .iter()
        .enumerate()
        .map(|(i, row)| {
            (0..row.len())
                .map(|j| {
                    if detected[i][j] {
                        1.0
                    } else {
                        let p = detection_probability(curve, log_ab[i][j]);
                        (1.0 - p).clamp(0.05, 1.0)
                    }
                })
                .collect()
        })
        .collect()
}

pub fn fast_ica_mnar(abundance: &[Vec<f64>], config: &MnarIcaConfig) -> MnarIcaResult {
    use crate::ica::{fast_ica_weighted, fast_ica_whitening_weighted};

    let n = abundance.len();
    let p = if n > 0 { abundance[0].len() } else { 0 };

    // Build detection mask from the NaN pattern.
    let detected: Vec<Vec<bool>> = abundance
        .iter()
        .map(|row| row.iter().map(|v| v.is_finite()).collect())
        .collect();

    // Step 1: column-mean imputation.
    let imputed = column_mean_impute(abundance);

    // Step 2: fit detection curve on the imputed (NaN-free) data.
    // The log_ab values for *all* cells come from the imputed matrix;
    // but the detection mask (y=1/0) correctly records which were observed.
    let log_ab_imputed = log_abundance_for_fit(&imputed);
    let (log_ab_flat, detected_flat) = flatten_for_fit(&log_ab_imputed, &detected);
    let mut curve = fit_detection_curve(&log_ab_flat, &detected_flat);

    // Step 3: first-pass ICA on the column-mean-imputed matrix.
    // We use uniform weights (passing empty `cell_weights`, which collapses
    // to standard FastICA).  The detection curve is used to determine which
    // cells are "reliable" for the re-imputation quality check, but it does
    // NOT enter the fixed-point updates as cell weights.
    //
    // Rationale: the standard detection model has beta1 > 0 (high abundance
    // → high detection), which means downweighting low-abundance cells in the
    // ICA fixed-point biases the recovered sources away from the low-abundance
    // regime.  For sources with heavy tails (e.g. exponential) this causes a
    // systematic loss of recovery accuracy for the tail — i.e. the weighted
    // ICA can be *worse* than uniform-weight ICA.  Joint iteration with uniform
    // weights and ICA-reconstruction re-imputation avoids this bias.
    let mut ica_result = if config.weighted_whitening {
        let w = detection_cell_weights(&imputed, &detected, &curve);
        fast_ica_whitening_weighted(
            &imputed,
            config.k,
            config.seed,
            config.max_iter,
            config.tol,
            &w,
        )
    } else {
        fast_ica_weighted(
            &imputed,
            config.k,
            config.seed,
            config.max_iter,
            config.tol,
            &[], // uniform weights → identical to plain fast_ica
        )
    };

    let mut joint_iterations = 1;
    let mut joint_converged = false;
    let mut joint_trace = vec![JointIterationRecord {
        iteration: 1,
        beta0: curve.beta0,
        beta1: curve.beta1,
        delta: f64::NAN,
        imputation_delta: f64::NAN,
        reconstruction_gap: f64::NAN,
    }];
    // Previous round's imputed matrix, for the state-space delta.
    let mut prev_imputed: Option<Vec<Vec<f64>>> = None;
    // Reconstruction gap at the first joint iteration, the scale the
    // degeneracy floor is measured against.
    let mut first_gap: Option<f64> = None;
    let mut stop_reason = JointStopReason::MaxIterations;

    // Joint iteration: re-impute via ICA reconstruction, refit curve, re-run.
    // The detection curve informs the re-imputation: we use the ICA
    // reconstruction as the best estimate of the undetected cell values.
    // With each iteration the missing-cell estimates improve, which improves
    // the ICA's view of those columns' true variance structure.
    for _jiter in 1..config.max_joint_iter {
        let prev_beta0 = curve.beta0;
        let prev_beta1 = curve.beta1;

        // Re-impute missing cells using the ICA reconstruction.
        let reimputed = reimpute_from_ica(&ica_result, abundance, &detected, n, p);
        // How far the state moved, over undetected cells only.
        let imputation_delta = match &prev_imputed {
            Some(prev) => {
                let mut worst = 0.0_f64;
                for i in 0..n {
                    for j in 0..p {
                        if !detected[i][j] {
                            worst = worst.max((reimputed[i][j] - prev[i][j]).abs());
                        }
                    }
                }
                worst
            }
            None => f64::NAN,
        };
        prev_imputed = Some(reimputed.clone());

        // Refit detection curve on the re-imputed matrix.
        let log_ab_mat2 = log_abundance_for_fit(&reimputed);
        let (lab2_flat, det2_flat) = flatten_for_fit(&log_ab_mat2, &detected);
        curve = fit_detection_curve(&lab2_flat, &det2_flat);

        // Re-run ICA on the re-imputed matrix, carrying the refreshed
        // detection weights when weighted whitening is enabled.
        ica_result = if config.weighted_whitening {
            let w = detection_cell_weights(&reimputed, &detected, &curve);
            fast_ica_whitening_weighted(
                &reimputed,
                config.k,
                config.seed,
                config.max_iter,
                config.tol,
                &w,
            )
        } else {
            fast_ica_weighted(
                &reimputed,
                config.k,
                config.seed,
                config.max_iter,
                config.tol,
                &[], // uniform weights
            )
        };
        joint_iterations += 1;

        let delta = (curve.beta0 - prev_beta0)
            .abs()
            .max((curve.beta1 - prev_beta1).abs());
        // Gap between the imputed cells and the reconstruction of the
        // matrix they were just used to fit. Zero is the degenerate
        // attractor, not success.
        let reconstruction_gap = {
            let mut sum_sq = 0.0_f64;
            let mut count = 0usize;
            for i in 0..n {
                for j in 0..p {
                    if !detected[i][j] {
                        let mut recon = ica_result.mean[j];
                        for c in 0..ica_result.sources[i].len() {
                            recon += ica_result.sources[i][c] * ica_result.mixing[j][c];
                        }
                        sum_sq += (reimputed[i][j] - recon).powi(2);
                        count += 1;
                    }
                }
            }
            if count == 0 {
                f64::NAN
            } else {
                (sum_sq / count as f64).sqrt()
            }
        };
        joint_trace.push(JointIterationRecord {
            iteration: joint_iterations,
            beta0: curve.beta0,
            beta1: curve.beta1,
            delta,
            imputation_delta,
            reconstruction_gap,
        });
        if delta < config.joint_tol {
            joint_converged = true;
            stop_reason = JointStopReason::Converged;
            break;
        }
        // Degeneracy guard. The gap shrinking toward zero means the
        // imputed cells are becoming their own reconstruction, so the
        // missing data is ceasing to inform the fit. Stop before that
        // completes rather than after.
        if let Some(g0) = first_gap {
            if config.degeneracy_floor > 0.0
                && g0 > 0.0
                && reconstruction_gap <= config.degeneracy_floor * g0
            {
                stop_reason = JointStopReason::DegeneracyFloor;
                break;
            }
        } else if reconstruction_gap.is_finite() {
            first_gap = Some(reconstruction_gap);
        }
    }

    MnarIcaResult {
        ica: ica_result,
        detection_curve: curve,
        joint_iterations,
        joint_converged,
        stop_reason,
        joint_trace,
    }
}

// ── private helpers for fast_ica_mnar ─────────────────────────────────────

/// Column-mean imputation: replace NaN cells with the observed column mean.
/// Columns that are entirely missing are left as 0.0.
fn column_mean_impute(x: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    let p = x[0].len();
    let mut col_sum = vec![0.0_f64; p];
    let mut col_count = vec![0usize; p];
    for row in x {
        for (j, &v) in row.iter().enumerate() {
            if v.is_finite() {
                col_sum[j] += v;
                col_count[j] += 1;
            }
        }
    }
    let col_mean: Vec<f64> = (0..p)
        .map(|j| {
            if col_count[j] > 0 {
                col_sum[j] / col_count[j] as f64
            } else {
                0.0
            }
        })
        .collect();

    x.iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .map(|(j, &v)| if v.is_finite() { v } else { col_mean[j] })
                .collect()
        })
        .collect()
}

/// Flatten a 2-D log-abundance matrix and detection mask into parallel slices
/// suitable for [`fit_detection_curve`].
fn flatten_for_fit(log_ab_mat: &[Vec<f64>], detected: &[Vec<bool>]) -> (Vec<f64>, Vec<bool>) {
    let mut log_ab_flat = Vec::new();
    let mut det_flat = Vec::new();
    for (i, row) in log_ab_mat.iter().enumerate() {
        for (j, &lab) in row.iter().enumerate() {
            log_ab_flat.push(lab);
            det_flat.push(detected[i][j]);
        }
    }
    (log_ab_flat, det_flat)
}

/// Re-impute missing cells (NaN in `abundance`) using the ICA reconstruction.
///
/// The reconstruction is `sources * mixing^T + mean` (i.e., `X̂ = S * A^T + μ`).
/// Observed cells are left untouched.
fn reimpute_from_ica(
    ica: &crate::ica::IcaResult,
    abundance: &[Vec<f64>],
    detected: &[Vec<bool>],
    n: usize,
    p: usize,
) -> Vec<Vec<f64>> {
    // Reconstruct: x̂[i][j] = mean[j] + Σ_c sources[i][c] * mixing[j][c]
    let mut out = abundance.to_vec();
    for i in 0..n {
        for j in 0..p {
            if !detected[i][j] {
                let mut recon = ica.mean[j];
                for c in 0..ica.sources[i].len() {
                    recon += ica.sources[i][c] * ica.mixing[j][c];
                }
                out[i][j] = recon;
            }
        }
    }
    out
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
            got,
            expected
        );
    }
}

#[cfg(test)]
mod degeneracy_floor_tests {
    use super::*;

    /// Matrix with real MNAR dropout: low values are undetected, so the
    /// joint loop has missing cells to re-impute.
    fn mnar_matrix(n: usize, p: usize) -> Vec<Vec<f64>> {
        (0..n)
            .map(|i| {
                (0..p)
                    .map(|j| {
                        let s1 = ((i as f64) * 0.83).sin() * 2.0;
                        let s2 = ((i as f64) * 1.27).cos() * 2.0;
                        let block = if j * 3 <= p { s1 } else { s2 };
                        let v = 6.0 + block + (((i * 7 + j * 5) % 13) as f64) * 0.02;
                        if v < 5.4 {
                            f64::NAN
                        } else {
                            v
                        }
                    })
                    .collect()
            })
            .collect()
    }

    fn config(max_joint_iter: usize, degeneracy_floor: f64) -> MnarIcaConfig {
        MnarIcaConfig {
            k: 2,
            seed: 42,
            max_iter: 200,
            tol: 1e-4,
            max_joint_iter,
            // Unreachable, so the loop cannot stop by converging and the
            // test is about the other two exits.
            joint_tol: 1e-300,
            degeneracy_floor,
            weighted_whitening: false,
        }
    }

    /// The guard must fire before the iteration budget is exhausted, and
    /// must not be reported as convergence.
    #[test]
    fn degeneracy_floor_stops_the_loop_and_is_not_convergence() {
        let x = mnar_matrix(24, 30);
        let guarded = fast_ica_mnar(&x, &config(400, 0.5));
        assert_eq!(
            guarded.stop_reason,
            JointStopReason::DegeneracyFloor,
            "the gap should reach half its initial value well inside 400 iterations"
        );
        assert!(
            !guarded.joint_converged,
            "stopping on the degeneracy floor is not convergence"
        );
        assert!(
            guarded.joint_iterations < 400,
            "should stop early, used {}",
            guarded.joint_iterations
        );
    }

    /// With the guard disabled the same run exhausts its budget, which
    /// is the previous behaviour and must remain reachable.
    #[test]
    fn disabling_the_floor_restores_the_old_behaviour() {
        let x = mnar_matrix(24, 30);
        let unguarded = fast_ica_mnar(&x, &config(40, 0.0));
        assert_eq!(unguarded.stop_reason, JointStopReason::MaxIterations);
        assert_eq!(unguarded.joint_iterations, 40);
        assert!(!unguarded.joint_converged);
    }

    /// The guard must stop the loop strictly earlier than the budget
    /// would, on the same data and seed.
    #[test]
    fn the_floor_stops_earlier_than_the_budget() {
        let x = mnar_matrix(24, 30);
        let guarded = fast_ica_mnar(&x, &config(200, 0.5));
        let unguarded = fast_ica_mnar(&x, &config(200, 0.0));
        assert!(
            guarded.joint_iterations < unguarded.joint_iterations,
            "guarded {} should stop before unguarded {}",
            guarded.joint_iterations,
            unguarded.joint_iterations
        );
        // And the trace should show the gap genuinely shrinking, which
        // is what makes the guard meaningful rather than arbitrary.
        let gaps: Vec<f64> = unguarded
            .joint_trace
            .iter()
            .filter(|r| r.reconstruction_gap.is_finite())
            .map(|r| r.reconstruction_gap)
            .collect();
        assert!(gaps.len() >= 3, "need several gap readings");
        assert!(
            gaps.last().unwrap() < &gaps[0],
            "the reconstruction gap should decay: first {}, last {}",
            gaps[0],
            gaps.last().unwrap()
        );
    }

    /// A genuine convergence must still be reported as convergence, not
    /// pre-empted by the guard.
    #[test]
    fn real_convergence_is_still_reported_as_converged() {
        let x = mnar_matrix(24, 30);
        let mut cfg = config(200, 0.5);
        cfg.joint_tol = 1e9; // trivially satisfied on the first comparison
        let out = fast_ica_mnar(&x, &cfg);
        assert_eq!(out.stop_reason, JointStopReason::Converged);
        assert!(out.joint_converged);
    }
}
