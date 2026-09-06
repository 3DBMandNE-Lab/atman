use super::{k_selection_bound_hit, warn_if_k_selection_hit_bound};
use anyhow::{bail, Context, Result};
use atman_core::compositional::{apply_transform, Transform};
use atman_core::ica::{
    fast_ica, jaccard_top_n, pca_whiten, select_k_cumulative_variance, IcaResult,
};
use atman_core::ica_mnar::{fast_ica_mnar, MnarIcaConfig};
use clap::{Args as ClapArgs, ValueEnum};
use serde_json::json;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use super::program_name;
use crate::io::{
    atomic_write, format_float, hash_canonical_inputs, read_measurements_long, read_samples,
    sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct IcaArgs {
    /// Canonical Atman input directory (contains `measurements.tsv`).
    #[arg(long)]
    pub(super) input_dir: PathBuf,

    /// Explicit number of components. Overrides `--k-selection` if set.
    #[arg(long)]
    pub(super) k: Option<usize>,

    /// K-selection method; currently supports `cumulative-variance=<target>`.
    #[arg(long, default_value = "cumulative-variance=0.80")]
    pub(super) k_selection: String,

    /// Minimum allowed K (lower clamp for `--k-selection`).
    #[arg(long, default_value_t = 3)]
    pub(super) k_min: usize,

    /// Maximum allowed K (upper clamp for `--k-selection`).
    #[arg(long, default_value_t = 30)]
    pub(super) k_max: usize,

    /// Number of seeds (must be >= 1; first seed is the reference).
    #[arg(long, default_value_t = 50)]
    pub(super) n_seeds: usize,

    /// First seed value; subsequent seeds are sequential (seed, seed+1, ...).
    #[arg(long, default_value_t = 20260418)]
    pub(super) seed: u64,

    /// Jaccard threshold for a program to be considered "recovered" in another seed.
    #[arg(long, default_value_t = 0.9)]
    pub(super) seed_stability_threshold: f64,

    /// A program is flagged unstable when the fraction of alternative
    /// seeds that recover it (best-Jaccard >= `--seed-stability-threshold`)
    /// falls below this value. Default 0.9 — i.e. a program must recover
    /// in at least 90% of alternative seeds to pass.
    #[arg(long, default_value_t = 0.9)]
    pub(super) min_stable_seed_fraction: f64,

    /// Stability metric; currently `jaccard-top20` (i.e. top-N |loading| overlap).
    #[arg(long, value_enum, default_value_t = StabilityMetric::JaccardTop20)]
    pub(super) stability_metric: StabilityMetric,

    /// Top-N loadings used by the stability metric.
    #[arg(long, default_value_t = 20)]
    pub(super) stability_top_n: usize,

    /// Maximum FastICA iterations per seed.
    #[arg(long, default_value_t = 300)]
    pub(super) max_iter: usize,

    /// FastICA convergence tolerance.
    #[arg(long, default_value_t = 1e-4)]
    pub(super) tol: f64,

    /// Drop assays where the fraction of missing samples exceeds this value. 0 = strict complete-case.
    #[arg(long, default_value_t = 0.0)]
    pub(super) max_missing_fraction: f64,

    /// Imputation for residual missingness after the assay filter: `none` (fail) or `mean` (per-assay mean).
    #[arg(long, default_value = "none")]
    pub(super) impute: String,

    /// Loadings output TSV path (`program\tassay_id\tgene_symbol\tloading`).
    #[arg(long)]
    pub(super) output_loadings: PathBuf,

    /// Activations output TSV path (`cohort\tsubject_id\tsample_id\tprogram\tactivation`).
    #[arg(long)]
    pub(super) output_activations: PathBuf,

    /// Stability output TSV path.
    #[arg(long)]
    pub(super) output_stability: PathBuf,

    /// Cohort label written into the activations table (optional; defaults to directory name).
    #[arg(long)]
    pub(super) cohort: Option<String>,

    /// Compositional transform applied to the sample × protein matrix
    /// before ICA. Atman's canonical abundance is already on a log
    /// scale, so these are linear operations on log values. `clr`
    /// (per-sample mean centering) is the recommended default for
    /// closed-sum proteomics. `alr` / `ratio-anchor` require
    /// `--alr-reference` (a gene symbol). `ilr` emits loadings in
    /// the Helmert coordinate basis rather than raw protein space.
    #[arg(long, value_enum, default_value_t = TransformArg::None)]
    pub(super) transform: TransformArg,

    /// Missingness model. `none` (default) applies no MNAR weighting —
    /// all cells must be observed or already imputed. `abundance-conditional`
    /// fits a closed-form logistic detection curve and runs weighted FastICA
    /// where each cell's contribution is proportional to its estimated
    /// detection probability.
    #[arg(long, value_enum, default_value_t = MissingnessModel::None)]
    pub(super) missingness_model: MissingnessModel,

    /// Maximum outer (impute → fit curve → weighted ICA) iterations when
    /// `--missingness-model abundance-conditional` is active. 1 = single pass.
    #[arg(long, default_value_t = 5)]
    pub(super) max_joint_iter: usize,

    /// Convergence tolerance for the outer joint loop: the loop stops
    /// when `max(|Δbeta0|, |Δbeta1|)` of the detection curve falls
    /// below this. Only used with `--missingness-model
    /// abundance-conditional`. The per-iteration trace is written
    /// alongside the loadings as `mnar_joint_trace.tsv`, so a run that
    /// does not converge can be read as closing, oscillating, or
    /// diverging rather than only as a failure.
    #[arg(long, default_value_t = 1e-6)]
    pub(super) joint_tol: f64,

    /// Stop the joint loop when the imputed cells come within this
    /// fraction of their first-iteration distance to their own
    /// reconstruction. `0` disables the guard.
    ///
    /// The loop is iterated reconstruction-imputation and its fixed
    /// point is the state where imputed cells EQUAL their
    /// reconstruction, carrying no independent information while the
    /// fit explains them by construction. Running to convergence
    /// therefore produces a worse model than stopping earlier, so this
    /// is not a performance setting. Stopping here is reported as
    /// `degeneracy-floor` rather than as convergence.
    #[arg(long, default_value_t = 0.1)]
    pub(super) degeneracy_floor: f64,

    /// Gene symbol (matched against `samples` metadata `gene_symbol`)
    /// used as the reference for `--transform alr` or
    /// `--transform ratio-anchor`. Ignored for other transforms.
    #[arg(long)]
    pub(super) alr_reference: Option<String>,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum TransformArg {
    #[value(name = "none")]
    None,
    #[value(name = "clr")]
    Clr,
    #[value(name = "alr")]
    Alr,
    #[value(name = "ilr")]
    Ilr,
    #[value(name = "ratio-anchor")]
    RatioAnchor,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum StabilityMetric {
    #[value(name = "jaccard-top20")]
    JaccardTop20,
}

/// Missingness model for the ICA decomposition.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum MissingnessModel {
    /// No missingness modelling (default). Missing cells must have been
    /// already handled by `--impute` / `--max-missing-fraction`.
    #[value(name = "none")]
    None,
    /// Abundance-conditional MNAR model (Phase G). Fits a closed-form
    /// logistic detection curve and runs weighted FastICA where each cell
    /// contributes proportionally to its estimated detection probability.
    #[value(name = "abundance-conditional")]
    AbundanceConditional,
}

pub(super) fn run_ica(args: IcaArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if args.n_seeds == 0 {
        bail!("--n-seeds must be at least 1");
    }
    if args.k_min == 0 {
        bail!("--k-min must be at least 1");
    }
    if args.k_max < args.k_min {
        bail!("--k-max must be >= --k-min");
    }
    if !(0.0..=1.0).contains(&args.seed_stability_threshold) {
        bail!("--seed-stability-threshold must be in [0, 1]");
    }
    if !(0.0..=1.0).contains(&args.min_stable_seed_fraction) {
        bail!("--min-stable-seed-fraction must be in [0, 1]");
    }
    if args.stability_top_n == 0 {
        bail!("--stability-top-n must be at least 1");
    }

    // Load the matrix.  For `abundance-conditional`, we also load a NaN-
    // bearing copy for the MNAR ICA; the imputed matrix is still used for
    // k-selection and alternative-seed stability runs.
    let mnar_raw: Option<Vec<Vec<f64>>> =
        if args.missingness_model == MissingnessModel::AbundanceConditional {
            // load_matrix_mnar returns NaN for unobserved cells.
            let m = load_matrix_mnar(&args)?;
            Some(m.data)
        } else {
            None
        };
    // Standard load (imputes / bails as configured) used for k-selection,
    // stability alt-seeds, and as the plain-ICA path.
    let (matrix, transform_meta) = {
        let m = load_matrix(&args)?;
        apply_compositional(m, &args)?
    };

    let k = resolve_k(&matrix, &args)?;
    // Effective ceiling the rule actually searched under, for the
    // sidecar's bound-hit field below.
    let k_max_effective = args
        .k_max
        .min(matrix.samples.len().saturating_sub(1))
        .max(1);
    if k > matrix.samples.len() || k > matrix.assays.len() {
        bail!(
            "k={k} exceeds min(n_samples={}, n_assays={})",
            matrix.samples.len(),
            matrix.assays.len()
        );
    }

    eprintln!(
        "decompose ica: n_samples={} n_assays={} k={} n_seeds={} threshold={} transform={} missingness_model={:?}",
        matrix.samples.len(),
        matrix.assays.len(),
        k,
        args.n_seeds,
        args.seed_stability_threshold,
        transform_meta.name,
        args.missingness_model,
    );

    // Reference run: use MNAR-aware ICA if requested, else plain FastICA.
    let mut mnar_trace: Vec<atman_core::ica_mnar::JointIterationRecord> = Vec::new();
    let mut mnar_stop_reason: Option<&'static str> = None;
    let (ref_seed, ref_run, mnar_joint_iters) = if args.missingness_model
        == MissingnessModel::AbundanceConditional
    {
        let raw = mnar_raw.as_ref().expect("MNAR raw matrix");
        let mnar_config = MnarIcaConfig {
            k,
            seed: args.seed,
            max_iter: args.max_iter,
            tol: args.tol,
            max_joint_iter: args.max_joint_iter,
            joint_tol: args.joint_tol,
            degeneracy_floor: args.degeneracy_floor,
        };
        let mnar_result = fast_ica_mnar(raw, &mnar_config);
        eprintln!(
                "decompose ica: mnar joint_iterations={} joint_converged={} detection_curve beta0={:.4} beta1={:.4}",
                mnar_result.joint_iterations,
                mnar_result.joint_converged,
                mnar_result.detection_curve.beta0,
                mnar_result.detection_curve.beta1,
            );
        if mnar_result.stop_reason == atman_core::ica_mnar::JointStopReason::DegeneracyFloor {
            eprintln!(
                "decompose ica: the MNAR joint loop stopped at iteration {} on the degeneracy \
                 floor, not by converging: the imputed cells had come within {:.0}% of their \
                 first-iteration distance to their own reconstruction. Continuing would drive \
                 that to zero, at which point the missing cells carry no independent \
                 information and the fit explains them by construction. Raise \
                 --degeneracy-floor to stop earlier, or set it to 0 to disable the guard.",
                mnar_result.joint_iterations,
                100.0 * args.degeneracy_floor,
            );
        }
        if !mnar_result.joint_converged
            && mnar_result.stop_reason != atman_core::ica_mnar::JointStopReason::DegeneracyFloor
        {
            eprintln!(
                "decompose ica: warning: the MNAR joint loop did not converge — stopped at \
                 --max-joint-iter={} with final step {:.3e}, above --joint-tol={:.1e}. See \
                 mnar_joint_trace.tsv beside the loadings for the per-iteration betas.",
                args.max_joint_iter,
                mnar_result
                    .joint_trace
                    .last()
                    .map(|r| r.delta)
                    .unwrap_or(f64::NAN),
                args.joint_tol,
            );
        }
        mnar_trace = mnar_result.joint_trace.clone();
        mnar_stop_reason = Some(mnar_result.stop_reason.as_str());
        (
            args.seed,
            mnar_result.ica,
            Some(mnar_result.joint_iterations),
        )
    } else {
        let result = fast_ica(&matrix.data, k, args.seed, args.max_iter, args.tol);
        (args.seed, result, None)
    };

    // Alternative seeds run the SAME method as the reference.
    //
    // They used to always run plain FastICA, even when the reference was
    // the missingness-aware fit. The seed-stability column then compared
    // an MNAR reference against plain-FastICA alternatives at a Jaccard
    // threshold, which measures method disagreement rather than seed
    // stability, and in practice returned exactly zero for every program
    // — a column named `n_stable_runs` that could not report stability.
    let mut alt_runs: Vec<(u64, IcaResult)> = Vec::with_capacity(args.n_seeds.saturating_sub(1));
    for offset in 1..args.n_seeds {
        let seed = args.seed.wrapping_add(offset as u64);
        let result = if args.missingness_model == MissingnessModel::AbundanceConditional {
            let raw = mnar_raw.as_ref().expect("MNAR raw matrix");
            fast_ica_mnar(
                raw,
                &MnarIcaConfig {
                    k,
                    seed,
                    max_iter: args.max_iter,
                    tol: args.tol,
                    max_joint_iter: args.max_joint_iter,
                    joint_tol: args.joint_tol,
                    degeneracy_floor: args.degeneracy_floor,
                },
            )
            .ica
        } else {
            fast_ica(&matrix.data, k, seed, args.max_iter, args.tol)
        };
        alt_runs.push((seed, result));
    }

    // Convergence reporting. `IcaResult` has always carried
    // `n_iterations` and `final_tol`; nothing read them, so a fit that
    // stopped at the iteration cap was indistinguishable from a
    // converged one in both stderr and the sidecar. Same failure class
    // as a selection rule silently hitting its bound.
    let ref_converged = ref_run.final_tol < args.tol;
    let alt_not_converged = alt_runs
        .iter()
        .filter(|(_, r)| r.final_tol >= args.tol)
        .count();
    let seeds_not_converged = alt_not_converged + usize::from(!ref_converged);
    if !ref_converged {
        eprintln!(
            "decompose ica: warning: the reference fit did not converge — stopped at \
             --max-iter={} with final tolerance {:.3e}, above the requested --tol={:.3e}. The \
             loadings are the state of the solver when it ran out of iterations, not a \
             converged decomposition. Raise --max-iter or relax --tol.",
            args.max_iter, ref_run.final_tol, args.tol,
        );
    }
    if alt_not_converged > 0 {
        eprintln!(
            "decompose ica: warning: {alt_not_converged} of {} alternative seeds did not \
             converge; the seed-stability ranking is computed over unconverged fits.",
            alt_runs.len(),
        );
    }

    let reference = canonicalize(&ref_run);
    let cohort = args
        .cohort
        .clone()
        .unwrap_or_else(|| derive_cohort(&args.input_dir));

    write_loadings(&args.output_loadings, &matrix, &reference)?;
    write_activations(&args.output_activations, &matrix, &reference, &cohort)?;
    write_stability(
        &args.output_stability,
        &reference,
        &alt_runs,
        args.stability_top_n,
        args.seed_stability_threshold,
        args.min_stable_seed_fraction,
        ref_seed,
    )?;

    eprintln!(
        "decompose ica: reference_seed={} loadings={} activations={} stability={}",
        ref_seed,
        args.output_loadings.display(),
        args.output_activations.display(),
        args.output_stability.display()
    );

    // Emit `transform_applied.json` next to the loadings whenever a
    // non-`none` transform ran, so downstream tools can audit which
    // coordinate system the archetypes live in.
    // Per-iteration joint-loop trace, so a non-converged MNAR fit can be
    // characterised rather than only reported as a failure.
    let mut mnar_trace_path: Option<PathBuf> = None;
    if !mnar_trace.is_empty() {
        let trace_path = args
            .output_loadings
            .parent()
            .unwrap_or(Path::new("."))
            .join("mnar_joint_trace.tsv");
        let mut buf = String::from(
            "iteration\tbeta0\tbeta1\tdelta\timputation_delta\treconstruction_gap\tconverged_here\n",
        );
        for (i, r) in mnar_trace.iter().enumerate() {
            let is_last = i + 1 == mnar_trace.len();
            buf.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                r.iteration,
                format_float(r.beta0),
                format_float(r.beta1),
                if r.delta.is_finite() {
                    format_float(r.delta)
                } else {
                    "NA".to_string()
                },
                if r.imputation_delta.is_finite() {
                    format_float(r.imputation_delta)
                } else {
                    "NA".to_string()
                },
                if r.reconstruction_gap.is_finite() {
                    format_float(r.reconstruction_gap)
                } else {
                    "NA".to_string()
                },
                u8::from(is_last && r.delta.is_finite() && r.delta < args.joint_tol),
            ));
        }
        atomic_write(&trace_path, buf.as_bytes())?;
        eprintln!("decompose ica: mnar joint trace={}", trace_path.display());
        mnar_trace_path = Some(trace_path);
    }

    let transform_applied_path = if transform_meta.name != "none" {
        let p = args
            .output_loadings
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| PathBuf::from("."))
            .join("transform_applied.json");
        write_transform_applied(&p, &transform_meta)?;
        eprintln!("decompose ica: transform_applied={}", p.display());
        Some(p)
    } else {
        None
    };

    let finished_at = SystemTime::now();
    let canonical_inputs: &[&str] = &["measurements.tsv", "samples.tsv", "proteins.tsv"];
    let inputs_sha256 = hash_canonical_inputs(&args.input_dir, canonical_inputs)?;
    let sidecar = sidecar_path_for(&args.output_loadings);
    let stability_metric = match args.stability_metric {
        StabilityMetric::JaccardTop20 => "jaccard-top20",
    };
    let missingness_model_str = match args.missingness_model {
        MissingnessModel::None => "none",
        MissingnessModel::AbundanceConditional => "abundance-conditional",
    };
    let mut outputs = vec![
        args.output_loadings.clone(),
        args.output_activations.clone(),
        args.output_stability.clone(),
    ];
    if let Some(p) = &transform_applied_path {
        outputs.push(p.clone());
    }
    if let Some(p) = &mnar_trace_path {
        outputs.push(p.clone());
    }
    write_run_sidecar(
        &sidecar,
        "decompose ica",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "k": args.k,
            "k-selection": args.k_selection,
            "k-min": args.k_min,
            "k-max": args.k_max,
            "n-seeds": args.n_seeds,
            "seed": args.seed,
            "seed-stability-threshold": args.seed_stability_threshold,
            "min-stable-seed-fraction": args.min_stable_seed_fraction,
            "stability-metric": stability_metric,
            "stability-top-n": args.stability_top_n,
            "max-iter": args.max_iter,
            "tol": args.tol,
            "max-missing-fraction": args.max_missing_fraction,
            "impute": args.impute,
            "missingness-model": missingness_model_str,
            "max-joint-iter": args.max_joint_iter,
            "joint-tol": args.joint_tol,
            "mnar-joint-iterations": mnar_joint_iters,
            "output-loadings": args.output_loadings.display().to_string(),
            "output-activations": args.output_activations.display().to_string(),
            "output-stability": args.output_stability.display().to_string(),
            "cohort": args.cohort,
            "transform": transform_meta.name,
            "alr-reference": transform_meta.alr_reference,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        {
            // Resolved-value fields: record the integer `k` actually used
            // so a reviewer reading the sidecar doesn't have to count
            // columns in loadings.tsv to find out which `k` the
            // `--k-selection` rule picked. `k_resolution_source` names
            // the path: `cli` when --k was passed, `<rule>` otherwise.
            let mut extras = serde_json::Map::new();
            extras.insert("k_resolved".into(), serde_json::json!(k));
            let source = if args.k.is_some() {
                "cli".to_string()
            } else {
                args.k_selection.clone()
            };
            extras.insert("k_resolution_source".into(), serde_json::json!(source));
            // Whether the solver converged, and the evidence for it.
            // `n_iterations` equal to `--max-iter` with `final_tol`
            // above `--tol` is a fit that ran out of iterations.
            extras.insert("ica_converged".into(), serde_json::json!(ref_converged));
            if let Some(reason) = mnar_stop_reason {
                extras.insert("mnar_joint_stop_reason".into(), serde_json::json!(reason));
            }
            // Alternative seeds now run the same method as the
            // reference, so the seed-stability column measures seed
            // stability under both missingness models.
            extras.insert(
                "stability_alt_seed_method".into(),
                serde_json::json!(if args.missingness_model
                    == MissingnessModel::AbundanceConditional
                {
                    "abundance-conditional"
                } else {
                    "fastica"
                }),
            );
            extras.insert(
                "ica_n_iterations".into(),
                serde_json::json!(ref_run.n_iterations),
            );
            extras.insert("ica_final_tol".into(), serde_json::json!(ref_run.final_tol));
            extras.insert(
                "ica_n_seeds_not_converged".into(),
                serde_json::json!(seeds_not_converged),
            );
            // Whether the rule ran into its own search bound. `k_max`
            // means the criterion was never satisfied inside the sweep
            // and `k_resolved` is a truncation, not a selection.
            extras.insert(
                "k_selection_bound_hit".into(),
                serde_json::json!(k_selection_bound_hit(
                    args.k.is_some(),
                    k,
                    args.k_min,
                    k_max_effective,
                )),
            );
            Some(extras)
        },
    )?;
    eprintln!("decompose ica: sidecar={}", sidecar.display());
    Ok(())
}

fn resolve_k(matrix: &AbundanceMatrix, args: &IcaArgs) -> Result<usize> {
    if let Some(k) = args.k {
        if k == 0 {
            bail!("--k must be at least 1");
        }
        return Ok(k);
    }
    let target = parse_k_selection(&args.k_selection)?;
    let whitening = pca_whiten(&matrix.data, matrix.samples.len().min(matrix.assays.len()));
    let effective_k_max = args
        .k_max
        .min(matrix.samples.len().saturating_sub(1))
        .max(1);
    let chosen = select_k_cumulative_variance(
        &whitening.full_spectrum,
        target,
        args.k_min,
        effective_k_max,
    );
    warn_if_k_selection_hit_bound(
        "decompose ica",
        &args.k_selection,
        chosen,
        args.k_min,
        effective_k_max,
    );
    Ok(chosen)
}

fn parse_k_selection(spec: &str) -> Result<f64> {
    let trimmed = spec.trim();
    if let Some(rest) = trimmed.strip_prefix("cumulative-variance=") {
        let value: f64 = rest
            .parse()
            .with_context(|| format!("parsing k-selection target in {:?}", spec))?;
        if !(0.0..=1.0).contains(&value) {
            bail!("cumulative-variance target must be in (0, 1]");
        }
        Ok(value)
    } else {
        bail!(
            "unsupported --k-selection {:?}; use cumulative-variance=<target>",
            spec
        )
    }
}

fn derive_cohort(path: &Path) -> String {
    path.file_name()
        .and_then(|os| os.to_str())
        .unwrap_or("cohort")
        .to_string()
}

#[derive(Debug)]
pub(super) struct AbundanceMatrix {
    pub(super) samples: Vec<SampleInfo>,
    pub(super) assays: Vec<AssayInfo>,
    pub(super) data: Vec<Vec<f64>>,
}

#[derive(Debug, Clone)]
pub(super) struct SampleInfo {
    sample_id: String,
    subject_id: String,
}

#[derive(Debug, Clone)]
pub(super) struct AssayInfo {
    assay_id: String,
    gene_symbol: String,
}

pub(super) fn load_matrix(args: &IcaArgs) -> Result<AbundanceMatrix> {
    let tsv = args.input_dir.join("measurements.tsv");
    let records = read_measurements_long(&tsv)?;
    if records.is_empty() {
        bail!("no measurements in {:?}", tsv);
    }

    let mut abundance_by_key: BTreeMap<(String, String), (f64, Option<String>)> = BTreeMap::new();
    let mut sample_order: Vec<String> = Vec::new();
    let mut seen_samples: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    let mut assay_meta: BTreeMap<String, Option<String>> = BTreeMap::new();
    for r in &records {
        if r.dropped_by_qc {
            continue;
        }
        let assay = r.assay_id.0.clone();
        let sample = r.sample_id.clone();
        let abundance = r.abundance.as_f64();
        if !abundance.is_finite() {
            continue;
        }
        if seen_samples.insert(sample.clone()) {
            sample_order.push(sample.clone());
        }
        assay_meta
            .entry(assay.clone())
            .or_insert_with(|| r.gene_symbol.clone());
        abundance_by_key.insert((assay, sample), (abundance, r.gene_symbol.clone()));
    }

    // Load samples.tsv for subject IDs (if present).
    let samples_path = args.input_dir.join("samples.tsv");
    let subject_lookup: BTreeMap<String, String> = if samples_path.exists() {
        crate::io::read_samples(&samples_path)?
            .into_iter()
            .map(|s| (s.sample_id.clone(), s.subject_id.unwrap_or_default()))
            .collect()
    } else {
        BTreeMap::new()
    };

    if !(0.0..=1.0).contains(&args.max_missing_fraction) {
        bail!("--max-missing-fraction must be in [0, 1]");
    }
    let impute_mean = match args.impute.as_str() {
        "none" => false,
        "mean" => true,
        other => bail!("--impute {:?}; expected none or mean", other),
    };

    // Drop assays whose missing-sample fraction exceeds the configured threshold.
    let n_samples = sample_order.len();
    if n_samples < 4 {
        bail!("need at least 4 samples for ICA, got {}", n_samples);
    }
    let mut kept_assays: Vec<String> = Vec::new();
    let mut dropped_assays = 0usize;
    for assay in assay_meta.keys() {
        let present = sample_order
            .iter()
            .filter(|s| abundance_by_key.contains_key(&(assay.clone(), (*s).clone())))
            .count();
        let missing_fraction = 1.0 - present as f64 / n_samples as f64;
        if missing_fraction <= args.max_missing_fraction + 1e-12 {
            kept_assays.push(assay.clone());
        } else {
            dropped_assays += 1;
        }
    }
    if kept_assays.is_empty() {
        bail!(
            "no assays retained with --max-missing-fraction={} over {} samples in {:?}",
            args.max_missing_fraction,
            n_samples,
            tsv
        );
    }
    eprintln!(
        "decompose ica: retained {} assays, dropped {} for exceeding --max-missing-fraction={}",
        kept_assays.len(),
        dropped_assays,
        args.max_missing_fraction
    );

    let assays: Vec<AssayInfo> = kept_assays
        .iter()
        .map(|a| AssayInfo {
            assay_id: a.clone(),
            gene_symbol: assay_meta.get(a).cloned().flatten().unwrap_or_default(),
        })
        .collect();
    let samples: Vec<SampleInfo> = sample_order
        .iter()
        .map(|s| SampleInfo {
            sample_id: s.clone(),
            subject_id: subject_lookup.get(s).cloned().unwrap_or_default(),
        })
        .collect();

    let mut data = vec![vec![0.0_f64; assays.len()]; samples.len()];
    let mut missing_cells = 0usize;
    for (j, assay) in assays.iter().enumerate() {
        let mut values: Vec<Option<f64>> = Vec::with_capacity(samples.len());
        let mut sum = 0.0_f64;
        let mut count = 0usize;
        for sample in &samples {
            match abundance_by_key.get(&(assay.assay_id.clone(), sample.sample_id.clone())) {
                Some((v, _)) => {
                    values.push(Some(*v));
                    sum += *v;
                    count += 1;
                }
                None => values.push(None),
            }
        }
        let mean = if count > 0 { sum / count as f64 } else { 0.0 };
        for (i, v) in values.iter().enumerate() {
            match v {
                Some(x) => data[i][j] = *x,
                None => {
                    if impute_mean {
                        data[i][j] = mean;
                        missing_cells += 1;
                    } else {
                        bail!(
                            "assay {} is missing sample {} and --impute=none; rerun with --impute mean or lower --max-missing-fraction",
                            assay.assay_id,
                            samples[i].sample_id
                        );
                    }
                }
            }
        }
    }
    if missing_cells > 0 {
        eprintln!(
            "decompose ica: mean-imputed {} missing cells across {} assays",
            missing_cells,
            assays.len()
        );
    }
    Ok(AbundanceMatrix {
        samples,
        assays,
        data,
    })
}

/// Load the abundance matrix **preserving NaN for missing cells**.
///
/// Like [`load_matrix`] but never imputes: missing cells are left as
/// `f64::NAN` so that [`fast_ica_mnar`] can model the missingness pattern.
/// The `--max-missing-fraction` assay filter still applies; the `--impute`
/// flag is ignored (the MNAR model handles missingness internally).
fn load_matrix_mnar(args: &IcaArgs) -> Result<AbundanceMatrix> {
    let tsv = args.input_dir.join("measurements.tsv");
    let records = read_measurements_long(&tsv)?;
    if records.is_empty() {
        bail!("no measurements in {:?}", tsv);
    }

    let mut abundance_by_key: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut sample_order: Vec<String> = Vec::new();
    let mut seen_samples: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    let mut assay_meta: BTreeMap<String, Option<String>> = BTreeMap::new();

    for r in &records {
        if r.dropped_by_qc {
            continue;
        }
        let assay = r.assay_id.0.clone();
        let sample = r.sample_id.clone();
        let abundance = r.abundance.as_f64();
        if !abundance.is_finite() {
            // Non-finite in source data (e.g. NaN sentinel): mark as missing
            // but still register the sample and assay.
            if seen_samples.insert(sample.clone()) {
                sample_order.push(sample.clone());
            }
            assay_meta
                .entry(assay.clone())
                .or_insert_with(|| r.gene_symbol.clone());
            continue;
        }
        if seen_samples.insert(sample.clone()) {
            sample_order.push(sample.clone());
        }
        assay_meta
            .entry(assay.clone())
            .or_insert_with(|| r.gene_symbol.clone());
        abundance_by_key.insert((assay, sample), abundance);
    }

    // Also register samples/assays where the value is NaN (missing sentinel).
    // We need to know all sample_ids even for NaN rows so the matrix dimensions
    // are consistent. Re-scan for NaN rows to capture any missing sample/assay.
    for r in &records {
        if r.dropped_by_qc {
            continue;
        }
        let assay = r.assay_id.0.clone();
        let sample = r.sample_id.clone();
        if seen_samples.insert(sample.clone()) {
            sample_order.push(sample.clone());
        }
        assay_meta
            .entry(assay)
            .or_insert_with(|| r.gene_symbol.clone());
    }

    let samples_path = args.input_dir.join("samples.tsv");
    let subject_lookup: BTreeMap<String, String> = if samples_path.exists() {
        read_samples(&samples_path)?
            .into_iter()
            .map(|s| (s.sample_id.clone(), s.subject_id.unwrap_or_default()))
            .collect()
    } else {
        BTreeMap::new()
    };

    if !(0.0..=1.0).contains(&args.max_missing_fraction) {
        bail!("--max-missing-fraction must be in [0, 1]");
    }

    let n_samples = sample_order.len();
    if n_samples < 4 {
        bail!("need at least 4 samples for ICA, got {}", n_samples);
    }
    let mut kept_assays: Vec<String> = Vec::new();
    let mut dropped_assays = 0usize;
    for assay in assay_meta.keys() {
        let present = sample_order
            .iter()
            .filter(|s| abundance_by_key.contains_key(&(assay.clone(), (*s).clone())))
            .count();
        let missing_fraction = 1.0 - present as f64 / n_samples as f64;
        if missing_fraction <= args.max_missing_fraction + 1e-12 {
            kept_assays.push(assay.clone());
        } else {
            dropped_assays += 1;
        }
    }
    if kept_assays.is_empty() {
        bail!(
            "no assays retained with --max-missing-fraction={} over {} samples in {:?}",
            args.max_missing_fraction,
            n_samples,
            tsv
        );
    }
    eprintln!(
        "decompose ica (mnar): retained {} assays, dropped {} for exceeding --max-missing-fraction={}",
        kept_assays.len(),
        dropped_assays,
        args.max_missing_fraction
    );

    let assays: Vec<AssayInfo> = kept_assays
        .iter()
        .map(|a| AssayInfo {
            assay_id: a.clone(),
            gene_symbol: assay_meta.get(a).cloned().flatten().unwrap_or_default(),
        })
        .collect();
    let samples: Vec<SampleInfo> = sample_order
        .iter()
        .map(|s| SampleInfo {
            sample_id: s.clone(),
            subject_id: subject_lookup.get(s).cloned().unwrap_or_default(),
        })
        .collect();

    // Fill data: NaN for missing cells.
    let mut data = vec![vec![f64::NAN; assays.len()]; samples.len()];
    let mut missing_cells = 0usize;
    for (j, assay) in assays.iter().enumerate() {
        for (i, sample) in samples.iter().enumerate() {
            match abundance_by_key.get(&(assay.assay_id.clone(), sample.sample_id.clone())) {
                Some(&v) => data[i][j] = v,
                None => {
                    missing_cells += 1;
                    // data[i][j] is already NaN from initialisation.
                }
            }
        }
    }
    eprintln!(
        "decompose ica (mnar): {} NaN cells ({} total)",
        missing_cells,
        samples.len() * assays.len()
    );
    Ok(AbundanceMatrix {
        samples,
        assays,
        data,
    })
}

/// Canonicalize an ICA result: sign-flip each program so its max |loading| is positive,
/// then sort programs by descending max |loading|. This makes output deterministic
/// under sign ambiguity and yields stable program ordering for reporting.
/// Metadata written into `transform_applied.json` alongside
/// `loadings.tsv` whenever a non-`none` transform runs. Lets reviewers
/// audit which coordinate system the archetypes live in.
struct TransformMeta {
    name: &'static str,
    alr_reference: Option<String>,
    n_input_proteins: usize,
    n_output_coords: usize,
}

fn apply_compositional(
    matrix: AbundanceMatrix,
    args: &IcaArgs,
) -> Result<(AbundanceMatrix, TransformMeta)> {
    let transform_name = match args.transform {
        TransformArg::None => "none",
        TransformArg::Clr => "clr",
        TransformArg::Alr => "alr",
        TransformArg::Ilr => "ilr",
        TransformArg::RatioAnchor => "ratio-anchor",
    };
    if args.transform == TransformArg::None {
        let n = matrix.assays.len();
        return Ok((
            matrix,
            TransformMeta {
                name: transform_name,
                alr_reference: None,
                n_input_proteins: n,
                n_output_coords: n,
            },
        ));
    }
    if matches!(args.transform, TransformArg::Ilr) && args.alr_reference.is_some() {
        eprintln!("decompose ica: --alr-reference is ignored for --transform ilr");
    }
    let needs_reference = matches!(
        args.transform,
        TransformArg::Alr | TransformArg::RatioAnchor
    );
    let resolved_reference: Option<(usize, String)> = if needs_reference {
        let gene = match args.alr_reference.as_deref() {
            Some(g) if !g.is_empty() => g,
            _ => bail!("--transform {transform_name} requires --alr-reference <gene_symbol>"),
        };
        let idx = matrix
            .assays
            .iter()
            .position(|a| a.gene_symbol == gene)
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "--alr-reference {gene:?} not found among {} assays",
                    matrix.assays.len()
                )
            })?;
        Some((idx, gene.to_string()))
    } else {
        None
    };
    let core_transform = match args.transform {
        TransformArg::None => Transform::None,
        TransformArg::Clr => Transform::Clr,
        TransformArg::Ilr => Transform::Ilr,
        TransformArg::Alr => Transform::Alr {
            reference_index: resolved_reference.as_ref().unwrap().0,
        },
        TransformArg::RatioAnchor => Transform::RatioAnchor {
            reference_index: resolved_reference.as_ref().unwrap().0,
        },
    };
    let n_input = matrix.assays.len();
    let transformed = apply_transform(&matrix.data, core_transform)
        .map_err(|e| anyhow::anyhow!("transform {transform_name} failed: {e}"))?;
    let n_output = if transformed.is_empty() {
        0
    } else {
        transformed[0].len()
    };
    let assays = if matches!(args.transform, TransformArg::Ilr) {
        // ILR coordinates are in the Helmert basis, not raw
        // protein space. Emit synthetic assay IDs so downstream
        // writers don't try to match them back to a protein.
        (0..n_output)
            .map(|i| AssayInfo {
                assay_id: format!("ilr_coord_{:03}", i + 1),
                gene_symbol: format!("ILR{:03}", i + 1),
            })
            .collect()
    } else {
        matrix.assays
    };
    Ok((
        AbundanceMatrix {
            samples: matrix.samples,
            assays,
            data: transformed,
        },
        TransformMeta {
            name: transform_name,
            alr_reference: resolved_reference.map(|(_, g)| g),
            n_input_proteins: n_input,
            n_output_coords: n_output,
        },
    ))
}

fn write_transform_applied(path: &Path, meta: &TransformMeta) -> Result<()> {
    let mut obj = serde_json::Map::new();
    obj.insert("transform".into(), json!(meta.name));
    obj.insert(
        "alr_reference".into(),
        match &meta.alr_reference {
            Some(g) => json!(g),
            None => serde_json::Value::Null,
        },
    );
    obj.insert("n_input_proteins".into(), json!(meta.n_input_proteins));
    obj.insert("n_output_coords".into(), json!(meta.n_output_coords));
    atomic_write(
        path,
        serde_json::to_string_pretty(&serde_json::Value::Object(obj))?.as_bytes(),
    )?;
    Ok(())
}

fn canonicalize(result: &IcaResult) -> CanonicalIca {
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
        mb.partial_cmp(&ma).unwrap_or(Ordering::Equal)
    });
    let ordered_loadings: Vec<Vec<f64>> = order.iter().map(|&i| loadings[i].clone()).collect();
    let ordered_activations: Vec<Vec<f64>> =
        order.iter().map(|&i| activations[i].clone()).collect();
    CanonicalIca {
        loadings: ordered_loadings,
        activations: ordered_activations,
    }
}

struct CanonicalIca {
    /// `k` rows, each of length `p`.
    loadings: Vec<Vec<f64>>,
    /// `k` rows, each of length `n`.
    activations: Vec<Vec<f64>>,
}

fn write_loadings(path: &Path, matrix: &AbundanceMatrix, run: &CanonicalIca) -> Result<()> {
    let mut out = String::from("program\tassay_id\tgene_symbol\tloading\n");
    for (c, row) in run.loadings.iter().enumerate() {
        let program = program_name(c);
        for (j, load) in row.iter().enumerate() {
            let assay = &matrix.assays[j];
            out.push_str(&program);
            out.push('\t');
            out.push_str(&assay.assay_id);
            out.push('\t');
            out.push_str(&assay.gene_symbol);
            out.push('\t');
            out.push_str(&format_float(*load));
            out.push('\n');
        }
    }
    atomic_write(path, out.as_bytes())
}

fn write_activations(
    path: &Path,
    matrix: &AbundanceMatrix,
    run: &CanonicalIca,
    cohort: &str,
) -> Result<()> {
    let mut out = String::from("cohort\tsubject_id\tsample_id\tprogram\tactivation\n");
    for (c, row) in run.activations.iter().enumerate() {
        let program = program_name(c);
        for (i, v) in row.iter().enumerate() {
            let sample = &matrix.samples[i];
            out.push_str(cohort);
            out.push('\t');
            out.push_str(&sample.subject_id);
            out.push('\t');
            out.push_str(&sample.sample_id);
            out.push('\t');
            out.push_str(&program);
            out.push('\t');
            out.push_str(&format_float(*v));
            out.push('\n');
        }
    }
    atomic_write(path, out.as_bytes())
}

#[allow(clippy::too_many_arguments)]
fn write_stability(
    path: &Path,
    reference: &CanonicalIca,
    others: &[(u64, IcaResult)],
    top_n: usize,
    jaccard_threshold: f64,
    min_stable_seed_fraction: f64,
    ref_seed: u64,
) -> Result<()> {
    let mut out = String::from(
        "program\treference_seed\tn_alt_seeds\tn_stable_runs\tseed_stability_fraction\tmean_best_jaccard\tflag_below_threshold\n",
    );
    let canon_alts: Vec<CanonicalIca> = others.iter().map(|(_, r)| canonicalize(r)).collect();
    for (c, row) in reference.loadings.iter().enumerate() {
        let program = program_name(c);
        let mut n_stable = 0usize;
        let mut jaccards: Vec<f64> = Vec::with_capacity(canon_alts.len());
        for alt in &canon_alts {
            let best = alt
                .loadings
                .iter()
                .map(|alt_row| jaccard_top_n(row, alt_row, top_n))
                .fold(0.0_f64, f64::max);
            jaccards.push(best);
            if best >= jaccard_threshold {
                n_stable += 1;
            }
        }
        let n_alt = canon_alts.len();
        let fraction = if n_alt == 0 {
            1.0
        } else {
            n_stable as f64 / n_alt as f64
        };
        let mean = if jaccards.is_empty() {
            f64::NAN
        } else {
            jaccards.iter().sum::<f64>() / jaccards.len() as f64
        };
        let flag = if n_alt > 0 && fraction < min_stable_seed_fraction {
            "1"
        } else {
            "0"
        };
        out.push_str(&format!(
            "{program}\t{ref_seed}\t{n_alt}\t{n_stable}\t{}\t{}\t{flag}\n",
            format_float(fraction),
            format_float(mean),
        ));
    }
    atomic_write(path, out.as_bytes())
}
