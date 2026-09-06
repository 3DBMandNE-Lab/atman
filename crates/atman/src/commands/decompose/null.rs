use anyhow::{bail, Context, Result};
use atman_core::bh_fdr;
use atman_core::ica_null::{archetype_null, ArchetypeNullRow, NullMode, NullParams};
use clap::{Args as ClapArgs, ValueEnum};
use serde_json::json;
use std::path::PathBuf;
use std::time::SystemTime;

use super::ica::{load_matrix, IcaArgs, MissingnessModel, StabilityMetric, TransformArg};
use super::program_name;
use crate::io::{
    atomic_write, format_float, hash_canonical_inputs, sidecar_path_for, write_run_sidecar,
};

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum NullModeArg {
    #[value(name = "sample-shuffle")]
    SampleShuffle,
    #[value(name = "protein-shuffle")]
    ProteinShuffle,
    #[value(name = "gaussian-matched")]
    GaussianMatched,
}

impl NullModeArg {
    fn to_core(self) -> NullMode {
        match self {
            Self::SampleShuffle => NullMode::SampleShuffle,
            Self::ProteinShuffle => NullMode::ProteinShuffle,
            Self::GaussianMatched => NullMode::GaussianMatched,
        }
    }
}

#[derive(ClapArgs, Debug)]
pub struct NullArgs {
    /// Canonical Atman input directory (contains `measurements.tsv`).
    #[arg(long)]
    input_dir: PathBuf,

    /// Number of ICA components to fit. Required — no k-selection here
    /// so the null run exactly mirrors the K chosen in `decompose ica`.
    #[arg(long)]
    k: usize,

    /// Number of null permutation iterations. Higher = tighter p-value
    /// resolution. Floor on p is `1/(n_perm + 1)` so use ≥ 200 for
    /// `q < 0.05` claims to mean anything.
    #[arg(long, default_value_t = 200)]
    n_perm: usize,

    /// FastICA seeds per matrix (real and each null). Seed 0 is the
    /// reference; stability is mean best-Jaccard across the remaining
    /// `n_seeds − 1` alternative seeds. Must be ≥ 2.
    #[arg(long, default_value_t = 3)]
    n_seeds: usize,

    /// Top-level seed. Each null iteration's sub-seed derives from
    /// `(seed, iter)` deterministically.
    #[arg(long, default_value_t = 20260418)]
    seed: u64,

    /// Null-matrix construction mode.
    #[arg(long, value_enum, default_value_t = NullModeArg::ProteinShuffle)]
    null_mode: NullModeArg,

    /// Top-N loadings used by the Jaccard stability metric.
    #[arg(long, default_value_t = 20)]
    top_n: usize,

    /// FastICA max iterations per seed.
    #[arg(long, default_value_t = 300)]
    max_iter: usize,

    /// FastICA convergence tolerance.
    #[arg(long, default_value_t = 1e-4)]
    tol: f64,

    /// Drop assays with more than this fraction of missing samples.
    #[arg(long, default_value_t = 0.0)]
    max_missing_fraction: f64,

    /// Residual missingness imputation: `none` fails, `mean` uses
    /// per-assay mean.
    #[arg(long, default_value = "none")]
    impute: String,

    /// BH-q threshold used by the output `decision` column
    /// (`signal` when `null_q < threshold`, else `noise`).
    #[arg(long, default_value_t = 0.05)]
    q_threshold: f64,

    /// Output TSV: one row per program with observed_stability,
    /// null distribution summary, null_p, null_q, decision.
    #[arg(long)]
    output: PathBuf,
}

pub(super) fn run_null(args: NullArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if args.k == 0 {
        bail!("--k must be at least 1");
    }
    if args.n_perm == 0 {
        bail!("--n-perm must be at least 1");
    }
    if args.n_seeds < 2 {
        bail!("--n-seeds must be >= 2 for stability calibration");
    }
    if args.top_n == 0 {
        bail!("--top-n must be at least 1");
    }
    if !(0.0..=1.0).contains(&args.q_threshold) {
        bail!("--q-threshold must be in [0, 1]");
    }

    // Reuse the abundance-matrix loader from `run_ica` by reconstructing
    // an equivalent IcaArgs; `load_matrix` only reads the fields below.
    let adapter = IcaArgs {
        input_dir: args.input_dir.clone(),
        k: Some(args.k),
        k_selection: "cumulative-variance=0.80".to_string(),
        k_min: 1,
        k_max: args.k.max(1),
        n_seeds: args.n_seeds,
        seed: args.seed,
        seed_stability_threshold: 0.9,
        min_stable_seed_fraction: 0.9,
        stability_metric: StabilityMetric::JaccardTop20,
        stability_top_n: args.top_n,
        max_iter: args.max_iter,
        tol: args.tol,
        max_missing_fraction: args.max_missing_fraction,
        impute: args.impute.clone(),
        output_loadings: PathBuf::new(),
        output_activations: PathBuf::new(),
        output_stability: PathBuf::new(),
        cohort: None,
        transform: TransformArg::None,
        alr_reference: None,
        missingness_model: MissingnessModel::None,
        max_joint_iter: 5,
        joint_tol: 1e-6,
        degeneracy_floor: 0.1,
        weighted_whitening: false,
    };
    let matrix = load_matrix(&adapter)?;
    if args.k > matrix.samples.len() || args.k > matrix.assays.len() {
        bail!(
            "k={} exceeds min(n_samples={}, n_assays={})",
            args.k,
            matrix.samples.len(),
            matrix.assays.len()
        );
    }

    eprintln!(
        "decompose null: n_samples={} n_assays={} k={} n_perm={} n_seeds={} mode={}",
        matrix.samples.len(),
        matrix.assays.len(),
        args.k,
        args.n_perm,
        args.n_seeds,
        args.null_mode.to_core().as_str()
    );

    let params = NullParams {
        k: args.k,
        n_perm: args.n_perm,
        n_seeds: args.n_seeds,
        seed: args.seed,
        mode: args.null_mode.to_core(),
        top_n: args.top_n,
        max_iter: args.max_iter,
        tol: args.tol,
    };
    let rows: Vec<ArchetypeNullRow> =
        archetype_null(&matrix.data, &params).map_err(|e| anyhow::anyhow!(e))?;

    // BH on null_p across programs in this run.
    let ps: Vec<Option<f64>> = rows.iter().map(|r| Some(r.null_p)).collect();
    let qs = bh_fdr(&ps);

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output dir {:?}", parent))?;
    }
    let mut out = String::from(
        "program\tobserved_stability\tnull_stability_mean\tnull_stability_p95\tnull_p\tnull_q\tdecision\n",
    );
    for (row, q) in rows.iter().zip(qs.iter()) {
        let program = program_name(row.program - 1);
        let q_val = q.unwrap_or(f64::NAN);
        let decision = if q_val.is_finite() && q_val < args.q_threshold {
            "signal"
        } else {
            "noise"
        };
        out.push_str(&format!(
            "{program}\t{}\t{}\t{}\t{}\t{}\t{decision}\n",
            format_float(row.observed_stability),
            format_float(row.null_stability_mean),
            format_float(row.null_stability_p95),
            format_float(row.null_p),
            format_float(q_val),
        ));
    }
    atomic_write(&args.output, out.as_bytes())?;

    eprintln!(
        "decompose null: wrote {} programs to {}",
        rows.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let canonical_inputs: &[&str] = &["measurements.tsv", "samples.tsv", "proteins.tsv"];
    let inputs_sha256 = hash_canonical_inputs(&args.input_dir, canonical_inputs)?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "decompose null",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output": args.output.display().to_string(),
            "k": args.k,
            "n-perm": args.n_perm,
            "n-seeds": args.n_seeds,
            "seed": args.seed,
            "null-mode": args.null_mode.to_core().as_str(),
            "top-n": args.top_n,
            "max-iter": args.max_iter,
            "tol": args.tol,
            "max-missing-fraction": args.max_missing_fraction,
            "impute": args.impute,
            "q-threshold": args.q_threshold,
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("decompose null: sidecar={}", sidecar.display());
    Ok(())
}
