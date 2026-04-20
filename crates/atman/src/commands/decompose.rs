use anyhow::{bail, Context, Result};
use atman_core::bh_fdr;
use atman_core::compositional::{apply_transform, Transform};
use atman_core::ica::{fast_ica, jaccard_top_n, pca_whiten, select_k_cumulative_variance, IcaResult};
use atman_core::ica_null::{archetype_null, ArchetypeNullRow, NullMode, NullParams};
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use serde_json::json;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, format_float, hash_canonical_inputs, read_measurements_long, sidecar_path_for,
    write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Multi-seed FastICA decomposition with seed-stability reporting.
    Ica(IcaArgs),
    /// Permutation-null calibration of archetype stability.
    Null(NullArgs),
}

#[derive(ClapArgs, Debug)]
pub struct IcaArgs {
    /// Canonical Atman input directory (contains `qc_measurements.tsv` or `measurements.tsv`).
    #[arg(long)]
    input_dir: PathBuf,

    /// Explicit number of components. Overrides `--k-selection` if set.
    #[arg(long)]
    k: Option<usize>,

    /// K-selection method; currently supports `cumulative-variance=<target>`.
    #[arg(long, default_value = "cumulative-variance=0.80")]
    k_selection: String,

    /// Minimum allowed K (lower clamp for `--k-selection`).
    #[arg(long, default_value_t = 3)]
    k_min: usize,

    /// Maximum allowed K (upper clamp for `--k-selection`).
    #[arg(long, default_value_t = 30)]
    k_max: usize,

    /// Number of seeds (must be >= 1; first seed is the reference).
    #[arg(long, default_value_t = 50)]
    n_seeds: usize,

    /// First seed value; subsequent seeds are sequential (seed, seed+1, ...).
    #[arg(long, default_value_t = 20260418)]
    seed: u64,

    /// Jaccard threshold for a program to be considered "recovered" in another seed.
    #[arg(long, default_value_t = 0.9)]
    seed_stability_threshold: f64,

    /// Stability metric; currently `jaccard-top20` (i.e. top-N |loading| overlap).
    #[arg(long, value_enum, default_value_t = StabilityMetric::JaccardTop20)]
    stability_metric: StabilityMetric,

    /// Top-N loadings used by the stability metric.
    #[arg(long, default_value_t = 20)]
    stability_top_n: usize,

    /// Maximum FastICA iterations per seed.
    #[arg(long, default_value_t = 300)]
    max_iter: usize,

    /// FastICA convergence tolerance.
    #[arg(long, default_value_t = 1e-4)]
    tol: f64,

    /// Source for the abundance matrix: `qc` (default, `qc_measurements.tsv`) or `raw` (`measurements.tsv`).
    #[arg(long, default_value = "qc")]
    source: String,

    /// Drop assays where the fraction of missing samples exceeds this value. 0 = strict complete-case.
    #[arg(long, default_value_t = 0.0)]
    max_missing_fraction: f64,

    /// Imputation for residual missingness after the assay filter: `none` (fail) or `mean` (per-assay mean).
    #[arg(long, default_value = "none")]
    impute: String,

    /// Loadings output TSV path (`program\tassay_id\tgene_symbol\tloading`).
    #[arg(long)]
    output_loadings: PathBuf,

    /// Activations output TSV path (`cohort\tsubject_id\tsample_id\tprogram\tactivation`).
    #[arg(long)]
    output_activations: PathBuf,

    /// Stability output TSV path.
    #[arg(long)]
    output_stability: PathBuf,

    /// Cohort label written into the activations table (optional; defaults to directory name).
    #[arg(long)]
    cohort: Option<String>,

    /// Compositional transform applied to the sample × protein matrix
    /// before ICA. Atman's canonical abundance is already on a log
    /// scale, so these are linear operations on log values. `clr`
    /// (per-sample mean centering) is the recommended default for
    /// closed-sum proteomics. `alr` / `ratio-anchor` require
    /// `--alr-reference` (a gene symbol). `ilr` emits loadings in
    /// the Helmert coordinate basis rather than raw protein space.
    #[arg(long, value_enum, default_value_t = TransformArg::None)]
    transform: TransformArg,

    /// Gene symbol (matched against `samples` metadata `gene_symbol`)
    /// used as the reference for `--transform alr` or
    /// `--transform ratio-anchor`. Ignored for other transforms.
    #[arg(long)]
    alr_reference: Option<String>,
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

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Ica(args) => run_ica(args),
        Command::Null(args) => run_null(args),
    }
}

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
    /// Canonical Atman input directory (contains `qc_measurements.tsv` or `measurements.tsv`).
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

    /// Measurement source — `qc` reads `qc_measurements.tsv`, `raw`
    /// reads `measurements.tsv`. Defaults mirror `decompose ica`.
    #[arg(long, default_value = "qc")]
    source: String,

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

fn run_null(args: NullArgs) -> Result<()> {
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
        stability_metric: StabilityMetric::JaccardTop20,
        stability_top_n: args.top_n,
        max_iter: args.max_iter,
        tol: args.tol,
        source: args.source.clone(),
        max_missing_fraction: args.max_missing_fraction,
        impute: args.impute.clone(),
        output_loadings: PathBuf::new(),
        output_activations: PathBuf::new(),
        output_stability: PathBuf::new(),
        cohort: None,
        transform: TransformArg::None,
        alr_reference: None,
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
    let canonical_inputs: &[&str] = match args.source.as_str() {
        "qc" => &["qc_measurements.tsv", "samples.tsv", "proteins.tsv"],
        _ => &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    };
    let input_dir_sha256 = hash_canonical_inputs(&args.input_dir, canonical_inputs)?;
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
            "source": args.source,
            "max-missing-fraction": args.max_missing_fraction,
            "impute": args.impute,
            "q-threshold": args.q_threshold,
        }),
        &input_dir_sha256,
        &[args.output.clone()],
        started_at,
        finished_at,
    )?;
    eprintln!("decompose null: sidecar={}", sidecar.display());
    Ok(())
}

fn run_ica(args: IcaArgs) -> Result<()> {
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
    if args.stability_top_n == 0 {
        bail!("--stability-top-n must be at least 1");
    }

    let matrix = load_matrix(&args)?;
    let (matrix, transform_meta) = apply_compositional(matrix, &args)?;
    let k = resolve_k(&matrix, &args)?;
    if k > matrix.samples.len() || k > matrix.assays.len() {
        bail!(
            "k={k} exceeds min(n_samples={}, n_assays={})",
            matrix.samples.len(),
            matrix.assays.len()
        );
    }

    eprintln!(
        "decompose ica: n_samples={} n_assays={} k={} n_seeds={} threshold={} transform={}",
        matrix.samples.len(),
        matrix.assays.len(),
        k,
        args.n_seeds,
        args.seed_stability_threshold,
        transform_meta.name,
    );

    let mut runs = Vec::with_capacity(args.n_seeds);
    for offset in 0..args.n_seeds {
        let seed = args.seed.wrapping_add(offset as u64);
        let result = fast_ica(&matrix.data, k, seed, args.max_iter, args.tol);
        runs.push((seed, result));
    }

    let (ref_seed, ref_run) = runs.first().expect("at least one seed");
    let reference = canonicalize(ref_run);
    let cohort = args
        .cohort
        .clone()
        .unwrap_or_else(|| derive_cohort(&args.input_dir));

    write_loadings(&args.output_loadings, &matrix, &reference)?;
    write_activations(&args.output_activations, &matrix, &reference, &cohort)?;
    write_stability(
        &args.output_stability,
        &reference,
        &runs[1..],
        args.stability_top_n,
        args.seed_stability_threshold,
        *ref_seed,
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
    let canonical_inputs: &[&str] = match args.source.as_str() {
        "qc" => &["qc_measurements.tsv", "samples.tsv", "proteins.tsv"],
        _ => &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    };
    let input_dir_sha256 = hash_canonical_inputs(&args.input_dir, canonical_inputs)?;
    let sidecar = sidecar_path_for(&args.output_loadings);
    let stability_metric = match args.stability_metric {
        StabilityMetric::JaccardTop20 => "jaccard-top20",
    };
    let mut outputs = vec![
        args.output_loadings.clone(),
        args.output_activations.clone(),
        args.output_stability.clone(),
    ];
    if let Some(p) = &transform_applied_path {
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
            "stability-metric": stability_metric,
            "stability-top-n": args.stability_top_n,
            "max-iter": args.max_iter,
            "tol": args.tol,
            "source": args.source,
            "max-missing-fraction": args.max_missing_fraction,
            "impute": args.impute,
            "output-loadings": args.output_loadings.display().to_string(),
            "output-activations": args.output_activations.display().to_string(),
            "output-stability": args.output_stability.display().to_string(),
            "cohort": args.cohort,
            "transform": transform_meta.name,
            "alr-reference": transform_meta.alr_reference,
        }),
        &input_dir_sha256,
        &outputs,
        started_at,
        finished_at,
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
    let chosen = select_k_cumulative_variance(
        &whitening.full_spectrum,
        target,
        args.k_min,
        args.k_max.min(matrix.samples.len().saturating_sub(1)).max(1),
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
struct AbundanceMatrix {
    samples: Vec<SampleInfo>,
    assays: Vec<AssayInfo>,
    data: Vec<Vec<f64>>,
}

#[derive(Debug, Clone)]
struct SampleInfo {
    sample_id: String,
    subject_id: String,
}

#[derive(Debug, Clone)]
struct AssayInfo {
    assay_id: String,
    gene_symbol: String,
}

fn load_matrix(args: &IcaArgs) -> Result<AbundanceMatrix> {
    let filename = match args.source.as_str() {
        "qc" => "qc_measurements.tsv",
        "raw" => "measurements.tsv",
        other => bail!("--source {:?}; expected qc or raw", other),
    };
    let tsv = args.input_dir.join(filename);
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
        assay_meta.entry(assay.clone()).or_insert_with(|| r.gene_symbol.clone());
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
    if matches!(args.transform, TransformArg::Ilr)
        && args.alr_reference.is_some()
    {
        eprintln!("decompose ica: --alr-reference is ignored for --transform ilr");
    }
    let needs_reference = matches!(
        args.transform,
        TransformArg::Alr | TransformArg::RatioAnchor
    );
    let resolved_reference: Option<(usize, String)> = if needs_reference {
        let gene = match args.alr_reference.as_deref() {
            Some(g) if !g.is_empty() => g,
            _ => bail!(
                "--transform {transform_name} requires --alr-reference <gene_symbol>"
            ),
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
    let n = result.sources.len();

    let mut signs = vec![1.0_f64; k];
    let mut loadings: Vec<Vec<f64>> = vec![vec![0.0; p]; k];
    for c in 0..k {
        let mut max_abs = 0.0_f64;
        let mut argmax = 0usize;
        for (j, row) in result.mixing.iter().enumerate() {
            if row[c].abs() > max_abs {
                max_abs = row[c].abs();
                argmax = j;
            }
        }
        if result.mixing[argmax][c] < 0.0 {
            signs[c] = -1.0;
        }
        for j in 0..p {
            loadings[c][j] = signs[c] * result.mixing[j][c];
        }
    }
    let mut activations: Vec<Vec<f64>> = vec![vec![0.0; n]; k];
    for c in 0..k {
        for i in 0..n {
            activations[c][i] = signs[c] * result.sources[i][c];
        }
    }
    let mut order: Vec<usize> = (0..k).collect();
    order.sort_by(|&a, &b| {
        let ma = loadings[a]
            .iter()
            .fold(0.0_f64, |acc, v| acc.max(v.abs()));
        let mb = loadings[b]
            .iter()
            .fold(0.0_f64, |acc, v| acc.max(v.abs()));
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

fn program_name(index: usize) -> String {
    format!("program_{:02}", index + 1)
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

fn write_stability(
    path: &Path,
    reference: &CanonicalIca,
    others: &[(u64, IcaResult)],
    top_n: usize,
    threshold: f64,
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
            if best >= threshold {
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
        let flag = if n_alt > 0 && fraction < threshold_fraction(threshold) {
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

fn threshold_fraction(_threshold: f64) -> f64 {
    // A program is flagged unstable if fewer than this fraction of other-seed runs
    // achieve best-Jaccard >= threshold. We use 0.9 (matches the feature-request
    // default of "recover in at least 90% of seeds"); keeping it fixed avoids
    // introducing another user-facing knob for now.
    0.9
}

