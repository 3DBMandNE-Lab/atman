use anyhow::{bail, Context, Result};
use atman_core::bh_fdr;
use atman_core::compositional::{apply_transform, Transform};
use atman_core::decompose_unmix::{
    bootstrap_ci, select_k_auto, unmix, AbundanceCi, AbundanceMethod, EndmemberMethod,
    KSweepRow, LoadingCi, UnmixResult,
};
use atman_core::ica::{fast_ica, jaccard_top_n, pca_whiten, select_k_cumulative_variance, IcaResult};
use atman_core::ica_null::{archetype_null, ArchetypeNullRow, NullMode, NullParams};
use atman_core::variance_decomposition::{
    decompose_archetype_variance, FixedFactor, VarianceRow,
};
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use csv::ReaderBuilder;
use serde_json::json;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, format_float, hash_canonical_inputs, read_measurements_long,
    read_samples, sidecar_path_for, write_run_sidecar,
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
    /// Partition archetype activation variance into fixed-effect,
    /// random-intercept, and residual components per archetype.
    Variance(VarianceArgs),
    /// Geometric compartmental unmixing (VCA endmember extraction
    /// + FCLS/UCLS abundance estimation).
    Unmix(UnmixArgs),
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
        Command::Variance(args) => run_variance(args),
        Command::Unmix(args) => run_unmix(args),
    }
}

#[derive(ClapArgs, Debug)]
pub struct VarianceArgs {
    /// Per-sample activation TSV produced by `atman decompose ica`:
    /// columns `cohort, subject_id, sample_id, program, activation`.
    #[arg(long)]
    activations: PathBuf,

    /// samples.tsv providing per-sample metadata columns used as
    /// fixed factors or the random-intercept grouping.
    #[arg(long)]
    samples: PathBuf,

    /// Mixed-model formula: `<fixed_terms> [+ (1|group_col)]`.
    /// Fixed terms are comma- or plus-separated column names from
    /// `samples.tsv` (or the built-ins `cohort`, `condition`).
    /// Categorical columns are one-hot encoded (alphabetically
    /// first level dropped as reference); numeric columns are single
    /// coefficient columns. At most one random-intercept term is
    /// supported; the grouping column must exist in samples.tsv.
    /// Example: `--factors "cohort + condition + (1|subject_id)"`.
    #[arg(long)]
    factors: String,

    /// Minimum samples per archetype for the fit to run. Below this
    /// threshold, the archetype emits a skip row with NaN stats.
    #[arg(long, default_value_t = 5)]
    min_samples: usize,

    /// Output TSV path.
    #[arg(long)]
    output: PathBuf,
}

fn run_variance(args: VarianceArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let samples = read_samples(&args.samples)?;
    if samples.is_empty() {
        bail!("samples.tsv at {:?} is empty", args.samples);
    }
    let (fixed_terms, random_group) = parse_variance_formula(&args.factors)?;
    if fixed_terms.is_empty() {
        bail!("--factors must include at least one fixed term");
    }

    // Read per-sample extras from the raw samples.tsv file (headers
    // include user-defined covariates beyond the canonical cols).
    let extras = read_samples_extras(&args.samples)?;

    // Activations: one row per (sample_id, program). Re-shape to
    // per-program per-sample.
    let (program_by_sample, programs_order, samples_order) =
        read_activations(&args.activations)?;
    if programs_order.is_empty() {
        bail!("no programs found in {:?}", args.activations);
    }
    if samples_order.is_empty() {
        bail!("no samples found in {:?}", args.activations);
    }

    // Build the design matrix and fixed-factor column ranges. Only
    // samples that appear in both `activations` and `samples.tsv`
    // are kept, and only samples with finite covariates.
    let (design, fixed_factors, kept_sample_ids, group_labels) =
        build_variance_design(
            &samples,
            &extras,
            &fixed_terms,
            random_group.as_deref(),
            &samples_order,
        )?;
    if design.is_empty() {
        bail!("no samples remain after fixed-covariate complete-case filtering");
    }

    // Build per-program activation rows in the reduced sample order.
    let mut activations_mat: Vec<Vec<f64>> = Vec::with_capacity(programs_order.len());
    for program in &programs_order {
        let mut row = Vec::with_capacity(kept_sample_ids.len());
        for sid in &kept_sample_ids {
            let v = program_by_sample
                .get(&(program.clone(), sid.clone()))
                .copied()
                .unwrap_or(f64::NAN);
            row.push(v);
        }
        activations_mat.push(row);
    }

    let rows: Vec<VarianceRow> = decompose_archetype_variance(
        &programs_order,
        &activations_mat,
        &design,
        &fixed_factors,
        group_labels.as_deref(),
        args.min_samples,
    )
    .map_err(|e| anyhow::anyhow!(e))?;

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output dir {:?}", parent))?;
    }
    write_variance_tsv(&args.output, &rows, &fixed_factors, random_group.as_deref())?;
    eprintln!(
        "decompose variance: wrote {} archetypes to {}",
        rows.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let input_dir_sha256 = crate::io::hash_labeled_inputs(&[
        ("activations", args.activations.as_path()),
        ("samples", args.samples.as_path()),
    ])?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "decompose variance",
        json!({
            "activations": args.activations.display().to_string(),
            "samples": args.samples.display().to_string(),
            "factors": args.factors,
            "min-samples": args.min_samples,
            "output": args.output.display().to_string(),
        }),
        &input_dir_sha256,
        &[args.output.clone()],
        started_at,
        finished_at,
    )?;
    eprintln!("decompose variance: sidecar={}", sidecar.display());
    Ok(())
}

fn parse_variance_formula(formula: &str) -> Result<(Vec<String>, Option<String>)> {
    let mut fixed = Vec::new();
    let mut random: Option<String> = None;
    for raw in formula.split(|c: char| c == '+' || c == ',') {
        let term = raw.trim();
        if term.is_empty() {
            continue;
        }
        if let Some(rest) = term
            .strip_prefix('(')
            .and_then(|s| s.strip_suffix(')'))
        {
            // (1|group_col)
            let inner = rest.trim();
            if let Some(group) = inner.strip_prefix("1|") {
                let g = group.trim();
                if g.is_empty() {
                    bail!("empty random-intercept grouping column in formula");
                }
                if random.is_some() {
                    bail!("at most one random-intercept term is supported");
                }
                random = Some(g.to_string());
            } else {
                bail!(
                    "unsupported random-effect term {:?}; expected `(1|<column>)`",
                    term
                );
            }
        } else {
            fixed.push(term.to_string());
        }
    }
    Ok((fixed, random))
}

/// Load the user-visible samples.tsv extra columns (beyond the
/// canonical `sample_id, subject_id, condition, is_control,
/// sample_type, ingest_order` that `read_samples` returns). Needed
/// for variance-decomposition formulas that reference cohort / sex
/// / batch / etc.
fn read_samples_extras(path: &Path) -> Result<BTreeMap<String, BTreeMap<String, String>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers: Vec<String> = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let sample_col = headers
        .iter()
        .position(|h| h == "sample_id")
        .ok_or_else(|| anyhow::anyhow!("samples.tsv {:?} missing sample_id", path))?;
    let mut out: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
    for row in reader.records() {
        let row = row.with_context(|| format!("reading {:?}", path))?;
        let sample = row
            .get(sample_col)
            .unwrap_or_default()
            .to_string();
        if sample.is_empty() {
            continue;
        }
        let mut extras: BTreeMap<String, String> = BTreeMap::new();
        for (i, h) in headers.iter().enumerate() {
            if i == sample_col {
                continue;
            }
            let v = row.get(i).unwrap_or_default().to_string();
            extras.insert(h.clone(), v);
        }
        out.insert(sample, extras);
    }
    Ok(out)
}

fn read_activations(
    path: &Path,
) -> Result<(BTreeMap<(String, String), f64>, Vec<String>, Vec<String>)> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers: Vec<String> = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let col = |name: &str| -> Result<usize> {
        headers
            .iter()
            .position(|h| h == name)
            .ok_or_else(|| anyhow::anyhow!("activations.tsv {:?} missing column {:?}", path, name))
    };
    let c_sample = col("sample_id")?;
    let c_program = col("program")?;
    let c_activation = col("activation")?;
    let mut by_key: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut programs: BTreeSet<String> = BTreeSet::new();
    let mut samples_seen: BTreeSet<String> = BTreeSet::new();
    let mut samples_order: Vec<String> = Vec::new();
    for row in reader.records() {
        let row = row.with_context(|| format!("reading {:?}", path))?;
        let sample = row.get(c_sample).unwrap_or_default().to_string();
        let program = row.get(c_program).unwrap_or_default().to_string();
        let value: f64 = row
            .get(c_activation)
            .unwrap_or_default()
            .parse()
            .unwrap_or(f64::NAN);
        if sample.is_empty() || program.is_empty() {
            continue;
        }
        if samples_seen.insert(sample.clone()) {
            samples_order.push(sample.clone());
        }
        programs.insert(program.clone());
        by_key.insert((program, sample), value);
    }
    let programs_order: Vec<String> = programs.into_iter().collect();
    Ok((by_key, programs_order, samples_order))
}

/// Build the design matrix + per-factor column ranges + per-sample
/// random-group labels in the order of `samples_order`. Drops samples
/// whose fixed covariates can't be resolved.
fn build_variance_design(
    samples: &[atman_core::Sample],
    extras: &BTreeMap<String, BTreeMap<String, String>>,
    fixed_terms: &[String],
    random_group: Option<&str>,
    samples_order: &[String],
) -> Result<(Vec<Vec<f64>>, Vec<FixedFactor>, Vec<String>, Option<Vec<String>>)> {
    let sample_by_id: BTreeMap<&str, &atman_core::Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    // For each fixed term, determine whether it's numeric or
    // categorical by probing values across samples_order.
    fn value_of(
        s: &atman_core::Sample,
        extras: &BTreeMap<String, String>,
        term: &str,
    ) -> Option<String> {
        match term {
            "condition" => s.condition.clone(),
            "is_control" => Some(if s.is_control { "1".into() } else { "0".into() }),
            "sample_type" => s.sample_type.clone(),
            "ingest_order" => Some(s.ingest_order.to_string()),
            "subject_id" => s.subject_id.clone(),
            other => extras.get(other).cloned(),
        }
    }

    // Pass 1: collect all values per term; decide type (numeric if
    // every non-empty value parses as f64, else categorical).
    struct TermInfo {
        name: String,
        numeric: bool,
        // For categorical: alphabetically sorted non-reference levels.
        levels: Vec<String>,
    }
    let mut term_info: Vec<TermInfo> = Vec::new();
    for term in fixed_terms {
        let mut seen: BTreeSet<String> = BTreeSet::new();
        let mut all_numeric = true;
        for sid in samples_order {
            let s = match sample_by_id.get(sid.as_str()) {
                Some(s) => *s,
                None => continue,
            };
            let extras_for_s = extras.get(sid).cloned().unwrap_or_default();
            if let Some(v) = value_of(s, &extras_for_s, term) {
                if v.is_empty() {
                    continue;
                }
                if v.parse::<f64>().is_err() {
                    all_numeric = false;
                }
                seen.insert(v);
            }
        }
        let levels: Vec<String> = seen.into_iter().collect();
        if levels.is_empty() {
            bail!("fixed term {:?} has no observed values", term);
        }
        term_info.push(TermInfo {
            name: term.clone(),
            numeric: all_numeric,
            levels,
        });
    }

    // Pass 2: determine column layout. Column 0 is intercept. Each
    // numeric term is 1 column. Each categorical term drops the
    // alphabetically-first level as reference and uses `k-1` columns.
    let mut fixed_factors: Vec<FixedFactor> = Vec::new();
    let mut total_cols = 1; // intercept
    // Intercept isn't a named factor row in the output; factor rows
    // are for the user-requested terms.
    for t in &term_info {
        let width = if t.numeric { 1 } else { t.levels.len().saturating_sub(1) };
        if width == 0 {
            bail!(
                "fixed term {:?} has only one level; cannot contribute variance",
                t.name
            );
        }
        fixed_factors.push(FixedFactor {
            name: t.name.clone(),
            columns: total_cols..(total_cols + width),
        });
        total_cols += width;
    }

    let mut design: Vec<Vec<f64>> = Vec::new();
    let mut kept: Vec<String> = Vec::new();
    let mut groups: Vec<String> = Vec::new();
    for sid in samples_order {
        let s = match sample_by_id.get(sid.as_str()) {
            Some(s) => *s,
            None => continue,
        };
        let extras_for_s = extras.get(sid).cloned().unwrap_or_default();
        let mut row = vec![0.0_f64; total_cols];
        row[0] = 1.0; // intercept
        let mut complete = true;
        for (ti, info) in term_info.iter().enumerate() {
            let factor = &fixed_factors[ti];
            let v = value_of(s, &extras_for_s, &info.name).unwrap_or_default();
            if v.is_empty() {
                complete = false;
                break;
            }
            if info.numeric {
                let parsed: f64 = match v.parse() {
                    Ok(x) => x,
                    Err(_) => {
                        complete = false;
                        break;
                    }
                };
                row[factor.columns.start] = parsed;
            } else {
                // Reference level = levels[0]; other levels get a 1 in
                // their respective column.
                let ref_level = &info.levels[0];
                if v == *ref_level {
                    // All zeros is correct.
                } else {
                    let position = info.levels[1..].iter().position(|l| l == &v);
                    match position {
                        Some(pos) => row[factor.columns.start + pos] = 1.0,
                        None => {
                            complete = false;
                            break;
                        }
                    }
                }
            }
        }
        if !complete {
            continue;
        }
        if let Some(rg) = random_group {
            let v = value_of(s, &extras_for_s, rg).unwrap_or_default();
            if v.is_empty() {
                continue;
            }
            groups.push(v);
        }
        design.push(row);
        kept.push(sid.clone());
    }
    let group_labels = if random_group.is_some() { Some(groups) } else { None };
    Ok((design, fixed_factors, kept, group_labels))
}

/// 12-decimal float formatter for variance-decomposition outputs so
/// parity tests against R lm + Type III can assert sub-1e-8 drift.
/// The shared `format_float` truncates to 6 decimals which is too
/// coarse for parity on statistics that run into the hundreds (F,
/// SS) when scaled.
fn format_stat(v: f64) -> String {
    if !v.is_finite() {
        if v.is_nan() {
            "NaN".into()
        } else if v > 0.0 {
            "Inf".into()
        } else {
            "-Inf".into()
        }
    } else if v == 0.0 {
        "0".into()
    } else {
        format!("{v:.12}")
    }
}

fn write_variance_tsv(
    path: &Path,
    rows: &[VarianceRow],
    factors: &[FixedFactor],
    random_group: Option<&str>,
) -> Result<()> {
    let mut header = String::from("archetype_id\ttotal_var\tvar_residual");
    if let Some(g) = random_group {
        header.push_str(&format!("\tvar_random_{g}\ticc_random_{g}"));
    }
    for f in factors {
        // Type III suite per factor: SS_III, F, df_num, df_den, p,
        // plus the per-coefficient max|t|/min p summary for quick
        // inspection.
        header.push_str(&format!(
            "\tss_type3_{name}\tf_statistic_{name}\tdf_num_{name}\tdf_den_{name}\tp_value_{name}\tmax_abs_t_{name}\tmin_p_{name}",
            name = f.name
        ));
    }
    header.push_str("\tskip_reason\n");

    let mut buf = header;
    for r in rows {
        buf.push_str(&r.archetype_id);
        buf.push('\t');
        buf.push_str(&format_stat(r.total_var));
        buf.push('\t');
        buf.push_str(&format_stat(r.var_residual));
        if random_group.is_some() {
            buf.push('\t');
            buf.push_str(&format_stat(r.var_random.unwrap_or(f64::NAN)));
            buf.push('\t');
            buf.push_str(&format_stat(r.icc_random.unwrap_or(f64::NAN)));
        }
        for (i, _) in factors.iter().enumerate() {
            let f = &r.per_factor[i];
            buf.push('\t');
            buf.push_str(&format_stat(f.ss_type3));
            buf.push('\t');
            buf.push_str(&format_stat(f.f_statistic));
            buf.push('\t');
            buf.push_str(&f.df_num.to_string());
            buf.push('\t');
            buf.push_str(&format_stat(f.df_den));
            buf.push('\t');
            buf.push_str(&format_stat(f.p_value));
            buf.push('\t');
            buf.push_str(&format_stat(f.max_abs_t));
            buf.push('\t');
            buf.push_str(&format_stat(f.min_p));
        }
        buf.push('\t');
        buf.push_str(r.skip_reason.as_deref().unwrap_or(""));
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
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

#[derive(ClapArgs, Debug)]
pub struct UnmixArgs {
    /// Canonical Atman input directory (expects `qc_measurements.tsv`
    /// or `measurements.tsv`, plus `samples.tsv`, `proteins.tsv`).
    #[arg(long)]
    input_dir: PathBuf,

    /// Number of endmembers. Pass an integer for a fixed `k`, or
    /// `auto` to sweep `[--k-min, --k-max]` and pick the smallest
    /// `k` whose marginal reconstruction-residual improvement drops
    /// below `--k-elbow-threshold` times the sweep's peak
    /// improvement. Must satisfy `2 ≤ k ≤ n_samples / 2`.
    #[arg(long)]
    k: String,

    /// Lower bound for `--k auto` sweep.
    #[arg(long, default_value_t = 2)]
    k_min: usize,

    /// Upper bound for `--k auto` sweep.
    #[arg(long, default_value_t = 8)]
    k_max: usize,

    /// Elbow threshold for `--k auto`: marginal improvement / peak
    /// improvement ratio at which the sweep stops.
    #[arg(long, default_value_t = 0.10)]
    k_elbow_threshold: f64,

    /// Endmember extraction method. `vca` (default): deterministic
    /// projection-onto-complement. `nfindr`: iterative simplex-volume
    /// maximization initialized from VCA.
    #[arg(long, default_value = "vca")]
    method: String,

    /// Max N-FINDR passes (ignored for `--method vca`). A full pass
    /// tries every (endmember slot × candidate sample) swap.
    #[arg(long, default_value_t = 20)]
    nfindr_max_passes: usize,

    /// Abundance estimator. `fcls` (default) enforces simplex
    /// constraints `α ≥ 0 ∧ Σα = 1`. `ucls` drops them.
    #[arg(long, default_value = "fcls")]
    abundance: String,

    /// Pre-transform applied to the subject × protein matrix before
    /// unmixing. `none` passes atman's already-log-scale canonical
    /// data through; `clr` centers each subject row on its own mean
    /// (valid only with `--abundance ucls`, or with
    /// `--allow-unconstrained-simplex` for FCLS); `ilr` is refused
    /// because it changes the feature ordering.
    #[arg(long, default_value = "none")]
    transform: String,

    /// Reference gene symbol for `--transform alr` / `ratio-anchor`.
    #[arg(long)]
    alr_reference: Option<String>,

    /// Escape hatch: run FCLS even with a compositional transform
    /// whose simplex interpretation is debatable (CLR, ILR). Off by
    /// default — FCLS's sum-to-one constraint has no natural
    /// meaning on log-ratio coordinates.
    #[arg(long, default_value_t = false)]
    allow_unconstrained_simplex: bool,

    /// Canonical measurements file: `qc` (default) or `raw`.
    #[arg(long, default_value = "qc")]
    source: String,

    /// Seed for VCA's initial projection direction. Sub-seeds for
    /// each VCA step derive from this via SplitMix64.
    #[arg(long, default_value_t = 20260420)]
    seed: u64,

    /// FCLS maximum projected-gradient iterations.
    #[arg(long, default_value_t = 1000)]
    fcls_max_iter: usize,

    /// FCLS convergence tolerance (`max |Δα|` between iterations).
    #[arg(long, default_value_t = 1e-9)]
    fcls_tol: f64,

    /// Subject-level bootstrap iterations for loading + abundance CI.
    /// 0 (default) disables bootstrap — endmembers.tsv and
    /// abundances.tsv carry `NA` in the CI columns.
    #[arg(long, default_value_t = 0)]
    n_boot: usize,

    /// Drop proteins with more than this fraction of missing samples.
    #[arg(long, default_value_t = 0.0)]
    max_missing_fraction: f64,

    /// `--impute none|mean` for residual missingness.
    #[arg(long, default_value = "none")]
    impute: String,

    /// Output directory for `endmembers.tsv`, `abundances.tsv`,
    /// `unmix_diagnostics.tsv`.
    #[arg(long)]
    output_dir: PathBuf,
}

fn run_unmix(args: UnmixArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let endmember_method = match args.method.as_str() {
        "vca" => EndmemberMethod::Vca,
        "nfindr" => EndmemberMethod::Nfindr {
            max_passes: args.nfindr_max_passes,
        },
        other => bail!("--method {other:?}; expected `vca` or `nfindr`"),
    };
    let abundance_method = match args.abundance.as_str() {
        "fcls" => AbundanceMethod::Fcls,
        "ucls" => AbundanceMethod::Ucls,
        other => bail!("--abundance {other:?}; expected fcls or ucls"),
    };
    match args.transform.as_str() {
        "none" | "log" => {}
        "clr" | "alr" | "ratio-anchor" => {
            if matches!(abundance_method, AbundanceMethod::Fcls)
                && !args.allow_unconstrained_simplex
            {
                bail!(
                    "--transform {:?} combined with --abundance fcls is refused: \
                     log-ratio coordinates do not admit a `Σα = 1` interpretation. \
                     Pass --allow-unconstrained-simplex to override, or switch to \
                     --abundance ucls.",
                    args.transform
                );
            }
        }
        "ilr" => {
            bail!(
                "--transform ilr is refused in `decompose unmix`: ILR changes \
                 the feature ordering, so recovered endmember loadings would \
                 not line up with the input protein labels"
            );
        }
        other => bail!("--transform {other:?}; expected none|log|clr|alr|ratio-anchor"),
    }
    let impute_mean = match args.impute.as_str() {
        "none" => false,
        "mean" => true,
        other => bail!("--impute {other:?}; expected none or mean"),
    };
    if args.fcls_max_iter < 1 {
        bail!("--fcls-max-iter must be >= 1");
    }

    // Load subject × protein matrix, keyed by gene_symbol.
    let source_file = match args.source.as_str() {
        "qc" => "qc_measurements.tsv",
        "raw" => "measurements.tsv",
        other => bail!("--source {other:?}; expected qc or raw"),
    };
    let (subject_ids, protein_labels, data) = load_subject_protein_matrix(
        &args.input_dir.join(source_file),
        args.max_missing_fraction,
        impute_mean,
    )?;
    // Resolve --k: either an integer or "auto" (sweep + elbow).
    let (effective_k, k_sweep): (usize, Vec<KSweepRow>) = if args.k == "auto" {
        if args.k_min < 2 || args.k_max < args.k_min {
            bail!(
                "invalid --k auto sweep: k-min={} k-max={}",
                args.k_min, args.k_max
            );
        }
        if args.k_max * 2 > subject_ids.len() {
            bail!(
                "--k-max={} requires >= {} samples; got {}",
                args.k_max,
                args.k_max * 2,
                subject_ids.len()
            );
        }
        eprintln!(
            "decompose unmix: --k auto sweeping [{}, {}]",
            args.k_min, args.k_max
        );
        // Apply transform later; here we need it already for the
        // sweep's VCA+FCLS runs. Rebuild the transformed matrix
        // early.
        // (The transform below is still the canonical one — we
        // just need it before the sweep runs.)
        select_k_auto(
            // placeholder: applied below after transform
            &data, args.k_min, args.k_max, args.seed,
            EndmemberMethod::Vca,  // sweep uses VCA; caller can
                                   // re-run with nfindr post-hoc
            abundance_method,
            args.fcls_max_iter,
            args.fcls_tol,
            args.k_elbow_threshold,
        )
        .map_err(|e| anyhow::anyhow!(e))?
    } else {
        let k: usize = args.k.parse().with_context(|| {
            format!("--k {:?}: expected an integer or `auto`", args.k)
        })?;
        (k, Vec::new())
    };
    if subject_ids.len() < effective_k * 2 {
        bail!(
            "k={} requires >= {} samples (k > n/2 is under-determined); got {}",
            effective_k,
            effective_k * 2,
            subject_ids.len()
        );
    }
    // Apply transform (if any) — "log" is a no-op on atman canonical
    // data since it's already on a log scale.
    let transform = match args.transform.as_str() {
        "none" | "log" => Transform::None,
        "clr" => Transform::Clr,
        "alr" => {
            let reference = args
                .alr_reference
                .as_deref()
                .context("--transform alr requires --alr-reference <gene>")?;
            let idx = protein_labels
                .iter()
                .position(|l| l == reference)
                .with_context(|| {
                    format!("--alr-reference {reference:?} not in protein labels")
                })?;
            Transform::Alr { reference_index: idx }
        }
        "ratio-anchor" => {
            let reference = args
                .alr_reference
                .as_deref()
                .context("--transform ratio-anchor requires --alr-reference <gene>")?;
            let idx = protein_labels
                .iter()
                .position(|l| l == reference)
                .with_context(|| {
                    format!("--alr-reference {reference:?} not in protein labels")
                })?;
            Transform::RatioAnchor { reference_index: idx }
        }
        _ => unreachable!(),
    };
    let transformed = apply_transform(&data, transform)
        .map_err(|e| anyhow::anyhow!("transform failed: {e}"))?;

    let result = unmix(
        &transformed,
        effective_k,
        args.seed,
        endmember_method,
        abundance_method,
        args.fcls_max_iter,
        args.fcls_tol,
    )
    .map_err(|e| anyhow::anyhow!(e))?;
    let (loading_ci, abundance_ci) = bootstrap_ci(
        &result,
        &transformed,
        args.seed,
        endmember_method,
        abundance_method,
        args.fcls_max_iter,
        args.fcls_tol,
        args.n_boot,
    )
    .map_err(|e| anyhow::anyhow!(e))?;

    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating {:?}", args.output_dir))?;
    let endmembers_path = args.output_dir.join("endmembers.tsv");
    write_unmix_endmembers(
        &endmembers_path,
        &result,
        &protein_labels,
        &subject_ids,
        loading_ci.as_ref(),
    )?;
    let abundances_path = args.output_dir.join("abundances.tsv");
    write_unmix_abundances(&abundances_path, &result, &subject_ids, abundance_ci.as_ref())?;
    let diag_path = args.output_dir.join("unmix_diagnostics.tsv");
    write_unmix_diagnostics(&diag_path, &result, &subject_ids)?;
    let mut outputs = vec![endmembers_path.clone(), abundances_path.clone(), diag_path.clone()];
    if !k_sweep.is_empty() {
        let k_path = args.output_dir.join("k_selection.tsv");
        write_k_selection(&k_path, &k_sweep, effective_k)?;
        outputs.push(k_path);
    }

    eprintln!(
        "decompose unmix: k={} n={} p={} method={} abundance={} transform={}",
        effective_k,
        subject_ids.len(),
        protein_labels.len(),
        args.method,
        args.abundance,
        args.transform,
    );

    let finished_at = SystemTime::now();
    let input_dir_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &[
            "qc_measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
    )?;
    let sidecar = sidecar_path_for(&endmembers_path);
    write_run_sidecar(
        &sidecar,
        "decompose unmix",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "k": args.k,
            "effective-k": effective_k,
            "k-min": args.k_min,
            "k-max": args.k_max,
            "k-elbow-threshold": args.k_elbow_threshold,
            "method": args.method,
            "nfindr-max-passes": args.nfindr_max_passes,
            "abundance": args.abundance,
            "transform": args.transform,
            "alr-reference": args.alr_reference,
            "allow-unconstrained-simplex": args.allow_unconstrained_simplex,
            "source": args.source,
            "seed": args.seed,
            "fcls-max-iter": args.fcls_max_iter,
            "fcls-tol": args.fcls_tol,
            "n-boot": args.n_boot,
            "max-missing-fraction": args.max_missing_fraction,
            "impute": args.impute,
        }),
        &input_dir_sha256,
        &outputs,
        started_at,
        finished_at,
    )?;
    eprintln!("decompose unmix: sidecar={}", sidecar.display());
    Ok(())
}

fn write_k_selection(path: &Path, rows: &[KSweepRow], chosen_k: usize) -> Result<()> {
    let mut buf = String::from(
        "k\tmean_residual_norm\tmarginal_improvement\tchosen\n",
    );
    for r in rows {
        buf.push_str(&format!(
            "{}\t{:.6}\t{:.6}\t{}\n",
            r.k,
            r.mean_residual_norm,
            r.marginal_improvement,
            if r.k == chosen_k { 1 } else { 0 },
        ));
    }
    atomic_write(path, buf.as_bytes())
}

fn load_subject_protein_matrix(
    path: &Path,
    max_missing_fraction: f64,
    impute_mean: bool,
) -> Result<(Vec<String>, Vec<String>, Vec<Vec<f64>>)> {
    let records = read_measurements_long(path)?;
    let mut samples: BTreeSet<String> = BTreeSet::new();
    let mut genes: BTreeSet<String> = BTreeSet::new();
    let mut cells: BTreeMap<(String, String), f64> = BTreeMap::new();
    for r in &records {
        if r.dropped_by_qc {
            continue;
        }
        let v = r.abundance.as_f64();
        if !v.is_finite() {
            continue;
        }
        let gene = match &r.gene_symbol {
            Some(g) => g.clone(),
            None => continue,
        };
        samples.insert(r.sample_id.clone());
        genes.insert(gene.clone());
        cells.insert((gene, r.sample_id.clone()), v);
    }
    let sample_list: Vec<String> = samples.into_iter().collect();
    let gene_list: Vec<String> = genes.into_iter().collect();
    if sample_list.is_empty() || gene_list.is_empty() {
        bail!("no usable measurements in {:?}", path);
    }
    // Build with missingness map.
    let mut raw = vec![vec![f64::NAN; gene_list.len()]; sample_list.len()];
    for (si, sid) in sample_list.iter().enumerate() {
        for (gi, gene) in gene_list.iter().enumerate() {
            if let Some(&v) = cells.get(&(gene.clone(), sid.clone())) {
                raw[si][gi] = v;
            }
        }
    }
    // Drop features exceeding the missing-fraction budget.
    let n = sample_list.len() as f64;
    let keep: Vec<bool> = (0..gene_list.len())
        .map(|gi| {
            let present = (0..sample_list.len())
                .filter(|&si| raw[si][gi].is_finite())
                .count();
            let missing = 1.0 - (present as f64 / n);
            missing <= max_missing_fraction + 1e-12
        })
        .collect();
    let kept_genes: Vec<String> = gene_list
        .iter()
        .zip(keep.iter())
        .filter_map(|(g, k)| if *k { Some(g.clone()) } else { None })
        .collect();
    if kept_genes.is_empty() {
        bail!(
            "no proteins retained after --max-missing-fraction {} filter",
            max_missing_fraction
        );
    }
    let mut data = vec![vec![0.0_f64; kept_genes.len()]; sample_list.len()];
    for (new_gi, gi) in (0..gene_list.len()).filter(|i| keep[*i]).enumerate() {
        // Compute column mean for imputation.
        let mut sum = 0.0;
        let mut count = 0usize;
        for si in 0..sample_list.len() {
            if raw[si][gi].is_finite() {
                sum += raw[si][gi];
                count += 1;
            }
        }
        let mean = if count > 0 { sum / count as f64 } else { 0.0 };
        for si in 0..sample_list.len() {
            let v = raw[si][gi];
            if v.is_finite() {
                data[si][new_gi] = v;
            } else if impute_mean {
                data[si][new_gi] = mean;
            } else {
                bail!(
                    "missing value at sample {} protein {}; rerun with --impute mean",
                    sample_list[si],
                    kept_genes[new_gi]
                );
            }
        }
    }
    Ok((sample_list, kept_genes, data))
}

fn write_unmix_endmembers(
    path: &Path,
    result: &UnmixResult,
    protein_labels: &[String],
    subject_ids: &[String],
    loading_ci: Option<&LoadingCi>,
) -> Result<()> {
    let mut buf = String::from(
        "endmember_id\tsource_sample_id\tprotein\tloading\trank_in_endmember\t\
         loading_ci_lower\tloading_ci_upper\n",
    );
    for (ai, loading) in result.endmember_loadings.iter().enumerate() {
        let src_idx = result.endmember_sample_indices[ai];
        let source_sid = subject_ids
            .get(src_idx)
            .cloned()
            .unwrap_or_else(|| format!("sample{src_idx}"));
        let mut ranking: Vec<(usize, f64)> = loading
            .iter()
            .enumerate()
            .map(|(i, v)| (i, v.abs()))
            .collect();
        ranking.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
        let mut rank_for = vec![0usize; loading.len()];
        for (r, (i, _)) in ranking.iter().enumerate() {
            rank_for[*i] = r + 1;
        }
        for (pi, &v) in loading.iter().enumerate() {
            let (lo, hi) = loading_ci
                .map(|ci| (ci.lower[ai][pi], ci.upper[ai][pi]))
                .map(|(a, b)| (format!("{a:.6}"), format!("{b:.6}")))
                .unwrap_or_else(|| ("NA".into(), "NA".into()));
            buf.push_str(&format!(
                "E{:03}\t{}\t{}\t{:.6}\t{}\t{}\t{}\n",
                ai + 1,
                source_sid,
                protein_labels[pi],
                v,
                rank_for[pi],
                lo,
                hi,
            ));
        }
    }
    atomic_write(path, buf.as_bytes())
}

fn write_unmix_abundances(
    path: &Path,
    result: &UnmixResult,
    subject_ids: &[String],
    abundance_ci: Option<&AbundanceCi>,
) -> Result<()> {
    // Emit per-endmember `E001, E001_ci_lower, E001_ci_upper, E002, ...`
    // so downstream code can extract either the point estimate or the
    // interval by column name.
    let mut buf = String::from("sample_id");
    for ai in 0..result.endmember_loadings.len() {
        let lbl = format!("E{:03}", ai + 1);
        buf.push('\t');
        buf.push_str(&lbl);
        buf.push('\t');
        buf.push_str(&format!("{lbl}_ci_lower"));
        buf.push('\t');
        buf.push_str(&format!("{lbl}_ci_upper"));
    }
    buf.push('\n');
    for (si, sid) in subject_ids.iter().enumerate() {
        buf.push_str(sid);
        for (ai, &v) in result.abundances[si].iter().enumerate() {
            buf.push('\t');
            buf.push_str(&format!("{v:.6}"));
            let (lo, hi) = abundance_ci
                .map(|ci| (ci.lower[si][ai], ci.upper[si][ai]))
                .map(|(a, b)| (format!("{a:.6}"), format!("{b:.6}")))
                .unwrap_or_else(|| ("NA".into(), "NA".into()));
            buf.push('\t');
            buf.push_str(&lo);
            buf.push('\t');
            buf.push_str(&hi);
        }
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn write_unmix_diagnostics(
    path: &Path,
    result: &UnmixResult,
    subject_ids: &[String],
) -> Result<()> {
    let mut buf = String::from("sample_id\treconstruction_residual_norm\tabundance_sum\n");
    for (si, sid) in subject_ids.iter().enumerate() {
        let sum: f64 = result.abundances[si].iter().sum();
        buf.push_str(&format!(
            "{}\t{:.6}\t{:.6}\n",
            sid, result.residual_norms[si], sum,
        ));
    }
    atomic_write(path, buf.as_bytes())
}

