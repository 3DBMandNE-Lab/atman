use anyhow::{bail, Context, Result};
use atman_core::align::{
    build_archetypes, summarize_archetypes, AlignMetric, AlignedProgram,
};
use atman_core::align_bootstrap::{
    align_bootstrap, BootstrapParams, BootstrapRow, CohortMatrix,
};
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use csv::ReaderBuilder;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, hash_labeled_inputs, read_measurements_long, sidecar_path_for,
    write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Cross-cohort program alignment and sweep mode.
    Programs(ProgramsArgs),
    /// Subject-level bootstrap of cross-cohort archetype alignment
    /// producing per-archetype universality probabilities.
    Bootstrap(BootstrapArgs),
}

#[derive(ClapArgs, Debug)]
pub struct ProgramsArgs {
    /// Comma-separated per-cohort loadings TSV paths.
    #[arg(long)]
    loadings: String,

    /// Comma-separated cohort labels (same order as --loadings). Defaults to file stems.
    #[arg(long)]
    cohorts: Option<String>,

    /// Comma-separated per-cohort annotations TSV paths (same order as --loadings).
    #[arg(long)]
    annotations: Option<String>,

    /// Label column in the loadings TSV: `gene_symbol` (default), `assay_id`, or `protein`.
    #[arg(long, default_value = "gene_symbol")]
    label_col: String,

    /// Similarity metric (single-mode). Ignored when --sweep is set.
    #[arg(long, value_enum, default_value_t = SingleMetric::Jaccard)]
    metric: SingleMetric,

    /// Top-N loadings for Jaccard (single-mode).
    #[arg(long, default_value_t = 40)]
    top_n: usize,

    /// Similarity threshold (single-mode).
    #[arg(long, default_value_t = 0.15)]
    tau: f64,

    /// Require reciprocal-best matches (single-mode and sweep).
    #[arg(long, default_value_t = false)]
    reciprocal_best: bool,

    /// Require same annotation category on both sides (single-mode).
    #[arg(long, default_value_t = false)]
    category_constraint: bool,

    /// Archetype output TSV (single-mode).
    #[arg(long)]
    output: Option<PathBuf>,

    /// Enable sweep mode; expects --metrics and --output-matrix.
    #[arg(long, default_value_t = false)]
    sweep: bool,

    /// Sweep spec: e.g. `jaccard:top_n=[20,40]:tau=[0.10,0.15],cosine:tau=[0.20,0.30]`.
    #[arg(long)]
    metrics: Option<String>,

    /// Emit both category-constrained and unconstrained rows for each sweep cell.
    #[arg(long, default_value_t = false)]
    compare_constrained_vs_unconstrained: bool,

    /// Sweep output matrix TSV.
    #[arg(long)]
    output_matrix: Option<PathBuf>,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum SingleMetric {
    Jaccard,
    Cosine,
    Spearman,
}

impl From<SingleMetric> for AlignMetric {
    fn from(m: SingleMetric) -> Self {
        match m {
            SingleMetric::Jaccard => AlignMetric::Jaccard,
            SingleMetric::Cosine => AlignMetric::Cosine,
            SingleMetric::Spearman => AlignMetric::Spearman,
        }
    }
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Programs(args) => run_programs(args),
        Command::Bootstrap(args) => run_bootstrap(args),
    }
}

#[derive(ClapArgs, Debug)]
pub struct BootstrapArgs {
    /// Comma-separated list of canonical Atman directories, one per
    /// cohort. Each must contain `qc_measurements.tsv`, `samples.tsv`,
    /// `proteins.tsv`. All cohorts share the same protein universe
    /// (intersection across cohorts is enforced).
    #[arg(long)]
    cohorts: String,

    /// Comma-separated cohort labels (same order as --cohorts).
    /// Defaults to directory basenames.
    #[arg(long)]
    labels: Option<String>,

    /// Number of ICA components per cohort.
    #[arg(long)]
    k: usize,

    /// Number of bootstrap iterations.
    #[arg(long, default_value_t = 100)]
    n_boot: usize,

    /// Top-level seed. Per-iteration sub-seeds derive deterministically
    /// from SplitMix64`(seed, iter)`.
    #[arg(long, default_value_t = 20260418)]
    seed: u64,

    /// Cosine similarity threshold for the alignment step that groups
    /// bootstrap programs into archetypes.
    #[arg(long, default_value_t = 0.30)]
    cosine_tau: f64,

    /// Floor cosine similarity between a point-estimate archetype's
    /// representative loading and a bootstrap program for that
    /// bootstrap archetype to be counted as a match.
    #[arg(long, default_value_t = 0.50)]
    match_tau: f64,

    /// Top-N loadings used by the Jaccard branch of the alignment
    /// (unused for cosine metric; kept for API symmetry).
    #[arg(long, default_value_t = 40)]
    top_n: usize,

    /// FastICA max iterations per seed per cohort.
    #[arg(long, default_value_t = 300)]
    max_iter: usize,

    /// FastICA convergence tolerance.
    #[arg(long, default_value_t = 1e-4)]
    tol: f64,

    /// Refuse when any cohort has fewer than this many subjects.
    #[arg(long, default_value_t = 20)]
    min_subjects: usize,

    /// Drop assays with more than this fraction of missing samples
    /// per cohort (matches `decompose ica`).
    #[arg(long, default_value_t = 0.0)]
    max_missing_fraction: f64,

    /// Imputation for residual missingness: `none` (fail), `mean`.
    #[arg(long, default_value = "none")]
    impute: String,

    /// Canonical input filename to read. `qc` (default) or `raw`.
    #[arg(long, default_value = "qc")]
    source: String,

    /// Output TSV path. Summary with one row per point-estimate
    /// archetype.
    #[arg(long)]
    output: PathBuf,
}

fn run_bootstrap(args: BootstrapArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if args.k == 0 {
        bail!("--k must be >= 1");
    }
    if args.n_boot == 0 {
        bail!("--n-boot must be >= 1");
    }
    if args.min_subjects < 2 {
        bail!("--min-subjects must be >= 2");
    }
    if !(0.0..=1.0).contains(&args.cosine_tau) {
        bail!("--cosine-tau must be in [0, 1]");
    }
    if !(0.0..=1.0).contains(&args.match_tau) {
        bail!("--match-tau must be in [0, 1]");
    }
    let cohort_dirs: Vec<PathBuf> = args
        .cohorts
        .split(',')
        .map(|s| PathBuf::from(s.trim()))
        .filter(|p| !p.as_os_str().is_empty())
        .collect();
    if cohort_dirs.len() < 2 {
        bail!("--cohorts must list at least 2 cohort directories");
    }
    let labels: Vec<String> = match &args.labels {
        Some(s) => s.split(',').map(|s| s.trim().to_string()).collect(),
        None => cohort_dirs
            .iter()
            .map(|p| file_stem(p).to_string())
            .collect(),
    };
    if labels.len() != cohort_dirs.len() {
        bail!(
            "--labels has {} entries but --cohorts has {}",
            labels.len(),
            cohort_dirs.len()
        );
    }
    let impute_mean = match args.impute.as_str() {
        "none" => false,
        "mean" => true,
        other => bail!("--impute {other:?}; expected none or mean"),
    };

    let matrices: Vec<CohortMatrix> = cohort_dirs
        .iter()
        .zip(labels.iter())
        .map(|(dir, label)| {
            load_cohort_matrix(
                dir,
                label,
                &args.source,
                args.max_missing_fraction,
                impute_mean,
            )
        })
        .collect::<Result<Vec<_>>>()?;
    // Enforce common protein universe across cohorts by intersecting
    // their label sets and restricting each matrix to the
    // intersection in a canonical order.
    let matrices = intersect_cohorts(matrices)?;

    eprintln!(
        "align bootstrap: cohorts={} k={} n_boot={} cosine_tau={} match_tau={} seed={}",
        matrices.len(),
        args.k,
        args.n_boot,
        args.cosine_tau,
        args.match_tau,
        args.seed
    );

    let params = BootstrapParams {
        k: args.k,
        n_boot: args.n_boot,
        seed: args.seed,
        top_n: args.top_n,
        cosine_tau: args.cosine_tau,
        match_tau: args.match_tau,
        max_iter: args.max_iter,
        tol: args.tol,
        min_subjects: args.min_subjects,
    };
    let rows: Vec<BootstrapRow> =
        align_bootstrap(&matrices, params).map_err(|e| anyhow::anyhow!(e))?;

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent).with_context(|| {
            format!("creating output dir {:?}", parent)
        })?;
    }
    write_bootstrap_summary(&args.output, &rows)?;
    eprintln!(
        "align bootstrap: wrote {} archetypes to {}",
        rows.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    // Hash each cohort's canonical input directory separately.
    let mut labeled: Vec<(String, PathBuf)> = Vec::new();
    for (dir, label) in cohort_dirs.iter().zip(labels.iter()) {
        let file = match args.source.as_str() {
            "qc" => "qc_measurements.tsv",
            _ => "measurements.tsv",
        };
        labeled.push((format!("{}_{}", label, file), dir.join(file)));
        labeled.push((format!("{}_samples", label), dir.join("samples.tsv")));
        labeled.push((format!("{}_proteins", label), dir.join("proteins.tsv")));
    }
    let refs: Vec<(&str, &Path)> = labeled
        .iter()
        .map(|(l, p)| (l.as_str(), p.as_path()))
        .collect();
    let input_dir_sha256 = hash_labeled_inputs(&refs)?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "align bootstrap",
        json!({
            "cohorts": args.cohorts,
            "labels": args.labels,
            "k": args.k,
            "n-boot": args.n_boot,
            "seed": args.seed,
            "cosine-tau": args.cosine_tau,
            "match-tau": args.match_tau,
            "top-n": args.top_n,
            "max-iter": args.max_iter,
            "tol": args.tol,
            "min-subjects": args.min_subjects,
            "max-missing-fraction": args.max_missing_fraction,
            "impute": args.impute,
            "source": args.source,
            "output": args.output.display().to_string(),
        }),
        &input_dir_sha256,
        &[args.output.clone()],
        started_at,
        finished_at,
    )?;
    eprintln!("align bootstrap: sidecar={}", sidecar.display());
    Ok(())
}

fn load_cohort_matrix(
    dir: &Path,
    label: &str,
    source: &str,
    max_missing_fraction: f64,
    impute_mean: bool,
) -> Result<CohortMatrix> {
    let file = match source {
        "qc" => "qc_measurements.tsv",
        "raw" => "measurements.tsv",
        other => bail!("--source {other:?}; expected qc or raw"),
    };
    let records = read_measurements_long(&dir.join(file))?;
    if records.is_empty() {
        bail!("no measurements in {:?}", dir.join(file));
    }
    // Abundance by (assay, sample). Only rows with finite values not
    // dropped by QC contribute.
    let mut abundance: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut sample_order: Vec<String> = Vec::new();
    let mut seen_samples: BTreeSet<String> = BTreeSet::new();
    let mut assays: BTreeSet<String> = BTreeSet::new();
    for r in &records {
        if r.dropped_by_qc {
            continue;
        }
        let v = r.abundance.as_f64();
        if !v.is_finite() {
            continue;
        }
        let sid = r.sample_id.clone();
        let aid = r.assay_id.0.clone();
        if seen_samples.insert(sid.clone()) {
            sample_order.push(sid.clone());
        }
        assays.insert(aid.clone());
        abundance.insert((aid, sid), v);
    }
    if sample_order.len() < 2 {
        bail!("cohort {label:?}: need >= 2 samples");
    }
    let n = sample_order.len();
    let mut kept: Vec<String> = Vec::new();
    for a in &assays {
        let present = sample_order
            .iter()
            .filter(|s| abundance.contains_key(&(a.clone(), (*s).clone())))
            .count();
        let missing_frac = 1.0 - present as f64 / n as f64;
        if missing_frac <= max_missing_fraction + 1e-12 {
            kept.push(a.clone());
        }
    }
    if kept.is_empty() {
        bail!(
            "cohort {label:?}: no assays retained at --max-missing-fraction={max_missing_fraction}"
        );
    }

    let mut data = vec![vec![0.0_f64; kept.len()]; n];
    for (j, a) in kept.iter().enumerate() {
        let mut vals: Vec<Option<f64>> = Vec::with_capacity(n);
        let mut sum = 0.0_f64;
        let mut count = 0usize;
        for s in &sample_order {
            let v = abundance.get(&(a.clone(), s.clone())).copied();
            if let Some(x) = v {
                sum += x;
                count += 1;
            }
            vals.push(v);
        }
        let mean = if count > 0 { sum / count as f64 } else { 0.0 };
        for (i, v) in vals.iter().enumerate() {
            match v {
                Some(x) => data[i][j] = *x,
                None => {
                    if impute_mean {
                        data[i][j] = mean;
                    } else {
                        bail!(
                            "cohort {label:?}: assay {} missing sample {}; rerun with --impute mean",
                            a,
                            sample_order[i]
                        );
                    }
                }
            }
        }
    }
    Ok(CohortMatrix {
        label: label.to_string(),
        data,
        protein_labels: kept,
    })
}

/// Restrict each cohort's matrix to the intersection of protein label
/// sets, in a canonical (sorted) order.
fn intersect_cohorts(mut matrices: Vec<CohortMatrix>) -> Result<Vec<CohortMatrix>> {
    if matrices.is_empty() {
        return Ok(matrices);
    }
    let mut universe: BTreeSet<String> =
        matrices[0].protein_labels.iter().cloned().collect();
    for m in matrices.iter().skip(1) {
        let s: BTreeSet<String> = m.protein_labels.iter().cloned().collect();
        universe = universe.intersection(&s).cloned().collect();
    }
    if universe.is_empty() {
        bail!("cohorts share zero proteins after intersection");
    }
    let ordered: Vec<String> = {
        let mut v: Vec<String> = universe.into_iter().collect();
        v.sort();
        v
    };
    for m in &mut matrices {
        let index_by_label: BTreeMap<String, usize> = m
            .protein_labels
            .iter()
            .enumerate()
            .map(|(i, l)| (l.clone(), i))
            .collect();
        let mut new_data = vec![vec![0.0_f64; ordered.len()]; m.data.len()];
        for (j_new, label) in ordered.iter().enumerate() {
            let j_old = index_by_label[label];
            for (i, row) in m.data.iter().enumerate() {
                new_data[i][j_new] = row[j_old];
            }
        }
        m.data = new_data;
        m.protein_labels = ordered.clone();
    }
    Ok(matrices)
}

fn write_bootstrap_summary(path: &Path, rows: &[BootstrapRow]) -> Result<()> {
    let mut out = String::from(
        "archetype_id\tobserved_n_cohorts\tobserved_cohorts\t\
         bootstrap_mean_n_cohorts\tbootstrap_prob_universal\tbootstrap_prob_multi\t\
         ci_lower_n_cohorts\tci_upper_n_cohorts\tbootstrap_match_rate\t\
         alignment_entropy\tbca_lower_n_cohorts\tbca_upper_n_cohorts\t\
         bca_fallback_to_percentile\n",
    );
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t\
             {:.6}\t{:.6}\t{:.6}\t{}\n",
            r.archetype_id,
            r.observed_n_cohorts,
            r.observed_cohorts.join(","),
            r.bootstrap_mean_n_cohorts,
            r.bootstrap_prob_universal,
            r.bootstrap_prob_multi,
            r.ci_lower_n_cohorts,
            r.ci_upper_n_cohorts,
            r.bootstrap_match_rate,
            r.alignment_entropy,
            r.bca_lower_n_cohorts,
            r.bca_upper_n_cohorts,
            if r.bca_fallback_to_percentile { 1 } else { 0 },
        ));
    }
    atomic_write(path, out.as_bytes())
}

fn run_programs(args: ProgramsArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let paths: Vec<PathBuf> = args
        .loadings
        .split(',')
        .map(|s| PathBuf::from(s.trim()))
        .collect();
    if paths.is_empty() {
        bail!("--loadings must list at least one TSV");
    }
    let cohort_labels: Vec<String> = match &args.cohorts {
        Some(c) => c.split(',').map(|s| s.trim().to_string()).collect(),
        None => paths
            .iter()
            .map(|p| file_stem(p).to_string())
            .collect(),
    };
    if cohort_labels.len() != paths.len() {
        bail!(
            "--cohorts has {} labels but --loadings has {} paths",
            cohort_labels.len(),
            paths.len()
        );
    }
    let annotation_paths: Vec<Option<PathBuf>> = match &args.annotations {
        Some(a) => {
            let parts: Vec<PathBuf> = a.split(',').map(|s| PathBuf::from(s.trim())).collect();
            if parts.len() != paths.len() {
                bail!(
                    "--annotations has {} paths but --loadings has {} paths",
                    parts.len(),
                    paths.len()
                );
            }
            parts.into_iter().map(Some).collect()
        }
        None => vec![None; paths.len()],
    };

    // Pass 1: build union of labels across all cohorts.
    let mut loadings_raw: Vec<Vec<(String, String, f64)>> = Vec::with_capacity(paths.len());
    for path in &paths {
        loadings_raw.push(read_loadings(path, &args.label_col)?);
    }
    let mut label_index: BTreeMap<String, usize> = BTreeMap::new();
    for cohort in &loadings_raw {
        for (_prog, label, _val) in cohort {
            if !label_index.contains_key(label) {
                let idx = label_index.len();
                label_index.insert(label.clone(), idx);
            }
        }
    }
    let universe_size = label_index.len();
    if universe_size == 0 {
        bail!("no labels found across cohorts");
    }

    // Assemble AlignedProgram structs across cohorts.
    let mut programs: Vec<AlignedProgram> = Vec::new();
    for (ci, cohort_rows) in loadings_raw.iter().enumerate() {
        let cohort_label = &cohort_labels[ci];
        let annotation_map: BTreeMap<String, String> = match &annotation_paths[ci] {
            Some(ap) => read_annotations(ap)?,
            None => BTreeMap::new(),
        };
        let mut by_program: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        for (program, label, value) in cohort_rows {
            let vec = by_program
                .entry(program.clone())
                .or_insert_with(|| vec![0.0_f64; universe_size]);
            if let Some(idx) = label_index.get(label) {
                vec[*idx] = *value;
            }
        }
        for (program, values) in by_program {
            programs.push(AlignedProgram {
                cohort: cohort_label.clone(),
                program: program.clone(),
                category: annotation_map.get(&program).cloned(),
                values,
            });
        }
    }

    let primary_output = if args.sweep {
        let matrix_path = args
            .output_matrix
            .as_ref()
            .context("sweep mode requires --output-matrix")?;
        let metrics_spec = args
            .metrics
            .as_deref()
            .context("sweep mode requires --metrics")?;
        let grid = parse_sweep_spec(metrics_spec)?;
        run_sweep(
            &programs,
            &grid,
            args.reciprocal_best,
            args.compare_constrained_vs_unconstrained,
            matrix_path,
        )?;
        matrix_path.clone()
    } else {
        let output = args
            .output
            .as_ref()
            .context("single-mode requires --output")?;
        let labels = build_archetypes(
            &programs,
            args.metric.into(),
            args.top_n,
            args.tau,
            args.reciprocal_best,
            args.category_constraint,
        );
        write_archetypes(output, &programs, &labels)?;
        let summary = summarize_archetypes(&programs, &labels);
        eprintln!(
            "align programs: archetypes={} multi_cohort={} universal={} category_recovery={:.3}",
            summary.n_archetypes,
            summary.n_multi_cohort,
            summary.n_universal,
            summary.category_recovery
        );
        output.clone()
    };

    let finished_at = SystemTime::now();
    let mut labeled: Vec<(String, PathBuf)> = Vec::new();
    for (ci, path) in paths.iter().enumerate() {
        labeled.push((format!("loadings_{}", cohort_labels[ci]), path.clone()));
    }
    for (ci, ap) in annotation_paths.iter().enumerate() {
        if let Some(p) = ap {
            labeled.push((format!("annotations_{}", cohort_labels[ci]), p.clone()));
        }
    }
    let refs: Vec<(&str, &Path)> = labeled
        .iter()
        .map(|(l, p)| (l.as_str(), p.as_path()))
        .collect();
    let input_dir_sha256 = hash_labeled_inputs(&refs)?;
    let sidecar = sidecar_path_for(&primary_output);
    let metric = match args.metric {
        SingleMetric::Jaccard => "jaccard",
        SingleMetric::Cosine => "cosine",
        SingleMetric::Spearman => "spearman",
    };
    write_run_sidecar(
        &sidecar,
        "align programs",
        json!({
            "loadings": args.loadings,
            "cohorts": args.cohorts,
            "annotations": args.annotations,
            "label-col": args.label_col,
            "metric": metric,
            "top-n": args.top_n,
            "tau": args.tau,
            "reciprocal-best": args.reciprocal_best,
            "category-constraint": args.category_constraint,
            "output": args.output.as_ref().map(|p| p.display().to_string()),
            "sweep": args.sweep,
            "metrics": args.metrics,
            "compare-constrained-vs-unconstrained": args.compare_constrained_vs_unconstrained,
            "output-matrix": args.output_matrix.as_ref().map(|p| p.display().to_string()),
        }),
        &input_dir_sha256,
        &[primary_output.clone()],
        started_at,
        finished_at,
    )?;
    eprintln!("align programs: sidecar={}", sidecar.display());
    Ok(())
}

fn read_loadings(path: &Path, label_col: &str) -> Result<Vec<(String, String, f64)>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let program_col = headers
        .iter()
        .position(|h| h == "program")
        .with_context(|| format!("{:?} missing 'program' column", path))?;
    let loading_col = headers
        .iter()
        .position(|h| h == "loading")
        .with_context(|| format!("{:?} missing 'loading' column", path))?;
    let label_idx = headers
        .iter()
        .position(|h| h == label_col)
        .or_else(|| headers.iter().position(|h| h == "gene_symbol"))
        .or_else(|| headers.iter().position(|h| h == "assay_id"))
        .or_else(|| headers.iter().position(|h| h == "protein"))
        .with_context(|| format!("{:?} missing label column", path))?;
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        let value: f64 = row[loading_col]
            .parse()
            .with_context(|| format!("parsing loading in {:?}", path))?;
        if !value.is_finite() {
            continue;
        }
        let label = row[label_idx].trim();
        if label.is_empty() {
            continue;
        }
        out.push((
            row[program_col].to_string(),
            label.to_string(),
            value,
        ));
    }
    if out.is_empty() {
        bail!("no usable loadings in {:?}", path);
    }
    Ok(out)
}

fn read_annotations(path: &Path) -> Result<BTreeMap<String, String>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let program_col = headers
        .iter()
        .position(|h| h == "program")
        .with_context(|| format!("{:?} missing 'program' column", path))?;
    let category_col = headers
        .iter()
        .position(|h| h == "category")
        .with_context(|| format!("{:?} missing 'category' column", path))?;
    let mut out = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let program = row[program_col].to_string();
        let category = row[category_col].trim().to_string();
        if !category.is_empty() {
            out.insert(program, category);
        }
    }
    Ok(out)
}

fn file_stem(path: &Path) -> &str {
    path.file_stem()
        .and_then(|os| os.to_str())
        .unwrap_or("cohort")
}

fn write_archetypes(path: &Path, programs: &[AlignedProgram], labels: &[usize]) -> Result<()> {
    // Assign compact archetype ids: sort members by (cohort, program), group by label.
    let mut by_label: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, lab) in labels.iter().enumerate() {
        by_label.entry(*lab).or_default().push(i);
    }
    let mut archetype_id_map: BTreeMap<usize, usize> = BTreeMap::new();
    let mut next_id = 1usize;
    // Assign multi-member archetypes first in deterministic order.
    for (lab, members) in by_label.iter() {
        if members.len() >= 2 {
            archetype_id_map.insert(*lab, next_id);
            next_id += 1;
        }
    }
    for (lab, members) in by_label.iter() {
        if members.len() < 2 && !archetype_id_map.contains_key(lab) {
            archetype_id_map.insert(*lab, next_id);
            next_id += 1;
        }
    }
    let mut out = String::from(
        "archetype_id\tcohort\tprogram\tcategory\tn_members\tn_cohorts\tis_singleton\n",
    );
    let mut archetype_sizes: BTreeMap<usize, (usize, BTreeSet<String>)> = BTreeMap::new();
    for (lab, members) in by_label.iter() {
        let cohorts: BTreeSet<String> = members
            .iter()
            .map(|i| programs[*i].cohort.clone())
            .collect();
        archetype_sizes.insert(*lab, (members.len(), cohorts));
    }
    let mut rows: Vec<(usize, &AlignedProgram, usize, usize, bool)> = Vec::new();
    for (lab, members) in by_label.iter() {
        let id = archetype_id_map[lab];
        let (n_members, cohorts) = &archetype_sizes[lab];
        for &i in members {
            rows.push((id, &programs[i], *n_members, cohorts.len(), *n_members < 2));
        }
    }
    rows.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.cohort.cmp(&b.1.cohort))
            .then(a.1.program.cmp(&b.1.program))
    });
    for (id, prog, n_members, n_cohorts, singleton) in rows {
        out.push_str(&format!(
            "A{:04}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            id,
            prog.cohort,
            prog.program,
            prog.category.clone().unwrap_or_default(),
            n_members,
            n_cohorts,
            if singleton { 1 } else { 0 }
        ));
    }
    atomic_write(path, out.as_bytes())
}

#[derive(Debug, Clone)]
struct SweepCell {
    metric: AlignMetric,
    top_n: usize,
    tau: f64,
}

fn split_top_level(input: &str, separator: char) -> Vec<String> {
    let mut out = Vec::new();
    let mut depth = 0i32;
    let mut current = String::new();
    for c in input.chars() {
        match c {
            '[' | '(' | '{' => {
                depth += 1;
                current.push(c);
            }
            ']' | ')' | '}' => {
                depth -= 1;
                current.push(c);
            }
            other if other == separator && depth == 0 => {
                out.push(std::mem::take(&mut current));
            }
            other => current.push(other),
        }
    }
    if !current.is_empty() {
        out.push(current);
    }
    out
}

fn parse_sweep_spec(spec: &str) -> Result<Vec<SweepCell>> {
    let mut out = Vec::new();
    for block in split_top_level(spec, ',') {
        let block = block.trim();
        if block.is_empty() {
            continue;
        }
        let parts = split_top_level(&block, ':');
        let metric_name = parts.first().cloned().unwrap_or_default();
        let metric = AlignMetric::parse(&metric_name)
            .with_context(|| format!("unknown metric {:?} in --metrics", metric_name))?;
        let mut top_ns = vec![40usize];
        let mut taus: Vec<f64> = Vec::new();
        for kv in parts.iter().skip(1) {
            let kv = kv.trim();
            if kv.is_empty() {
                continue;
            }
            let (key, raw_values) = kv
                .split_once('=')
                .with_context(|| format!("expected key=value in {:?}", kv))?;
            let inner = raw_values
                .trim()
                .trim_start_matches('[')
                .trim_end_matches(']');
            match key.trim() {
                "top_n" => {
                    top_ns = inner
                        .split(',')
                        .map(|s| s.trim().parse())
                        .collect::<std::result::Result<Vec<usize>, _>>()
                        .with_context(|| format!("parsing top_n list in {:?}", kv))?;
                }
                "tau" => {
                    taus = inner
                        .split(',')
                        .map(|s| s.trim().parse())
                        .collect::<std::result::Result<Vec<f64>, _>>()
                        .with_context(|| format!("parsing tau list in {:?}", kv))?;
                }
                other => bail!("unknown sweep key {:?}", other),
            }
        }
        if taus.is_empty() {
            bail!("metric {:?} has no tau values", metric_name);
        }
        for &tn in &top_ns {
            for &tau in &taus {
                out.push(SweepCell {
                    metric,
                    top_n: tn,
                    tau,
                });
            }
        }
    }
    if out.is_empty() {
        bail!("--metrics produced no sweep cells");
    }
    Ok(out)
}

fn run_sweep(
    programs: &[AlignedProgram],
    grid: &[SweepCell],
    reciprocal_best: bool,
    compare_constrained: bool,
    output: &Path,
) -> Result<()> {
    let mut out = String::from(
        "metric\ttop_n\ttau\treciprocal_best\tcategory_constraint\tn_archetypes\tn_multi_cohort\tn_universal\tcategory_recovery\n",
    );
    let modes: Vec<(bool, &str)> = if compare_constrained {
        vec![(true, "1"), (false, "0")]
    } else {
        vec![(false, "0")]
    };
    for cell in grid {
        for (cat_flag, cat_flag_str) in &modes {
            let labels = build_archetypes(
                programs,
                cell.metric,
                cell.top_n,
                cell.tau,
                reciprocal_best,
                *cat_flag,
            );
            let summary = summarize_archetypes(programs, &labels);
            let cat_recovery = if summary.category_recovery.is_finite() {
                format!("{:.6}", summary.category_recovery)
            } else {
                "NaN".to_string()
            };
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                cell.metric.as_str(),
                cell.top_n,
                cell.tau,
                if reciprocal_best { 1 } else { 0 },
                cat_flag_str,
                summary.n_archetypes,
                summary.n_multi_cohort,
                summary.n_universal,
                cat_recovery,
            ));
        }
    }
    atomic_write(output, out.as_bytes())?;
    eprintln!(
        "align programs sweep: cells={} output={}",
        grid.len() * modes.len(),
        output.display()
    );
    Ok(())
}
