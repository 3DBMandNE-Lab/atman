//! `atman network influence` — feature-covariance hub scoring.
//!
//! Builds a feature × feature similarity graph over subjects and
//! scores each feature by its role as a hub via eigenvector
//! centrality × betweenness centrality. Optional stratification
//! (`--stratify <column>`) builds one graph per stratum.

use anyhow::{bail, Context, Result};
use atman_core::network::{
    adjacency, influence_scores, pairwise_similarity, AdjacencyPolicy, InfluenceRow,
    SimilarityMetric,
};
use atman_core::network_differential::{
    edge_pairwise_differential, edge_summary_differential, module_rewiring,
    signed_pairwise_correlations, CohortCorrelations, CohortData, SignedMetric,
};
use atman_core::MeasurementRecord;
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use csv::ReaderBuilder;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, format_float, hash_labeled_inputs, read_measurements_long, read_samples,
    sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Feature-covariance hub scoring via eigenvector ×
    /// betweenness centrality.
    Influence(InfluenceArgs),
    /// Cross-cohort differential coexpression. Three output modes:
    /// edge-pairwise (DGCA-style Fisher-z test), edge-summary
    /// (per-edge cross-cohort summary), module (per-module within-
    /// module connectivity rewiring).
    Differential(DifferentialArgs),
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum DifferentialMode {
    /// Per-edge × per-cohort-pair Fisher-z test of `H_0: corr_A = corr_B`.
    EdgePairwise,
    /// Per-edge cross-cohort summary (mean, sd, sign-flips, conservation,
    /// divergence). One row per edge.
    EdgeSummary,
    /// Per-module within-module connectivity per cohort, with
    /// cross-cohort rewiring score. Requires --gene-sets.
    Module,
}

#[derive(ClapArgs, Debug)]
pub struct DifferentialArgs {
    /// Per-cohort canonical input directories. Pass as
    /// `LABEL=path,LABEL=path,...` to override the default cohort label
    /// (which is the directory basename). Each directory must contain
    /// `measurements.tsv` and `samples.tsv`.
    #[arg(long, value_delimiter = ',')]
    inputs: Vec<String>,

    /// Output TSV path.
    #[arg(long)]
    output: PathBuf,

    /// Differential mode (edge-pairwise / edge-summary / module).
    #[arg(long, value_enum)]
    mode: DifferentialMode,

    /// Gene-set TSV (`set_name`, `gene_symbol`). Required when
    /// `--mode module`.
    #[arg(long)]
    gene_sets: Option<PathBuf>,

    /// Feature column to key on (`gene_symbol` or `assay_id`).
    #[arg(long, default_value = "gene_symbol")]
    feature_col: String,

    /// Signed correlation metric.
    #[arg(long, value_enum, default_value_t = MetricArg::Pearson)]
    method: MetricArg,

    /// Minimum per-pair sample overlap within each cohort. Edges with
    /// fewer overlapping subjects in any cohort are dropped from the
    /// output.
    #[arg(long, default_value_t = 5)]
    min_overlap: usize,

    /// Cap on emitted rows for `edge-pairwise` and `edge-summary` modes
    /// (sorted by `|z_diff|` and `divergence_score` respectively, then
    /// by lexicographic edge key as the deterministic tie-break). 0
    /// means no cap.
    #[arg(long, default_value_t = 0)]
    top_rows: usize,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum MetricArg {
    #[value(name = "pearson")]
    Pearson,
    #[value(name = "spearman")]
    Spearman,
    #[value(name = "covariance")]
    Covariance,
}

impl MetricArg {
    fn to_core(self) -> SimilarityMetric {
        match self {
            Self::Pearson => SimilarityMetric::Pearson,
            Self::Spearman => SimilarityMetric::Spearman,
            Self::Covariance => SimilarityMetric::Covariance,
        }
    }
}

#[derive(ClapArgs, Debug)]
pub struct InfluenceArgs {
    /// Canonical long-format measurements TSV (usually
    /// `measurements.tsv`) with `sample_id, assay_id, gene_symbol,
    /// abundance, dropped_by_qc` columns.
    #[arg(long)]
    input: PathBuf,

    /// `samples.tsv` for stratification metadata.
    #[arg(long)]
    samples: PathBuf,

    /// Feature column to key nodes by. `gene_symbol` (default) or
    /// `assay_id`; anything else errors.
    #[arg(long, default_value = "gene_symbol")]
    feature_col: String,

    /// Similarity metric.
    #[arg(long, value_enum, default_value_t = MetricArg::Spearman)]
    method: MetricArg,

    /// Hard-threshold on `|similarity|` — edges below this get
    /// weight 0. Mutually exclusive with `--soft-power`.
    #[arg(long)]
    threshold: Option<f64>,

    /// Soft-power β in the WGCNA `|r|^β` construction. Mutually
    /// exclusive with `--threshold`.
    #[arg(long)]
    soft_power: Option<u32>,

    /// Optional column in `samples.tsv` — compute one graph per
    /// stratum. `stratum` column in the output carries the stratum
    /// name; `""` when stratification is not set.
    #[arg(long)]
    stratify: Option<String>,

    /// Output TSV.
    #[arg(long)]
    output: PathBuf,

    /// When set, also write `<output_dir>/network_edges.tsv` with
    /// every retained edge. Can be large on dense graphs.
    #[arg(long, default_value_t = false)]
    emit_edges: bool,

    /// Max power-iteration steps for eigenvector centrality.
    #[arg(long, default_value_t = 500)]
    eigen_max_iter: usize,

    /// Convergence tolerance (max elementwise change) for the
    /// eigenvector iteration.
    #[arg(long, default_value_t = 1e-8)]
    eigen_tol: f64,

    /// Minimum subject count per stratum required for a graph to be
    /// built. Below this, the stratum is skipped with a stderr note.
    #[arg(long, default_value_t = 5)]
    min_subjects: usize,
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Influence(args) => run_influence(args),
        Command::Differential(args) => run_differential(args),
    }
}

fn metric_arg_to_signed(metric: MetricArg) -> Result<SignedMetric> {
    match metric {
        MetricArg::Pearson => Ok(SignedMetric::Pearson),
        MetricArg::Spearman => Ok(SignedMetric::Spearman),
        MetricArg::Covariance => bail!(
            "network differential requires a signed correlation metric; pass --method pearson or spearman"
        ),
    }
}

fn parse_cohort_inputs(inputs: &[String]) -> Result<Vec<(String, PathBuf)>> {
    if inputs.len() < 2 {
        bail!("network differential requires at least 2 cohorts via --inputs");
    }
    let mut out = Vec::with_capacity(inputs.len());
    let mut seen_labels: BTreeSet<String> = BTreeSet::new();
    for raw in inputs {
        let token = raw.trim();
        let (label, path) = if let Some((lhs, rhs)) = token.split_once('=') {
            (lhs.trim().to_string(), PathBuf::from(rhs.trim()))
        } else {
            let path = PathBuf::from(token);
            let label = path
                .file_name()
                .and_then(|s| s.to_str())
                .map(|s| s.to_string())
                .unwrap_or_else(|| token.to_string());
            (label, path)
        };
        if label.is_empty() {
            bail!("empty cohort label in --inputs entry {token:?}");
        }
        if !seen_labels.insert(label.clone()) {
            bail!("duplicate cohort label {label:?} in --inputs");
        }
        if !path.is_dir() {
            bail!("cohort directory does not exist: {path:?}");
        }
        out.push((label, path));
    }
    Ok(out)
}

struct LoadedCohort {
    label: String,
    feature_labels: Vec<String>,
    /// `feature_labels[i]` → values across the kept subjects (NaN if missing
    /// in that subject). Length per row equals the number of subjects in
    /// the cohort.
    data: Vec<Vec<f64>>,
}

fn load_cohort(
    label: &str,
    dir: &Path,
    shared_features: &[String],
    feature_col: &str,
) -> Result<LoadedCohort> {
    let measurements_path = dir.join("measurements.tsv");
    let samples_path = dir.join("samples.tsv");
    let records = read_measurements_long(&measurements_path)
        .with_context(|| format!("reading {measurements_path:?}"))?;
    let _samples = read_samples(&samples_path)
        .with_context(|| format!("reading {samples_path:?}"))?;
    let FeatureMatrix {
        features: _,
        per_sample,
    } = collect_feature_matrix(&records, feature_col)?;
    if per_sample.is_empty() {
        bail!("no subjects with measurements in {measurements_path:?}");
    }
    // Build the rectangular feature × subject matrix on the shared
    // feature universe. Subjects appear in deterministic sorted order.
    let subject_ids: Vec<String> = per_sample.keys().cloned().collect();
    let mut data = Vec::with_capacity(shared_features.len());
    for feat in shared_features {
        let row: Vec<f64> = subject_ids
            .iter()
            .map(|sid| {
                per_sample
                    .get(sid)
                    .and_then(|m| m.get(feat))
                    .copied()
                    .unwrap_or(f64::NAN)
            })
            .collect();
        data.push(row);
    }
    Ok(LoadedCohort {
        label: label.to_string(),
        feature_labels: shared_features.to_vec(),
        data,
    })
}

fn shared_feature_universe(cohort_dirs: &[(String, PathBuf)], feature_col: &str) -> Result<Vec<String>> {
    let mut intersection: Option<BTreeSet<String>> = None;
    for (label, dir) in cohort_dirs {
        let measurements_path = dir.join("measurements.tsv");
        let records = read_measurements_long(&measurements_path)
            .with_context(|| format!("reading {measurements_path:?}"))?;
        let FeatureMatrix { features, .. } = collect_feature_matrix(&records, feature_col)?;
        let here: BTreeSet<String> = features.into_iter().collect();
        intersection = Some(match intersection {
            None => here,
            Some(prev) => prev.intersection(&here).cloned().collect(),
        });
        eprintln!(
            "network differential: {} contributed {} features",
            label,
            intersection.as_ref().unwrap().len()
        );
    }
    let shared = intersection.unwrap_or_default();
    if shared.len() < 2 {
        bail!(
            "shared feature universe has fewer than 2 features ({})",
            shared.len()
        );
    }
    Ok(shared.into_iter().collect())
}

fn read_gene_sets_for_modules(path: &Path) -> Result<BTreeMap<String, BTreeSet<String>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {path:?}"))?;
    let headers = reader.headers()?.clone();
    let set_col = headers
        .iter()
        .position(|h| h == "set_name")
        .with_context(|| format!("missing column `set_name` in {path:?}"))?;
    let gene_col = headers
        .iter()
        .position(|h| h == "gene_symbol")
        .with_context(|| format!("missing column `gene_symbol` in {path:?}"))?;
    let mut out: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let s = row[set_col].trim();
        let g = row[gene_col].trim();
        if s.is_empty() || g.is_empty() {
            continue;
        }
        out.entry(s.to_string()).or_default().insert(g.to_string());
    }
    Ok(out)
}

fn run_differential(args: DifferentialArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let metric = metric_arg_to_signed(args.method)?;
    if args.feature_col != "gene_symbol" && args.feature_col != "assay_id" {
        bail!(
            "--feature-col must be `gene_symbol` or `assay_id`; got {:?}",
            args.feature_col
        );
    }
    if matches!(args.mode, DifferentialMode::Module) && args.gene_sets.is_none() {
        bail!("--mode module requires --gene-sets");
    }
    let cohort_dirs = parse_cohort_inputs(&args.inputs)?;
    let shared_features = shared_feature_universe(&cohort_dirs, &args.feature_col)?;
    eprintln!(
        "network differential: shared feature universe = {} features across {} cohorts",
        shared_features.len(),
        cohort_dirs.len()
    );

    let mut loaded: Vec<LoadedCohort> = Vec::with_capacity(cohort_dirs.len());
    for (label, dir) in &cohort_dirs {
        loaded.push(load_cohort(label, dir, &shared_features, &args.feature_col)?);
    }
    let mut correlations: Vec<CohortCorrelations> = Vec::with_capacity(loaded.len());
    for c in &loaded {
        let cd = CohortData {
            label: c.label.clone(),
            feature_labels: &c.feature_labels,
            data: &c.data,
        };
        let cc = signed_pairwise_correlations(&cd, metric, args.min_overlap)
            .with_context(|| format!("computing correlations for cohort {}", c.label))?;
        correlations.push(cc);
    }
    let cohort_refs: Vec<&CohortCorrelations> = correlations.iter().collect();

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output dir {parent:?}"))?;
    }

    let n_rows = match args.mode {
        DifferentialMode::EdgePairwise => {
            let mut rows = edge_pairwise_differential(&cohort_refs);
            rows.sort_by(|a, b| {
                b.z_diff
                    .abs()
                    .partial_cmp(&a.z_diff.abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.feature_a.cmp(&b.feature_a))
                    .then_with(|| a.feature_b.cmp(&b.feature_b))
                    .then_with(|| a.cohort_a.cmp(&b.cohort_a))
                    .then_with(|| a.cohort_b.cmp(&b.cohort_b))
            });
            if args.top_rows > 0 && rows.len() > args.top_rows {
                rows.truncate(args.top_rows);
            }
            write_edge_pairwise(&args.output, &rows)?;
            rows.len()
        }
        DifferentialMode::EdgeSummary => {
            let mut rows = edge_summary_differential(&cohort_refs, args.min_overlap);
            rows.sort_by(|a, b| {
                b.divergence_score
                    .partial_cmp(&a.divergence_score)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.feature_a.cmp(&b.feature_a))
                    .then_with(|| a.feature_b.cmp(&b.feature_b))
            });
            if args.top_rows > 0 && rows.len() > args.top_rows {
                rows.truncate(args.top_rows);
            }
            write_edge_summary(&args.output, &rows, &cohort_refs)?;
            rows.len()
        }
        DifferentialMode::Module => {
            let gene_sets_path = args.gene_sets.as_ref().expect("guarded above");
            let modules = read_gene_sets_for_modules(gene_sets_path)?;
            if modules.is_empty() {
                bail!("no gene sets found in {gene_sets_path:?}");
            }
            let mut rows = module_rewiring(&cohort_refs, &modules);
            rows.sort_by(|a, b| {
                b.rewiring_score
                    .partial_cmp(&a.rewiring_score)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.set_name.cmp(&b.set_name))
            });
            write_module_rewiring(&args.output, &rows, &cohort_refs)?;
            rows.len()
        }
    };
    eprintln!(
        "network differential: mode={:?} rows={n_rows}",
        args.mode
    );

    let finished_at = SystemTime::now();
    let mut labeled: Vec<(&str, &Path)> = Vec::with_capacity(2 * cohort_dirs.len() + 1);
    let measurements_paths: Vec<PathBuf> = cohort_dirs
        .iter()
        .map(|(_, d)| d.join("measurements.tsv"))
        .collect();
    let samples_paths: Vec<PathBuf> = cohort_dirs
        .iter()
        .map(|(_, d)| d.join("samples.tsv"))
        .collect();
    for (i, (label, _)) in cohort_dirs.iter().enumerate() {
        labeled.push((
            Box::leak(format!("{label}_measurements").into_boxed_str()),
            measurements_paths[i].as_path(),
        ));
        labeled.push((
            Box::leak(format!("{label}_samples").into_boxed_str()),
            samples_paths[i].as_path(),
        ));
    }
    let gene_sets_path_storage: PathBuf;
    if let Some(p) = args.gene_sets.as_ref() {
        gene_sets_path_storage = p.clone();
        labeled.push(("gene_sets", gene_sets_path_storage.as_path()));
    }
    let inputs_sha256 = hash_labeled_inputs(&labeled)?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "network differential",
        json!({
            "inputs": cohort_dirs
                .iter()
                .map(|(l, p)| format!("{l}={}", p.display()))
                .collect::<Vec<_>>(),
            "output": args.output.display().to_string(),
            "mode": format!("{:?}", args.mode).to_ascii_lowercase().replace("::", "-"),
            "feature-col": args.feature_col,
            "method": match metric {
                SignedMetric::Pearson => "pearson",
                SignedMetric::Spearman => "spearman",
            },
            "min-overlap": args.min_overlap,
            "top-rows": args.top_rows,
            "gene-sets": args.gene_sets.as_ref().map(|p| p.display().to_string()),
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("network differential: sidecar={}", sidecar.display());
    Ok(())
}

fn write_edge_pairwise(
    path: &Path,
    rows: &[atman_core::network_differential::EdgePairwiseRow],
) -> Result<()> {
    let mut buf = String::from(
        "feature_a\tfeature_b\tcohort_a\tcohort_b\tn_a\tn_b\tcorr_a\tcorr_b\tz_diff\tp_value\n",
    );
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.feature_a,
            r.feature_b,
            r.cohort_a,
            r.cohort_b,
            r.n_a,
            r.n_b,
            format_float(r.corr_a),
            format_float(r.corr_b),
            format_float(r.z_diff),
            format_float(r.p_value),
        ));
    }
    atomic_write(path, buf.as_bytes())
}

fn write_edge_summary(
    path: &Path,
    rows: &[atman_core::network_differential::EdgeSummaryRow],
    cohorts: &[&CohortCorrelations],
) -> Result<()> {
    let mut header = String::from(
        "feature_a\tfeature_b\tn_cohorts\tmean_corr\tsd_corr\tmin_abs_corr\tmax_abs_corr\trange_corr\tn_sign_flips\tconservation_score\tdivergence_score",
    );
    for c in cohorts {
        header.push('\t');
        header.push_str(&format!("{}_corr", c.label));
    }
    header.push('\n');
    let mut buf = header;
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.feature_a,
            r.feature_b,
            r.n_cohorts,
            format_float(r.mean_corr),
            format_float(r.sd_corr),
            format_float(r.min_abs_corr),
            format_float(r.max_abs_corr),
            format_float(r.range_corr),
            r.n_sign_flips,
            format_float(r.conservation_score),
            format_float(r.divergence_score),
        ));
        for (_, c) in &r.per_cohort_corrs {
            buf.push('\t');
            buf.push_str(&format_float(*c));
        }
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn write_module_rewiring(
    path: &Path,
    rows: &[atman_core::network_differential::ModuleRewiringRow],
    cohorts: &[&CohortCorrelations],
) -> Result<()> {
    let mut header = String::from(
        "set_name\tset_size_declared\tset_size_observed\tmean_connectivity\tsd_connectivity\trange_connectivity\trewiring_score",
    );
    for c in cohorts {
        header.push('\t');
        header.push_str(&format!("{}_connectivity", c.label));
    }
    header.push('\n');
    let mut buf = header;
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.set_name,
            r.set_size_declared,
            r.set_size_observed,
            format_float(r.mean_connectivity),
            format_float(r.sd_connectivity),
            format_float(r.range_connectivity),
            format_float(r.rewiring_score),
        ));
        for (_, c) in &r.per_cohort_connectivity {
            buf.push('\t');
            buf.push_str(&format_float(*c));
        }
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn run_influence(args: InfluenceArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if args.threshold.is_some() && args.soft_power.is_some() {
        bail!("--threshold and --soft-power are mutually exclusive");
    }
    let policy = match (args.threshold, args.soft_power) {
        (Some(t), None) => {
            if !(0.0..=1.0).contains(&t) {
                bail!("--threshold must be in [0, 1]");
            }
            AdjacencyPolicy::HardThreshold { threshold: t }
        }
        (None, Some(p)) => {
            if p == 0 {
                bail!("--soft-power must be >= 1");
            }
            AdjacencyPolicy::SoftPower { power: p }
        }
        (None, None) => AdjacencyPolicy::HardThreshold { threshold: 0.30 },
        (Some(_), Some(_)) => unreachable!(),
    };
    if args.feature_col != "gene_symbol" && args.feature_col != "assay_id" {
        bail!(
            "--feature-col must be `gene_symbol` or `assay_id`; got {:?}",
            args.feature_col
        );
    }

    let records = read_measurements_long(&args.input)?;
    if records.is_empty() {
        bail!("no measurements in {:?}", args.input);
    }
    let samples = read_samples(&args.samples)?;
    if samples.is_empty() {
        bail!("samples.tsv at {:?} is empty", args.samples);
    }

    // Optional stratum map; when stratifying, look up the column
    // value in the samples.tsv extras (since `read_samples` returns
    // a fixed-shape Sample struct).
    let stratum_by_sample: BTreeMap<String, String> = match args.stratify.as_deref() {
        Some(col) => read_stratum(&args.samples, col)?,
        None => BTreeMap::new(),
    };

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output dir {:?}", parent))?;
    }

    // Build feature universe + per-sample per-feature values.
    let FeatureMatrix {
        features,
        per_sample,
    } = collect_feature_matrix(&records, &args.feature_col)?;
    if features.is_empty() {
        bail!("no features extracted from measurements");
    }

    // Partition sample IDs by stratum (or a single "" key).
    let mut samples_by_stratum: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let sample_universe: BTreeSet<&String> = per_sample.keys().collect();
    for sid in &sample_universe {
        let stratum = if args.stratify.is_some() {
            stratum_by_sample.get(*sid).cloned().unwrap_or_default()
        } else {
            String::new()
        };
        samples_by_stratum
            .entry(stratum)
            .or_default()
            .push((*sid).clone());
    }

    let mut all_rows: Vec<(String, InfluenceRow, usize)> = Vec::new(); // (stratum, row, n_subjects)
    let mut all_edges: Vec<EdgeRow> = Vec::new();
    let mut strata_summary: Vec<(String, usize, usize)> = Vec::new(); // (stratum, n_subjects, n_features)
    let mut audit_total = atman_core::network::SimilarityAudit::default();
    let mut per_stratum_audit: Vec<(String, atman_core::network::SimilarityAudit)> = Vec::new();

    for (stratum, sample_ids) in &samples_by_stratum {
        if sample_ids.len() < args.min_subjects {
            eprintln!(
                "network influence: stratum {:?} skipped ({} subjects < --min-subjects={})",
                stratum,
                sample_ids.len(),
                args.min_subjects
            );
            continue;
        }
        // Build feature × subject matrix (complete-case rows).
        let mut data: Vec<Vec<f64>> = Vec::with_capacity(features.len());
        for feature in &features {
            let mut row = Vec::with_capacity(sample_ids.len());
            for sid in sample_ids {
                let v = per_sample
                    .get(sid)
                    .and_then(|m| m.get(feature))
                    .copied()
                    .unwrap_or(f64::NAN);
                row.push(v);
            }
            data.push(row);
        }
        let Some((similarity, audit)) =
            pairwise_similarity(&features, &data, args.method.to_core())
        else {
            eprintln!(
                "network influence: stratum {:?} produced no similarity matrix (insufficient data)",
                stratum
            );
            continue;
        };
        if audit.n_pairs_insufficient_overlap > 0 || audit.n_pairs_undefined_metric > 0 {
            eprintln!(
                "network influence: stratum {:?} similarity audit: {}/{} pairs fell back to 0.0 \
                 ({} insufficient overlap, {} undefined metric)",
                stratum,
                audit.n_pairs_insufficient_overlap + audit.n_pairs_undefined_metric,
                audit.n_pairs_total,
                audit.n_pairs_insufficient_overlap,
                audit.n_pairs_undefined_metric,
            );
        }
        audit_total.n_pairs_total += audit.n_pairs_total;
        audit_total.n_pairs_insufficient_overlap += audit.n_pairs_insufficient_overlap;
        audit_total.n_pairs_undefined_metric += audit.n_pairs_undefined_metric;
        per_stratum_audit.push((stratum.clone(), audit));
        let adj = adjacency(&similarity, policy);
        let rows = influence_scores(&features, &adj, args.eigen_max_iter, args.eigen_tol);
        let n_nontrivial = rows.iter().filter(|r| r.degree > 0).count();
        if n_nontrivial == 0 {
            bail!(
                "network influence: stratum {:?} yielded an empty graph under policy {:?}; \
                 lower --threshold or switch to --soft-power",
                stratum,
                policy
            );
        }
        let n_subjects = sample_ids.len();
        strata_summary.push((stratum.clone(), n_subjects, features.len()));
        for row in rows {
            all_rows.push((stratum.clone(), row, n_subjects));
        }
        if args.emit_edges {
            for i in 0..features.len() {
                for j in (i + 1)..features.len() {
                    let w = adj[i][j];
                    if w > 0.0 {
                        all_edges.push(EdgeRow {
                            feature_i: features[i].clone(),
                            feature_j: features[j].clone(),
                            edge_weight: w,
                            stratum: stratum.clone(),
                        });
                    }
                }
            }
        }
    }

    if all_rows.is_empty() {
        bail!(
            "network influence: no strata produced a graph at --min-subjects={}",
            args.min_subjects
        );
    }

    write_influence_tsv(&args.output, &all_rows)?;
    let mut outputs = vec![args.output.clone()];
    let edges_path = args.output.with_file_name("network_edges.tsv");
    if args.emit_edges {
        write_edges_tsv(&edges_path, &all_edges)?;
        outputs.push(edges_path.clone());
    }

    eprintln!(
        "network influence: wrote {} influence rows over {} stratum(a)",
        all_rows.len(),
        strata_summary.len()
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_labeled_inputs(&[
        ("input", args.input.as_path()),
        ("samples", args.samples.as_path()),
    ])?;
    let sidecar = sidecar_path_for(&args.output);
    let policy_str = match policy {
        AdjacencyPolicy::HardThreshold { threshold } => format!("hard:{threshold}"),
        AdjacencyPolicy::SoftPower { power } => format!("soft:{power}"),
    };
    write_run_sidecar(
        &sidecar,
        "network influence",
        json!({
            "input": args.input.display().to_string(),
            "samples": args.samples.display().to_string(),
            "feature-col": args.feature_col,
            "method": args.method.to_core().as_str(),
            "adjacency-policy": policy_str,
            "stratify": args.stratify,
            "min-subjects": args.min_subjects,
            "eigen-max-iter": args.eigen_max_iter,
            "eigen-tol": args.eigen_tol,
            "emit-edges": args.emit_edges,
            "output": args.output.display().to_string(),
            "similarity_audit": {
                "n_pairs_total": audit_total.n_pairs_total,
                "n_pairs_insufficient_overlap": audit_total.n_pairs_insufficient_overlap,
                "n_pairs_undefined_metric": audit_total.n_pairs_undefined_metric,
                "per_stratum": per_stratum_audit
                    .iter()
                    .map(|(s, a)| {
                        json!({
                            "stratum": s,
                            "n_pairs_total": a.n_pairs_total,
                            "n_pairs_insufficient_overlap": a.n_pairs_insufficient_overlap,
                            "n_pairs_undefined_metric": a.n_pairs_undefined_metric,
                        })
                    })
                    .collect::<Vec<_>>(),
            },
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("network influence: sidecar={}", sidecar.display());
    Ok(())
}

struct FeatureMatrix {
    /// All distinct feature IDs seen, in sorted order.
    features: Vec<String>,
    /// `sample_id → (feature_id → abundance)`.
    per_sample: BTreeMap<String, BTreeMap<String, f64>>,
}

fn collect_feature_matrix(
    records: &[MeasurementRecord],
    feature_col: &str,
) -> Result<FeatureMatrix> {
    let mut features: BTreeSet<String> = BTreeSet::new();
    let mut per_sample: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    for r in records {
        if r.dropped_by_qc {
            continue;
        }
        let v = r.abundance.as_f64();
        if !v.is_finite() {
            continue;
        }
        let key = match feature_col {
            "gene_symbol" => r.gene_symbol.clone(),
            "assay_id" => Some(r.assay_id.0.clone()),
            _ => unreachable!(),
        };
        let Some(feature) = key else { continue };
        if feature.is_empty() {
            continue;
        }
        features.insert(feature.clone());
        per_sample
            .entry(r.sample_id.clone())
            .or_default()
            .insert(feature, v);
    }
    Ok(FeatureMatrix {
        features: features.into_iter().collect(),
        per_sample,
    })
}

fn read_stratum(path: &Path, column: &str) -> Result<BTreeMap<String, String>> {
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
    let stratum_col = headers.iter().position(|h| h == column).ok_or_else(|| {
        anyhow::anyhow!(
            "samples.tsv {:?} missing stratification column {:?}",
            path,
            column
        )
    })?;
    let mut out: BTreeMap<String, String> = BTreeMap::new();
    for row in reader.records() {
        let row = row.with_context(|| format!("reading {:?}", path))?;
        let sample = row.get(sample_col).unwrap_or_default().to_string();
        let stratum = row.get(stratum_col).unwrap_or_default().to_string();
        if sample.is_empty() || stratum.is_empty() {
            continue;
        }
        out.insert(sample, stratum);
    }
    Ok(out)
}

fn write_influence_tsv(path: &Path, rows: &[(String, InfluenceRow, usize)]) -> Result<()> {
    let mut buf = String::from(
        "feature_id\tstratum\teigenvector_centrality\tbetweenness_centrality\t\
         influence_score\tdegree\tn_subjects_used\n",
    );
    for (stratum, row, n_subjects) in rows {
        buf.push_str(&row.feature_id);
        buf.push('\t');
        buf.push_str(stratum);
        buf.push('\t');
        buf.push_str(&format_float(row.eigenvector_centrality));
        buf.push('\t');
        buf.push_str(&format_float(row.betweenness_centrality));
        buf.push('\t');
        buf.push_str(&format_float(row.influence_score));
        buf.push('\t');
        buf.push_str(&row.degree.to_string());
        buf.push('\t');
        buf.push_str(&n_subjects.to_string());
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

struct EdgeRow {
    feature_i: String,
    feature_j: String,
    edge_weight: f64,
    stratum: String,
}

fn write_edges_tsv(path: &Path, rows: &[EdgeRow]) -> Result<()> {
    let mut buf = String::from("feature_i\tfeature_j\tedge_weight\tstratum\n");
    for r in rows {
        buf.push_str(&r.feature_i);
        buf.push('\t');
        buf.push_str(&r.feature_j);
        buf.push('\t');
        buf.push_str(&format_float(r.edge_weight));
        buf.push('\t');
        buf.push_str(&r.stratum);
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}
