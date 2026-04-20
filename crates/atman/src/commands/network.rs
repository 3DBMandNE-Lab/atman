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
    /// `qc_measurements.tsv`) with `sample_id, assay_id, gene_symbol,
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
    }
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
    let (features, per_sample) = collect_feature_matrix(&records, &args.feature_col)?;
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
        let Some(similarity) = pairwise_similarity(&features, &data, args.method.to_core()) else {
            eprintln!(
                "network influence: stratum {:?} produced no similarity matrix (insufficient data)",
                stratum
            );
            continue;
        };
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
    let input_dir_sha256 = hash_labeled_inputs(&[
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
        }),
        &input_dir_sha256,
        &outputs,
        started_at,
        finished_at,
    )?;
    eprintln!("network influence: sidecar={}", sidecar.display());
    Ok(())
}

fn collect_feature_matrix(
    records: &[MeasurementRecord],
    feature_col: &str,
) -> Result<(Vec<String>, BTreeMap<String, BTreeMap<String, f64>>)> {
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
    Ok((features.into_iter().collect(), per_sample))
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
    let stratum_col = headers
        .iter()
        .position(|h| h == column)
        .ok_or_else(|| anyhow::anyhow!(
            "samples.tsv {:?} missing stratification column {:?}",
            path,
            column
        ))?;
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

fn write_influence_tsv(
    path: &Path,
    rows: &[(String, InfluenceRow, usize)],
) -> Result<()> {
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
