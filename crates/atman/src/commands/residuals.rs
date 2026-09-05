//! `atman residuals`: per-protein OLS residuals on nuisance covariates.
//!
//! Covariates come from `samples.tsv` plus any `--covariates-tsv` (for
//! example axis scores), and `--design` accepts covariate expressions
//! (`~ axis1_raw`, `~ log10(QAlb) + sex`). Each protein is fitted on the
//! complete-case design rows where it is observed (no imputation); the
//! optional variance tables report per-protein `R²` and the cohort-level
//! `Σss_model / Σss_total`, and `--output-canonical-dir` writes the
//! residual matrix as a canonical directory so `de` or `score weighted` can
//! run on it directly.

use anyhow::{bail, Context, Result};
use atman_core::contrast::{median, percentile};
use atman_core::de::{ols, OlsOutcome};
use atman_core::stats::mean;
use atman_core::{Abundance, MeasurementRecord, Sample};
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use crate::design::{build_design, parse_design, CollapseGenes, CovariateFrame, Term};
use crate::io::{
    atomic_write, format_f64, hash_canonical_inputs, hash_labeled_inputs, read_measurements_long,
    read_samples, sidecar_path_for, write_measurements_long, write_run_sidecar,
};
use std::time::SystemTime;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Nuisance-covariate formula, e.g. "~ age + sex" or "~ axis1_raw".
    /// Accepts covariate expressions (`log10()`, `z()`, arithmetic).
    #[arg(long)]
    design: String,

    /// Long residual output TSV.
    #[arg(long)]
    output: PathBuf,

    /// Optional wide residual matrix with proteins as rows and subjects as columns.
    #[arg(long)]
    output_wide: Option<PathBuf>,

    /// Minimum complete samples required per protein.
    #[arg(long, default_value_t = 3)]
    min_samples: usize,

    /// Extra `sample_id`-keyed covariate TSVs joined to samples.tsv (repeatable).
    #[arg(long, action = clap::ArgAction::Append)]
    covariates_tsv: Vec<PathBuf>,

    /// Drop proteins missing in more than this fraction of the design
    /// samples (1.0 = keep all).
    #[arg(long, default_value_t = 1.0)]
    max_missing_fraction: f64,

    /// Per-protein `n, ss_model, ss_total, r2` TSV.
    #[arg(long)]
    output_variance: Option<PathBuf>,

    /// One-row cohort summary: `frac_variance = sum(ss_model)/sum(ss_total)`,
    /// median/q75 R², fraction of proteins with R² > 0.25.
    #[arg(long)]
    output_variance_summary: Option<PathBuf>,

    /// Canonical directory whose `measurements.tsv` carries the residual as
    /// abundance (observed cells only) with `samples.tsv`/`proteins.tsv` copied.
    #[arg(long)]
    output_canonical_dir: Option<PathBuf>,

    /// How protein groups sharing a gene symbol are reduced: `none` keeps
    /// every assay as its own row; `mean` and `max-observed` collapse to one
    /// row per gene (assay_id = the representative assay). The missingness
    /// filter runs on the collapsed gene.
    #[arg(long, value_enum, default_value_t = CollapseGenes::None)]
    collapse_genes: CollapseGenes,
}

#[derive(Debug, Clone)]
struct SampleDesign {
    sample: Sample,
    row: Vec<f64>,
}

#[derive(Debug, Clone)]
struct ResidualRow {
    sample_id: String,
    subject_id: String,
    condition: String,
    assay_id: String,
    gene_symbol: String,
    residual: f64,
}

#[derive(Debug, Clone)]
struct VarianceRow {
    assay_id: String,
    gene_symbol: String,
    n: usize,
    ss_model: f64,
    ss_total: f64,
    r2: Option<f64>,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    if !(0.0..=1.0).contains(&args.max_missing_fraction) {
        bail!("--max-missing-fraction must be in [0, 1]");
    }
    let samples_path = args.input_dir.join("samples.tsv");
    let samples = read_samples(&samples_path)?;
    let mut frame = CovariateFrame::new();
    frame.load_tsv(&samples_path)?;
    for p in &args.covariates_tsv {
        frame.load_tsv(p)?;
    }
    let spec = parse_design(&args.design, None)?;
    if spec
        .terms
        .iter()
        .any(|t| matches!(t, Term::Column(c) if c == "condition"))
    {
        bail!("residuals expects nuisance covariates only; omit `condition` from --design");
    }
    if spec.terms.is_empty() {
        bail!("--design must contain at least one covariate besides the intercept");
    }
    let dm = build_design(
        &spec,
        samples.len(),
        &|i, col| frame.get(&samples[i].sample_id, col).map(str::to_string),
        &|_| None,
    )
    .with_context(|| format!("building the design for {}", args.design))?;
    let sample_design: Vec<SampleDesign> = dm
        .kept
        .iter()
        .enumerate()
        .map(|(k, &i)| SampleDesign {
            sample: samples[i].clone(),
            row: dm.rows[k].clone(),
        })
        .collect();
    if sample_design.is_empty() {
        bail!(
            "no samples have complete design metadata for {}",
            args.design
        );
    }

    let measurements_path = args.input_dir.join("measurements.tsv");
    let measurements = read_measurements_long(&measurements_path)?;
    let (rows, variance, n_dropped_missingness, collapse_stats) = compute_residuals(
        &sample_design,
        &measurements,
        args.min_samples,
        args.max_missing_fraction,
        args.collapse_genes,
    )?;
    if collapse_stats.n_genes_multi_assay > 0 {
        eprintln!(
            "residuals: {} gene symbols are carried by more than one assay; --collapse-genes {}",
            collapse_stats.n_genes_multi_assay,
            args.collapse_genes.as_str()
        );
    }
    write_long(&args.output, &rows)?;
    let mut outputs = vec![args.output.clone()];
    if let Some(path) = &args.output_wide {
        write_wide(path, &rows)?;
        outputs.push(path.clone());
    }
    if let Some(path) = &args.output_variance {
        write_variance(path, &variance)?;
        outputs.push(path.clone());
    }
    if let Some(path) = &args.output_variance_summary {
        write_variance_summary(path, &args.design, sample_design.len(), &variance)?;
        outputs.push(path.clone());
    }
    if let Some(dir) = &args.output_canonical_dir {
        outputs.extend(write_canonical_dir(
            dir,
            &args.input_dir,
            &measurements,
            &rows,
        )?);
    }
    eprintln!(
        "residuals: design={} proteins={} rows={} dropped_missingness={} output={}",
        args.design,
        variance.len(),
        rows.len(),
        n_dropped_missingness,
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let mut inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    for path in &args.covariates_tsv {
        let label = format!(
            "covariates_tsv/{}",
            path.file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("covariates")
        );
        inputs_sha256.extend(hash_labeled_inputs(&[(label.as_str(), path.as_path())])?);
    }
    let sidecar = sidecar_path_for(&args.output);
    let mut extras = serde_json::Map::new();
    extras.insert("n_design_samples".into(), json!(sample_design.len()));
    extras.insert("n_dropped_missingness".into(), json!(n_dropped_missingness));
    extras.insert(
        "gene_symbol_collapse".into(),
        json!({
            "rule": args.collapse_genes.as_str(),
            "n_assays": collapse_stats.n_assays,
            "n_genes_with_multiple_assays": collapse_stats.n_genes_multi_assay,
            "n_genes_after_collapse": collapse_stats.n_genes_after_collapse,
            "n_rows_after_collapse": collapse_stats.n_rows_after_collapse,
        }),
    );
    write_run_sidecar(
        &sidecar,
        "residuals",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "design": args.design,
            "output": args.output.display().to_string(),
            "output-wide": args.output_wide.as_ref().map(|p| p.display().to_string()),
            "min-samples": args.min_samples,
            "covariates-tsv": args.covariates_tsv.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
            "max-missing-fraction": args.max_missing_fraction,
            "output-variance": args.output_variance.as_ref().map(|p| p.display().to_string()),
            "output-variance-summary": args.output_variance_summary.as_ref().map(|p| p.display().to_string()),
            "output-canonical-dir": args.output_canonical_dir.as_ref().map(|p| p.display().to_string()),
            "collapse-genes": args.collapse_genes.as_str(),
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        Some(extras),
    )?;
    eprintln!("residuals: sidecar={}", sidecar.display());
    Ok(())
}

type ProteinValues<'a> = BTreeMap<(String, String), Vec<(&'a SampleDesign, f64)>>;

/// Gene-symbol collapse bookkeeping for the sidecar.
#[derive(Debug, Clone, Copy, Default)]
struct CollapseStats {
    /// Distinct assays measured on the design samples.
    n_assays: usize,
    /// Gene symbols carried by more than one assay.
    n_genes_multi_assay: usize,
    /// Distinct gene symbols after the collapse.
    n_genes_after_collapse: usize,
    /// Protein rows entering the fit (equals genes unless `none` keeps
    /// several assays per symbol), before the missingness filter.
    n_rows_after_collapse: usize,
}

fn compute_residuals(
    sample_design: &[SampleDesign],
    measurements: &[MeasurementRecord],
    min_samples: usize,
    max_missing_fraction: f64,
    collapse: CollapseGenes,
) -> Result<(Vec<ResidualRow>, Vec<VarianceRow>, usize, CollapseStats)> {
    let design_by_sample: HashMap<&str, usize> = sample_design
        .iter()
        .enumerate()
        .map(|(i, row)| (row.sample.sample_id.as_str(), i))
        .collect();
    // gene → design index → [(assay_id, value)]
    let mut raw: BTreeMap<String, BTreeMap<usize, Vec<(String, f64)>>> = BTreeMap::new();
    for measurement in measurements {
        let Some(value) = measurement.effective_abundance() else {
            continue;
        };
        let Some(&idx) = design_by_sample.get(measurement.sample_id.as_str()) else {
            continue;
        };
        let gene = measurement
            .gene_symbol
            .clone()
            .unwrap_or_else(|| measurement.assay_id.0.clone());
        raw.entry(gene)
            .or_default()
            .entry(idx)
            .or_default()
            .push((measurement.assay_id.0.clone(), value));
    }
    let mut by_protein: ProteinValues = BTreeMap::new();
    let mut stats = CollapseStats {
        n_genes_after_collapse: raw.len(),
        ..CollapseStats::default()
    };
    for (gene, by_sample) in raw {
        let mut assay_counts: BTreeMap<String, usize> = BTreeMap::new();
        for values in by_sample.values() {
            for (assay, _) in values {
                *assay_counts.entry(assay.clone()).or_default() += 1;
            }
        }
        stats.n_assays += assay_counts.len();
        if assay_counts.len() > 1 {
            stats.n_genes_multi_assay += 1;
        }
        if collapse == CollapseGenes::None {
            for (idx, values) in &by_sample {
                for (assay, value) in values {
                    by_protein
                        .entry((assay.clone(), gene.clone()))
                        .or_default()
                        .push((&sample_design[*idx], *value));
                }
            }
            continue;
        }
        let representative = collapse
            .representative(&assay_counts)
            .cloned()
            .unwrap_or_default();
        for (idx, values) in &by_sample {
            if let Some(value) = collapse.collapse(values, &representative) {
                by_protein
                    .entry((representative.clone(), gene.clone()))
                    .or_default()
                    .push((&sample_design[*idx], value));
            }
        }
    }

    stats.n_rows_after_collapse = by_protein.len();
    let n_design = sample_design.len();
    let mut out = Vec::new();
    let mut variance = Vec::new();
    let mut n_dropped = 0usize;
    for ((assay_id, gene_symbol), values) in by_protein {
        let missing = 1.0 - values.len() as f64 / n_design.max(1) as f64;
        if missing > max_missing_fraction + 1e-12 {
            n_dropped += 1;
            continue;
        }
        let design: Vec<Vec<f64>> = values
            .iter()
            .map(|(sample, _)| sample.row.clone())
            .collect();
        let y: Vec<f64> = values.iter().map(|(_, value)| *value).collect();
        let fit = match ols(&design, &y, min_samples) {
            OlsOutcome::Computed(fit) => fit,
            OlsOutcome::Skipped { .. } => continue,
        };
        let ybar = mean(&y);
        let mut ss_total = 0.0;
        let mut ss_model = 0.0;
        for ((sample, observed), row) in values.iter().zip(design.iter()) {
            let fitted: f64 = row
                .iter()
                .zip(fit.beta.iter())
                .map(|(x, beta)| x * beta)
                .sum();
            ss_total += (observed - ybar).powi(2);
            ss_model += (fitted - ybar).powi(2);
            out.push(ResidualRow {
                sample_id: sample.sample.sample_id.clone(),
                subject_id: sample.sample.subject_id.clone().unwrap_or_default(),
                condition: sample.sample.condition.clone().unwrap_or_default(),
                assay_id: assay_id.clone(),
                gene_symbol: gene_symbol.clone(),
                residual: observed - fitted,
            });
        }
        variance.push(VarianceRow {
            assay_id,
            gene_symbol,
            n: y.len(),
            ss_model,
            ss_total,
            r2: if ss_total > 0.0 {
                Some(ss_model / ss_total)
            } else {
                None
            },
        });
    }
    out.sort_by(|a, b| {
        a.gene_symbol
            .cmp(&b.gene_symbol)
            .then_with(|| a.assay_id.cmp(&b.assay_id))
            .then_with(|| a.sample_id.cmp(&b.sample_id))
    });
    variance.sort_by(|a, b| {
        a.gene_symbol
            .cmp(&b.gene_symbol)
            .then_with(|| a.assay_id.cmp(&b.assay_id))
    });
    Ok((out, variance, n_dropped, stats))
}

fn write_long(path: &Path, rows: &[ResidualRow]) -> Result<()> {
    let mut out =
        String::from("sample_id\tsubject_id\tcondition\tassay_id\tgene_symbol\tresidual\n");
    for row in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            row.sample_id,
            row.subject_id,
            row.condition,
            row.assay_id,
            row.gene_symbol,
            fmt(row.residual)
        ));
    }
    atomic_write(path, out.as_bytes())
}

fn write_wide(path: &Path, rows: &[ResidualRow]) -> Result<()> {
    let samples: Vec<String> = rows
        .iter()
        .map(|row| row.sample_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut by_protein: BTreeMap<(String, String), BTreeMap<String, Vec<f64>>> = BTreeMap::new();
    for row in rows {
        by_protein
            .entry((row.assay_id.clone(), row.gene_symbol.clone()))
            .or_default()
            .entry(row.sample_id.clone())
            .or_default()
            .push(row.residual);
    }
    let mut out = String::from("assay_id\tgene_symbol");
    for sample in &samples {
        out.push('\t');
        out.push_str(sample);
    }
    out.push('\n');
    for ((assay_id, gene_symbol), per_sample) in by_protein {
        out.push_str(&assay_id);
        out.push('\t');
        out.push_str(&gene_symbol);
        for sample in &samples {
            out.push('\t');
            if let Some(values) = per_sample.get(sample) {
                out.push_str(&fmt(mean(values)));
            } else {
                out.push('0');
            }
        }
        out.push('\n');
    }
    atomic_write(path, out.as_bytes())
}

fn write_variance(path: &Path, rows: &[VarianceRow]) -> Result<()> {
    let mut out = String::from("assay_id\tgene_symbol\tn\tss_model\tss_total\tr2\n");
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            r.assay_id,
            r.gene_symbol,
            r.n,
            format_f64(r.ss_model),
            format_f64(r.ss_total),
            r.r2.map(format_f64).unwrap_or_default()
        ));
    }
    atomic_write(path, out.as_bytes())
}

fn write_variance_summary(
    path: &Path,
    design: &str,
    n_samples: usize,
    rows: &[VarianceRow],
) -> Result<()> {
    let with_r2: Vec<&VarianceRow> = rows.iter().filter(|r| r.r2.is_some()).collect();
    let n_proteins = with_r2.len();
    let ss_model: f64 = with_r2.iter().map(|r| r.ss_model).sum();
    let ss_total: f64 = with_r2.iter().map(|r| r.ss_total).sum();
    let frac_variance = if ss_total > 0.0 {
        Some(ss_model / ss_total)
    } else {
        None
    };
    let mut r2s: Vec<f64> = with_r2.iter().filter_map(|r| r.r2).collect();
    r2s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_r2 = median(&r2s);
    let q75_r2 = percentile(&r2s, 0.75);
    let frac_gt = if n_proteins > 0 {
        Some(r2s.iter().filter(|v| **v > 0.25).count() as f64 / n_proteins as f64)
    } else {
        None
    };
    let f = |v: Option<f64>| v.map(format_f64).unwrap_or_default();
    let out = format!(
        "design\tn_samples\tn_proteins\tfrac_variance\tmedian_r2\tq75_r2\tfrac_r2_gt_0_25\n{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        design,
        n_samples,
        n_proteins,
        f(frac_variance),
        f(median_r2),
        f(q75_r2),
        f(frac_gt)
    );
    atomic_write(path, out.as_bytes())
}

fn with_value(abundance: &Abundance, value: f64) -> Abundance {
    match abundance {
        Abundance::Log2Npx(_) => Abundance::Log2Npx(value),
        Abundance::Log2Intensity(_) => Abundance::Log2Intensity(value),
        Abundance::Ibaq(_) => Abundance::Ibaq(value),
        Abundance::Raw(_) => Abundance::Raw(value),
    }
}

/// Write the residual matrix as a canonical directory: `measurements.tsv`
/// carries the residual as abundance for every fitted (sample, assay)
/// cell; `samples.tsv` and `proteins.tsv` are byte copies of the input.
fn write_canonical_dir(
    dir: &Path,
    input_dir: &Path,
    measurements: &[MeasurementRecord],
    rows: &[ResidualRow],
) -> Result<Vec<PathBuf>> {
    std::fs::create_dir_all(dir).with_context(|| format!("creating {:?}", dir))?;
    // Template records: exact (sample, assay) match first; otherwise any
    // record of that (sample, gene) when the residual is a collapsed gene
    // whose representative assay was not measured in that sample.
    let mut by_sample_assay: HashMap<(&str, &str), &MeasurementRecord> = HashMap::new();
    let mut by_sample_gene: HashMap<(&str, &str), &MeasurementRecord> = HashMap::new();
    for m in measurements {
        by_sample_assay
            .entry((m.sample_id.as_str(), m.assay_id.0.as_str()))
            .or_insert(m);
        let gene = m.gene_symbol.as_deref().unwrap_or(m.assay_id.0.as_str());
        by_sample_gene
            .entry((m.sample_id.as_str(), gene))
            .or_insert(m);
    }
    let mut records: Vec<MeasurementRecord> = Vec::with_capacity(rows.len());
    for row in rows {
        let template = by_sample_assay
            .get(&(row.sample_id.as_str(), row.assay_id.as_str()))
            .or_else(|| by_sample_gene.get(&(row.sample_id.as_str(), row.gene_symbol.as_str())));
        let Some(&m) = template else {
            continue;
        };
        let mut r = m.clone();
        r.assay_id = atman_core::AssayId(row.assay_id.clone());
        r.gene_symbol = Some(row.gene_symbol.clone());
        r.abundance = with_value(&m.abundance, row.residual);
        r.abundance_raw = with_value(&m.abundance_raw, row.residual);
        r.npx_source_str = format_f64(row.residual);
        r.dropped_by_qc = false;
        records.push(r);
    }
    records.sort_by(|a, b| {
        a.ingest_order
            .cmp(&b.ingest_order)
            .then_with(|| a.assay_id.0.cmp(&b.assay_id.0))
    });
    let measurements_path = dir.join("measurements.tsv");
    write_measurements_long(&measurements_path, &records)?;
    let mut written = vec![measurements_path];
    for name in ["samples.tsv", "proteins.tsv"] {
        let src = input_dir.join(name);
        let dst = dir.join(name);
        std::fs::copy(&src, &dst).with_context(|| format!("copying {:?} to {:?}", src, dst))?;
        written.push(dst);
    }
    Ok(written)
}

fn fmt(value: f64) -> String {
    if value.abs() < 5e-13 {
        "0".to_string()
    } else {
        format!("{value:.12}")
    }
}
