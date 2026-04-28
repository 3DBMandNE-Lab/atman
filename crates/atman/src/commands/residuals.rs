use anyhow::{anyhow, bail, Context, Result};
use atman_core::de::{ols, OlsOutcome};
use atman_core::stats::mean;
use atman_core::{MeasurementRecord, Sample};
use clap::Args as ClapArgs;
use csv::ReaderBuilder;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use crate::io::{
    atomic_write, hash_canonical_inputs, need_col, read_measurements_long, read_samples,
    sidecar_path_for, write_run_sidecar,
};
use serde_json::json;
use std::time::SystemTime;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Nuisance-covariate formula, e.g. "~ age + sex".
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
enum CovariateSpec {
    Numeric { name: String },
    Categorical { name: String, levels: Vec<String> },
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    let terms = parse_design_terms(&args.design)?;
    let samples_path = args.input_dir.join("samples.tsv");
    let samples = read_samples(&samples_path)?;
    let metadata = read_sample_metadata(&samples_path, &terms)?;
    let specs = infer_covariates(&terms, &metadata)?;
    let sample_design = build_sample_design(samples, &metadata, &specs)?;
    if sample_design.is_empty() {
        bail!(
            "no samples have complete design metadata for {}",
            args.design
        );
    }

    let measurements_path = args.input_dir.join("measurements.tsv");
    let measurements = read_measurements_long(&measurements_path)?;
    let rows = compute_residuals(&sample_design, &measurements, args.min_samples)?;
    write_long(&args.output, &rows)?;
    if let Some(path) = &args.output_wide {
        write_wide(path, &rows)?;
    }
    eprintln!(
        "residuals: design={} proteins={} rows={} output={}",
        args.design,
        rows.iter()
            .map(|r| (r.assay_id.as_str(), r.gene_symbol.as_str()))
            .collect::<BTreeSet<_>>()
            .len(),
        rows.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    let mut outputs = vec![args.output.clone()];
    if let Some(p) = &args.output_wide {
        outputs.push(p.clone());
    }
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "residuals",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "design": args.design,
            "output": args.output.display().to_string(),
            "output-wide": args.output_wide.as_ref().map(|p| p.display().to_string()),
            "min-samples": args.min_samples,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("residuals: sidecar={}", sidecar.display());
    Ok(())
}

fn parse_design_terms(design: &str) -> Result<Vec<String>> {
    let rhs = design
        .trim()
        .strip_prefix('~')
        .ok_or_else(|| anyhow!("--design must start with `~`"))?
        .trim();
    if rhs.is_empty() {
        bail!("--design must contain at least one covariate");
    }
    let mut terms = Vec::new();
    let mut seen = BTreeSet::new();
    for raw in rhs.split('+') {
        let term = raw.trim();
        if term.is_empty() {
            bail!("empty term in --design");
        }
        if term == "1" {
            continue;
        }
        if term == "0" || term == "-1" {
            bail!("intercept removal is not supported; Atman includes an intercept");
        }
        if term == "condition" {
            bail!("residuals expects nuisance covariates only; omit `condition` from --design");
        }
        if !seen.insert(term.to_string()) {
            bail!("duplicate term {:?} in --design", term);
        }
        terms.push(term.to_string());
    }
    if terms.is_empty() {
        bail!("--design must contain at least one covariate besides the intercept");
    }
    Ok(terms)
}

fn read_sample_metadata(
    path: &Path,
    terms: &[String],
) -> Result<BTreeMap<String, BTreeMap<String, String>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let sample_col = need_col(&headers, "sample_id", path)?;
    let mut term_cols = Vec::new();
    for term in terms {
        term_cols.push((
            term.clone(),
            need_col(&headers, term, path)
                .with_context(|| format!("required by residuals --design {}", term))?,
        ));
    }
    let mut out = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let mut values = BTreeMap::new();
        for (term, col) in &term_cols {
            values.insert(term.clone(), row[*col].trim().to_string());
        }
        out.insert(row[sample_col].to_string(), values);
    }
    Ok(out)
}

fn infer_covariates(
    terms: &[String],
    metadata: &BTreeMap<String, BTreeMap<String, String>>,
) -> Result<Vec<CovariateSpec>> {
    let mut specs = Vec::new();
    for term in terms {
        let values: Vec<&str> = metadata
            .values()
            .filter_map(|row| row.get(term).map(String::as_str))
            .filter(|value| !value.is_empty())
            .collect();
        if values.is_empty() {
            bail!("covariate {:?} has no non-empty values", term);
        }
        if values.iter().all(|value| value.parse::<f64>().is_ok()) {
            specs.push(CovariateSpec::Numeric { name: term.clone() });
        } else {
            let levels: Vec<String> = values
                .into_iter()
                .map(str::to_string)
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            if levels.len() < 2 {
                bail!("categorical covariate {:?} has fewer than two levels", term);
            }
            specs.push(CovariateSpec::Categorical {
                name: term.clone(),
                levels,
            });
        }
    }
    Ok(specs)
}

fn build_sample_design(
    samples: Vec<Sample>,
    metadata: &BTreeMap<String, BTreeMap<String, String>>,
    specs: &[CovariateSpec],
) -> Result<Vec<SampleDesign>> {
    let mut out = Vec::new();
    for sample in samples {
        let Some(meta) = metadata.get(&sample.sample_id) else {
            continue;
        };
        let mut row = vec![1.0];
        let mut complete = true;
        for spec in specs {
            match spec {
                CovariateSpec::Numeric { name } => {
                    let value = meta.get(name).map(String::as_str).unwrap_or("").trim();
                    if value.is_empty() {
                        complete = false;
                        break;
                    }
                    row.push(value.parse()?);
                }
                CovariateSpec::Categorical { name, levels } => {
                    let value = meta.get(name).map(String::as_str).unwrap_or("").trim();
                    if value.is_empty() {
                        complete = false;
                        break;
                    }
                    if !levels.iter().any(|level| level == value) {
                        bail!("unknown level {:?} for covariate {:?}", value, name);
                    }
                    for level in levels.iter().skip(1) {
                        row.push((value == level) as u8 as f64);
                    }
                }
            }
        }
        if complete {
            out.push(SampleDesign { sample, row });
        }
    }
    Ok(out)
}

fn compute_residuals(
    sample_design: &[SampleDesign],
    measurements: &[MeasurementRecord],
    min_samples: usize,
) -> Result<Vec<ResidualRow>> {
    let design_by_sample: HashMap<&str, &SampleDesign> = sample_design
        .iter()
        .map(|row| (row.sample.sample_id.as_str(), row))
        .collect();
    let mut by_protein: BTreeMap<(String, String), Vec<(&SampleDesign, f64)>> = BTreeMap::new();
    for measurement in measurements {
        let Some(value) = measurement.effective_abundance() else {
            continue;
        };
        let Some(sample) = design_by_sample
            .get(measurement.sample_id.as_str())
            .copied()
        else {
            continue;
        };
        let gene = measurement
            .gene_symbol
            .clone()
            .unwrap_or_else(|| measurement.assay_id.0.clone());
        by_protein
            .entry((measurement.assay_id.0.clone(), gene))
            .or_default()
            .push((sample, value));
    }

    let mut out = Vec::new();
    for ((assay_id, gene_symbol), values) in by_protein {
        let design: Vec<Vec<f64>> = values
            .iter()
            .map(|(sample, _)| sample.row.clone())
            .collect();
        let y: Vec<f64> = values.iter().map(|(_, value)| *value).collect();
        let fit = match ols(&design, &y, min_samples) {
            OlsOutcome::Computed(fit) => fit,
            OlsOutcome::Skipped { .. } => continue,
        };
        for ((sample, observed), row) in values.iter().zip(design.iter()) {
            let fitted: f64 = row
                .iter()
                .zip(fit.beta.iter())
                .map(|(x, beta)| x * beta)
                .sum();
            out.push(ResidualRow {
                sample_id: sample.sample.sample_id.clone(),
                subject_id: sample.sample.subject_id.clone().unwrap_or_default(),
                condition: sample.sample.condition.clone().unwrap_or_default(),
                assay_id: assay_id.clone(),
                gene_symbol: gene_symbol.clone(),
                residual: observed - fitted,
            });
        }
    }
    out.sort_by(|a, b| {
        a.gene_symbol
            .cmp(&b.gene_symbol)
            .then_with(|| a.assay_id.cmp(&b.assay_id))
            .then_with(|| a.sample_id.cmp(&b.sample_id))
    });
    Ok(out)
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

fn fmt(value: f64) -> String {
    if value.abs() < 5e-13 {
        "0".to_string()
    } else {
        format!("{value:.12}")
    }
}
