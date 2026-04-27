use anyhow::{bail, Context, Result};
use atman_core::stats::mean;
use atman_core::Sample;
use clap::{Args as ClapArgs, Subcommand};
use csv::ReaderBuilder;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use crate::io::{atomic_write, read_measurements_long, read_samples};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Score user-defined modules per sample.
    Modules(ModulesArgs),
}

#[derive(ClapArgs, Debug)]
pub struct ModulesArgs {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Modules TSV with columns `module` and `gene_symbol`.
    #[arg(long, alias = "modules")]
    modules_tsv: PathBuf,

    /// Output TSV path for per-sample module scores.
    #[arg(long)]
    output: PathBuf,

    /// Scoring method: mean, median, zscore, or pc1.
    #[arg(long, default_value = "mean")]
    method: String,

    /// Optional output directory for canonical TSVs usable by `atman de`.
    #[arg(long)]
    canonical_output_dir: Option<PathBuf>,
}

#[derive(Debug, Clone)]
struct ScoreRow {
    sample_id: String,
    subject_id: String,
    condition: String,
    module: String,
    method: String,
    score: f64,
    n_genes_declared: usize,
    n_genes_observed: usize,
    coverage: f64,
}

#[derive(Default)]
struct Acc {
    sum: f64,
    n: usize,
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Modules(args) => run_modules(args),
    }
}

fn run_modules(args: ModulesArgs) -> Result<()> {
    if !matches!(args.method.as_str(), "mean" | "median" | "zscore" | "pc1") {
        bail!("score modules supports --method mean, median, zscore, or pc1");
    }
    let modules = read_modules(&args.modules_tsv)?;
    if modules.is_empty() {
        bail!("no module definitions found in {:?}", args.modules_tsv);
    }
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();
    let measurements_path = args.input_dir.join("measurements.tsv");
    let measurements = read_measurements_long(&measurements_path)?;

    let mut gene_to_modules: HashMap<&str, Vec<&str>> = HashMap::new();
    for (module, genes) in &modules {
        for gene in genes {
            gene_to_modules
                .entry(gene.as_str())
                .or_default()
                .push(module.as_str());
        }
    }

    let mut values: BTreeMap<(String, String, String), Acc> = BTreeMap::new();
    let mut observed: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for m in &measurements {
        let Some(value) = m.effective_abundance() else {
            continue;
        };
        let Some(gene) = m.gene_symbol.as_deref() else {
            continue;
        };
        let Some(module_names) = gene_to_modules.get(gene) else {
            continue;
        };
        for module in module_names {
            observed
                .entry((*module).to_string())
                .or_default()
                .insert(gene.to_string());
            let acc = values
                .entry(((*module).to_string(), m.sample_id.clone(), gene.to_string()))
                .or_default();
            acc.sum += value;
            acc.n += 1;
        }
    }

    let mut rows = Vec::new();
    for (module, genes) in &modules {
        let scores = score_module(module, genes, &samples, &values, &args.method);
        let n_observed = observed.get(module).map(BTreeSet::len).unwrap_or(0);
        for (sample_id, score, n_sample_genes) in scores {
            let Some(sample) = sample_by_id.get(sample_id.as_str()) else {
                continue;
            };
            rows.push(ScoreRow {
                sample_id,
                subject_id: sample.subject_id.clone().unwrap_or_default(),
                condition: sample.condition.clone().unwrap_or_default(),
                module: module.clone(),
                method: args.method.clone(),
                score,
                n_genes_declared: genes.len(),
                n_genes_observed: n_observed,
                coverage: n_sample_genes as f64 / genes.len() as f64,
            });
        }
    }
    rows.sort_by(|a, b| {
        a.module
            .cmp(&b.module)
            .then_with(|| a.sample_id.cmp(&b.sample_id))
    });
    write_scores(&args.output, &rows)?;
    if let Some(dir) = &args.canonical_output_dir {
        write_canonical(dir, &args.input_dir, &rows, &modules)?;
    }
    eprintln!(
        "score modules: method={} modules={} rows={}",
        args.method,
        modules.len(),
        rows.len()
    );
    Ok(())
}

fn score_module(
    module: &str,
    genes: &BTreeSet<String>,
    samples: &[Sample],
    values: &BTreeMap<(String, String, String), Acc>,
    method: &str,
) -> Vec<(String, f64, usize)> {
    let sample_values = |sample_id: &str| -> Vec<(String, f64)> {
        genes
            .iter()
            .filter_map(|gene| {
                values
                    .get(&(module.to_string(), sample_id.to_string(), gene.clone()))
                    .map(|acc| (gene.clone(), acc.sum / acc.n as f64))
            })
            .collect()
    };
    if method == "zscore" || method == "pc1" {
        let stats = gene_stats(module, genes, samples, values);
        if method == "pc1" {
            return pc1_scores(module, genes, samples, values, &stats);
        }
        return samples
            .iter()
            .filter_map(|sample| {
                let zs: Vec<f64> = sample_values(&sample.sample_id)
                    .into_iter()
                    .filter_map(|(gene, value)| {
                        let (mu, sd) = stats.get(&gene)?;
                        if *sd > 0.0 {
                            Some((value - mu) / sd)
                        } else {
                            None
                        }
                    })
                    .collect();
                if zs.is_empty() {
                    None
                } else {
                    Some((sample.sample_id.clone(), mean(&zs), zs.len()))
                }
            })
            .collect();
    }
    samples
        .iter()
        .filter_map(|sample| {
            let vals: Vec<f64> = sample_values(&sample.sample_id)
                .into_iter()
                .map(|(_, v)| v)
                .collect();
            if vals.is_empty() {
                None
            } else if method == "median" {
                median(vals.clone()).map(|score| (sample.sample_id.clone(), score, vals.len()))
            } else {
                Some((sample.sample_id.clone(), mean(&vals), vals.len()))
            }
        })
        .collect()
}

fn gene_stats(
    module: &str,
    genes: &BTreeSet<String>,
    samples: &[Sample],
    values: &BTreeMap<(String, String, String), Acc>,
) -> BTreeMap<String, (f64, f64)> {
    let mut out = BTreeMap::new();
    for gene in genes {
        let vals: Vec<f64> = samples
            .iter()
            .filter_map(|sample| {
                values
                    .get(&(module.to_string(), sample.sample_id.clone(), gene.clone()))
                    .map(|acc| acc.sum / acc.n as f64)
            })
            .collect();
        if vals.len() >= 2 {
            let mu = mean(&vals);
            let sd = sample_sd(&vals, mu);
            out.insert(gene.clone(), (mu, sd));
        }
    }
    out
}

fn pc1_scores(
    module: &str,
    genes: &BTreeSet<String>,
    samples: &[Sample],
    values: &BTreeMap<(String, String, String), Acc>,
    stats: &BTreeMap<String, (f64, f64)>,
) -> Vec<(String, f64, usize)> {
    let genes_used: Vec<&String> = genes
        .iter()
        .filter(|gene| stats.get(*gene).map(|(_, sd)| *sd > 0.0).unwrap_or(false))
        .collect();
    if genes_used.is_empty() {
        return Vec::new();
    }
    let mut matrix = Vec::new();
    let mut sample_ids = Vec::new();
    let mut coverage = Vec::new();
    for sample in samples {
        let mut row = Vec::with_capacity(genes_used.len());
        let mut n_seen = 0;
        for gene in &genes_used {
            let (mu, sd) = stats.get(*gene).expect("filtered above");
            let z = values
                .get(&(
                    module.to_string(),
                    sample.sample_id.clone(),
                    (*gene).clone(),
                ))
                .map(|acc| {
                    n_seen += 1;
                    ((acc.sum / acc.n as f64) - mu) / sd
                })
                .unwrap_or(0.0);
            row.push(z);
        }
        if n_seen > 0 {
            matrix.push(row);
            sample_ids.push(sample.sample_id.clone());
            coverage.push(n_seen);
        }
    }
    if matrix.is_empty() {
        return Vec::new();
    }
    let mut loading = vec![1.0 / (genes_used.len() as f64).sqrt(); genes_used.len()];
    for _ in 0..50 {
        let mut next = vec![0.0; genes_used.len()];
        for row in &matrix {
            let projection: f64 = row.iter().zip(loading.iter()).map(|(x, l)| x * l).sum();
            for (idx, value) in row.iter().enumerate() {
                next[idx] += value * projection;
            }
        }
        let norm = next.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm == 0.0 {
            break;
        }
        for value in &mut next {
            *value /= norm;
        }
        loading = next;
    }
    if loading.iter().sum::<f64>() < 0.0 {
        for value in &mut loading {
            *value = -*value;
        }
    }
    sample_ids
        .into_iter()
        .zip(matrix)
        .zip(coverage)
        .map(|((sample_id, row), n)| {
            let score = row.iter().zip(loading.iter()).map(|(x, l)| x * l).sum();
            (sample_id, score, n)
        })
        .collect()
}

fn write_scores(path: &Path, rows: &[ScoreRow]) -> Result<()> {
    let mut out = String::from(
        "sample_id\tsubject_id\tcondition\tmodule\tmethod\tscore\tn_genes_declared\tn_genes_observed\tcoverage\n",
    );
    for row in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.sample_id,
            row.subject_id,
            row.condition,
            row.module,
            row.method,
            row.score,
            row.n_genes_declared,
            row.n_genes_observed,
            row.coverage,
        ));
    }
    atomic_write(path, out.as_bytes())
}

fn write_canonical(
    output_dir: &Path,
    input_dir: &Path,
    rows: &[ScoreRow],
    modules: &BTreeMap<String, BTreeSet<String>>,
) -> Result<()> {
    std::fs::create_dir_all(output_dir).with_context(|| format!("creating {:?}", output_dir))?;
    let samples = std::fs::read(input_dir.join("samples.tsv"))
        .with_context(|| format!("reading {:?}", input_dir.join("samples.tsv")))?;
    atomic_write(&output_dir.join("samples.tsv"), &samples)?;

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for module in modules.keys() {
        proteins.push_str(&format!("somascan\t{module}\t\t{module}\tmodule_score\t\n"));
    }
    atomic_write(&output_dir.join("proteins.tsv"), proteins.as_bytes())?;

    let mut measurements = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    for (idx, row) in rows.iter().enumerate() {
        measurements.push_str(&format!(
            "somascan\t{}\t{}\t{}\tmodule_score\t{}\t{}\t{}\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            row.sample_id,
            row.module,
            row.module,
            row.score,
            row.score,
            row.score,
            idx + 1
        ));
    }
    atomic_write(
        &output_dir.join("measurements.tsv"),
        measurements.as_bytes(),
    )
}

fn read_modules(path: &Path) -> Result<BTreeMap<String, BTreeSet<String>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let module_col = headers
        .iter()
        .position(|h| h == "module")
        .with_context(|| format!("missing column `module` in {:?}", path))?;
    let gene_col = headers
        .iter()
        .position(|h| h == "gene_symbol")
        .with_context(|| format!("missing column `gene_symbol` in {:?}", path))?;
    let mut out: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let module = row[module_col].trim();
        let gene = row[gene_col].trim();
        if module.is_empty() || gene.is_empty() {
            continue;
        }
        out.entry(module.to_string())
            .or_default()
            .insert(gene.to_string());
    }
    Ok(out)
}

fn sample_sd(values: &[f64], mean: f64) -> f64 {
    (values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64).sqrt()
}

fn median(mut values: Vec<f64>) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        Some(values[mid])
    } else {
        Some((values[mid - 1] + values[mid]) / 2.0)
    }
}
