use anyhow::{bail, Context, Result};
use atman_core::stats::mean;
use atman_core::Sample;
use clap::{Args as ClapArgs, Subcommand};
use csv::ReaderBuilder;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use super::parse_comparisons;
use crate::io::{
    atomic_write, hash_canonical_inputs, hash_labeled_inputs, read_measurements_long,
    read_proteins, read_samples, sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Subject-level bootstrap confidence intervals for protein effects.
    Protein(ProteinArgs),
    /// Subject-level bootstrap confidence intervals for module effects.
    Module(ModuleArgs),
    /// Subject-level bootstrap confidence intervals for signed-loading program effects.
    Program(ProgramArgs),
}

#[derive(ClapArgs, Debug)]
pub struct ProteinArgs {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output TSV path.
    #[arg(long)]
    output: PathBuf,

    /// Comma-separated comparisons in A-B form.
    #[arg(long)]
    groups: String,

    /// Bootstrap mode: welch-t for unpaired or paired-t for matched subjects.
    #[arg(long, default_value = "welch-t")]
    test: String,

    /// Number of bootstrap replicates.
    #[arg(long, default_value_t = 2000)]
    n: usize,

    /// Minimum subjects per group, or matched subjects in paired mode.
    #[arg(long, default_value_t = 2)]
    min_pairs: usize,

    /// Deterministic RNG seed.
    #[arg(long, default_value_t = 1)]
    seed: u64,

    /// Lower CI quantile.
    #[arg(long, default_value_t = 0.025)]
    ci_low: f64,

    /// Upper CI quantile.
    #[arg(long, default_value_t = 0.975)]
    ci_high: f64,
}

#[derive(ClapArgs, Debug)]
pub struct ModuleArgs {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Modules TSV with columns `module` and `gene_symbol`.
    #[arg(long)]
    modules_tsv: PathBuf,

    /// Output TSV path.
    #[arg(long)]
    output: PathBuf,

    /// Comma-separated comparisons in A-B form.
    #[arg(long)]
    groups: String,

    /// Bootstrap mode: welch-t for unpaired or paired-t for matched subjects.
    #[arg(long, default_value = "welch-t")]
    test: String,

    /// Number of bootstrap replicates.
    #[arg(long, default_value_t = 2000)]
    n: usize,

    /// Minimum subjects per group, or matched subjects in paired mode.
    #[arg(long, default_value_t = 2)]
    min_pairs: usize,

    /// Deterministic RNG seed.
    #[arg(long, default_value_t = 1)]
    seed: u64,

    /// Lower CI quantile.
    #[arg(long, default_value_t = 0.025)]
    ci_low: f64,

    /// Upper CI quantile.
    #[arg(long, default_value_t = 0.975)]
    ci_high: f64,
}

#[derive(ClapArgs, Debug)]
pub struct ProgramArgs {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Signed loading TSV with columns `program`, `gene_symbol`, and `loading`.
    #[arg(long)]
    loadings: PathBuf,

    /// Output TSV path.
    #[arg(long)]
    output: PathBuf,

    /// Comma-separated comparisons in A-B form.
    #[arg(long)]
    groups: String,

    /// Bootstrap mode: welch-t for unpaired or paired-t for matched subjects.
    #[arg(long, default_value = "welch-t")]
    test: String,

    /// Number of bootstrap replicates.
    #[arg(long, default_value_t = 2000)]
    n: usize,

    /// Minimum subjects per group, or matched subjects in paired mode.
    #[arg(long, default_value_t = 2)]
    min_pairs: usize,

    /// Deterministic RNG seed.
    #[arg(long, default_value_t = 20260418)]
    seed: u64,

    /// Lower CI quantile.
    #[arg(long, default_value_t = 0.025)]
    ci_low: f64,

    /// Upper CI quantile.
    #[arg(long, default_value_t = 0.975)]
    ci_high: f64,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ProteinKey {
    platform: String,
    panel: String,
    assay_id: String,
    gene_symbol: String,
}

#[derive(Debug)]
struct ProteinMeta {
    uniprot: String,
}

#[derive(Debug)]
struct BootstrapRow {
    comparison: String,
    platform: String,
    panel: String,
    assay_id: String,
    gene_symbol: String,
    uniprot: String,
    n_a: usize,
    n_b: usize,
    n_pairs: usize,
    point_effect: Option<f64>,
    bootstrap_mean: Option<f64>,
    ci_low: Option<f64>,
    ci_high: Option<f64>,
    sign_stability: Option<f64>,
    n_bootstrap: usize,
    skip_reason: String,
}

#[derive(Debug)]
struct ModuleBootstrapRow {
    comparison: String,
    module: String,
    n_a: usize,
    n_b: usize,
    n_pairs: usize,
    n_genes_declared: usize,
    n_genes_observed: usize,
    point_effect: Option<f64>,
    bootstrap_mean: Option<f64>,
    ci_low: Option<f64>,
    ci_high: Option<f64>,
    sign_stability: Option<f64>,
    n_bootstrap: usize,
    skip_reason: String,
}

#[derive(Debug)]
struct ProgramBootstrapRow {
    comparison: String,
    program: String,
    n_a: usize,
    n_b: usize,
    n_pairs: usize,
    n_loadings_declared: usize,
    n_loadings_observed: usize,
    point_effect: Option<f64>,
    bootstrap_mean: Option<f64>,
    ci_low: Option<f64>,
    ci_high: Option<f64>,
    sign_stability: Option<f64>,
    n_bootstrap: usize,
    skip_reason: String,
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Protein(args) => run_protein(args),
        Command::Module(args) => run_module(args),
        Command::Program(args) => run_program(args),
    }
}

fn run_protein(args: ProteinArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if args.test != "welch-t" && args.test != "paired-t" {
        bail!("bootstrap protein supports --test welch-t or paired-t");
    }
    if args.n == 0 {
        bail!("--n must be > 0");
    }
    if args.min_pairs < 2 {
        bail!("--min-pairs must be >= 2");
    }
    if !(0.0..=1.0).contains(&args.ci_low)
        || !(0.0..=1.0).contains(&args.ci_high)
        || args.ci_low >= args.ci_high
    {
        bail!("CI quantiles must satisfy 0 <= ci-low < ci-high <= 1");
    }

    let comparisons = parse_comparisons(&args.groups)?;
    let measurements_path = if args.input_dir.join("measurements.tsv").exists() {
        args.input_dir.join("measurements.tsv")
    } else {
        args.input_dir.join("measurements.tsv")
    };
    let measurements = read_measurements_long(&measurements_path)?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let proteins = read_proteins(&args.input_dir.join("proteins.tsv"))?;
    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    let mut meta: HashMap<(String, String), ProteinMeta> = HashMap::new();
    for p in proteins {
        meta.insert(
            (p.platform.as_str().to_string(), p.assay_id.0),
            ProteinMeta {
                uniprot: p.uniprot.join(","),
            },
        );
    }

    let mut cells: BTreeMap<ProteinKey, BTreeMap<String, BTreeMap<String, Vec<f64>>>> =
        BTreeMap::new();
    for m in &measurements {
        let Some(value) = m.effective_abundance() else {
            continue;
        };
        let Some(sample) = sample_by_id.get(m.sample_id.as_str()) else {
            continue;
        };
        if sample.is_control {
            continue;
        }
        let Some(condition) = &sample.condition else {
            continue;
        };
        let subject = sample.subject_id.as_ref().unwrap_or(&sample.sample_id);
        let key = ProteinKey {
            platform: m.platform.as_str().to_string(),
            panel: m.panel.clone().unwrap_or_default(),
            assay_id: m.assay_id.0.clone(),
            gene_symbol: m.gene_symbol.clone().unwrap_or_default(),
        };
        cells
            .entry(key)
            .or_default()
            .entry(condition.clone())
            .or_default()
            .entry(subject.clone())
            .or_default()
            .push(value);
    }

    let mut rng = Rng64::new(args.seed);
    let mut rows = Vec::new();
    for (a, b) in &comparisons {
        let comparison = format!("{a}-{b}");
        for (key, by_condition) in &cells {
            let a_values = subject_means(by_condition.get(a));
            let b_values = subject_means(by_condition.get(b));
            let uniprot = meta
                .get(&(key.platform.clone(), key.assay_id.clone()))
                .map(|m| m.uniprot.clone())
                .unwrap_or_default();
            let row = if args.test == "paired-t" {
                bootstrap_paired(
                    &args,
                    &mut rng,
                    &comparison,
                    key,
                    &uniprot,
                    &a_values,
                    &b_values,
                )
            } else {
                bootstrap_unpaired(
                    &args,
                    &mut rng,
                    &comparison,
                    key,
                    &uniprot,
                    &a_values,
                    &b_values,
                )
            };
            rows.push(row);
        }
    }

    rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| a.panel.cmp(&b.panel))
            .then_with(|| a.gene_symbol.cmp(&b.gene_symbol))
            .then_with(|| a.assay_id.cmp(&b.assay_id))
    });
    write_rows(&args.output, &rows)?;
    let computed = rows.iter().filter(|r| r.point_effect.is_some()).count();
    eprintln!(
        "bootstrap protein: test={} comparisons={} rows={} computed={} n={} seed={}",
        args.test,
        comparisons.len(),
        rows.len(),
        computed,
        args.n,
        args.seed
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &[
            "measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
    )?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "bootstrap protein",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output": args.output.display().to_string(),
            "groups": args.groups,
            "test": args.test,
            "n": args.n,
            "min-pairs": args.min_pairs,
            "seed": args.seed,
            "ci-low": args.ci_low,
            "ci-high": args.ci_high,
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("bootstrap protein: sidecar={}", sidecar.display());
    Ok(())
}

fn run_module(args: ModuleArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if args.test != "welch-t" && args.test != "paired-t" {
        bail!("bootstrap module supports --test welch-t or paired-t");
    }
    if args.n == 0 {
        bail!("--n must be > 0");
    }
    if args.min_pairs < 2 {
        bail!("--min-pairs must be >= 2");
    }
    if !(0.0..=1.0).contains(&args.ci_low)
        || !(0.0..=1.0).contains(&args.ci_high)
        || args.ci_low >= args.ci_high
    {
        bail!("CI quantiles must satisfy 0 <= ci-low < ci-high <= 1");
    }

    let comparisons = parse_comparisons(&args.groups)?;
    let modules = read_modules(&args.modules_tsv)?;
    let measurements_path = if args.input_dir.join("measurements.tsv").exists() {
        args.input_dir.join("measurements.tsv")
    } else {
        args.input_dir.join("measurements.tsv")
    };
    let measurements = read_measurements_long(&measurements_path)?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    let mut gene_to_modules: HashMap<&str, Vec<&str>> = HashMap::new();
    for (module, genes) in &modules {
        for gene in genes {
            gene_to_modules
                .entry(gene.as_str())
                .or_default()
                .push(module.as_str());
        }
    }

    #[derive(Default)]
    struct Acc {
        sum: f64,
        n: usize,
    }

    let mut observed_genes: BTreeMap<String, BTreeMap<String, ()>> = BTreeMap::new();
    let mut per_module_sample: HashMap<(String, String), Acc> = HashMap::new();
    for m in &measurements {
        let Some(gene) = m.gene_symbol.as_deref() else {
            continue;
        };
        let Some(module_names) = gene_to_modules.get(gene) else {
            continue;
        };
        let Some(value) = m.effective_abundance() else {
            continue;
        };
        for module in module_names {
            observed_genes
                .entry((*module).to_string())
                .or_default()
                .insert(gene.to_string(), ());
            let acc = per_module_sample
                .entry(((*module).to_string(), m.sample_id.clone()))
                .or_default();
            acc.sum += value;
            acc.n += 1;
        }
    }

    let mut cells: BTreeMap<String, BTreeMap<String, BTreeMap<String, Vec<f64>>>> = BTreeMap::new();
    for ((module, sample_id), acc) in per_module_sample {
        if acc.n == 0 {
            continue;
        }
        let Some(sample) = sample_by_id.get(sample_id.as_str()) else {
            continue;
        };
        if sample.is_control {
            continue;
        }
        let Some(condition) = &sample.condition else {
            continue;
        };
        let subject = sample.subject_id.as_ref().unwrap_or(&sample.sample_id);
        cells
            .entry(module)
            .or_default()
            .entry(condition.clone())
            .or_default()
            .entry(subject.clone())
            .or_default()
            .push(acc.sum / acc.n as f64);
    }

    let mut rng = Rng64::new(args.seed);
    let mut rows = Vec::new();
    for (a, b) in &comparisons {
        let comparison = format!("{a}-{b}");
        for (module, genes) in &modules {
            let by_condition = cells.get(module);
            let a_values = subject_means(by_condition.and_then(|m| m.get(a)));
            let b_values = subject_means(by_condition.and_then(|m| m.get(b)));
            let n_observed = observed_genes.get(module).map(BTreeMap::len).unwrap_or(0);
            let ctx = ModuleRowContext {
                comparison: &comparison,
                module,
                n_a: a_values.len(),
                n_b: b_values.len(),
                n_pairs: 0,
                n_genes_declared: genes.len(),
                n_genes_observed: n_observed,
            };
            let row = if args.test == "paired-t" {
                bootstrap_module_paired(&args, &mut rng, ctx, &a_values, &b_values)
            } else {
                bootstrap_module_unpaired(&args, &mut rng, ctx, &a_values, &b_values)
            };
            rows.push(row);
        }
    }

    rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| a.module.cmp(&b.module))
    });
    write_module_rows(&args.output, &rows)?;
    let computed = rows.iter().filter(|r| r.point_effect.is_some()).count();
    eprintln!(
        "bootstrap module: test={} comparisons={} rows={} computed={} n={} seed={}",
        args.test,
        comparisons.len(),
        rows.len(),
        computed,
        args.n,
        args.seed
    );

    let finished_at = SystemTime::now();
    let qc_path = args.input_dir.join("measurements.tsv");
    let raw_path = args.input_dir.join("measurements.tsv");
    let samples_path = args.input_dir.join("samples.tsv");
    let proteins_path = args.input_dir.join("proteins.tsv");
    let inputs_sha256 = hash_labeled_inputs(&[
        ("measurements.tsv", qc_path.as_path()),
        ("measurements.tsv", raw_path.as_path()),
        ("samples.tsv", samples_path.as_path()),
        ("proteins.tsv", proteins_path.as_path()),
        ("modules_tsv", args.modules_tsv.as_path()),
    ])?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "bootstrap module",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "modules-tsv": args.modules_tsv.display().to_string(),
            "output": args.output.display().to_string(),
            "groups": args.groups,
            "test": args.test,
            "n": args.n,
            "min-pairs": args.min_pairs,
            "seed": args.seed,
            "ci-low": args.ci_low,
            "ci-high": args.ci_high,
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("bootstrap module: sidecar={}", sidecar.display());
    Ok(())
}

fn run_program(args: ProgramArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if args.test != "welch-t" && args.test != "paired-t" {
        bail!("bootstrap program supports --test welch-t or paired-t");
    }
    if args.n == 0 {
        bail!("--n must be > 0");
    }
    if args.min_pairs < 2 {
        bail!("--min-pairs must be >= 2");
    }
    if !(0.0..=1.0).contains(&args.ci_low)
        || !(0.0..=1.0).contains(&args.ci_high)
        || args.ci_low >= args.ci_high
    {
        bail!("CI quantiles must satisfy 0 <= ci-low < ci-high <= 1");
    }

    let comparisons = parse_comparisons(&args.groups)?;
    let loadings = read_program_loadings(&args.loadings)?;
    let measurements_path = if args.input_dir.join("measurements.tsv").exists() {
        args.input_dir.join("measurements.tsv")
    } else {
        args.input_dir.join("measurements.tsv")
    };
    let measurements = read_measurements_long(&measurements_path)?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    let mut gene_to_programs: HashMap<&str, Vec<(&str, f64)>> = HashMap::new();
    for (program, genes) in &loadings {
        for (gene, loading) in genes {
            gene_to_programs
                .entry(gene.as_str())
                .or_default()
                .push((program.as_str(), *loading));
        }
    }

    #[derive(Default)]
    struct WeightedAcc {
        weighted_sum: f64,
        weight_abs_sum: f64,
        n: usize,
    }

    let mut observed_genes: BTreeMap<String, BTreeMap<String, ()>> = BTreeMap::new();
    let mut per_program_sample: HashMap<(String, String), WeightedAcc> = HashMap::new();
    for m in &measurements {
        let Some(gene) = m.gene_symbol.as_deref() else {
            continue;
        };
        let Some(program_names) = gene_to_programs.get(gene) else {
            continue;
        };
        let Some(value) = m.effective_abundance() else {
            continue;
        };
        for (program, loading) in program_names {
            observed_genes
                .entry((*program).to_string())
                .or_default()
                .insert(gene.to_string(), ());
            let acc = per_program_sample
                .entry(((*program).to_string(), m.sample_id.clone()))
                .or_default();
            acc.weighted_sum += value * loading;
            acc.weight_abs_sum += loading.abs();
            acc.n += 1;
        }
    }

    let mut cells: BTreeMap<String, BTreeMap<String, BTreeMap<String, Vec<f64>>>> = BTreeMap::new();
    for ((program, sample_id), acc) in per_program_sample {
        if acc.n == 0 || acc.weight_abs_sum == 0.0 {
            continue;
        }
        let Some(sample) = sample_by_id.get(sample_id.as_str()) else {
            continue;
        };
        if sample.is_control {
            continue;
        }
        let Some(condition) = &sample.condition else {
            continue;
        };
        let subject = sample.subject_id.as_ref().unwrap_or(&sample.sample_id);
        cells
            .entry(program)
            .or_default()
            .entry(condition.clone())
            .or_default()
            .entry(subject.clone())
            .or_default()
            .push(acc.weighted_sum / acc.weight_abs_sum);
    }

    let mut rng = Rng64::new(args.seed);
    let mut rows = Vec::new();
    for (a, b) in &comparisons {
        let comparison = format!("{a}-{b}");
        for (program, genes) in &loadings {
            let by_condition = cells.get(program);
            let a_values = subject_means(by_condition.and_then(|m| m.get(a)));
            let b_values = subject_means(by_condition.and_then(|m| m.get(b)));
            let n_observed = observed_genes.get(program).map(BTreeMap::len).unwrap_or(0);
            let mut ctx = ProgramRowContext {
                comparison: &comparison,
                program,
                n_a: a_values.len(),
                n_b: b_values.len(),
                n_pairs: 0,
                n_loadings_declared: genes.len(),
                n_loadings_observed: n_observed,
            };
            let row = if args.test == "paired-t" {
                let diffs: Vec<f64> = a_values
                    .iter()
                    .filter_map(|(subject, a)| b_values.get(subject).map(|b| a - b))
                    .collect();
                ctx.n_pairs = diffs.len();
                bootstrap_program_from_values(&args, &mut rng, ctx, None, Some(&diffs))
            } else {
                bootstrap_program_from_values(
                    &args,
                    &mut rng,
                    ctx,
                    Some((&a_values, &b_values)),
                    None,
                )
            };
            rows.push(row);
        }
    }

    rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| a.program.cmp(&b.program))
    });
    write_program_rows(&args.output, &rows)?;
    let computed = rows.iter().filter(|r| r.point_effect.is_some()).count();
    eprintln!(
        "bootstrap program: test={} comparisons={} rows={} computed={} n={} seed={}",
        args.test,
        comparisons.len(),
        rows.len(),
        computed,
        args.n,
        args.seed
    );

    let finished_at = SystemTime::now();
    let qc_path = args.input_dir.join("measurements.tsv");
    let raw_path = args.input_dir.join("measurements.tsv");
    let samples_path = args.input_dir.join("samples.tsv");
    let proteins_path = args.input_dir.join("proteins.tsv");
    let inputs_sha256 = hash_labeled_inputs(&[
        ("measurements.tsv", qc_path.as_path()),
        ("measurements.tsv", raw_path.as_path()),
        ("samples.tsv", samples_path.as_path()),
        ("proteins.tsv", proteins_path.as_path()),
        ("loadings", args.loadings.as_path()),
    ])?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "bootstrap program",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "loadings": args.loadings.display().to_string(),
            "output": args.output.display().to_string(),
            "groups": args.groups,
            "test": args.test,
            "n": args.n,
            "min-pairs": args.min_pairs,
            "seed": args.seed,
            "ci-low": args.ci_low,
            "ci-high": args.ci_high,
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("bootstrap program: sidecar={}", sidecar.display());
    Ok(())
}

fn bootstrap_unpaired(
    args: &ProteinArgs,
    rng: &mut Rng64,
    comparison: &str,
    key: &ProteinKey,
    uniprot: &str,
    a_values: &BTreeMap<String, f64>,
    b_values: &BTreeMap<String, f64>,
) -> BootstrapRow {
    let a: Vec<f64> = a_values.values().copied().collect();
    let b: Vec<f64> = b_values.values().copied().collect();
    if a.len() < args.min_pairs || b.len() < args.min_pairs {
        return skipped_row(
            comparison,
            key,
            uniprot,
            a.len(),
            b.len(),
            0,
            "insufficient_subjects",
        );
    }
    let point = mean(&a) - mean(&b);
    let mut effects = Vec::with_capacity(args.n);
    for _ in 0..args.n {
        let ma = bootstrap_mean(&a, rng);
        let mb = bootstrap_mean(&b, rng);
        effects.push(ma - mb);
    }
    computed_row(
        RowContext {
            comparison,
            key,
            uniprot,
            n_a: a.len(),
            n_b: b.len(),
            n_pairs: 0,
        },
        point,
        effects,
        args,
    )
}

fn bootstrap_paired(
    args: &ProteinArgs,
    rng: &mut Rng64,
    comparison: &str,
    key: &ProteinKey,
    uniprot: &str,
    a_values: &BTreeMap<String, f64>,
    b_values: &BTreeMap<String, f64>,
) -> BootstrapRow {
    let diffs: Vec<f64> = a_values
        .iter()
        .filter_map(|(subject, a)| b_values.get(subject).map(|b| a - b))
        .collect();
    if diffs.len() < args.min_pairs {
        return skipped_row(
            comparison,
            key,
            uniprot,
            a_values.len(),
            b_values.len(),
            diffs.len(),
            "insufficient_pairs",
        );
    }
    let point = mean(&diffs);
    let mut effects = Vec::with_capacity(args.n);
    for _ in 0..args.n {
        effects.push(bootstrap_mean(&diffs, rng));
    }
    computed_row(
        RowContext {
            comparison,
            key,
            uniprot,
            n_a: a_values.len(),
            n_b: b_values.len(),
            n_pairs: diffs.len(),
        },
        point,
        effects,
        args,
    )
}

fn bootstrap_module_unpaired(
    args: &ModuleArgs,
    rng: &mut Rng64,
    ctx: ModuleRowContext<'_>,
    a_values: &BTreeMap<String, f64>,
    b_values: &BTreeMap<String, f64>,
) -> ModuleBootstrapRow {
    let a: Vec<f64> = a_values.values().copied().collect();
    let b: Vec<f64> = b_values.values().copied().collect();
    if ctx.n_genes_observed == 0 {
        return skipped_module_row(ctx, "no_observed_genes");
    }
    if a.len() < args.min_pairs || b.len() < args.min_pairs {
        return skipped_module_row(ctx, "insufficient_subjects");
    }
    let point = mean(&a) - mean(&b);
    let mut effects = Vec::with_capacity(args.n);
    for _ in 0..args.n {
        effects.push(bootstrap_mean(&a, rng) - bootstrap_mean(&b, rng));
    }
    computed_module_row(ctx, point, effects, args.ci_low, args.ci_high)
}

fn bootstrap_module_paired(
    args: &ModuleArgs,
    rng: &mut Rng64,
    mut ctx: ModuleRowContext<'_>,
    a_values: &BTreeMap<String, f64>,
    b_values: &BTreeMap<String, f64>,
) -> ModuleBootstrapRow {
    let diffs: Vec<f64> = a_values
        .iter()
        .filter_map(|(subject, a)| b_values.get(subject).map(|b| a - b))
        .collect();
    ctx.n_pairs = diffs.len();
    if ctx.n_genes_observed == 0 {
        return skipped_module_row(ctx, "no_observed_genes");
    }
    if diffs.len() < args.min_pairs {
        return skipped_module_row(ctx, "insufficient_pairs");
    }
    let point = mean(&diffs);
    let mut effects = Vec::with_capacity(args.n);
    for _ in 0..args.n {
        effects.push(bootstrap_mean(&diffs, rng));
    }
    computed_module_row(ctx, point, effects, args.ci_low, args.ci_high)
}

struct RowContext<'a> {
    comparison: &'a str,
    key: &'a ProteinKey,
    uniprot: &'a str,
    n_a: usize,
    n_b: usize,
    n_pairs: usize,
}

fn computed_row(
    ctx: RowContext<'_>,
    point: f64,
    mut effects: Vec<f64>,
    args: &ProteinArgs,
) -> BootstrapRow {
    let boot_mean = mean(&effects);
    let sign_stability = sign_stability(point, &effects);
    effects.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    BootstrapRow {
        comparison: ctx.comparison.to_string(),
        platform: ctx.key.platform.clone(),
        panel: ctx.key.panel.clone(),
        assay_id: ctx.key.assay_id.clone(),
        gene_symbol: ctx.key.gene_symbol.clone(),
        uniprot: ctx.uniprot.to_string(),
        n_a: ctx.n_a,
        n_b: ctx.n_b,
        n_pairs: ctx.n_pairs,
        point_effect: Some(point),
        bootstrap_mean: Some(boot_mean),
        ci_low: Some(quantile_sorted(&effects, args.ci_low)),
        ci_high: Some(quantile_sorted(&effects, args.ci_high)),
        sign_stability: Some(sign_stability),
        n_bootstrap: effects.len(),
        skip_reason: String::new(),
    }
}

fn skipped_row(
    comparison: &str,
    key: &ProteinKey,
    uniprot: &str,
    n_a: usize,
    n_b: usize,
    n_pairs: usize,
    reason: &str,
) -> BootstrapRow {
    BootstrapRow {
        comparison: comparison.to_string(),
        platform: key.platform.clone(),
        panel: key.panel.clone(),
        assay_id: key.assay_id.clone(),
        gene_symbol: key.gene_symbol.clone(),
        uniprot: uniprot.to_string(),
        n_a,
        n_b,
        n_pairs,
        point_effect: None,
        bootstrap_mean: None,
        ci_low: None,
        ci_high: None,
        sign_stability: None,
        n_bootstrap: 0,
        skip_reason: reason.to_string(),
    }
}

#[derive(Clone, Copy)]
struct ModuleRowContext<'a> {
    comparison: &'a str,
    module: &'a str,
    n_a: usize,
    n_b: usize,
    n_pairs: usize,
    n_genes_declared: usize,
    n_genes_observed: usize,
}

#[derive(Clone, Copy)]
struct ProgramRowContext<'a> {
    comparison: &'a str,
    program: &'a str,
    n_a: usize,
    n_b: usize,
    n_pairs: usize,
    n_loadings_declared: usize,
    n_loadings_observed: usize,
}

type SubjectValues<'a> = (&'a BTreeMap<String, f64>, &'a BTreeMap<String, f64>);

fn computed_module_row(
    ctx: ModuleRowContext<'_>,
    point: f64,
    mut effects: Vec<f64>,
    ci_low: f64,
    ci_high: f64,
) -> ModuleBootstrapRow {
    let boot_mean = mean(&effects);
    let sign_stability = sign_stability(point, &effects);
    effects.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    ModuleBootstrapRow {
        comparison: ctx.comparison.to_string(),
        module: ctx.module.to_string(),
        n_a: ctx.n_a,
        n_b: ctx.n_b,
        n_pairs: ctx.n_pairs,
        n_genes_declared: ctx.n_genes_declared,
        n_genes_observed: ctx.n_genes_observed,
        point_effect: Some(point),
        bootstrap_mean: Some(boot_mean),
        ci_low: Some(quantile_sorted(&effects, ci_low)),
        ci_high: Some(quantile_sorted(&effects, ci_high)),
        sign_stability: Some(sign_stability),
        n_bootstrap: effects.len(),
        skip_reason: String::new(),
    }
}

fn skipped_module_row(ctx: ModuleRowContext<'_>, reason: &str) -> ModuleBootstrapRow {
    ModuleBootstrapRow {
        comparison: ctx.comparison.to_string(),
        module: ctx.module.to_string(),
        n_a: ctx.n_a,
        n_b: ctx.n_b,
        n_pairs: ctx.n_pairs,
        n_genes_declared: ctx.n_genes_declared,
        n_genes_observed: ctx.n_genes_observed,
        point_effect: None,
        bootstrap_mean: None,
        ci_low: None,
        ci_high: None,
        sign_stability: None,
        n_bootstrap: 0,
        skip_reason: reason.to_string(),
    }
}

fn bootstrap_program_from_values(
    args: &ProgramArgs,
    rng: &mut Rng64,
    ctx: ProgramRowContext<'_>,
    unpaired: Option<SubjectValues<'_>>,
    paired_diffs: Option<&[f64]>,
) -> ProgramBootstrapRow {
    if ctx.n_loadings_observed == 0 {
        return skipped_program_row(ctx, "no_observed_loadings");
    }
    let mut effects = Vec::with_capacity(args.n);
    let point = if let Some(diffs) = paired_diffs {
        if diffs.len() < args.min_pairs {
            return skipped_program_row(ctx, "insufficient_pairs");
        }
        for _ in 0..args.n {
            effects.push(bootstrap_mean(diffs, rng));
        }
        mean(diffs)
    } else if let Some((a_values, b_values)) = unpaired {
        let a: Vec<f64> = a_values.values().copied().collect();
        let b: Vec<f64> = b_values.values().copied().collect();
        if a.len() < args.min_pairs || b.len() < args.min_pairs {
            return skipped_program_row(ctx, "insufficient_subjects");
        }
        for _ in 0..args.n {
            effects.push(bootstrap_mean(&a, rng) - bootstrap_mean(&b, rng));
        }
        mean(&a) - mean(&b)
    } else {
        return skipped_program_row(ctx, "internal_error");
    };
    computed_program_row(ctx, point, effects, args.ci_low, args.ci_high)
}

fn computed_program_row(
    ctx: ProgramRowContext<'_>,
    point: f64,
    mut effects: Vec<f64>,
    ci_low: f64,
    ci_high: f64,
) -> ProgramBootstrapRow {
    let boot_mean = mean(&effects);
    let sign_stability = sign_stability(point, &effects);
    effects.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    ProgramBootstrapRow {
        comparison: ctx.comparison.to_string(),
        program: ctx.program.to_string(),
        n_a: ctx.n_a,
        n_b: ctx.n_b,
        n_pairs: ctx.n_pairs,
        n_loadings_declared: ctx.n_loadings_declared,
        n_loadings_observed: ctx.n_loadings_observed,
        point_effect: Some(point),
        bootstrap_mean: Some(boot_mean),
        ci_low: Some(quantile_sorted(&effects, ci_low)),
        ci_high: Some(quantile_sorted(&effects, ci_high)),
        sign_stability: Some(sign_stability),
        n_bootstrap: effects.len(),
        skip_reason: String::new(),
    }
}

fn skipped_program_row(ctx: ProgramRowContext<'_>, reason: &str) -> ProgramBootstrapRow {
    ProgramBootstrapRow {
        comparison: ctx.comparison.to_string(),
        program: ctx.program.to_string(),
        n_a: ctx.n_a,
        n_b: ctx.n_b,
        n_pairs: ctx.n_pairs,
        n_loadings_declared: ctx.n_loadings_declared,
        n_loadings_observed: ctx.n_loadings_observed,
        point_effect: None,
        bootstrap_mean: None,
        ci_low: None,
        ci_high: None,
        sign_stability: None,
        n_bootstrap: 0,
        skip_reason: reason.to_string(),
    }
}

fn subject_means(input: Option<&BTreeMap<String, Vec<f64>>>) -> BTreeMap<String, f64> {
    input
        .into_iter()
        .flat_map(|m| m.iter())
        .map(|(subject, values)| (subject.clone(), mean(values)))
        .collect()
}

fn bootstrap_mean(values: &[f64], rng: &mut Rng64) -> f64 {
    let mut total = 0.0;
    for _ in 0..values.len() {
        total += values[rng.gen_range(values.len())];
    }
    total / values.len() as f64
}

fn quantile_sorted(values: &[f64], q: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let idx = ((values.len() - 1) as f64 * q).round() as usize;
    values[idx.min(values.len() - 1)]
}

fn sign_stability(point: f64, effects: &[f64]) -> f64 {
    if point > 0.0 {
        effects.iter().filter(|v| **v > 0.0).count() as f64 / effects.len() as f64
    } else if point < 0.0 {
        effects.iter().filter(|v| **v < 0.0).count() as f64 / effects.len() as f64
    } else {
        effects.iter().filter(|v| **v == 0.0).count() as f64 / effects.len() as f64
    }
}

fn write_rows(path: &Path, rows: &[BootstrapRow]) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output dir {:?}", parent))?;
    }
    let mut buf = String::from(
        "comparison\tplatform\tpanel\tassay_id\tgene_symbol\tuniprot\tn_a\tn_b\tn_pairs\tpoint_effect\tbootstrap_mean\tci_low\tci_high\tsign_stability\tn_bootstrap\tskip_reason\n",
    );
    for r in rows {
        buf.push_str(&r.comparison);
        buf.push('\t');
        buf.push_str(&r.platform);
        buf.push('\t');
        buf.push_str(&r.panel);
        buf.push('\t');
        buf.push_str(&r.assay_id);
        buf.push('\t');
        buf.push_str(&r.gene_symbol);
        buf.push('\t');
        buf.push_str(&r.uniprot);
        buf.push('\t');
        buf.push_str(&r.n_a.to_string());
        buf.push('\t');
        buf.push_str(&r.n_b.to_string());
        buf.push('\t');
        buf.push_str(&r.n_pairs.to_string());
        buf.push('\t');
        push_opt(&mut buf, r.point_effect);
        buf.push('\t');
        push_opt(&mut buf, r.bootstrap_mean);
        buf.push('\t');
        push_opt(&mut buf, r.ci_low);
        buf.push('\t');
        push_opt(&mut buf, r.ci_high);
        buf.push('\t');
        push_opt(&mut buf, r.sign_stability);
        buf.push('\t');
        buf.push_str(&r.n_bootstrap.to_string());
        buf.push('\t');
        buf.push_str(&r.skip_reason);
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn write_module_rows(path: &Path, rows: &[ModuleBootstrapRow]) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output dir {:?}", parent))?;
    }
    let mut buf = String::from(
        "comparison\tmodule\tn_a\tn_b\tn_pairs\tn_genes_declared\tn_genes_observed\tpoint_effect\tbootstrap_mean\tci_low\tci_high\tsign_stability\tn_bootstrap\tskip_reason\n",
    );
    for r in rows {
        buf.push_str(&r.comparison);
        buf.push('\t');
        buf.push_str(&r.module);
        buf.push('\t');
        buf.push_str(&r.n_a.to_string());
        buf.push('\t');
        buf.push_str(&r.n_b.to_string());
        buf.push('\t');
        buf.push_str(&r.n_pairs.to_string());
        buf.push('\t');
        buf.push_str(&r.n_genes_declared.to_string());
        buf.push('\t');
        buf.push_str(&r.n_genes_observed.to_string());
        buf.push('\t');
        push_opt(&mut buf, r.point_effect);
        buf.push('\t');
        push_opt(&mut buf, r.bootstrap_mean);
        buf.push('\t');
        push_opt(&mut buf, r.ci_low);
        buf.push('\t');
        push_opt(&mut buf, r.ci_high);
        buf.push('\t');
        push_opt(&mut buf, r.sign_stability);
        buf.push('\t');
        buf.push_str(&r.n_bootstrap.to_string());
        buf.push('\t');
        buf.push_str(&r.skip_reason);
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn write_program_rows(path: &Path, rows: &[ProgramBootstrapRow]) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output dir {:?}", parent))?;
    }
    let mut buf = String::from(
        "comparison\tprogram\tn_a\tn_b\tn_pairs\tn_loadings_declared\tn_loadings_observed\tpoint_effect\tbootstrap_mean\tci_low\tci_high\tsign_stability\tn_bootstrap\tskip_reason\n",
    );
    for r in rows {
        buf.push_str(&r.comparison);
        buf.push('\t');
        buf.push_str(&r.program);
        buf.push('\t');
        buf.push_str(&r.n_a.to_string());
        buf.push('\t');
        buf.push_str(&r.n_b.to_string());
        buf.push('\t');
        buf.push_str(&r.n_pairs.to_string());
        buf.push('\t');
        buf.push_str(&r.n_loadings_declared.to_string());
        buf.push('\t');
        buf.push_str(&r.n_loadings_observed.to_string());
        buf.push('\t');
        push_opt(&mut buf, r.point_effect);
        buf.push('\t');
        push_opt(&mut buf, r.bootstrap_mean);
        buf.push('\t');
        push_opt(&mut buf, r.ci_low);
        buf.push('\t');
        push_opt(&mut buf, r.ci_high);
        buf.push('\t');
        push_opt(&mut buf, r.sign_stability);
        buf.push('\t');
        buf.push_str(&r.n_bootstrap.to_string());
        buf.push('\t');
        buf.push_str(&r.skip_reason);
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn read_modules(path: &Path) -> Result<BTreeMap<String, Vec<String>>> {
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
    let mut out: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let module = row[module_col].trim();
        let gene = row[gene_col].trim();
        if module.is_empty() || gene.is_empty() {
            continue;
        }
        out.entry(module.to_string())
            .or_default()
            .push(gene.to_string());
    }
    Ok(out)
}

fn read_program_loadings(path: &Path) -> Result<BTreeMap<String, Vec<(String, f64)>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let program_col = headers
        .iter()
        .position(|h| h == "program")
        .with_context(|| format!("missing column `program` in {:?}", path))?;
    let gene_col = headers
        .iter()
        .position(|h| h == "gene_symbol")
        .with_context(|| format!("missing column `gene_symbol` in {:?}", path))?;
    let loading_col = headers
        .iter()
        .position(|h| h == "loading")
        .with_context(|| format!("missing column `loading` in {:?}", path))?;
    let mut out: BTreeMap<String, Vec<(String, f64)>> = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let program = row[program_col].trim();
        let gene = row[gene_col].trim();
        let loading: f64 = row[loading_col]
            .trim()
            .parse()
            .with_context(|| format!("parsing loading in {:?}", path))?;
        if program.is_empty() || gene.is_empty() || !loading.is_finite() || loading == 0.0 {
            continue;
        }
        out.entry(program.to_string())
            .or_default()
            .push((gene.to_string(), loading));
    }
    if out.is_empty() {
        bail!("no finite non-zero program loadings found in {:?}", path);
    }
    Ok(out)
}

fn push_opt(buf: &mut String, value: Option<f64>) {
    if let Some(value) = value {
        buf.push_str(&format!("{value}"));
    }
}

struct Rng64 {
    state: u64,
}

impl Rng64 {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 0x9e3779b97f4a7c15 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545f4914f6cdd1d)
    }

    fn gen_range(&mut self, upper: usize) -> usize {
        (self.next_u64() as usize) % upper
    }
}
