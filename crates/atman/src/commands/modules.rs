//! `atman modules` — data-driven module discovery surface.
//!
//! Currently a single subcommand `discover`. Produces a canonical
//! `modules_discovered.tsv` whose schema matches the `modules.tsv`
//! contract consumed by `atman score modules` / `bootstrap module` /
//! `module-de` (`module, gene_symbol`) so the discovery output drops
//! straight into the existing scoring surface.

use anyhow::{bail, Context, Result};
use atman_core::modules_discover::{
    discover, DiscoveryMethod, ModuleAssignment, ModuleReportRow, Similarity,
    SoftPowerSweepRow,
};
use clap::{Args as ClapArgs, Subcommand};
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, hash_canonical_inputs, read_measurements_long, sidecar_path_for,
    write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// WGCNA-style soft-power adjacency + TOM + UPGMA clustering with a
    /// fixed height cut — produces a canonical `modules_discovered.tsv`
    /// keyed by `gene_symbol` that drops into `atman score modules`.
    Discover(DiscoverArgs),
}

#[derive(ClapArgs, Debug)]
pub struct DiscoverArgs {
    /// Canonical Atman directory (expects `qc_measurements.tsv`,
    /// `samples.tsv`, `proteins.tsv`).
    #[arg(long)]
    input_dir: PathBuf,

    /// Similarity metric used when building the adjacency: `pearson`
    /// (default) or `spearman`.
    #[arg(long, default_value = "pearson")]
    similarity: String,

    /// Adjacency construction: `wgcna-soft` (default, soft-power on |r|)
    /// or `hard-threshold` (binary adjacency at `|r| ≥ --threshold`).
    #[arg(long, default_value = "wgcna-soft")]
    method: String,

    /// Soft power β. `auto` picks the smallest integer in
    /// `{1..--max-beta}` whose scale-free topology `R² ≥ --r2-target`
    /// and slope < 0. Integer values set β explicitly.
    #[arg(long, default_value = "auto")]
    soft_power: String,

    /// Target scale-free `R²` for `--soft-power auto`.
    #[arg(long, default_value_t = 0.8)]
    r2_target: f64,

    /// Upper bound on the β sweep.
    #[arg(long, default_value_t = 20)]
    max_beta: usize,

    /// Number of log-bins in the scale-free topology histogram.
    #[arg(long, default_value_t = 10)]
    n_bins: usize,

    /// Hard-threshold value for `--method hard-threshold`.
    #[arg(long, default_value_t = 0.3)]
    threshold: f64,

    /// Minimum number of features per retained module. Features below
    /// threshold go to the `grey` catch-all.
    #[arg(long, default_value_t = 5)]
    min_module_size: usize,

    /// UPGMA tree-cut height. Clusters whose deepest merge ≤ this
    /// height form one module.
    #[arg(long, default_value_t = 0.5)]
    cut_height: f64,

    /// Canonical measurements file to read.
    #[arg(long, default_value = "qc")]
    source: String,

    /// Output directory (`modules_discovered.tsv`,
    /// `module_discovery_report.tsv`, and `soft_power_diagnostics.tsv`
    /// for wgcna-soft).
    #[arg(long)]
    output_dir: PathBuf,
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Discover(args) => run_discover(args),
    }
}

fn run_discover(args: DiscoverArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let similarity = match args.similarity.as_str() {
        "pearson" => Similarity::Pearson,
        "spearman" => Similarity::Spearman,
        other => bail!("--similarity {other:?}; expected pearson or spearman"),
    };
    let method = match args.method.as_str() {
        "wgcna-soft" => {
            let beta = if args.soft_power == "auto" {
                None
            } else {
                Some(
                    args.soft_power
                        .parse::<usize>()
                        .with_context(|| format!("--soft-power {:?}", args.soft_power))?,
                )
            };
            if let Some(b) = beta {
                if b < 1 {
                    bail!("--soft-power must be >= 1");
                }
            }
            if args.max_beta < 1 {
                bail!("--max-beta must be >= 1");
            }
            if !(0.0..=1.0).contains(&args.r2_target) {
                bail!("--r2-target must be in [0, 1]");
            }
            if args.n_bins < 3 {
                bail!("--n-bins must be >= 3");
            }
            DiscoveryMethod::WgcnaSoft {
                beta,
                r2_target: args.r2_target,
                n_bins: args.n_bins,
                max_beta: args.max_beta,
            }
        }
        "hard-threshold" => {
            if !(0.0..=1.0).contains(&args.threshold) {
                bail!("--threshold must be in [0, 1]");
            }
            DiscoveryMethod::HardThreshold {
                threshold: args.threshold,
            }
        }
        other => bail!("--method {other:?}; expected wgcna-soft or hard-threshold"),
    };
    if args.min_module_size < 2 {
        bail!("--min-module-size must be >= 2");
    }
    if !(0.0..=2.0).contains(&args.cut_height) {
        bail!("--cut-height must be in [0, 2]");
    }

    // Load the canonical measurements into a subject × feature matrix
    // keyed by gene_symbol.
    let source_file = match args.source.as_str() {
        "qc" => "qc_measurements.tsv",
        "raw" => "measurements.tsv",
        other => bail!("--source {other:?}; expected qc or raw"),
    };
    let records = read_measurements_long(&args.input_dir.join(source_file))?;
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
    if samples.len() < 3 {
        bail!(
            "need >= 3 samples to build a meaningful similarity graph; got {}",
            samples.len()
        );
    }
    if genes.len() < args.min_module_size {
        bail!(
            "only {} features retained; less than --min-module-size",
            genes.len()
        );
    }
    let sample_list: Vec<String> = samples.into_iter().collect();
    let feature_list: Vec<String> = genes.into_iter().collect();
    let mut data = vec![vec![f64::NAN; feature_list.len()]; sample_list.len()];
    for (si, sid) in sample_list.iter().enumerate() {
        for (gi, gene) in feature_list.iter().enumerate() {
            if let Some(&v) = cells.get(&(gene.clone(), sid.clone())) {
                data[si][gi] = v;
            }
        }
    }
    // Drop features with any missing values to keep the similarity
    // computation clean — WGCNA is not missingness-aware and neither
    // is v1 here. Impute is a deliberate follow-on.
    let keep_feature: Vec<bool> = (0..feature_list.len())
        .map(|gi| (0..sample_list.len()).all(|si| data[si][gi].is_finite()))
        .collect();
    let kept_features: Vec<String> = feature_list
        .iter()
        .zip(keep_feature.iter())
        .filter_map(|(f, k)| if *k { Some(f.clone()) } else { None })
        .collect();
    if kept_features.len() < args.min_module_size {
        bail!(
            "only {} fully-observed features remain; less than --min-module-size",
            kept_features.len()
        );
    }
    let filtered_data: Vec<Vec<f64>> = data
        .iter()
        .map(|row| {
            row.iter()
                .zip(keep_feature.iter())
                .filter_map(|(v, k)| if *k { Some(*v) } else { None })
                .collect()
        })
        .collect();

    let result = discover(
        &filtered_data,
        &kept_features,
        similarity,
        method,
        args.cut_height,
        args.min_module_size,
    )
    .map_err(|e| anyhow::anyhow!(e))?;

    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating {:?}", args.output_dir))?;

    let modules_path = args.output_dir.join("modules_discovered.tsv");
    write_modules(&modules_path, &result.modules)?;
    let report_path = args.output_dir.join("module_discovery_report.tsv");
    write_report(&report_path, &result.report)?;
    let mut outputs = vec![modules_path.clone(), report_path.clone()];
    if !result.soft_power_sweep.is_empty() {
        let diag_path = args.output_dir.join("soft_power_diagnostics.tsv");
        write_soft_power_diagnostics(&diag_path, &result.soft_power_sweep)?;
        outputs.push(diag_path);
    }

    eprintln!(
        "modules discover: features={} non_grey_modules={} chosen_beta={} method={:?}",
        kept_features.len(),
        result.report.iter().filter(|r| r.module != "grey").count(),
        result.soft_power_chosen,
        args.method,
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &[
            "qc_measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
    )?;
    let sidecar = sidecar_path_for(&modules_path);
    write_run_sidecar(
        &sidecar,
        "modules discover",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "similarity": args.similarity,
            "method": args.method,
            "soft-power": args.soft_power,
            "r2-target": args.r2_target,
            "max-beta": args.max_beta,
            "n-bins": args.n_bins,
            "threshold": args.threshold,
            "min-module-size": args.min_module_size,
            "cut-height": args.cut_height,
            "source": args.source,
            "chosen-beta": result.soft_power_chosen,
            "n-features-retained": kept_features.len(),
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        None,
)?;
    eprintln!("modules discover: sidecar={}", sidecar.display());
    Ok(())
}

fn write_modules(path: &Path, rows: &[ModuleAssignment]) -> Result<()> {
    let mut buf = String::from("module\tgene_symbol\n");
    for r in rows {
        buf.push_str(&r.module);
        buf.push('\t');
        buf.push_str(&r.feature);
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn write_report(path: &Path, rows: &[ModuleReportRow]) -> Result<()> {
    let mut buf = String::from(
        "module\tsize\tmean_within_abs_correlation\thub_feature\t\
         eigenprotein_pc1_variance_explained\n",
    );
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{:.6}\t{}\t{:.6}\n",
            r.module,
            r.size,
            r.mean_within_abs_correlation,
            r.hub_feature,
            r.eigenprotein_pc1_variance_explained,
        ));
    }
    atomic_write(path, buf.as_bytes())
}

fn write_soft_power_diagnostics(path: &Path, rows: &[SoftPowerSweepRow]) -> Result<()> {
    let mut buf = String::from("beta\tscale_free_r_squared\tslope\tmean_k\n");
    for r in rows {
        buf.push_str(&format!(
            "{}\t{:.6}\t{:.6}\t{:.6}\n",
            r.beta, r.scale_free_r_squared, r.slope, r.mean_k,
        ));
    }
    atomic_write(path, buf.as_bytes())
}
