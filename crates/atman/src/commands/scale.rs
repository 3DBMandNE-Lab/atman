//! `atman scale absolute`: put a canonical log2 matrix on an absolute
//! per-sample scale for Reiber-style analyses.
//!
//! For every sample with a total-protein value in `samples.tsv`:
//!
//! ```text
//! a' = a − log2( Σ_assays 2^a ) + log2(total_protein)
//! ```
//!
//! where `a` is the canonical log2 abundance after QC masks and the sum runs
//! over that sample's measured assays. Samples without a total-protein value
//! are dropped from the output measurements (their count is recorded in the
//! sidecar); `samples.tsv` and `proteins.tsv` are copied verbatim.

use anyhow::{bail, Context, Result};
use atman_core::Abundance;
use clap::{Args as ClapArgs, Subcommand};
use serde_json::json;
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use std::time::SystemTime;

use crate::design::CovariateFrame;
use crate::io::{
    format_f64, hash_canonical_inputs, read_measurements_long, sidecar_path_for,
    write_measurements_long, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Rescale each sample to an absolute scale using a total-protein column.
    Absolute(AbsoluteArgs),
}

#[derive(ClapArgs, Debug)]
pub struct AbsoluteArgs {
    /// Canonical input directory (log2 abundances).
    #[arg(long)]
    input_dir: PathBuf,

    /// samples.tsv column holding the per-sample total protein (linear scale).
    #[arg(long, default_value = "total_protein")]
    total_col: String,

    /// Output canonical directory (`measurements.tsv` rescaled; `samples.tsv`
    /// and `proteins.tsv` copied).
    #[arg(long)]
    output_canonical_dir: PathBuf,
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Absolute(a) => run_absolute(a),
    }
}

fn with_value(abundance: &Abundance, value: f64) -> Abundance {
    match abundance {
        Abundance::Log2Npx(_) => Abundance::Log2Npx(value),
        Abundance::Log2Intensity(_) => Abundance::Log2Intensity(value),
        Abundance::Ibaq(_) => Abundance::Ibaq(value),
        Abundance::Raw(_) => Abundance::Raw(value),
    }
}

fn run_absolute(args: AbsoluteArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let samples_path = args.input_dir.join("samples.tsv");
    let mut frame = CovariateFrame::new();
    frame.load_tsv(&samples_path)?;
    if !frame.has_column(&args.total_col) {
        bail!(
            "--total-col {:?} not found in {:?}",
            args.total_col,
            samples_path
        );
    }
    let measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;

    // Per-sample linear sum over measured assays.
    let mut linear_sum: BTreeMap<String, f64> = BTreeMap::new();
    for m in &measurements {
        if let Some(a) = m.effective_abundance() {
            *linear_sum.entry(m.sample_id.clone()).or_default() += a.exp2();
        }
    }
    let mut total_of: HashMap<String, f64> = HashMap::new();
    let mut n_samples_dropped = 0usize;
    for sample in linear_sum.keys() {
        match frame
            .get(sample, &args.total_col)
            .and_then(|v| v.parse::<f64>().ok())
            .filter(|v| v.is_finite() && *v > 0.0)
        {
            Some(t) => {
                total_of.insert(sample.clone(), t);
            }
            None => n_samples_dropped += 1,
        }
    }
    if total_of.is_empty() {
        bail!(
            "no sample has a positive {:?} value in {:?}",
            args.total_col,
            samples_path
        );
    }

    let mut records = Vec::with_capacity(measurements.len());
    let mut n_cells = 0usize;
    for m in &measurements {
        let Some(&total) = total_of.get(&m.sample_id) else {
            continue;
        };
        let Some(a) = m.effective_abundance() else {
            continue;
        };
        let scaled = a - linear_sum[&m.sample_id].log2() + total.log2();
        let mut r = m.clone();
        r.abundance = with_value(&m.abundance, scaled);
        r.abundance_raw = with_value(&m.abundance_raw, scaled);
        r.npx_source_str = format_f64(scaled);
        records.push(r);
        n_cells += 1;
    }
    std::fs::create_dir_all(&args.output_canonical_dir)
        .with_context(|| format!("creating {:?}", args.output_canonical_dir))?;
    let measurements_out = args.output_canonical_dir.join("measurements.tsv");
    write_measurements_long(&measurements_out, &records)?;
    let mut outputs = vec![measurements_out.clone()];
    for name in ["samples.tsv", "proteins.tsv"] {
        let src = args.input_dir.join(name);
        let dst = args.output_canonical_dir.join(name);
        std::fs::copy(&src, &dst).with_context(|| format!("copying {:?} to {:?}", src, dst))?;
        outputs.push(dst);
    }
    eprintln!(
        "scale absolute: samples={} dropped_without_total={} cells={} output={}",
        total_of.len(),
        n_samples_dropped,
        n_cells,
        args.output_canonical_dir.display()
    );

    let finished_at = SystemTime::now();
    let inputs = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    let sidecar = sidecar_path_for(&measurements_out);
    let mut extras = serde_json::Map::new();
    extras.insert("n_samples_scaled".into(), json!(total_of.len()));
    extras.insert(
        "n_samples_dropped_without_total".into(),
        json!(n_samples_dropped),
    );
    extras.insert("n_cells".into(), json!(n_cells));
    write_run_sidecar(
        &sidecar,
        "scale absolute",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "total-col": args.total_col,
            "output-canonical-dir": args.output_canonical_dir.display().to_string(),
        }),
        &inputs,
        &outputs,
        started_at,
        finished_at,
        Some(extras),
    )?;
    eprintln!("scale absolute: sidecar={}", sidecar.display());
    Ok(())
}
