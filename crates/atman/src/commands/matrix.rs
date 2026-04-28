use anyhow::{Context, Result};
use atman_core::matrix::pivot_wide_panels;
use clap::Args as ClapArgs;
use serde_json::json;
use std::path::PathBuf;
use std::time::SystemTime;

use crate::io::{
    hash_canonical_inputs, read_measurements_long, read_samples, sidecar_path_for,
    write_run_sidecar, write_wide_panel,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing `measurements.tsv` and `samples.tsv`.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for per-panel CSVs.
    #[arg(long)]
    output_dir: PathBuf,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;

    let panels = pivot_wide_panels(&measurements, &samples);
    let mut output_paths: Vec<PathBuf> = Vec::with_capacity(panels.len());
    for panel in &panels {
        let p = write_wide_panel(&args.output_dir, panel)?;
        eprintln!(
            "matrix: wrote {:?} ({} assays, {} rows)",
            p,
            panel.assays.len(),
            panel.rows.len()
        );
        output_paths.push(p);
    }

    let finished_at = SystemTime::now();
    let inputs_sha256 =
        hash_canonical_inputs(&args.input_dir, &["measurements.tsv", "samples.tsv"])?;
    let primary_output = output_paths
        .first()
        .cloned()
        .unwrap_or_else(|| args.output_dir.join("matrix.run"));
    let sidecar = sidecar_path_for(&primary_output);
    write_run_sidecar(
        &sidecar,
        "matrix",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
        }),
        &inputs_sha256,
        &output_paths,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("matrix: sidecar={}", sidecar.display());
    Ok(())
}
