use anyhow::{Context, Result};
use atman_core::matrix::pivot_wide_panels;
use clap::Args as ClapArgs;
use std::path::PathBuf;

use crate::io::{read_measurements_long, read_samples, write_wide_panel};

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
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;

    let panels = pivot_wide_panels(&measurements, &samples);
    for panel in &panels {
        let p = write_wide_panel(&args.output_dir, panel)?;
        eprintln!(
            "matrix: wrote {:?} ({} assays, {} rows)",
            p,
            panel.assays.len(),
            panel.rows.len()
        );
    }
    Ok(())
}
