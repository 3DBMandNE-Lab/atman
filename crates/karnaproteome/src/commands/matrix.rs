use anyhow::{Context, Result};
use clap::Args as ClapArgs;
use proteome_core::matrix::dube_wide_pivot;
use std::path::PathBuf;

use crate::io::{read_measurements_long, read_samples, write_dube_wide_panel};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing `qc_measurements.tsv` and `samples.tsv`.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for per-panel CSVs.
    #[arg(long)]
    output_dir: PathBuf,

    /// Output format. v0.1 supports only `dube-wide`.
    #[arg(long, default_value = "dube-wide")]
    format: String,

    /// Split dimension. v0.1 supports only `panel`.
    #[arg(long, default_value = "panel")]
    split_by: String,
}

pub fn run(args: Args) -> Result<()> {
    if args.format != "dube-wide" {
        anyhow::bail!("format {:?} not supported in v0.1", args.format);
    }
    if args.split_by != "panel" {
        anyhow::bail!("split-by {:?} not supported in v0.1", args.split_by);
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let measurements = read_measurements_long(&args.input_dir.join("qc_measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;

    let panels = dube_wide_pivot(&measurements, &samples);
    for panel in &panels {
        let p = write_dube_wide_panel(&args.output_dir, panel)?;
        eprintln!(
            "matrix: wrote {:?} ({} assays, {} rows)",
            p,
            panel.assays.len(),
            panel.rows.len()
        );
    }
    Ok(())
}
