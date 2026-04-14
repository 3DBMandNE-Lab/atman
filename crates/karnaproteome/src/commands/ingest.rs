use anyhow::{Context, Result};
use clap::Args as ClapArgs;
use proteome_core::{
    ingest::olink_explore::OlinkExploreLongCsv, DubeSampleIdParser, ProteomeIngest,
};
use std::path::PathBuf;

use crate::io::{write_measurements_long, write_proteins, write_samples};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Proteomics platform. v0.1 supports only `olink-explore-ngs`.
    #[arg(long, default_value = "olink-explore-ngs")]
    platform: String,

    /// Sample ID parser. v0.1 supports only `dube`.
    #[arg(long, default_value = "dube")]
    parser: String,

    /// Output directory for measurements.tsv, proteins.tsv, samples.tsv.
    #[arg(long)]
    output_dir: PathBuf,

    /// Input NPX files (one or more).
    inputs: Vec<PathBuf>,
}

pub fn run(args: Args) -> Result<()> {
    if args.platform != "olink-explore-ngs" {
        anyhow::bail!("platform {:?} not supported in v0.1", args.platform);
    }
    if args.parser != "dube" {
        anyhow::bail!("parser {:?} not supported in v0.1", args.parser);
    }
    if args.inputs.is_empty() {
        anyhow::bail!("no input files given");
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let adapter = OlinkExploreLongCsv;
    let parser = DubeSampleIdParser;
    let out = adapter.read(&args.inputs, &parser)?;

    write_measurements_long(&args.output_dir.join("measurements.tsv"), &out.measurements)?;
    write_proteins(&args.output_dir.join("proteins.tsv"), &out.proteins)?;
    write_samples(&args.output_dir.join("samples.tsv"), &out.samples)?;

    eprintln!(
        "ingest: {} rows, {} assays, {} samples, {} inputs",
        out.measurements.len(),
        out.proteins.len(),
        out.samples.len(),
        args.inputs.len()
    );
    Ok(())
}
