use anyhow::{Context, Result};
use clap::Args as ClapArgs;
use proteome_core::qc::apply_dube_rule_all;
use std::path::PathBuf;

use crate::io::{read_measurements_long, write_measurements_long};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing `measurements.tsv` from `ingest`.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for `qc_measurements.tsv`.
    #[arg(long)]
    output_dir: PathBuf,

    /// QC rule name. v0.1 supports only `dube`.
    #[arg(long, default_value = "dube")]
    rule: String,
}

pub fn run(args: Args) -> Result<()> {
    if args.rule != "dube" {
        anyhow::bail!("rule {:?} not supported in v0.1", args.rule);
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let in_path = args.input_dir.join("measurements.tsv");
    let mut records = read_measurements_long(&in_path)?;
    let before = records.iter().filter(|r| !r.dropped_by_qc).count();
    apply_dube_rule_all(&mut records);
    let after = records.iter().filter(|r| !r.dropped_by_qc).count();
    let masked = before - after;

    write_measurements_long(
        &args.output_dir.join("qc_measurements.tsv"),
        &records,
    )?;
    eprintln!(
        "qc: rule={} total={} masked={} passed={}",
        args.rule,
        records.len(),
        masked,
        after,
    );
    Ok(())
}
