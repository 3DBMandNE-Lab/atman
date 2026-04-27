use anyhow::{Context, Result};
use atman_core::qc::apply_mask_warn_fail_all;
use atman_core::MeasurementRecord;
use clap::Args as ClapArgs;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::io::{read_measurements_long, write_measurements_long};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing `measurements.tsv` from `ingest`.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for the rewritten `measurements.tsv` and the
    /// `qc_report.tsv` summary side-channel. The common case is
    /// `--output-dir` equal to `--input-dir` (in-place update).
    #[arg(long)]
    output_dir: PathBuf,

    /// QC rule name. Supports `mask-warn-fail` (default): mark
    /// `dropped_by_qc=true` when either `qc_sample` or `qc_assay` is not
    /// `Pass`. This is the Olink-recommended default.
    #[arg(long, default_value = "mask-warn-fail")]
    rule: String,
}

pub fn run(args: Args) -> Result<()> {
    if args.rule != "mask-warn-fail" {
        anyhow::bail!("rule {:?} not supported", args.rule);
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let in_path = args.input_dir.join("measurements.tsv");
    let mut records = read_measurements_long(&in_path)?;
    let before = records.iter().filter(|r| !r.dropped_by_qc).count();
    apply_mask_warn_fail_all(&mut records);
    let after = records.iter().filter(|r| !r.dropped_by_qc).count();
    let masked = before - after;

    write_measurements_long(&args.output_dir.join("measurements.tsv"), &records)?;
    write_qc_report(&args.output_dir.join("qc_report.tsv"), &records, &args.rule)?;
    eprintln!(
        "qc: rule={} total={} masked={} passed={}",
        args.rule,
        records.len(),
        masked,
        after,
    );
    Ok(())
}

fn write_qc_report(path: &Path, records: &[MeasurementRecord], rule: &str) -> Result<()> {
    let mut per_sample: BTreeMap<&str, (usize, usize)> = BTreeMap::new();
    for r in records {
        let entry = per_sample.entry(r.sample_id.as_str()).or_insert((0, 0));
        entry.0 += 1;
        if r.dropped_by_qc {
            entry.1 += 1;
        }
    }

    let f = File::create(path).with_context(|| format!("creating {:?}", path))?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "sample_id\tn_measurements\tn_masked\tmask_rate\trule_applied"
    )?;
    for (sample_id, (n, masked)) in &per_sample {
        let rate = if *n > 0 {
            *masked as f64 / *n as f64
        } else {
            0.0
        };
        writeln!(w, "{}\t{}\t{}\t{:.6}\t{}", sample_id, n, masked, rate, rule)?;
    }
    w.flush()?;
    Ok(())
}
