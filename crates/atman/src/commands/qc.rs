use anyhow::{Context, Result};
use atman_core::qc::apply_mask_warn_fail_all;
use atman_core::MeasurementRecord;
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    hash_canonical_inputs, read_measurements_long, sidecar_path_for, write_measurements_long,
    write_run_sidecar,
};

const RULE_NAME: &str = "mask-warn-fail";

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing `measurements.tsv` from upstream ingest.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for the rewritten `measurements.tsv` and the
    /// `qc_report.tsv` summary side-channel. The common case is
    /// `--output-dir` equal to `--input-dir` (in-place update).
    #[arg(long)]
    output_dir: PathBuf,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let in_path = args.input_dir.join("measurements.tsv");
    let mut records = read_measurements_long(&in_path)?;
    let before = records.iter().filter(|r| !r.dropped_by_qc).count();
    apply_mask_warn_fail_all(&mut records);
    let after = records.iter().filter(|r| !r.dropped_by_qc).count();
    let masked = before - after;

    let measurements_out = args.output_dir.join("measurements.tsv");
    let report_out = args.output_dir.join("qc_report.tsv");
    write_measurements_long(&measurements_out, &records)?;
    write_qc_report(&report_out, &records)?;
    eprintln!(
        "qc: rule={} total={} masked={} passed={}",
        RULE_NAME,
        records.len(),
        masked,
        after,
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(&args.input_dir, &["measurements.tsv"])?;
    let outputs = [measurements_out.clone(), report_out.clone()];
    let sidecar = sidecar_path_for(&measurements_out);
    write_run_sidecar(
        &sidecar,
        "qc",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
            "rule": RULE_NAME,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("qc: sidecar={}", sidecar.display());
    Ok(())
}

fn write_qc_report(path: &Path, records: &[MeasurementRecord]) -> Result<()> {
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
        writeln!(
            w,
            "{}\t{}\t{}\t{:.6}\t{}",
            sample_id, n, masked, rate, RULE_NAME
        )?;
    }
    w.flush()?;
    Ok(())
}
