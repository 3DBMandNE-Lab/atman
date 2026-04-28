use anyhow::{Context, Result};
use atman_core::stats::ranks;
use atman_core::{Abundance, AssayId, MeasurementRecord, Platform};
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;
use std::time::SystemTime;

use crate::io::{
    hash_canonical_inputs, read_measurements_long, read_proteins, read_samples, sidecar_path_for,
    write_measurements_long, write_proteins, write_run_sidecar, write_samples,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Canonical Atman input directory.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for rank-transformed canonical TSVs.
    #[arg(long)]
    output_dir: PathBuf,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating {:?}", args.output_dir))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let proteins = read_proteins(&args.input_dir.join("proteins.tsv"))?;
    let samples_out = args.output_dir.join("samples.tsv");
    let proteins_out = args.output_dir.join("proteins.tsv");
    let measurements_out = args.output_dir.join("measurements.tsv");
    write_samples(&samples_out, &samples)?;
    write_proteins(&proteins_out, &proteins)?;

    let measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;
    let ranked_measurements = rank_records(measurements);
    write_measurements_long(&measurements_out, &ranked_measurements)?;

    eprintln!(
        "within-cohort-rank: samples={} proteins={} measurements={}",
        samples.len(),
        proteins.len(),
        ranked_measurements.len()
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    let outputs = [
        samples_out.clone(),
        proteins_out.clone(),
        measurements_out.clone(),
    ];
    let sidecar = sidecar_path_for(&measurements_out);
    write_run_sidecar(
        &sidecar,
        "within-cohort-rank",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("within-cohort-rank: sidecar={}", sidecar.display());
    Ok(())
}

fn rank_records(mut records: Vec<MeasurementRecord>) -> Vec<MeasurementRecord> {
    let mut by_assay: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for (idx, record) in records.iter().enumerate() {
        if record.effective_abundance().is_some_and(f64::is_finite) {
            by_assay
                .entry(assay_key(record.platform, &record.assay_id))
                .or_default()
                .push(idx);
        }
    }

    for indexes in by_assay.values() {
        let values: Vec<f64> = indexes
            .iter()
            .map(|idx| records[*idx].abundance.as_f64())
            .collect();
        let ranks = ranks(&values);
        for (idx, rank) in indexes.iter().zip(ranks) {
            records[*idx].abundance = Abundance::Raw(rank);
            records[*idx].abundance_raw = Abundance::Raw(rank);
            records[*idx].npx_source_str = rank.to_string();
            records[*idx].detection_limit = atman_core::DetectionLimit(None);
            records[*idx].below_lod = false;
        }
    }
    records
}

fn assay_key(platform: Platform, assay_id: &AssayId) -> (String, String) {
    (platform.as_str().to_string(), assay_id.0.clone())
}
