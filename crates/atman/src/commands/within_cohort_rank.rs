use anyhow::{Context, Result};
use atman_core::stats::ranks;
use atman_core::{Abundance, AssayId, MeasurementRecord, Platform};
use clap::Args as ClapArgs;
use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

use crate::io::{
    read_measurements_long, read_proteins, read_samples, write_measurements_long, write_proteins,
    write_samples,
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
    fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating {:?}", args.output_dir))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let proteins = read_proteins(&args.input_dir.join("proteins.tsv"))?;
    write_samples(&args.output_dir.join("samples.tsv"), &samples)?;
    write_proteins(&args.output_dir.join("proteins.tsv"), &proteins)?;

    let measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;
    let ranked_measurements = rank_records(measurements);
    write_measurements_long(
        &args.output_dir.join("measurements.tsv"),
        &ranked_measurements,
    )?;

    let qc_measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;
    let ranked_qc = rank_records(qc_measurements);
    write_measurements_long(&args.output_dir.join("measurements.tsv"), &ranked_qc)?;

    eprintln!(
        "within-cohort-rank: samples={} proteins={} measurements={}",
        samples.len(),
        proteins.len(),
        ranked_measurements.len()
    );
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
