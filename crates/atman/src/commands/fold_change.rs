use anyhow::{Context, Result};
use atman_core::fold_change::{compute_log2_fc, Comparison, FoldChangeInput};
use clap::Args as ClapArgs;
use serde_json::json;
use std::path::PathBuf;
use std::time::SystemTime;

use super::parse_comparisons;
use crate::io::{
    hash_canonical_inputs, read_measurements_long, read_samples, sidecar_path_for,
    write_fold_change_panel, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for per-panel fold-change CSVs and the
    /// `fold_change.run.json` sidecar.
    #[arg(long)]
    output_dir: PathBuf,

    /// Comma-separated comparisons in `A-B` form. Example: "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1".
    #[arg(long)]
    groups: String,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let comparisons: Vec<Comparison> = parse_comparisons(&args.groups)?
        .into_iter()
        .map(|(a, b)| Comparison { a, b })
        .collect();

    let measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;

    use std::collections::HashMap;
    let sample_by_id: HashMap<&str, &atman_core::Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    // First: collect the universe of (panel, gene) pairs from ALL measurements
    // regardless of QC state. The published reference emits a row per assay in
    // each panel even if some have zero non-missing values after QC.
    let mut universe: std::collections::BTreeSet<(String, String)> =
        std::collections::BTreeSet::new();
    for m in &measurements {
        if let (Some(panel), Some(gene)) = (m.panel.clone(), m.gene_symbol.clone()) {
            universe.insert((panel, gene));
        }
    }

    // Second: collect the actual data cells, filtered for QC + biological.
    let cells: Vec<(String, String, String, String, Option<f64>)> = measurements
        .iter()
        .filter_map(|m| {
            let abundance = m.effective_abundance()?;
            let s = sample_by_id.get(m.sample_id.as_str())?;
            if s.is_control {
                return None;
            }
            let participant = s.subject_id.clone()?;
            let exposure = s.condition.clone()?;
            let panel = m.panel.clone()?;
            let gene = m.gene_symbol.clone()?;
            Some((panel, gene, participant, exposure, Some(abundance)))
        })
        .collect();
    let mut input = FoldChangeInput::from_cells(cells);
    for (panel, gene) in &universe {
        input.ensure_assay(panel, gene);
    }

    let output = compute_log2_fc(&input, &comparisons);
    let mut output_paths: Vec<PathBuf> = Vec::with_capacity(output.panels.len());
    for panel in &output.panels {
        let p = write_fold_change_panel(&args.output_dir, panel, &comparisons)?;
        eprintln!("fold-change: wrote {:?} ({} assays)", p, panel.assays.len());
        output_paths.push(p);
    }

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    let primary_output = output_paths
        .first()
        .cloned()
        .unwrap_or_else(|| args.output_dir.join("fold_change.run"));
    let sidecar = sidecar_path_for(&primary_output);
    write_run_sidecar(
        &sidecar,
        "fold-change",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
            "groups": args.groups,
        }),
        &inputs_sha256,
        &output_paths,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("fold-change: sidecar={}", sidecar.display());
    Ok(())
}
