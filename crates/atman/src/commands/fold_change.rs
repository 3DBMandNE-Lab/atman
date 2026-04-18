use anyhow::{Context, Result};
use atman_core::fold_change::{compute_log2_fc, Comparison, FoldChangeInput};
use clap::Args as ClapArgs;
use std::path::PathBuf;

use super::parse_comparisons;
use crate::io::{read_measurements_long, read_samples, write_fold_change_panel};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[arg(long)]
    input_dir: PathBuf,

    #[arg(long)]
    output_dir: PathBuf,

    /// Comma-separated comparisons in `A-B` form. Example: "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1".
    #[arg(long)]
    groups: String,

    /// Split dimension. Supports `panel`.
    #[arg(long, default_value = "panel")]
    split_by: String,
}

pub fn run(args: Args) -> Result<()> {
    if args.split_by != "panel" {
        anyhow::bail!("split-by {:?} not supported", args.split_by);
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let comparisons: Vec<Comparison> = parse_comparisons(&args.groups)?
        .into_iter()
        .map(|(a, b)| Comparison { a, b })
        .collect();

    let measurements = read_measurements_long(&args.input_dir.join("qc_measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;

    use std::collections::HashMap;
    let sample_by_id: HashMap<&str, &atman_core::Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    // First: collect the universe of (panel, gene) pairs from ALL measurements
    // regardless of QC state. Dube emits a row per assay in each panel even if
    // some have zero non-missing values after QC.
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
    for panel in &output.panels {
        let p = write_fold_change_panel(&args.output_dir, panel, &comparisons)?;
        eprintln!("fold-change: wrote {:?} ({} assays)", p, panel.assays.len());
    }
    Ok(())
}
