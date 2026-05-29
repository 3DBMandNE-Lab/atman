use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use csv::ReaderBuilder;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use std::time::SystemTime;

use crate::io::{atomic_write, format_float, sidecar_path_for, write_run_sidecar};

#[derive(ClapArgs, Debug)]
pub struct CounterfactualArgs {
    /// Loadings TSV from `atman decompose ica` (columns: program, assay_id,
    /// gene_symbol, loading).
    #[arg(long)]
    loadings: PathBuf,

    /// Activations TSV from `atman decompose ica` (columns: cohort,
    /// subject_id, sample_id, program, activation).
    #[arg(long)]
    activations: PathBuf,

    /// Program to manipulate (e.g. `program_04`).
    #[arg(long)]
    set_program: String,

    /// Value to set the target program's activation to for all samples.
    #[arg(long, default_value_t = 0.0)]
    to: f64,

    /// Output TSV for the delta matrix (protein × sample contribution of
    /// the manipulated program). Wide format: assay_id, gene_symbol,
    /// then one column per sample_id.
    #[arg(long)]
    output: PathBuf,
}

// ── counterfactual ──────────────────────────────────────────────────

pub(super) fn run_counterfactual(args: CounterfactualArgs) -> Result<()> {
    let started_at = SystemTime::now();

    // Read loadings: program, assay_id, gene_symbol, loading
    let mut loadings_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&args.loadings)
        .with_context(|| format!("opening {:?}", args.loadings))?;
    let mut loadings_map: BTreeMap<(String, String), BTreeMap<String, f64>> = BTreeMap::new();
    let mut all_programs: BTreeSet<String> = BTreeSet::new();
    for record in loadings_reader.records() {
        let r = record?;
        let program = r.get(0).unwrap_or("").to_string();
        let assay_id = r.get(1).unwrap_or("").to_string();
        let gene_symbol = r.get(2).unwrap_or("").to_string();
        let loading: f64 = r
            .get(3)
            .unwrap_or("")
            .parse()
            .with_context(|| format!("parsing loading for {}/{}", program, assay_id))?;
        all_programs.insert(program.clone());
        loadings_map
            .entry((assay_id, gene_symbol))
            .or_default()
            .insert(program, loading);
    }

    if !all_programs.contains(&args.set_program) {
        bail!(
            "--set-program {:?} not found in loadings; available: {}",
            args.set_program,
            all_programs.iter().cloned().collect::<Vec<_>>().join(", ")
        );
    }

    // Read activations: cohort, subject_id, sample_id, program, activation
    let mut act_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&args.activations)
        .with_context(|| format!("opening {:?}", args.activations))?;
    let mut activations_map: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut sample_ids: BTreeSet<String> = BTreeSet::new();
    for record in act_reader.records() {
        let r = record?;
        let sample_id = r.get(2).unwrap_or("").to_string();
        let program = r.get(3).unwrap_or("").to_string();
        let activation: f64 = r
            .get(4)
            .unwrap_or("")
            .parse()
            .with_context(|| format!("parsing activation for {}/{}", sample_id, program))?;
        sample_ids.insert(sample_id.clone());
        activations_map.insert((sample_id, program), activation);
    }
    let sample_ids: Vec<String> = sample_ids.into_iter().collect();

    // Compute delta: for each protein × sample, the contribution of the
    // target program is loading[protein, program] × (activation_original - to).
    let mut buf = String::from("assay_id\tgene_symbol");
    for sid in &sample_ids {
        buf.push('\t');
        buf.push_str(sid);
    }
    buf.push('\n');

    let mut n_proteins = 0u64;
    for ((assay_id, gene_symbol), per_program) in &loadings_map {
        let loading = per_program.get(&args.set_program).copied().unwrap_or(0.0);
        buf.push_str(assay_id);
        buf.push('\t');
        buf.push_str(gene_symbol);
        for sid in &sample_ids {
            let original_activation = activations_map
                .get(&(sid.clone(), args.set_program.clone()))
                .copied()
                .unwrap_or(0.0);
            let delta = loading * (original_activation - args.to);
            buf.push('\t');
            buf.push_str(&format_float(delta));
        }
        buf.push('\n');
        n_proteins += 1;
    }

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)?;
    }
    atomic_write(&args.output, buf.as_bytes())?;

    eprintln!(
        "decompose counterfactual: {} proteins × {} samples, program={} set to {}, output={}",
        n_proteins,
        sample_ids.len(),
        args.set_program,
        args.to,
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = crate::io::hash_labeled_inputs(&[
        ("loadings", args.loadings.as_path()),
        ("activations", args.activations.as_path()),
    ])?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "decompose counterfactual",
        json!({
            "loadings": args.loadings.display().to_string(),
            "activations": args.activations.display().to_string(),
            "set-program": args.set_program,
            "to": args.to,
            "output": args.output.display().to_string(),
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("decompose counterfactual: sidecar={}", sidecar.display());
    Ok(())
}
