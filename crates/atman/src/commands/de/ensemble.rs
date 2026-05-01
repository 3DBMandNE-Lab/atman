//! Cross-method consensus dispatcher (`atman de --test ensemble`).
//!
//! Fans out to every requested DE method in a tempdir, aggregates
//! per (comparison, protein) into a Stouffer-combined ensemble
//! p / q, and assigns a VALIDATED / PROVISIONAL / INSUFFICIENT grade.
//! Methods that fail (e.g. `msqrob` without `--peptide-measurements`)
//! are collected into `methods_skipped` rather than aborting the
//! whole run.

use anyhow::{Context, Result};
use atman_core::ensemble::{aggregate_per_protein, assign_grade, EnsembleInput, GradeThresholds};
use serde_json::json;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::SystemTime;
use tempfile::TempDir;

use super::{run, Args};
use crate::io::{
    hash_canonical_inputs, read_de_results, sidecar_path_for, write_de_ensemble, write_de_results,
    write_run_sidecar, DeResultRow, EnsembleRow,
};

pub(super) fn run_ensemble(args: Args, started_at: SystemTime) -> Result<()> {
    let requested: Vec<String> = args
        .ensemble_methods
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if requested.is_empty() {
        anyhow::bail!("--ensemble-methods resolved to zero methods");
    }
    for m in &requested {
        if !matches!(
            m.as_str(),
            "paired-t" | "moderated" | "welch-t" | "ols" | "mixed" | "limma" | "msqrob"
        ) {
            anyhow::bail!("unsupported ensemble method {:?}", m);
        }
    }

    let thresholds = GradeThresholds {
        q_threshold: args.ensemble_q_threshold,
        validated_sign_fraction: args.ensemble_sign_fraction,
        provisional_sign_fraction: args.ensemble_provisional_fraction,
    };

    let tmp = TempDir::new().with_context(|| "creating ensemble tempdir")?;
    let mut merged_rows: Vec<DeResultRow> = Vec::new();
    let mut skipped: Vec<String> = Vec::new();
    let mut applied: Vec<String> = Vec::new();

    for method in &requested {
        let mut sub = args.clone();
        sub.test = method.clone();
        sub.output_dir = tmp.path().join(method);
        // Strip peptide flags for methods that don't consume them.
        match method.as_str() {
            "msqrob" => { /* keeps both */ }
            "limma" => {
                sub.peptide_measurements = None;
            }
            _ => {
                sub.peptide_measurements = None;
                sub.peptide_metadata = None;
            }
        }
        // Strip formula / mixed flags for methods that can't use them.
        match method.as_str() {
            "ols" => {
                sub.fixed = None;
                sub.random = None;
                sub.adjust_for = vec![];
            }
            "mixed" => {
                sub.covariates = None;
                sub.design = None;
                sub.contrast = None;
                sub.per_subject_proxy = None;
                sub.adjust_for = vec![];
            }
            "limma" => {
                // --adjust-for is limma-specific; keep it on the limma sub-call.
                sub.fixed = None;
                sub.random = None;
                sub.covariates = None;
                sub.design = None;
                sub.contrast = None;
                sub.per_subject_proxy = None;
            }
            _ => {
                sub.covariates = None;
                sub.design = None;
                sub.contrast = None;
                sub.per_subject_proxy = None;
                sub.fixed = None;
                sub.random = None;
                sub.adjust_for = vec![];
            }
        }
        match run(sub.clone()) {
            Ok(()) => {
                let results_path = sub.output_dir.join("de_results.tsv");
                let rows = read_de_results(&results_path)
                    .with_context(|| format!("reading sub-method de_results for {method:?}"))?;
                merged_rows.extend(rows);
                applied.push(method.clone());
            }
            Err(e) => {
                let short = e.to_string().lines().next().unwrap_or("").to_string();
                eprintln!("ensemble: method {method:?} skipped — {short}");
                skipped.push(format!("{method}:{short}"));
            }
        }
    }

    // Aggregate per (comparison, assay_id).
    type Key = (String, String);
    let mut grouped: BTreeMap<
        Key,
        (Vec<EnsembleInput>, String, String, String), // (inputs, panel, gene, uniprot)
    > = BTreeMap::new();
    for r in &merged_rows {
        let key = (r.comparison.clone(), r.assay_id.clone());
        let entry = grouped.entry(key).or_insert_with(|| {
            (
                Vec::new(),
                r.panel.clone(),
                r.gene_symbol.clone(),
                r.uniprot.clone(),
            )
        });
        entry.0.push(EnsembleInput {
            method: r.method.clone(),
            mean_diff: r.mean_diff,
            p_value: r.p_value,
            bh_q: r.bh_q,
        });
    }

    let mut ensemble_rows: Vec<EnsembleRow> = Vec::new();
    let mut by_comparison: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    let methods_skipped_str = if skipped.is_empty() {
        String::new()
    } else {
        skipped.join("|")
    };
    for ((comparison, assay_id), (inputs, panel, gene, uniprot)) in &grouped {
        let agg = aggregate_per_protein(inputs, thresholds);
        let idx = ensemble_rows.len();
        by_comparison
            .entry(comparison.clone())
            .or_default()
            .push(idx);
        ensemble_rows.push(EnsembleRow {
            comparison: comparison.clone(),
            panel: panel.clone(),
            assay_id: assay_id.clone(),
            gene_symbol: gene.clone(),
            uniprot: uniprot.clone(),
            n_applied: agg.n_applied,
            n_significant: agg.n_significant,
            n_sign_consistent: agg.n_sign_consistent,
            majority_sign: agg.majority_sign,
            ensemble_p: agg.ensemble_p,
            ensemble_q: None,
            grade: String::new(), // filled after BH
            methods_applied: agg.methods_applied.join(","),
            methods_skipped: methods_skipped_str.clone(),
        });
    }
    // BH on ensemble_p within each comparison, then assign grade per
    // protein from (ensemble_q, sign consistency).
    for indices in by_comparison.values() {
        let ps: Vec<Option<f64>> = indices
            .iter()
            .map(|&i| ensemble_rows[i].ensemble_p)
            .collect();
        let qs = atman_core::bh_fdr(&ps);
        for (j, &i) in indices.iter().enumerate() {
            ensemble_rows[i].ensemble_q = qs[j];
            let row = &ensemble_rows[i];
            let grade = assign_grade(
                row.ensemble_q,
                row.n_applied,
                row.n_sign_consistent,
                thresholds,
            );
            ensemble_rows[i].grade = grade.as_str().to_string();
        }
    }

    // Stable output ordering: (comparison, grade-priority, ensemble_q asc).
    ensemble_rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| grade_rank(&a.grade).cmp(&grade_rank(&b.grade)))
            .then_with(|| match (a.ensemble_q, b.ensemble_q) {
                (Some(qa), Some(qb)) => qa.partial_cmp(&qb).unwrap_or(std::cmp::Ordering::Equal),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            })
            .then_with(|| a.assay_id.cmp(&b.assay_id))
    });
    merged_rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| a.assay_id.cmp(&b.assay_id))
            .then_with(|| a.method.cmp(&b.method))
    });

    let results_path = args.output_dir.join("de_results.tsv");
    let ensemble_path = args.output_dir.join("de_ensemble.tsv");
    write_de_results(&results_path, &merged_rows)?;
    write_de_ensemble(&ensemble_path, &ensemble_rows)?;

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &[
            "measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
    )?;
    let outputs: Vec<PathBuf> = vec![results_path.clone(), ensemble_path.clone()];
    let sidecar = sidecar_path_for(&results_path);
    write_run_sidecar(
        &sidecar,
        "de",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
            "test": "ensemble",
            "groups": args.groups,
            "min-pairs": args.min_pairs,
            "peptide-measurements": args.peptide_measurements.as_ref().map(|p| p.display().to_string()),
            "peptide-metadata": args.peptide_metadata.as_ref().map(|p| p.display().to_string()),
            "ensemble-methods": args.ensemble_methods,
            "ensemble-q-threshold": args.ensemble_q_threshold,
            "ensemble-provisional-fraction": args.ensemble_provisional_fraction,
            "ensemble-sign-fraction": args.ensemble_sign_fraction,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        {
            // Resolved-value fields: which methods ran vs were
            // auto-skipped (e.g. msqrob without --peptide-measurements).
            // Kept outside `args` so user-intent vs runtime-outcome
            // don't conflate in the sidecar.
            let mut extras = serde_json::Map::new();
            extras.insert(
                "ensemble_methods_applied".into(),
                serde_json::json!(applied),
            );
            extras.insert(
                "ensemble_methods_skipped".into(),
                serde_json::json!(skipped),
            );
            Some(extras)
        },
    )?;
    eprintln!(
        "de ensemble: {} rows, {} applied: [{}], {} skipped: [{}]",
        ensemble_rows.len(),
        applied.len(),
        applied.join(", "),
        skipped.len(),
        skipped.join(", ")
    );
    eprintln!("de: sidecar={}", sidecar.display());
    Ok(())
}

fn grade_rank(grade: &str) -> u8 {
    match grade {
        "VALIDATED" => 0,
        "PROVISIONAL" => 1,
        _ => 2,
    }
}
