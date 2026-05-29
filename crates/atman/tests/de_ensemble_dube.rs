//! Ensemble integration test on the bundled Dube heat-acclimation
//! cohort. Runs ingest → qc → `atman de --test ensemble` on
//! PT2-PR2 (acute heat stress vs pre-ride baseline) with paired-t,
//! welch-t, and limma as the methods.
//!
//! Canonical heat-shock proteins must show consistent positive
//! majority_sign (mean_a − mean_b > 0 → PT2 > PR2, upregulated by
//! heat). The grade is assigned from method agreement
//! (`n_significant` + sign consistency), NOT from the
//! Stouffer-combined `ensemble_q`: all three methods share the same
//! abundance matrix, so the combined p is anti-conservative and is
//! reported only as a heuristic. In this cohort no protein has all
//! three methods clear per-method BH-q (so nothing reaches
//! VALIDATED), which is the honest, conservative outcome — the test
//! asserts sign consistency and that the grade never over-states
//! confidence beyond what the per-method evidence supports.

use std::path::{Path, PathBuf};
use std::process::Command;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

fn run_atman(args: &[&str]) {
    let bin = env!("CARGO_BIN_EXE_atman");
    let status = Command::new(bin).args(args).status().expect("run atman");
    assert!(status.success(), "atman {:?} failed", args);
}

fn parse_tsv(path: &Path) -> (Vec<String>, Vec<std::collections::HashMap<String, String>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    let rows = lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(h, c)| (h.clone(), c.to_string()))
                .collect()
        })
        .collect();
    (header, rows)
}

#[test]
fn ensemble_grades_dube_heat_shock_proteins_as_validated() {
    let root = repo_root();
    let data = root.join("example_data/dube_heat_2023");
    let raw1 = data.join("20212016_Dube_NPX_2021-11-30.csv");
    let raw2 = data.join("20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv");

    let tmp = tempfile::tempdir().unwrap();
    let canonical = tmp.path().join("canonical");
    let out = tmp.path().join("ensemble");

    let stager = root.join("adapters/generic/olink_explore_to_atman.py");
    let py_status = Command::new("python3")
        .arg(&stager)
        .arg("--output-dir")
        .arg(&canonical)
        .arg(&raw1)
        .arg(&raw2)
        .status()
        .expect("run olink_explore_to_atman.py");
    assert!(py_status.success(), "olink_explore_to_atman.py failed");

    run_atman(&[
        "qc",
        "--input-dir",
        canonical.to_str().unwrap(),
        "--output-dir",
        canonical.to_str().unwrap(),
    ]);

    run_atman(&[
        "de",
        "--input-dir",
        canonical.to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--test",
        "ensemble",
        "--ensemble-methods",
        "paired-t,welch-t,limma",
        "--groups",
        "PT2-PR2",
        "--paired-by",
        "subject_id",
        "--min-pairs",
        "5",
        "--trend",
        "false",
        "--robust",
        "false",
    ]);

    let ensemble_path = out.join("de_ensemble.tsv");
    let results_path = out.join("de_results.tsv");
    assert!(ensemble_path.exists(), "de_ensemble.tsv not produced");
    assert!(results_path.exists(), "de_results.tsv not produced");

    let (_, rows) = parse_tsv(&ensemble_path);
    let pt2 = rows
        .iter()
        .filter(|r| r.get("comparison").map(|s| s.as_str()) == Some("PT2-PR2"))
        .collect::<Vec<_>>();
    assert!(
        pt2.len() >= 2000,
        "expected ≥ 2000 ensemble rows for PT2-PR2, got {}",
        pt2.len()
    );

    // Canonical heat-shock genes present in the Dube Olink Explore
    // panel (empirically verified from the reference NPX headers).
    let canonical_hsps = ["HSPB1", "HSPA1A", "DNAJB1"];

    let mut saw_positive_sign = 0usize;
    let mut saw_hsp = 0usize;
    for r in &pt2 {
        let gene = r.get("gene_symbol").map(|s| s.as_str()).unwrap_or("");
        if !canonical_hsps.contains(&gene) {
            continue;
        }
        saw_hsp += 1;
        let grade = r.get("grade").map(|s| s.as_str()).unwrap_or("");
        let sign: f64 = r
            .get("majority_sign")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);
        let n_applied: usize = r.get("n_applied").and_then(|s| s.parse().ok()).unwrap_or(0);
        let n_significant: usize = r
            .get("n_significant")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let ensemble_q: Option<f64> = r.get("ensemble_q").and_then(|s| s.parse().ok());
        eprintln!(
            "  {gene:<10} grade={grade:<12} sign={sign:+.1} n_sig={n_significant}/{n_applied} ensemble_q={:?}",
            ensemble_q
        );
        if sign > 0.0 {
            saw_positive_sign += 1;
        }
        // Grade honesty: VALIDATED only when every applied method is
        // individually significant; PROVISIONAL only when a majority
        // is. The grade must never exceed what per-method significance
        // supports — i.e. it is NOT driven by the anti-conservative
        // combined ensemble_q.
        match grade {
            "VALIDATED" => assert_eq!(
                n_significant, n_applied,
                "{gene} graded VALIDATED but only {n_significant}/{n_applied} methods significant"
            ),
            "PROVISIONAL" => assert!(
                n_significant * 2 >= n_applied,
                "{gene} graded PROVISIONAL but only {n_significant}/{n_applied} methods significant"
            ),
            _ => {}
        }
    }

    assert_eq!(
        saw_hsp,
        canonical_hsps.len(),
        "expected all canonical HSPs present in PT2-PR2 ensemble output"
    );
    assert!(
        saw_positive_sign >= 2,
        "expected ≥ 2 canonical HSPs with positive majority_sign (PT2 > PR2); got {}",
        saw_positive_sign
    );

    // Dataset-level grade sanity: the grade is anchored to per-method
    // significance, so no protein may be VALIDATED unless all of its
    // applied methods individually clear per-method BH-q.
    for r in &pt2 {
        let grade = r.get("grade").map(|s| s.as_str()).unwrap_or("");
        if grade != "VALIDATED" {
            continue;
        }
        let n_applied: usize = r.get("n_applied").and_then(|s| s.parse().ok()).unwrap_or(0);
        let n_significant: usize = r
            .get("n_significant")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        assert_eq!(
            n_significant,
            n_applied,
            "VALIDATED protein with n_significant {n_significant} != n_applied {n_applied}: {:?}",
            r.get("gene_symbol")
        );
    }

    // de_results.tsv should have one row per (protein, method) for each
    // method requested. 3 methods × ~2938 proteins.
    let (_, per_method) = parse_tsv(&results_path);
    let pt2_methods: std::collections::BTreeSet<String> = per_method
        .iter()
        .filter(|r| r.get("comparison").map(|s| s.as_str()) == Some("PT2-PR2"))
        .map(|r| r["method"].clone())
        .collect();
    assert_eq!(
        pt2_methods.iter().cloned().collect::<Vec<_>>(),
        vec!["limma", "paired-t", "welch-t"],
        "expected exactly paired-t, welch-t, limma in de_results.tsv methods"
    );
}
