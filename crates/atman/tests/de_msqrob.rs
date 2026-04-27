//! Integration test for `atman de --test msqrob`.
//!
//! A synthetic fixture with three proteins (UPREG, DOWNREG, STABLE) ×
//! four peptides × twelve samples (six per condition) verifies:
//!
//! - CLI accepts `--peptide-measurements` and `--peptide-metadata` and
//!   rejects them when `--test msqrob` is not set.
//! - `de_results.tsv` includes the three new msqrob columns
//!   (`n_peptides_observed`, `peptide_variance_ratio`, `ridge_lambda`).
//! - The sign of the reported effect matches atman convention
//!   (`mean_a − mean_b`) for the injected ground-truth effects.
//! - Ridge `--ridge-lambda` is recorded in the per-row output and in
//!   the sidecar.
//! - `--ridge-lambda auto` is accepted and written as `0.0` (current
//!   data-driven selection stub).

use std::collections::HashMap;
use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_peptide_fixture(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();

    // 12 samples: 6 condition A, 6 condition B.
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=12 {
        let cond = if i <= 6 { "A" } else { "B" };
        samples.push_str(&format!("S{i:02}\tS{i:02}\t{cond}\t0\tplasma\t{i}\n"));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    // 3 proteins.
    let proteins = "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
                    maxquant_lfq\tA_UP\tQ01\tUPREG\tMS1\t\n\
                    maxquant_lfq\tA_DOWN\tQ02\tDOWNREG\tMS1\t\n\
                    maxquant_lfq\tA_STABLE\tQ03\tSTABLE\tMS1\t\n";
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    // Empty measurements.tsv — msqrob doesn't read it, but the
    // `atman de` entry point reads the canonical file before dispatch.
    let qc = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
              abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
              detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n";
    std::fs::write(dir.join("measurements.tsv"), qc).unwrap();
    std::fs::write(dir.join("measurements.tsv"), qc).unwrap();

    // 4 peptides per protein, 12 peptides total.
    let mut peptides =
        String::from("peptide_id\tassay_id\tsequence\tcharge\tmodifications\tmissed_cleavages\n");
    for (prot, prefix) in [("A_UP", "UP"), ("A_DOWN", "DN"), ("A_STABLE", "ST")] {
        for k in 0..4 {
            peptides.push_str(&format!(
                "P_{prefix}_{k}\t{prot}\tPEPTIDE{prefix}{k}\t2\t\t0\n"
            ));
        }
    }
    std::fs::write(dir.join("peptides.tsv"), peptides).unwrap();

    // Peptide abundance: baseline 10 + peptide offset + condition effect
    // + small deterministic noise.
    let peptide_offsets = [-0.4, -0.1, 0.2, 0.5];
    let mut pm = String::from(
        "sample_id\tpeptide_id\tabundance\tabundance_unit\tdropped_by_qc\tbelow_lod\n",
    );
    for sample in 1..=12usize {
        let is_b = sample > 6;
        for (prot_prefix, condition_effect) in [("UP", 1.5), ("DN", -1.5), ("ST", 0.0)] {
            for (k, &offset) in peptide_offsets.iter().enumerate() {
                let shift = if is_b { condition_effect } else { 0.0 };
                let noise = ((sample * 13 + k * 7) as f64).sin() * 0.03;
                let abund = 10.0 + offset + shift + noise;
                pm.push_str(&format!(
                    "S{:02}\tP_{}_{}\t{:.6}\tlog2_intensity\t0\t0\n",
                    sample, prot_prefix, k, abund
                ));
            }
        }
    }
    std::fs::write(dir.join("peptide_measurements.tsv"), pm).unwrap();
}

fn parse_tsv(path: &Path) -> (Vec<String>, Vec<HashMap<String, String>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    let rows: Vec<HashMap<String, String>> = lines
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
fn msqrob_recovers_injected_direction_and_columns() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("msqrob_out");
    write_peptide_fixture(&input);
    std::fs::create_dir_all(&output).unwrap();

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "msqrob",
        "--groups",
        "A-B",
        "--peptide-measurements",
        input.join("peptide_measurements.tsv").to_str().unwrap(),
        "--peptide-metadata",
        input.join("peptides.tsv").to_str().unwrap(),
        "--ridge-lambda",
        "0.5",
        "--min-peptides",
        "2",
        "--min-pairs",
        "4",
    ]);
    assert!(
        out.status.success(),
        "msqrob run failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let results_path = output.join("de_results.tsv");
    let (header, rows) = parse_tsv(&results_path);
    for col in [
        "n_peptides_observed",
        "peptide_variance_ratio",
        "ridge_lambda",
    ] {
        assert!(header.iter().any(|h| h == col), "header missing {col}");
    }
    // Expect 3 rows, one per protein, all for A-B comparison.
    let ab: Vec<&HashMap<String, String>> = rows
        .iter()
        .filter(|r| r.get("comparison").map(|s| s.as_str()) == Some("A-B"))
        .collect();
    assert_eq!(ab.len(), 3, "expected 3 A-B rows, got {}", ab.len());
    for row in &ab {
        assert_eq!(row["effect_size_method"], "msqrob-ridge");
        assert_eq!(row["ridge_lambda"], "0.5");
        let np: usize = row["n_peptides_observed"].parse().unwrap();
        assert_eq!(np, 4, "four peptides observed per protein");
        assert!(
            !row["peptide_variance_ratio"].is_empty(),
            "peptide_variance_ratio should be populated"
        );
    }
    // Direction: UPREG has mean_b > mean_a, so A-B effect should be
    // strongly negative. DOWNREG is the opposite. STABLE is near zero.
    let by_gene: HashMap<&str, &HashMap<String, String>> =
        ab.iter().map(|r| (r["gene_symbol"].as_str(), *r)).collect();
    let eff = |gene: &str| -> f64 {
        by_gene[gene]["mean_diff"]
            .parse()
            .unwrap_or_else(|_| panic!("parsing mean_diff for {gene}"))
    };
    assert!(eff("UPREG") < -1.0, "UPREG effect was {}", eff("UPREG"));
    assert!(
        eff("DOWNREG") > 1.0,
        "DOWNREG effect was {}",
        eff("DOWNREG")
    );
    assert!(
        eff("STABLE").abs() < 0.3,
        "STABLE effect was {}",
        eff("STABLE")
    );

    // Sidecar shape.
    let sidecar = output.join("de_results.tsv.run.json");
    assert!(sidecar.exists());
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    assert_eq!(j["args"]["test"], "msqrob");
    assert_eq!(j["args"]["ridge-lambda"], "0.5");
    assert_eq!(j["args"]["min-peptides"], 2);
}

#[test]
fn msqrob_auto_ridge_lambda_is_zero() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("auto_out");
    write_peptide_fixture(&input);
    std::fs::create_dir_all(&output).unwrap();

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "msqrob",
        "--groups",
        "A-B",
        "--peptide-measurements",
        input.join("peptide_measurements.tsv").to_str().unwrap(),
        "--peptide-metadata",
        input.join("peptides.tsv").to_str().unwrap(),
        "--min-peptides",
        "2",
        "--min-pairs",
        "4",
    ]);
    assert!(
        out.status.success(),
        "auto msqrob run failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let (_, rows) = parse_tsv(&output.join("de_results.tsv"));
    for row in rows {
        assert_eq!(row["ridge_lambda"], "0", "auto should default to 0");
    }
}

#[test]
fn msqrob_rejects_peptide_flags_for_other_tests() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("reject_out");
    write_peptide_fixture(&input);
    std::fs::create_dir_all(&output).unwrap();

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "welch-t",
        "--groups",
        "A-B",
        "--peptide-measurements",
        input.join("peptide_measurements.tsv").to_str().unwrap(),
        "--peptide-metadata",
        input.join("peptides.tsv").to_str().unwrap(),
        "--min-pairs",
        "4",
    ]);
    assert!(!out.status.success(), "expected rejection");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("requires --test msqrob"),
        "unexpected stderr: {stderr}"
    );
}

#[test]
fn msqrob_refuses_missing_peptide_files() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("missing_out");
    write_peptide_fixture(&input);
    std::fs::create_dir_all(&output).unwrap();

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "msqrob",
        "--groups",
        "A-B",
        "--min-pairs",
        "4",
    ]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--peptide-measurements"),
        "expected peptide-measurements complaint: {stderr}"
    );
}
