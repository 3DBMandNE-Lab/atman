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
//! - The `--ridge-lambda` default applies no shrinkage, and the
//!   unimplemented `auto` is refused under msqrob and reported as
//!   ignored elsewhere, rather than silently meaning `0.0`.

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

/// The default penalty is no shrinkage. This used to be spelled `auto`,
/// which read as data-driven selection while doing nothing; the default
/// is now literally `0.0` and the numeric behaviour is unchanged.
#[test]
fn msqrob_default_ridge_lambda_applies_no_shrinkage() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("default_out");
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
        "default msqrob run failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let (_, rows) = parse_tsv(&output.join("de_results.tsv"));
    for row in rows {
        assert_eq!(
            row["ridge_lambda"], "0",
            "default should apply no shrinkage"
        );
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

// ---- --ridge-lambda auto is not implemented ---------------------------
//
// `auto` used to be the flag's DEFAULT and silently resolved to 0.0, so
// every `de` run recorded `ridge-lambda: auto` in its sidecar whether or
// not msqrob ran, and a user who chose `auto` got no shrinkage and no
// notice. The default is now the honest `0.0`; `auto` is only ever
// present because someone asked for it, and asking gets an answer.

/// Minimal two-group canonical dir with real abundances, for exercising
/// non-msqrob paths (the peptide fixture writes an empty
/// `measurements.tsv` on purpose).
fn write_two_group_fixture(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=12 {
        let cond = if i <= 6 { "A" } else { "B" };
        samples.push_str(&format!("S{i:02}\tS{i:02}\t{cond}\t0\tplasma\t{i}\n"));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    std::fs::write(
        dir.join("proteins.tsv"),
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
         olink_explore_ngs\tA001\tQ00001\tG001\tP1\t\n\
         olink_explore_ngs\tA002\tQ00002\tG002\tP1\t\n",
    )
    .unwrap();
    let mut m = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0;
    for i in 1..=12 {
        for (j, assay) in ["A001", "A002"].iter().enumerate() {
            order += 1;
            // Group B shifted up on A001; A002 is flat.
            let base = if j == 0 && i > 6 { 3.0 } else { 1.0 };
            let v = base + (i as f64) * 0.01 + (j as f64) * 0.5;
            m.push_str(&format!(
                "olink_explore_ngs\tS{i:02}\t{assay}\tG{n:03}\tP1\t{v:.6}\t{v:.6}\t{v:.6}\t\
                 log2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n",
                n = j + 1
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), m).unwrap();
}

fn read_sidecar(output_dir: &Path) -> serde_json::Value {
    let p = output_dir.join("de_results.tsv.run.json");
    let text = std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()));
    serde_json::from_str(&text).expect("parse sidecar")
}

#[test]
fn msqrob_refuses_ridge_lambda_auto_as_unimplemented() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("auto_reject");
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
        "auto",
        "--min-peptides",
        "2",
        "--min-pairs",
        "4",
    ]);
    assert!(
        !out.status.success(),
        "--ridge-lambda auto must be refused under --test msqrob, not silently treated as 0.0"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("not implemented") && stderr.contains("0.0"),
        "refusal must say auto is unimplemented and name the no-shrinkage value; got: {stderr}"
    );
}

#[test]
fn non_msqrob_run_warns_but_succeeds_on_explicit_ridge_lambda_auto() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("auto_warn");
    write_two_group_fixture(&input);
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
        "--ridge-lambda",
        "auto",
        "--min-pairs",
        "4",
    ]);
    assert!(
        out.status.success(),
        "the flag is inert outside msqrob and must not fail the run:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--ridge-lambda") && stderr.contains("ignored"),
        "an inert --ridge-lambda should say so on stderr; got: {stderr}"
    );
}

#[test]
fn sidecar_records_resolved_ridge_lambda_only_when_msqrob_ran() {
    // msqrob run: the resolved numeric penalty is recorded.
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("resolved");
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
    let j = read_sidecar(&output);
    assert_eq!(j["args"]["ridge-lambda"], "0.5");
    assert_eq!(
        j["args"]["ridge-lambda-resolved"], 0.5,
        "the executed penalty belongs in the sidecar as a number"
    );

    // Non-msqrob run: the shrinkage path never executed, so the resolved
    // value must be null rather than implying msqrob ran.
    let tmp2 = tempfile::tempdir().unwrap();
    let input2 = tmp2.path().join("canonical");
    let output2 = tmp2.path().join("unresolved");
    write_two_group_fixture(&input2);
    std::fs::create_dir_all(&output2).unwrap();
    let out2 = run_atman(&[
        "de",
        "--input-dir",
        input2.to_str().unwrap(),
        "--output-dir",
        output2.to_str().unwrap(),
        "--test",
        "welch-t",
        "--groups",
        "A-B",
        "--min-pairs",
        "4",
    ]);
    assert!(
        out2.status.success(),
        "welch-t run failed:\n{}",
        String::from_utf8_lossy(&out2.stderr)
    );
    let j2 = read_sidecar(&output2);
    assert!(
        j2["args"]["ridge-lambda-resolved"].is_null(),
        "a run that never entered the shrinkage path must record null, got {}",
        j2["args"]["ridge-lambda-resolved"]
    );
}
