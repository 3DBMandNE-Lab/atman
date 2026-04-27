//! Integration test for `atman de --test ensemble` on a 3-protein
//! synthetic fixture. UP is upregulated in condition B, DN is
//! downregulated in B, ST is null. Ensemble assigns:
//!   UP → VALIDATED (negative mean_diff because atman reports
//!        mean_a − mean_b)
//!   DN → VALIDATED (positive mean_diff)
//!   ST → INSUFFICIENT

use std::collections::HashMap;
use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
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

fn write_synthetic(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();

    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=12 {
        let cond = if i <= 6 { "A" } else { "B" };
        samples.push_str(&format!("S{i:02}\tS{i:02}\t{cond}\t0\tplasma\t{i}\n"));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    let proteins = "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
                    olink_explore_ngs\tA001\tQ01\tUP\tP1\t\n\
                    olink_explore_ngs\tA002\tQ02\tDN\tP1\t\n\
                    olink_explore_ngs\tA003\tQ03\tST\tP1\t\n";
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0;
    for i in 1..=12 {
        let is_b = i > 6;
        for (assay, gene, effect) in [
            ("A001", "UP", 1.5_f64),
            ("A002", "DN", -1.5),
            ("A003", "ST", 0.0),
        ] {
            order += 1;
            let v = 10.0 + if is_b { effect } else { 0.0 } + ((i * 7 + order) as f64).sin() * 0.05;
            qc.push_str(&format!(
                "olink_explore_ngs\tS{i:02}\t{assay}\t{gene}\tP1\t{v:.6}\t\
                 {v:.6}\t{v:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

#[test]
fn ensemble_validates_known_effects_and_insufficients_null() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("ensemble_out");
    write_synthetic(&input);
    std::fs::create_dir_all(&output).unwrap();

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "ensemble",
        "--groups",
        "A-B",
        "--ensemble-methods",
        "welch-t,limma",
        "--min-pairs",
        "3",
        "--trend",
        "false",
        "--robust",
        "false",
    ]);
    assert!(
        out.status.success(),
        "atman de --test ensemble failed:\nstderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let (header, rows) = parse_tsv(&output.join("de_ensemble.tsv"));
    for col in [
        "comparison",
        "panel",
        "assay_id",
        "gene_symbol",
        "uniprot",
        "n_applied",
        "n_significant",
        "n_sign_consistent",
        "majority_sign",
        "ensemble_p",
        "ensemble_q",
        "grade",
        "methods_applied",
        "methods_skipped",
    ] {
        assert!(header.iter().any(|h| h == col), "missing column {col}");
    }

    let by_gene: HashMap<&str, &HashMap<String, String>> = rows
        .iter()
        .map(|r| (r["gene_symbol"].as_str(), r))
        .collect();

    assert_eq!(
        by_gene["UP"]["grade"], "VALIDATED",
        "UP should be VALIDATED; got row {:?}",
        by_gene["UP"]
    );
    assert_eq!(by_gene["DN"]["grade"], "VALIDATED");
    assert_eq!(by_gene["ST"]["grade"], "INSUFFICIENT");
    // atman convention: mean_diff = mean_a − mean_b. UP is higher in B
    // (positive injected effect), so atman reports negative majority_sign.
    let up_sign: f64 = by_gene["UP"]["majority_sign"].parse().unwrap();
    let dn_sign: f64 = by_gene["DN"]["majority_sign"].parse().unwrap();
    assert!(
        up_sign < 0.0,
        "UP majority_sign should be negative (B > A), got {up_sign}"
    );
    assert!(
        dn_sign > 0.0,
        "DN majority_sign should be positive (A > B), got {dn_sign}"
    );

    // methods_applied covers both requested methods on every row.
    for gene in ["UP", "DN", "ST"] {
        let methods = &by_gene[gene]["methods_applied"];
        assert!(
            methods.contains("welch-t"),
            "missing welch-t for {gene}: {methods}"
        );
        assert!(
            methods.contains("limma"),
            "missing limma for {gene}: {methods}"
        );
    }

    // Per-method rows in de_results.tsv: expect 3 proteins × 2 methods = 6.
    let (_, per_method) = parse_tsv(&output.join("de_results.tsv"));
    assert_eq!(per_method.len(), 6, "expected 6 per-method rows");
    let method_values: std::collections::BTreeSet<String> =
        per_method.iter().map(|r| r["method"].clone()).collect();
    assert_eq!(
        method_values.iter().cloned().collect::<Vec<_>>(),
        vec!["limma", "welch-t"]
    );

    // Sidecar records ensemble metadata.
    let sidecar = output.join("de_results.tsv.run.json");
    assert!(sidecar.exists());
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    assert_eq!(j["args"]["test"], "ensemble");
    assert_eq!(j["args"]["ensemble-methods"], "welch-t,limma");
    // Resolved-value fields moved out of `args` into top-level extras
    // in schema v1 so user-intent vs runtime-outcome don't conflate.
    let applied: Vec<&str> = j["ensemble_methods_applied"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    assert_eq!(applied, vec!["welch-t", "limma"]);
    // `output_files` maps each written file's absolute path to its
    // sha256; both de_results.tsv and de_ensemble.tsv must appear.
    let outputs = j["output_files"]
        .as_object()
        .expect("output_files object in sidecar");
    let keys: Vec<&str> = outputs.keys().map(|s| s.as_str()).collect();
    assert!(
        keys.iter().any(|k| k.ends_with("de_results.tsv")),
        "de_results.tsv missing from output_files: {keys:?}"
    );
    assert!(
        keys.iter().any(|k| k.ends_with("de_ensemble.tsv")),
        "de_ensemble.tsv missing from output_files: {keys:?}"
    );
}

#[test]
fn ensemble_skips_msqrob_without_peptide_inputs_without_aborting() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("ensemble_skip");
    write_synthetic(&input);
    std::fs::create_dir_all(&output).unwrap();

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "ensemble",
        "--groups",
        "A-B",
        "--ensemble-methods",
        "welch-t,msqrob",
        "--min-pairs",
        "3",
        "--trend",
        "false",
        "--robust",
        "false",
    ]);
    assert!(
        out.status.success(),
        "ensemble should not abort on msqrob missing inputs; stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let sidecar = output.join("de_results.tsv.run.json");
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    // Resolved-value fields under top-level extras (v1 schema).
    let applied: Vec<&str> = j["ensemble_methods_applied"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    let skipped: Vec<&str> = j["ensemble_methods_skipped"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    assert_eq!(applied, vec!["welch-t"], "only welch-t should apply");
    assert!(
        skipped.iter().any(|s| s.starts_with("msqrob:")),
        "msqrob should be in skipped: {skipped:?}"
    );
}
