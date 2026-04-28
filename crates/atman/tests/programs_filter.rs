use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn programs_filter_flags_interpretable_and_failed_programs() {
    let tmp = tempfile::tempdir().unwrap();
    let loadings = tmp.path().join("ica_loadings.tsv");
    let annotations = tmp.path().join("program_annotations.tsv");
    let output = tmp.path().join("programs_interpretable.tsv");

    std::fs::write(
        &loadings,
        "program\tgene_symbol\tloading\n\
         P1\tAQP4\t0.80\n\
         P1\tGFAP\t-0.30\n\
         P1\tKRT18\t0.10\n\
         P2\tKRT8\t0.70\n\
         P2\tKRT18\t0.60\n\
         P2\tKRT19\t0.50\n\
         P3\tSYN1\t0.20\n\
         P3\tSYP\t0.10\n\
         P4\tIGHG1\t0.90\n",
    )
    .unwrap();
    std::fs::write(
        &annotations,
        "program\tcategory\ttop_annotation\ttop_annotation_p_value\n\
         P1\tastrocyte\tGliosis\t0.001\n\
         P2\tcontaminant\tKeratinization\t0.001\n\
         P3\tsynaptic\tSynapse\t0.20\n",
    )
    .unwrap();

    let result = run_atman(&[
        "programs",
        "filter",
        "--loadings",
        loadings.to_str().unwrap(),
        "--annotations",
        annotations.to_str().unwrap(),
        "--min-annotation-pvalue",
        "0.05",
        "--min-top-loading",
        "0.5",
        "--max-contamination-fraction",
        "0.5",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let text = std::fs::read_to_string(output).unwrap();
    assert!(text.starts_with("program\tn_loadings\tmax_abs_loading\tcontamination_fraction_top_n"));
    let rows: Vec<Vec<&str>> = text
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 4);

    assert_eq!(rows[0][0], "P1");
    assert_eq!(rows[0][10], "1");
    assert_eq!(rows[0][11], "none");

    assert_eq!(rows[1][0], "P2");
    assert_eq!(rows[1][9], "0");
    assert_eq!(rows[1][10], "0");
    assert_eq!(rows[1][11], "contamination_signature");

    assert_eq!(rows[2][0], "P3");
    assert_eq!(rows[2][7], "0");
    assert_eq!(rows[2][8], "0");
    assert_eq!(rows[2][10], "0");
    assert_eq!(rows[2][11], "annotation_p_value;diffuse_loading");

    assert_eq!(rows[3][0], "P4");
    assert_eq!(rows[3][7], "0");
    assert_eq!(rows[3][10], "0");
    assert_eq!(rows[3][11], "missing_annotation");
}

#[test]
fn programs_filter_honors_custom_contamination_pattern() {
    let tmp = tempfile::tempdir().unwrap();
    let loadings = tmp.path().join("ica_loadings.tsv");
    let annotations = tmp.path().join("program_annotations.tsv");
    let output = tmp.path().join("programs_interpretable.tsv");

    std::fs::write(
        &loadings,
        "program\tgene_symbol\tloading\n\
         P1\tIGHG1\t0.90\n\
         P1\tIGHM\t0.80\n\
         P1\tIGKC\t0.70\n\
         P2\tKRT8\t0.90\n\
         P2\tKRT18\t0.80\n\
         P3\tAQP4\t0.80\n\
         P3\tGFAP\t0.60\n",
    )
    .unwrap();
    std::fs::write(
        &annotations,
        "program\tcategory\ttop_annotation\ttop_annotation_p_value\n\
         P1\tplasma_cell\tImmunoglobulin\t0.001\n\
         P2\tcontaminant\tKeratinization\t0.001\n\
         P3\tastrocyte\tGliosis\t0.001\n",
    )
    .unwrap();

    // Override default pattern with immunoglobulin regex. Keratin programs
    // (P2) should pass; immunoglobulin program (P1) should now fail.
    let result = run_atman(&[
        "programs",
        "filter",
        "--loadings",
        loadings.to_str().unwrap(),
        "--annotations",
        annotations.to_str().unwrap(),
        "--min-annotation-pvalue",
        "0.05",
        "--min-top-loading",
        "0.5",
        "--max-contamination-fraction",
        "0.5",
        "--contamination-pattern",
        "(?i)^IG[HKL]",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let text = std::fs::read_to_string(output).unwrap();
    let rows: Vec<Vec<&str>> = text
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 3);

    // P1 (IGHG1/IGHM/IGKC) is now flagged by the immunoglobulin pattern.
    assert_eq!(rows[0][0], "P1");
    assert_eq!(rows[0][9], "0");
    assert_eq!(rows[0][11], "contamination_signature");

    // P2 (KRT8/KRT18) is no longer flagged — keratin doesn't match `^IG[HKL]`.
    assert_eq!(rows[1][0], "P2");
    assert_eq!(rows[1][9], "1");
    assert_eq!(rows[1][10], "1");
    assert_eq!(rows[1][11], "none");

    assert_eq!(rows[2][0], "P3");
    assert_eq!(rows[2][10], "1");
    assert_eq!(rows[2][11], "none");
}

#[test]
fn programs_filter_rejects_invalid_contamination_pattern() {
    let tmp = tempfile::tempdir().unwrap();
    let loadings = tmp.path().join("ica_loadings.tsv");
    let annotations = tmp.path().join("program_annotations.tsv");
    let output = tmp.path().join("programs_interpretable.tsv");

    std::fs::write(&loadings, "program\tgene_symbol\tloading\nP1\tAQP4\t0.8\n").unwrap();
    std::fs::write(
        &annotations,
        "program\tcategory\ttop_annotation\ttop_annotation_p_value\nP1\tx\ty\t0.001\n",
    )
    .unwrap();

    let result = run_atman(&[
        "programs",
        "filter",
        "--loadings",
        loadings.to_str().unwrap(),
        "--annotations",
        annotations.to_str().unwrap(),
        "--contamination-pattern",
        "[unterminated",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("invalid --contamination-pattern"),
        "stderr was:\n{stderr}"
    );
}
