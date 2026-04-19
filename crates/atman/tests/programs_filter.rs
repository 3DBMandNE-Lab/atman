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
        "--max-keratin-fraction",
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
    assert!(text.starts_with("program\tn_loadings\tmax_abs_loading"));
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
    assert_eq!(rows[1][11], "keratin_contamination");

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
