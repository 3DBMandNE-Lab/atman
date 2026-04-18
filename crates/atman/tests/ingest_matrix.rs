use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn ingest_matrix_supports_proteins_rows_and_preserves_sample_covariates() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path();
    let out = input.join("out");
    std::fs::write(
        input.join("samples.tsv"),
        "\
sample_id\tsubject\tgroup\tage\n\
S1\tA\tControl\t30\n\
S2\tB\tCase\t42\n",
    )
    .unwrap();
    std::fs::write(
        input.join("matrix.tsv"),
        "\
Protein.Group\tGenes\tPanel\tS1\tS2\n\
P00001\tGENE1\tms\t10\t12\n\
P00002\tGENE2\tms\t8\t\n",
    )
    .unwrap();

    let output = run_atman(&[
        "ingest-matrix",
        "--matrix",
        input.join("matrix.tsv").to_str().unwrap(),
        "--samples",
        input.join("samples.tsv").to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--orientation",
        "proteins-rows",
        "--platform",
        "spectronaut_report",
        "--abundance-unit",
        "log2_intensity",
        "--sample-id-col",
        "sample_id",
        "--subject-id-col",
        "subject",
        "--condition-col",
        "group",
        "--assay-id-col",
        "Protein.Group",
        "--gene-col",
        "Genes",
        "--panel-col",
        "Panel",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let samples = std::fs::read_to_string(out.join("samples.tsv")).unwrap();
    assert!(samples.starts_with(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tage\n"
    ));
    assert!(samples.contains("S2\tB\tCase\t0\t\t2\t42\n"));

    let proteins = std::fs::read_to_string(out.join("proteins.tsv")).unwrap();
    assert!(proteins.contains("spectronaut_report\tP00001\t\tGENE1\tms\t\n"));

    let measurements = std::fs::read_to_string(out.join("measurements.tsv")).unwrap();
    assert!(measurements
        .contains("spectronaut_report\tS1\tP00001\tGENE1\tms\t10\t10\t10\tlog2_intensity"));
    assert!(measurements.contains(
        "spectronaut_report\tS2\tP00002\tGENE2\tms\t\t\t\tlog2_intensity\tPASS\tPASS\t\t0\t1"
    ));
    assert!(out.join("qc_measurements.tsv").exists());
}

#[test]
fn ingest_matrix_supports_samples_rows_and_log2_transform() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path();
    let out = input.join("out");
    std::fs::write(
        input.join("samples.tsv"),
        "\
sample_id\tcondition\n\
S1\tControl\n\
S2\tCase\n",
    )
    .unwrap();
    std::fs::write(
        input.join("proteins.tsv"),
        "\
assay\tgene\n\
A1\tGENE1\n\
A2\tGENE2\n",
    )
    .unwrap();
    std::fs::write(
        input.join("matrix.tsv"),
        "\
sample_id\tA1\tA2\n\
S1\t8\t0\n\
S2\t16\t4\n",
    )
    .unwrap();

    let output = run_atman(&[
        "ingest-matrix",
        "--matrix",
        input.join("matrix.tsv").to_str().unwrap(),
        "--samples",
        input.join("samples.tsv").to_str().unwrap(),
        "--proteins",
        input.join("proteins.tsv").to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--orientation",
        "samples-rows",
        "--platform",
        "diann_report",
        "--abundance-unit",
        "log2_diann_pg_quantity",
        "--condition-col",
        "condition",
        "--assay-id-col",
        "assay",
        "--gene-col",
        "gene",
        "--log2-transform",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let measurements = std::fs::read_to_string(out.join("measurements.tsv")).unwrap();
    assert!(measurements.contains("diann_report\tS1\tA1\tGENE1\t\t8\t3\t3\tlog2_diann_pg_quantity"));
    assert!(measurements.contains(
        "diann_report\tS1\tA2\tGENE2\t\t0\t\t\tlog2_diann_pg_quantity\tPASS\tPASS\t\t0\t1"
    ));

    let validate = run_atman(&["validate", "--input-dir", out.to_str().unwrap()]);
    assert!(
        validate.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&validate.stderr)
    );
}

#[test]
fn ingest_matrix_fails_when_no_sample_columns_match() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path();
    let out = input.join("out");
    std::fs::write(
        input.join("samples.tsv"),
        "\
sample_id\tcondition\n\
S1\tControl\n",
    )
    .unwrap();
    std::fs::write(
        input.join("matrix.tsv"),
        "\
assay\tgene\tX1\n\
A1\tGENE1\t1\n",
    )
    .unwrap();

    let output = run_atman(&[
        "ingest-matrix",
        "--matrix",
        input.join("matrix.tsv").to_str().unwrap(),
        "--samples",
        input.join("samples.tsv").to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--platform",
        "maxquant_lfq",
        "--abundance-unit",
        "log2_lfq",
        "--condition-col",
        "condition",
        "--assay-id-col",
        "assay",
        "--gene-col",
        "gene",
    ]);

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("no matrix columns matched sample IDs")
    );
}
