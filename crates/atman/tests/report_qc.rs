use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn report_qc_writes_summary_sample_protein_and_condition_tables() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_dir = tmp.path().join("report");
    std::fs::create_dir_all(&input).unwrap();

    std::fs::write(
        input.join("samples.tsv"),
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
C1\tC1\tControl\t0\tcsf\t1\n\
C2\tC2\tControl\t0\tcsf\t2\n\
K1\tK1\tCase\t0\tcsf\t3\n\
QC1\t\t\t1\tcontrol\t4\n",
    )
    .unwrap();
    std::fs::write(
        input.join("proteins.tsv"),
        "\
platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n\
spectronaut_report\tP00002\tP00002\tGENE2\tspectronaut\t\n",
    )
    .unwrap();
    std::fs::write(
        input.join("measurements.tsv"),
        "\
platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n\
spectronaut_report\tC1\tP00001\tGENE1\tspectronaut\t10.0\t10.0\t10.0\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t1\n\
spectronaut_report\tC2\tP00001\tGENE1\tspectronaut\t10.2\t10.2\t10.2\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t2\n\
spectronaut_report\tK1\tP00001\tGENE1\tspectronaut\t\t12.0\t12.0\tlog2_intensity\tWARN\tPASS\t\t0\t1\t\t\t3\n\
spectronaut_report\tC1\tP00002\tGENE2\tspectronaut\t8.0\t8.0\t8.0\tlog2_intensity\tPASS\tPASS\t9.0\t1\t0\t\t\t4\n",
    )
    .unwrap();

    let output = run_atman(&[
        "report",
        "qc",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--min-subjects",
        "1",
    ]);

    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let summary = std::fs::read_to_string(output_dir.join("qc_summary.tsv")).unwrap();
    assert!(summary.contains("sample_count\t4\n"));
    assert!(summary.contains("protein_count\t2\n"));
    assert!(summary.contains("measurement_count\t4\n"));
    assert!(summary.contains("effective_measurement_count\t3\n"));
    assert!(summary.contains("qc_masked_count\t1\n"));
    assert!(summary.contains("below_lod_count\t1\n"));

    let sample_qc = std::fs::read_to_string(output_dir.join("sample_qc.tsv")).unwrap();
    assert!(sample_qc.contains("K1\tK1\tCase\t0\t1\t0\t1\t1\t1\t0\n"));

    let protein_qc = std::fs::read_to_string(output_dir.join("protein_qc.tsv")).unwrap();
    assert!(protein_qc.contains(
        "spectronaut_report\tP00001\tGENE1\tspectronaut\t3\t2\t1\t0.3333333333333333\t1\t0\t1\n"
    ));
    assert!(protein_qc
        .contains("spectronaut_report\tP00002\tGENE2\tspectronaut\t1\t1\t0\t0\t0\t1\t1\n"));

    let conditions = std::fs::read_to_string(output_dir.join("condition_counts.tsv")).unwrap();
    assert!(conditions.contains("Control\t2\t2\t2\t2\t3\n"));

    let missingness = std::fs::read_to_string(output_dir.join("missingness_summary.tsv")).unwrap();
    assert!(missingness.contains("all\tall\t4\t2\t8\t4\t4\t0.5\t1\t1\t3\n"));
    assert!(missingness.contains("condition\tControl\t2\t2\t4\t3\t1\t0.25\t0\t1\t3\n"));
    assert!(missingness.contains("condition\tCase\t1\t2\t2\t1\t1\t0.5\t1\t0\t0\n"));
}
