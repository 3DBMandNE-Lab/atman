use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_fixture(dir: &std::path::Path) {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tage\n\
K1\tK1\tCase\t0\tcsf\t1\t60\n\
K2\tK2\tCase\t0\tcsf\t2\t61\n\
K3\tK3\tCase\t0\tcsf\t3\t62\n\
K4\tK4\tCase\t0\tcsf\t4\t63\n\
K5\tK5\tCase\t0\tcsf\t5\t64\n\
K6\tK6\tCase\t0\tcsf\t6\t65\n\
C1\tC1\tControl\t0\tcsf\t7\t60\n\
C2\tC2\tControl\t0\tcsf\t8\t61\n\
C3\tC3\tControl\t0\tcsf\t9\t62\n\
C4\tC4\tControl\t0\tcsf\t10\t63\n\
C5\tC5\tControl\t0\tcsf\t11\t64\n\
C6\tC6\tControl\t0\tcsf\t12\t65\n",
    )
    .unwrap();
    std::fs::write(
        dir.join("proteins.tsv"),
        "\
platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
diann_report\tP1\t\tDET1\tms\t\n\
diann_report\tP2\t\tABUND1\tms\t\n",
    )
    .unwrap();
    let mut buf = String::from("platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n");
    let mut order = 1u64;
    for sample in [
        "K1", "K2", "K3", "K4", "K5", "K6", "C1", "C2", "C3", "C4", "C5", "C6",
    ] {
        let case = sample.starts_with('K');
        let below = if sample == "C1" || (case && sample != "K6") {
            0
        } else {
            1
        };
        let value = if case { 10.0 } else { 9.8 };
        buf.push_str(&format!(
            "diann_report\t{sample}\tP1\tDET1\tms\t\t{value}\t{value}\tlog2_intensity\tPASS\tPASS\t\t{below}\t0\t\t\t{order}\n"
        ));
        order += 1;
        let value = if case { 12.0 } else { 10.0 };
        buf.push_str(&format!(
            "diann_report\t{sample}\tP2\tABUND1\tms\t\t{value}\t{value}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
        ));
        order += 1;
    }
    std::fs::write(dir.join("measurements.tsv"), buf).unwrap();
}

#[test]
fn detectability_writes_detection_and_abundance_layers() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_path = tmp.path().join("detectability.tsv");
    write_fixture(&input);

    let output = run_atman(&[
        "detectability",
        "--input-dir",
        input.to_str().unwrap(),
        "--groups",
        "Case-Control",
        "--design",
        "~ 1",
        "--min-samples",
        "5",
        "--min-detected",
        "3",
        "--min-missing",
        "3",
        "--output",
        output_path.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let table = std::fs::read_to_string(output_path).unwrap();
    assert!(table.contains("detection_or"));
    assert!(table.contains("DET1"));
    assert!(table.contains("ABUND1"));
    let summary = std::fs::read_to_string(tmp.path().join("detectability_summary.tsv")).unwrap();
    assert!(summary.contains("diagnostic"));
    assert!(summary.contains("detection_informative"));
}
