use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_mixed_fixture(dir: &std::path::Path, repeated: bool) {
    std::fs::create_dir_all(dir).unwrap();
    let samples = if repeated {
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
C1\tS1\tControl\t0\tcsf\t1\n\
K1\tS1\tCase\t0\tcsf\t2\n\
C2\tS2\tControl\t0\tcsf\t3\n\
K2\tS2\tCase\t0\tcsf\t4\n\
C3\tS3\tControl\t0\tcsf\t5\n\
K3\tS3\tCase\t0\tcsf\t6\n\
C4\tS4\tControl\t0\tcsf\t7\n\
K4\tS4\tCase\t0\tcsf\t8\n"
    } else {
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
C1\tC1\tControl\t0\tcsf\t1\n\
C2\tC2\tControl\t0\tcsf\t2\n\
K1\tK1\tCase\t0\tcsf\t3\n\
K2\tK2\tCase\t0\tcsf\t4\n"
    };
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    std::fs::write(
        dir.join("proteins.tsv"),
        "\
platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n",
    )
    .unwrap();

    let rows: Vec<(&str, f64)> = if repeated {
        vec![
            ("C1", 10.0),
            ("K1", 12.0),
            ("C2", 11.0),
            ("K2", 13.0),
            ("C3", 9.5),
            ("K3", 11.5),
            ("C4", 12.0),
            ("K4", 14.0),
        ]
    } else {
        vec![("C1", 10.0), ("C2", 11.0), ("K1", 12.0), ("K2", 13.0)]
    };
    let mut measurements = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    for (order, (sample, value)) in rows.into_iter().enumerate() {
        measurements.push_str(&format!(
            "spectronaut_report\t{sample}\tP00001\tGENE1\tspectronaut\t{value}\t{value}\t{value}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            order + 1
        ));
    }
    std::fs::write(dir.join("qc_measurements.tsv"), measurements).unwrap();
}

#[test]
fn mixed_random_intercept_runs_on_repeated_subjects() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_dir = tmp.path().join("out");
    write_mixed_fixture(&input, true);

    let output = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--test",
        "mixed",
        "--groups",
        "Case-Control",
        "--fixed",
        "condition",
        "--random",
        "1|subject_id",
        "--min-pairs",
        "2",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let de_results = std::fs::read_to_string(output_dir.join("de_results.tsv")).unwrap();
    assert!(de_results.contains("Case-Control"));
    let row = de_results.lines().nth(1).expect("one DE row");
    let fields: Vec<&str> = row.split('\t').collect();
    let mean_diff: f64 = fields[8].parse().expect("mean_diff");
    assert!(
        (mean_diff - 2.0).abs() < 1e-9,
        "unexpected mixed effect in:\n{de_results}"
    );
    let design = std::fs::read_to_string(output_dir.join("de_design.tsv")).unwrap();
    assert!(design.contains("conditionCase"));
}

#[test]
fn mixed_random_intercept_refuses_unrepeated_subjects() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_dir = tmp.path().join("out");
    write_mixed_fixture(&input, false);

    let output = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--test",
        "mixed",
        "--groups",
        "Case-Control",
        "--fixed",
        "condition",
        "--random",
        "1|subject_id",
        "--min-pairs",
        "2",
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("repeated subjects"));
}
