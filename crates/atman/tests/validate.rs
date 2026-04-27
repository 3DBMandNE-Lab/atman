use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_valid_input(dir: &std::path::Path, abundance_unit: &str, duplicate_measurement: bool) {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
C1\tC1\tControl\t0\tcsf\t1\n\
C2\tC2\tControl\t0\tcsf\t2\n\
K1\tK1\tCase\t0\tcsf\t3\n\
K2\tK2\tCase\t0\tcsf\t4\n",
    )
    .unwrap();
    std::fs::write(
        dir.join("proteins.tsv"),
        "\
platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n",
    )
    .unwrap();

    let mut measurements = format!(
        "\
platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n\
spectronaut_report\tC1\tP00001\tGENE1\tspectronaut\t10.0\t10.0\t10.0\t{abundance_unit}\tPASS\tPASS\t\t0\t0\t\t\t1\n\
spectronaut_report\tC2\tP00001\tGENE1\tspectronaut\t10.2\t10.2\t10.2\t{abundance_unit}\tPASS\tPASS\t\t0\t0\t\t\t2\n\
spectronaut_report\tK1\tP00001\tGENE1\tspectronaut\t12.0\t12.0\t12.0\t{abundance_unit}\tPASS\tPASS\t\t0\t0\t\t\t3\n\
spectronaut_report\tK2\tP00001\tGENE1\tspectronaut\t12.1\t12.1\t12.1\t{abundance_unit}\tPASS\tPASS\t\t0\t0\t\t\t4\n",
    );
    if duplicate_measurement {
        measurements.push_str(&format!(
            "spectronaut_report\tK2\tP00001\tGENE1\tspectronaut\t12.1\t12.1\t12.1\t{abundance_unit}\tPASS\tPASS\t\t0\t0\t\t\t5\n"
        ));
    }
    std::fs::write(dir.join("measurements.tsv"), measurements).unwrap();
}

#[test]
fn validate_accepts_clean_non_olink_dataset_and_writes_report() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let report = tmp.path().join("validate_report.tsv");
    write_valid_input(&input, "log2_intensity", false);

    let output = run_atman(&[
        "validate",
        "--input-dir",
        input.to_str().unwrap(),
        "--groups",
        "Case-Control",
        "--min-pairs",
        "2",
        "--report",
        report.to_str().unwrap(),
    ]);

    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("validate: 0 error(s), 0 warning(s)"));
    assert_eq!(
        std::fs::read_to_string(report).unwrap(),
        "severity\tcode\tmessage\n"
    );
}

#[test]
fn validate_accepts_reordered_samples_columns_and_de_still_runs() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_dir = tmp.path().join("out");
    std::fs::create_dir_all(&output_dir).unwrap();
    write_valid_input(&input, "log2_intensity", false);
    std::fs::write(
        input.join("samples.tsv"),
        "\
condition\tsample_type\tingest_order\tsample_id\tis_control\tsubject_id\n\
Control\tcsf\t1\tC1\t0\tC1\n\
Control\tcsf\t2\tC2\t0\tC2\n\
Case\tcsf\t3\tK1\t0\tK1\n\
Case\tcsf\t4\tK2\t0\tK2\n",
    )
    .unwrap();

    let validate = run_atman(&[
        "validate",
        "--input-dir",
        input.to_str().unwrap(),
        "--groups",
        "Case-Control",
        "--min-pairs",
        "2",
    ]);
    assert!(
        validate.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&validate.stderr)
    );

    let de = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--test",
        "welch-t",
        "--groups",
        "Case-Control",
        "--min-pairs",
        "2",
    ]);
    assert!(
        de.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&de.stdout),
        String::from_utf8_lossy(&de.stderr)
    );
    let results = std::fs::read_to_string(output_dir.join("de_results.tsv")).unwrap();
    assert!(results.lines().count() > 1, "expected computed DE rows");
}

#[test]
fn validate_rejects_missing_required_column() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    write_valid_input(&input, "log2_intensity", false);
    std::fs::write(
        input.join("proteins.tsv"),
        "\
assay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
P00001\tP00001\tGENE1\tspectronaut\t\n",
    )
    .unwrap();

    let output = run_atman(&["validate", "--input-dir", input.to_str().unwrap()]);

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("missing_column"));
    assert!(stderr.contains("platform"));
}

#[test]
fn validate_rejects_duplicate_measurement_keys() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    write_valid_input(&input, "log2_intensity", true);

    let output = run_atman(&["validate", "--input-dir", input.to_str().unwrap()]);

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("duplicate_measurement_key"));
}

#[test]
fn validate_strict_treats_warnings_as_failures() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    write_valid_input(&input, "vendor_area", false);

    let non_strict = run_atman(&["validate", "--input-dir", input.to_str().unwrap()]);
    assert!(
        non_strict.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&non_strict.stderr)
    );
    assert!(String::from_utf8_lossy(&non_strict.stderr).contains("unknown_abundance_unit"));

    let strict = run_atman(&[
        "validate",
        "--input-dir",
        input.to_str().unwrap(),
        "--strict",
    ]);
    assert!(!strict.status.success());
    assert!(String::from_utf8_lossy(&strict.stderr).contains("unknown_abundance_unit"));
}
