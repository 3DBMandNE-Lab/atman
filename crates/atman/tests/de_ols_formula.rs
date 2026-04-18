use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_ols_fixture(dir: &std::path::Path, confounded_sex: bool) {
    std::fs::create_dir_all(dir).unwrap();
    let samples = if confounded_sex {
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tage\tsex\n\
C1\tC1\tControl\t0\tcsf\t1\t30\tF\n\
C2\tC2\tControl\t0\tcsf\t2\t40\tF\n\
K1\tK1\tCase\t0\tcsf\t3\t31\tM\n\
K2\tK2\tCase\t0\tcsf\t4\t41\tM\n"
    } else {
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tage\tsex\n\
C1\tC1\tControl\t0\tcsf\t1\t30\tF\n\
C2\tC2\tControl\t0\tcsf\t2\t42\tM\n\
C3\tC3\tControl\t0\tcsf\t3\t35\tF\n\
C4\tC4\tControl\t0\tcsf\t4\t47\tM\n\
K1\tK1\tCase\t0\tcsf\t5\t32\tF\n\
K2\tK2\tCase\t0\tcsf\t6\t44\tM\n\
K3\tK3\tCase\t0\tcsf\t7\t37\tF\n\
K4\tK4\tCase\t0\tcsf\t8\t49\tM\n"
    };
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    std::fs::write(
        dir.join("proteins.tsv"),
        "\
platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n\
spectronaut_report\tP00002\tP00002\tGENE2\tspectronaut\t\n",
    )
    .unwrap();

    let sample_rows: Vec<(&str, &str, f64, &str)> = if confounded_sex {
        vec![
            ("C1", "Control", 30.0, "F"),
            ("C2", "Control", 40.0, "F"),
            ("K1", "Case", 31.0, "M"),
            ("K2", "Case", 41.0, "M"),
        ]
    } else {
        vec![
            ("C1", "Control", 30.0, "F"),
            ("C2", "Control", 42.0, "M"),
            ("C3", "Control", 35.0, "F"),
            ("C4", "Control", 47.0, "M"),
            ("K1", "Case", 32.0, "F"),
            ("K2", "Case", 44.0, "M"),
            ("K3", "Case", 37.0, "F"),
            ("K4", "Case", 49.0, "M"),
        ]
    };

    let mut measurements = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 1;
    for (sample, condition, age, sex) in sample_rows {
        let case = if condition == "Case" { 1.0 } else { 0.0 };
        let male = if sex == "M" { 1.0 } else { 0.0 };
        let gene1 = 5.0 + 2.0 * case + 0.10 * age + 0.4 * male;
        let gene2 = 2.0 - 1.0 * case + 0.05 * age - 0.2 * male;
        for (assay, gene, value) in [("P00001", "GENE1", gene1), ("P00002", "GENE2", gene2)] {
            measurements.push_str(&format!(
                "spectronaut_report\t{sample}\t{assay}\t{gene}\tspectronaut\t{value}\t{value}\t{value}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
            order += 1;
        }
    }
    std::fs::write(dir.join("qc_measurements.tsv"), measurements).unwrap();
}

#[test]
fn ols_design_reproduces_covariates_shortcut() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let legacy = tmp.path().join("legacy");
    let formula = tmp.path().join("formula");
    write_ols_fixture(&input, false);

    let legacy_output = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        legacy.to_str().unwrap(),
        "--test",
        "ols",
        "--groups",
        "Case-Control",
        "--covariates",
        "age,sex",
        "--min-pairs",
        "2",
    ]);
    assert!(
        legacy_output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&legacy_output.stderr)
    );

    let formula_output = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        formula.to_str().unwrap(),
        "--test",
        "ols",
        "--groups",
        "Case-Control",
        "--design",
        "~ condition + age + sex",
        "--contrast",
        "conditionCase",
        "--min-pairs",
        "2",
    ]);
    assert!(
        formula_output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&formula_output.stderr)
    );

    assert_eq!(
        std::fs::read_to_string(legacy.join("de_results.tsv")).unwrap(),
        std::fs::read_to_string(formula.join("de_results.tsv")).unwrap()
    );
    let design_report = std::fs::read_to_string(formula.join("de_design.tsv")).unwrap();
    assert!(design_report.contains("conditionCase,age,sexM"));
}

#[test]
fn ols_design_can_infer_two_condition_groups_from_contrast() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_dir = tmp.path().join("out");
    write_ols_fixture(&input, false);

    let output = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--test",
        "ols",
        "--design",
        "~ condition + age + sex",
        "--contrast",
        "conditionCase",
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
}

#[test]
fn ols_design_refuses_singular_model() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_dir = tmp.path().join("out");
    write_ols_fixture(&input, true);

    let output = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--test",
        "ols",
        "--groups",
        "Case-Control",
        "--design",
        "~ condition + sex",
        "--contrast",
        "conditionCase",
        "--min-pairs",
        "2",
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("singular design matrix"));
}
