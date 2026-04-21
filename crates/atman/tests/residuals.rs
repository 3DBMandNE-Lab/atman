use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn residuals_exports_long_and_wide_covariate_adjusted_values() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = tmp.path().join("canonical");
    std::fs::create_dir(&input_dir).unwrap();
    let long = tmp.path().join("residuals_long.tsv");
    let wide = tmp.path().join("residuals_wide.tsv");

    std::fs::write(
        input_dir.join("samples.tsv"),
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tage\tsex\n\
         S1\tS1\tControl\t1\tbio\t1\t10\tF\n\
         S2\tS2\tControl\t1\tbio\t2\t20\tM\n\
         S3\tS3\tCase\t0\tbio\t3\t30\tF\n\
         S4\tS4\tCase\t0\tbio\t4\t40\tM\n",
    )
    .unwrap();
    std::fs::write(
        input_dir.join("measurements.tsv"),
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n\
         somascan\tS1\tA1\tGENEA\tpanel\t11\t11\t11\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t1\n\
         somascan\tS2\tA1\tGENEA\tpanel\t19\t19\t19\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t2\n\
         somascan\tS3\tA1\tGENEA\tpanel\t29\t29\t29\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t3\n\
         somascan\tS4\tA1\tGENEA\tpanel\t41\t41\t41\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t4\n\
         somascan\tB1\tB1\tGENEB\tpanel\t21\t21\t21\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t5\n\
         somascan\tS1\tB1\tGENEB\tpanel\t21\t21\t21\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t6\n\
         somascan\tS2\tB1\tGENEB\tpanel\t41\t41\t41\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t7\n\
         somascan\tS3\tB1\tGENEB\tpanel\t61\t61\t61\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t8\n\
         somascan\tS4\tB1\tGENEB\tpanel\t81\t81\t81\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t9\n",
    )
    .unwrap();

    let result = run_atman(&[
        "residuals",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--design",
        "~ age",
        "--output",
        long.to_str().unwrap(),
        "--output-wide",
        wide.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let long_text = std::fs::read_to_string(long).unwrap();
    assert!(
        long_text.starts_with("sample_id\tsubject_id\tcondition\tassay_id\tgene_symbol\tresidual")
    );
    assert_eq!(long_text.lines().count(), 9);
    assert!(long_text.contains("S1\tS1\tControl\tB1\tGENEB\t0\n"));
    assert!(long_text.contains("S4\tS4\tCase\tB1\tGENEB\t0\n"));

    let wide_text = std::fs::read_to_string(wide).unwrap();
    assert!(wide_text.starts_with("assay_id\tgene_symbol\tS1\tS2\tS3\tS4\n"));
    assert!(wide_text.contains("B1\tGENEB\t0\t0\t0\t0\n"));
}

#[test]
fn residuals_wide_uses_sample_id_not_subject_id() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = tmp.path().join("canonical");
    std::fs::create_dir(&input_dir).unwrap();
    let wide = tmp.path().join("wide.tsv");

    std::fs::write(
        input_dir.join("samples.tsv"),
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tage\n\
         run_001\tSUBJ_A\tCtrl\t1\tbio\t1\t25\n\
         run_002\tSUBJ_B\tCtrl\t1\tbio\t2\t35\n\
         run_003\tSUBJ_C\tCase\t0\tbio\t3\t45\n\
         run_004\tSUBJ_D\tCase\t0\tbio\t4\t55\n",
    )
    .unwrap();
    std::fs::write(
        input_dir.join("measurements.tsv"),
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n\
         somascan\trun_001\tP1\tTP53\tpanel\t10\t10\t10\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t1\n\
         somascan\trun_002\tP1\tTP53\tpanel\t20\t20\t20\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t2\n\
         somascan\trun_003\tP1\tTP53\tpanel\t30\t30\t30\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t3\n\
         somascan\trun_004\tP1\tTP53\tpanel\t40\t40\t40\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t4\n",
    )
    .unwrap();

    let result = run_atman(&[
        "residuals",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--design",
        "~ age",
        "--output",
        tmp.path().join("long.tsv").to_str().unwrap(),
        "--output-wide",
        wide.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let header = std::fs::read_to_string(&wide)
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_string();
    let cols: Vec<&str> = header.split('\t').collect();
    assert_eq!(
        cols,
        vec!["assay_id", "gene_symbol", "run_001", "run_002", "run_003", "run_004"],
        "wide output must use sample_id columns, not subject_id"
    );
}
