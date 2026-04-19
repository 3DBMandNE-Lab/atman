use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn ratio_tests_module_score_log_ratios() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = tmp.path().join("canonical");
    std::fs::create_dir(&input_dir).unwrap();
    let output = tmp.path().join("ratio.tsv");

    std::fs::write(
        input_dir.join("samples.tsv"),
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
         C1\tC1\tControl\t1\tbio\t1\n\
         C2\tC2\tControl\t1\tbio\t2\n\
         C3\tC3\tControl\t1\tbio\t3\n\
         C4\tC4\tControl\t1\tbio\t4\n\
         M1\tM1\tMS\t0\tbio\t5\n\
         M2\tM2\tMS\t0\tbio\t6\n\
         M3\tM3\tMS\t0\tbio\t7\n\
         M4\tM4\tMS\t0\tbio\t8\n",
    )
    .unwrap();
    std::fs::write(
        input_dir.join("measurements.tsv"),
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n\
         somascan\tC1\tintrathecal_IgV\tintrathecal_IgV\tmodule_score\t1.0\t1.0\t1.0\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t1\n\
         somascan\tC1\tplasma_IgG\tplasma_IgG\tmodule_score\t0.8\t0.8\t0.8\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t2\n\
         somascan\tC2\tintrathecal_IgV\tintrathecal_IgV\tmodule_score\t1.1\t1.1\t1.1\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t3\n\
         somascan\tC2\tplasma_IgG\tplasma_IgG\tmodule_score\t0.8\t0.8\t0.8\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t4\n\
         somascan\tC3\tintrathecal_IgV\tintrathecal_IgV\tmodule_score\t1.2\t1.2\t1.2\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t5\n\
         somascan\tC3\tplasma_IgG\tplasma_IgG\tmodule_score\t0.8\t0.8\t0.8\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t6\n\
         somascan\tC4\tintrathecal_IgV\tintrathecal_IgV\tmodule_score\t1.3\t1.3\t1.3\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t7\n\
         somascan\tC4\tplasma_IgG\tplasma_IgG\tmodule_score\t0.8\t0.8\t0.8\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t8\n\
         somascan\tM1\tintrathecal_IgV\tintrathecal_IgV\tmodule_score\t2.0\t2.0\t2.0\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t9\n\
         somascan\tM1\tplasma_IgG\tplasma_IgG\tmodule_score\t1.0\t1.0\t1.0\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t10\n\
         somascan\tM2\tintrathecal_IgV\tintrathecal_IgV\tmodule_score\t2.1\t2.1\t2.1\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t11\n\
         somascan\tM2\tplasma_IgG\tplasma_IgG\tmodule_score\t1.0\t1.0\t1.0\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t12\n\
         somascan\tM3\tintrathecal_IgV\tintrathecal_IgV\tmodule_score\t2.2\t2.2\t2.2\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t13\n\
         somascan\tM3\tplasma_IgG\tplasma_IgG\tmodule_score\t1.0\t1.0\t1.0\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t14\n\
         somascan\tM4\tintrathecal_IgV\tintrathecal_IgV\tmodule_score\t2.3\t2.3\t2.3\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t15\n\
         somascan\tM4\tplasma_IgG\tplasma_IgG\tmodule_score\t1.0\t1.0\t1.0\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t16\n",
    )
    .unwrap();

    let result = run_atman(&[
        "ratio",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--numerator",
        "modules.intrathecal_IgV",
        "--denominator",
        "modules.plasma_IgG",
        "--groups",
        "MS-Control",
        "--test",
        "wilcoxon",
        "--n-bootstrap",
        "200",
        "--seed",
        "7",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let text = std::fs::read_to_string(output).unwrap();
    assert!(text.starts_with("comparison\tnumerator\tdenominator\ttest"));
    let row: Vec<&str> = text.lines().nth(1).unwrap().split('\t').collect();
    assert_eq!(row[0], "MS-Control");
    assert_eq!(row[1], "intrathecal_IgV");
    assert_eq!(row[2], "plasma_IgG");
    assert_eq!(row[3], "wilcoxon");
    assert_eq!(row[5], "4");
    assert_eq!(row[6], "4");
    assert_eq!(row[9], "0.800000");
    assert_eq!(row[14], "200");
    assert!(row[12].parse::<f64>().unwrap() < 0.05);
    assert!(row[13].parse::<f64>().unwrap() > 0.95);
}
