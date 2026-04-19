use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn bootstrap_program_uses_signed_loadings_for_subject_scores() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    std::fs::create_dir(&input).unwrap();
    let loadings = tmp.path().join("ica_loadings.tsv");
    let output = tmp.path().join("program_bootstrap.tsv");

    std::fs::write(
        input.join("samples.tsv"),
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
         C1\tC1\tControl\t0\tbio\t1\n\
         C2\tC2\tControl\t0\tbio\t2\n\
         K1\tK1\tCase\t0\tbio\t3\n\
         K2\tK2\tCase\t0\tbio\t4\n",
    )
    .unwrap();
    std::fs::write(
        input.join("proteins.tsv"),
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
         somascan\tA1\tP1\tG1\tpanel\t\n\
         somascan\tA2\tP2\tG2\tpanel\t\n",
    )
    .unwrap();
    std::fs::write(
        input.join("measurements.tsv"),
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n\
         somascan\tC1\tA1\tG1\tpanel\t1\t1\t1\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t1\n\
         somascan\tC1\tA2\tG2\tpanel\t1\t1\t1\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t2\n\
         somascan\tC2\tA1\tG1\tpanel\t2\t2\t2\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t3\n\
         somascan\tC2\tA2\tG2\tpanel\t2\t2\t2\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t4\n\
         somascan\tK1\tA1\tG1\tpanel\t3\t3\t3\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t5\n\
         somascan\tK1\tA2\tG2\tpanel\t-1\t-1\t-1\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t6\n\
         somascan\tK2\tA1\tG1\tpanel\t4\t4\t4\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t7\n\
         somascan\tK2\tA2\tG2\tpanel\t0\t0\t0\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t8\n",
    )
    .unwrap();
    std::fs::write(
        &loadings,
        "program\tgene_symbol\tloading\n\
         P1\tG1\t1\n\
         P1\tG2\t-1\n",
    )
    .unwrap();

    let result = run_atman(&[
        "bootstrap",
        "program",
        "--input-dir",
        input.to_str().unwrap(),
        "--loadings",
        loadings.to_str().unwrap(),
        "--groups",
        "Case-Control",
        "--n",
        "100",
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
    assert!(text.starts_with("comparison\tprogram\tn_a\tn_b"));
    let row: Vec<&str> = text.lines().nth(1).unwrap().split('\t').collect();
    assert_eq!(row[0], "Case-Control");
    assert_eq!(row[1], "P1");
    assert_eq!(row[2], "2");
    assert_eq!(row[3], "2");
    assert_eq!(row[5], "2");
    assert_eq!(row[6], "2");
    assert_eq!(row[7], "2");
    assert_eq!(row[12], "100");
    assert!(row[11].parse::<f64>().unwrap() > 0.95);
}
