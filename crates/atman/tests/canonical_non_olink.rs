use std::process::Command;

fn run_atman(args: &[&str]) {
    let bin = env!("CARGO_BIN_EXE_atman");
    let status = Command::new(bin).args(args).status().expect("run atman");
    assert!(status.success(), "atman {:?} failed", args);
}

#[test]
fn de_accepts_non_olink_canonical_tsv() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output = tmp.path().join("output");
    std::fs::create_dir_all(&input).unwrap();

    std::fs::write(
        input.join("samples.tsv"),
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
C1\tC1\tControl\t0\tcsf\t1\n\
C2\tC2\tControl\t0\tcsf\t2\n\
C3\tC3\tControl\t0\tcsf\t3\n\
K1\tK1\tCase\t0\tcsf\t4\n\
K2\tK2\tCase\t0\tcsf\t5\n\
K3\tK3\tCase\t0\tcsf\t6\n",
    )
    .unwrap();

    std::fs::write(
        input.join("proteins.tsv"),
        "\
platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n",
    )
    .unwrap();

    std::fs::write(
        input.join("measurements.tsv"),
        "\
platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n\
spectronaut_report\tC1\tP00001\tGENE1\tspectronaut\t10.0\t10.0\t10.0\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t1\n\
spectronaut_report\tC2\tP00001\tGENE1\tspectronaut\t10.2\t10.2\t10.2\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t2\n\
spectronaut_report\tC3\tP00001\tGENE1\tspectronaut\t9.9\t9.9\t9.9\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t3\n\
spectronaut_report\tK1\tP00001\tGENE1\tspectronaut\t12.0\t12.0\t12.0\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t4\n\
spectronaut_report\tK2\tP00001\tGENE1\tspectronaut\t12.1\t12.1\t12.1\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t5\n\
spectronaut_report\tK3\tP00001\tGENE1\tspectronaut\t11.8\t11.8\t11.8\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t6\n",
    )
    .unwrap();

    run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "welch-t",
        "--groups",
        "Case-Control",
        "--min-pairs",
        "2",
    ]);

    let de_results = std::fs::read_to_string(output.join("de_results.tsv")).unwrap();
    assert!(de_results.contains("spectronaut\tP00001\tGENE1"));
    assert!(de_results.contains("Case-Control"));
}
