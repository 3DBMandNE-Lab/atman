use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_bootstrap_fixture(dir: &std::path::Path, paired: bool) {
    std::fs::create_dir_all(dir).unwrap();
    let samples = if paired {
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
A1\tP1\tCase\t0\tcsf\t1\n\
B1\tP1\tControl\t0\tcsf\t2\n\
A2\tP2\tCase\t0\tcsf\t3\n\
B2\tP2\tControl\t0\tcsf\t4\n\
A3\tP3\tCase\t0\tcsf\t5\n\
B3\tP3\tControl\t0\tcsf\t6\n"
    } else {
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
C1\tC1\tControl\t0\tcsf\t1\n\
C2\tC2\tControl\t0\tcsf\t2\n\
C3\tC3\tControl\t0\tcsf\t3\n\
K1\tK1\tCase\t0\tcsf\t4\n\
K2\tK2\tCase\t0\tcsf\t5\n\
K3\tK3\tCase\t0\tcsf\t6\n"
    };
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    std::fs::write(
        dir.join("proteins.tsv"),
        "\
platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n",
    )
    .unwrap();

    let rows = if paired {
        vec![
            ("A1", "12.0", 1),
            ("B1", "10.0", 2),
            ("A2", "13.0", 3),
            ("B2", "11.0", 4),
            ("A3", "14.0", 5),
            ("B3", "12.0", 6),
        ]
    } else {
        vec![
            ("C1", "10.0", 1),
            ("C2", "11.0", 2),
            ("C3", "12.0", 3),
            ("K1", "13.0", 4),
            ("K2", "14.0", 5),
            ("K3", "15.0", 6),
        ]
    };
    let mut measurements = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    for (sample, value, order) in rows {
        measurements.push_str(&format!(
            "spectronaut_report\t{sample}\tP00001\tGENE1\tspectronaut\t{value}\t{value}\t{value}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
        ));
    }
    std::fs::write(dir.join("qc_measurements.tsv"), measurements).unwrap();
}

#[test]
fn bootstrap_protein_unpaired_is_seeded_and_reports_effect_interval() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let out1 = tmp.path().join("boot1.tsv");
    let out2 = tmp.path().join("boot2.tsv");
    write_bootstrap_fixture(&input, false);

    for out in [&out1, &out2] {
        let output = run_atman(&[
            "bootstrap",
            "protein",
            "--input-dir",
            input.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--groups",
            "Case-Control",
            "--test",
            "welch-t",
            "--n",
            "50",
            "--seed",
            "123",
            "--min-pairs",
            "2",
        ]);
        assert!(
            output.status.success(),
            "stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let text1 = std::fs::read_to_string(out1).unwrap();
    let text2 = std::fs::read_to_string(out2).unwrap();
    assert_eq!(text1, text2);
    assert!(text1.contains("Case-Control\tspectronaut_report\tspectronaut\tP00001\tGENE1"));
    assert!(text1.contains("\t3\t3\t0\t3\t"));
    assert!(text1.contains("\t50\t\n"));
}

#[test]
fn bootstrap_protein_paired_resamples_matched_subjects() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let out = tmp.path().join("paired.tsv");
    write_bootstrap_fixture(&input, true);

    let output = run_atman(&[
        "bootstrap",
        "protein",
        "--input-dir",
        input.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--groups",
        "Case-Control",
        "--test",
        "paired-t",
        "--n",
        "25",
        "--seed",
        "7",
        "--min-pairs",
        "2",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("\t3\t3\t3\t2\t"));
    assert!(text.contains("\t1\t25\t\n"));
}
