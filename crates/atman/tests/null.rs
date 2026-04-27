use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_null_fixture(dir: &std::path::Path, paired: bool) {
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
spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n\
spectronaut_report\tP00002\tP00002\tGENE2\tspectronaut\t\n",
    )
    .unwrap();

    let rows = if paired {
        vec![
            ("A1", "P00001", "GENE1", "13.0", 1),
            ("B1", "P00001", "GENE1", "10.0", 2),
            ("A2", "P00001", "GENE1", "15.0", 3),
            ("B2", "P00001", "GENE1", "11.0", 4),
            ("A3", "P00001", "GENE1", "16.0", 5),
            ("B3", "P00001", "GENE1", "14.0", 6),
            ("A1", "P00002", "GENE2", "8.0", 7),
            ("B1", "P00002", "GENE2", "8.5", 8),
            ("A2", "P00002", "GENE2", "9.0", 9),
            ("B2", "P00002", "GENE2", "9.5", 10),
            ("A3", "P00002", "GENE2", "10.0", 11),
            ("B3", "P00002", "GENE2", "10.5", 12),
        ]
    } else {
        vec![
            ("C1", "P00001", "GENE1", "10.0", 1),
            ("C2", "P00001", "GENE1", "11.0", 2),
            ("C3", "P00001", "GENE1", "12.0", 3),
            ("K1", "P00001", "GENE1", "13.0", 4),
            ("K2", "P00001", "GENE1", "14.5", 5),
            ("K3", "P00001", "GENE1", "16.0", 6),
            ("C1", "P00002", "GENE2", "8.0", 7),
            ("C2", "P00002", "GENE2", "9.0", 8),
            ("C3", "P00002", "GENE2", "10.0", 9),
            ("K1", "P00002", "GENE2", "8.5", 10),
            ("K2", "P00002", "GENE2", "9.5", 11),
            ("K3", "P00002", "GENE2", "10.5", 12),
        ]
    };
    let mut measurements = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    for (sample, assay, gene, value, order) in rows {
        measurements.push_str(&format!(
            "spectronaut_report\t{sample}\t{assay}\t{gene}\tspectronaut\t{value}\t{value}\t{value}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
        ));
    }
    std::fs::write(dir.join("measurements.tsv"), measurements).unwrap();
}

#[test]
fn null_unpaired_permutation_is_seeded() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let out1 = tmp.path().join("null1");
    let out2 = tmp.path().join("null2");
    write_null_fixture(&input, false);

    for out in [&out1, &out2] {
        let output = run_atman(&[
            "null",
            "--input-dir",
            input.to_str().unwrap(),
            "--output-dir",
            out.to_str().unwrap(),
            "--groups",
            "Case-Control",
            "--test",
            "welch-t",
            "--n",
            "30",
            "--seed",
            "99",
            "--min-pairs",
            "2",
        ]);
        assert!(
            output.status.success(),
            "stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let summary1 = std::fs::read_to_string(out1.join("null_summary.tsv")).unwrap();
    let summary2 = std::fs::read_to_string(out2.join("null_summary.tsv")).unwrap();
    let empirical1 = std::fs::read_to_string(out1.join("empirical_p.tsv")).unwrap();
    let empirical2 = std::fs::read_to_string(out2.join("empirical_p.tsv")).unwrap();
    assert_eq!(summary1, summary2);
    assert_eq!(empirical1, empirical2);
    assert!(summary1.contains("Case-Control\twelch-t\t0.05\t"));
    assert!(empirical1.contains("Case-Control\tspectronaut\tP00001\tGENE1"));
}

#[test]
fn null_paired_sign_flip_reports_matched_subjects() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let out = tmp.path().join("null_paired");
    write_null_fixture(&input, true);

    let output = run_atman(&[
        "null",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--groups",
        "Case-Control",
        "--test",
        "paired-t",
        "--n",
        "20",
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
    let summary = std::fs::read_to_string(out.join("null_summary.tsv")).unwrap();
    let empirical = std::fs::read_to_string(out.join("empirical_p.tsv")).unwrap();
    assert!(summary.contains("Case-Control\tpaired-t\t0.05\t"));
    assert!(empirical.contains("Case-Control\tspectronaut\tP00001\tGENE1\tP00001\t3\t"));
}

#[test]
fn null_paired_refuses_unmatched_design() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    write_null_fixture(&input, false);

    let output = run_atman(&[
        "null",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        tmp.path().join("out").to_str().unwrap(),
        "--groups",
        "Case-Control",
        "--test",
        "paired-t",
        "--n",
        "10",
        "--min-pairs",
        "2",
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("insufficient matched subjects"));
}
