use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_module_fixture(dir: &std::path::Path, paired: bool) -> std::path::PathBuf {
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
            ("A1", "GENE1", "12.0", 1),
            ("A1", "GENE2", "14.0", 2),
            ("B1", "GENE1", "10.0", 3),
            ("B1", "GENE2", "12.0", 4),
            ("A2", "GENE1", "13.0", 5),
            ("A2", "GENE2", "15.0", 6),
            ("B2", "GENE1", "11.0", 7),
            ("B2", "GENE2", "13.0", 8),
            ("A3", "GENE1", "14.0", 9),
            ("A3", "GENE2", "16.0", 10),
            ("B3", "GENE1", "12.0", 11),
            ("B3", "GENE2", "14.0", 12),
        ]
    } else {
        vec![
            ("C1", "GENE1", "10.0", 1),
            ("C1", "GENE2", "12.0", 2),
            ("C2", "GENE1", "11.0", 3),
            ("C2", "GENE2", "13.0", 4),
            ("C3", "GENE1", "12.0", 5),
            ("C3", "GENE2", "14.0", 6),
            ("K1", "GENE1", "13.0", 7),
            ("K1", "GENE2", "15.0", 8),
            ("K2", "GENE1", "14.0", 9),
            ("K2", "GENE2", "16.0", 10),
            ("K3", "GENE1", "15.0", 11),
            ("K3", "GENE2", "17.0", 12),
        ]
    };
    let mut measurements = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    for (sample, gene, value, order) in rows {
        let assay = if gene == "GENE1" { "P00001" } else { "P00002" };
        measurements.push_str(&format!(
            "spectronaut_report\t{sample}\t{assay}\t{gene}\tspectronaut\t{value}\t{value}\t{value}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
        ));
    }
    std::fs::write(dir.join("measurements.tsv"), measurements).unwrap();

    let modules_tsv = dir.join("modules.tsv");
    std::fs::write(
        &modules_tsv,
        "\
module\tgene_symbol\n\
heat\tGENE1\n\
heat\tGENE2\n\
missing\tNOPE\n",
    )
    .unwrap();
    modules_tsv
}

#[test]
fn bootstrap_module_unpaired_is_seeded_and_reports_gene_coverage() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let out1 = tmp.path().join("module_boot1.tsv");
    let out2 = tmp.path().join("module_boot2.tsv");
    let modules = write_module_fixture(&input, false);

    for out in [&out1, &out2] {
        let output = run_atman(&[
            "bootstrap",
            "module",
            "--input-dir",
            input.to_str().unwrap(),
            "--modules-tsv",
            modules.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--groups",
            "Case-Control",
            "--test",
            "welch-t",
            "--n",
            "40",
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
    assert!(text1.contains("Case-Control\theat\t3\t3\t0\t2\t2\t"));
    assert!(text1.contains("\t40\t\n"));
    assert!(
        text1.contains("Case-Control\tmissing\t0\t0\t0\t1\t0\t\t\t\t\t\t0\tno_observed_genes\n")
    );
}

#[test]
fn bootstrap_module_paired_resamples_matched_subject_scores() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let out = tmp.path().join("module_paired.tsv");
    let modules = write_module_fixture(&input, true);

    let output = run_atman(&[
        "bootstrap",
        "module",
        "--input-dir",
        input.to_str().unwrap(),
        "--modules-tsv",
        modules.to_str().unwrap(),
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
    assert!(text.contains("Case-Control\theat\t3\t3\t3\t2\t2\t"));
    assert!(text.contains("\t1\t25\t\n"));
}
