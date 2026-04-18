use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_paired_fixture(dir: &std::path::Path) {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
A1\tP1\tCase\t0\tcsf\t1\n\
B1\tP1\tControl\t0\tcsf\t2\n\
A2\tP2\tCase\t0\tcsf\t3\n\
B2\tP2\tControl\t0\tcsf\t4\n\
A3\tP3\tCase\t0\tcsf\t5\n\
B3\tP3\tControl\t0\tcsf\t6\n\
A4\tP4\tCase\t0\tcsf\t7\n\
B4\tP4\tControl\t0\tcsf\t8\n",
    )
    .unwrap();
    write_proteins_and_measurements(
        dir,
        &[
            ("A1", 12.0, 1),
            ("B1", 10.0, 2),
            ("A2", 13.0, 3),
            ("B2", 10.0, 4),
            ("A3", 15.0, 5),
            ("B3", 11.0, 6),
            ("A4", 14.0, 7),
            ("B4", 11.0, 8),
        ],
    );
}

fn write_unpaired_fixture(dir: &std::path::Path) {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
C1\tC1\tControl\t0\tcsf\t1\n\
C2\tC2\tControl\t0\tcsf\t2\n\
C3\tC3\tControl\t0\tcsf\t3\n\
C4\tC4\tControl\t0\tcsf\t4\n\
K1\tK1\tCase\t0\tcsf\t5\n\
K2\tK2\tCase\t0\tcsf\t6\n\
K3\tK3\tCase\t0\tcsf\t7\n\
K4\tK4\tCase\t0\tcsf\t8\n",
    )
    .unwrap();
    write_proteins_and_measurements(
        dir,
        &[
            ("C1", 10.0, 1),
            ("C2", 11.0, 2),
            ("C3", 12.0, 3),
            ("C4", 13.0, 4),
            ("K1", 12.0, 5),
            ("K2", 14.0, 6),
            ("K3", 15.0, 7),
            ("K4", 17.0, 8),
        ],
    );
}

fn write_proteins_and_measurements(dir: &std::path::Path, rows: &[(&str, f64, usize)]) {
    std::fs::write(
        dir.join("proteins.tsv"),
        "\
platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n",
    )
    .unwrap();
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

fn read_first_row(path: &std::path::Path) -> (Vec<String>, Vec<String>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_string)
        .collect();
    let row = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_string)
        .collect();
    (header, row)
}

fn value(header: &[String], row: &[String], name: &str) -> String {
    let idx = header.iter().position(|h| h == name).unwrap();
    row[idx].clone()
}

#[test]
fn paired_de_appends_cohen_dz_ci_and_signed_rank() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_dir = tmp.path().join("out");
    write_paired_fixture(&input);

    let output = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--test",
        "paired-t",
        "--groups",
        "Case-Control",
        "--min-pairs",
        "2",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (header, row) = read_first_row(&output_dir.join("de_results.tsv"));
    assert_eq!(value(&header, &row, "effect_size_method"), "cohen_dz");
    assert_eq!(value(&header, &row, "wilcoxon_method"), "signed_rank");
    let expected_dz = 3.0 / (2.0_f64 / 3.0).sqrt();
    assert!(
        (value(&header, &row, "effect_size").parse::<f64>().unwrap() - expected_dz).abs() < 1e-9
    );
    assert!(!value(&header, &row, "ci_low").is_empty());
    assert!(!value(&header, &row, "wilcoxon_p").is_empty());
    assert_eq!(value(&header, &row, "median_diff"), "3");
    assert_eq!(value(&header, &row, "trimmed_mean_diff"), "3");
}

#[test]
fn welch_de_appends_hedges_g_ci_and_rank_sum() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_dir = tmp.path().join("out");
    write_unpaired_fixture(&input);

    let output = run_atman(&[
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
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (header, row) = read_first_row(&output_dir.join("de_results.tsv"));
    assert_eq!(value(&header, &row, "effect_size_method"), "hedges_g");
    assert_eq!(value(&header, &row, "wilcoxon_method"), "rank_sum");
    assert!(!value(&header, &row, "effect_size").is_empty());
    assert!(!value(&header, &row, "ci_high").is_empty());
    assert!(!value(&header, &row, "wilcoxon_p").is_empty());
    assert_eq!(value(&header, &row, "median_diff"), "3");
    assert_eq!(value(&header, &row, "trimmed_mean_diff"), "3");
}
