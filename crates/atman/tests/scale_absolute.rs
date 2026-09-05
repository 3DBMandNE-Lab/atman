use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

const MEAS_HEADER: &str = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n";

fn measurement(sid: &str, gene: &str, v: f64, order: usize) -> String {
    format!("maxquant_lfq\t{sid}\t{gene}\t{gene}\tP\t{v}\t{v}\t{v}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n")
}

#[test]
fn scale_absolute_rescales_per_sample_and_drops_samples_without_total() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    std::fs::create_dir_all(&input).unwrap();
    // S1: total 400 ; S2: total 800 ; S3: no total ⇒ dropped.
    std::fs::write(input.join("samples.tsv"), "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\ttotal_protein\nS1\tS1\tC\t1\tbio\t1\t400\nS2\tS2\tC\t1\tbio\t2\t800\nS3\tS3\tC\t1\tbio\t3\t\n").unwrap();
    std::fs::write(input.join("proteins.tsv"), "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\nmaxquant_lfq\tGA\tGA\tGA\tP\t\nmaxquant_lfq\tGB\tGB\tGB\tP\t\n").unwrap();
    let mut m = String::from(MEAS_HEADER);
    // S1: GA = 3 (8), GB = 3 (8) ⇒ sum 16 ⇒ a' = a − 4 + log2(400)
    // S2: GA = 5 (32), GB = 4 (16) ⇒ sum 48 ⇒ a' = a − log2(48) + log2(800)
    m.push_str(&measurement("S1", "GA", 3.0, 1));
    m.push_str(&measurement("S1", "GB", 3.0, 2));
    m.push_str(&measurement("S2", "GA", 5.0, 3));
    m.push_str(&measurement("S2", "GB", 4.0, 4));
    m.push_str(&measurement("S3", "GA", 1.0, 5));
    std::fs::write(input.join("measurements.tsv"), m).unwrap();
    let out = tmp.path().join("abs");
    let r = run_atman(&[
        "scale",
        "absolute",
        "--input-dir",
        input.to_str().unwrap(),
        "--total-col",
        "total_protein",
        "--output-canonical-dir",
        out.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let text = std::fs::read_to_string(out.join("measurements.tsv")).unwrap();
    let rows: Vec<Vec<&str>> = text
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 4, "S3 dropped:\n{text}");
    let value = |sid: &str, gene: &str| -> f64 {
        rows.iter()
            .find(|r| r[1] == sid && r[3] == gene)
            .map(|r| r[6].parse::<f64>().unwrap())
            .unwrap()
    };
    let expect_s1 = 3.0 - 4.0 + 400f64.log2();
    assert!((value("S1", "GA") - expect_s1).abs() < 1e-12);
    assert!((value("S1", "GB") - expect_s1).abs() < 1e-12);
    assert!((value("S2", "GA") - (5.0 - 48f64.log2() + 800f64.log2())).abs() < 1e-12);
    assert!((value("S2", "GB") - (4.0 - 48f64.log2() + 800f64.log2())).abs() < 1e-12);
    // Linear sums after rescaling equal the totals.
    let s1_sum = value("S1", "GA").exp2() + value("S1", "GB").exp2();
    assert!((s1_sum - 400.0).abs() < 1e-9);
    assert!(out.join("samples.tsv").exists() && out.join("proteins.tsv").exists());
    let sidecar = std::fs::read_to_string(out.join("measurements.tsv.run.json")).unwrap();
    assert!(sidecar.contains("\"scale absolute\""));
    assert!(sidecar.contains("\"n_samples_dropped_without_total\": 1"));
    assert!(sidecar.contains("\"n_samples_scaled\": 2"));
}

#[test]
fn scale_absolute_rejects_missing_total_column() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    std::fs::create_dir_all(&input).unwrap();
    std::fs::write(input.join("samples.tsv"), "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\nS1\tS1\tC\t1\tbio\t1\n").unwrap();
    std::fs::write(
        input.join("proteins.tsv"),
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n",
    )
    .unwrap();
    std::fs::write(input.join("measurements.tsv"), MEAS_HEADER).unwrap();
    let r = run_atman(&[
        "scale",
        "absolute",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-canonical-dir",
        tmp.path().join("o").to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("total_protein"));
}
