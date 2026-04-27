use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_fixture(dir: &std::path::Path, below_lod_fraction: f64) {
    // 4 samples x 4 proteins = 16 usable rows. Stamp `below_lod=1` on the
    // first ceil(16 * fraction) rows; the rest carry `below_lod=0`.
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
S1\tP1\tCase\t0\t\t1\n\
S2\tP1\tControl\t0\t\t2\n\
S3\tP2\tCase\t0\t\t3\n\
S4\tP2\tControl\t0\t\t4\n",
    )
    .unwrap();
    std::fs::write(
        dir.join("proteins.tsv"),
        "\
platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
diann_report\tA1\t\tG1\tms\t\n\
diann_report\tA2\t\tG2\tms\t\n\
diann_report\tA3\t\tG3\tms\t\n\
diann_report\tA4\t\tG4\tms\t\n",
    )
    .unwrap();
    let header = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n";
    let mut buf = String::from(header);
    let samples = ["S1", "S2", "S3", "S4"];
    let assays = ["A1", "A2", "A3", "A4"];
    let total_usable: usize = samples.len() * assays.len();
    let below_count: usize = (total_usable as f64 * below_lod_fraction).ceil() as usize;
    let mut k = 0usize;
    let mut order = 1u64;
    for s in samples {
        for a in assays {
            let val = 10.0 + (k as f64) * 0.1;
            let below = if k < below_count { 1 } else { 0 };
            buf.push_str(&format!(
                "diann_report\t{s}\t{a}\tG{}\tms\t\t{val}\t{val}\tlog2_diann_pg_quantity\tPASS\tPASS\t\t{below}\t0\t\t\t{order}\n",
                a.trim_start_matches('A'),
            ));
            k += 1;
            order += 1;
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &buf).unwrap();
    std::fs::write(dir.join("measurements.tsv"), &buf).unwrap();
}

#[test]
fn de_refuses_when_below_lod_fraction_exceeds_default_threshold() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    write_fixture(dir, 0.70);
    let output = run_atman(&[
        "de",
        "--input-dir",
        dir.to_str().unwrap(),
        "--output-dir",
        dir.join("out").to_str().unwrap(),
        "--test",
        "welch-t",
        "--groups",
        "Case-Control",
        "--min-pairs",
        "2",
    ]);
    assert!(!output.status.success(), "expected failure");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("below-LOD fraction"),
        "expected below-LOD error; stderr:\n{stderr}"
    );
}

#[test]
fn de_runs_when_below_lod_fraction_below_threshold() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    write_fixture(dir, 0.20);
    let output = run_atman(&[
        "de",
        "--input-dir",
        dir.to_str().unwrap(),
        "--output-dir",
        dir.join("out").to_str().unwrap(),
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
}

#[test]
fn de_runs_when_below_lod_fraction_high_but_allow_censored_set() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    write_fixture(dir, 0.70);
    let output = run_atman(&[
        "de",
        "--input-dir",
        dir.to_str().unwrap(),
        "--output-dir",
        dir.join("out").to_str().unwrap(),
        "--test",
        "welch-t",
        "--groups",
        "Case-Control",
        "--min-pairs",
        "2",
        "--allow-censored",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("WARNING below-LOD fraction"),
        "expected allow-censored warning; stderr:\n{stderr}"
    );
}
