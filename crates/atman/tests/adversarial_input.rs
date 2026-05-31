//! Adversarial / malformed input robustness.
//!
//! `atman` is a local analysis CLI that reads untrusted scientific data files.
//! Malformed input must produce a clean error and a non-zero exit — never a
//! panic, out-of-bounds, or hang. These tests codify the robustness verified in
//! the 1.1.0 security review (ragged rows rejected by the strict CSV reader; no
//! header-claimed-dimension allocation; non-finite handled as MNAR, not a
//! crash). They are the stable-toolchain complement to the `fuzz/` targets.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_valid_input(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
         C1\tC1\tControl\t0\tcsf\t1\n\
         C2\tC2\tControl\t0\tcsf\t2\n\
         K1\tK1\tCase\t0\tcsf\t3\n\
         K2\tK2\tCase\t0\tcsf\t4\n",
    )
    .unwrap();
    std::fs::write(
        dir.join("proteins.tsv"),
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
         spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n",
    )
    .unwrap();
    let hdr = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n";
    let mut m = String::from(hdr);
    for (s, a) in [
        ("C1", "10.0"),
        ("C2", "10.2"),
        ("K1", "12.0"),
        ("K2", "12.1"),
    ] {
        m.push_str(&format!(
            "spectronaut_report\t{s}\tP00001\tGENE1\tspectronaut\t{a}\t{a}\t{a}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t1\n"
        ));
    }
    std::fs::write(dir.join("measurements.tsv"), m).unwrap();
}

/// The load-bearing invariant: whatever the input, the process must not panic.
fn assert_no_panic(out: &Output, label: &str) {
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !stderr.contains("panicked") && !stderr.contains("RUST_BACKTRACE"),
        "{label}: process PANICKED on malformed input:\n{stderr}"
    );
}

fn validate(dir: &Path) -> Output {
    run_atman(&["validate", "--input-dir", dir.to_str().unwrap()])
}

#[test]
fn ragged_measurements_row_errors_without_panic() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    write_valid_input(dir);
    // Append a row with too few fields.
    let mut m = std::fs::read_to_string(dir.join("measurements.tsv")).unwrap();
    m.push_str("spectronaut_report\tC1\tP00001\tonly\tthree\textra\n");
    std::fs::write(dir.join("measurements.tsv"), m).unwrap();

    let out = validate(dir);
    assert_no_panic(&out, "ragged row");
    assert!(!out.status.success(), "ragged row should be rejected");
}

#[test]
fn missing_required_column_errors_without_panic() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    write_valid_input(dir);
    // Drop the `abundance` column from the header only (now a missing column).
    let m = std::fs::read_to_string(dir.join("measurements.tsv")).unwrap();
    let stripped = m.replacen("\tabundance\t", "\t", 1);
    std::fs::write(dir.join("measurements.tsv"), stripped).unwrap();

    let out = validate(dir);
    assert_no_panic(&out, "missing column");
    assert!(
        !out.status.success(),
        "missing required column should error"
    );
}

#[test]
fn binary_garbage_measurements_errors_without_panic() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    write_valid_input(dir);
    std::fs::write(
        dir.join("measurements.tsv"),
        [0u8, 159, 146, 150, 255, 0, 1, 2, 9],
    )
    .unwrap();

    let out = validate(dir);
    assert_no_panic(&out, "binary garbage");
    assert!(!out.status.success(), "binary garbage should error");
}

#[test]
fn truncated_midrow_measurements_errors_without_panic() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    write_valid_input(dir);
    let m = std::fs::read_to_string(dir.join("measurements.tsv")).unwrap();
    // Cut off mid-file (keep header + a partial row).
    let cut = &m[..m.len().saturating_sub(15)];
    std::fs::write(dir.join("measurements.tsv"), cut).unwrap();

    let out = validate(dir);
    assert_no_panic(&out, "truncated");
    // May or may not error depending on where the cut lands; the invariant is
    // simply: no panic.
}

#[test]
fn empty_and_header_only_files_do_not_panic() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    write_valid_input(dir);
    // Header-only measurements (zero data rows).
    let hdr = std::fs::read_to_string(dir.join("measurements.tsv"))
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_string();
    std::fs::write(dir.join("measurements.tsv"), format!("{hdr}\n")).unwrap();
    assert_no_panic(&validate(dir), "header-only");

    // Completely empty file.
    std::fs::write(dir.join("measurements.tsv"), "").unwrap();
    assert_no_panic(&validate(dir), "empty file");

    // Missing file entirely.
    std::fs::remove_file(dir.join("measurements.tsv")).unwrap();
    assert_no_panic(&validate(dir), "missing file");
}
