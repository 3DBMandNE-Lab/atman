//! Regression: `atman de --test paired-t` must accept both the new canonical
//! `--paired-by subject_id` (default) and the legacy `--paired-by participant`
//! alias, and reject other values with an informative error.
//!
//! Prior to the fix in this commit, the CLI hardcoded acceptance of only
//! `participant`, even though the canonical sample schema's pairing field is
//! `subject_id` and the actual pairing logic reads `subject_id` regardless of
//! the flag's value. The check was disconnected from what the code did.

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
    for (sample, value, order) in [
        ("A1", 12.0_f64, 1usize),
        ("B1", 10.0, 2),
        ("A2", 13.0, 3),
        ("B2", 10.0, 4),
        ("A3", 15.0, 5),
        ("B3", 11.0, 6),
        ("A4", 14.0, 7),
        ("B4", 11.0, 8),
    ] {
        measurements.push_str(&format!(
            "spectronaut_report\t{sample}\tP00001\tGENE1\tspectronaut\t{value}\t{value}\t{value}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
        ));
    }
    std::fs::write(dir.join("qc_measurements.tsv"), measurements).unwrap();
}

fn run_paired_t(input_dir: &std::path::Path, paired_by: Option<&str>) -> Output {
    let out = input_dir.join("out");
    let _ = std::fs::remove_dir_all(&out);
    let in_path = input_dir.to_str().unwrap().to_string();
    let out_path = out.to_str().unwrap().to_string();
    let mut args = vec![
        "de",
        "--input-dir",
        &in_path,
        "--output-dir",
        &out_path,
        "--test",
        "paired-t",
        "--groups",
        "Case-Control",
        "--min-pairs",
        "3",
    ];
    if let Some(pb) = paired_by {
        args.push("--paired-by");
        args.push(pb);
    }
    run_atman(&args)
}

#[test]
fn paired_t_works_with_default_paired_by() {
    let tmp = tempfile::tempdir().unwrap();
    write_paired_fixture(tmp.path());
    let out = run_paired_t(tmp.path(), None);
    assert!(
        out.status.success(),
        "default --paired-by must succeed; stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(tmp.path().join("out/de_results.tsv").exists());
}

#[test]
fn paired_t_works_with_explicit_subject_id() {
    let tmp = tempfile::tempdir().unwrap();
    write_paired_fixture(tmp.path());
    let out = run_paired_t(tmp.path(), Some("subject_id"));
    assert!(
        out.status.success(),
        "--paired-by subject_id must succeed; stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn paired_t_accepts_legacy_participant_alias() {
    let tmp = tempfile::tempdir().unwrap();
    write_paired_fixture(tmp.path());
    let out = run_paired_t(tmp.path(), Some("participant"));
    assert!(
        out.status.success(),
        "--paired-by participant (legacy alias) must still succeed; stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn paired_t_rejects_arbitrary_paired_by_with_helpful_error() {
    let tmp = tempfile::tempdir().unwrap();
    write_paired_fixture(tmp.path());
    let out = run_paired_t(tmp.path(), Some("not_a_real_column"));
    assert!(
        !out.status.success(),
        "unsupported --paired-by must be rejected; stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("subject_id"),
        "rejection message should mention `subject_id` so the user knows what to use; stderr:\n{stderr}"
    );
}
