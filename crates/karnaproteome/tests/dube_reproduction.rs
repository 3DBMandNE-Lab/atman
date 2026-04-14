//! End-to-end reproduction of Dube et al. Scientific Data 2023 Olink Explore
//! published filtered NPX and fold change files. Passing this test is the
//! definition of v0.1 done.
//!
//! - Filtered NPX: strict byte-exact cell-for-cell match (spec §8).
//! - Fold change:  |delta| ≤ 1e-4 numerical match (spec §8). Observed max
//!   delta in practice is ~1e-15 (IEEE 754 last-bit drift).

mod common;

use common::{diff_numeric, diff_strict};
use std::path::PathBuf;
use std::process::Command;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

fn run_karnaproteome(args: &[&str]) {
    let bin = env!("CARGO_BIN_EXE_karnaproteome");
    let status = Command::new(bin)
        .args(args)
        .status()
        .expect("run karnaproteome");
    assert!(status.success(), "karnaproteome {:?} failed", args);
}

const PANELS: &[&str] = &[
    "cardiometabolic",
    "cardiometabolic_ii",
    "inflammation",
    "inflammation_ii",
    "neurology",
    "neurology_ii",
    "oncology",
    "oncology_ii",
];

/// Numeric-diff tolerance for fold change. Spec §8 calls for 1e-4; we
/// observe ~1e-15 in practice. Keep the spec-level slack so this test
/// remains green against small numerical drift in future refactors, but
/// print the actual observed max delta so regressions are visible.
const FC_TOLERANCE: f64 = 1e-4;

#[test]
fn dube_reproduction_end_to_end() {
    let root = repo_root();
    let data = root.join("example_data/dube_heat_2023");
    let raw1 = data.join("20212016_Dube_NPX_2021-11-30.csv");
    let raw2 = data.join("20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv");
    assert!(raw1.exists(), "raw file 1 missing: {:?}", raw1);
    assert!(raw2.exists(), "raw file 2 missing: {:?}", raw2);

    let tmp = tempfile::tempdir().unwrap();
    let tmp_path = tmp.path();

    // Stage 1: ingest
    run_karnaproteome(&[
        "ingest",
        "--platform", "olink-explore-ngs",
        "--parser", "dube",
        "--output-dir", tmp_path.to_str().unwrap(),
        raw1.to_str().unwrap(),
        raw2.to_str().unwrap(),
    ]);

    // Stage 2: qc
    run_karnaproteome(&[
        "qc",
        "--input-dir", tmp_path.to_str().unwrap(),
        "--output-dir", tmp_path.to_str().unwrap(),
        "--rule", "dube",
    ]);

    // Stage 3: matrix (Dube-wide per-panel CSVs)
    run_karnaproteome(&[
        "matrix",
        "--input-dir", tmp_path.to_str().unwrap(),
        "--output-dir", tmp_path.to_str().unwrap(),
        "--format", "dube-wide",
        "--split-by", "panel",
    ]);

    // Stage 4: fold change
    run_karnaproteome(&[
        "fold-change",
        "--input-dir", tmp_path.to_str().unwrap(),
        "--output-dir", tmp_path.to_str().unwrap(),
        "--groups", "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1",
    ]);

    // Stage 5a: strict diff of filtered NPX files (byte-exact cell match)
    for panel in PANELS {
        let ours = tmp_path.join(format!("{}_npx.csv", panel));
        let reference = data.join(format!("filtered_npx/npx/{}_npx.csv", panel));
        assert!(ours.exists(), "missing ours: {:?}", ours);
        assert!(reference.exists(), "missing reference: {:?}", reference);
        let report = diff_strict(&ours, &reference);
        report.assert_clean();
    }

    // Stage 5b: numeric diff of fold change files (|delta| ≤ FC_TOLERANCE)
    let mut worst_fc_delta = 0.0_f64;
    for panel in PANELS {
        let ours = tmp_path.join(format!("{}_log2_fc.csv", panel));
        let reference = data.join(format!(
            "fold_changes/fold_changes/{}_log2_fc.csv",
            panel
        ));
        assert!(ours.exists(), "missing ours: {:?}", ours);
        assert!(reference.exists(), "missing reference: {:?}", reference);
        let report = diff_numeric(&ours, &reference, FC_TOLERANCE);
        if report.max_numeric_delta > worst_fc_delta {
            worst_fc_delta = report.max_numeric_delta;
        }
        report.assert_clean();
    }

    eprintln!(
        "dube_reproduction: all 16 reference files diff-clean. \
         Worst FC numeric delta: {:.2e} (tolerance: {:.0e})",
        worst_fc_delta, FC_TOLERANCE,
    );
}
