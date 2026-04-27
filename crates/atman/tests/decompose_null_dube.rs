//! Smoke test for `atman decompose null` against the bundled Dube
//! heat-acclimation cohort. Runs a short null calibration
//! (`n_perm=15`, `n_seeds=2`, `k=4`) and asserts output schema +
//! finite statistics. Not a claims test — it guards against regressions
//! in the end-to-end path on real multi-panel Olink data.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn parse_tsv(path: &Path) -> (Vec<String>, Vec<std::collections::HashMap<String, String>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    let rows = lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(h, c)| (h.clone(), c.to_string()))
                .collect()
        })
        .collect();
    (header, rows)
}

#[test]
fn decompose_null_runs_end_to_end_on_dube_cohort() {
    let root = repo_root();
    let data = root.join("example_data/dube_heat_2023");
    let raw1 = data.join("20212016_Dube_NPX_2021-11-30.csv");
    let raw2 = data.join("20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv");

    let tmp = tempfile::tempdir().unwrap();
    let canonical = tmp.path().join("canonical");
    let out = tmp.path().join("null_dube").join("archetype_null.tsv");

    let stager = root.join("adapters/generic/olink_explore_to_atman.py");
    let py_status = Command::new("python3")
        .arg(&stager)
        .arg("--output-dir")
        .arg(&canonical)
        .arg(&raw1)
        .arg(&raw2)
        .status()
        .expect("run olink_explore_to_atman.py");
    assert!(py_status.success(), "olink_explore_to_atman.py failed");
    assert!(run_atman(&[
        "qc",
        "--input-dir",
        canonical.to_str().unwrap(),
        "--output-dir",
        canonical.to_str().unwrap(),
    ])
    .status
    .success());

    let status = run_atman(&[
        "decompose",
        "null",
        "--input-dir",
        canonical.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--k",
        "4",
        "--n-perm",
        "15",
        "--n-seeds",
        "2",
        "--seed",
        "20260418",
        "--null-mode",
        "protein-shuffle",
        "--top-n",
        "5",
        "--max-iter",
        "120",
        "--tol",
        "1e-3",
        "--max-missing-fraction",
        "0.5",
        "--impute",
        "mean",
    ]);
    assert!(
        status.status.success(),
        "decompose null failed on Dube:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&status.stdout),
        String::from_utf8_lossy(&status.stderr)
    );

    let (header, rows) = parse_tsv(&out);
    for col in [
        "program",
        "observed_stability",
        "null_stability_mean",
        "null_stability_p95",
        "null_p",
        "null_q",
        "decision",
    ] {
        assert!(header.iter().any(|h| h == col), "missing column {col}");
    }
    assert_eq!(rows.len(), 4, "expected k=4 Dube programs");

    let n_perm_f = 15_f64;
    let min_p = 1.0 / (n_perm_f + 1.0);
    for r in &rows {
        let p: f64 = r["null_p"].parse().unwrap();
        let q: f64 = r["null_q"].parse().unwrap();
        let obs: f64 = r["observed_stability"].parse().unwrap();
        let mean: f64 = r["null_stability_mean"].parse().unwrap();
        assert!(
            p >= min_p - 1e-12 && p <= 1.0 + 1e-12,
            "null_p {p} on {} outside [{min_p},1]",
            r["program"]
        );
        assert!(q.is_finite(), "null_q non-finite on {}", r["program"]);
        assert!(
            obs.is_finite() && (0.0..=1.0 + 1e-9).contains(&obs),
            "observed_stability {obs} outside [0,1] on {}",
            r["program"]
        );
        assert!(
            mean.is_finite() && (0.0..=1.0 + 1e-9).contains(&mean),
            "null_stability_mean {mean} outside [0,1] on {}",
            r["program"]
        );
    }
}
