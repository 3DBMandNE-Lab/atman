//! Instrument sanity check for the `de` subcommand: on the Dube heat
//! acclimation dataset, canonical heat-shock proteins must go UP in the
//! acute heat-stress comparisons (PT1-PR1 and PT2-PR2). If this fails,
//! the statistical pipeline is wrong in some fundamental way and no
//! downstream findings can be trusted.

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

fn run_atman(args: &[&str]) {
    let bin = env!("CARGO_BIN_EXE_atman");
    let status = Command::new(bin).args(args).status().expect("run atman");
    assert!(status.success(), "atman {:?} failed", args);
}

/// One row from de_results.tsv: (gene, mean_diff, p_value, bh_q).
type DeRow = (String, Option<f64>, Option<f64>, Option<f64>);

/// Read de_results.tsv rows for a given comparison.
fn load_de_rows(path: &std::path::Path, comparison: &str) -> Vec<DeRow> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .unwrap();
    let headers: Vec<String> = reader
        .headers()
        .unwrap()
        .iter()
        .map(|s| s.to_string())
        .collect();
    let col = |n: &str| headers.iter().position(|h| h == n).unwrap();
    let c_gene = col("gene_symbol");
    let c_comp = col("comparison");
    let c_mean_diff = col("mean_diff");
    let c_p = col("p_value");
    let c_q = col("bh_q");
    let mut out = Vec::new();
    for r in reader.records() {
        let r = r.unwrap();
        if &r[c_comp] != comparison {
            continue;
        }
        let parse_opt = |s: &str| s.parse::<f64>().ok();
        out.push((
            r[c_gene].to_string(),
            parse_opt(&r[c_mean_diff]),
            parse_opt(&r[c_p]),
            parse_opt(&r[c_q]),
        ));
    }
    out
}

#[test]
fn dube_de_heat_shock_sanity() {
    let root = repo_root();
    let data = root.join("example_data/dube_heat_2023");
    let raw1 = data.join("20212016_Dube_NPX_2021-11-30.csv");
    let raw2 = data.join("20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv");

    let tmp = tempfile::tempdir().unwrap();
    let tmp_path = tmp.path();

    let stager = root.join("adapters/generic/olink_explore_to_atman.py");
    let py_status = Command::new("python3")
        .arg(&stager)
        .arg("--output-dir")
        .arg(tmp_path)
        .arg(&raw1)
        .arg(&raw2)
        .status()
        .expect("run olink_explore_to_atman.py");
    assert!(py_status.success(), "olink_explore_to_atman.py failed");

    run_atman(&[
        "qc",
        "--input-dir",
        tmp_path.to_str().unwrap(),
        "--output-dir",
        tmp_path.to_str().unwrap(),
        "--rule",
        "mask-warn-fail",
    ]);

    run_atman(&[
        "de",
        "--input-dir",
        tmp_path.to_str().unwrap(),
        "--output-dir",
        tmp_path.to_str().unwrap(),
        "--test",
        "paired-t",
        "--paired-by",
        "participant",
        "--groups",
        "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2",
    ]);

    let de_path = tmp_path.join("de_results.tsv");
    assert!(de_path.exists(), "de_results.tsv not produced");

    // Canonical heat-shock proteins present on Olink Explore panels
    // (verified empirically from the reference filtered NPX headers).
    let canonical_hsps = ["HSPB1", "HSPA1A", "DNAJB1", "HSPG2"];

    for comparison in &["PT1-PR1", "PT2-PR2"] {
        let rows = load_de_rows(&de_path, comparison);
        let mut found_up = 0;
        let mut best_q: Option<f64> = None;
        for (gene, mean_diff, _p, q) in &rows {
            if canonical_hsps.iter().any(|h| h == gene) {
                // Direction check: mean_diff must be > 0 (heat-stress goes up).
                if let Some(d) = mean_diff {
                    if *d > 0.0 {
                        found_up += 1;
                    } else {
                        panic!(
                            "sanity violation: {} has mean_diff={} in {} (expected > 0 under heat stress)",
                            gene, d, comparison
                        );
                    }
                }
                if let Some(qv) = q {
                    best_q = Some(match best_q {
                        Some(m) => m.min(*qv),
                        None => *qv,
                    });
                }
            }
        }
        assert!(
            found_up >= 3,
            "only {} of {} canonical heat-shock proteins went up in {} (expected >= 3)",
            found_up,
            canonical_hsps.len(),
            comparison
        );
        assert!(
            best_q.map(|q| q < 0.2).unwrap_or(false),
            "best heat-shock protein q-value in {} is {:?}; expected at least one below 0.2",
            comparison,
            best_q
        );
        eprintln!(
            "sanity: {} — {}/{} canonical HSPs up, best q = {:.4}",
            comparison,
            found_up,
            canonical_hsps.len(),
            best_q.unwrap()
        );
    }
}
