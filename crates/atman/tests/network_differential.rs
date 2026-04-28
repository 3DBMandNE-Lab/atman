//! End-to-end tests for `atman network differential` across all three
//! modes (edge-pairwise, edge-summary, module). Builds a synthetic
//! two-cohort fixture with planted signal:
//!
//! - 6 proteins (P0..P5)
//! - Cohort C1: P0/P1/P2 are tightly co-expressed; P3/P4/P5 are noise.
//! - Cohort C2: P0/P1/P2 are decoupled (P0 randomized); P3/P4/P5 are
//!   tightly co-expressed.
//!
//! That gives a clean planted signal: the (P0, P1) and (P0, P2) edges
//! should sign-flip / decouple between cohorts; the (P3, P4) edges
//! should conserve high correlation in C2 but be noise in C1.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_cohort(dir: &Path, planted_block: usize, n_subjects: usize) {
    // planted_block = 0: P0/P1/P2 tight, P3/P4/P5 noise.
    // planted_block = 1: P3/P4/P5 tight, P0/P1/P2 noise.
    std::fs::create_dir_all(dir).unwrap();
    let mut samples = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n",
    );
    for i in 1..=n_subjects {
        samples.push_str(&format!("S{i:03}\tS{i:03}\tCase\t0\ttissue\t{i}\n"));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 0..6 {
        proteins.push_str(&format!("cptac_tmt_proteome\tP{j}\t\tP{j}\tsynth\t\n"));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    let mut measurements = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0u64;
    for i in 1..=n_subjects {
        // Deterministic pseudo-random per-subject signal driver.
        let s = ((i as f64) * 1.7).sin();
        let n = ((i as f64) * 0.31).cos() * 0.05;
        for j in 0..6 {
            order += 1;
            let value = if (planted_block == 0 && j < 3) || (planted_block == 1 && j >= 3) {
                // Tight block: all proteins follow the per-subject driver
                // with tiny per-protein noise.
                s + (j as f64) * 0.001 + n
            } else {
                // Noise block: per-protein independent oscillation, no
                // shared driver.
                ((i as f64) * (j as f64 + 1.7)).cos()
            };
            measurements.push_str(&format!(
                "cptac_tmt_proteome\tS{i:03}\tP{j}\tP{j}\tsynth\t{v:.6}\t{v:.6}\t{v:.6}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n",
                v = value
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), measurements).unwrap();
}

#[test]
fn network_differential_edge_summary_flags_planted_divergence() {
    let tmp = tempfile::tempdir().unwrap();
    let c1 = tmp.path().join("C1");
    let c2 = tmp.path().join("C2");
    write_cohort(&c1, 0, 30);
    write_cohort(&c2, 1, 30);
    let out = tmp.path().join("summary.tsv");

    let r = run_atman(&[
        "network",
        "differential",
        "--inputs",
        &format!("C1={},C2={}", c1.display(), c2.display()),
        "--output",
        out.to_str().unwrap(),
        "--mode",
        "edge-summary",
        "--min-overlap",
        "5",
    ]);
    assert!(
        r.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );

    let body = std::fs::read_to_string(&out).unwrap();
    let mut lines = body.lines();
    let header = lines.next().unwrap();
    assert!(header.starts_with("feature_a\tfeature_b\tn_cohorts"));
    assert!(header.ends_with("C1_corr\tC2_corr"), "per-cohort columns at end: {header}");
    let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
    assert!(!rows.is_empty(), "expected non-empty edge summary");

    // Find the (P0, P1) edge — should diverge between cohorts.
    let p0_p1 = rows
        .iter()
        .find(|r| (r[0] == "P0" && r[1] == "P1") || (r[0] == "P1" && r[1] == "P0"))
        .expect("P0-P1 edge missing");
    let c1_corr: f64 = p0_p1[p0_p1.len() - 2].parse().unwrap();
    let c2_corr: f64 = p0_p1[p0_p1.len() - 1].parse().unwrap();
    assert!(
        c1_corr.abs() > 0.7,
        "P0-P1 should be strongly correlated in C1, got {c1_corr}"
    );
    assert!(
        c2_corr.abs() < 0.5,
        "P0-P1 should decouple in C2, got {c2_corr}"
    );

    // Sidecar with input hashes.
    let sidecar = tmp.path().join("summary.tsv.run.json");
    assert!(sidecar.exists());
    let s = std::fs::read_to_string(&sidecar).unwrap();
    assert!(s.contains("network differential"));
    assert!(s.contains("\"C1_measurements\""));
    assert!(s.contains("\"C2_samples\""));
}

#[test]
fn network_differential_edge_pairwise_emits_significant_z_for_planted_flip() {
    let tmp = tempfile::tempdir().unwrap();
    let c1 = tmp.path().join("C1");
    let c2 = tmp.path().join("C2");
    write_cohort(&c1, 0, 30);
    write_cohort(&c2, 1, 30);
    let out = tmp.path().join("pairwise.tsv");

    let r = run_atman(&[
        "network",
        "differential",
        "--inputs",
        &format!("C1={},C2={}", c1.display(), c2.display()),
        "--output",
        out.to_str().unwrap(),
        "--mode",
        "edge-pairwise",
        "--min-overlap",
        "5",
    ]);
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));

    let body = std::fs::read_to_string(&out).unwrap();
    let mut lines = body.lines();
    let header = lines.next().unwrap();
    assert_eq!(
        header,
        "feature_a\tfeature_b\tcohort_a\tcohort_b\tn_a\tn_b\tcorr_a\tcorr_b\tz_diff\tp_value"
    );
    let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
    let p0_p1 = rows
        .iter()
        .find(|r| (r[0] == "P0" && r[1] == "P1") || (r[0] == "P1" && r[1] == "P0"))
        .expect("P0-P1 edge missing");
    let z_diff: f64 = p0_p1[8].parse().unwrap();
    let p_val: f64 = p0_p1[9].parse().unwrap();
    assert!(
        z_diff.abs() > 2.0,
        "Fisher-z should be |Z| > 2 for planted decoupling, got {z_diff}"
    );
    assert!(p_val < 0.05, "p-value should be < 0.05 for planted flip, got {p_val}");

    // Top row (sorted by |z_diff| desc) should be one of the planted
    // divergent edges (anything within the planted blocks).
    let top = &rows[0];
    let planted_edges: std::collections::HashSet<(&str, &str)> =
        [("P0", "P1"), ("P0", "P2"), ("P1", "P2"), ("P3", "P4"), ("P3", "P5"), ("P4", "P5")]
            .into_iter()
            .collect();
    let key = if top[0] < top[1] { (top[0], top[1]) } else { (top[1], top[0]) };
    assert!(
        planted_edges.contains(&key),
        "top |z_diff| edge should be a planted divergence, got {key:?}"
    );
}

#[test]
fn network_differential_module_mode_scores_rewiring() {
    let tmp = tempfile::tempdir().unwrap();
    let c1 = tmp.path().join("C1");
    let c2 = tmp.path().join("C2");
    write_cohort(&c1, 0, 30);
    write_cohort(&c2, 1, 30);
    let out = tmp.path().join("module.tsv");

    // Two modules: one over the C1-tight block, one over the C2-tight
    // block. Expect both to have non-zero rewiring (within-module
    // connectivity flips between cohorts).
    let sets = tmp.path().join("modules.tsv");
    std::fs::write(
        &sets,
        "set_name\tgene_symbol\nblock_lo\tP0\nblock_lo\tP1\nblock_lo\tP2\nblock_hi\tP3\nblock_hi\tP4\nblock_hi\tP5\n",
    )
    .unwrap();

    let r = run_atman(&[
        "network",
        "differential",
        "--inputs",
        &format!("C1={},C2={}", c1.display(), c2.display()),
        "--output",
        out.to_str().unwrap(),
        "--mode",
        "module",
        "--gene-sets",
        sets.to_str().unwrap(),
    ]);
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));

    let body = std::fs::read_to_string(&out).unwrap();
    let mut lines = body.lines();
    let header = lines.next().unwrap();
    assert!(header.starts_with("set_name\tset_size_declared"));
    assert!(header.ends_with("C1_connectivity\tC2_connectivity"));
    let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows.len(), 2);
    let by_set: std::collections::BTreeMap<&str, Vec<&str>> =
        rows.iter().map(|r| (r[0], r.clone())).collect();
    let lo = &by_set["block_lo"];
    let hi = &by_set["block_hi"];
    let lo_c1: f64 = lo[lo.len() - 2].parse().unwrap();
    let lo_c2: f64 = lo[lo.len() - 1].parse().unwrap();
    let hi_c1: f64 = hi[hi.len() - 2].parse().unwrap();
    let hi_c2: f64 = hi[hi.len() - 1].parse().unwrap();
    assert!(
        lo_c1 > lo_c2,
        "block_lo should be more connected in C1 (planted), got C1={lo_c1} C2={lo_c2}"
    );
    assert!(
        hi_c2 > hi_c1,
        "block_hi should be more connected in C2 (planted), got C1={hi_c1} C2={hi_c2}"
    );
    let lo_rewire: f64 = lo[6].parse().unwrap();
    assert!(lo_rewire > 0.0, "block_lo rewiring should be > 0, got {lo_rewire}");
}

#[test]
fn network_differential_module_mode_requires_gene_sets() {
    let tmp = tempfile::tempdir().unwrap();
    let c1 = tmp.path().join("C1");
    let c2 = tmp.path().join("C2");
    write_cohort(&c1, 0, 20);
    write_cohort(&c2, 1, 20);
    let out = tmp.path().join("module.tsv");
    let r = run_atman(&[
        "network",
        "differential",
        "--inputs",
        &format!("C1={},C2={}", c1.display(), c2.display()),
        "--output",
        out.to_str().unwrap(),
        "--mode",
        "module",
    ]);
    assert!(!r.status.success(), "missing --gene-sets must error");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(stderr.contains("--gene-sets"), "stderr: {stderr}");
}

#[test]
fn network_differential_requires_two_cohorts() {
    let tmp = tempfile::tempdir().unwrap();
    let c1 = tmp.path().join("C1");
    write_cohort(&c1, 0, 20);
    let out = tmp.path().join("summary.tsv");
    let r = run_atman(&[
        "network",
        "differential",
        "--inputs",
        c1.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--mode",
        "edge-summary",
    ]);
    assert!(!r.status.success(), "single-cohort input must error");
}
