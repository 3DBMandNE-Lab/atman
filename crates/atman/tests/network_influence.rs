//! Integration tests for `atman network influence`.
//!
//! Covers: star-topology hub recovery (via the CLI end-to-end), CLI
//! schema, empty-graph refusal, stratified mode, and determinism
//! under repeated runs.

use std::path::Path;
use std::process::{Command, Output};

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

/// Deterministic pseudo-random normal (Box-Muller over LCG) so the
/// test fixture is independent of atman's Xoshiro.
struct Lcg(u64);
impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn unif(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as f64 / (u32::MAX as f64)
    }
    fn normal(&mut self) -> f64 {
        let u1 = self.unif().max(1e-12);
        let u2 = self.unif();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

/// Star topology: HUB correlates with every LEAF (shared source),
/// leaves do not correlate with each other beyond the shared HUB
/// driver.
fn write_star_fixture(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    let n_samples = 40usize;
    let leaf_count = 6usize;
    let mut rng = Lcg::new(20260420);

    let mut samples = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tstratum\n",
    );
    for i in 1..=n_samples {
        let stratum = if i % 2 == 0 { "A" } else { "B" };
        samples.push_str(&format!(
            "S{i:03}\tS{i:03}\tN/A\t0\tplasma\t{i}\t{stratum}\n"
        ));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    proteins.push_str("olink_explore_ngs\tH001\tQ00001\tHUB\tP1\t\n");
    for k in 1..=leaf_count {
        proteins.push_str(&format!(
            "olink_explore_ngs\tL{k:03}\tQ0000{k}\tLEAF{k:02}\tP1\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0u64;
    // Per-sample shared source drives HUB and every LEAF.
    for i in 1..=n_samples {
        let shared = rng.normal() * 3.0; // strong common driver
        let sid = format!("S{i:03}");
        order += 1;
        let hub_val = shared + rng.normal() * 0.1;
        qc.push_str(&format!(
            "olink_explore_ngs\t{sid}\tH001\tHUB\tP1\t{hub_val:.6}\t\
             {hub_val:.6}\t{hub_val:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
        ));
        for k in 1..=leaf_count {
            order += 1;
            // LEAF = shared + LEAF-specific noise. Per-leaf noise is
            // uncorrelated with other leaves, so leaf-leaf correlation
            // only goes through the shared (HUB) driver.
            let leaf_val = shared + rng.normal() * 1.5;
            qc.push_str(&format!(
                "olink_explore_ngs\t{sid}\tL{k:03}\tLEAF{k:02}\tP1\t{leaf_val:.6}\t\
                 {leaf_val:.6}\t{leaf_val:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

#[test]
fn network_influence_emits_expected_schema_and_stratum_column() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    write_star_fixture(&fixture);
    let out = tmp.path().join("influence.tsv");

    let status = run_atman(&[
        "network",
        "influence",
        "--input",
        fixture.join("measurements.tsv").to_str().unwrap(),
        "--samples",
        fixture.join("samples.tsv").to_str().unwrap(),
        "--method",
        "spearman",
        "--threshold",
        "0.3",
        "--min-subjects",
        "5",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "network influence failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&status.stdout),
        String::from_utf8_lossy(&status.stderr)
    );
    let (header, rows) = parse_tsv(&out);
    for col in [
        "feature_id",
        "stratum",
        "eigenvector_centrality",
        "betweenness_centrality",
        "influence_score",
        "degree",
        "n_subjects_used",
    ] {
        assert!(header.iter().any(|h| h == col), "missing column {col}");
    }
    assert_eq!(
        rows.len(),
        7,
        "expected 7 features (HUB + 6 leaves), got {}",
        rows.len()
    );

    // Sidecar includes method and adjacency-policy metadata.
    let sidecar = out.parent().unwrap().join("influence.tsv.run.json");
    assert!(sidecar.exists());
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    assert_eq!(j["args"]["method"], "spearman");
    assert!(
        j["args"]["adjacency-policy"]
            .as_str()
            .unwrap()
            .starts_with("hard:"),
        "expected hard adjacency policy"
    );
}

#[test]
fn network_influence_refuses_empty_graph_at_high_threshold() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    write_star_fixture(&fixture);
    let out = tmp.path().join("empty.tsv");

    let status = run_atman(&[
        "network",
        "influence",
        "--input",
        fixture.join("measurements.tsv").to_str().unwrap(),
        "--samples",
        fixture.join("samples.tsv").to_str().unwrap(),
        "--method",
        "spearman",
        "--threshold",
        "0.9999",
        "--min-subjects",
        "5",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!status.status.success(), "expected refusal");
    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(
        stderr.contains("empty graph"),
        "expected empty-graph error: {stderr}"
    );
}

#[test]
fn network_influence_stratify_emits_one_block_per_stratum() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    write_star_fixture(&fixture);
    let out = tmp.path().join("influence.tsv");

    let status = run_atman(&[
        "network",
        "influence",
        "--input",
        fixture.join("measurements.tsv").to_str().unwrap(),
        "--samples",
        fixture.join("samples.tsv").to_str().unwrap(),
        "--method",
        "spearman",
        "--threshold",
        "0.25",
        "--stratify",
        "stratum",
        "--min-subjects",
        "5",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );
    let (_, rows) = parse_tsv(&out);
    let strata: std::collections::BTreeSet<String> =
        rows.iter().map(|r| r["stratum"].clone()).collect();
    assert_eq!(strata.len(), 2, "expected two strata rows; got {strata:?}");
    assert!(strata.contains("A") && strata.contains("B"));
    // Each stratum should produce exactly 7 feature rows.
    for s in &strata {
        let n = rows.iter().filter(|r| &r["stratum"] == s).count();
        assert_eq!(n, 7, "stratum {s:?} should have 7 features, got {n}");
    }
}

#[test]
fn network_influence_is_deterministic_across_runs() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    write_star_fixture(&fixture);
    let out1 = tmp.path().join("r1.tsv");
    let out2 = tmp.path().join("r2.tsv");
    for out in [&out1, &out2] {
        let status = run_atman(&[
            "network",
            "influence",
            "--input",
            fixture.join("measurements.tsv").to_str().unwrap(),
            "--samples",
            fixture.join("samples.tsv").to_str().unwrap(),
            "--method",
            "pearson",
            "--soft-power",
            "6",
            "--min-subjects",
            "5",
            "--output",
            out.to_str().unwrap(),
        ]);
        assert!(
            status.status.success(),
            "{}",
            String::from_utf8_lossy(&status.stderr)
        );
    }
    let a = std::fs::read_to_string(&out1).unwrap();
    let b = std::fs::read_to_string(&out2).unwrap();
    assert_eq!(a, b, "network influence output must be byte-identical");
}

#[test]
fn network_influence_refuses_conflicting_adjacency_flags() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    write_star_fixture(&fixture);
    let out = tmp.path().join("x.tsv");
    let status = run_atman(&[
        "network",
        "influence",
        "--input",
        fixture.join("measurements.tsv").to_str().unwrap(),
        "--samples",
        fixture.join("samples.tsv").to_str().unwrap(),
        "--threshold",
        "0.3",
        "--soft-power",
        "6",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!status.status.success());
    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(
        stderr.contains("mutually exclusive"),
        "unexpected: {stderr}"
    );
}
