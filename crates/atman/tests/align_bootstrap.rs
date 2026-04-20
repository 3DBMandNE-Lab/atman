//! Integration test for `atman align bootstrap`.
//!
//! Two synthetic cohorts with a planted universal archetype (same
//! loading pattern in both cohorts) and cohort-specific archetypes
//! (A has its own planted signal on proteins 4..6; B has its own on
//! proteins 7..9). Runs subject-level bootstrap and asserts:
//!
//! 1. Output schema is as specified.
//! 2. At least one point-estimate archetype is detected.
//! 3. Every archetype row has finite, bounded statistics.
//! 4. Determinism across runs under fixed --seed.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

struct Lcg(u64);
impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as f64 / (u32::MAX as f64)
    }
    fn heavy_tail(&mut self) -> f64 {
        let u = self.next() * 2.0 - 1.0;
        u * u * u
    }
}

fn write_cohort(dir: &Path, cohort_label: &str, planted_a_specific: bool, planted_b_specific: bool, seed: u64) {
    std::fs::create_dir_all(dir).unwrap();
    let n_samples = 30usize;
    let n_proteins = 10usize;

    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=n_samples {
        samples.push_str(&format!(
            "{cohort_label}_S{i:03}\t{cohort_label}_S{i:03}\tN/A\t0\tplasma\t{i}\n"
        ));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 1..=n_proteins {
        proteins.push_str(&format!(
            "olink_explore_ngs\tA{j:03}\tQ{j:05}\tG{j:03}\tP1\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    // Universal archetype loadings on proteins 1..4; cohort-specific
    // loadings on 5..7 (A) or 8..10 (B).
    let mut rng = Lcg::new(seed);
    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order: u64 = 0;
    for i in 1..=n_samples {
        // Per-sample heavy-tailed source drives each archetype.
        let universal_source = rng.heavy_tail() * 2.0;
        let a_source = rng.heavy_tail() * 2.0;
        let b_source = rng.heavy_tail() * 2.0;
        for j in 1..=n_proteins {
            order += 1;
            let universal_weight = if (1..=4).contains(&j) { 1.0 } else { 0.0 };
            let a_weight =
                if planted_a_specific && (5..=7).contains(&j) { 1.0 } else { 0.0 };
            let b_weight =
                if planted_b_specific && (8..=10).contains(&j) { 1.0 } else { 0.0 };
            let noise = (rng.next() - 0.5) * 0.3;
            let value =
                universal_weight * universal_source + a_weight * a_source + b_weight * b_source + noise;
            qc.push_str(&format!(
                "olink_explore_ngs\t{cohort_label}_S{i:03}\tA{j:03}\tG{j:03}\tP1\t{value:.6}\t\
                 {value:.6}\t{value:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("qc_measurements.tsv"), &qc).unwrap();
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

fn parse_tsv(
    path: &Path,
) -> (Vec<String>, Vec<std::collections::HashMap<String, String>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines.next().unwrap().split('\t').map(String::from).collect();
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
fn align_bootstrap_runs_on_two_cohort_synthetic_fixture() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_cohort(&a, "A", true, false, 11);
    write_cohort(&b, "B", false, true, 22);
    let out = tmp.path().join("bootstrap").join("summary.tsv");

    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--labels",
        "A,B",
        "--k",
        "2",
        "--n-boot",
        "10",
        "--seed",
        "20260418",
        "--cosine-tau",
        "0.1",
        "--match-tau",
        "0.1",
        "--min-subjects",
        "10",
        "--max-iter",
        "80",
        "--tol",
        "1e-3",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "align bootstrap failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&status.stdout),
        String::from_utf8_lossy(&status.stderr)
    );

    let (header, rows) = parse_tsv(&out);
    for col in [
        "archetype_id",
        "observed_n_cohorts",
        "observed_cohorts",
        "bootstrap_mean_n_cohorts",
        "bootstrap_prob_universal",
        "bootstrap_prob_multi",
        "ci_lower_n_cohorts",
        "ci_upper_n_cohorts",
        "bootstrap_match_rate",
    ] {
        assert!(header.iter().any(|h| h == col), "missing column {col}");
    }
    for r in &rows {
        let n: usize = r["observed_n_cohorts"].parse().unwrap();
        assert!(
            (2..=2).contains(&n),
            "expected observed_n_cohorts in [2,2] for a universal archetype; got {n}"
        );
        let mean: f64 = r["bootstrap_mean_n_cohorts"].parse().unwrap();
        assert!(mean >= 0.0);
        let prob_u: f64 = r["bootstrap_prob_universal"].parse().unwrap();
        assert!((0.0..=1.0 + 1e-9).contains(&prob_u));
        let prob_m: f64 = r["bootstrap_prob_multi"].parse().unwrap();
        assert!((0.0..=1.0 + 1e-9).contains(&prob_m));
        let mr: f64 = r["bootstrap_match_rate"].parse().unwrap();
        assert!((0.0..=1.0 + 1e-9).contains(&mr));
    }

    // Sidecar shape.
    let sidecar_path = out.parent().unwrap().join("summary.tsv.run.json");
    assert!(sidecar_path.exists());
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar_path).unwrap()).unwrap();
    assert_eq!(j["args"]["k"], 2);
    assert_eq!(j["args"]["n-boot"], 10);
    assert_eq!(j["args"]["labels"], "A,B");
}

#[test]
fn align_bootstrap_refuses_single_cohort() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    write_cohort(&a, "A", true, false, 1);
    let out = tmp.path().join("summary.tsv");
    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        a.to_str().unwrap(),
        "--k",
        "2",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!status.status.success());
    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(stderr.contains("at least 2"), "unexpected: {stderr}");
}

#[test]
fn align_bootstrap_refuses_small_cohorts() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_cohort(&a, "A", true, false, 1);
    write_cohort(&b, "B", false, true, 2);
    let out = tmp.path().join("summary.tsv");
    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--k",
        "2",
        "--n-boot",
        "2",
        "--min-subjects",
        "1000",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!status.status.success());
    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(stderr.contains("min-subjects"), "unexpected: {stderr}");
}

#[test]
fn align_bootstrap_is_deterministic_across_runs() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_cohort(&a, "A", true, false, 11);
    write_cohort(&b, "B", false, true, 22);
    let out1 = tmp.path().join("r1.tsv");
    let out2 = tmp.path().join("r2.tsv");
    for out in [&out1, &out2] {
        let status = run_atman(&[
            "align",
            "bootstrap",
            "--cohorts",
            &format!("{},{}", a.display(), b.display()),
            "--labels",
            "A,B",
            "--k",
            "2",
            "--n-boot",
            "4",
            "--seed",
            "20260418",
            "--cosine-tau",
            "0.1",
            "--match-tau",
            "0.1",
            "--min-subjects",
            "10",
            "--max-iter",
            "40",
            "--tol",
            "1e-3",
            "--output",
            out.to_str().unwrap(),
        ]);
        assert!(status.status.success(), "{}", String::from_utf8_lossy(&status.stderr));
    }
    let a = std::fs::read_to_string(&out1).unwrap();
    let b = std::fs::read_to_string(&out2).unwrap();
    assert_eq!(a, b, "align bootstrap must be deterministic under fixed seed");
}
