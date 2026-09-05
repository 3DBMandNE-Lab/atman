//! Integration tests for `atman align bootstrap --decomposition nmf`.
//!
//! Mirrors the fixture style used by `align_bootstrap.rs`: per-sample
//! correlated "source" intensities drive planted programs shared across
//! both cohorts. [`write_nonneg_cohort`] keeps every value `>= 1.0` by
//! construction (a positive baseline plus non-negative source/noise
//! terms) so the default `--transform none` satisfies NMF's
//! non-negativity requirement without needing a transform.
//! [`write_signed_cohort`] produces ordinary signed data (can go
//! negative) to exercise the `--transform none` rejection path and the
//! `--transform exp2-clip` restoration path.

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
}

fn write_samples_and_proteins(dir: &Path, cohort_label: &str, n_samples: usize, n_proteins: usize) {
    std::fs::create_dir_all(dir).unwrap();
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
}

/// 12 samples × 20 assays. Proteins 1..=6 carry a "program A" source,
/// 7..=12 carry a "program B" source; both blocks are shared identically
/// across cohorts (independently drawn per-cohort per-sample source
/// magnitudes on the same protein blocks), so a correctly wired
/// alignment recovers both as universal archetypes. Proteins 13..=20 are
/// noise-only filler. Every value is `>= 1.0` by construction.
fn write_nonneg_cohort(dir: &Path, cohort_label: &str, seed: u64) {
    let n_samples = 12usize;
    let n_proteins = 20usize;
    write_samples_and_proteins(dir, cohort_label, n_samples, n_proteins);

    let mut rng = Lcg::new(seed);
    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order: u64 = 0;
    for i in 1..=n_samples {
        let source_a = rng.next() * 2.0;
        let source_b = rng.next() * 2.0;
        for j in 1..=n_proteins {
            order += 1;
            let a_weight = if (1..=6).contains(&j) { 1.0 } else { 0.0 };
            let b_weight = if (7..=12).contains(&j) { 1.0 } else { 0.0 };
            let noise = rng.next() * 0.1;
            let value = 1.0 + a_weight * source_a + b_weight * source_b + noise;
            qc.push_str(&format!(
                "olink_explore_ngs\t{cohort_label}_S{i:03}\tA{j:03}\tG{j:03}\tP1\t{value:.6}\t\
                 {value:.6}\t{value:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

/// Same 12 × 20 shape, but centered noise (`[-1, 1)`) with no positive
/// baseline: guaranteed to contain negative entries, for exercising the
/// non-negativity rejection and the `exp2-clip` restoration path.
fn write_signed_cohort(dir: &Path, cohort_label: &str, seed: u64) {
    let n_samples = 12usize;
    let n_proteins = 20usize;
    write_samples_and_proteins(dir, cohort_label, n_samples, n_proteins);

    let mut rng = Lcg::new(seed);
    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order: u64 = 0;
    for i in 1..=n_samples {
        let source_a = (rng.next() - 0.5) * 2.0;
        let source_b = (rng.next() - 0.5) * 2.0;
        for j in 1..=n_proteins {
            order += 1;
            let a_weight = if (1..=6).contains(&j) { 1.0 } else { 0.0 };
            let b_weight = if (7..=12).contains(&j) { 1.0 } else { 0.0 };
            let noise = (rng.next() - 0.5) * 0.2;
            let value = a_weight * source_a + b_weight * source_b + noise;
            qc.push_str(&format!(
                "olink_explore_ngs\t{cohort_label}_S{i:03}\tA{j:03}\tG{j:03}\tP1\t{value:.6}\t\
                 {value:.6}\t{value:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
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

fn read_sidecar(primary_output: &Path) -> serde_json::Value {
    let mut s = primary_output.as_os_str().to_os_string();
    s.push(".run.json");
    let sidecar_path = std::path::PathBuf::from(s);
    serde_json::from_str(&std::fs::read_to_string(&sidecar_path).unwrap()).unwrap()
}

/// Step 1: `--decomposition nmf` runs end to end on a non-negative
/// two-cohort fixture with two shared planted programs.
#[test]
fn align_bootstrap_nmf_runs_on_two_cohort_synthetic_fixture() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_nonneg_cohort(&a, "A", 11);
    write_nonneg_cohort(&b, "B", 22);
    let out = tmp.path().join("bootstrap").join("summary.tsv");

    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--labels",
        "A,B",
        "--decomposition",
        "nmf",
        "--k",
        "2",
        "--n-boot",
        "8",
        "--seed",
        "42",
        "--cosine-tau",
        "0.1",
        "--match-tau",
        "0.1",
        "--min-subjects",
        "8",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "align bootstrap --decomposition nmf failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&status.stdout),
        String::from_utf8_lossy(&status.stderr)
    );

    let (header, rows) = parse_tsv(&out);
    for col in [
        "archetype_id",
        "bootstrap_prob_universal",
        "bootstrap_prob_multi",
    ] {
        assert!(header.iter().any(|h| h == col), "missing column {col}");
    }
    assert!(!rows.is_empty(), "expected >= 1 archetype row");
    let mut any_positive = false;
    for r in &rows {
        let pu: f64 = r["bootstrap_prob_universal"].parse().unwrap();
        let pm: f64 = r["bootstrap_prob_multi"].parse().unwrap();
        assert!(
            (0.0..=1.0 + 1e-9).contains(&pu),
            "prob_universal out of range: {pu}"
        );
        assert!(
            (0.0..=1.0 + 1e-9).contains(&pm),
            "prob_multi out of range: {pm}"
        );
        if pm > 0.0 {
            any_positive = true;
        }
    }
    assert!(
        any_positive,
        "expected at least one archetype with bootstrap_prob_multi > 0; rows:\n{}",
        rows.iter()
            .map(|r| format!(
                "  id={} prob_universal={} prob_multi={}",
                r["archetype_id"], r["bootstrap_prob_universal"], r["bootstrap_prob_multi"]
            ))
            .collect::<Vec<_>>()
            .join("\n")
    );

    let sidecar = read_sidecar(&out);
    assert_eq!(sidecar["args"]["decomposition_method"], "nmf");
}

/// Step 2: two identical `--decomposition nmf` invocations produce
/// byte-identical output TSVs.
#[test]
fn align_bootstrap_nmf_is_deterministic_across_runs() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_nonneg_cohort(&a, "A", 11);
    write_nonneg_cohort(&b, "B", 22);
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
            "--decomposition",
            "nmf",
            "--beta-loss",
            "frobenius",
            "--init",
            "random",
            "--nmf-max-iter",
            "500",
            "--nmf-tol",
            "1e-5",
            "--transform",
            "none",
            "--k",
            "2",
            "--n-boot",
            "8",
            "--seed",
            "42",
            "--cosine-tau",
            "0.1",
            "--match-tau",
            "0.1",
            "--min-subjects",
            "8",
            "--output",
            out.to_str().unwrap(),
        ]);
        assert!(
            status.status.success(),
            "{}",
            String::from_utf8_lossy(&status.stderr)
        );
    }
    let x = std::fs::read_to_string(&out1).unwrap();
    let y = std::fs::read_to_string(&out2).unwrap();
    assert_eq!(
        x, y,
        "align bootstrap --decomposition nmf must be deterministic under a fixed seed"
    );
}

/// Step 3 (direct variant): `--decomposition ica` explicit must be
/// byte-identical to the flag omitted entirely — adding the new flag
/// with its default must not perturb today's ICA behavior. (The
/// existing `align_bootstrap.rs` suite is the other half of this
/// contract and is run unmodified alongside this file.)
#[test]
fn align_bootstrap_decomposition_ica_explicit_matches_omitted_flag() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_nonneg_cohort(&a, "A", 11);
    write_nonneg_cohort(&b, "B", 22);
    let out_explicit = tmp.path().join("explicit.tsv");
    let out_omitted = tmp.path().join("omitted.tsv");

    let base_args = |out: &Path| -> Vec<String> {
        vec![
            "align".into(),
            "bootstrap".into(),
            "--cohorts".into(),
            format!("{},{}", a.display(), b.display()),
            "--labels".into(),
            "A,B".into(),
            "--k".into(),
            "2".into(),
            "--n-boot".into(),
            "6".into(),
            "--seed".into(),
            "20260418".into(),
            "--cosine-tau".into(),
            "0.1".into(),
            "--match-tau".into(),
            "0.1".into(),
            "--min-subjects".into(),
            "8".into(),
            "--max-iter".into(),
            "60".into(),
            "--tol".into(),
            "1e-3".into(),
            "--output".into(),
            out.to_str().unwrap().into(),
        ]
    };

    let mut explicit_args = base_args(&out_explicit);
    explicit_args.push("--decomposition".into());
    explicit_args.push("ica".into());
    let explicit_refs: Vec<&str> = explicit_args.iter().map(|s| s.as_str()).collect();
    let status = run_atman(&explicit_refs);
    assert!(
        status.status.success(),
        "{}",
        String::from_utf8_lossy(&status.stderr)
    );

    let omitted_args = base_args(&out_omitted);
    let omitted_refs: Vec<&str> = omitted_args.iter().map(|s| s.as_str()).collect();
    let status = run_atman(&omitted_refs);
    assert!(
        status.status.success(),
        "{}",
        String::from_utf8_lossy(&status.stderr)
    );

    let x = std::fs::read_to_string(&out_explicit).unwrap();
    let y = std::fs::read_to_string(&out_omitted).unwrap();
    assert_eq!(
        x, y,
        "--decomposition ica must be byte-identical to the flag omitted"
    );

    for out in [&out_explicit, &out_omitted] {
        let sidecar = read_sidecar(out);
        assert_eq!(sidecar["args"]["decomposition_method"], "ica");
    }
}

/// NMF non-negativity constraint: `--transform none` (the default) on
/// signed input must fail loudly, the same way `decompose nmf` does.
#[test]
fn align_bootstrap_nmf_transform_none_rejects_negative_input() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_signed_cohort(&a, "A", 1);
    write_signed_cohort(&b, "B", 2);
    let out = tmp.path().join("summary.tsv");
    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--labels",
        "A,B",
        "--decomposition",
        "nmf",
        "--k",
        "2",
        "--n-boot",
        "2",
        "--seed",
        "1",
        "--min-subjects",
        "8",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!status.status.success());
    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(
        stderr.to_lowercase().contains("non-negative"),
        "expected a non-negative-input error, got: {stderr}"
    );
}

/// `--transform exp2-clip` restores non-negativity from signed input so
/// the same fixture that fails under `--transform none` succeeds here,
/// exercising every NMF-specific flag together.
#[test]
fn align_bootstrap_nmf_exp2_clip_transform_handles_signed_input() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_signed_cohort(&a, "A", 1);
    write_signed_cohort(&b, "B", 2);
    let out = tmp.path().join("summary.tsv");
    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--labels",
        "A,B",
        "--decomposition",
        "nmf",
        "--beta-loss",
        "frobenius",
        "--init",
        "random",
        "--nmf-max-iter",
        "500",
        "--nmf-tol",
        "1e-5",
        "--transform",
        "exp2-clip",
        "--transform-clamp",
        "6.0",
        "--k",
        "2",
        "--n-boot",
        "4",
        "--seed",
        "7",
        "--cosine-tau",
        "0.1",
        "--match-tau",
        "0.1",
        "--min-subjects",
        "8",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "align bootstrap with exp2-clip failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&status.stdout),
        String::from_utf8_lossy(&status.stderr)
    );
    let sidecar = read_sidecar(&out);
    assert_eq!(sidecar["args"]["decomposition_method"], "nmf");
    assert_eq!(sidecar["args"]["transform"], "exp2-clip");
    assert_eq!(sidecar["args"]["transform-clamp"], 6.0);
}

/// `--transform-clamp` is refused unless `--transform exp2-clip`,
/// mirroring `decompose nmf`.
#[test]
fn align_bootstrap_transform_clamp_requires_exp2_clip() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_nonneg_cohort(&a, "A", 1);
    write_nonneg_cohort(&b, "B", 2);
    let out = tmp.path().join("summary.tsv");
    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--decomposition",
        "nmf",
        "--transform",
        "none",
        "--transform-clamp",
        "3.0",
        "--k",
        "2",
        "--min-subjects",
        "8",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!status.status.success());
    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(stderr.contains("transform-clamp"), "unexpected: {stderr}");
}

/// Unknown `--decomposition` values are rejected with a clear error.
#[test]
fn align_bootstrap_rejects_unknown_decomposition() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_nonneg_cohort(&a, "A", 1);
    write_nonneg_cohort(&b, "B", 2);
    let out = tmp.path().join("summary.tsv");
    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--decomposition",
        "pca",
        "--k",
        "2",
        "--min-subjects",
        "8",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!status.status.success());
    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(stderr.contains("decomposition"), "unexpected: {stderr}");
}
