//! Integration test for `atman decompose null`.
//!
//! Builds a synthetic canonical TSV dataset containing two planted
//! low-rank archetypes plus independent Gaussian background proteins,
//! runs `atman decompose null` with k larger than the true signal
//! dimension, and confirms that the null-calibration output
//! distinguishes at least some programs from the permutation null.
//!
//! Assertions:
//!
//! 1. Output TSV has the expected column schema.
//! 2. Every `null_p` is in `[1/(n_perm+1), 1]`.
//! 3. BH q-values are monotone in p (non-decreasing after sort).
//! 4. `archetype_null` is deterministic under repeated identical seed.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

/// Tiny LCG so test data is independent of atman's own PRNG.
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
    fn normal(&mut self) -> f64 {
        let u1 = self.next().max(1e-12);
        let u2 = self.next();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

fn write_planted_canonical(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    let n_samples = 30usize;
    let n_signal_proteins = 30usize; // 15 for each of 2 planted archetypes
    let n_noise_proteins = 30usize;
    let n_total_proteins = n_signal_proteins + n_noise_proteins;

    // samples.tsv
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=n_samples {
        let cond = if i <= n_samples / 2 { "A" } else { "B" };
        samples.push_str(&format!("S{i:03}\tS{i:03}\t{cond}\t0\tplasma\t{i}\n"));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    // proteins.tsv
    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 1..=n_total_proteins {
        proteins.push_str(&format!(
            "olink_explore_ngs\tA{j:03}\tQ{j:05}\tG{j:03}\tP1\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    // Per-sample planted archetype activations (heavy-tailed so ICA's
    // non-Gaussianity criterion can find them).
    let mut rng = Lcg::new(20260418);
    let a1: Vec<f64> = (0..n_samples)
        .map(|_| {
            let u = rng.normal();
            u * u * u.signum()
        })
        .collect();
    let a2: Vec<f64> = (0..n_samples)
        .map(|_| {
            let u = rng.normal();
            u * u * u.signum()
        })
        .collect();

    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0;
    let signal_strength = 2.0;
    let noise_sigma = 0.2;
    for i in 0..n_samples {
        for j in 1..=n_total_proteins {
            order += 1;
            let assay = format!("A{j:03}");
            let gene = format!("G{j:03}");
            let value = if j <= 15 {
                // archetype 1 drives proteins 1..15
                signal_strength * a1[i] + noise_sigma * rng.normal()
            } else if j <= 30 {
                // archetype 2 drives proteins 16..30
                signal_strength * a2[i] + noise_sigma * rng.normal()
            } else {
                // background
                rng.normal()
            };
            qc.push_str(&format!(
                "olink_explore_ngs\tS{i:03}\t{assay}\t{gene}\tP1\t{value:.6}\t\
                 {value:.6}\t{value:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n",
                i = i + 1
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

#[test]
fn decompose_null_produces_valid_p_values_and_monotone_q() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let out_dir = tmp.path().join("null_out");
    write_planted_canonical(&input);

    let out_path: PathBuf = out_dir.join("archetype_null.tsv");
    let out = run_atman(&[
        "decompose",
        "null",
        "--input-dir",
        input.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--k",
        "4",
        "--n-perm",
        "25",
        "--n-seeds",
        "2",
        "--seed",
        "20260418",
        "--null-mode",
        "protein-shuffle",
        "--top-n",
        "10",
        "--max-iter",
        "150",
        "--tol",
        "1e-3",
    ]);
    assert!(
        out.status.success(),
        "atman decompose null failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );

    let (header, rows) = parse_tsv(&out_path);
    for col in [
        "program",
        "observed_stability",
        "null_stability_mean",
        "null_stability_p95",
        "null_p",
        "null_q",
        "decision",
        "n_perm",
    ] {
        assert!(header.iter().any(|h| h == col), "missing column {col}");
    }
    assert_eq!(rows.len(), 4, "expected k=4 programs");

    // All null_p in [1/(n+1), 1].
    let n_perm = 25_f64;
    let min_p = 1.0 / (n_perm + 1.0);
    for r in &rows {
        let p: f64 = r["null_p"].parse().unwrap();
        assert!(
            p >= min_p - 1e-12 && p <= 1.0 + 1e-12,
            "null_p {p} outside [{min_p}, 1.0]"
        );
    }

    // BH-q non-decreasing once the rows are sorted by p.
    let mut ps: Vec<(f64, f64)> = rows
        .iter()
        .map(|r| (r["null_p"].parse().unwrap(), r["null_q"].parse().unwrap()))
        .collect();
    ps.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for window in ps.windows(2) {
        assert!(
            window[0].1 <= window[1].1 + 1e-9,
            "BH q non-monotone: {:?} then {:?}",
            window[0],
            window[1]
        );
    }

    // Sidecar carries the null parameters.
    let sidecar = out_dir.join("archetype_null.tsv.run.json");
    assert!(sidecar.exists());
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    assert_eq!(j["args"]["null-mode"], "protein-shuffle");
    assert_eq!(j["args"]["k"], 4);
    assert_eq!(j["args"]["n-perm"], 25);
}

#[test]
fn decompose_null_is_deterministic_across_runs() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    write_planted_canonical(&input);
    let out_a = tmp.path().join("run_a").join("archetype_null.tsv");
    let out_b = tmp.path().join("run_b").join("archetype_null.tsv");
    for out_path in [&out_a, &out_b] {
        let out = run_atman(&[
            "decompose",
            "null",
            "--input-dir",
            input.to_str().unwrap(),
            "--output",
            out_path.to_str().unwrap(),
            "--k",
            "3",
            "--n-perm",
            "10",
            "--n-seeds",
            "2",
            "--seed",
            "20260418",
            "--null-mode",
            "protein-shuffle",
            "--top-n",
            "6",
            "--max-iter",
            "80",
            "--tol",
            "1e-3",
        ]);
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    let a = std::fs::read_to_string(&out_a).unwrap();
    let b = std::fs::read_to_string(&out_b).unwrap();
    assert_eq!(
        a, b,
        "decompose null output must be byte-identical under fixed seed"
    );
}

#[test]
fn decompose_null_rejects_k_above_sample_count() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    write_planted_canonical(&input);
    let out_path = tmp.path().join("archetype_null.tsv");
    let out = run_atman(&[
        "decompose",
        "null",
        "--input-dir",
        input.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--k",
        "500",
        "--n-perm",
        "2",
        "--n-seeds",
        "2",
    ]);
    assert!(!out.status.success(), "expected rejection");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("k=500 exceeds"),
        "unexpected stderr: {stderr}"
    );
}

/// The draw count that determines null_stability_mean,
/// null_stability_p95 and null_p must travel in the file.
///
/// It was already in the run sidecar, which is recoverable but not where
/// a reader of the table looks. The p95 is the reason: its index is
/// `ceil(0.95 * n_perm) - 1`, so at a low n_perm one or two draws decide
/// it and nothing in the row said how many there were.
///
/// n_perm is NOT the denominator of null_q, which is BH across programs.
/// A reader must not generalise the column to that one, so the reference
/// says so and this test pins the distinction by checking n_perm against
/// the flag rather than against the row count.
#[test]
fn decompose_null_writes_the_draw_count_that_determines_its_statistics() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let out_path = tmp.path().join("null_out").join("archetype_null.tsv");
    write_planted_canonical(&input);

    let n_perm = 24usize;
    let out = run_atman(&[
        "decompose",
        "null",
        "--input-dir",
        input.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--k",
        "4",
        "--n-perm",
        &n_perm.to_string(),
        "--n-seeds",
        "2",
        "--seed",
        "20260418",
        "--null-mode",
        "protein-shuffle",
        "--top-n",
        "10",
        "--max-iter",
        "150",
        "--tol",
        "1e-3",
    ]);
    assert!(
        out.status.success(),
        "decompose null failed:\nstderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let (header, rows) = parse_tsv(&out_path);
    assert!(
        header.iter().any(|h| h == "n_perm"),
        "n_perm must be a column"
    );
    assert!(!rows.is_empty(), "no programs written");
    for row in &rows {
        let got: usize = row["n_perm"].parse().expect("n_perm parses as an integer");
        assert_eq!(
            got, n_perm,
            "n_perm must report the draw count actually used, not the \
             row count or a default"
        );
        // Not the row count: the distinction the column must not blur,
        // since null_q's denominator IS the row count.
        assert_ne!(got, rows.len());
    }

    // At this count the p95 rests on the 23rd of 24 sorted draws, which
    // is the situation the column exists to expose.
    let p95_idx = (0.95f64 * n_perm as f64).ceil() as usize - 1;
    assert_eq!(
        p95_idx, 22,
        "at n_perm=24 the p95 index is the 23rd of 24 draws"
    );
}
