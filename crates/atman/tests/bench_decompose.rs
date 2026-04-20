//! Priority 6 integration test: `atman bench decompose`.
//!
//! Builds a tiny planted-archetype fixture on the fly (2 archetypes ×
//! 30 proteins × 40 subjects, well-separated heavy-tailed sources),
//! runs `atman bench decompose --tools atman` against it, and verifies:
//!
//! 1. Output TSV schema includes every expected column.
//! 2. Atman's `recovery_jaccard` and `archetype_correlation` are both
//!    high (≥ 0.6 and ≥ 0.9 respectively) on both planted archetypes.
//! 3. Determinism score is 1.0 under fixed seed.
//! 4. Requesting a non-existent external tool produces a single
//!    `tool_not_available = 1` row rather than aborting.

use std::collections::HashMap;
use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn parse_tsv(
    path: &Path,
) -> (Vec<String>, Vec<HashMap<String, String>>) {
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

/// Deterministic Gaussian draws via LCG + Box-Muller.
fn make_gaussian_source(seed: u64, n: usize) -> Vec<f64> {
    let mut state = std::num::Wrapping(seed);
    let mut next_u = || {
        state = state * std::num::Wrapping(6364136223846793005_u64)
            + std::num::Wrapping(1442695040888963407_u64);
        ((state.0 >> 33) as f64) / (u32::MAX as f64)
    };
    let mut out = Vec::with_capacity(n);
    while out.len() < n {
        let u1 = next_u().max(1e-12);
        let u2 = next_u();
        // Box-Muller; take only one of the two to keep the LCG stream simple.
        let z = (-2.0_f64 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        out.push(z);
    }
    out
}

fn write_fixture(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    // 2 archetypes on 30 proteins. A1 loads on P01..P10; A2 on P11..P20.
    // Background proteins P21..P30 unloaded.
    let p = 30usize;
    let n = 40usize;
    let archs = [
        // A1
        {
            let mut v = vec![0.0_f64; p];
            for (i, slot) in v.iter_mut().enumerate().take(10) {
                *slot = 1.0 - (i as f64) * 0.03;
            }
            v
        },
        // A2
        {
            let mut v = vec![0.0_f64; p];
            for (offset, slot) in v.iter_mut().skip(10).take(10).enumerate() {
                *slot = 1.0 - (offset as f64) * 0.03;
            }
            v
        },
    ];

    // Write planted_loadings.tsv.
    let mut planted = String::from("archetype_id\tprotein\tloading\n");
    for (ai, row) in archs.iter().enumerate() {
        for (pi, &v) in row.iter().enumerate() {
            planted.push_str(&format!("P_A{:02}\tG{:03}\t{v:.6}\n", ai + 1, pi + 1));
        }
    }
    std::fs::write(dir.join("planted_loadings.tsv"), planted).unwrap();

    // Sources: per-subject heavy-tailed activations for each archetype.
    let s1 = make_gaussian_source(77, n);
    let s2 = make_gaussian_source(88, n);
    let noise = make_gaussian_source(99, n * p);

    let mut abundance = String::from("sample_id");
    for pi in 0..p {
        abundance.push('\t');
        abundance.push_str(&format!("G{:03}", pi + 1));
    }
    abundance.push('\n');
    for si in 0..n {
        abundance.push_str(&format!("S{si:03}"));
        for pi in 0..p {
            // Heavy-tailed-ish: cube the source to break Gaussianity.
            let v = s1[si].powi(3) * archs[0][pi]
                + s2[si].powi(3) * archs[1][pi]
                + 0.1 * noise[si * p + pi];
            abundance.push('\t');
            abundance.push_str(&format!("{v:.6}"));
        }
        abundance.push('\n');
    }
    std::fs::write(dir.join("abundance.tsv"), abundance).unwrap();
}

#[test]
fn bench_decompose_atman_recovers_planted_archetypes() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    let output = tmp.path().join("bench.tsv");
    write_fixture(&fixture);

    let status = run_atman(&[
        "bench", "decompose",
        "--fixture", fixture.to_str().unwrap(),
        "--tools", "atman",
        "--seed", "20260420",
        "--top-n", "10",
        "--output", output.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "bench decompose failed:\nstderr:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );
    let (header, rows) = parse_tsv(&output);
    for col in [
        "tool",
        "planted_archetype",
        "archetype_correlation",
        "recovery_jaccard",
        "runtime_seconds",
        "determinism_score",
        "tool_not_available",
    ] {
        assert!(header.iter().any(|h| h == col), "missing column {col}");
    }
    let atman_rows: Vec<&HashMap<String, String>> = rows
        .iter()
        .filter(|r| r.get("tool").map(|t| t == "atman").unwrap_or(false))
        .collect();
    assert_eq!(atman_rows.len(), 2, "expected one row per planted archetype");
    for r in &atman_rows {
        let corr: f64 = r["archetype_correlation"].parse().unwrap();
        let jacc: f64 = r["recovery_jaccard"].parse().unwrap();
        let det: f64 = r["determinism_score"].parse().unwrap_or(f64::NAN);
        assert!(
            corr >= 0.9,
            "atman should recover archetype {} to |pearson| ≥ 0.9; got {corr}",
            r["planted_archetype"]
        );
        assert!(
            jacc >= 0.6,
            "atman top-10 jaccard should be ≥ 0.6; got {jacc}",

        );
        assert!(
            (det - 1.0).abs() < 1e-12,
            "atman should be deterministic under fixed seed; got {det}"
        );
    }
}

#[test]
fn bench_decompose_missing_adapter_emits_tool_not_available() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    let output = tmp.path().join("bench.tsv");
    write_fixture(&fixture);

    let status = run_atman(&[
        "bench", "decompose",
        "--fixture", fixture.to_str().unwrap(),
        "--tools", "atman,nonexistent-tool",
        "--seed", "20260420",
        "--top-n", "10",
        "--output", output.to_str().unwrap(),
    ]);
    assert!(status.status.success());
    let (_, rows) = parse_tsv(&output);
    let nonexistent: Vec<&HashMap<String, String>> = rows
        .iter()
        .filter(|r| {
            r.get("tool")
                .map(|t| t == "nonexistent-tool")
                .unwrap_or(false)
        })
        .collect();
    assert_eq!(nonexistent.len(), 1, "expected one row for the missing tool");
    let flag: u8 = nonexistent[0]["tool_not_available"].parse().unwrap();
    assert_eq!(flag, 1);
    assert!(!nonexistent[0]["unavailable_reason"].is_empty());
}
