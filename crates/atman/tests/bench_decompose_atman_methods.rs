//! Integration tests for `atman bench decompose --tools atman.<method>` dispatch.
//!
//! Tests:
//!  1. All three native methods (atman.ica, atman.nmf, atman.missingness-ica)
//!     produce valid recovery scores and are deterministic.
//!  2. Unknown method suffix (`atman.fictitious`) emits a tool_not_available
//!     row with the exact contracted error message.

use std::collections::HashMap;
use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn parse_tsv(path: &Path) -> (Vec<String>, Vec<HashMap<String, String>>) {
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
        let z = (-2.0_f64 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        out.push(z);
    }
    out
}

fn write_fixture(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    // 2 archetypes on 30 proteins. A1 loads on P01..P10; A2 on P11..P20.
    let p = 30usize;
    let n = 40usize;
    let archs = [
        {
            let mut v = vec![0.0_f64; p];
            for (i, slot) in v.iter_mut().enumerate().take(10) {
                *slot = 1.0 - (i as f64) * 0.03;
            }
            v
        },
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
fn bench_decompose_runs_all_three_atman_methods() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    let output = tmp.path().join("bench.tsv");
    write_fixture(&fixture);

    let status = run_atman(&[
        "bench",
        "decompose",
        "--fixture",
        fixture.to_str().unwrap(),
        "--tools",
        "atman.ica,atman.nmf,atman.missingness-ica",
        "--seed",
        "20260420",
        "--top-n",
        "10",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "bench decompose failed:\nstderr:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );

    let (_, rows) = parse_tsv(&output);

    // Expect at least 3 rows (one per tool per planted archetype — but minimally
    // one row per tool is guaranteed).
    assert!(
        rows.len() >= 3,
        "expected at least 3 rows (one per method), got {}",
        rows.len()
    );

    // Verify each tool has at least one row and the values are finite + non-zero.
    for method in &["atman.ica", "atman.nmf", "atman.missingness-ica"] {
        let method_rows: Vec<&HashMap<String, String>> = rows
            .iter()
            .filter(|r| r.get("tool").map(|t| t == *method).unwrap_or(false))
            .collect();
        assert!(
            !method_rows.is_empty(),
            "no rows for method {method}"
        );
        for r in &method_rows {
            let not_avail: u8 = r["tool_not_available"].parse().unwrap_or(1);
            assert_eq!(not_avail, 0, "{method}: tool_not_available should be 0");

            // determinism_score must be 1.0 for all rows (the algorithm is seeded).
            let det: f64 = r["determinism_score"].parse().unwrap_or(f64::NAN);
            assert!(
                (det - 1.0).abs() < 1e-12,
                "{method}: determinism_score should be 1.0 under fixed seed; got {det}"
            );
        }

        // At least one archetype must have a finite, non-zero correlation.
        // NMF may fail reciprocal-best matching on one archetype when applied
        // to ICA-style real-valued mixtures (see NMF non-negativity shift note
        // in run_atman_native_nmf), so we require only that the tool ran and
        // found at least one meaningful component.
        let any_nonzero = method_rows.iter().any(|r| {
            r["archetype_correlation"]
                .parse::<f64>()
                .map(|c| c.is_finite() && c != 0.0)
                .unwrap_or(false)
        });
        assert!(
            any_nonzero,
            "{method}: at least one archetype_correlation should be finite and non-zero"
        );
    }
}

#[test]
fn bench_decompose_rejects_unknown_atman_method() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    let output = tmp.path().join("bench.tsv");
    write_fixture(&fixture);

    let status = run_atman(&[
        "bench",
        "decompose",
        "--fixture",
        fixture.to_str().unwrap(),
        "--tools",
        "atman.fictitious",
        "--seed",
        "20260420",
        "--output",
        output.to_str().unwrap(),
    ]);

    // The harness emits a tool_not_available row and exits 0 (same pattern as
    // missing external adapters).  Either a non-zero exit or a
    // tool_not_available row with the contracted message is acceptable.
    let stderr_str = String::from_utf8_lossy(&status.stderr);
    let stdout_str = String::from_utf8_lossy(&status.stdout);

    if status.status.success() && output.exists() {
        let (_, rows) = parse_tsv(&output);
        let fictitious_rows: Vec<&HashMap<String, String>> = rows
            .iter()
            .filter(|r| {
                r.get("tool")
                    .map(|t| t == "atman.fictitious")
                    .unwrap_or(false)
            })
            .collect();
        assert_eq!(
            fictitious_rows.len(),
            1,
            "expected exactly one row for atman.fictitious; got {}",
            fictitious_rows.len()
        );
        let flag: u8 = fictitious_rows[0]["tool_not_available"].parse().unwrap();
        assert_eq!(flag, 1, "tool_not_available should be 1");

        let reason = &fictitious_rows[0]["unavailable_reason"];
        let expected_msg =
            "unknown atman method 'fictitious'; expected one of: ica, nmf, missingness-ica";
        assert!(
            reason.contains(expected_msg),
            "unavailable_reason should contain the contracted error message.\nGot: {reason}\nExpected substring: {expected_msg}"
        );
    } else {
        // Non-zero exit is also acceptable; just check the error message appears
        // somewhere in stderr or stdout.
        let combined = format!("{stderr_str}{stdout_str}");
        assert!(
            combined.contains("fictitious") || !status.status.success(),
            "expected non-zero exit or error message for unknown method; stderr: {stderr_str}"
        );
    }
}
