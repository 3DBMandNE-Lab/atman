//! DEBT-6 integration test: `atman de --post-hoc dunnett` on the
//! shared 3-level stage fixture used by sidak and tukey.
//!
//! Verifies the dispatch end-to-end:
//!
//! 1. Self-consistency — recompute the expected `1 − pdunnett(|t|, m,
//!    df, ρ=0.5)` via `atman_core::multivariate_t` and check the
//!    CLI-emitted `posthoc_adj_p` matches to floating-point precision.
//! 2. Magnitude on planted effects — `RESPONDER` has a strong stage
//!    gradient, so its two non-reference contrasts (`CN vs AD`,
//!    `MCI vs AD` with reference level AD) must hit `adj_p < 0.05`;
//!    `NULLP`'s must stay `adj_p > 0.10`.
//! 3. Family size — exactly `levels − 1 = 2` rows per protein.

use std::collections::HashMap;
use std::path::Path;
use std::process::{Command, Output};

use atman_core::multivariate_t::pdunnett;

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
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

fn copy_fixture_to_canonical(tmp: &Path) -> std::path::PathBuf {
    let input = tmp.join("canonical");
    std::fs::create_dir_all(&input).unwrap();
    let fixtures = Path::new("tests/fixtures");
    for (src, dst) in [
        ("posthoc_sidak_samples.tsv", "samples.tsv"),
        ("posthoc_sidak_proteins.tsv", "proteins.tsv"),
        ("posthoc_sidak_qc_measurements.tsv", "qc_measurements.tsv"),
        ("posthoc_sidak_measurements.tsv", "measurements.tsv"),
    ] {
        std::fs::copy(fixtures.join(src), input.join(dst)).unwrap();
    }
    input
}

#[test]
fn posthoc_dunnett_adjustment_is_self_consistent_with_core_pdunnett() {
    let tmp = tempfile::tempdir().unwrap();
    let input = copy_fixture_to_canonical(tmp.path());
    let output = tmp.path().join("dunnett_out");
    let out = run_atman(&[
        "de",
        "--input-dir", input.to_str().unwrap(),
        "--output-dir", output.to_str().unwrap(),
        "--test", "ols",
        "--design", "~ stage + age",
        "--post-hoc", "dunnett",
        "--post-hoc-factor", "stage",
        "--min-pairs", "4",
    ]);
    assert!(
        out.status.success(),
        "atman failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    let (_, rows) = parse_tsv(&output.join("de_results.tsv"));
    assert!(!rows.is_empty());

    // The stage factor has 3 levels ⇒ 2 non-reference contrasts per protein.
    let mut per_protein: HashMap<String, usize> = HashMap::new();
    for r in &rows {
        *per_protein.entry(r["gene_symbol"].clone()).or_default() += 1;
    }
    for (gene, n) in &per_protein {
        assert_eq!(*n, 2, "{gene} should have 2 Dunnett contrasts, got {n}");
    }

    let m_contrasts = 2usize;
    let mut max_drift = 0.0_f64;
    for r in &rows {
        assert_eq!(r["posthoc_method"], "dunnett");
        let t_stat: f64 = r["t"].parse().unwrap_or(f64::NAN);
        let df: f64 = r["df"].parse().unwrap_or(f64::NAN);
        let atman_adj: f64 = r["posthoc_adj_p"].parse().unwrap_or(f64::NAN);
        if !t_stat.is_finite() || !df.is_finite() || !atman_adj.is_finite() {
            continue;
        }
        let expected = (1.0 - pdunnett(t_stat.abs(), m_contrasts, df, 0.5)).clamp(0.0, 1.0);
        max_drift = max_drift.max((atman_adj - expected).abs());
    }
    assert!(
        max_drift < 1e-10,
        "atman's dunnett adj_p must match atman_core::pdunnett to floating-point; \
         max drift = {max_drift:e}"
    );
}

#[test]
fn posthoc_dunnett_magnitude_detects_planted_stage_effect() {
    let tmp = tempfile::tempdir().unwrap();
    let input = copy_fixture_to_canonical(tmp.path());
    let output = tmp.path().join("dunnett_magnitude");
    let out = run_atman(&[
        "de",
        "--input-dir", input.to_str().unwrap(),
        "--output-dir", output.to_str().unwrap(),
        "--test", "ols",
        "--design", "~ stage + age",
        "--post-hoc", "dunnett",
        "--post-hoc-factor", "stage",
        "--min-pairs", "4",
    ]);
    assert!(out.status.success());
    let (_, rows) = parse_tsv(&output.join("de_results.tsv"));

    let responder: Vec<&HashMap<String, String>> = rows
        .iter()
        .filter(|r| r.get("gene_symbol").map(|g| g == "RESPONDER").unwrap_or(false))
        .collect();
    let nullp: Vec<&HashMap<String, String>> = rows
        .iter()
        .filter(|r| r.get("gene_symbol").map(|g| g == "NULLP").unwrap_or(false))
        .collect();
    assert_eq!(responder.len(), 2);
    assert_eq!(nullp.len(), 2);
    for r in &responder {
        let adj: f64 = r["posthoc_adj_p"].parse().unwrap();
        assert!(
            adj < 0.05,
            "planted effect should reject Dunnett: {:?} adj={adj}",
            r["comparison"]
        );
    }
    for r in &nullp {
        let adj: f64 = r["posthoc_adj_p"].parse().unwrap();
        assert!(
            adj > 0.10,
            "planted null should not reject Dunnett: {:?} adj={adj}",
            r["comparison"]
        );
    }
}
