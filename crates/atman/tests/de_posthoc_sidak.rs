//! DEBT-1 parity test: `atman de --post-hoc sidak` must match the
//! R reference `fixtures/posthoc_sidak_reference.R` (lm + pairwise
//! linear contrasts + Sidak adjustment) to ≥ 6 decimal places on a
//! 36-sample 3-level fixture.

use statrs::distribution::{ContinuousCDF, StudentsT};
use std::collections::HashMap;
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

fn copy_fixture_to_canonical(tmp: &Path) -> std::path::PathBuf {
    let input = tmp.join("canonical");
    std::fs::create_dir_all(&input).unwrap();
    let fixtures = Path::new("tests/fixtures");
    std::fs::copy(
        fixtures.join("posthoc_sidak_samples.tsv"),
        input.join("samples.tsv"),
    )
    .unwrap();
    std::fs::copy(
        fixtures.join("posthoc_sidak_proteins.tsv"),
        input.join("proteins.tsv"),
    )
    .unwrap();
    std::fs::copy(
        fixtures.join("posthoc_sidak_qc_measurements.tsv"),
        input.join("measurements.tsv"),
    )
    .unwrap();
    std::fs::copy(
        fixtures.join("posthoc_sidak_measurements.tsv"),
        input.join("measurements.tsv"),
    )
    .unwrap();
    input
}

#[test]
fn posthoc_sidak_matches_r_lm_contrasts_to_6_decimals() {
    let ref_path = Path::new("tests/fixtures/posthoc_sidak_reference.tsv");
    assert!(
        ref_path.exists(),
        "missing R reference; run `Rscript tests/fixtures/posthoc_sidak_reference.R` first"
    );
    let tmp = tempfile::tempdir().unwrap();
    let input = copy_fixture_to_canonical(tmp.path());
    let output = tmp.path().join("posthoc_out");

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "ols",
        "--design",
        "~ stage + age",
        "--post-hoc",
        "sidak",
        "--contrast-list",
        "MCI-CN,AD-CN,AD-MCI",
        "--post-hoc-factor",
        "stage",
        "--min-pairs",
        "4",
    ]);
    assert!(
        out.status.success(),
        "atman failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );

    let (_, atman_rows) = parse_tsv(&output.join("de_results.tsv"));
    let (_, ref_rows) = parse_tsv(ref_path);
    // Key: (gene_symbol, comparison)
    type Key = (String, String);
    let atman_by_key: HashMap<Key, &HashMap<String, String>> = atman_rows
        .iter()
        .map(|r| ((r["gene_symbol"].clone(), r["comparison"].clone()), r))
        .collect();
    let mut max_p_err = 0.0_f64;
    let mut max_adj_err = 0.0_f64;
    let mut max_est_err = 0.0_f64;
    let mut max_se_err = 0.0_f64;
    for ref_row in &ref_rows {
        let key = (
            ref_row["gene_symbol"].clone(),
            ref_row["comparison"].clone(),
        );
        let atman_row = atman_by_key
            .get(&key)
            .unwrap_or_else(|| panic!("atman missing row for {key:?}"));
        let ref_p: f64 = ref_row["posthoc_p"].parse().unwrap();
        let ref_adj: f64 = ref_row["posthoc_adj_p"].parse().unwrap();
        let ref_est: f64 = ref_row["estimate"].parse().unwrap();
        let ref_se: f64 = ref_row["se"].parse().unwrap();
        let atman_p: f64 = atman_row["posthoc_p"].parse().unwrap();
        let atman_adj: f64 = atman_row["posthoc_adj_p"].parse().unwrap();
        let atman_est: f64 = atman_row["mean_diff"].parse().unwrap();
        // SE isn't emitted on DeResultRow, but the CI half-width is
        // `t_crit * se` where `t_crit` is the two-sided 97.5% Student-t
        // quantile at the fit's df (TASK-023: matches the t-based p-value,
        // replacing the prior fixed 1.96 z-multiplier). Recover SE via
        // (mean_diff - ci_low) / t_crit using the df emitted on the row.
        let ci_low: f64 = atman_row["ci_low"].parse().unwrap();
        let df: f64 = atman_row["df"].parse().unwrap();
        let t_crit = StudentsT::new(0.0, 1.0, df).unwrap().inverse_cdf(0.975);
        let atman_se = (atman_est - ci_low) / t_crit;
        max_p_err = max_p_err.max((ref_p - atman_p).abs());
        max_adj_err = max_adj_err.max((ref_adj - atman_adj).abs());
        max_est_err = max_est_err.max((ref_est - atman_est).abs());
        max_se_err = max_se_err.max((ref_se - atman_se).abs());
    }
    eprintln!(
        "sidak parity: max |Δp|={max_p_err:.3e}, \
         max |Δp_adj|={max_adj_err:.3e}, \
         max |Δest|={max_est_err:.3e}, max |Δse|={max_se_err:.3e}"
    );
    assert!(max_p_err < 1e-6, "posthoc_p drift {max_p_err:e}");
    assert!(max_adj_err < 1e-6, "posthoc_adj_p drift {max_adj_err:e}");
    assert!(max_est_err < 1e-6, "estimate drift {max_est_err:e}");
    assert!(max_se_err < 1e-4, "se drift {max_se_err:e}");
}

#[test]
fn posthoc_sidak_refuses_unknown_factor_level() {
    let tmp = tempfile::tempdir().unwrap();
    let input = copy_fixture_to_canonical(tmp.path());
    let output = tmp.path().join("bad_level");
    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "ols",
        "--design",
        "~ stage + age",
        "--post-hoc",
        "sidak",
        "--contrast-list",
        "MCI-BOGUS",
        "--post-hoc-factor",
        "stage",
    ]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("BOGUS"), "unexpected stderr: {stderr}");
}
