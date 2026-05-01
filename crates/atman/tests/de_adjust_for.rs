//! Parity tests: `atman de --adjust-for <covariates-tsv>` must reproduce R
//! reference implementations within 1e-6 absolute tolerance on log_fc and
//! the test statistic (moderated t for limma; t from msqrob2 robust LMM).
//!
//! J1 / K1 fixtures (committed under tests/fixtures/):
//!   - `de_adjust_for_input/`               — 12 samples × 30 proteins
//!   - `de_adjust_for_input/peptides.tsv`        — 3 peptides per protein
//!   - `de_adjust_for_input/peptide_measurements.tsv` — peptide abundances
//!   - `de_adjust_for_covariates.tsv`        — 2 NMF covariates
//!   - `de_adjust_for_limma_reference.tsv`   — R limma ground truth
//!   - `de_adjust_for_msqrob2_reference.tsv` — R msqrob2 ground truth
//!
//! Sign convention note (limma path):
//!   R's "conditionb" coefficient is mean_b − mean_a (condition b is the
//!   second level). atman with --groups b-a: A="b", B="a".
//!   build_two_group_design uses group_b_indicator=1 for samples in the B
//!   group (condition "a"). β_group_b = mean_a − mean_b.
//!   contrast [0, -1, 0, 0] negates that, giving mean_b − mean_a. Same sign
//!   as R. No sign correction needed when comparing atman b-a to R conditionb.
//!
//! Sign convention note (msqrob path):
//!   atman reports effect = mean_a − mean_b = −β_group_b. The R msqrob2
//!   reference is fit with condition levels c("a","b") so conditionb =
//!   mean_b − mean_a = +β_group_b. To compare: atman_effect ≈ −ref_log_fc.
//!   atman_t ≈ −ref_t_stat (same magnitude, opposite sign).
//!   The test compares absolute values, so no sign correction in the assertion.

use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

fn atman_bin() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_atman"))
}

/// Read a TSV that may start with a `# ...` comment line.
/// Returns rows-as-hashmaps keyed by column name.
fn read_tsv_skip_comment(path: &Path) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {:?}: {e}", path));
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .peekable();
    let header: Vec<&str> = match lines.next() {
        Some(h) => h.split('\t').collect(),
        None => return vec![],
    };
    let mut rows = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let mut m = HashMap::new();
        for (h, c) in header.iter().zip(cells.iter()) {
            m.insert(h.to_string(), c.to_string());
        }
        rows.push(m);
    }
    rows
}

/// Shared parity test runner.
///
/// Runs `atman de` with `extra_args`, then for each gene in the atman output
/// finds the matching row in `reference_tsv` and checks that:
///   - |atman_effect − ref_log_fc| ≤ `tol`
///   - |atman_t − ref_t_stat| ≤ `tol`
///
/// `sign_flip`: set `true` when atman's sign convention for the effect is
/// opposite to the R reference (e.g. msqrob path: atman reports mean_a−mean_b,
/// R reports mean_b−mean_a). When `true`, `ref_log_fc` and `ref_t_stat` are
/// negated before the absolute-delta check. The tolerance check uses absolute
/// delta, so only the negation matters for comparing at the assertion boundary.
///
/// Returns `(max_lfc_delta, max_t_delta)` for reporting.
fn run_adjusted_de_parity_test(
    extra_args: &[&str],
    output_dir: &Path,
    reference_tsv: &Path,
    sign_flip: bool,
    tol: f64,
) -> (f64, f64) {
    let fixture_dir = Path::new("tests/fixtures");
    let input_dir   = fixture_dir.join("de_adjust_for_input");
    let cov_tsv     = fixture_dir.join("de_adjust_for_covariates.tsv");

    std::fs::create_dir_all(output_dir).unwrap();

    let mut args = vec![
        "de",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--adjust-for",
        cov_tsv.to_str().unwrap(),
        "--groups",
        "b-a",
        "--min-pairs",
        "3",
    ];
    args.extend_from_slice(extra_args);

    let out = Command::new(atman_bin())
        .args(&args)
        .output()
        .expect("run atman");

    assert!(
        out.status.success(),
        "atman de failed.\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr),
    );

    let atman_rows = read_tsv_skip_comment(&output_dir.join("de_results.tsv"));
    let ref_rows   = read_tsv_skip_comment(reference_tsv);

    assert_eq!(
        atman_rows.len(),
        ref_rows.len(),
        "row count mismatch: atman={} ref={}",
        atman_rows.len(),
        ref_rows.len()
    );

    // Index reference by gene_symbol.
    let ref_by_gene: HashMap<String, HashMap<String, String>> = ref_rows
        .into_iter()
        .map(|r| (r["gene_symbol"].clone(), r))
        .collect();

    let mut max_lfc_delta: f64 = 0.0;
    let mut max_t_delta:   f64 = 0.0;

    for row in &atman_rows {
        let gene = row.get("gene_symbol").expect("gene_symbol column");
        let r = ref_by_gene
            .get(gene)
            .unwrap_or_else(|| panic!("no R reference for gene {gene}"));

        let atman_lfc: f64 = row
            .get("mean_diff")
            .or_else(|| row.get("effect_size"))
            .unwrap_or_else(|| panic!("no effect column for gene {gene}"))
            .parse()
            .unwrap_or_else(|_| panic!("cannot parse effect for gene {gene}"));

        let atman_t: f64 = row
            .get("t")
            .unwrap_or_else(|| panic!("no t column for gene {gene}"))
            .parse()
            .unwrap_or_else(|_| panic!("cannot parse t for gene {gene}"));

        let ref_lfc_raw: f64 = r["log_fc"].parse().expect("ref log_fc");
        let ref_t_raw:   f64 = r["t_stat"].parse().expect("ref t_stat");

        // When sign conventions differ, negate the R values before delta.
        let ref_lfc = if sign_flip { -ref_lfc_raw } else { ref_lfc_raw };
        let ref_t   = if sign_flip { -ref_t_raw   } else { ref_t_raw   };

        let lfc_delta = (atman_lfc - ref_lfc).abs();
        let t_delta   = (atman_t   - ref_t).abs();

        if lfc_delta > max_lfc_delta { max_lfc_delta = lfc_delta; }
        if t_delta   > max_t_delta   { max_t_delta   = t_delta;   }

        assert!(
            lfc_delta <= tol,
            "gene={gene}: log_fc delta {lfc_delta:.2e} > {tol:.0e}\n  atman={atman_lfc}  R={ref_lfc} (sign_flip={sign_flip})"
        );
        assert!(
            t_delta <= tol,
            "gene={gene}: t_stat delta {t_delta:.2e} > {tol:.0e}\n  atman={atman_t}  R={ref_t} (sign_flip={sign_flip})"
        );
    }

    (max_lfc_delta, max_t_delta)
}

#[test]
fn limma_adjusted_de_matches_r_reference() {
    let tmp = tempfile::tempdir().expect("tempdir");
    let output_dir = tmp.path().join("de_out");

    let fixture_dir = Path::new("tests/fixtures");
    let ref_tsv = fixture_dir.join("de_adjust_for_limma_reference.tsv");

    let (max_lfc_delta, max_t_delta) = run_adjusted_de_parity_test(
        &[
            "--test", "limma",
            "--trend", "false",
            "--robust", "false",
        ],
        &output_dir,
        &ref_tsv,
        // limma path: atman b-a means effect = mean_b − mean_a,
        // same sign as R's conditionb. No flip needed.
        false,
        1e-6,
    );

    eprintln!(
        "limma_adjusted_de_matches_r_reference: max |Δlog_fc|={:.3e}  max |Δt|={:.3e}",
        max_lfc_delta, max_t_delta
    );

    // Verify sidecar records the --adjust-for paths.
    let sidecar_path = output_dir.join("de_results.tsv.run.json");
    assert!(sidecar_path.exists(), "sidecar missing");
    let sidecar: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar_path).unwrap()).unwrap();
    let adj = sidecar["args"]["adjust-for"].as_array();
    assert!(
        adj.is_some() && !adj.unwrap().is_empty(),
        "sidecar does not record adjust-for paths: {}",
        sidecar["args"]["adjust-for"]
    );
}

#[test]
fn msqrob2_adjusted_de_matches_r_reference() {
    let tmp = tempfile::tempdir().expect("tempdir");
    let output_dir = tmp.path().join("de_out");

    let fixture_dir = Path::new("tests/fixtures");
    let input_dir   = fixture_dir.join("de_adjust_for_input");
    let ref_tsv     = fixture_dir.join("de_adjust_for_msqrob2_reference.tsv");

    let (max_lfc_delta, max_t_delta) = run_adjusted_de_parity_test(
        &[
            "--test", "msqrob",
            "--peptide-metadata",
            input_dir.join("peptides.tsv").to_str().unwrap(),
            "--peptide-measurements",
            input_dir.join("peptide_measurements.tsv").to_str().unwrap(),
            "--min-peptides", "2",
            "--ridge-lambda", "0.0",
            "--robust", "false",
        ],
        &output_dir,
        &ref_tsv,
        // msqrob path: atman reports mean_a − mean_b; R msqrob2 conditionb
        // is mean_b − mean_a. Negate R values before delta check.
        true,
        1e-6,
    );

    eprintln!(
        "msqrob2_adjusted_de_matches_r_reference: max |Δlog_fc|={:.3e}  max |Δt|={:.3e}",
        max_lfc_delta, max_t_delta
    );

    // Verify sidecar records the --adjust-for paths.
    let sidecar_path = output_dir.join("de_results.tsv.run.json");
    assert!(sidecar_path.exists(), "sidecar missing");
    let sidecar: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar_path).unwrap()).unwrap();
    let adj = sidecar["args"]["adjust-for"].as_array();
    assert!(
        adj.is_some() && !adj.unwrap().is_empty(),
        "sidecar does not record adjust-for paths: {}",
        sidecar["args"]["adjust-for"]
    );
}
