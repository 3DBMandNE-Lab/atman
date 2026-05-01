//! Parity test: `atman de --test limma --adjust-for <covariates-tsv>` must
//! reproduce R limma + covariates within 1e-6 absolute tolerance on log_fc
//! and moderated t-statistic.
//!
//! J1 fixture (committed):
//!   - `tests/fixtures/de_adjust_for_input/`  — 12 samples × 30 proteins
//!   - `tests/fixtures/de_adjust_for_covariates.tsv`  — 2 NMF covariates
//!   - `tests/fixtures/de_adjust_for_limma_reference.tsv`  — R ground truth
//!
//! Sign convention note:
//!   R's "conditionb" coefficient is mean_b − mean_a.
//!   atman's A-B comparison (contrast [0, -1, 0, 0]) is mean_a − mean_b.
//!   Both the reference t_stat and log_fc must be sign-negated before
//!   comparison with atman's output.

use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

fn atman_bin() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_atman"))
}

/// Read a TSV that may start with a `# ...` comment line.
/// Returns (headers, rows-as-hashmaps).
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

#[test]
fn limma_adjusted_de_matches_r_reference() {
    let tmp = tempfile::tempdir().expect("tempdir");
    let output_dir = tmp.path().join("de_out");
    std::fs::create_dir_all(&output_dir).unwrap();

    let fixture_dir = Path::new("tests/fixtures");
    let input_dir   = fixture_dir.join("de_adjust_for_input");
    let cov_tsv     = fixture_dir.join("de_adjust_for_covariates.tsv");
    let ref_tsv     = fixture_dir.join("de_adjust_for_limma_reference.tsv");

    // Run atman de with --adjust-for.
    let out = Command::new(atman_bin())
        .args([
            "de",
            "--input-dir",
            input_dir.to_str().unwrap(),
            "--output-dir",
            output_dir.to_str().unwrap(),
            "--test",
            "limma",
            "--adjust-for",
            cov_tsv.to_str().unwrap(),
            "--groups",
            "b-a",
            "--min-pairs",
            "3",
            "--trend",
            "false",
            "--robust",
            "false",
        ])
        .output()
        .expect("run atman");

    assert!(
        out.status.success(),
        "atman de failed.\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr),
    );

    let atman_rows = read_tsv_skip_comment(&output_dir.join("de_results.tsv"));
    let ref_rows   = read_tsv_skip_comment(&ref_tsv);

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
    const TOL: f64 = 1e-6;

    for row in &atman_rows {
        let gene = row.get("gene_symbol").expect("gene_symbol column");
        let r = ref_by_gene
            .get(gene)
            .unwrap_or_else(|| panic!("no R reference for gene {gene}"));

        // atman: mean_a − mean_b  (contrast [0, -1, 0, 0])
        // R ref: mean_b − mean_a  (conditionb coefficient)
        // → negate R values before comparing.
        let atman_lfc: f64 = row
            .get("mean_diff")
            .or_else(|| row.get("effect_size"))
            .unwrap_or_else(|| panic!("no effect column for gene {gene}"))
            .parse()
            .unwrap_or_else(|_| panic!("cannot parse effect for gene {gene}"));
        let ref_lfc: f64 = r["log_fc"].parse().expect("ref log_fc");
        let ref_lfc_negated = -ref_lfc;

        let atman_t: f64 = row
            .get("t")
            .unwrap_or_else(|| panic!("no t column for gene {gene}"))
            .parse()
            .unwrap_or_else(|_| panic!("cannot parse t for gene {gene}"));
        let ref_t: f64 = r["t_stat"].parse().expect("ref t_stat");
        let ref_t_negated = -ref_t;

        let lfc_delta = (atman_lfc - ref_lfc_negated).abs();
        let t_delta   = (atman_t   - ref_t_negated).abs();

        if lfc_delta > max_lfc_delta { max_lfc_delta = lfc_delta; }
        if t_delta   > max_t_delta   { max_t_delta   = t_delta;   }

        assert!(
            lfc_delta <= TOL,
            "gene={gene}: log_fc delta {lfc_delta:.2e} > {TOL:.0e}\n  atman={atman_lfc}  R={ref_lfc_negated}"
        );
        assert!(
            t_delta <= TOL,
            "gene={gene}: t_stat delta {t_delta:.2e} > {TOL:.0e}\n  atman={atman_t}  R={ref_t_negated}"
        );
    }

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
