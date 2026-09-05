use std::collections::HashMap;
use std::path::PathBuf;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/de_csf")
}

fn read_rows(path: &std::path::Path) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| {
            header
                .iter()
                .zip(l.split('\t'))
                .map(|(h, v)| (h.to_string(), v.to_string()))
                .collect()
        })
        .collect()
}

fn num(r: &HashMap<String, String>, c: &str) -> f64 {
    r[c].parse().unwrap_or_else(|_| panic!("{c} = {:?}", r[c]))
}

#[test]
fn de_ols_without_include_controls_has_no_control_samples() {
    let tmp = tempfile::tempdir().unwrap();
    let r = run_atman(&[
        "de",
        "--input-dir",
        fixture().to_str().unwrap(),
        "--output-dir",
        tmp.path().to_str().unwrap(),
        "--test",
        "ols",
        "--design",
        "~ condition + age + sex",
        "--groups",
        "Case-Ctrl",
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("singular"));
}

#[test]
fn de_ols_include_controls_matches_reference_and_emits_cohen_d() {
    let tmp = tempfile::tempdir().unwrap();
    let r = run_atman(&[
        "de",
        "--input-dir",
        fixture().to_str().unwrap(),
        "--output-dir",
        tmp.path().to_str().unwrap(),
        "--test",
        "ols",
        "--design",
        "~ condition + age + sex",
        "--groups",
        "Case-Ctrl",
        "--include-controls",
        "--min-pairs",
        "5",
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&tmp.path().join("de_results.tsv"));
    let reference = read_rows(&fixture().join("reference.tsv"));
    assert_eq!(rows.len(), 3);
    for rr in &reference {
        let g = rows
            .iter()
            .find(|x| x["gene_symbol"] == rr["gene_symbol"])
            .unwrap();
        assert!(
            (num(g, "mean_diff") - num(rr, "beta")).abs() < 1e-8,
            "beta {}",
            rr["gene_symbol"]
        );
        assert!(
            (num(g, "p_value") - num(rr, "p_value")).abs() < 1e-8,
            "p {}",
            rr["gene_symbol"]
        );
        if rr["cohen_d"].is_empty() {
            // A single observed case: d is undefined on both sides.
            assert!(g["effect_size"].is_empty());
            assert_eq!(g["effect_size_method"], "");
        } else {
            assert!(
                (num(g, "effect_size") - num(rr, "cohen_d")).abs() < 1e-8,
                "d {}",
                rr["gene_symbol"]
            );
            assert_eq!(g["effect_size_method"], "cohen_d");
        }
        assert_eq!(g["n_a"], rr["n_a"]);
        assert_eq!(g["n_b"], rr["n_b"]);
        assert_eq!(g["n_pairs"], rr["n_obs"]);
        if !rr["cohen_d"].is_empty() {
            assert!(!g["ci_low"].is_empty() && !g["ci_high"].is_empty());
        }
    }
    let sidecar = std::fs::read_to_string(tmp.path().join("de_results.tsv.run.json")).unwrap();
    assert!(sidecar.contains("\"include-controls\": true"));
}

#[test]
fn de_max_missing_fraction_drops_sparse_protein_and_records_counts() {
    let tmp = tempfile::tempdir().unwrap();
    let r = run_atman(&[
        "de",
        "--input-dir",
        fixture().to_str().unwrap(),
        "--output-dir",
        tmp.path().to_str().unwrap(),
        "--test",
        "ols",
        "--design",
        "~ condition",
        "--groups",
        "Case-Ctrl",
        "--include-controls",
        "--max-missing-fraction",
        "0.3",
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&tmp.path().join("de_results.tsv"));
    assert_eq!(rows.len(), 2);
    assert!(rows.iter().all(|r| r["gene_symbol"] != "G3"));
    let sidecar = std::fs::read_to_string(tmp.path().join("de_results.tsv.run.json")).unwrap();
    assert!(sidecar.contains("\"missingness_filter\""));
    assert!(sidecar.contains("\"n_dropped\": 1"));
}

#[test]
fn de_condition_col_subset_and_collapse_others() {
    let tmp = tempfile::tempdir().unwrap();
    let r = run_atman(&[
        "de",
        "--input-dir",
        fixture().to_str().unwrap(),
        "--output-dir",
        tmp.path().to_str().unwrap(),
        "--test",
        "welch-t",
        "--groups",
        "A-Other",
        "--include-controls",
        "--condition-col",
        "grp",
        "--subset",
        "site!=Y",
        "--collapse-others",
        "Other",
        "--min-pairs",
        "2",
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&tmp.path().join("de_results.tsv"));
    let g1 = rows.iter().find(|r| r["gene_symbol"] == "G1").unwrap();
    assert_eq!(g1["comparison"], "A-Other");
    assert_eq!(g1["n_a"], "2");
    assert_eq!(g1["n_b"], "4");
    let sidecar = std::fs::read_to_string(tmp.path().join("de_results.tsv.run.json")).unwrap();
    assert!(sidecar.contains("\"condition-col\": \"grp\""));
    assert!(sidecar.contains("site!=Y"));
}

#[test]
fn de_continuous_contrast_with_expression_term() {
    let tmp = tempfile::tempdir().unwrap();
    let r = run_atman(&[
        "de",
        "--input-dir",
        fixture().to_str().unwrap(),
        "--output-dir",
        tmp.path().to_str().unwrap(),
        "--test",
        "ols",
        "--design",
        "~ log10(age) + sex",
        "--contrast",
        "log10(age)",
        "--include-controls",
        "--min-pairs",
        "5",
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&tmp.path().join("de_results.tsv"));
    let reference = read_rows(&fixture().join("reference_continuous.tsv"));
    assert_eq!(rows.len(), 3);
    for rr in &reference {
        let g = rows
            .iter()
            .find(|x| x["gene_symbol"] == rr["gene_symbol"])
            .unwrap();
        assert_eq!(g["comparison"], "log10(age)");
        assert!(g["mean_a"].is_empty() && g["mean_b"].is_empty());
        assert!(g["effect_size"].is_empty());
        assert!((num(g, "mean_diff") - num(rr, "beta")).abs() < 1e-8);
        assert!((num(g, "t") - num(rr, "t")).abs() < 1e-8);
        assert!((num(g, "p_value") - num(rr, "p_value")).abs() < 1e-8);
        assert_eq!(g["n_pairs"], rr["n_obs"]);
    }
    let covs = read_rows(&tmp.path().join("de_covariates.tsv"));
    assert!(covs.iter().any(|c| c["covariate"] == "sexM"));
}

#[test]
fn de_z_age_leaves_condition_effect_unchanged() {
    let run = |design: &str, dir: &std::path::Path| {
        let r = run_atman(&[
            "de",
            "--input-dir",
            fixture().to_str().unwrap(),
            "--output-dir",
            dir.to_str().unwrap(),
            "--test",
            "ols",
            "--design",
            design,
            "--groups",
            "Case-Ctrl",
            "--include-controls",
            "--min-pairs",
            "5",
        ]);
        assert!(
            r.status.success(),
            "stderr:\n{}",
            String::from_utf8_lossy(&r.stderr)
        );
        read_rows(&dir.join("de_results.tsv"))
    };
    let t1 = tempfile::tempdir().unwrap();
    let t2 = tempfile::tempdir().unwrap();
    let raw = run("~ condition + age + sex", t1.path());
    let z = run("~ condition + z(age) + sex", t2.path());
    for (a, b) in raw.iter().zip(z.iter()) {
        assert_eq!(a["gene_symbol"], b["gene_symbol"]);
        assert!((num(a, "mean_diff") - num(b, "mean_diff")).abs() < 1e-9);
        assert!((num(a, "p_value") - num(b, "p_value")).abs() < 1e-9);
    }
    let covs = read_rows(&t2.path().join("de_covariates.tsv"));
    assert!(covs.iter().any(|c| c["covariate"] == "z(age)"));
}

#[test]
fn de_continuous_contrast_rejects_groups_and_condition_term() {
    let tmp = tempfile::tempdir().unwrap();
    let r = run_atman(&[
        "de",
        "--input-dir",
        fixture().to_str().unwrap(),
        "--output-dir",
        tmp.path().to_str().unwrap(),
        "--test",
        "ols",
        "--design",
        "~ condition + log10(age)",
        "--contrast",
        "log10(age)",
        "--include-controls",
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("--groups"));
}
