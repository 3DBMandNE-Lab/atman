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

fn dup_fixture(dir: &std::path::Path) {
    std::fs::write(dir.join("samples.tsv"), "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\nS1\tS1\tCase\t0\tbio\t1\nS2\tS2\tCase\t0\tbio\t2\nS3\tS3\tCtrl\t1\tbio\t3\nS4\tS4\tCtrl\t1\tbio\t4\n").unwrap();
    std::fs::write(dir.join("proteins.tsv"), "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\nmaxquant_lfq\tA1\tA1\tG\tP\t\nmaxquant_lfq\tA2\tA2\tG\tP\t\nmaxquant_lfq\tB1\tB1\tH\tP\t\n").unwrap();
    let mut m = String::from("platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n");
    let mut order = 1;
    let mut push = |sid: &str, assay: &str, gene: &str, v: f64| {
        m.push_str(&format!("maxquant_lfq\t{sid}\t{assay}\t{gene}\tP\t{v}\t{v}\t{v}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"));
        order += 1;
    };
    // G / A1 on every sample; G / A2 on the two cases only; H / B1 everywhere.
    for (sid, v) in [("S1", 1.0), ("S2", 2.0), ("S3", 3.0), ("S4", 4.0)] {
        push(sid, "A1", "G", v);
        push(sid, "B1", "H", 10.0 + v);
    }
    push("S1", "A2", "G", 5.0);
    push("S2", "A2", "G", 6.0);
    std::fs::write(dir.join("measurements.tsv"), m).unwrap();
}

#[test]
fn de_collapse_genes_rules() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    std::fs::create_dir_all(&input).unwrap();
    dup_fixture(&input);
    let run = |rule: &str| {
        let out = tmp.path().join(format!("out_{rule}"));
        let r = run_atman(&[
            "de",
            "--input-dir",
            input.to_str().unwrap(),
            "--output-dir",
            out.to_str().unwrap(),
            "--test",
            "welch-t",
            "--groups",
            "Case-Ctrl",
            "--include-controls",
            "--min-pairs",
            "2",
            "--collapse-genes",
            rule,
        ]);
        assert!(
            r.status.success(),
            "stderr:\n{}",
            String::from_utf8_lossy(&r.stderr)
        );
        let rows = read_rows(&out.join("de_results.tsv"));
        let g = rows
            .iter()
            .find(|x| x["gene_symbol"] == "G")
            .unwrap()
            .clone();
        let sidecar = std::fs::read_to_string(out.join("de_results.tsv.run.json")).unwrap();
        (g, sidecar)
    };
    // none: lexically first assay A1 ⇒ Case [1,2] vs Ctrl [3,4] ⇒ mean_diff -2
    let (g, sidecar) = run("none");
    assert!((num(&g, "mean_diff") + 2.0).abs() < 1e-12);
    assert_eq!(g["assay_id"], "A1");
    assert!(sidecar.contains("\"n_genes_with_multiple_assays\": 1"));
    assert!(sidecar.contains("\"n_extra_assays\": 1"));
    assert!(sidecar.contains("\"collapse-genes\": \"none\""));
    // max-observed: A1 (4 obs) beats A2 (2 obs) ⇒ same as none
    let (g, _) = run("max-observed");
    assert!((num(&g, "mean_diff") + 2.0).abs() < 1e-12);
    // mean: Case S1 = (1+5)/2 = 3, S2 = (2+6)/2 = 4 ⇒ 3.5 - 3.5 = 0
    let (g, sidecar) = run("mean");
    assert!(num(&g, "mean_diff").abs() < 1e-12);
    assert_eq!(g["n_a"], "2");
    assert!(sidecar.contains("\"rule\": \"mean\""));
}
