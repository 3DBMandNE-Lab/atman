use std::collections::HashMap;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
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

fn setup(tmp: &std::path::Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let scores = tmp.join("scores.tsv");
    let mut s = String::from("sample_id\tcohort\tcondition\tis_control\ts1\n");
    let mut m = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tQAlb\tage\n",
    );
    for i in 0..12 {
        let cond = if i % 2 == 0 { "MS" } else { "nonMS" };
        let is_control = if i % 2 == 0 { 0 } else { 1 };
        s.push_str(&format!("P{i}\tms\t{cond}\t{is_control}\t{}\n", i as f64));
        m.push_str(&format!(
            "P{i}\tP{i}\t{cond}\t{is_control}\tbio\t{}\t{}\t{}\n",
            i + 1,
            10f64.powf(i as f64 / 4.0),
            100 - i
        ));
    }
    std::fs::write(&scores, s).unwrap();
    let dir = tmp.join("ms");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("samples.tsv"), m).unwrap();
    (scores, dir)
}

#[test]
fn axes_anchor_computes_spearman_per_scope_with_bootstrap_ci() {
    let tmp = tempfile::tempdir().unwrap();
    let (scores, dir) = setup(tmp.path());
    let out = tmp.path().join("anchors.tsv");
    let r = run_atman(&[
        "axes",
        "anchor",
        "--scores",
        scores.to_str().unwrap(),
        "--cohort-dirs",
        &format!("ms={}", dir.display()),
        "--cohorts",
        "ms",
        "--score-cols",
        "s1",
        "--anchors",
        "log10(QAlb),age",
        "--scope",
        "all",
        "--scope",
        "controls",
        "--scope",
        "condition=MS",
        "--n-bootstrap",
        "100",
        "--seed",
        "3",
        "--min-n",
        "5",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&out);
    assert_eq!(rows.len(), 6);
    let all_q = rows
        .iter()
        .find(|r| r["scope"] == "all" && r["anchor"] == "log10(QAlb)")
        .unwrap();
    assert_eq!(all_q["n"], "12");
    assert_eq!(all_q["rho"], "1");
    assert_eq!(all_q["p"], "0");
    assert!((all_q["ci_lo"].parse::<f64>().unwrap() - 1.0).abs() < 1e-12);
    let ctrl_age = rows
        .iter()
        .find(|r| r["scope"] == "controls" && r["anchor"] == "age")
        .unwrap();
    assert_eq!(ctrl_age["n"], "6");
    assert_eq!(ctrl_age["rho"], "-1");
    let ms_age = rows
        .iter()
        .find(|r| r["scope"] == "condition=MS" && r["anchor"] == "age")
        .unwrap();
    assert_eq!(ms_age["n"], "6");
    assert!(rows.iter().all(|r| r["given"].is_empty()));
    assert!(tmp.path().join("anchors.tsv.run.json").exists());
}

#[test]
fn axes_anchor_partial_and_min_n() {
    let tmp = tempfile::tempdir().unwrap();
    let (scores, dir) = setup(tmp.path());
    let out = tmp.path().join("anchors.tsv");
    let r = run_atman(&[
        "axes",
        "anchor",
        "--scores",
        scores.to_str().unwrap(),
        "--cohort-dirs",
        &format!("ms={}", dir.display()),
        "--score-cols",
        "s1",
        "--anchors",
        "age",
        "--scope",
        "all",
        "--partial",
        "log10(QAlb)",
        "--n-bootstrap",
        "0",
        "--min-n",
        "20",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&out);
    assert_eq!(rows.len(), 2);
    let plain = rows.iter().find(|r| r["given"].is_empty()).unwrap();
    assert_eq!(plain["n"], "12");
    assert!(plain["rho"].is_empty());
    let partial = rows.iter().find(|r| !r["given"].is_empty()).unwrap();
    assert_eq!(partial["scope"], "all_partial");
    assert_eq!(partial["given"], "log10(QAlb)");
    assert!(partial["rho"].is_empty());
}
