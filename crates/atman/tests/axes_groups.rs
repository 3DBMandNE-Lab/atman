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

#[test]
fn axes_groups_reports_means_ci_medians_kruskal_and_reference_effects() {
    let tmp = tempfile::tempdir().unwrap();
    let scores = tmp.path().join("scores.tsv");
    std::fs::write(
        &scores,
        "sample_id\tcohort\tcondition\tis_control\ts1\n\
         A1\tms\tMS\t0\t1\nA2\tms\tMS\t0\t2\nA3\tms\tMS\t0\t3\n\
         B1\tms\tnonMS\t1\t4\nB2\tms\tnonMS\t1\t5\nB3\tms\tnonMS\t1\t6\n\
         C1\tms\tnonMS\t1\t7\nC2\tms\tnonMS\t1\t8\nC3\tms\tnonMS\t1\t9\n\
         Z1\tother\tX\t0\t100\n",
    )
    .unwrap();
    let dir = tmp.path().join("ms");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tdiagnosis_group\tQAlb\tQIgG\n\
         A1\tA1\tMS\t0\tbio\t1\tMS\t4\t2\nA2\tA2\tMS\t0\tbio\t2\tMS\t6\t3\nA3\tA3\tMS\t0\tbio\t3\tMS\t8\t2\n\
         B1\tB1\tnonMS\t1\tbio\t4\tHeadache\t10\t5\nB2\tB2\tnonMS\t1\tbio\t5\tHeadache\t\t5\nB3\tB3\tnonMS\t1\tbio\t6\tHeadache\t20\t5\n\
         C1\tC1\tnonMS\t1\tbio\t7\tInfectious\t1\t1\nC2\tC2\tnonMS\t1\tbio\t8\tInfectious\t2\t1\nC3\tC3\tnonMS\t1\tbio\t9\tInfectious\t3\t1\n",
    )
    .unwrap();
    let out = tmp.path().join("groups.tsv");
    let tests = tmp.path().join("tests.tsv");
    let r = run_atman(&[
        "axes",
        "groups",
        "--scores",
        scores.to_str().unwrap(),
        "--cohort-dirs",
        &format!("ms={}", dir.display()),
        "--cohort",
        "ms",
        "--group-by",
        "diagnosis_group",
        "--groups",
        "MS,Headache,Infectious",
        "--score-cols",
        "s1",
        "--median-cols",
        "QAlb,QIgG/QAlb",
        "--reference",
        "MS",
        "--output",
        out.to_str().unwrap(),
        "--output-tests",
        tests.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&out);
    assert_eq!(
        rows.iter().map(|r| r["group"].as_str()).collect::<Vec<_>>(),
        vec!["MS", "Headache", "Infectious"]
    );
    let ms = &rows[0];
    assert_eq!(ms["n"], "3");
    assert_eq!(ms["s1_mean"], "2");
    assert!((ms["s1_ci_lo"].parse::<f64>().unwrap() - (2.0 - 2.4841377)).abs() < 1e-6);
    assert!((ms["s1_ci_hi"].parse::<f64>().unwrap() - (2.0 + 2.4841377)).abs() < 1e-6);
    assert_eq!(ms["QAlb_median"], "6");
    assert_eq!(ms["QIgG/QAlb_median"], "0.5");
    assert_eq!(rows[1]["QAlb_median"], "15");
    let t = read_rows(&tests);
    let kw = t
        .iter()
        .find(|r| r["kind"] == "kruskal" && r["score"] == "s1")
        .unwrap();
    assert!((kw["statistic"].parse::<f64>().unwrap() - 7.2).abs() < 1e-9);
    assert!((kw["p"].parse::<f64>().unwrap() - 0.02732372244729256).abs() < 1e-9);
    let vs = t
        .iter()
        .find(|r| r["kind"] == "reference_vs_group" && r["group"] == "Infectious")
        .unwrap();
    assert_eq!(vs["reference"], "MS");
    assert_eq!(vs["n_ref"], "3");
    assert_eq!(vs["n_group"], "3");
    assert!((vs["statistic"].parse::<f64>().unwrap() + 6.0).abs() < 1e-9);
    assert!(vs["p"].parse::<f64>().unwrap() < 0.01);
    assert!(tmp.path().join("groups.tsv.run.json").exists());
}

#[test]
fn axes_groups_rejects_unknown_group() {
    let tmp = tempfile::tempdir().unwrap();
    let scores = tmp.path().join("scores.tsv");
    std::fs::write(
        &scores,
        "sample_id\tcohort\tcondition\tis_control\ts1\nA1\tms\tMS\t0\t1\nA2\tms\tnonMS\t1\t2\n",
    )
    .unwrap();
    let r = run_atman(&[
        "axes",
        "groups",
        "--scores",
        scores.to_str().unwrap(),
        "--cohort",
        "ms",
        "--group-by",
        "condition",
        "--groups",
        "MS,Nope",
        "--score-cols",
        "s1",
        "--output",
        tmp.path().join("g.tsv").to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("Nope"));
}
