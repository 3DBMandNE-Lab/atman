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

fn num(r: &HashMap<String, String>, c: &str) -> f64 {
    r[c].parse().unwrap_or_else(|_| panic!("{c} = {:?}", r[c]))
}

#[test]
fn axes_displacement_cosines_norms_and_bootstrap() {
    let tmp = tempfile::tempdir().unwrap();
    let scores = tmp.path().join("scores.tsv");
    let mut s = String::from("sample_id\tcohort\tcondition\tis_control\ts1\ts2\n");
    for (cohort, d1, d2) in [("c1", 2.0, 0.0), ("c2", 0.0, 2.0), ("c3", 1.0, 1.0)] {
        for (i, base) in [(0.0, 0.0), (1.0, -1.0), (-1.0, 1.0), (0.5, 0.5)]
            .iter()
            .enumerate()
        {
            s.push_str(&format!(
                "{cohort}_case{i}\t{cohort}\tCase\t0\t{}\t{}\n",
                base.0 + d1,
                base.1 + d2
            ));
            s.push_str(&format!(
                "{cohort}_ctrl{i}\t{cohort}\tCtrl\t1\t{}\t{}\n",
                base.0, base.1
            ));
        }
    }
    std::fs::write(&scores, s).unwrap();
    let manifest = tmp.path().join("contrasts.tsv");
    std::fs::write(
        &manifest,
        "label\tcohort\tcase\tcontrol\tfamily\nA\tc1\tCase\tCtrl\tf\nB\tc2\tCase\tCtrl\tf\nC\tc3\tCase\tCtrl\tf\n",
    )
    .unwrap();
    let out = tmp.path().join("cos.tsv");
    let vec = tmp.path().join("vec.tsv");
    let r = run_atman(&[
        "axes",
        "displacement",
        "--scores",
        scores.to_str().unwrap(),
        "--manifest",
        manifest.to_str().unwrap(),
        "--score-sets",
        "two=s1,s2;one=s1",
        "--standardize",
        "none",
        "--n-bootstrap",
        "100",
        "--seed",
        "9",
        "--output",
        out.to_str().unwrap(),
        "--output-vectors",
        vec.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&out);
    assert_eq!(rows.len(), 6);
    let ab = rows
        .iter()
        .find(|r| r["a"] == "A" && r["b"] == "B" && r["space"] == "two")
        .unwrap();
    assert!(num(ab, "cosine").abs() < 1e-12);
    assert!((num(ab, "norm_a") - 2.0).abs() < 1e-12);
    assert!((num(ab, "norm_b") - 2.0).abs() < 1e-12);
    let ac = rows
        .iter()
        .find(|r| r["a"] == "A" && r["b"] == "C" && r["space"] == "two")
        .unwrap();
    assert!((num(ac, "cosine") - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-12);
    assert!((num(ac, "norm_b") - std::f64::consts::SQRT_2).abs() < 1e-12);
    assert!(
        num(ac, "ci_lo") <= num(ac, "cosine") + 1e-9
            && num(ac, "cosine") - 1e-9 <= num(ac, "ci_hi")
    );
    let v = read_rows(&vec);
    let a_s1 = v
        .iter()
        .find(|r| r["label"] == "A" && r["space"] == "two" && r["score"] == "s1")
        .unwrap();
    assert!((num(a_s1, "displacement") - 2.0).abs() < 1e-12);
    assert_eq!(a_s1["n_case"], "4");
    assert!(tmp.path().join("cos.tsv.run.json").exists());
}

#[test]
fn axes_displacement_global_standardization_equalizes_norms() {
    let tmp = tempfile::tempdir().unwrap();
    let scores = tmp.path().join("scores.tsv");
    std::fs::write(&scores, "sample_id\tcohort\tcondition\tis_control\ts1\ts2\nA1\tc\tCase\t0\t10\t1\nA2\tc\tCase\t0\t12\t1.2\nB1\tc\tCtrl\t1\t0\t0\nB2\tc\tCtrl\t1\t2\t0.2\n").unwrap();
    let manifest = tmp.path().join("m.tsv");
    std::fs::write(
        &manifest,
        "label\tcohort\tcase\tcontrol\tfamily\nX\tc\tCase\tCtrl\tf\nY\tc\tCase\tCtrl\tf\n",
    )
    .unwrap();
    let out = tmp.path().join("cos.tsv");
    let r = run_atman(&[
        "axes",
        "displacement",
        "--scores",
        scores.to_str().unwrap(),
        "--manifest",
        manifest.to_str().unwrap(),
        "--score-sets",
        "two=s1,s2",
        "--n-bootstrap",
        "0",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&out);
    assert!((num(&rows[0], "cosine") - 1.0).abs() < 1e-12);
    assert!((num(&rows[0], "norm_a") - num(&rows[0], "norm_b")).abs() < 1e-12);
    assert!(rows[0]["ci_lo"].is_empty());
}
