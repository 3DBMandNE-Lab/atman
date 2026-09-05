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
fn axes_tree_builds_upgma_with_bootstrap_support() {
    let tmp = tempfile::tempdir().unwrap();
    let scores = tmp.path().join("scores.tsv");
    let mut s = String::from("sample_id\tcohort\tcondition\tis_control\ts1\ts2\n");
    // A: +2 on s1 ; B: +2 on s2 ; C: +1,+1 (45° from both).
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
    let linkage = tmp.path().join("linkage.tsv");
    let support = tmp.path().join("support.tsv");
    let newick = tmp.path().join("tree.newick");
    let r = run_atman(&[
        "axes",
        "tree",
        "--scores",
        scores.to_str().unwrap(),
        "--manifest",
        manifest.to_str().unwrap(),
        "--score-cols",
        "s1,s2",
        "--standardize",
        "none",
        "--n-bootstrap",
        "100",
        "--seed",
        "3",
        "--output-linkage",
        linkage.to_str().unwrap(),
        "--output-support",
        support.to_str().unwrap(),
        "--output-newick",
        newick.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let l = read_rows(&linkage);
    assert_eq!(l.len(), 2);
    // First merge: A (0) with C (2) at 1 − cos 45° = 0.2929; then B joins at the average of 1.0 and 0.2929.
    assert_eq!(l[0]["left"], "0");
    assert_eq!(l[0]["right"], "2");
    assert!(
        (l[0]["distance"].parse::<f64>().unwrap() - (1.0 - std::f64::consts::FRAC_1_SQRT_2)).abs()
            < 1e-9
    );
    assert_eq!(l[0]["n"], "2");
    assert_eq!(l[1]["left"], "1");
    assert_eq!(l[1]["right"], "3");
    assert!(
        (l[1]["distance"].parse::<f64>().unwrap()
            - (1.0 + (1.0 - std::f64::consts::FRAC_1_SQRT_2)) / 2.0)
            .abs()
            < 1e-9
    );
    let sp = read_rows(&support);
    assert_eq!(sp.len(), 2);
    assert_eq!(sp[0]["node"], "A | C");
    assert_eq!(sp[0]["n_leaves"], "2");
    let s0: f64 = sp[0]["support"].parse().unwrap();
    assert!((0.0..=1.0).contains(&s0));
    assert_eq!(sp[1]["support"], "1");
    assert_eq!(sp[0]["n_boot"], "100");
    let nw = std::fs::read_to_string(&newick).unwrap();
    let label = (s0 * 100.0).round() as i64;
    assert_eq!(
        nw.trim_end(),
        format!("(B:0.646447,(A:0.292893,C:0.292893){label}:0.353553)100;")
    );
    assert!(tmp.path().join("linkage.tsv.run.json").exists());
}
