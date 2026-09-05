use std::collections::HashMap;
use std::path::PathBuf;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/axes_contrast")
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

fn num(row: &HashMap<String, String>, col: &str) -> f64 {
    row[col]
        .parse::<f64>()
        .unwrap_or_else(|_| panic!("column {col} = {:?} not numeric", row[col]))
}

fn int_text(value: &str) -> String {
    value.trim_end_matches(".0000000000").to_string()
}

#[test]
fn axes_contrast_matches_python_reference() {
    let f = fixtures();
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("modulation.tsv");
    let cov = tmp.path().join("covariates.tsv");
    let result = run_atman(&[
        "axes",
        "contrast",
        "--scores",
        f.join("scores.tsv").to_str().unwrap(),
        "--cohort-dirs",
        &format!(
            "c1={},c2={}",
            f.join("c1").display(),
            f.join("c2").display()
        ),
        "--manifest",
        f.join("contrasts.tsv").to_str().unwrap(),
        "--score-cols",
        "s1,s2",
        "--designs",
        "~ case; ~ case + z(age) + sex",
        "--output",
        out.to_str().unwrap(),
        "--output-covariates",
        cov.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let got = read_rows(&out);
    let reference = read_rows(&f.join("reference.tsv"));
    assert_eq!(got.len(), reference.len(), "row count");
    for r in &reference {
        let g = got
            .iter()
            .find(|g| {
                g["label"] == r["label"] && g["score"] == r["score"] && g["design"] == r["design"]
            })
            .unwrap_or_else(|| panic!("missing row {:?}", r));
        for col in ["n_case", "n_control", "n_fit"] {
            assert_eq!(g[col], int_text(&r[col]), "{col} for {:?}", r);
        }
        for col in [
            "mean_diff",
            "cohen_d",
            "d_ci_lo",
            "d_ci_hi",
            "welch_t",
            "welch_p",
            "auc",
            "beta",
            "se",
            "ci_lo",
            "ci_hi",
            "p",
            "beta_over_resid_sd",
        ] {
            let (a, b) = (num(g, col), num(r, col));
            assert!(
                (a - b).abs() < 1e-6,
                "{col}: atman {a} vs python {b} for {:?}",
                r
            );
        }
    }
    for g in &got {
        assert!(num(g, "q") >= num(g, "p") - 1e-12);
        assert!(num(g, "welch_q") >= num(g, "welch_p") - 1e-12);
        assert!(g["boot_n"].is_empty());
    }
    for g in got.iter().filter(|g| g["design"] == "~ case") {
        assert!((num(g, "beta") - num(g, "mean_diff")).abs() < 1e-9);
    }
    let covs = read_rows(&cov);
    assert!(covs.iter().any(|c| c["term"] == "z(age)"
        && c["design"] == "~ case + z(age) + sex"
        && c["kind"] == "coef"));
    assert!(covs.iter().any(|c| c["term"] == "sexM"));
    assert!(
        std::fs::read_to_string(tmp.path().join("modulation.tsv.run.json"))
            .unwrap()
            .contains("\"axes contrast\"")
    );
}

#[test]
fn axes_contrast_multilevel_factor_with_omnibus_f_test() {
    let f = fixtures();
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("m.tsv");
    let cov = tmp.path().join("c.tsv");
    let r = run_atman(&[
        "axes",
        "contrast",
        "--scores",
        f.join("scores.tsv").to_str().unwrap(),
        "--cohort-dirs",
        &format!(
            "c1={},c2={}",
            f.join("c1").display(),
            f.join("c2").display()
        ),
        "--manifest",
        f.join("contrasts.tsv").to_str().unwrap(),
        "--score-cols",
        "s1",
        "--designs",
        "~ case + z(age) + sex + site",
        "--omnibus-factor",
        "site",
        "--output",
        out.to_str().unwrap(),
        "--output-covariates",
        cov.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let reference = &read_rows(&f.join("reference_site.tsv"))[0];
    let row = read_rows(&out)
        .into_iter()
        .find(|g| g["label"] == "c1_case")
        .unwrap();
    assert!((num(&row, "beta") - num(reference, "beta_case")).abs() < 1e-6);
    assert!((num(&row, "p") - num(reference, "p_case")).abs() < 1e-6);
    assert_eq!(row["n_fit"], int_text(&reference["n_fit"]));
    let covs = read_rows(&cov);
    let coef_terms: Vec<&str> = covs
        .iter()
        .filter(|c| c["label"] == "c1_case" && c["kind"] == "coef")
        .map(|c| c["term"].as_str())
        .collect();
    assert_eq!(coef_terms, vec!["z(age)", "sexM", "siteKiel", "siteSweden"]);
    let omni = covs
        .iter()
        .find(|c| c["label"] == "c1_case" && c["kind"] == "omnibus_f")
        .unwrap();
    assert_eq!(omni["term"], "site");
    assert!((num(omni, "f") - num(reference, "f")).abs() < 1e-6);
    assert_eq!(num(omni, "df_num"), num(reference, "df_num"));
    assert_eq!(num(omni, "df_den"), num(reference, "df_den"));
    assert!((num(omni, "p") - num(reference, "p_f")).abs() < 1e-6);
    assert!(omni["beta"].is_empty());
}

#[test]
fn axes_contrast_bootstrap_is_deterministic_and_brackets_beta() {
    let f = fixtures();
    let run = |dir: &std::path::Path| {
        let out = dir.join("m.tsv");
        let boot = dir.join("b.tsv");
        let r = run_atman(&[
            "axes",
            "contrast",
            "--scores",
            f.join("scores.tsv").to_str().unwrap(),
            "--cohort-dirs",
            &format!(
                "c1={},c2={}",
                f.join("c1").display(),
                f.join("c2").display()
            ),
            "--manifest",
            f.join("contrasts.tsv").to_str().unwrap(),
            "--score-cols",
            "s1",
            "--designs",
            "~ case + z(age) + sex",
            "--n-bootstrap",
            "200",
            "--seed",
            "20260905",
            "--output",
            out.to_str().unwrap(),
            "--output-bootstrap",
            boot.to_str().unwrap(),
        ]);
        assert!(
            r.status.success(),
            "stderr:\n{}",
            String::from_utf8_lossy(&r.stderr)
        );
        (
            std::fs::read_to_string(&out).unwrap(),
            std::fs::read_to_string(&boot).unwrap(),
        )
    };
    let t1 = tempfile::tempdir().unwrap();
    let t2 = tempfile::tempdir().unwrap();
    let (m1, b1) = run(t1.path());
    let (m2, b2) = run(t2.path());
    assert_eq!(m1, m2);
    assert_eq!(b1, b2);
    for g in read_rows(&t1.path().join("m.tsv")) {
        assert_eq!(g["boot_n"], "200");
        let (lo, hi, beta) = (
            num(&g, "boot_ci_lo"),
            num(&g, "boot_ci_hi"),
            num(&g, "beta"),
        );
        assert!(
            lo < beta && beta < hi,
            "boot CI [{lo}, {hi}] must bracket beta {beta}"
        );
        let frac = num(&g, "boot_frac_gt0");
        assert!((0.0..=1.0).contains(&frac));
    }
    assert_eq!(read_rows(&t1.path().join("b.tsv")).len(), 400);
}

#[test]
fn axes_contrast_rejects_empty_group_and_missing_case_term() {
    let f = fixtures();
    let tmp = tempfile::tempdir().unwrap();
    let manifest = tmp.path().join("bad.tsv");
    std::fs::write(
        &manifest,
        "label\tcohort\tcase\tcontrol\tfamily\nx\tc1\tNope\tCtrl\tfam\n",
    )
    .unwrap();
    let r = run_atman(&[
        "axes",
        "contrast",
        "--scores",
        f.join("scores.tsv").to_str().unwrap(),
        "--manifest",
        manifest.to_str().unwrap(),
        "--score-cols",
        "s1",
        "--output",
        tmp.path().join("o.tsv").to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("no case rows"));
    let r = run_atman(&[
        "axes",
        "contrast",
        "--scores",
        f.join("scores.tsv").to_str().unwrap(),
        "--manifest",
        f.join("contrasts.tsv").to_str().unwrap(),
        "--score-cols",
        "s1",
        "--designs",
        "~ age",
        "--output",
        tmp.path().join("o.tsv").to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("`case`"));
}

#[test]
fn axes_contrast_bootstrap_skips_replicates_with_a_degenerate_factor() {
    let tmp = tempfile::tempdir().unwrap();
    let scores = tmp.path().join("scores.tsv");
    let dir = tmp.path().join("c");
    std::fs::create_dir_all(&dir).unwrap();
    let mut s = String::from("sample_id\tcohort\tcondition\tis_control\ts1\n");
    let mut m = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tsex\n",
    );
    // Three cases and three controls; exactly one case is male.
    for (i, (cond, ctrl, sex, v)) in [
        ("Case", 0, "M", 2.0),
        ("Case", 0, "F", 1.5),
        ("Case", 0, "F", 2.5),
        ("Ctrl", 1, "F", 0.0),
        ("Ctrl", 1, "F", 0.5),
        ("Ctrl", 1, "F", -0.5),
    ]
    .iter()
    .enumerate()
    {
        s.push_str(&format!("P{i}\tc\t{cond}\t{ctrl}\t{v}\n"));
        m.push_str(&format!(
            "P{i}\tP{i}\t{cond}\t{ctrl}\tbio\t{}\t{sex}\n",
            i + 1
        ));
    }
    std::fs::write(&scores, s).unwrap();
    std::fs::write(dir.join("samples.tsv"), m).unwrap();
    let manifest = tmp.path().join("m.tsv");
    std::fs::write(
        &manifest,
        "label\tcohort\tcase\tcontrol\tfamily\nX\tc\tCase\tCtrl\tf\n",
    )
    .unwrap();
    let out = tmp.path().join("o.tsv");
    let r = run_atman(&[
        "axes",
        "contrast",
        "--scores",
        scores.to_str().unwrap(),
        "--cohort-dirs",
        &format!("c={}", dir.display()),
        "--manifest",
        manifest.to_str().unwrap(),
        "--score-cols",
        "s1",
        "--designs",
        "~ case + sex",
        "--n-bootstrap",
        "50",
        "--seed",
        "1",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let row = &read_rows(&out)[0];
    let kept: usize = row["boot_n"].parse().unwrap();
    let skipped: usize = row["boot_n_skipped"].parse().unwrap();
    assert_eq!(kept + skipped, 50);
    assert!(skipped >= 1, "expected at least one degenerate replicate");
    assert!(kept >= 1);
}
