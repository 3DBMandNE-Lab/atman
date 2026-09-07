//! `--k-selection` must be validated, and must say when it is ignored.
//!
//! `resolve_k` returns early when `--k` is set, so the selection string
//! was never parsed on that path. Two consequences, both found by the
//! GBM manuscript session hitting them:
//!
//! 1. `--k 10 --k-selection fixed=30` ran happily with k=10 and wrote
//!    `"k-selection": "fixed=30"` into the run sidecar — a provenance
//!    record naming a rule that was never validated and never applied.
//!
//! 2. Without `--k`, the same string failed inside `resolve_k`, which
//!    runs AFTER the matrix is loaded. On a 123 MB cohort that is a
//!    3.5 s failure, and in a timing loop a 3.5 s failure is
//!    indistinguishable from a 3.5 s success. Four such runs were read
//!    as a flat max-iter slope and drawn a conclusion from. The
//!    conclusion was right by luck; the measurement was of a command
//!    that wrote no output.
//!
//! Validation now happens before the input is read, so a bad flag costs
//! milliseconds and any timing that includes one is obviously wrong.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_atman"))
        .args(args)
        .output()
        .expect("run atman")
}

fn write_cohort(dir: &Path, n: usize, p: usize) {
    std::fs::create_dir_all(dir).unwrap();
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 0..n {
        samples.push_str(&format!("S{i:04}\tS{i:04}\tN/A\t0\tplasma\t{}\n", i + 1));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 0..p {
        proteins.push_str(&format!(
            "olink_explore_ngs\tP{j:05}\tP{j:05}\tGP{j:05}\tP1\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    let mut st = 0x9E3779B97F4A7C15_u64;
    let mut rnd = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut m = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\t\
         abundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\t\
         dropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0usize;
    for i in 0..n {
        for j in 0..p {
            order += 1;
            let v = 5.0 + rnd();
            m.push_str(&format!(
                "olink_explore_ngs\tS{i:04}\tP{j:05}\tGP{j:05}\tP1\t{v:.6}\t{v:.6}\t{v:.6}\t\
                 npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), m).unwrap();
}

struct Run {
    _tmp: tempfile::TempDir,
    dir: std::path::PathBuf,
    out: std::path::PathBuf,
}

impl Run {
    fn new() -> Self {
        let tmp = tempfile::tempdir().unwrap();
        let dir = tmp.path().join("cohort");
        let out = tmp.path().join("out");
        write_cohort(&dir, 24, 120);
        std::fs::create_dir_all(&out).unwrap();
        Self {
            _tmp: tmp,
            dir,
            out,
        }
    }
    fn ica(&self, extra: &[&str]) -> Output {
        let l = self.out.join("l.tsv");
        let a = self.out.join("a.tsv");
        let s = self.out.join("s.tsv");
        let mut args = vec![
            "decompose",
            "ica",
            "--input-dir",
            self.dir.to_str().unwrap(),
            "--max-iter",
            "15",
            "--tol",
            "1e-3",
            "--seed",
            "1",
            "--n-seeds",
            "1",
            "--output-loadings",
            l.to_str().unwrap(),
            "--output-activations",
            a.to_str().unwrap(),
            "--output-stability",
            s.to_str().unwrap(),
        ];
        args.extend_from_slice(extra);
        run_atman(&args)
    }
    fn sidecar(&self) -> serde_json::Value {
        let t = std::fs::read_to_string(self.out.join("l.tsv.run.json")).unwrap();
        serde_json::from_str(&t).unwrap()
    }
}

/// An unsupported `--k-selection` must be refused even when `--k` makes
/// it irrelevant. Accepting it silently is what let a sidecar record a
/// rule that had never been parsed.
#[test]
fn an_unsupported_k_selection_is_refused_even_when_k_overrides_it() {
    let r = Run::new();
    let out = r.ica(&["--k", "4", "--k-selection", "fixed=30"]);
    assert!(
        !out.status.success(),
        "an unsupported --k-selection must be refused even with --k set"
    );
    let err = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(err.contains("unsupported --k-selection"), "got:\n{err}");
    assert!(err.contains("fixed=30"), "the error must quote the value");
}

/// The same value without `--k` must also be refused, and must NOT
/// require the input to be read first.
#[test]
fn an_unsupported_k_selection_is_refused_without_reading_the_input() {
    let r = Run::new();
    // Point at a directory with no measurements.tsv at all: if the flag
    // is validated first, the run still fails on the FLAG, not on the
    // missing input. That pins the ordering without timing anything.
    let empty = r.out.join("no_such_cohort");
    std::fs::create_dir_all(&empty).unwrap();
    let l = r.out.join("l.tsv");
    let a = r.out.join("a.tsv");
    let s = r.out.join("s.tsv");
    let out = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        empty.to_str().unwrap(),
        "--k-selection",
        "fixed=30",
        "--max-iter",
        "15",
        "--tol",
        "1e-3",
        "--seed",
        "1",
        "--n-seeds",
        "1",
        "--output-loadings",
        l.to_str().unwrap(),
        "--output-activations",
        a.to_str().unwrap(),
        "--output-stability",
        s.to_str().unwrap(),
    ]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(
        err.contains("unsupported --k-selection"),
        "the flag must be rejected BEFORE the input is opened, so the \
         error names the flag and not the missing file; got:\n{err}"
    );
}

/// When `--k` overrides a non-default rule the run must say so, and the
/// sidecar must record which one set k.
#[test]
fn k_overriding_a_rule_warns_and_is_recorded_in_the_sidecar() {
    let r = Run::new();
    let out = r.ica(&["--k", "4", "--k-selection", "cumulative-variance=0.5"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let err = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(
        err.contains("is IGNORED because --k"),
        "the run must say the rule was ignored; got:\n{err}"
    );
    let sc = r.sidecar();
    assert_eq!(sc["k_selection_applied"], serde_json::json!(false));
    assert_eq!(sc["k_resolved"], serde_json::json!(4));
}

/// And when the rule DOES set k, the sidecar must say that too --
/// otherwise the flag only ever records one of its two states.
#[test]
fn a_rule_that_actually_selects_k_is_recorded_as_applied() {
    let r = Run::new();
    let out = r.ica(&["--k-selection", "cumulative-variance=0.5"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let err = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(
        !err.contains("is IGNORED because --k"),
        "no --k was given, so nothing was ignored; got:\n{err}"
    );
    let sc = r.sidecar();
    assert_eq!(sc["k_selection_applied"], serde_json::json!(true));
}
