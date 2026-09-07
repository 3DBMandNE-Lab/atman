//! `--weighted-whitening` has a scale ceiling, and it must refuse rather
//! than run overnight.
//!
//! The weighted path cannot use the sample-space Gram shortcut, because a
//! weighted inner product between two SAMPLES is not the weighted
//! covariance between two FEATURES. It therefore builds a dense p x p
//! covariance and eigendecomposes it at O(p^3).
//!
//! Measured on an M4, n=60, k=5, uniform weights: p=200 0.46 s, p=400
//! 5.22 s, p=800 26.5 s, p=1200 88.0 s. The fitted exponent between the
//! last two points is 2.96. Extrapolated to real cohort widths that is
//! about 16 hours at p=10,491 and 19 hours at p=10,977, on a covariance
//! matrix of roughly 1 GB.
//!
//! A flag that allocates a gigabyte and runs overnight with no warning is
//! the same defect class as a well-formed but meaningless output: it
//! looks like it is working. The GBM manuscript session confirmed the
//! flag has never been used -- zero of 23 decompose sidecars in the paper
//! tree pass any weighting -- so refusing costs no existing analysis.

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

    let mut st = 0x2545F491_4F6CDD1D_u64;
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

fn ica(dir: &Path, out: &Path, k: &str, weighted: bool) -> Output {
    let l = out.join("l.tsv");
    let a = out.join("a.tsv");
    let s = out.join("s.tsv");
    let mut args = vec![
        "decompose",
        "ica",
        "--input-dir",
        dir.to_str().unwrap(),
        "--k",
        k,
        "--max-iter",
        "20",
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
    if weighted {
        args.extend_from_slice(&["--missingness-model", "abundance-conditional"]);
        args.push("--weighted-whitening");
    }
    run_atman(&args)
}

/// Above the ceiling the command must refuse, and the message must carry
/// the numbers a user needs to act: the matrix size, the cost, and the
/// limit.
#[test]
fn weighted_whitening_refuses_above_the_assay_ceiling() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().join("cohort");
    // Just over the 2000-assay limit, with few samples so the fixture is
    // cheap to build. The refusal fires before any decomposition work.
    write_cohort(&dir, 8, 2001);
    let out = tmp.path().join("out");
    std::fs::create_dir_all(&out).unwrap();

    let res = ica(&dir, &out, "3", true);
    assert!(
        !res.status.success(),
        "2001 assays must be refused with --weighted-whitening"
    );
    let err = String::from_utf8_lossy(&res.stderr).to_string();
    for needle in ["--weighted-whitening", "2001", "2000", "p^3"] {
        assert!(
            err.contains(needle),
            "the refusal must name {needle:?} so the user can act on it; got:\n{err}"
        );
    }
}

/// Below the ceiling it must still run. A gate that refuses everything
/// would pass the test above while removing the feature.
#[test]
fn weighted_whitening_still_runs_below_the_ceiling() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().join("cohort");
    write_cohort(&dir, 20, 120);
    let out = tmp.path().join("out");
    std::fs::create_dir_all(&out).unwrap();

    let res = ica(&dir, &out, "3", true);
    assert!(
        res.status.success(),
        "120 assays is far below the ceiling and must run:\n{}",
        String::from_utf8_lossy(&res.stderr)
    );
    assert!(out.join("l.tsv").exists(), "loadings were not written");
}

/// The ceiling applies ONLY to the weighted path. Unweighted whitening
/// uses the sample-space Gram and is unaffected, so a wide cohort must
/// still decompose.
#[test]
fn the_ceiling_does_not_apply_to_unweighted_whitening() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().join("cohort");
    write_cohort(&dir, 8, 2001);
    let out = tmp.path().join("out");
    std::fs::create_dir_all(&out).unwrap();

    let res = ica(&dir, &out, "3", false);
    assert!(
        res.status.success(),
        "2001 assays without --weighted-whitening must still run:\n{}",
        String::from_utf8_lossy(&res.stderr)
    );
}
