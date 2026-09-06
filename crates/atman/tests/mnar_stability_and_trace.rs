//! Two MNAR-ICA reporting gaps, both reported from a real run tree.
//!
//! 1. The seed-stability column ran its alternative seeds through plain
//!    FastICA even when the reference was the missingness-aware fit, so
//!    `n_stable_runs` compared two different methods at a Jaccard
//!    threshold. In practice it was exactly zero for every program in
//!    every cohort — a column named for a measurement that could not
//!    report one. Alternative seeds now run the reference's method.
//! 2. The joint loop reported only `joint_converged`, with the
//!    tolerance hard-coded. A run that stopped at `--max-joint-iter`
//!    could not be described as closing, oscillating, or diverging.
//!    `--joint-tol` is now exposed and the per-iteration betas are
//!    written to `mnar_joint_trace.tsv`.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_atman"))
        .args(args)
        .output()
        .expect("run atman")
}

/// Fixture with genuine below-LOD dropout, so the MNAR path has
/// something to model, and with block structure so programs are
/// recoverable across seeds.
fn write_mnar_fixture(dir: &Path, n_samples: usize, n_genes: usize) {
    std::fs::create_dir_all(dir).unwrap();
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=n_samples {
        samples.push_str(&format!("S{i:03}\tS{i:03}\tCase\t0\tplasma\t{i}\n"));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 1..=n_genes {
        proteins.push_str(&format!(
            "olink_explore_ngs\tA{j:04}\tQ{j:05}\tG{j:04}\tP1\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    let mut m = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0u64;
    for i in 1..=n_samples {
        // Two planted source signals drive two blocks of genes.
        let s1 = ((i as f64) * 0.83).sin() * 2.0;
        let s2 = ((i as f64) * 1.27).cos() * 2.0;
        for j in 1..=n_genes {
            order += 1;
            let block = if j * 3 <= n_genes {
                s1
            } else if j * 3 <= 2 * n_genes {
                s2
            } else {
                0.0
            };
            let v = 6.0 + block + (((i * 7 + j * 5) % 13) as f64) * 0.02;
            // Low values drop out, which is what makes the missingness
            // abundance-conditional rather than random.
            let below = v < 5.4;
            // Undetected cells are written as NaN, which is how
            // `load_matrix_mnar` represents a non-detection.
            let ab = if below {
                "NaN".to_string()
            } else {
                format!("{v:.6}")
            };
            let lod_flag = u8::from(below);
            m.push_str(&format!(
                "olink_explore_ngs\tS{i:03}\tA{j:04}\tG{j:04}\tP1\t{ab}\t{ab}\t{ab}\t\
                 log2_npx\tPASS\tPASS\t5.4\t{lod_flag}\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), m).unwrap();
}

fn run_mnar(input: &Path, output: &Path, extra: &[&str]) -> Output {
    std::fs::create_dir_all(output).unwrap();
    let mut args: Vec<String> = vec![
        "decompose".into(),
        "ica".into(),
        "--input-dir".into(),
        input.display().to_string(),
        "--output-loadings".into(),
        output.join("loadings.tsv").display().to_string(),
        "--output-activations".into(),
        output.join("activations.tsv").display().to_string(),
        "--output-stability".into(),
        output.join("stability.tsv").display().to_string(),
        "--k".into(),
        "2".into(),
        "--seed".into(),
        "42".into(),
        "--missingness-model".into(),
        "abundance-conditional".into(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    let refs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
    run_atman(&refs)
}

fn column(path: &Path, name: &str) -> Vec<String> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let idx = lines
        .next()
        .unwrap()
        .split('\t')
        .position(|h| h == name)
        .unwrap_or_else(|| panic!("missing column {name} in {}", path.display()));
    lines
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').nth(idx).unwrap_or("").to_string())
        .collect()
}

#[test]
fn mnar_seed_stability_can_be_non_zero() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("out");
    write_mnar_fixture(&input, 24, 30);

    let out = run_mnar(
        &input,
        &output,
        &["--n-seeds", "3", "--max-joint-iter", "3"],
    );
    assert!(
        out.status.success(),
        "mnar ica failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let n_stable = column(&output.join("stability.tsv"), "n_stable_runs");
    assert!(!n_stable.is_empty(), "stability table has no rows");
    let total: usize = n_stable.iter().map(|v| v.parse::<usize>().unwrap()).sum();
    assert!(
        total > 0,
        "under MNAR every program reported n_stable_runs = 0, which is the bug: the \
         alternative seeds must run the same method as the reference. Rows: {n_stable:?}"
    );

    // And the sidecar must say which method the alt seeds used, so the
    // number is interpretable without reading the source.
    let j: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(output.join("loadings.tsv.run.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(j["stability_alt_seed_method"], "abundance-conditional");
}

#[test]
fn mnar_writes_a_per_iteration_joint_trace() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("out");
    write_mnar_fixture(&input, 24, 30);

    let out = run_mnar(
        &input,
        &output,
        &[
            "--n-seeds",
            "1",
            "--max-joint-iter",
            "4",
            "--joint-tol",
            "1e-12",
        ],
    );
    assert!(
        out.status.success(),
        "mnar ica failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let trace = output.join("mnar_joint_trace.tsv");
    assert!(trace.is_file(), "no mnar_joint_trace.tsv was written");
    let text = std::fs::read_to_string(&trace).unwrap();
    let header = text.lines().next().unwrap();
    assert_eq!(header, "iteration\tbeta0\tbeta1\tdelta\tconverged_here");
    let rows: Vec<&str> = text.lines().skip(1).filter(|l| !l.is_empty()).collect();
    assert!(!rows.is_empty(), "trace has no iterations");
    // First iteration has no predecessor, so no step size.
    assert!(
        rows[0].contains("\tNA"),
        "the first iteration should record delta = NA, got {:?}",
        rows[0]
    );
    // The trace is listed as an output so its hash is recorded.
    let j: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(output.join("loadings.tsv.run.json")).unwrap(),
    )
    .unwrap();
    let named: Vec<String> = j["output_files"]
        .as_object()
        .unwrap()
        .keys()
        .map(|k| {
            Path::new(k)
                .file_name()
                .unwrap()
                .to_string_lossy()
                .to_string()
        })
        .collect();
    assert!(
        named.iter().any(|n| n == "mnar_joint_trace.tsv"),
        "the trace must be a recorded output; got {named:?}"
    );
}

#[test]
fn mnar_joint_tol_is_honoured() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    write_mnar_fixture(&input, 24, 30);

    // A tolerance so loose that the loop must stop at iteration 2.
    let loose = tmp.path().join("loose");
    let out = run_mnar(
        &input,
        &loose,
        &[
            "--n-seeds",
            "1",
            "--max-joint-iter",
            "6",
            "--joint-tol",
            "1e9",
        ],
    );
    assert!(out.status.success());
    let loose_rows = column(&loose.join("mnar_joint_trace.tsv"), "iteration").len();
    assert_eq!(
        loose_rows, 2,
        "a huge --joint-tol should converge at the first comparable step"
    );
}

/// The alternative seeds exist only to rank stability. Changing how many
/// there are, or which method they use, must not move the reference
/// fit's loadings or activations by a single byte.
///
/// This is the invariant that bounds a caller's exposure when the
/// stability machinery changes: if it holds, a pinned tree keeps its
/// loadings and only its stability column comes from the newer commit.
/// Verified across this session's changes against the pre-session
/// binary; pinned here so it stays true.
#[test]
fn mnar_reference_output_does_not_depend_on_n_seeds() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    write_mnar_fixture(&input, 24, 30);

    let one = tmp.path().join("one");
    let many = tmp.path().join("many");
    for (out_dir, n_seeds) in [(&one, "1"), (&many, "4")] {
        let out = run_mnar(
            &input,
            out_dir,
            &["--n-seeds", n_seeds, "--max-joint-iter", "3"],
        );
        assert!(
            out.status.success(),
            "mnar ica failed at --n-seeds {n_seeds}:\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    for f in ["loadings.tsv", "activations.tsv"] {
        let a = std::fs::read(one.join(f)).unwrap();
        let b = std::fs::read(many.join(f)).unwrap();
        assert_eq!(
            a, b,
            "{f} changed with --n-seeds; the reference fit must be independent of the \
             stability machinery"
        );
    }
}
