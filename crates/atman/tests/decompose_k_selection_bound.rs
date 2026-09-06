//! An automatic `k`-selection rule that returns its own search ceiling
//! has not selected anything — it ran out of room, and the criterion it
//! was asked to satisfy may never have been met.
//!
//! This was found in the wild: a six-cohort `decompose ica` run with
//! `--k-selection cumulative-variance=0.80 --k-max 30` resolved to
//! exactly 30 in five of six cohorts. The real requirements were 28 to
//! 56, so those five captured 70-79% variance while the run
//! configuration claimed 80%. The sidecar recorded `k_resolved: 30`
//! beside `k-max: 30`, which is indistinguishable from a genuine
//! selection landing on 30.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_atman"))
        .args(args)
        .output()
        .expect("run atman")
}

/// Correlated multi-factor matrix whose variance spectrum is spread
/// widely enough that an 80% cumulative-variance rule needs many
/// components.
fn write_wide_spectrum_fixture(dir: &Path, n_samples: usize, n_assays: usize) {
    std::fs::create_dir_all(dir).unwrap();
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=n_samples {
        samples.push_str(&format!("S{i:03}\tS{i:03}\tCase\t0\tplasma\t{i}\n"));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 1..=n_assays {
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
        for j in 1..=n_assays {
            order += 1;
            // Many weak, near-equal factors: no small set of components
            // captures 80% of the variance.
            let mut v = 0.0_f64;
            for f in 1..=24u64 {
                let phase = ((i as f64) * (f as f64) * 0.7).sin();
                let load = ((j as f64) * (f as f64) * 1.3).cos();
                v += phase * load;
            }
            v = 6.0 + v / 6.0;
            m.push_str(&format!(
                "olink_explore_ngs\tS{i:03}\tA{j:04}\tG{j:04}\tP1\t{v:.6}\t{v:.6}\t{v:.6}\t\
                 log2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), m).unwrap();
}

fn sidecar(output_dir: &Path, file: &str) -> serde_json::Value {
    let p = output_dir.join(format!("{file}.run.json"));
    serde_json::from_str(
        &std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display())),
    )
    .expect("parse sidecar")
}

#[test]
fn ica_k_selection_that_hits_k_max_warns_and_is_recorded() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("out");
    std::fs::create_dir_all(&output).unwrap();
    write_wide_spectrum_fixture(&input, 40, 60);

    // A ceiling far below what an 80% rule needs on this spectrum.
    let out = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-loadings",
        output.join("loadings.tsv").to_str().unwrap(),
        "--output-activations",
        output.join("activations.tsv").to_str().unwrap(),
        "--output-stability",
        output.join("stability.tsv").to_str().unwrap(),
        "--k-selection",
        "cumulative-variance=0.80",
        "--k-min",
        "2",
        "--k-max",
        "3",
        "--seed",
        "42",
    ]);
    assert!(
        out.status.success(),
        "ica run failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--k-max") && stderr.contains("truncation, not a selection"),
        "hitting --k-max must warn that the criterion was never met; got:\n{stderr}"
    );

    let j = sidecar(&output, "loadings.tsv");
    assert_eq!(j["k_resolved"], 3, "fixture should saturate the ceiling");
    assert_eq!(
        j["k_selection_bound_hit"], "k_max",
        "the sidecar must distinguish a truncation from a selection"
    );
}

#[test]
fn ica_k_selection_inside_its_range_is_silent_and_records_null() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("out");
    std::fs::create_dir_all(&output).unwrap();
    write_wide_spectrum_fixture(&input, 40, 60);

    let out = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-loadings",
        output.join("loadings.tsv").to_str().unwrap(),
        "--output-activations",
        output.join("activations.tsv").to_str().unwrap(),
        "--output-stability",
        output.join("stability.tsv").to_str().unwrap(),
        // 0.70 resolves to k=3 on this fixture: strictly inside
        // [2, 30], so neither bound is touched.
        "--k-selection",
        "cumulative-variance=0.70",
        "--k-min",
        "2",
        "--k-max",
        "30",
        "--seed",
        "42",
    ]);
    assert!(
        out.status.success(),
        "ica run failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !stderr.contains("truncation, not a selection"),
        "a rule that resolved inside its range must not warn; got:\n{stderr}"
    );
    let j = sidecar(&output, "loadings.tsv");
    // The field must be present-and-null, not absent: an absent field
    // would make this assertion pass for the wrong reason.
    assert!(
        j.as_object().unwrap().contains_key("k_selection_bound_hit"),
        "sidecar is missing k_selection_bound_hit entirely"
    );
    assert!(
        j["k_selection_bound_hit"].is_null(),
        "expected null bound-hit, got {}",
        j["k_selection_bound_hit"]
    );
}

#[test]
fn explicit_k_is_not_reported_as_a_bound_hit() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("out");
    std::fs::create_dir_all(&output).unwrap();
    write_wide_spectrum_fixture(&input, 40, 60);

    let out = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-loadings",
        output.join("loadings.tsv").to_str().unwrap(),
        "--output-activations",
        output.join("activations.tsv").to_str().unwrap(),
        "--output-stability",
        output.join("stability.tsv").to_str().unwrap(),
        "--k",
        "3",
        "--seed",
        "42",
    ]);
    assert!(
        out.status.success(),
        "ica run failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let j = sidecar(&output, "loadings.tsv");
    assert!(
        j.as_object().unwrap().contains_key("k_selection_bound_hit"),
        "sidecar is missing k_selection_bound_hit entirely"
    );
    assert!(
        j["k_selection_bound_hit"].is_null(),
        "an explicit --k is not a selection and must record null"
    );
}
