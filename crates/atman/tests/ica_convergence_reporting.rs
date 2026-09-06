//! FastICA that exhausts `--max-iter` without converging must say so.
//!
//! Fourth instance of a pattern found three times in real run trees
//! (`--ridge-lambda auto`, k-selection hitting its ceiling, soft-power
//! falling back to beta = 1): the output is well-formed, the recorded
//! parameters are faithful to what was requested, and the failure is
//! invisible. `IcaResult` has carried `n_iterations` and `final_tol`
//! all along; nothing read them, so a run that stopped at the iteration
//! cap produced loadings indistinguishable from a converged one.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_atman"))
        .args(args)
        .output()
        .expect("run atman")
}

fn write_fixture(dir: &Path, n_samples: usize, n_genes: usize) {
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
        for j in 1..=n_genes {
            order += 1;
            let v = 6.0
                + ((i as f64) * 0.7).sin() * ((j as f64) * 0.3).cos()
                + ((i as f64) * (j as f64) * 0.11).sin() * 0.5;
            m.push_str(&format!(
                "olink_explore_ngs\tS{i:03}\tA{j:04}\tG{j:04}\tP1\t{v:.6}\t{v:.6}\t{v:.6}\t\
                 log2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), m).unwrap();
}

fn run_ica(input: &Path, output: &Path, extra: &[&str]) -> Output {
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
        "3".into(),
        "--seed".into(),
        "42".into(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    let refs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
    run_atman(&refs)
}

fn sidecar(output: &Path) -> serde_json::Value {
    let p = output.join("loadings.tsv.run.json");
    serde_json::from_str(
        &std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display())),
    )
    .expect("parse sidecar")
}

#[test]
fn ica_stopping_at_max_iter_warns_and_is_recorded() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("out");
    write_fixture(&input, 30, 40);

    // One iteration and an impossible tolerance: cannot converge.
    let out = run_ica(
        &input,
        &output,
        &["--max-iter", "1", "--tol", "1e-300", "--n-seeds", "1"],
    );
    assert!(
        out.status.success(),
        "the run should still produce output:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("did not converge"),
        "hitting --max-iter must warn; got:\n{stderr}"
    );

    let j = sidecar(&output);
    assert_eq!(
        j["ica_converged"], false,
        "the sidecar must record that the fit did not converge"
    );
    assert_eq!(j["ica_n_iterations"], 1);
    assert!(
        j["ica_final_tol"].as_f64().unwrap() > 1e-300,
        "final tolerance should be recorded and above the requested tol"
    );
}

#[test]
fn ica_that_converges_is_silent_and_records_true() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("out");
    write_fixture(&input, 30, 40);

    let out = run_ica(
        &input,
        &output,
        &["--max-iter", "500", "--tol", "1e-4", "--n-seeds", "1"],
    );
    assert!(
        out.status.success(),
        "ica failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !stderr.contains("did not converge"),
        "a converged fit must not warn; got:\n{stderr}"
    );
    let j = sidecar(&output);
    assert!(
        j.as_object().unwrap().contains_key("ica_converged"),
        "sidecar is missing ica_converged entirely"
    );
    assert_eq!(j["ica_converged"], true);
}

#[test]
fn ica_reports_non_convergence_among_alternative_seeds() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("out");
    write_fixture(&input, 30, 40);

    // Stability seeds run the same solver; if they fail to converge the
    // stability ranking is built on unconverged fits and the user should
    // hear about it.
    let out = run_ica(
        &input,
        &output,
        &["--max-iter", "1", "--tol", "1e-300", "--n-seeds", "3"],
    );
    assert!(out.status.success());
    let j = sidecar(&output);
    assert_eq!(
        j["ica_n_seeds_not_converged"], 3,
        "all three seeds hit the iteration cap"
    );
}
