//! Guards against a degenerate module discovery that looks well-formed.
//!
//! Found in the wild by the atman methods-paper session: on a 93-sample
//! bulk-proteome contrast the scale-free criterion was met at no power
//! in 1..=20 (R² declined monotonically from 0.505 to 0.133), the sweep
//! fell back to β = 1 — no soft-thresholding at all — and every one of
//! 9,376 features landed in the grey catch-all. `module-de` then ran on
//! that K=1 "module set" and emitted a single row testing the mean of
//! all 9,376 proteins between subtypes, at q = 0.04. Exit code 0, a
//! well-formed output file, and a significant-looking number that means
//! nothing.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_atman"))
        .args(args)
        .output()
        .expect("run atman")
}

/// Dense, near-uniform correlation structure: no soft power yields a
/// scale-free degree distribution, and nothing separates into modules.
fn write_unstructured_fixture(dir: &Path, n_samples: usize, n_genes: usize) {
    std::fs::create_dir_all(dir).unwrap();
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=n_samples {
        let cond = if i * 2 <= n_samples { "A" } else { "B" };
        samples.push_str(&format!("S{i:03}\tS{i:03}\t{cond}\t0\tplasma\t{i}\n"));
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
            // Each gene follows its own incommensurate oscillation, so
            // no pair is reliably correlated and nothing clusters into a
            // module of the required size.
            let v = 6.0 + ((i as f64) * (j as f64) * 1.618).sin();
            m.push_str(&format!(
                "olink_explore_ngs\tS{i:03}\tA{j:04}\tG{j:04}\tP1\t{v:.6}\t{v:.6}\t{v:.6}\t\
                 log2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), m).unwrap();
}

#[test]
fn modules_discover_refuses_to_emit_an_all_grey_module_set() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("out");
    write_unstructured_fixture(&input, 24, 40);

    let out = run_atman(&[
        "modules",
        "discover",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--method",
        "wgcna-soft",
        "--soft-power",
        "auto",
        "--r2-target",
        "0.8",
        "--max-beta",
        "20",
        "--cut-height",
        "0.02",
        "--min-module-size",
        "8",
    ]);
    assert!(
        !out.status.success(),
        "an all-grey result must be refused, not written with exit 0"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("grey catch-all"),
        "the refusal must explain that no module was found; got:\n{stderr}"
    );
    assert!(
        !output.join("modules_discovered.tsv").exists(),
        "no module file should be written when the result is degenerate"
    );
    assert!(
        !output.join("module_discovery_report.tsv").exists(),
        "no report should be written when the result is degenerate"
    );

    // The refusal must leave its evidence behind: the beta sweep that
    // shows the criterion could not be met, with a sidecar saying so.
    let diag = output.join("soft_power_diagnostics.tsv");
    assert!(
        diag.exists(),
        "the soft-power sweep must survive the refusal"
    );
    let body = std::fs::read_to_string(&diag).unwrap();
    let n_rows = body.lines().count().saturating_sub(1);
    assert_eq!(n_rows, 20, "one sweep row per beta in 1..=20; got:\n{body}");
    assert!(
        stderr.contains("soft_power_diagnostics.tsv"),
        "the refusal must name the evidence file; got:\n{stderr}"
    );
    let sidecar = output.join("soft_power_diagnostics.tsv.run.json");
    assert!(sidecar.exists(), "the sweep needs a sidecar to be citable");
    let sc: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    assert_eq!(sc["outcome"], "refused");
    assert_eq!(sc["command"], "modules discover");
    assert_eq!(sc["args"]["scale-free-fit-achieved"], false);
    assert!(
        sc["output_files"]
            .as_object()
            .unwrap()
            .keys()
            .any(|k| k.ends_with("soft_power_diagnostics.tsv")),
        "the sidecar must hash the sweep it describes; got {}",
        sc["output_files"]
    );
    let written: Vec<String> = std::fs::read_dir(&output)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    let mut sorted = written.clone();
    sorted.sort();
    assert_eq!(
        sorted,
        vec![
            "soft_power_diagnostics.tsv".to_string(),
            "soft_power_diagnostics.tsv.run.json".to_string()
        ],
        "the refusal writes the evidence and nothing else"
    );
}

#[test]
fn module_de_refuses_a_grey_only_module_set() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    write_unstructured_fixture(&input, 24, 40);

    // A modules file of exactly the shape the old discover path emitted.
    let modules_tsv = tmp.path().join("modules_discovered.tsv");
    let mut body = String::from("module\tgene_symbol\n");
    for j in 1..=40 {
        body.push_str(&format!("grey\tG{j:04}\n"));
    }
    std::fs::write(&modules_tsv, body).unwrap();

    let out = run_atman(&[
        "module-de",
        "--input-dir",
        input.to_str().unwrap(),
        "--modules-tsv",
        modules_tsv.to_str().unwrap(),
        "--groups",
        "A-B",
        "--output-dir",
        tmp.path().join("mde_out").to_str().unwrap(),
    ]);
    assert!(
        !out.status.success(),
        "module-de must refuse a module set whose only member is grey"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("no module other than the grey catch-all"),
        "unexpected stderr:\n{stderr}"
    );
}
