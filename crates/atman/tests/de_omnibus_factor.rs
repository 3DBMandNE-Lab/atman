//! Integration test for `atman de --test ols --omnibus-factor`.
//!
//! Fixture has 24 samples split across 3 stages (CN, MCI, AD) with
//! 8 samples per stage. Two proteins are planted:
//!   `RESPONDER`  — activation depends on stage (CN=0, MCI=1, AD=2)
//!   `NONRESPONDER` — activation independent of stage
//!
//! Runs `atman de --test ols --design "~ condition + stage" --contrast
//! conditionCase --omnibus-factor stage` and asserts that:
//!   * `de_omnibus.tsv` contains a row per protein with the factor
//!     `stage`, a finite F-statistic, and p-values in (0, 1).
//!   * `RESPONDER` has omnibus `bh_q < 0.05`.
//!   * `NONRESPONDER` has omnibus `bh_q >= 0.05`.
//!   * Sidecar captures the `omnibus-factor` arg.
//!
//! Also a refusal test for `--omnibus-factor` without `--design`.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn parse_tsv(
    path: &Path,
) -> (Vec<String>, Vec<std::collections::HashMap<String, String>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines.next().unwrap().split('\t').map(String::from).collect();
    let rows = lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(h, c)| (h.clone(), c.to_string()))
                .collect()
        })
        .collect();
    (header, rows)
}

fn write_three_stage_canonical(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    let mut samples = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tstage\n",
    );
    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    proteins.push_str("olink_explore_ngs\tR001\tQ00001\tRESPONDER\tP1\t\n");
    proteins.push_str("olink_explore_ngs\tN001\tQ00002\tNONRESPONDER\tP1\t\n");
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let stages = ["CN", "MCI", "AD"];
    let mut order = 0u64;
    let mut qc_order = 0u64;
    // Deterministic LCG so noise looks random but fixed under the
    // test seed (no dependency on atman's Xoshiro).
    let mut state: u64 = 20260420;
    let mut noise = || -> f64 {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let u1 = (state >> 33) as f64 / (u32::MAX as f64);
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let u2 = (state >> 33) as f64 / (u32::MAX as f64);
        // Box-Muller standard normal.
        let u1 = u1.max(1e-12);
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    };
    for (stage_idx, stage) in stages.iter().enumerate() {
        for i in 0..8 {
            order += 1;
            let sid = format!("S{order:03}");
            // Alternate condition: half Case, half Control in each stage.
            let condition = if i % 2 == 0 { "Case" } else { "Control" };
            samples.push_str(&format!(
                "{sid}\t{sid}\t{condition}\t0\tplasma\t{order}\t{stage}\n"
            ));

            // RESPONDER: stage_idx * 1.0 + N(0, 0.5).
            let responder_val = stage_idx as f64 * 1.0 + 0.5 * noise();
            qc_order += 1;
            qc.push_str(&format!(
                "olink_explore_ngs\t{sid}\tR001\tRESPONDER\tP1\t{responder_val:.6}\t\
                 {responder_val:.6}\t{responder_val:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{qc_order}\n"
            ));

            // NONRESPONDER: stage-independent N(0, 1.0) noise.
            let non_val = noise();
            qc_order += 1;
            qc.push_str(&format!(
                "olink_explore_ngs\t{sid}\tN001\tNONRESPONDER\tP1\t{non_val:.6}\t\
                 {non_val:.6}\t{non_val:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{qc_order}\n"
            ));
        }
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    std::fs::write(dir.join("qc_measurements.tsv"), &qc).unwrap();
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

#[test]
fn omnibus_factor_identifies_responder_vs_nonresponder() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("de_omnibus_out");
    write_three_stage_canonical(&input);
    std::fs::create_dir_all(&output).unwrap();

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "ols",
        "--design",
        "~ condition + stage",
        "--contrast",
        "conditionCase",
        "--groups",
        "Case-Control",
        "--min-pairs",
        "4",
        "--omnibus-factor",
        "stage",
    ]);
    assert!(
        out.status.success(),
        "atman de failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );

    let omnibus_path = output.join("de_omnibus.tsv");
    assert!(omnibus_path.exists(), "de_omnibus.tsv missing");
    let (header, rows) = parse_tsv(&omnibus_path);
    for col in [
        "panel",
        "assay_id",
        "gene_symbol",
        "factor",
        "comparison",
        "f_statistic",
        "df_num",
        "df_den",
        "p_value",
        "bh_q",
    ] {
        assert!(header.iter().any(|h| h == col), "missing col {col}");
    }
    assert_eq!(rows.len(), 2, "expected 2 omnibus rows, got {}", rows.len());
    let by_gene: std::collections::HashMap<&str, &std::collections::HashMap<String, String>> = rows
        .iter()
        .map(|r| (r["gene_symbol"].as_str(), r))
        .collect();
    let responder = by_gene["RESPONDER"];
    let nonresponder = by_gene["NONRESPONDER"];
    let r_f: f64 = responder["f_statistic"].parse().unwrap();
    let r_p: f64 = responder["p_value"].parse().unwrap();
    let r_q: f64 = responder["bh_q"].parse().unwrap();
    let n_f: f64 = nonresponder["f_statistic"].parse().unwrap();
    let n_p: f64 = nonresponder["p_value"].parse().unwrap();
    let n_q: f64 = nonresponder["bh_q"].parse().unwrap();
    eprintln!(
        "RESPONDER   F={r_f:.3} p={r_p:.4} q={r_q:.4}\n\
         NONRESPONDER F={n_f:.3} p={n_p:.4} q={n_q:.4}"
    );
    assert!(r_f > n_f, "RESPONDER F should exceed NONRESPONDER F");
    assert!(r_q < 0.05, "RESPONDER bh_q should clear 0.05: {r_q}");
    assert!(n_q >= 0.05, "NONRESPONDER bh_q should fail 0.05: {n_q}");
    for r in &rows {
        assert_eq!(r["factor"], "stage");
        let df_num: usize = r["df_num"].parse().unwrap();
        assert_eq!(df_num, 2, "stage has 3 levels → df_num = 2");
    }

    // Sidecar records the flag.
    let sidecar = output.join("de_results.tsv.run.json");
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    assert_eq!(j["args"]["omnibus-factor"], "stage");
    let keys: Vec<&str> = j["output_files"]
        .as_object()
        .unwrap()
        .keys()
        .map(|s| s.as_str())
        .collect();
    assert!(keys.iter().any(|k| k.ends_with("de_omnibus.tsv")));
}

#[test]
fn omnibus_factor_requires_design_and_ols() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("omnibus_refuse");
    write_three_stage_canonical(&input);

    // Without --design, --omnibus-factor must refuse.
    let no_design = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "ols",
        "--groups",
        "Case-Control",
        "--omnibus-factor",
        "stage",
        "--min-pairs",
        "4",
    ]);
    assert!(!no_design.status.success());
    let stderr = String::from_utf8_lossy(&no_design.stderr);
    assert!(stderr.contains("--design"), "expected --design error: {stderr}");

    // With welch-t, --omnibus-factor must refuse.
    let no_ols = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "welch-t",
        "--groups",
        "Case-Control",
        "--omnibus-factor",
        "stage",
        "--min-pairs",
        "4",
    ]);
    assert!(!no_ols.status.success());
    let stderr = String::from_utf8_lossy(&no_ols.stderr);
    assert!(
        stderr.contains("--test ols"),
        "expected ols-only error: {stderr}"
    );
}
