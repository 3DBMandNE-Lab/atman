//! Priority 8 integration test: `atman modules discover`.
//!
//! Two-block synthetic fixture: 40 subjects × (20 + 20 + 40) proteins,
//! where proteins 1-20 share a latent Gaussian (block 1), proteins
//! 21-40 share a different latent (block 2), and 41-80 are
//! independent noise. Runs `atman modules discover --method
//! wgcna-soft --soft-power auto` and verifies:
//!
//! 1. `modules_discovered.tsv` schema is (`module`, `gene_symbol`).
//! 2. At least two non-grey modules of size ≥ 15 are recovered.
//! 3. The soft-power diagnostics TSV is emitted with exactly
//!    `--max-beta` rows when auto-selecting.
//! 4. Refuses cleanly when the retained feature count falls below
//!    `--min-module-size`.

use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn parse_tsv(path: &Path) -> (Vec<String>, Vec<HashMap<String, String>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
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

fn box_muller(u1: f64, u2: f64) -> f64 {
    let u1 = u1.max(1e-12);
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

/// Write a canonical Atman directory with 40 subjects × 80 proteins
/// containing two planted correlation blocks + noise background.
fn write_two_block_cohort(dir: &Path, seed: u64) {
    std::fs::create_dir_all(dir).unwrap();
    let n = 40usize;
    let block = 20usize;
    let noise = 40usize;
    let p = 2 * block + noise;
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=n {
        samples.push_str(&format!("S{i:03}\tS{i:03}\tN/A\t0\tplasma\t{i}\n"));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 1..=p {
        proteins.push_str(&format!(
            "olink_explore_ngs\tA{j:03}\tQ{j:05}\tG{j:03}\tP1\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    // Deterministic LCG → Box-Muller Gaussians.
    let mut state = std::num::Wrapping(seed);
    let mut next_u = || {
        state = state * std::num::Wrapping(6364136223846793005_u64)
            + std::num::Wrapping(1442695040888963407_u64);
        ((state.0 >> 33) as f64) / (u32::MAX as f64)
    };
    let mut next_z = || box_muller(next_u(), next_u());

    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order: u64 = 0;
    for i in 1..=n {
        let z1 = next_z();
        let z2 = next_z();
        for j in 1..=p {
            order += 1;
            let val = if j <= block {
                2.0 * z1 + 0.3 * next_z()
            } else if j <= 2 * block {
                2.0 * z2 + 0.3 * next_z()
            } else {
                next_z()
            };
            qc.push_str(&format!(
                "olink_explore_ngs\tS{i:03}\tA{j:03}\tG{j:03}\tP1\t{val:.6}\t\
                 {val:.6}\t{val:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

#[test]
fn modules_discover_recovers_two_planted_blocks() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("cohort");
    let output = tmp.path().join("modules_out");
    write_two_block_cohort(&input, 20260420);

    let status = run_atman(&[
        "modules",
        "discover",
        "--input-dir",
        input.to_str().unwrap(),
        "--method",
        "wgcna-soft",
        "--soft-power",
        "auto",
        "--r2-target",
        "0.8",
        "--max-beta",
        "20",
        "--n-bins",
        "10",
        "--cut-height",
        "0.5",
        "--min-module-size",
        "10",
        "--output-dir",
        output.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "modules discover failed:\nstderr:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );

    let modules_path = output.join("modules_discovered.tsv");
    let (header, rows) = parse_tsv(&modules_path);
    assert_eq!(
        header,
        vec!["module".to_string(), "gene_symbol".to_string()]
    );
    assert!(!rows.is_empty());

    // Count members per module (ignoring grey).
    let mut per_module: HashMap<String, usize> = HashMap::new();
    for r in &rows {
        *per_module.entry(r["module"].clone()).or_insert(0) += 1;
    }
    let non_grey_large: Vec<(&String, &usize)> = per_module
        .iter()
        .filter(|(m, n)| *m != "grey" && **n >= 15)
        .collect();
    assert!(
        non_grey_large.len() >= 2,
        "expected ≥ 2 non-grey modules with size ≥ 15; got {:?}",
        per_module
    );

    // Verify block separation: block 1 genes (G001..G020) should
    // concentrate into one module; block 2 genes (G021..G040) into
    // another. Count the modal module per block.
    let mut by_gene: HashMap<String, String> = HashMap::new();
    for r in &rows {
        by_gene.insert(r["gene_symbol"].clone(), r["module"].clone());
    }
    let mut block1_mods: HashMap<String, usize> = HashMap::new();
    let mut block2_mods: HashMap<String, usize> = HashMap::new();
    for j in 1..=20 {
        let gene = format!("G{j:03}");
        *block1_mods.entry(by_gene[&gene].clone()).or_insert(0) += 1;
    }
    for j in 21..=40 {
        let gene = format!("G{j:03}");
        *block2_mods.entry(by_gene[&gene].clone()).or_insert(0) += 1;
    }
    let block1_mode = block1_mods.iter().max_by_key(|(_, n)| *n).unwrap();
    let block2_mode = block2_mods.iter().max_by_key(|(_, n)| *n).unwrap();
    assert!(
        *block1_mode.1 >= 15,
        "block 1 should concentrate into one module (≥15 members); got {:?}",
        block1_mods
    );
    assert!(
        *block2_mode.1 >= 15,
        "block 2 should concentrate into one module (≥15 members); got {:?}",
        block2_mods
    );
    // And the two block modes should be DIFFERENT modules.
    assert_ne!(
        block1_mode.0, block2_mode.0,
        "blocks should land in different modules; both got {}",
        block1_mode.0
    );

    // Soft-power diagnostics TSV should have --max-beta rows.
    let diag_path = output.join("soft_power_diagnostics.tsv");
    assert!(diag_path.exists());
    let (_, diag_rows) = parse_tsv(&diag_path);
    assert_eq!(diag_rows.len(), 20);

    // Per-module report exists and records hub + eigenprotein variance.
    let report_path = output.join("module_discovery_report.tsv");
    let (report_header, report_rows) = parse_tsv(&report_path);
    for col in [
        "module",
        "size",
        "mean_within_abs_correlation",
        "hub_feature",
        "eigenprotein_pc1_variance_explained",
    ] {
        assert!(
            report_header.iter().any(|h| h == col),
            "missing {col} column"
        );
    }
    assert!(!report_rows.is_empty());

    // Sidecar exists.
    let sidecar = output.join("modules_discovered.tsv.run.json");
    assert!(sidecar.exists());
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    assert!(j["args"]["chosen-beta"].as_u64().unwrap() >= 1);

    // Output is consumable by `score modules` / module-de: the set of
    // unique modules in the output file is small and each gene_symbol
    // appears exactly once.
    let genes_in_output: HashSet<String> = rows.iter().map(|r| r["gene_symbol"].clone()).collect();
    assert_eq!(
        genes_in_output.len(),
        rows.len(),
        "each gene must appear once"
    );
}

#[test]
fn modules_discover_refuses_when_features_below_min_module_size() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("small");
    // 2 proteins, 5 samples — below any reasonable min_module_size.
    std::fs::create_dir_all(&input).unwrap();
    std::fs::write(
        input.join("samples.tsv"),
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
         S001\tS001\tN/A\t0\tplasma\t1\n\
         S002\tS002\tN/A\t0\tplasma\t2\n\
         S003\tS003\tN/A\t0\tplasma\t3\n\
         S004\tS004\tN/A\t0\tplasma\t4\n\
         S005\tS005\tN/A\t0\tplasma\t5\n",
    )
    .unwrap();
    std::fs::write(
        input.join("proteins.tsv"),
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
         olink_explore_ngs\tA001\tQ00001\tG1\tP1\t\n\
         olink_explore_ngs\tA002\tQ00002\tG2\tP1\t\n",
    )
    .unwrap();
    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0;
    for i in 1..=5 {
        for (j, g) in ["G1", "G2"].iter().enumerate() {
            order += 1;
            qc.push_str(&format!(
                "olink_explore_ngs\tS{i:03}\tA{:03}\t{g}\tP1\t1.0\t1.0\t1.0\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n",
                j + 1
            ));
        }
    }
    std::fs::write(input.join("measurements.tsv"), &qc).unwrap();
    std::fs::write(input.join("measurements.tsv"), &qc).unwrap();

    let out = run_atman(&[
        "modules",
        "discover",
        "--input-dir",
        input.to_str().unwrap(),
        "--method",
        "wgcna-soft",
        "--min-module-size",
        "5",
        "--output-dir",
        tmp.path().join("out").to_str().unwrap(),
    ]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("min-module-size"),
        "expected min-module-size error, got: {stderr}"
    );
}
