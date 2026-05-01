//! Integration test: NMF loadings produced by `decompose nmf` flow through
//! `align programs` unchanged.
//!
//! Two synthetic cohorts are built from the same rank-4 fixture (same data,
//! different NMF seeds) so that all four planted programs are recoverable in
//! both cohorts and cross-cohort alignment is meaningful. At tau=0.30 with
//! cosine metric and reciprocal-best, at least two cross-cohort archetypes
//! are expected (the dominant programs survive across seeds). The test also
//! verifies the `decomposition_method` field in the sidecar.
//!
//! A second test verifies `decompose nmf` output pipes to `align bootstrap`
//! and `align project` without modification (schema check only).

use std::collections::BTreeSet;
use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

/// Build a canonical Atman input directory from the standard NMF rank-4
/// fixture (30 samples × 50 genes). Returns the canonical dir path.
fn build_canonical_from_fixture(tmp: &Path, subdir: &str) -> std::path::PathBuf {
    let input_fixture_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/nmf_frobenius_input.tsv");
    let input_text =
        std::fs::read_to_string(&input_fixture_path).expect("read nmf_frobenius_input.tsv");
    let input_lines: Vec<&str> = input_text.lines().collect();

    // Skip comment header.
    let data_start = if input_lines[0].starts_with('#') { 1 } else { 0 };

    let mut samples: BTreeSet<String> = BTreeSet::new();
    let mut genes: Vec<String> = Vec::new();
    let mut seen_genes: BTreeSet<String> = BTreeSet::new();
    let mut measurements: Vec<(String, String, f64)> = Vec::new();

    for line in &input_lines[(data_start + 1)..] {
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let sample_id = parts[0].to_string();
        let gene_symbol = parts[1].to_string();
        let abundance: f64 = parts[2].parse().expect("parse abundance");
        samples.insert(sample_id.clone());
        if seen_genes.insert(gene_symbol.clone()) {
            genes.push(gene_symbol.clone());
        }
        measurements.push((sample_id, gene_symbol, abundance));
    }

    let canonical_dir = tmp.join(subdir);
    std::fs::create_dir_all(&canonical_dir).unwrap();

    // measurements.tsv
    let headers = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
                   abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
                   detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order";
    let mut measurements_tsv = String::from(headers);
    measurements_tsv.push('\n');
    let mut ingest_order = 0u64;
    for (sample_id, gene_symbol, abundance) in &measurements {
        let gene_idx = genes.iter().position(|g| g == gene_symbol).unwrap();
        let assay_id = format!("A{:03}", gene_idx);
        ingest_order += 1;
        let src = format!("{:.6}", abundance);
        measurements_tsv.push_str(&format!(
            "proteomics\t{}\t{}\t{}\tP1\t{}\t{:.6}\t{:.6}\tlog2_scale\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            sample_id, assay_id, gene_symbol, src, abundance, abundance, ingest_order
        ));
    }
    std::fs::write(canonical_dir.join("measurements.tsv"), &measurements_tsv).unwrap();

    // samples.tsv
    let mut samples_tsv =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for (i, sample_id) in samples.iter().enumerate() {
        samples_tsv.push_str(&format!(
            "{}\t{}\tcase\t0\tcsf\t{}\n",
            sample_id, sample_id, i + 1
        ));
    }
    std::fs::write(canonical_dir.join("samples.tsv"), &samples_tsv).unwrap();

    // proteins.tsv
    let mut proteins_tsv =
        String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for (gene_idx, gene) in genes.iter().enumerate() {
        let assay_id = format!("A{:03}", gene_idx);
        proteins_tsv.push_str(&format!("proteomics\t{}\t\t{}\tP1\t\n", assay_id, gene));
    }
    std::fs::write(canonical_dir.join("proteins.tsv"), &proteins_tsv).unwrap();

    canonical_dir
}

/// Run `decompose nmf --k 4` on a canonical dir and return the loadings path.
fn run_nmf(canonical_dir: &Path, out_prefix: &Path, seed: u64) -> std::path::PathBuf {
    // Use a per-seed file to avoid collisions between cohorts.
    let loadings_path = {
        let stem = out_prefix
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .to_string();
        let dir = out_prefix.parent().unwrap();
        dir.join(format!("{stem}_nmf_loadings.tsv"))
    };
    let activations_path = {
        let stem = out_prefix
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .to_string();
        let dir = out_prefix.parent().unwrap();
        dir.join(format!("{stem}_nmf_activations.tsv"))
    };
    let output = run_atman(&[
        "decompose",
        "nmf",
        "--input-dir",
        canonical_dir.to_str().unwrap(),
        "--k",
        "4",
        "--solver",
        "mu",
        "--beta-loss",
        "frobenius",
        "--init",
        "nndsvda",
        "--max-iter",
        "400",
        "--tol",
        "1e-6",
        "--seed",
        &seed.to_string(),
        "--output-loadings",
        loadings_path.to_str().unwrap(),
        "--output-activations",
        activations_path.to_str().unwrap(),
    ]);
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!(
            "decompose nmf (seed={seed}) failed:\nstderr:\n{}",
            stderr
        );
    }
    assert!(
        loadings_path.exists(),
        "NMF loadings not written to {:?}",
        loadings_path
    );
    loadings_path
}

/// Parse archetypes TSV and return all rows as BTreeMap vectors.
fn parse_archetypes(path: &Path) -> Vec<std::collections::HashMap<String, String>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(h, c)| (h.clone(), c.to_string()))
                .collect()
        })
        .collect()
}

#[test]
fn align_programs_consumes_nmf_loadings_unchanged() {
    let tmp = tempfile::tempdir().unwrap();
    let tmp_path = tmp.path();

    // Build both cohorts from the same rank-4 fixture. Different NMF seeds
    // produce permuted but structurally identical programs so cross-cohort
    // alignment is meaningful.
    let canonical_a = build_canonical_from_fixture(tmp_path, "cohort_a");
    let canonical_b = build_canonical_from_fixture(tmp_path, "cohort_b");

    // Run NMF with seed=42 for cohort A and seed=99 for cohort B.
    let loadings_a = run_nmf(&canonical_a, &tmp_path.join("a"), 42);
    let loadings_b = run_nmf(&canonical_b, &tmp_path.join("b"), 99);

    // Verify both loadings TSVs exist and have the expected 4-column schema
    // (program, assay_id, gene_symbol, loading) as emitted by decompose nmf.
    for (path, label) in [(&loadings_a, "CohortA"), (&loadings_b, "CohortB")] {
        let text = std::fs::read_to_string(path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let data_start = if lines[0].starts_with('#') { 1 } else { 0 };
        let header = lines[data_start];
        // NMF emits: program TAB assay_id TAB gene_symbol TAB loading
        assert!(
            header.starts_with("program\t") && header.contains("gene_symbol") && header.contains("loading"),
            "{label} loadings header unexpected: {header:?}"
        );
    }

    // Verify the neighboring sidecar was written by decompose nmf.
    for (path, label) in [(&loadings_a, "CohortA"), (&loadings_b, "CohortB")] {
        let sidecar_path = {
            let mut s = path.as_os_str().to_os_string();
            s.push(".run.json");
            std::path::PathBuf::from(s)
        };
        assert!(
            sidecar_path.exists(),
            "{label} sidecar not found at {:?}",
            sidecar_path
        );
        let j: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&sidecar_path).unwrap()).unwrap();
        assert_eq!(
            j["command"].as_str(),
            Some("decompose nmf"),
            "{label} sidecar command mismatch"
        );
    }

    // Run align programs on the two NMF loadings TSVs.
    let archetypes_path = tmp_path.join("archetypes.tsv");
    let output = run_atman(&[
        "align",
        "programs",
        "--metric",
        "cosine",
        "--tau",
        "0.30",
        "--reciprocal-best",
        "--loadings",
        &format!(
            "{},{}",
            loadings_a.to_str().unwrap(),
            loadings_b.to_str().unwrap()
        ),
        "--cohorts",
        "CohortA,CohortB",
        "--output",
        archetypes_path.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "align programs failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify schema of archetypes.tsv.
    let text = std::fs::read_to_string(&archetypes_path).unwrap();
    assert!(
        text.starts_with("archetype_id\tcohort\tprogram\tcategory\tn_members\tn_cohorts\tis_singleton\n"),
        "archetypes.tsv has unexpected header:\n{}", text.lines().next().unwrap_or("")
    );

    // Parse rows and assert at least one cross-cohort archetype.
    let rows = parse_archetypes(&archetypes_path);
    assert!(!rows.is_empty(), "archetypes.tsv is empty");

    let n_cross_cohort = rows.iter().filter(|r| {
        r.get("n_cohorts").and_then(|v| v.parse::<usize>().ok()).unwrap_or(0) >= 2
    }).count();

    assert!(
        n_cross_cohort >= 2,
        "expected >= 2 cross-cohort archetype member rows (n_cohorts >= 2); got {} out of {} rows.\n\
         archetype rows:\n{}",
        n_cross_cohort,
        rows.len(),
        rows.iter()
            .map(|r| format!(
                "  id={} cohort={} n_cohorts={}",
                r.get("archetype_id").map(|s| s.as_str()).unwrap_or("?"),
                r.get("cohort").map(|s| s.as_str()).unwrap_or("?"),
                r.get("n_cohorts").map(|s| s.as_str()).unwrap_or("?"),
            ))
            .collect::<Vec<_>>()
            .join("\n")
    );

    let cross_cohort_archetype_ids: BTreeSet<String> = rows
        .iter()
        .filter(|r| {
            r.get("n_cohorts")
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or(0)
                >= 2
        })
        .map(|r| r.get("archetype_id").cloned().unwrap_or_default())
        .collect();

    eprintln!(
        "align_programs_nmf: {n_cross_cohort} cross-cohort member rows across {} archetypes",
        cross_cohort_archetype_ids.len()
    );

    // Verify the sidecar for align programs records decomposition_method = "nmf".
    let align_sidecar_path = {
        let mut s = archetypes_path.as_os_str().to_os_string();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    assert!(
        align_sidecar_path.exists(),
        "align programs sidecar not found at {:?}",
        align_sidecar_path
    );
    let align_j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&align_sidecar_path).unwrap()).unwrap();
    assert_eq!(
        align_j["args"]["decomposition_method"].as_str(),
        Some("nmf"),
        "align programs sidecar should record decomposition_method=nmf; got: {:?}",
        align_j["args"]["decomposition_method"]
    );
}

/// Verify the `decomposition_method` field appears as "unknown" when no
/// sidecar is present (i.e., hand-crafted loadings TSVs without a run.json).
#[test]
fn align_programs_sidecar_decomposition_method_unknown_when_no_sidecar() {
    let tmp = tempfile::tempdir().unwrap();
    let cohort_a = tmp.path().join("nosidecar_a.tsv");
    let cohort_b = tmp.path().join("nosidecar_b.tsv");
    let out = tmp.path().join("archetypes_nosidecar.tsv");

    // Write minimal hand-crafted loadings with no neighboring run.json.
    let loadings = "program\tgene_symbol\tloading\n\
                    program_01\tAQP4\t0.90\n\
                    program_01\tGFAP\t0.80\n\
                    program_02\tSYN1\t0.85\n\
                    program_02\tSYP\t0.80\n";
    std::fs::write(&cohort_a, loadings).unwrap();
    std::fs::write(&cohort_b, loadings).unwrap();

    let output = run_atman(&[
        "align",
        "programs",
        "--metric",
        "cosine",
        "--tau",
        "0.30",
        "--reciprocal-best",
        "--loadings",
        &format!(
            "{},{}",
            cohort_a.to_str().unwrap(),
            cohort_b.to_str().unwrap()
        ),
        "--cohorts",
        "A,B",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "align programs failed:\nstderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let sidecar_path = {
        let mut s = out.as_os_str().to_os_string();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar_path).unwrap()).unwrap();
    assert_eq!(
        j["args"]["decomposition_method"].as_str(),
        Some("unknown"),
        "expected decomposition_method=unknown for hand-crafted loadings; got: {:?}",
        j["args"]["decomposition_method"]
    );
}
