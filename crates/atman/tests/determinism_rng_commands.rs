//! TASK-022 — byte-identity determinism across all RNG-bearing / iterative commands.
//!
//! Post-RNG-unification (TASK-017, `atman_core::rng`), every command that draws
//! random numbers or runs an iterative solver must produce byte-identical output
//! across two runs with the same `--seed`. This file extends the coverage in
//! `determinism_new_methods.rs` (which handles nmf / missingness-ica / --adjust-for)
//! to:
//!
//!   - `null`            (permutation null)
//!   - `bootstrap protein` (subject-level bootstrap)
//!   - `ratio`           (bootstrap CI on module-score log-ratios)
//!   - `align bootstrap` (subject-level archetype-alignment bootstrap)
//!   - `de --test ensemble` (multi-method ensemble DE)
//!   - `decompose ica`   (standard FastICA, multi-seed stability)
//!   - `robust-paired`   (leave-one-pair-out diagnostic; deterministic by
//!                        construction — no seed — but covered for byte-identity)
//!
//! Each test runs the command TWICE into two separate output dirs and asserts the
//! produced output files are byte-identical. The `.run.json` sidecar carries
//! non-deterministic fields (run_uuid, started_at, finished_at); rather than diff
//! the whole sidecar, we diff its SHA fields (inputs_sha256 + output_files),
//! mirroring the exclusion approach in `determinism_new_methods.rs`. Identical
//! output-file SHAs is the byte-identity claim for the data products themselves;
//! we additionally byte-compare the TSVs directly.
//!
//! BLAS/OMP threads are pinned to 1 (as in `determinism_new_methods.rs`) so any
//! BLAS-backed path (ICA, ensemble/limma) stays deterministic.

use std::path::Path;
use std::process::Command;

fn run_atman(args: &[&str]) -> std::process::Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    let mut cmd = Command::new(bin);
    // Pin BLAS threading to 1 to ensure determinism.
    cmd.env("BLAS_NUM_THREADS", "1")
        .env("OMP_NUM_THREADS", "1")
        .env("OPENBLAS_NUM_THREADS", "1");
    cmd.args(args).output().expect("run atman")
}

/// Compare two files byte-for-byte.
fn assert_byte_identical(path_a: &Path, path_b: &Path, label: &str) {
    let bytes_a = std::fs::read(path_a).unwrap_or_else(|e| panic!("read {}: {e}", path_a.display()));
    let bytes_b = std::fs::read(path_b).unwrap_or_else(|e| panic!("read {}: {e}", path_b.display()));
    assert_eq!(
        bytes_a,
        bytes_b,
        "{label} byte mismatch:\n  A: {}\n  B: {}",
        path_a.display(),
        path_b.display()
    );
}

/// Parse a `.run.json` sidecar and compare only its deterministic SHA fields
/// (`inputs_sha256` and `output_files`), excluding non-deterministic fields
/// (run_uuid, started_at, finished_at) and path-bearing keys. This mirrors the
/// exclusion approach used by `determinism_new_methods.rs`.
fn assert_sidecar_shas_match(path_a: &Path, path_b: &Path, label: &str) {
    let read_json = |p: &Path| -> serde_json::Value {
        let text =
            std::fs::read_to_string(p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()));
        serde_json::from_str(&text).unwrap_or_else(|e| panic!("parse {}: {e}", p.display()))
    };
    let ja = read_json(path_a);
    let jb = read_json(path_b);

    // inputs_sha256: per-file dict {path: "sha256:..."}; paths are stable across
    // runs that share an input dir, but for runs with distinct input dirs we
    // compare only the sorted set of SHA values.
    let sha_values = |v: &serde_json::Value| -> Vec<String> {
        let mut out: Vec<String> = v
            .as_object()
            .map(|m| m.values().filter_map(|x| x.as_str().map(String::from)).collect())
            .unwrap_or_default();
        out.sort();
        out
    };
    assert_eq!(
        sha_values(&ja["inputs_sha256"]),
        sha_values(&jb["inputs_sha256"]),
        "{label} sidecar inputs_sha256 mismatch"
    );

    // output_files: {abs_path: "sha256:..."}. Paths differ (distinct output
    // dirs), so compare (filename, sha) pairs.
    let filename_shas = |v: &serde_json::Value| -> Vec<(String, String)> {
        let mut out = Vec::new();
        if let Some(m) = v.as_object() {
            for (full_path, sha) in m {
                let filename = std::path::Path::new(full_path)
                    .file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_else(|| full_path.clone());
                out.push((filename, sha.as_str().unwrap_or("").to_string()));
            }
        }
        out.sort();
        out
    };
    let out_a = filename_shas(&ja["output_files"]);
    // Guard against a trivial pass: an empty (or absent) output_files map would
    // make the equality below vacuously true. Every RNG command emits at least
    // one output file, so assert the set is non-empty.
    assert!(
        !out_a.is_empty(),
        "{label} sidecar output_files is empty — comparison would be vacuous"
    );
    assert_eq!(
        out_a,
        filename_shas(&jb["output_files"]),
        "{label} sidecar output_files SHA mismatch"
    );
}

/// The derived `<stem>_summary.<ext>` companion path a command writes next to
/// its primary `--output` file (used by robust-paired).
fn summary_sibling(out: &Path) -> std::path::PathBuf {
    let stem = out.file_stem().unwrap().to_string_lossy();
    let ext = out.extension().map(|e| e.to_string_lossy().to_string()).unwrap_or_default();
    out.with_file_name(format!("{stem}_summary.{ext}"))
}

fn sidecar_for(output: &Path) -> std::path::PathBuf {
    let mut s = output.as_os_str().to_os_string();
    s.push(".run.json");
    std::path::PathBuf::from(s)
}

// ---------------------------------------------------------------------------
// Fixture builders (mirrored from the commands' own integration tests).
// ---------------------------------------------------------------------------

/// Two-group canonical fixture (mirrors `tests/null.rs` / `tests/bootstrap_protein.rs`).
fn write_two_group_fixture(dir: &Path, two_proteins: bool) {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
         C1\tC1\tControl\t0\tcsf\t1\n\
         C2\tC2\tControl\t0\tcsf\t2\n\
         C3\tC3\tControl\t0\tcsf\t3\n\
         K1\tK1\tCase\t0\tcsf\t4\n\
         K2\tK2\tCase\t0\tcsf\t5\n\
         K3\tK3\tCase\t0\tcsf\t6\n",
    )
    .unwrap();
    let proteins = if two_proteins {
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
         spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n\
         spectronaut_report\tP00002\tP00002\tGENE2\tspectronaut\t\n"
    } else {
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
         spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n"
    };
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    let mut rows = vec![
        ("C1", "P00001", "GENE1", "10.0", 1),
        ("C2", "P00001", "GENE1", "11.0", 2),
        ("C3", "P00001", "GENE1", "12.0", 3),
        ("K1", "P00001", "GENE1", "13.0", 4),
        ("K2", "P00001", "GENE1", "14.5", 5),
        ("K3", "P00001", "GENE1", "16.0", 6),
    ];
    if two_proteins {
        rows.extend([
            ("C1", "P00002", "GENE2", "8.0", 7),
            ("C2", "P00002", "GENE2", "9.0", 8),
            ("C3", "P00002", "GENE2", "10.0", 9),
            ("K1", "P00002", "GENE2", "8.5", 10),
            ("K2", "P00002", "GENE2", "9.5", 11),
            ("K3", "P00002", "GENE2", "10.5", 12),
        ]);
    }
    let mut measurements = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    for (sample, assay, gene, value, order) in rows {
        measurements.push_str(&format!(
            "spectronaut_report\t{sample}\t{assay}\t{gene}\tspectronaut\t{value}\t{value}\t{value}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
        ));
    }
    std::fs::write(dir.join("measurements.tsv"), measurements).unwrap();
}

/// Module-score fixture for `ratio` (mirrors `tests/ratio.rs`).
fn write_ratio_fixture(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
         C1\tC1\tControl\t1\tbio\t1\n\
         C2\tC2\tControl\t1\tbio\t2\n\
         C3\tC3\tControl\t1\tbio\t3\n\
         C4\tC4\tControl\t1\tbio\t4\n\
         M1\tM1\tMS\t0\tbio\t5\n\
         M2\tM2\tMS\t0\tbio\t6\n\
         M3\tM3\tMS\t0\tbio\t7\n\
         M4\tM4\tMS\t0\tbio\t8\n",
    )
    .unwrap();
    let mut meas = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let rows: &[(&str, f64, f64)] = &[
        ("C1", 1.0, 0.8),
        ("C2", 1.1, 0.8),
        ("C3", 1.2, 0.8),
        ("C4", 1.3, 0.8),
        ("M1", 2.0, 1.0),
        ("M2", 2.1, 1.0),
        ("M3", 2.2, 1.0),
        ("M4", 2.3, 1.0),
    ];
    let mut order = 1;
    for (sample, num, den) in rows {
        meas.push_str(&format!(
            "somascan\t{sample}\tintrathecal_IgV\tintrathecal_IgV\tmodule_score\t{num}\t{num}\t{num}\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
        ));
        order += 1;
        meas.push_str(&format!(
            "somascan\t{sample}\tplasma_IgG\tplasma_IgG\tmodule_score\t{den}\t{den}\t{den}\tmodule_score\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
        ));
        order += 1;
    }
    std::fs::write(dir.join("measurements.tsv"), meas).unwrap();
}

/// Two-condition fixture for `de --test ensemble` (mirrors `tests/de_ensemble.rs`).
fn write_de_ensemble_fixture(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=12 {
        let cond = if i <= 6 { "A" } else { "B" };
        samples.push_str(&format!("S{i:02}\tS{i:02}\t{cond}\t0\tplasma\t{i}\n"));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    std::fs::write(
        dir.join("proteins.tsv"),
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
         olink_explore_ngs\tA001\tQ01\tUP\tP1\t\n\
         olink_explore_ngs\tA002\tQ02\tDN\tP1\t\n\
         olink_explore_ngs\tA003\tQ03\tST\tP1\t\n",
    )
    .unwrap();
    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0;
    for i in 1..=12 {
        let is_b = i > 6;
        for (assay, gene, effect) in [("A001", "UP", 1.5_f64), ("A002", "DN", -1.5), ("A003", "ST", 0.0)]
        {
            order += 1;
            let v = 10.0 + if is_b { effect } else { 0.0 } + ((i * 7 + order) as f64).sin() * 0.05;
            qc.push_str(&format!(
                "olink_explore_ngs\tS{i:02}\t{assay}\t{gene}\tP1\t{v:.6}\t\
                 {v:.6}\t{v:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

/// Synthetic linear-mix fixture for `decompose ica` (mirrors `tests/decompose_ica.rs`).
fn write_ica_fixture(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    let n_samples = 30;
    let n_assays = 40;
    let mut xs = Vec::with_capacity(n_samples);
    let mut state: u64 = 1234567;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) as f64) / (u32::MAX as f64)
    };
    for _ in 0..n_samples {
        let s1 = (next() - 0.5) * 2.0;
        let s2 = (next() - 0.5) * 2.0;
        let mut row = Vec::with_capacity(n_assays);
        for j in 0..n_assays {
            let a1 = if j < n_assays / 2 { 1.0 } else { 0.1 };
            let a2 = if j >= n_assays / 2 { 1.0 } else { 0.1 };
            let noise = (next() - 0.5) * 0.05;
            row.push(a1 * s1 + a2 * s2 + noise);
        }
        xs.push(row);
    }
    let mut measurements = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut ingest = 0u64;
    for (i, row) in xs.iter().enumerate() {
        let sample_id = format!("S{i:02}");
        for (j, v) in row.iter().enumerate() {
            let assay_id = format!("A{j:03}");
            let gene = format!("GENE{j:03}");
            ingest += 1;
            let src = format!("{v:.6}");
            measurements.push_str(&format!(
                "olink_explore_ngs\t{sample_id}\t{assay_id}\t{gene}\tP1\t{src}\t{v:.6}\t{v:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{ingest}\n",
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &measurements).unwrap();
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 0..n_samples {
        samples.push_str(&format!("S{i:02}\tS{i:02}\tcase\t0\tcsf\t{}\n", i + 1));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 0..n_assays {
        proteins.push_str(&format!("olink_explore_ngs\tA{j:03}\t\tGENE{j:03}\tP1\t\n"));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();
}

/// Paired fixture for `robust-paired` (mirrors `tests/robust_paired.rs`).
fn write_robust_paired_fixture(dir: &Path, n_pairs: usize) {
    std::fs::create_dir_all(dir).unwrap();
    let mut samples_buf = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tpatient_id\n",
    );
    let mut order = 1u64;
    for k in 0..n_pairs {
        let pid = format!("pair_{:03}", k);
        samples_buf.push_str(&format!("T{:03}\t{}\ttumor\t0\ttumor\t{}\t{}\n", k, k, order, pid));
        order += 1;
        samples_buf.push_str(&format!(
            "P{:03}\t{}\tpaired_non_tumor\t0\tnormal\t{}\t{}\n",
            k, k, order, pid
        ));
        order += 1;
    }
    std::fs::write(dir.join("samples.tsv"), samples_buf).unwrap();
    std::fs::write(
        dir.join("proteins.tsv"),
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
         diann_report\tA_STRONG\t\tSTRONG\tms\t\n\
         diann_report\tA_FRAGILE\t\tFRAGILE\tms\t\n",
    )
    .unwrap();
    let mut buf = String::from("platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n");
    let mut row = 1u64;
    for k in 0..n_pairs {
        let jitter = ((k as f64) * 0.013).sin() * 0.05;
        let t_val = 11.0 + jitter;
        let n_val = 10.0 - jitter;
        buf.push_str(&format!(
            "diann_report\tT{:03}\tA_STRONG\tSTRONG\tms\t\t{}\t{}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            k, t_val, t_val, row
        ));
        row += 1;
        buf.push_str(&format!(
            "diann_report\tP{:03}\tA_STRONG\tSTRONG\tms\t\t{}\t{}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            k, n_val, n_val, row
        ));
        row += 1;
        let f_jitter = ((k as f64) * 0.27).cos() * 0.05;
        let (tf, nf) = if k == 0 { (30.0, 10.0) } else { (10.0 + f_jitter, 10.0 - f_jitter) };
        buf.push_str(&format!(
            "diann_report\tT{:03}\tA_FRAGILE\tFRAGILE\tms\t\t{}\t{}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            k, tf, tf, row
        ));
        row += 1;
        buf.push_str(&format!(
            "diann_report\tP{:03}\tA_FRAGILE\tFRAGILE\tms\t\t{}\t{}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            k, nf, nf, row
        ));
        row += 1;
    }
    std::fs::write(dir.join("measurements.tsv"), buf).unwrap();
}

/// Two-cohort planted-archetype fixture for `align bootstrap`
/// (mirrors `tests/align_bootstrap.rs`).
struct Lcg(u64);
impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as f64 / (u32::MAX as f64)
    }
    fn heavy_tail(&mut self) -> f64 {
        let u = self.next() * 2.0 - 1.0;
        u * u * u
    }
}

fn write_align_cohort(
    dir: &Path,
    cohort_label: &str,
    planted_a_specific: bool,
    planted_b_specific: bool,
    seed: u64,
) {
    std::fs::create_dir_all(dir).unwrap();
    let n_samples = 30usize;
    let n_proteins = 10usize;
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=n_samples {
        samples.push_str(&format!(
            "{cohort_label}_S{i:03}\t{cohort_label}_S{i:03}\tN/A\t0\tplasma\t{i}\n"
        ));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 1..=n_proteins {
        proteins.push_str(&format!("olink_explore_ngs\tA{j:03}\tQ{j:05}\tG{j:03}\tP1\t\n"));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();
    let mut rng = Lcg::new(seed);
    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order: u64 = 0;
    for i in 1..=n_samples {
        let universal_source = rng.heavy_tail() * 2.0;
        let a_source = rng.heavy_tail() * 2.0;
        let b_source = rng.heavy_tail() * 2.0;
        for j in 1..=n_proteins {
            order += 1;
            let universal_weight = if (1..=4).contains(&j) { 1.0 } else { 0.0 };
            let a_weight = if planted_a_specific && (5..=7).contains(&j) { 1.0 } else { 0.0 };
            let b_weight = if planted_b_specific && (8..=10).contains(&j) { 1.0 } else { 0.0 };
            let noise = (rng.next() - 0.5) * 0.3;
            let value =
                universal_weight * universal_source + a_weight * a_source + b_weight * b_source + noise;
            qc.push_str(&format!(
                "olink_explore_ngs\t{cohort_label}_S{i:03}\tA{j:03}\tG{j:03}\tP1\t{value:.6}\t\
                 {value:.6}\t{value:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

#[test]
fn null_permutation_is_byte_deterministic() {
    let tmp1 = tempfile::tempdir().unwrap();
    let tmp2 = tempfile::tempdir().unwrap();
    let input1 = tmp1.path().join("input");
    let input2 = tmp2.path().join("input");
    let out1 = tmp1.path().join("out");
    let out2 = tmp2.path().join("out");
    write_two_group_fixture(&input1, true);
    write_two_group_fixture(&input2, true);

    for (input, out) in [(&input1, &out1), (&input2, &out2)] {
        let r = run_atman(&[
            "null",
            "--input-dir",
            input.to_str().unwrap(),
            "--output-dir",
            out.to_str().unwrap(),
            "--groups",
            "Case-Control",
            "--test",
            "welch-t",
            "--n",
            "100",
            "--seed",
            "99",
            "--min-pairs",
            "2",
        ]);
        assert!(r.status.success(), "null failed:\n{}", String::from_utf8_lossy(&r.stderr));
    }

    assert_byte_identical(
        &out1.join("null_summary.tsv"),
        &out2.join("null_summary.tsv"),
        "null_summary",
    );
    assert_byte_identical(
        &out1.join("empirical_p.tsv"),
        &out2.join("empirical_p.tsv"),
        "empirical_p",
    );
    assert_sidecar_shas_match(
        &sidecar_for(&out1.join("null_summary.tsv")),
        &sidecar_for(&out2.join("null_summary.tsv")),
        "null sidecar",
    );
}

#[test]
fn bootstrap_protein_is_byte_deterministic() {
    let tmp1 = tempfile::tempdir().unwrap();
    let tmp2 = tempfile::tempdir().unwrap();
    let input1 = tmp1.path().join("input");
    let input2 = tmp2.path().join("input");
    let out1 = tmp1.path().join("boot.tsv");
    let out2 = tmp2.path().join("boot.tsv");
    write_two_group_fixture(&input1, true);
    write_two_group_fixture(&input2, true);

    for (input, out) in [(&input1, &out1), (&input2, &out2)] {
        let r = run_atman(&[
            "bootstrap",
            "protein",
            "--input-dir",
            input.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--groups",
            "Case-Control",
            "--test",
            "welch-t",
            "--n",
            "200",
            "--seed",
            "123",
            "--min-pairs",
            "2",
        ]);
        assert!(r.status.success(), "bootstrap failed:\n{}", String::from_utf8_lossy(&r.stderr));
    }

    assert_byte_identical(&out1, &out2, "bootstrap protein");
    assert_sidecar_shas_match(&sidecar_for(&out1), &sidecar_for(&out2), "bootstrap sidecar");
}

#[test]
fn ratio_bootstrap_is_byte_deterministic() {
    let tmp1 = tempfile::tempdir().unwrap();
    let tmp2 = tempfile::tempdir().unwrap();
    let input1 = tmp1.path().join("input");
    let input2 = tmp2.path().join("input");
    let out1 = tmp1.path().join("ratio.tsv");
    let out2 = tmp2.path().join("ratio.tsv");
    write_ratio_fixture(&input1);
    write_ratio_fixture(&input2);

    for (input, out) in [(&input1, &out1), (&input2, &out2)] {
        let r = run_atman(&[
            "ratio",
            "--input-dir",
            input.to_str().unwrap(),
            "--numerator",
            "modules.intrathecal_IgV",
            "--denominator",
            "modules.plasma_IgG",
            "--groups",
            "MS-Control",
            "--test",
            "wilcoxon",
            "--n-bootstrap",
            "500",
            "--seed",
            "7",
            "--output",
            out.to_str().unwrap(),
        ]);
        assert!(r.status.success(), "ratio failed:\n{}", String::from_utf8_lossy(&r.stderr));
    }

    assert_byte_identical(&out1, &out2, "ratio");
    assert_sidecar_shas_match(&sidecar_for(&out1), &sidecar_for(&out2), "ratio sidecar");
}

#[test]
fn align_bootstrap_is_byte_deterministic() {
    let tmp1 = tempfile::tempdir().unwrap();
    let tmp2 = tempfile::tempdir().unwrap();
    let a1 = tmp1.path().join("cohort_a");
    let b1 = tmp1.path().join("cohort_b");
    let a2 = tmp2.path().join("cohort_a");
    let b2 = tmp2.path().join("cohort_b");
    write_align_cohort(&a1, "A", true, false, 11);
    write_align_cohort(&b1, "B", false, true, 22);
    write_align_cohort(&a2, "A", true, false, 11);
    write_align_cohort(&b2, "B", false, true, 22);
    let out1 = tmp1.path().join("summary.tsv");
    let out2 = tmp2.path().join("summary.tsv");

    for (a, b, out) in [(&a1, &b1, &out1), (&a2, &b2, &out2)] {
        let r = run_atman(&[
            "align",
            "bootstrap",
            "--cohorts",
            &format!("{},{}", a.display(), b.display()),
            "--labels",
            "A,B",
            "--k",
            "2",
            "--n-boot",
            "10",
            "--seed",
            "20260418",
            "--cosine-tau",
            "0.1",
            "--match-tau",
            "0.1",
            "--min-subjects",
            "10",
            "--max-iter",
            "80",
            "--tol",
            "1e-3",
            "--output",
            out.to_str().unwrap(),
        ]);
        assert!(
            r.status.success(),
            "align bootstrap failed:\n{}",
            String::from_utf8_lossy(&r.stderr)
        );
    }

    assert_byte_identical(&out1, &out2, "align bootstrap");
    assert_sidecar_shas_match(&sidecar_for(&out1), &sidecar_for(&out2), "align bootstrap sidecar");
}

#[test]
fn de_ensemble_is_byte_deterministic() {
    let tmp1 = tempfile::tempdir().unwrap();
    let tmp2 = tempfile::tempdir().unwrap();
    let input1 = tmp1.path().join("input");
    let input2 = tmp2.path().join("input");
    let out1 = tmp1.path().join("out");
    let out2 = tmp2.path().join("out");
    write_de_ensemble_fixture(&input1);
    write_de_ensemble_fixture(&input2);
    std::fs::create_dir_all(&out1).unwrap();
    std::fs::create_dir_all(&out2).unwrap();

    for (input, out) in [(&input1, &out1), (&input2, &out2)] {
        let r = run_atman(&[
            "de",
            "--input-dir",
            input.to_str().unwrap(),
            "--output-dir",
            out.to_str().unwrap(),
            "--test",
            "ensemble",
            "--groups",
            "A-B",
            "--ensemble-methods",
            "welch-t,limma",
            "--min-pairs",
            "3",
            "--trend",
            "false",
            "--robust",
            "false",
        ]);
        assert!(
            r.status.success(),
            "de ensemble failed:\n{}",
            String::from_utf8_lossy(&r.stderr)
        );
    }

    assert_byte_identical(
        &out1.join("de_results.tsv"),
        &out2.join("de_results.tsv"),
        "de_results",
    );
    assert_byte_identical(
        &out1.join("de_ensemble.tsv"),
        &out2.join("de_ensemble.tsv"),
        "de_ensemble",
    );
    assert_sidecar_shas_match(
        &sidecar_for(&out1.join("de_results.tsv")),
        &sidecar_for(&out2.join("de_results.tsv")),
        "de ensemble sidecar",
    );
}

#[test]
fn decompose_ica_is_byte_deterministic() {
    let tmp1 = tempfile::tempdir().unwrap();
    let tmp2 = tempfile::tempdir().unwrap();
    let input1 = tmp1.path().join("input");
    let input2 = tmp2.path().join("input");
    write_ica_fixture(&input1);
    write_ica_fixture(&input2);
    let loadings1 = tmp1.path().join("loadings.tsv");
    let activations1 = tmp1.path().join("activations.tsv");
    let stability1 = tmp1.path().join("stability.tsv");
    let loadings2 = tmp2.path().join("loadings.tsv");
    let activations2 = tmp2.path().join("activations.tsv");
    let stability2 = tmp2.path().join("stability.tsv");

    for (input, l, a, s) in [
        (&input1, &loadings1, &activations1, &stability1),
        (&input2, &loadings2, &activations2, &stability2),
    ] {
        let r = run_atman(&[
            "decompose",
            "ica",
            "--input-dir",
            input.to_str().unwrap(),
            "--k",
            "2",
            "--n-seeds",
            "5",
            "--seed",
            "20260418",
            "--stability-top-n",
            "10",
            "--output-loadings",
            l.to_str().unwrap(),
            "--output-activations",
            a.to_str().unwrap(),
            "--output-stability",
            s.to_str().unwrap(),
        ]);
        assert!(
            r.status.success(),
            "decompose ica failed:\n{}",
            String::from_utf8_lossy(&r.stderr)
        );
    }

    assert_byte_identical(&loadings1, &loadings2, "ica loadings");
    assert_byte_identical(&activations1, &activations2, "ica activations");
    assert_byte_identical(&stability1, &stability2, "ica stability");
    assert_sidecar_shas_match(
        &sidecar_for(&loadings1),
        &sidecar_for(&loadings2),
        "ica sidecar",
    );
}

#[test]
fn robust_paired_is_byte_deterministic() {
    // robust-paired is a leave-one-pair-out diagnostic: deterministic by
    // construction (no RNG, no seed). It is covered here for byte-identity to
    // guard against accidental nondeterminism (e.g. HashMap iteration order).
    let tmp1 = tempfile::tempdir().unwrap();
    let tmp2 = tempfile::tempdir().unwrap();
    let input1 = tmp1.path().join("input");
    let input2 = tmp2.path().join("input");
    write_robust_paired_fixture(&input1, 30);
    write_robust_paired_fixture(&input2, 30);
    let out1 = tmp1.path().join("robust.tsv");
    let out2 = tmp2.path().join("robust.tsv");

    for (input, out) in [(&input1, &out1), (&input2, &out2)] {
        let r = run_atman(&[
            "robust-paired",
            "--input-dir",
            input.to_str().unwrap(),
            "--groups",
            "tumor-paired_non_tumor",
            "--paired-by",
            "patient_id",
            "--min-pairs",
            "5",
            "--output",
            out.to_str().unwrap(),
        ]);
        assert!(
            r.status.success(),
            "robust-paired failed:\n{}",
            String::from_utf8_lossy(&r.stderr)
        );
    }

    assert_byte_identical(&out1, &out2, "robust-paired");
    // robust-paired also emits a <stem>_summary.<ext> companion; compare it
    // directly rather than relying solely on the sidecar SHA helper.
    assert_byte_identical(
        &summary_sibling(&out1),
        &summary_sibling(&out2),
        "robust-paired summary",
    );
    assert_sidecar_shas_match(&sidecar_for(&out1), &sidecar_for(&out2), "robust-paired sidecar");
}
