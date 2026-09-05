// Matrix-math variables use uppercase names (X, H, W) to mirror the standard
// NMF notation X ≈ W·H; keep them rather than forcing snake_case.
#![allow(non_snake_case)]

use std::collections::{BTreeMap, BTreeSet};
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

/// Reads an NMF loadings TSV (4 columns: program, assay_id, gene_symbol, loading)
/// and returns Vec<(program_id, BTreeMap<gene_symbol, loading_value>)> sorted by program.
fn read_loadings_tsv(path: &std::path::Path) -> Vec<(String, BTreeMap<String, f64>)> {
    let text = std::fs::read_to_string(path).expect("read loadings TSV");
    let lines: Vec<&str> = text.lines().collect();

    // Skip #-prefixed version comment line if present.
    let data_start = if lines.first().map(|l| l.starts_with('#')).unwrap_or(false) {
        1
    } else {
        0
    };

    // Expect header "program\tassay_id\tgene_symbol\tloading"
    assert!(
        lines
            .get(data_start)
            .map(|l| l.starts_with("program\t"))
            .unwrap_or(false),
        "invalid loadings header"
    );

    let mut programs: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();

    for line in &lines[(data_start + 1)..] {
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            parts.len(),
            4,
            "expected 4 columns in loadings row: {}",
            line
        );

        let program = parts[0].to_string();
        let gene_symbol = parts[2].to_string();
        let loading: f64 = parts[3].parse().expect("parse loading as f64");

        programs
            .entry(program)
            .or_default()
            .insert(gene_symbol, loading);
    }

    programs.into_iter().collect()
}

/// Computes Jaccard similarity of top-N genes (by absolute loading) between two gene-loading maps.
fn jaccard_top_n(map_a: &BTreeMap<String, f64>, map_b: &BTreeMap<String, f64>, n: usize) -> f64 {
    let mut genes_a: Vec<_> = map_a.iter().map(|(g, l)| (g.clone(), l.abs())).collect();
    let mut genes_b: Vec<_> = map_b.iter().map(|(g, l)| (g.clone(), l.abs())).collect();

    genes_a.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    genes_b.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let top_a: BTreeSet<String> = genes_a.iter().take(n).map(|(g, _)| g.clone()).collect();
    let top_b: BTreeSet<String> = genes_b.iter().take(n).map(|(g, _)| g.clone()).collect();

    let intersection = top_a.intersection(&top_b).collect::<BTreeSet<_>>().len();
    let union = top_a.union(&top_b).collect::<BTreeSet<_>>().len();

    if union == 0 {
        1.0
    } else {
        intersection as f64 / union as f64
    }
}

/// Best-match alignment of programs by Jaccard top-20.
/// Returns a mapping of reference_program_idx -> actual_program_idx.
fn align_programs_by_jaccard(
    reference: &[(String, BTreeMap<String, f64>)],
    actual: &[(String, BTreeMap<String, f64>)],
) -> Vec<usize> {
    let n_ref = reference.len();
    let n_actual = actual.len();
    assert_eq!(
        n_ref, n_actual,
        "reference and actual must have same number of programs"
    );

    // For each reference program, find best-matching actual program.
    let mut alignment = vec![0; n_ref];
    let mut used = vec![false; n_actual];

    for ref_idx in 0..n_ref {
        let ref_genes = &reference[ref_idx].1;
        let mut best_actual_idx = 0;
        let mut best_jaccard = -1.0;

        for act_idx in 0..n_actual {
            if used[act_idx] {
                continue;
            }
            let act_genes = &actual[act_idx].1;
            let j = jaccard_top_n(ref_genes, act_genes, 20);
            if j > best_jaccard {
                best_jaccard = j;
                best_actual_idx = act_idx;
            }
        }

        alignment[ref_idx] = best_actual_idx;
        used[best_actual_idx] = true;
    }

    alignment
}

/// Reads an NMF activations TSV (3 columns: sample_id, program, activation)
/// and returns BTreeMap<sample_id, BTreeMap<program_id, activation_value>>.
fn read_activations_tsv(path: &std::path::Path) -> BTreeMap<String, BTreeMap<String, f64>> {
    let text = std::fs::read_to_string(path).expect("read activations TSV");
    let lines: Vec<&str> = text.lines().collect();

    // Skip #-prefixed version comment line if present.
    let data_start = if lines.first().map(|l| l.starts_with('#')).unwrap_or(false) {
        1
    } else {
        0
    };

    // Expect header "sample_id\tprogram\tactivation"
    assert!(
        lines
            .get(data_start)
            .map(|l| l.starts_with("sample_id\t"))
            .unwrap_or(false),
        "invalid activations header"
    );

    let mut samples: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();

    for line in &lines[(data_start + 1)..] {
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }

        let sample_id = parts[0].to_string();
        let program = parts[1].to_string();
        let activation: f64 = parts[2].parse().expect("parse activation as f64");

        samples
            .entry(sample_id)
            .or_default()
            .insert(program, activation);
    }

    samples
}

/// Computes actual relative Frobenius reconstruction error: ||X - WH||_F / ||X||_F
/// where W is sample x program (activations), H is program x gene (loadings).
fn frobenius_reconstruction_error(
    data_tsv: &std::path::Path,
    activations: &BTreeMap<String, BTreeMap<String, f64>>,
    loadings: &[(String, BTreeMap<String, f64>)],
) -> f64 {
    // Read data: sample_id\tgene_symbol\tabundance (skip comment line).
    let text = std::fs::read_to_string(data_tsv).expect("read input data");
    let lines: Vec<&str> = text.lines().collect();

    // Skip #-prefixed version comment line if present.
    let data_start = if lines.first().map(|l| l.starts_with('#')).unwrap_or(false) {
        1
    } else {
        0
    };

    // Parse into a matrix: sample -> gene -> abundance.
    let mut X: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();

    for line in &lines[(data_start + 1)..] {
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }

        let sample_id = parts[0].to_string();
        let gene_symbol = parts[1].to_string();
        let abundance: f64 = parts[2].parse().unwrap_or(0.0);

        X.entry(sample_id)
            .or_default()
            .insert(gene_symbol, abundance);
    }

    // Build H as program_id -> gene_symbol -> loading.
    let mut H: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    for (program_id, genes) in loadings {
        H.insert(program_id.clone(), genes.clone());
    }

    // Compute X_hat[sample][gene] = sum_program W[sample][program] * H[program][gene].
    let mut X_hat: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();

    for (sample_id, programs) in activations {
        for (program_id, W_val) in programs {
            if let Some(gene_loadings) = H.get(program_id) {
                for (gene_symbol, H_val) in gene_loadings {
                    let reconstructed = W_val * H_val;
                    X_hat
                        .entry(sample_id.clone())
                        .or_default()
                        .entry(gene_symbol.clone())
                        .and_modify(|v| *v += reconstructed)
                        .or_insert(reconstructed);
                }
            }
        }
    }

    // Compute ||X - X_hat||_F.
    let mut error_sq = 0.0;
    let mut x_sq = 0.0;

    for (sample_id, genes) in &X {
        for (gene_symbol, x_val) in genes {
            x_sq += x_val * x_val;
            let x_hat_val = X_hat
                .get(sample_id)
                .and_then(|g| g.get(gene_symbol))
                .copied()
                .unwrap_or(0.0);
            let diff = x_val - x_hat_val;
            error_sq += diff * diff;
        }
    }

    // Handle entries in X_hat that were not in X (should be rare, but be complete).
    for (sample_id, genes) in &X_hat {
        for (gene_symbol, x_hat_val) in genes {
            if !X
                .get(sample_id)
                .map(|g| g.contains_key(gene_symbol))
                .unwrap_or(false)
            {
                error_sq += x_hat_val * x_hat_val;
            }
        }
    }

    let x_norm = x_sq.sqrt();
    if x_norm < 1e-10 {
        return 0.0;
    }

    error_sq.sqrt() / x_norm
}

/// Parameterized NMF parity test helper.
/// Runs atman with `--beta-loss <beta_loss>`, checks Jaccard top-20 >= 0.9 per program,
/// and checks relative Frobenius reconstruction error <= 5%.
/// Input fixture is always `nmf_frobenius_input.tsv` (loss-independent).
fn run_nmf_parity_test(beta_loss: &str, reference_tsv_path: &str) {
    let tmp = tempfile::tempdir().unwrap();
    let tmp_path = tmp.path();

    // Step 1: Read the long-format input fixture.
    let input_fixture_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/nmf_frobenius_input.tsv");
    let input_text =
        std::fs::read_to_string(&input_fixture_path).expect("read nmf_frobenius_input.tsv fixture");
    let input_lines: Vec<&str> = input_text.lines().collect();

    // First line is comment; skip it.
    let data_start = if input_lines[0].starts_with('#') {
        1
    } else {
        0
    };

    // Collect all unique samples and genes.
    let mut samples: BTreeSet<String> = BTreeSet::new();
    let mut genes: BTreeSet<String> = BTreeSet::new();
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
        genes.insert(gene_symbol.clone());
        measurements.push((sample_id, gene_symbol, abundance));
    }

    // Step 2: Build canonical Atman input directory.
    let canonical_dir = tmp_path.join("canonical");
    std::fs::create_dir(&canonical_dir).unwrap();

    // Write measurements.tsv in canonical long format.
    let headers = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
                    abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
                    detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order";
    let mut measurements_tsv = String::from(headers);
    measurements_tsv.push('\n');
    let mut ingest_order = 0u64;

    let genes_vec: Vec<String> = genes.iter().cloned().collect();
    for (sample_id, gene_symbol, abundance) in measurements.iter() {
        let gene_idx = genes_vec.iter().position(|g| g == gene_symbol).unwrap();
        let assay_id = format!("A{:03}", gene_idx);
        ingest_order += 1;
        let src = format!("{:.6}", abundance);
        measurements_tsv.push_str(&format!(
            "proteomics\t{}\t{}\t{}\tP1\t{}\t{:.6}\t{:.6}\tlog2_scale\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            sample_id, assay_id, gene_symbol, src, abundance, abundance, ingest_order
        ));
    }
    std::fs::write(canonical_dir.join("measurements.tsv"), &measurements_tsv).unwrap();

    // Write samples.tsv.
    let mut samples_tsv =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for (i, sample_id) in samples.iter().enumerate() {
        samples_tsv.push_str(&format!(
            "{}\t{}\tcase\t0\tcsf\t{}\n",
            sample_id,
            sample_id,
            i + 1
        ));
    }
    std::fs::write(canonical_dir.join("samples.tsv"), &samples_tsv).unwrap();

    // Write proteins.tsv.
    let mut proteins_tsv =
        String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for (gene_idx, gene) in genes_vec.iter().enumerate() {
        let assay_id = format!("A{:03}", gene_idx);
        proteins_tsv.push_str(&format!("proteomics\t{}\t\t{}\tP1\t\n", assay_id, gene));
    }
    std::fs::write(canonical_dir.join("proteins.tsv"), &proteins_tsv).unwrap();

    // Step 3: Run `atman decompose nmf`.
    let loadings_path = tmp_path.join("nmf_loadings.tsv");
    let activations_path = tmp_path.join("nmf_activations.tsv");
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
        beta_loss,
        "--init",
        "nndsvda",
        "--max-iter",
        "400",
        "--tol",
        "1e-6",
        "--seed",
        "42",
        "--output-loadings",
        loadings_path.to_str().unwrap(),
        "--output-activations",
        activations_path.to_str().unwrap(),
    ]);

    // Step 4: Check that atman succeeded.
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!(
            "atman decompose nmf --beta-loss {} failed:\nstderr:\n{}",
            beta_loss, stderr
        );
    }

    // Step 5: Read actual loadings and reference loadings.
    let actual_loadings = read_loadings_tsv(&loadings_path);
    let reference_fixture_path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(reference_tsv_path);
    let reference_loadings = read_loadings_tsv(&reference_fixture_path);

    // Step 6: Align programs by Jaccard top-20.
    let alignment = align_programs_by_jaccard(&reference_loadings, &actual_loadings);

    // Step 7: Check per-aligned-program Jaccard top-20 >= 0.9.
    for (ref_idx, act_idx) in alignment.iter().enumerate() {
        let ref_genes = &reference_loadings[ref_idx].1;
        let act_genes = &actual_loadings[*act_idx].1;
        let j = jaccard_top_n(ref_genes, act_genes, 20);
        assert!(
            j >= 0.9,
            "beta-loss={}: program {} Jaccard top-20 = {}, expected >= 0.9",
            beta_loss,
            ref_idx + 1,
            j
        );
    }

    // Step 8: Check relative Frobenius reconstruction error within 5%.
    // Read atman's activations and compute ||X - W_atman·H_atman||_F / ||X||_F.
    let actual_activations = read_activations_tsv(&activations_path);
    let act_error =
        frobenius_reconstruction_error(&input_fixture_path, &actual_activations, &actual_loadings);

    assert!(
        act_error <= 0.05,
        "beta-loss={}: Frobenius reconstruction error = {}, expected <= 0.05",
        beta_loss,
        act_error
    );
}

#[test]
fn nmf_frobenius_matches_sklearn_reference() {
    run_nmf_parity_test("frobenius", "tests/fixtures/nmf_frobenius_reference.tsv");
}

#[test]
fn nmf_kl_matches_sklearn_reference() {
    run_nmf_parity_test("kullback-leibler", "tests/fixtures/nmf_kl_reference.tsv");
}

/// Reads a stability TSV (3 columns: program, stable_seed_fraction, n_seeds_present).
fn read_stability_tsv(path: &std::path::Path) -> Vec<(String, f64, usize)> {
    let text = std::fs::read_to_string(path).expect("read stability TSV");
    let lines: Vec<&str> = text.lines().collect();

    // Skip #-prefixed version comment if present.
    let data_start = if lines.first().map(|l| l.starts_with('#')).unwrap_or(false) {
        1
    } else {
        0
    };

    assert!(
        lines
            .get(data_start)
            .map(|l| l.starts_with("program\t"))
            .unwrap_or(false),
        "invalid stability TSV header: {:?}",
        lines.get(data_start)
    );

    let mut rows = Vec::new();
    for line in &lines[(data_start + 1)..] {
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            parts.len(),
            3,
            "expected 3 columns in stability row: {}",
            line
        );
        let program = parts[0].to_string();
        let frac: f64 = parts[1].parse().expect("parse stable_seed_fraction");
        let n: usize = parts[2].parse().expect("parse n_seeds_present");
        rows.push((program, frac, n));
    }
    rows
}

/// Build a canonical Atman input directory from the standard NMF fixture and
/// return the path to the directory + a vec of (sample_id, gene_symbol, abundance).
fn build_canonical_dir(
    tmp_path: &std::path::Path,
) -> (std::path::PathBuf, Vec<(String, String, f64)>) {
    let input_fixture_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/nmf_frobenius_input.tsv");
    let input_text =
        std::fs::read_to_string(&input_fixture_path).expect("read nmf_frobenius_input.tsv");
    let input_lines: Vec<&str> = input_text.lines().collect();
    let data_start = if input_lines[0].starts_with('#') {
        1
    } else {
        0
    };

    let mut samples: BTreeSet<String> = BTreeSet::new();
    let mut genes: BTreeSet<String> = BTreeSet::new();
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
        genes.insert(gene_symbol.clone());
        measurements.push((sample_id, gene_symbol, abundance));
    }

    let canonical_dir = tmp_path.join("canonical");
    std::fs::create_dir(&canonical_dir).unwrap();

    let genes_vec: Vec<String> = genes.iter().cloned().collect();

    let headers = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
                   abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
                   detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order";
    let mut measurements_tsv = String::from(headers);
    measurements_tsv.push('\n');
    let mut ingest_order = 0u64;
    for (sample_id, gene_symbol, abundance) in &measurements {
        let gene_idx = genes_vec.iter().position(|g| g == gene_symbol).unwrap();
        let assay_id = format!("A{:03}", gene_idx);
        ingest_order += 1;
        let src = format!("{:.6}", abundance);
        measurements_tsv.push_str(&format!(
            "proteomics\t{}\t{}\t{}\tP1\t{}\t{:.6}\t{:.6}\tlog2_scale\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            sample_id, assay_id, gene_symbol, src, abundance, abundance, ingest_order
        ));
    }
    std::fs::write(canonical_dir.join("measurements.tsv"), &measurements_tsv).unwrap();

    let mut samples_tsv =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for (i, sample_id) in samples.iter().enumerate() {
        samples_tsv.push_str(&format!(
            "{}\t{}\tcase\t0\tcsf\t{}\n",
            sample_id,
            sample_id,
            i + 1
        ));
    }
    std::fs::write(canonical_dir.join("samples.tsv"), &samples_tsv).unwrap();

    let mut proteins_tsv =
        String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for (gene_idx, gene) in genes_vec.iter().enumerate() {
        let assay_id = format!("A{:03}", gene_idx);
        proteins_tsv.push_str(&format!("proteomics\t{}\t\t{}\tP1\t\n", assay_id, gene));
    }
    std::fs::write(canonical_dir.join("proteins.tsv"), &proteins_tsv).unwrap();

    (canonical_dir, measurements)
}

#[test]
fn nmf_multi_seed_emits_stability_and_filters_unstable_programs() {
    let tmp = tempfile::tempdir().unwrap();
    let tmp_path = tmp.path();

    let (canonical_dir, _measurements) = build_canonical_dir(tmp_path);

    let loadings_path = tmp_path.join("ms_loadings.tsv");
    let activations_path = tmp_path.join("ms_activations.tsv");
    let stability_path = tmp_path.join("ms_stability.tsv");

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
        "random",
        "--max-iter",
        "400",
        "--tol",
        "1e-6",
        "--n-seeds",
        "5",
        "--seed-base",
        "100",
        "--min-stable-seed-fraction",
        "0.9",
        "--stability-metric",
        "jaccard-top20",
        "--stability-top-n",
        "20",
        "--output-loadings",
        loadings_path.to_str().unwrap(),
        "--output-activations",
        activations_path.to_str().unwrap(),
        "--output-stability",
        stability_path.to_str().unwrap(),
    ]);

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!(
            "atman decompose nmf --n-seeds 5 failed:\nstderr:\n{}",
            stderr
        );
    }

    // ── 1. stability TSV must exist with the right columns ───────────────────
    assert!(
        stability_path.exists(),
        "stability TSV was not written to {:?}",
        stability_path
    );
    let stability_rows = read_stability_tsv(&stability_path);
    // The stability TSV reports ALL k=4 programs (pre-filter).
    assert_eq!(
        stability_rows.len(),
        4,
        "stability TSV should have 4 rows (one per original program)"
    );
    // Each n_seeds_present must be between 1 and 5.
    for (prog, frac, n_present) in &stability_rows {
        assert!(
            *n_present >= 1 && *n_present <= 5,
            "program {}: n_seeds_present={} out of range [1,5]",
            prog,
            n_present
        );
        assert!(
            *frac >= 0.0 && *frac <= 1.0,
            "program {}: stable_seed_fraction={} out of [0,1]",
            prog,
            frac
        );
    }

    // ── 2. loadings TSV: at most k=4 programs (filter never expands) ─────────
    let surviving_loadings = read_loadings_tsv(&loadings_path);
    let n_surviving = surviving_loadings.len();
    assert!(
        n_surviving <= 4,
        "loadings TSV has {} programs, expected <= 4 (k)",
        n_surviving
    );

    // ── 3. every program in loadings has stable_seed_fraction >= 0.9 ─────────
    // (The loadings TSV uses re-indexed program labels; match by position.)
    for (stab_prog, stab_frac, _) in stability_rows.iter().filter(|(_, frac, _)| *frac >= 0.9) {
        // Verify the fraction is non-NaN and non-negative.
        assert!(
            stab_frac.is_finite() && *stab_frac >= 0.0,
            "stability fraction for {} is invalid: {}",
            stab_prog,
            stab_frac
        );
    }
    // Surviving programs in loadings must all have been stable.
    let n_stable_in_stability = stability_rows
        .iter()
        .filter(|(_, frac, _)| *frac >= 0.9)
        .count();
    assert_eq!(
        n_surviving, n_stable_in_stability,
        "loadings TSV has {} programs but {} pass stability >= 0.9",
        n_surviving, n_stable_in_stability
    );

    // ── 4. activations TSV: samples × surviving programs ────────────────────
    let activations = read_activations_tsv(&activations_path);
    // Every surviving program should appear in activations for every sample.
    for (sample_id, prog_map) in &activations {
        assert_eq!(
            prog_map.len(),
            n_surviving,
            "sample {} has {} program activations, expected {}",
            sample_id,
            prog_map.len(),
            n_surviving
        );
    }
}

/// Reads a k-sweep TSV and returns Vec<(k, cophenetic, mean_rss, mean_kl_or_na, selected)>.
fn read_k_sweep_tsv(path: &std::path::Path) -> Vec<(usize, String, f64, String, u8)> {
    let text = std::fs::read_to_string(path).expect("read k-sweep TSV");
    let lines: Vec<&str> = text.lines().collect();

    // Skip #-prefixed comment lines.
    let data_start = if lines.first().map(|l| l.starts_with('#')).unwrap_or(false) {
        1
    } else {
        0
    };

    assert!(
        lines
            .get(data_start)
            .map(|l| l.starts_with("k\t"))
            .unwrap_or(false),
        "invalid k-sweep header: {:?}",
        lines.get(data_start)
    );

    let mut rows = Vec::new();
    for line in &lines[(data_start + 1)..] {
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            parts.len(),
            5,
            "expected 5 columns in k-sweep row: {}",
            line
        );
        let k: usize = parts[0].parse().expect("parse k");
        let cophenetic = parts[1].to_string();
        let mean_rss: f64 = parts[2].parse().expect("parse mean_rss");
        let mean_kl = parts[3].to_string();
        let selected: u8 = parts[4].parse().expect("parse selected");
        rows.push((k, cophenetic, mean_rss, mean_kl, selected));
    }
    rows
}

#[test]
fn nmf_k_selection_cophenetic_picks_planted_k() {
    // The existing fixture is rank-4 (generated by sklearn with n_components=4,
    // n_samples=30, n_features=50). We sweep k_min=2..k_max=7 and assert
    // selected_k == 4 and the sweep TSV has the right shape.
    let tmp = tempfile::tempdir().unwrap();
    let tmp_path = tmp.path();

    let (canonical_dir, _measurements) = build_canonical_dir(tmp_path);

    let loadings_path = tmp_path.join("ksel_loadings.tsv");
    let activations_path = tmp_path.join("ksel_activations.tsv");
    let sweep_path = tmp_path.join("ksel_sweep.tsv");

    let output = run_atman(&[
        "decompose",
        "nmf",
        "--input-dir",
        canonical_dir.to_str().unwrap(),
        "--k-selection",
        "cophenetic-knee",
        "--k-min",
        "2",
        "--k-max",
        "7",
        "--solver",
        "mu",
        "--beta-loss",
        "frobenius",
        "--init",
        "random",
        "--max-iter",
        "400",
        "--tol",
        "1e-6",
        "--n-seeds",
        "10",
        "--seed-base",
        "42",
        "--output-loadings",
        loadings_path.to_str().unwrap(),
        "--output-activations",
        activations_path.to_str().unwrap(),
        "--output-k-sweep",
        sweep_path.to_str().unwrap(),
    ]);

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!(
            "atman decompose nmf --k-selection cophenetic-knee failed:\nstderr:\n{}",
            stderr
        );
    }

    // ── 1. sweep TSV must exist and have k_max - k_min + 1 = 6 rows ─────────
    assert!(
        sweep_path.exists(),
        "k-sweep TSV was not written to {:?}",
        sweep_path
    );
    let sweep_rows = read_k_sweep_tsv(&sweep_path);
    assert_eq!(
        sweep_rows.len(),
        6, // k=2..=7
        "k-sweep TSV should have 6 rows (k_min=2..=k_max=7)"
    );

    // ── 2. exactly one row has selected=1 ───────────────────────────────────
    let selected_count = sweep_rows.iter().filter(|r| r.4 == 1).count();
    assert_eq!(selected_count, 1, "exactly one row should have selected=1");

    // ── 3. selected_k is 4 (ground truth of the rank-4 fixture) ─────────────
    let selected_k = sweep_rows
        .iter()
        .find(|r| r.4 == 1)
        .map(|r| r.0)
        .expect("no selected row found");

    assert_eq!(
        selected_k,
        4,
        "cophenetic-knee should select k=4 on the rank-4 fixture; selected k={}. \
         sweep: {:?}",
        selected_k,
        sweep_rows
            .iter()
            .map(|r| (r.0, r.1.as_str(), r.2))
            .collect::<Vec<_>>()
    );

    // ── 4. loadings TSV must exist with k=4 programs ────────────────────────
    assert!(
        loadings_path.exists(),
        "loadings TSV was not written to {:?}",
        loadings_path
    );
    let loadings = read_loadings_tsv(&loadings_path);
    assert_eq!(
        loadings.len(),
        4,
        "loadings TSV should have 4 programs (the selected k)"
    );
}

// ── `--max-missing-fraction` assay filter (mirrors `decompose ica`) ──────────

/// Builds a minimal canonical dir with `n_samples` samples and `n_assays`
/// assays (assay ids `A000`, `A001`, ... in row order), omitting the single
/// measurement at `(missing_assay_idx, missing_sample_idx)` to create exactly
/// one missing cell. All abundances are non-negative (required by NMF's
/// default `--transform none`).
fn build_canonical_dir_with_one_missing_cell(
    tmp_path: &std::path::Path,
    n_samples: usize,
    n_assays: usize,
    missing_assay_idx: usize,
    missing_sample_idx: usize,
) -> std::path::PathBuf {
    let input_dir = tmp_path.join("canonical_missing");
    std::fs::create_dir_all(&input_dir).unwrap();

    let headers = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
                   abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
                   detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order";
    let mut measurements_tsv = String::from(headers);
    measurements_tsv.push('\n');
    let mut ingest_order = 0u64;
    for i in 0..n_samples {
        let sample_id = format!("S{i:02}");
        for j in 0..n_assays {
            if j == missing_assay_idx && i == missing_sample_idx {
                continue; // intentional hole
            }
            let assay_id = format!("A{j:03}");
            let gene_symbol = format!("GENE{j:03}");
            let value = 1.0 + (i as f64) * 0.1 + (j as f64) * 0.05;
            ingest_order += 1;
            measurements_tsv.push_str(&format!(
                "proteomics\t{sample_id}\t{assay_id}\t{gene_symbol}\tP1\t{value:.6}\t{value:.6}\t{value:.6}\tlog2_scale\tPASS\tPASS\t\t0\t0\t\t\t{ingest_order}\n",
            ));
        }
    }
    std::fs::write(input_dir.join("measurements.tsv"), &measurements_tsv).unwrap();

    let mut samples_tsv =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 0..n_samples {
        samples_tsv.push_str(&format!("S{i:02}\tS{i:02}\tcase\t0\tcsf\t{}\n", i + 1));
    }
    std::fs::write(input_dir.join("samples.tsv"), samples_tsv).unwrap();

    let mut proteins_tsv =
        String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 0..n_assays {
        proteins_tsv.push_str(&format!("proteomics\tA{j:03}\t\tGENE{j:03}\tP1\t\n"));
    }
    std::fs::write(input_dir.join("proteins.tsv"), proteins_tsv).unwrap();

    input_dir
}

#[test]
fn nmf_max_missing_fraction_drops_incomplete_assays() {
    // Build a canonical dir where assay A001 is missing in 1 of 10 samples
    // (missing fraction 0.1) and all other assays are complete.
    let tmp = tempfile::tempdir().unwrap();
    let n_samples = 10;
    let n_assays = 5;
    let missing_assay_idx = 1; // A001
    let missing_sample_idx = 0; // S00
    let input_dir = build_canonical_dir_with_one_missing_cell(
        tmp.path(),
        n_samples,
        n_assays,
        missing_assay_idx,
        missing_sample_idx,
    );

    // With --max-missing-fraction 0.0 (the default), the incomplete assay
    // (missing fraction 0.1 > 0.0) must be dropped and the run must succeed.
    let loadings_default = tmp.path().join("loadings_default.tsv");
    let out_default = run_atman(&[
        "decompose",
        "nmf",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--k",
        "2",
        "--output-loadings",
        loadings_default.to_str().unwrap(),
    ]);
    assert!(
        out_default.status.success(),
        "default-threshold NMF should succeed by dropping the incomplete assay A001; stderr:\n{}",
        String::from_utf8_lossy(&out_default.stderr)
    );
    let loadings_text = std::fs::read_to_string(&loadings_default).unwrap();
    assert!(
        !loadings_text.contains("A001"),
        "loadings should not contain the dropped assay A001:\n{}",
        loadings_text
    );

    // With --max-missing-fraction 0.2, the assay's missing fraction (0.1) is
    // under the threshold so it is retained -- but NMF cannot take holes, so
    // the residual missing cell must cause a loud failure mentioning
    // "complete matrix".
    let loadings_lenient = tmp.path().join("loadings_lenient.tsv");
    let out_lenient = run_atman(&[
        "decompose",
        "nmf",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--k",
        "2",
        "--max-missing-fraction",
        "0.2",
        "--output-loadings",
        loadings_lenient.to_str().unwrap(),
    ]);
    assert!(
        !out_lenient.status.success(),
        "NMF with a retained-but-incomplete assay must fail, not silently succeed"
    );
    let stderr_lenient = String::from_utf8_lossy(&out_lenient.stderr);
    assert!(
        stderr_lenient.contains("complete matrix"),
        "expected stderr to mention 'complete matrix', got:\n{}",
        stderr_lenient
    );
}

// ── `--transform exp2-clip` / `shift-min` (Task 3) ───────────────────────────

/// Builds a minimal canonical dir with `n_samples` samples and `n_assays`
/// assays whose abundances span negative and positive values (simulating a
/// log2-ratio-to-reference scale, e.g. CPTAC). No missing cells. Used to
/// exercise `--transform exp2-clip` / `--transform shift-min`, which restore
/// non-negativity from signed input; `--transform none` must keep rejecting
/// this fixture.
fn build_canonical_dir_with_negative_values(
    tmp_path: &std::path::Path,
    n_samples: usize,
    n_assays: usize,
) -> std::path::PathBuf {
    let input_dir = tmp_path.join("canonical_negative");
    std::fs::create_dir_all(&input_dir).unwrap();

    let headers = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
                   abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
                   detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order";
    let mut measurements_tsv = String::from(headers);
    measurements_tsv.push('\n');
    let mut ingest_order = 0u64;
    for i in 0..n_samples {
        let sample_id = format!("S{i:02}");
        for j in 0..n_assays {
            let assay_id = format!("A{j:03}");
            let gene_symbol = format!("GENE{j:03}");
            // Log2-ratio-like: centered near 0, spans negative and positive.
            let value = (i as f64) * 0.6 - 3.0 + (j as f64) * 0.4;
            ingest_order += 1;
            measurements_tsv.push_str(&format!(
                "proteomics\t{sample_id}\t{assay_id}\t{gene_symbol}\tP1\t{value:.6}\t{value:.6}\t{value:.6}\tlog2_scale\tPASS\tPASS\t\t0\t0\t\t\t{ingest_order}\n",
            ));
        }
    }
    std::fs::write(input_dir.join("measurements.tsv"), &measurements_tsv).unwrap();

    let mut samples_tsv =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 0..n_samples {
        samples_tsv.push_str(&format!("S{i:02}\tS{i:02}\tcase\t0\tcsf\t{}\n", i + 1));
    }
    std::fs::write(input_dir.join("samples.tsv"), samples_tsv).unwrap();

    let mut proteins_tsv =
        String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 0..n_assays {
        proteins_tsv.push_str(&format!("proteomics\tA{j:03}\t\tGENE{j:03}\tP1\t\n"));
    }
    std::fs::write(input_dir.join("proteins.tsv"), proteins_tsv).unwrap();

    input_dir
}

/// Reads `<primary>.run.json` (the sidecar path convention from
/// `sidecar_path_for`) and parses it as JSON.
fn read_sidecar_json(primary_output: &std::path::Path) -> serde_json::Value {
    let sidecar_path = format!("{}.run.json", primary_output.display());
    let text = std::fs::read_to_string(&sidecar_path)
        .unwrap_or_else(|e| panic!("reading sidecar {sidecar_path:?}: {e}"));
    serde_json::from_str(&text).expect("parse sidecar JSON")
}

#[test]
fn nmf_transform_exp2_clip_succeeds_on_negative_input_and_records_sidecar() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = build_canonical_dir_with_negative_values(tmp.path(), 8, 6);

    let loadings_path = tmp.path().join("exp2clip_loadings.tsv");
    let output = run_atman(&[
        "decompose",
        "nmf",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--transform",
        "exp2-clip",
        "--transform-clamp",
        "6",
        "--k",
        "2",
        "--output-loadings",
        loadings_path.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "decompose nmf --transform exp2-clip should succeed on negative input; stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let sidecar = read_sidecar_json(&loadings_path);
    assert_eq!(sidecar["args"]["transform"], "exp2-clip");
    assert_eq!(
        sidecar["args"]["transform_clamp"].as_f64(),
        Some(6.0),
        "sidecar args: {}",
        sidecar["args"]
    );
    // exp2-clip never records a shift (that's shift-min's field).
    assert!(sidecar["args"]["transform_shift"].is_null());

    // Loadings must actually be non-negative-input-derived (i.e. NMF ran on
    // a transformed, non-negative matrix) — the run succeeding at all is
    // the load-bearing assertion here, but sanity-check the file is non-empty.
    let loadings_text = std::fs::read_to_string(&loadings_path).unwrap();
    assert!(loadings_text.lines().count() > 1, "loadings TSV is empty");
}

#[test]
fn nmf_transform_shift_min_succeeds_on_negative_input_and_records_shift() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = build_canonical_dir_with_negative_values(tmp.path(), 8, 6);

    let loadings_path = tmp.path().join("shiftmin_loadings.tsv");
    let output = run_atman(&[
        "decompose",
        "nmf",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--transform",
        "shift-min",
        "--k",
        "2",
        "--output-loadings",
        loadings_path.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "decompose nmf --transform shift-min should succeed on negative input; stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let sidecar = read_sidecar_json(&loadings_path);
    assert_eq!(sidecar["args"]["transform"], "shift-min");
    assert!(sidecar["args"]["transform_clamp"].is_null());
    // The fixture's minimum cell is at (sample=0, assay=0): -3.0.
    let shift = sidecar["args"]["transform_shift"]
        .as_f64()
        .expect("transform_shift should be a number");
    assert!(
        (shift - (-3.0)).abs() < 1e-9,
        "expected transform_shift == -3.0 (fixture's global min), got {shift}"
    );
}

#[test]
fn nmf_transform_none_still_rejects_negative_input() {
    // Same negative-valued fixture; `--transform none` (the default) must
    // keep failing with the existing non-negativity bail, unaffected by the
    // new exp2-clip / shift-min transforms.
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = build_canonical_dir_with_negative_values(tmp.path(), 8, 6);

    let loadings_path = tmp.path().join("none_loadings.tsv");
    let output = run_atman(&[
        "decompose",
        "nmf",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--transform",
        "none",
        "--k",
        "2",
        "--output-loadings",
        loadings_path.to_str().unwrap(),
    ]);
    assert!(
        !output.status.success(),
        "decompose nmf --transform none must keep rejecting negative input"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("requires non-negative input"),
        "expected the existing non-negativity bail message, got:\n{}",
        stderr
    );
}

#[test]
fn nmf_transform_clamp_rejected_unless_transform_is_exp2_clip() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = build_canonical_dir_with_negative_values(tmp.path(), 8, 6);

    let loadings_path = tmp.path().join("rejected_clamp_loadings.tsv");
    let output = run_atman(&[
        "decompose",
        "nmf",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--transform",
        "shift-min",
        "--transform-clamp",
        "3",
        "--k",
        "2",
        "--output-loadings",
        loadings_path.to_str().unwrap(),
    ]);
    assert!(
        !output.status.success(),
        "--transform-clamp with --transform shift-min must be rejected"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--transform-clamp") && stderr.contains("exp2-clip"),
        "expected an error naming --transform-clamp and exp2-clip, got:\n{}",
        stderr
    );
}
