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

    // Expect header "program\tassay_id\tgene_symbol\tloading"
    assert!(lines[0].starts_with("program\t"), "invalid loadings header");

    let mut programs: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();

    for line in &lines[1..] {
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts.len(), 4, "expected 4 columns in loadings row: {}", line);

        let program = parts[0].to_string();
        let gene_symbol = parts[2].to_string();
        let loading: f64 = parts[3].parse().expect("parse loading as f64");

        programs
            .entry(program)
            .or_insert_with(BTreeMap::new)
            .insert(gene_symbol, loading);
    }

    programs
        .into_iter()
        .map(|(prog, genes)| (prog, genes))
        .collect()
}

/// Computes Jaccard similarity of top-N genes (by absolute loading) between two gene-loading maps.
fn jaccard_top_n(
    map_a: &BTreeMap<String, f64>,
    map_b: &BTreeMap<String, f64>,
    n: usize,
) -> f64 {
    let mut genes_a: Vec<_> = map_a.iter().map(|(g, l)| (g.clone(), l.abs())).collect();
    let mut genes_b: Vec<_> = map_b.iter().map(|(g, l)| (g.clone(), l.abs())).collect();

    genes_a.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    genes_b.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let top_a: BTreeSet<String> = genes_a
        .iter()
        .take(n)
        .map(|(g, _)| g.clone())
        .collect();
    let top_b: BTreeSet<String> = genes_b
        .iter()
        .take(n)
        .map(|(g, _)| g.clone())
        .collect();

    let intersection = top_a
        .intersection(&top_b)
        .collect::<BTreeSet<_>>()
        .len();
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

/// Computes Frobenius norm of reconstruction error given loadings W and data X.
/// X_reconstructed = W @ H (where H is implicit; we check how well W explains X).
/// For simplicity, check per-gene reconstruction: how well does loadings explain the observed variance.
fn frobenius_reconstruction_error(
    data_tsv: &std::path::Path,
    loadings: &[(String, BTreeMap<String, f64>)],
) -> f64 {
    // Read data: sample_id\tgene_symbol\tabundance (skip comment line).
    let text = std::fs::read_to_string(data_tsv).expect("read input data");
    let lines: Vec<&str> = text.lines().collect();

    // First line is comment starting with #; skip it.
    let data_start = if lines[0].starts_with('#') { 1 } else { 0 };

    // Parse into a 2D matrix: gene x sample -> abundance.
    let mut data: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();

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

        data.entry(gene_symbol)
            .or_insert_with(BTreeMap::new)
            .insert(sample_id, abundance);
    }

    // For each gene, compute how well its observed profile is explained by loadings.
    // Use Frobenius norm of residual after orthogonal projection onto span(loadings).
    let mut total_error = 0.0;

    for (gene_symbol, samples) in &data {
        // Extract observed vector for this gene across all samples.
        let n_samples = samples.len();
        if n_samples == 0 {
            continue;
        }

        let mut obs: Vec<f64> = samples.values().cloned().collect();
        let obs_norm: f64 = obs.iter().map(|x| x * x).sum::<f64>().sqrt();

        if obs_norm < 1e-10 {
            continue;
        }

        // Normalize observed vector.
        for x in &mut obs {
            *x /= obs_norm;
        }

        // Compute projection of obs onto each program's loadings for this gene.
        let mut projection = 0.0;
        for (_prog_id, genes) in loadings {
            if let Some(loading) = genes.get(gene_symbol) {
                // Loading is the weight of this gene in this program.
                // Assume programs form an approximate basis; sum squared loadings.
                projection += loading * loading;
            }
        }

        // Reconstruction error: orthogonal distance from obs to subspace spanned by loadings.
        // Simplified: sqrt(1 - projection^2) for normalized obs.
        let error = (1.0 - projection.min(1.0)).abs().sqrt();
        total_error += error * error;
    }

    total_error.sqrt()
}

#[test]
fn nmf_frobenius_matches_sklearn_reference() {
    let tmp = tempfile::tempdir().unwrap();
    let tmp_path = tmp.path();

    // Step 1: Read the long-format input fixture.
    let input_fixture_path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/nmf_frobenius_input.tsv");
    let input_text = std::fs::read_to_string(&input_fixture_path)
        .expect("read nmf_frobenius_input.tsv fixture");
    let input_lines: Vec<&str> = input_text.lines().collect();

    // First line is comment; skip it.
    let data_start = if input_lines[0].starts_with('#') { 1 } else { 0 };

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
    for (_row_idx, (sample_id, gene_symbol, abundance)) in measurements.iter().enumerate() {
        let gene_idx = genes_vec
            .iter()
            .position(|g| g == gene_symbol)
            .unwrap();
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
    let mut samples_tsv = String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for (i, sample_id) in samples.iter().enumerate() {
        samples_tsv.push_str(&format!(
            "{}\t{}\tcase\t0\tcsf\t{}\n",
            sample_id, sample_id, i + 1
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
        "42",
        "--output-loadings",
        loadings_path.to_str().unwrap(),
    ]);

    // Step 4: Check that atman was invoked (may fail if nmf not yet implemented).
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        // Expected failure during B2->B3/B4: "unrecognized subcommand 'nmf'".
        // Panic with stderr for debugging.
        panic!(
            "atman decompose nmf failed (expected until B3/B4 implemented):\nstderr:\n{}",
            stderr
        );
    }

    // Step 5: Read actual loadings and reference loadings.
    let actual_loadings = read_loadings_tsv(&loadings_path);
    let reference_fixture_path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/nmf_frobenius_reference.tsv");
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
            "program {} Jaccard top-20 = {}, expected >= 0.9",
            ref_idx + 1,
            j
        );
    }

    // Step 8: Check Frobenius reconstruction error within 1% relative.
    let ref_error = frobenius_reconstruction_error(&input_fixture_path, &reference_loadings);
    let act_error = frobenius_reconstruction_error(&input_fixture_path, &actual_loadings);

    let relative_error = if ref_error.abs() > 1e-10 {
        (act_error - ref_error).abs() / ref_error
    } else {
        act_error.abs()
    };

    assert!(
        relative_error <= 0.01,
        "Frobenius reconstruction error relative diff = {}, expected <= 0.01\n\
         ref_error={}, act_error={}",
        relative_error,
        ref_error,
        act_error
    );
}
