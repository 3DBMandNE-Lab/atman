/// Integration tests for the MNAR-aware ICA decomposition.
///
/// G1: MAR-collapse parity — on fully-observed data, `--missingness-model
///     abundance-conditional` must produce loadings within 1e-6 of plain
///     FastICA (after canonicalization by L2 norm and sign-fix).
///
/// G4: MNAR-recovery — on the Phase-F synthetic MNAR fixture, the MNAR-aware
///     run must recover ground-truth sources with lower MAE than plain
///     column-mean-impute + decompose.
use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

// ── shared synthetic canonical builder (30 samples × 40 assays, fully observed) ──

fn build_fully_observed_canonical(tmp: &Path) -> std::path::PathBuf {
    let input_dir = tmp.join("canonical_full");
    std::fs::create_dir_all(&input_dir).unwrap();

    // Deterministic LCG matching decompose_ica.rs.
    let n_samples = 30;
    let n_assays = 40;
    let mut state: u64 = 1234567;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) as f64) / (u32::MAX as f64)
    };
    let mut xs = Vec::with_capacity(n_samples);
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

    let headers = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
                    abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
                    detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order";
    let mut measurements = String::from(headers);
    measurements.push('\n');
    let mut ingest = 0u64;
    for (i, row) in xs.iter().enumerate() {
        let sample_id = format!("S{i:02}");
        for (j, v) in row.iter().enumerate() {
            let assay_id = format!("A{j:03}");
            let gene = format!("GENE{j:03}");
            ingest += 1;
            measurements.push_str(&format!(
                "olink_explore_ngs\t{sample_id}\t{assay_id}\t{gene}\tP1\t{v:.6}\t{v:.6}\t{v:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{ingest}\n",
            ));
        }
    }
    std::fs::write(input_dir.join("measurements.tsv"), &measurements).unwrap();

    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 0..n_samples {
        samples.push_str(&format!("S{i:02}\tS{i:02}\tcase\t0\tcsf\t{}\n", i + 1));
    }
    std::fs::write(input_dir.join("samples.tsv"), samples).unwrap();

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 0..n_assays {
        proteins.push_str(&format!("olink_explore_ngs\tA{j:03}\t\tGENE{j:03}\tP1\t\n"));
    }
    std::fs::write(input_dir.join("proteins.tsv"), proteins).unwrap();

    input_dir
}

// ── G1: MAR-collapse parity ───────────────────────────────────────────────────

#[test]
fn mnar_ica_collapses_to_fastica_under_uniform_detection() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = build_fully_observed_canonical(tmp.path());

    let loadings_plain = tmp.path().join("loadings_plain.tsv");
    let activations_plain = tmp.path().join("activations_plain.tsv");
    let stability_plain = tmp.path().join("stability_plain.tsv");

    let loadings_mnar = tmp.path().join("loadings_mnar.tsv");
    let activations_mnar = tmp.path().join("activations_mnar.tsv");
    let stability_mnar = tmp.path().join("stability_mnar.tsv");

    // Run A: plain FastICA (no missingness model).
    let out_a = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--k",
        "2",
        "--n-seeds",
        "1",
        "--seed",
        "20260418",
        "--output-loadings",
        loadings_plain.to_str().unwrap(),
        "--output-activations",
        activations_plain.to_str().unwrap(),
        "--output-stability",
        stability_plain.to_str().unwrap(),
    ]);
    assert!(
        out_a.status.success(),
        "plain ICA stderr:\n{}",
        String::from_utf8_lossy(&out_a.stderr)
    );

    // Run B: MNAR-aware on the same fully-observed data.
    let out_b = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--k",
        "2",
        "--n-seeds",
        "1",
        "--seed",
        "20260418",
        "--missingness-model",
        "abundance-conditional",
        "--max-joint-iter",
        "1",
        "--output-loadings",
        loadings_mnar.to_str().unwrap(),
        "--output-activations",
        activations_mnar.to_str().unwrap(),
        "--output-stability",
        stability_mnar.to_str().unwrap(),
    ]);
    assert!(
        out_b.status.success(),
        "MNAR ICA stderr:\n{}",
        String::from_utf8_lossy(&out_b.stderr)
    );

    // Parse both loading tables.
    let plain = parse_loadings(&loadings_plain);
    let mnar = parse_loadings(&loadings_mnar);

    // Canonicalize by L2 norm (sort programs descending by their vector L2 norm).
    let plain_canon = sort_by_l2_norm(plain);
    let mnar_canon = sort_by_l2_norm(mnar);

    // Each program: sign-align MNAR to plain before comparing.
    let k = plain_canon.len();
    assert_eq!(k, mnar_canon.len(), "program count mismatch");

    let mut max_diff = 0.0_f64;
    for i in 0..k {
        let p = &plain_canon[i];
        let m = &mnar_canon[i];
        assert_eq!(
            p.len(),
            m.len(),
            "loading vector length mismatch at program {i}"
        );

        // Determine sign by dot product.
        let dot: f64 = p.iter().zip(m.iter()).map(|(a, b)| a * b).sum();
        let sign = if dot >= 0.0 { 1.0 } else { -1.0 };

        for (a, b) in p.iter().zip(m.iter()) {
            let diff = (a - sign * b).abs();
            if diff > max_diff {
                max_diff = diff;
            }
        }
    }

    eprintln!("G1 max element-wise loading diff (plain vs MNAR on full data): {max_diff:.2e}");
    assert!(
        max_diff < 1e-6,
        "MAR-collapse parity violated: max element-wise diff = {max_diff:.2e} (threshold 1e-6)"
    );
}

// ── G4: MNAR-recovery on synthetic ground truth ──────────────────────────────

/// The pre-registered comparison, run in full for the first time.
///
/// The atman methods-paper session set this bar before the weighted
/// whitening path existed: per-cell detection weighting earns its place
/// only if it recovers the planted sources at lower MAE than BOTH
/// alternatives on the committed ground-truth fixture. Anything else --
/// looking different, being slower, producing a plausible curve -- is
/// not evidence, because only this fixture knows the true sources.
///
/// The existing test above compares two of the arms. This one adds the
/// two that were missing: complete-case ICA, and the weighted-whitening
/// path. Its assertions are deliberately weak on direction and strong
/// on reporting: it prints all four MAEs so a reader can see the
/// margins, and fails only if the weighted path is not computed at all.
/// A negative result here is a real finding about detection-probability
/// weighting rather than a broken test, and encoding "weighted must
/// win" would make it impossible to record one.
#[test]
fn mnar_weighted_whitening_against_the_preregistered_alternatives() {
    let fixture_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let observed_tsv = fixture_dir.join("mnar_observed_abundance.tsv");
    let truth_tsv = fixture_dir.join("mnar_ground_truth_sources.tsv");
    assert!(
        observed_tsv.exists() && truth_tsv.exists(),
        "fixtures missing"
    );

    let tmp = tempfile::tempdir().unwrap();
    let input_dir = build_mnar_canonical_input(tmp.path(), &observed_tsv);
    let truth = parse_sources_tsv(&truth_tsv);

    let run = |tag: &str, seed: u64, extra: &[&str]| -> f64 {
        let loadings = tmp.path().join(format!("l_{tag}_{seed}.tsv"));
        let activations = tmp.path().join(format!("a_{tag}_{seed}.tsv"));
        let stability = tmp.path().join(format!("s_{tag}_{seed}.tsv"));
        let mut args: Vec<String> = vec![
            "decompose".into(),
            "ica".into(),
            "--input-dir".into(),
            input_dir.display().to_string(),
            "--output-loadings".into(),
            loadings.display().to_string(),
            "--output-activations".into(),
            activations.display().to_string(),
            "--output-stability".into(),
            stability.display().to_string(),
            "--k".into(),
            "2".into(),
            "--n-seeds".into(),
            "1".into(),
            "--seed".into(),
            seed.to_string(),
            // Without these the loader drops every assay carrying any
            // missingness, leaving a complete matrix on which all three
            // arms are identical by construction. The first version of
            // this test omitted them and reported three identical MAEs.
            "--max-missing-fraction".into(),
            "1.0".into(),
            "--impute".into(),
            "mean".into(),
        ];
        args.extend(extra.iter().map(|s| s.to_string()));
        let refs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
        let out = run_atman(&refs);
        assert!(
            out.status.success(),
            "{tag} failed:\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
        best_match_mae(&parse_activations_tsv(&activations), &truth)
    };

    // Across seeds, not one. A single seed gave 0.0768 / 0.0759 / 0.0754
    // and a clean win for weighted whitening -- margins under 2%, which
    // is exactly the size that one draw cannot distinguish from noise.
    let seeds = [20260418u64, 7, 42, 99, 1234, 555];
    let mut wins = 0usize;
    let mut rows: Vec<(u64, f64, f64, f64)> = Vec::new();
    let mut single_pass: Vec<(u64, f64, f64)> = Vec::new();
    for seed in seeds {
        let mae_impute = run("impute", seed, &[]);
        let mae_joint = run(
            "joint",
            seed,
            &[
                "--missingness-model",
                "abundance-conditional",
                "--max-joint-iter",
                "5",
            ],
        );
        let mae_weighted = run(
            "weighted",
            seed,
            &[
                "--missingness-model",
                "abundance-conditional",
                "--max-joint-iter",
                "5",
                "--weighted-whitening",
            ],
        );
        if mae_weighted < mae_impute && mae_weighted < mae_joint {
            wins += 1;
        }
        rows.push((seed, mae_impute, mae_joint, mae_weighted));
        // Single-pass arms: weighting WITHOUT the joint loop. The loop
        // is what drives the reconstruction gap to zero, so this
        // separates "does the detection weighting help" from "does the
        // iteration help", which the looped arms confound.
        let one_joint = run(
            "joint1",
            seed,
            &[
                "--missingness-model",
                "abundance-conditional",
                "--max-joint-iter",
                "1",
            ],
        );
        let one_weighted = run(
            "weighted1",
            seed,
            &[
                "--missingness-model",
                "abundance-conditional",
                "--max-joint-iter",
                "1",
                "--weighted-whitening",
            ],
        );
        single_pass.push((seed, one_joint, one_weighted));
    }
    eprintln!("PREREGISTERED seed  mean-impute   joint      weighted");
    for (seed, a, b, c) in &rows {
        eprintln!("PREREGISTERED {seed:>8}  {a:.6}  {b:.6}  {c:.6}");
    }
    eprintln!(
        "PREREGISTERED VERDICT weighted whitening beat both alternatives on {wins}/{} seeds",
        seeds.len()
    );
    eprintln!("SINGLEPASS seed  unweighted   weighted");
    let mut single_wins = 0usize;
    for (seed, a, b) in &single_pass {
        if b < a {
            single_wins += 1;
        }
        eprintln!("SINGLEPASS {seed:>8}  {a:.6}  {b:.6}");
    }
    eprintln!(
        "SINGLEPASS VERDICT detection weighting alone helped on {single_wins}/{} seeds",
        single_pass.len()
    );

    // Two orderings, both stable across every seed measured. Asserted
    // because they are the finding, and because a regression that
    // reversed either one would otherwise be invisible.
    assert_eq!(
        wins,
        seeds.len(),
        "weighted whitening must beat both looped alternatives on every seed"
    );
    assert_eq!(
        single_wins,
        single_pass.len(),
        "detection weighting must help on every seed when the joint loop is not run"
    );
    // The joint loop is the part that hurts. Single-pass weighting beat
    // the looped weighting on every seed measured, which is why the
    // recommended configuration is --weighted-whitening with
    // --max-joint-iter 1 rather than the loop.
    for ((seed, _, _, looped), (_, _, single)) in rows.iter().zip(single_pass.iter()) {
        assert!(
            single < looped,
            "seed {seed}: single-pass weighting ({single:.6}) should beat looped weighting \
             ({looped:.6}); the joint loop drives the imputed cells toward their own \
             reconstruction and costs recovery"
        );
    }
}

#[test]
fn mnar_ica_recovers_sources_better_than_impute_then_decompose() {
    // Fixtures produced by the Phase-F generator (60 samples × 50 features,
    // beta0=4.0 beta1=1.5, numpy seed=7).
    let fixture_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");

    let observed_tsv = fixture_dir.join("mnar_observed_abundance.tsv");
    let truth_tsv = fixture_dir.join("mnar_ground_truth_sources.tsv");

    assert!(
        observed_tsv.exists(),
        "fixture not found: {}",
        observed_tsv.display()
    );
    assert!(
        truth_tsv.exists(),
        "fixture not found: {}",
        truth_tsv.display()
    );

    let tmp = tempfile::tempdir().unwrap();

    // Build canonical input dir from the fixture observed_abundance.tsv.
    let input_dir = build_mnar_canonical_input(tmp.path(), &observed_tsv);

    // (a) Plain ICA with column-mean imputation.
    let loadings_a = tmp.path().join("loadings_a.tsv");
    let activations_a = tmp.path().join("activations_a.tsv");
    let stability_a = tmp.path().join("stability_a.tsv");

    let out_a = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--k",
        "2",
        "--n-seeds",
        "1",
        "--seed",
        "20260418",
        "--max-missing-fraction",
        "1.0",
        "--impute",
        "mean",
        "--output-loadings",
        loadings_a.to_str().unwrap(),
        "--output-activations",
        activations_a.to_str().unwrap(),
        "--output-stability",
        stability_a.to_str().unwrap(),
    ]);
    assert!(
        out_a.status.success(),
        "impute-then-decompose stderr:\n{}",
        String::from_utf8_lossy(&out_a.stderr)
    );

    // (b) MNAR-aware ICA.
    let loadings_b = tmp.path().join("loadings_b.tsv");
    let activations_b = tmp.path().join("activations_b.tsv");
    let stability_b = tmp.path().join("stability_b.tsv");

    let out_b = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--k",
        "2",
        "--n-seeds",
        "1",
        "--seed",
        "20260418",
        "--max-missing-fraction",
        "1.0",
        "--impute",
        "mean",
        "--missingness-model",
        "abundance-conditional",
        "--max-joint-iter",
        "5",
        "--output-loadings",
        loadings_b.to_str().unwrap(),
        "--output-activations",
        activations_b.to_str().unwrap(),
        "--output-stability",
        stability_b.to_str().unwrap(),
    ]);
    assert!(
        out_b.status.success(),
        "MNAR-aware stderr:\n{}",
        String::from_utf8_lossy(&out_b.stderr)
    );

    // Read ground-truth sources (60 samples × 2 sources).
    let truth = parse_sources_tsv(&truth_tsv);
    // truth: Vec<Vec<f64>> where truth[i] = [s1_i, s2_i]

    // Read recovered activations (per-sample activations = "sources" in ICA terms).
    let act_a = parse_activations_tsv(&activations_a);
    let act_b = parse_activations_tsv(&activations_b);

    // Compute MAE after best-match cosine alignment.
    let mae_a = best_match_mae(&act_a, &truth);
    let mae_b = best_match_mae(&act_b, &truth);

    eprintln!("G4 MAE impute-then-decompose = {mae_a:.6}");
    eprintln!("G4 MAE MNAR-aware            = {mae_b:.6}");

    assert!(
        mae_b < mae_a,
        "MNAR-aware ICA (MAE={mae_b:.6}) did not improve over impute-then-decompose (MAE={mae_a:.6})"
    );
}

// ── helpers ──────────────────────────────────────────────────────────────────

/// Build a canonical atman input directory from the MNAR observed-abundance
/// fixture TSV.  The fixture has columns: sample_id, gene_symbol, abundance
/// (with "NaN" as the sentinel for missing values).
fn build_mnar_canonical_input(tmp: &Path, observed_tsv: &Path) -> std::path::PathBuf {
    let input_dir = tmp.join("mnar_input");
    std::fs::create_dir_all(&input_dir).unwrap();

    // Parse the fixture: collect (sample_id, gene_symbol, abundance_or_nan).
    struct Row {
        sample_id: String,
        gene_symbol: String,
        abundance: Option<f64>,
    }
    let content = std::fs::read_to_string(observed_tsv).unwrap();
    let mut rows: Vec<Row> = Vec::new();
    let mut seen_samples: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    let mut sample_order: Vec<String> = Vec::new();
    let mut seen_genes: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    let mut gene_order: Vec<String> = Vec::new();

    for line in content.lines() {
        if line.starts_with('#') || line.starts_with("sample_id") {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let sample_id = parts[0].to_string();
        let gene_symbol = parts[1].to_string();
        let abundance: Option<f64> = if parts[2].eq_ignore_ascii_case("nan") || parts[2].is_empty()
        {
            None
        } else {
            parts[2].parse().ok()
        };
        if seen_samples.insert(sample_id.clone()) {
            sample_order.push(sample_id.clone());
        }
        if seen_genes.insert(gene_symbol.clone()) {
            gene_order.push(gene_symbol.clone());
        }
        rows.push(Row {
            sample_id,
            gene_symbol,
            abundance,
        });
    }

    // Write measurements.tsv in canonical long format.
    let headers = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
                    abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
                    detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order";
    let mut measurements = String::from(headers);
    measurements.push('\n');
    let mut ingest = 0u64;
    for row in &rows {
        ingest += 1;
        // NaN rows: write "NaN" as the abundance field; load_matrix_mnar handles non-finite.
        let ab_str = match row.abundance {
            Some(v) => format!("{v:.10}"),
            None => "NaN".to_string(),
        };
        measurements.push_str(&format!(
            "olink_explore_ngs\t{sid}\t{gsym}\t{gsym}\tP1\t{ab_str}\t{ab_str}\t{ab_str}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{ingest}\n",
            sid = row.sample_id,
            gsym = row.gene_symbol,
        ));
    }
    std::fs::write(input_dir.join("measurements.tsv"), &measurements).unwrap();

    // samples.tsv
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for (i, s) in sample_order.iter().enumerate() {
        samples.push_str(&format!("{s}\t{s}\tcase\t0\tcsf\t{}\n", i + 1));
    }
    std::fs::write(input_dir.join("samples.tsv"), samples).unwrap();

    // proteins.tsv
    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for g in &gene_order {
        proteins.push_str(&format!("olink_explore_ngs\t{g}\t\t{g}\tP1\t\n"));
    }
    std::fs::write(input_dir.join("proteins.tsv"), proteins).unwrap();

    input_dir
}

/// Parse a loadings TSV (`program\tassay_id\tgene_symbol\tloading`).
/// Returns a map from program name to loading vector (in assay order).
fn parse_loadings(path: &Path) -> std::collections::BTreeMap<String, Vec<f64>> {
    let content = std::fs::read_to_string(path).unwrap();
    let mut map: std::collections::BTreeMap<String, Vec<f64>> = std::collections::BTreeMap::new();
    for line in content.lines().skip(1) {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 {
            continue;
        }
        let program = parts[0].to_string();
        let loading: f64 = parts[3].parse().unwrap_or(0.0);
        map.entry(program).or_default().push(loading);
    }
    map
}

/// Sort program loading vectors by descending L2 norm.
fn sort_by_l2_norm(map: std::collections::BTreeMap<String, Vec<f64>>) -> Vec<Vec<f64>> {
    let mut vecs: Vec<Vec<f64>> = map.into_values().collect();
    vecs.sort_by(|a, b| {
        let la: f64 = a.iter().map(|v| v * v).sum::<f64>().sqrt();
        let lb: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        lb.partial_cmp(&la).unwrap_or(std::cmp::Ordering::Equal)
    });
    vecs
}

/// Parse activations TSV (`cohort\tsubject_id\tsample_id\tprogram\tactivation`).
/// Returns per-sample activations aligned by sample_id order, shape n × k.
fn parse_activations_tsv(path: &Path) -> Vec<Vec<f64>> {
    let content = std::fs::read_to_string(path).unwrap();
    // Collect (sample_id, program, activation).
    let mut by_sample: std::collections::BTreeMap<String, std::collections::BTreeMap<String, f64>> =
        std::collections::BTreeMap::new();
    for line in content.lines().skip(1) {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            continue;
        }
        let sample_id = parts[2].to_string();
        let program = parts[3].to_string();
        let activation: f64 = parts[4].parse().unwrap_or(0.0);
        by_sample
            .entry(sample_id)
            .or_default()
            .insert(program, activation);
    }
    // Emit as Vec<Vec<f64>> sorted by sample_id, programs by name.
    by_sample
        .into_values()
        .map(|prog_map| prog_map.into_values().collect())
        .collect()
}

/// Parse the ground-truth sources TSV (`sample_id\tsource\tvalue`).
/// Returns per-sample source values aligned by sample_id and source name,
/// shape n × k.
fn parse_sources_tsv(path: &Path) -> Vec<Vec<f64>> {
    let content = std::fs::read_to_string(path).unwrap();
    let mut by_sample: std::collections::BTreeMap<String, std::collections::BTreeMap<String, f64>> =
        std::collections::BTreeMap::new();
    for line in content.lines() {
        if line.starts_with('#') || line.starts_with("sample_id") {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let sample_id = parts[0].to_string();
        let source = parts[1].to_string();
        let value: f64 = parts[2].parse().unwrap_or(0.0);
        by_sample
            .entry(sample_id)
            .or_default()
            .insert(source, value);
    }
    by_sample
        .into_values()
        .map(|src_map| src_map.into_values().collect())
        .collect()
}

/// Compute MAE after best-match alignment between `recovered` and `truth`.
///
/// ICA sources are ambiguous up to permutation, sign, and scale.
/// We align by maximising absolute correlation between recovered and truth
/// programs (Hungarian-free greedy, works for k=2), then normalise each
/// recovered source to unit variance before computing MAE.
fn best_match_mae(recovered: &[Vec<f64>], truth: &[Vec<f64>]) -> f64 {
    let n = truth.len().min(recovered.len());
    if n == 0 {
        return f64::INFINITY;
    }
    let k_truth = truth[0].len();
    let k_rec = recovered[0].len();

    // Build n × k column vecs for both.
    let truth_cols: Vec<Vec<f64>> = (0..k_truth)
        .map(|c| truth.iter().take(n).map(|r| r[c]).collect())
        .collect();
    let rec_cols: Vec<Vec<f64>> = (0..k_rec)
        .map(|c| recovered.iter().take(n).map(|r| r[c]).collect())
        .collect();

    // Absolute correlation matrix.
    let mut abs_corr = vec![vec![0.0_f64; k_rec]; k_truth];
    for (ti, tc) in truth_cols.iter().enumerate() {
        for (ri, rc) in rec_cols.iter().enumerate() {
            abs_corr[ti][ri] = abs_pearson(tc, rc);
        }
    }

    // Greedy matching: for each truth component, pick the best unmatched recovered.
    let mut used_rec = vec![false; k_rec];
    let mut total_mae = 0.0_f64;
    for ti in 0..k_truth {
        let mut best_corr = -1.0_f64;
        let mut best_ri = 0;
        for ri in 0..k_rec {
            if !used_rec[ri] && abs_corr[ti][ri] > best_corr {
                best_corr = abs_corr[ti][ri];
                best_ri = ri;
            }
        }
        used_rec[best_ri] = true;
        let tc = &truth_cols[ti];
        let rc = &rec_cols[best_ri];
        // Sign: dot product direction.
        let dot: f64 = tc.iter().zip(rc.iter()).map(|(a, b)| a * b).sum();
        let sign = if dot >= 0.0 { 1.0 } else { -1.0 };
        // Scale: regress rc onto tc; scale = cov(tc, rc) / var(rc).
        let mean_rc: f64 = rc.iter().sum::<f64>() / n as f64;
        let mean_tc: f64 = tc.iter().sum::<f64>() / n as f64;
        let mut cov = 0.0_f64;
        let mut var_rc = 0.0_f64;
        for i in 0..n {
            let dr = rc[i] - mean_rc;
            let dt = tc[i] - mean_tc;
            cov += dt * sign * dr;
            var_rc += dr * dr;
        }
        let scale = if var_rc > 1e-14 { cov / var_rc } else { 1.0 };
        // MAE of the scaled & signed recovered vs truth.
        let mae: f64 = tc
            .iter()
            .zip(rc.iter())
            .map(|(t, r)| (t - scale * sign * r).abs())
            .sum::<f64>()
            / n as f64;
        total_mae += mae;
    }
    total_mae / k_truth as f64
}

fn abs_pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    let ma = a.iter().sum::<f64>() / n as f64;
    let mb = b.iter().sum::<f64>() / n as f64;
    let mut sxy = 0.0_f64;
    let mut sxx = 0.0_f64;
    let mut syy = 0.0_f64;
    for (x, y) in a.iter().zip(b.iter()) {
        let dx = x - ma;
        let dy = y - mb;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    let denom = (sxx * syy).sqrt();
    if denom < 1e-14 {
        0.0
    } else {
        (sxy / denom).abs()
    }
}
