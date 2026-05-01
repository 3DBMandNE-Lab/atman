//! End-to-end smoke test: Fig 3 demo-paper chain.
//!
//! Chain: NMF deconvolution → admixture-adjusted DE → GSEA → null calibration.
//! Each step's output feeds the next step's input. The load-bearing assertion
//! is chained-input-hash provenance in each `*.run.json` sidecar.
//!
//! Design (30 samples × 50 genes, 15 mesenchymal vs 15 other):
//!   1. `atman decompose nmf` → nmf_loadings.tsv + nmf_activations.tsv
//!   2. Pivot activations long→wide → nmf_activations_wide.tsv
//!   3. `atman de --test welch-t --adjust-for nmf_activations_wide.tsv`
//!      → de_output/de_results.tsv
//!   4. `atman enrich gsea` → gsea_results.tsv
//!   5. `atman null --test welch-t` → null_output/null_summary.tsv

use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use std::process::Command;

fn atman() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_atman"))
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(atman())
        .args(args)
        .output()
        .expect("spawn atman")
}

/// Assert a command succeeded, printing both streams on failure.
fn assert_success(output: &std::process::Output, label: &str) {
    assert!(
        output.status.success(),
        "{label} failed (exit {:?})\nstdout:\n{}\nstderr:\n{}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// Read a TSV (skip `#`-prefixed comment lines) into rows-as-hashmaps.
fn read_tsv(path: &Path) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {:?}: {e}", path));
    let mut lines = text.lines().filter(|l| !l.starts_with('#')).peekable();
    let header: Vec<&str> = match lines.next() {
        Some(h) => h.split('\t').collect(),
        None => return vec![],
    };
    let mut rows = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let mut m = HashMap::new();
        for (h, c) in header.iter().zip(cells.iter()) {
            m.insert(h.to_string(), c.to_string());
        }
        rows.push(m);
    }
    rows
}

/// Compute SHA-256 hex of a file — used to cross-check sidecar hashes.
fn sha256_of_file(path: &Path) -> String {
    use std::io::Read;
    let mut bytes = Vec::new();
    std::fs::File::open(path)
        .unwrap_or_else(|e| panic!("open {:?}: {e}", path))
        .read_to_end(&mut bytes)
        .unwrap();
    // Use the same approach atman uses: SHA-256 via sha2.
    // We don't have sha2 here as a direct dep, so we shell out.
    let output = Command::new("shasum")
        .args(["-a", "256", path.to_str().unwrap()])
        .output()
        .expect("shasum");
    let out_str = String::from_utf8_lossy(&output.stdout);
    out_str.split_whitespace().next().unwrap_or("").to_string()
}

/// Parse a run sidecar JSON and return the `inputs_sha256` object
/// as a map from key → value string (the "sha256:<hex>" strings).
fn parse_sidecar_inputs(sidecar_path: &Path) -> HashMap<String, String> {
    let text = std::fs::read_to_string(sidecar_path)
        .unwrap_or_else(|e| panic!("read sidecar {:?}: {e}", sidecar_path));
    // Parse the `inputs_sha256` block manually — no serde_json dep in tests.
    // Strategy: find `"inputs_sha256"` then extract the next `{...}` block.
    let marker = "\"inputs_sha256\"";
    let start = text
        .find(marker)
        .unwrap_or_else(|| panic!("no inputs_sha256 in {:?}", sidecar_path));
    let after = &text[start + marker.len()..];
    // Skip whitespace and `:`.
    let brace_start = after
        .find('{')
        .unwrap_or_else(|| panic!("no {{ after inputs_sha256 in {:?}", sidecar_path));
    let obj_str = &after[brace_start..];
    // Find matching close brace (depth-1 search; values are all strings so
    // no nested objects expected).
    let mut depth = 0usize;
    let mut end = 0usize;
    for (i, c) in obj_str.char_indices() {
        match c {
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 {
                    end = i + 1;
                    break;
                }
            }
            _ => {}
        }
    }
    let obj_body = &obj_str[1..end - 1]; // strip outer { }
    // Parse key:"value" pairs.
    let mut map = HashMap::new();
    // Split on `","` to get raw pairs, then parse "key": "value".
    // We iterate over comma-separated entries. Simple state machine.
    let mut remaining = obj_body.trim();
    loop {
        remaining = remaining.trim();
        if remaining.is_empty() {
            break;
        }
        // Expect: "key": "value"
        if !remaining.starts_with('"') {
            break;
        }
        let key_end = remaining[1..].find('"').expect("key close quote") + 1;
        let key = remaining[1..key_end].to_string();
        remaining = remaining[key_end + 1..].trim();
        // Expect: `: "value"`
        assert!(remaining.starts_with(':'), "expected colon after key {key:?}");
        remaining = remaining[1..].trim();
        assert!(remaining.starts_with('"'), "expected string value for key {key:?}");
        let val_end = remaining[1..].find('"').expect("value close quote") + 1;
        let val = remaining[1..val_end].to_string();
        remaining = &remaining[val_end + 1..];
        map.insert(key, val);
        // Skip past optional comma.
        remaining = remaining.trim();
        if remaining.starts_with(',') {
            remaining = &remaining[1..];
        }
    }
    map
}

/// Build a synthetic GBM-shaped canonical input directory.
///
/// 30 samples (15 mesenchymal, 15 other) × 50 genes.
/// Abundances are non-negative (raw values — atman's measurements schema
/// stores log-scale values in practice but NMF only requires non-negative
/// input; using raw positive floats satisfies both).
///
/// Returns the path to the canonical input dir.
fn build_synthetic_fixture(tmp: &Path) -> std::path::PathBuf {
    use std::fmt::Write as FmtWrite;

    let dir = tmp.join("canonical");
    std::fs::create_dir_all(&dir).unwrap();

    let n_samples = 30usize;
    let n_genes = 50usize;
    let n_mes = 15usize; // mesenchymal

    // ── samples.tsv ─────────────────────────────────────────────────────────
    let mut samples_tsv = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n",
    );
    for i in 0..n_samples {
        let sid = format!("S{:03}", i + 1);
        let cond = if i < n_mes { "mesenchymal" } else { "other" };
        writeln!(
            samples_tsv,
            "{sid}\t{sid}\t{cond}\t0\ttumour\t{}",
            i + 1
        )
        .unwrap();
    }
    std::fs::write(dir.join("samples.tsv"), &samples_tsv).unwrap();

    // ── proteins.tsv ────────────────────────────────────────────────────────
    let mut proteins_tsv =
        String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for g in 0..n_genes {
        let gene = format!("GENE{:03}", g + 1);
        let assay = format!("A{:04}", g + 1);
        writeln!(
            proteins_tsv,
            "custom\t{assay}\t\t{gene}\tgbm_panel\t"
        )
        .unwrap();
    }
    std::fs::write(dir.join("proteins.tsv"), &proteins_tsv).unwrap();

    // ── measurements.tsv ────────────────────────────────────────────────────
    // Non-negative gamma-mixed abundances: each gene has a base level + a
    // condition-specific offset so DE has something to find.
    // Seed with a simple LCG for determinism (no external RNG dep).
    let mut lcg_state: u64 = 0xdeadbeef_cafebabe;
    let mut lcg = move || -> f64 {
        lcg_state = lcg_state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        // u64 in [0, 2^64) → (0, 1)
        (lcg_state >> 11) as f64 / (1u64 << 53) as f64
    };

    let header = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n";
    let mut meas = String::from(header);
    let mut order = 1usize;
    for i in 0..n_samples {
        let sid = format!("S{:03}", i + 1);
        let is_mes = i < n_mes;
        for g in 0..n_genes {
            let gene = format!("GENE{:03}", g + 1);
            let assay = format!("A{:04}", g + 1);
            // Base level + small noise (all positive).
            let base = 8.0 + (g as f64) * 0.1;
            // Condition offset: first 10 genes elevated in mesenchymal.
            let cond_offset = if g < 10 && is_mes { 2.0 } else { 0.0 };
            let noise = lcg() * 0.5; // [0, 0.5)
            let abundance = base + cond_offset + noise;
            writeln!(
                meas,
                "custom\t{sid}\t{assay}\t{gene}\tgbm_panel\t{abundance:.6}\t{abundance:.6}\t{abundance:.6}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}"
            )
            .unwrap();
            order += 1;
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &meas).unwrap();

    dir
}

/// Build a small gene-set TSV exercising GSEA across the fixture gene universe.
/// Three sets: one loaded with the first 8 genes (should be enriched in
/// mesenchymal), one with the last 8 genes, one scattered set.
fn build_gene_sets(tmp: &Path) {
    let path = tmp.join("gene_sets.tsv");
    let mut s = String::from("set_name\tgene_symbol\n");
    for g in 1..=8 {
        s.push_str(&format!("MES_PROGRAM\tGENE{g:03}\n"));
    }
    for g in 43..=50 {
        s.push_str(&format!("OTHER_PROGRAM\tGENE{g:03}\n"));
    }
    for g in (10..=40).step_by(5) {
        s.push_str(&format!("SCATTERED\tGENE{g:03}\n"));
    }
    std::fs::write(path, s).unwrap();
}

/// Pivot NMF activations from long (sample_id, program, activation) to wide
/// (sample_id, program_0, program_1, ...) and write as a TSV file.
/// Returns the path written.
fn pivot_activations(activations_path: &Path, out_path: &Path) {
    let rows = read_tsv(activations_path);

    // Collect programs in sorted order.
    let mut programs: Vec<String> = rows
        .iter()
        .map(|r| r["program"].clone())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    programs.sort();

    // Build sample → program → activation map.
    let mut map: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    for r in &rows {
        let sid = r["sample_id"].clone();
        let prog = r["program"].clone();
        let act: f64 = r["activation"].parse().unwrap_or(0.0);
        map.entry(sid).or_default().insert(prog, act);
    }

    // Write wide format.
    let mut s = String::from("sample_id");
    for p in &programs {
        s.push('\t');
        s.push_str(p);
    }
    s.push('\n');
    for (sid, prog_map) in &map {
        s.push_str(sid);
        for p in &programs {
            s.push('\t');
            let v = prog_map.get(p).copied().unwrap_or(0.0);
            s.push_str(&format!("{v:.8}"));
        }
        s.push('\n');
    }
    std::fs::write(out_path, s).unwrap();
}

#[test]
fn e2e_demo_paper_chain() {
    let tmp = tempfile::tempdir().unwrap();
    let tmp = tmp.path();

    // ── Build synthetic fixture ──────────────────────────────────────────────
    let canonical_dir = build_synthetic_fixture(tmp);
    build_gene_sets(tmp);

    let nmf_loadings = tmp.join("nmf_loadings.tsv");
    let nmf_activations = tmp.join("nmf_activations.tsv");
    let nmf_stability = tmp.join("nmf_stability.tsv");
    let nmf_activations_wide = tmp.join("nmf_activations_wide.tsv");
    let de_output_dir = tmp.join("de_output");
    let gsea_output = tmp.join("gsea_results.tsv");
    let null_output_dir = tmp.join("null_output");
    let gene_sets = tmp.join("gene_sets.tsv");

    // ── Step 1: NMF decomposition ────────────────────────────────────────────
    let nmf_out = run(&[
        "decompose",
        "nmf",
        "--input-dir",
        canonical_dir.to_str().unwrap(),
        "--k",
        "3",
        "--n-seeds",
        "5",
        "--beta-loss",
        "frobenius",
        "--init",
        "nndsvda",
        "--max-iter",
        "200",
        "--tol",
        "1e-6",
        "--seed-base",
        "42",
        "--output-loadings",
        nmf_loadings.to_str().unwrap(),
        "--output-activations",
        nmf_activations.to_str().unwrap(),
        "--output-stability",
        nmf_stability.to_str().unwrap(),
    ]);
    assert_success(&nmf_out, "Step 1: decompose nmf");

    // Assert outputs exist and are non-empty.
    assert!(nmf_loadings.exists(), "nmf_loadings.tsv missing");
    assert!(nmf_activations.exists(), "nmf_activations.tsv missing");
    let loadings_rows = read_tsv(&nmf_loadings);
    assert!(!loadings_rows.is_empty(), "nmf_loadings.tsv has no rows");
    assert!(
        loadings_rows[0].contains_key("program"),
        "nmf_loadings.tsv missing 'program' column"
    );
    let activations_rows = read_tsv(&nmf_activations);
    assert!(!activations_rows.is_empty(), "nmf_activations.tsv has no rows");
    assert!(
        activations_rows[0].contains_key("sample_id"),
        "nmf_activations.tsv missing 'sample_id' column"
    );

    // Sidecar for NMF (keyed on loadings output path).
    let nmf_sidecar_path = {
        let p = nmf_loadings.clone();
        let mut s = p.as_os_str().to_owned();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    assert!(nmf_sidecar_path.exists(), "NMF sidecar missing: {nmf_sidecar_path:?}");
    let nmf_sidecar_text = std::fs::read_to_string(&nmf_sidecar_path).unwrap();
    assert!(
        nmf_sidecar_text.contains("\"command\""),
        "NMF sidecar missing 'command' field"
    );
    assert!(
        nmf_sidecar_text.contains("decompose nmf"),
        "NMF sidecar command is not 'decompose nmf'"
    );

    // ── Step 2: Pivot activations long → wide ────────────────────────────────
    // `--adjust-for` expects wide format: sample_id | program_0 | program_1 | ...
    pivot_activations(&nmf_activations, &nmf_activations_wide);
    assert!(
        nmf_activations_wide.exists(),
        "nmf_activations_wide.tsv was not created"
    );

    // ── Step 3: Admixture-adjusted DE (welch-t + --adjust-for) ───────────────
    // Condition comes from samples.tsv `condition` column: "mesenchymal" vs "other".
    // --adjust-for takes the wide activations as external covariates.
    let de_out = run(&[
        "de",
        "--test",
        "welch-t",
        "--input-dir",
        canonical_dir.to_str().unwrap(),
        "--output-dir",
        de_output_dir.to_str().unwrap(),
        "--groups",
        "mesenchymal-other",
        "--adjust-for",
        nmf_activations_wide.to_str().unwrap(),
        "--min-pairs",
        "5",
    ]);
    assert_success(&de_out, "Step 3: de --test welch-t --adjust-for");

    let de_results_path = de_output_dir.join("de_results.tsv");
    assert!(de_results_path.exists(), "de_results.tsv missing");
    let de_rows = read_tsv(&de_results_path);
    assert!(!de_rows.is_empty(), "de_results.tsv has no rows");
    assert!(
        de_rows[0].contains_key("gene_symbol"),
        "de_results.tsv missing 'gene_symbol' column"
    );
    assert!(
        de_rows[0].contains_key("p_value"),
        "de_results.tsv missing 'p_value' column"
    );

    // Sidecar for DE (keyed on de_results.tsv).
    let de_sidecar_path = {
        let mut s = de_results_path.as_os_str().to_owned();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    assert!(de_sidecar_path.exists(), "DE sidecar missing: {de_sidecar_path:?}");
    let de_sidecar_text = std::fs::read_to_string(&de_sidecar_path).unwrap();
    assert!(
        de_sidecar_text.contains("\"command\""),
        "DE sidecar missing 'command' field"
    );

    // ── Load-bearing assertion A: DE sidecar records nmf_activations_wide.tsv hash ──
    // The basename of the --adjust-for file is used as the key in inputs_sha256.
    let de_inputs = parse_sidecar_inputs(&de_sidecar_path);
    assert!(
        de_inputs.contains_key("nmf_activations_wide.tsv"),
        "DE sidecar inputs_sha256 is missing 'nmf_activations_wide.tsv' key.\n\
         This key is required for the auditability-under-loose-coupling claim.\n\
         Got keys: {:?}",
        de_inputs.keys().collect::<Vec<_>>()
    );

    // Cross-check: the hash in the sidecar must match the actual file.
    let expected_hash = sha256_of_file(&nmf_activations_wide);
    let sidecar_hash = de_inputs["nmf_activations_wide.tsv"]
        .strip_prefix("sha256:")
        .unwrap_or(&de_inputs["nmf_activations_wide.tsv"]);
    assert_eq!(
        sidecar_hash, expected_hash,
        "DE sidecar hash for nmf_activations_wide.tsv does not match actual file SHA-256.\n\
         sidecar={sidecar_hash}\n  actual={expected_hash}"
    );

    // ── Step 4: GSEA ────────────────────────────────────────────────────────
    let gsea_out = run(&[
        "enrich",
        "gsea",
        "--de-results",
        de_results_path.to_str().unwrap(),
        "--gene-sets",
        gene_sets.to_str().unwrap(),
        "--output",
        gsea_output.to_str().unwrap(),
        "--n-permutations",
        "500",
        "--seed",
        "42",
        "--min-set-size",
        "5",
        "--comparison",
        "mesenchymal-other",
    ]);
    assert_success(&gsea_out, "Step 4: enrich gsea");

    assert!(gsea_output.exists(), "gsea_results.tsv missing");
    let gsea_rows = read_tsv(&gsea_output);
    assert!(!gsea_rows.is_empty(), "gsea_results.tsv has no rows");
    assert!(
        gsea_rows[0].contains_key("set_name"),
        "gsea_results.tsv missing 'set_name' column"
    );
    assert!(
        gsea_rows[0].contains_key("p_value"),
        "gsea_results.tsv missing 'p_value' column"
    );
    assert!(
        gsea_rows[0].contains_key("nes"),
        "gsea_results.tsv missing 'nes' column"
    );

    // Sidecar for GSEA.
    let gsea_sidecar_path = {
        let mut s = gsea_output.as_os_str().to_owned();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    assert!(gsea_sidecar_path.exists(), "GSEA sidecar missing: {gsea_sidecar_path:?}");
    let gsea_sidecar_text = std::fs::read_to_string(&gsea_sidecar_path).unwrap();
    assert!(
        gsea_sidecar_text.contains("enrich gsea"),
        "GSEA sidecar command field is unexpected"
    );

    // ── Load-bearing assertion B: GSEA sidecar records de_results.tsv hash ──
    let gsea_inputs = parse_sidecar_inputs(&gsea_sidecar_path);
    assert!(
        gsea_inputs.contains_key("de_results"),
        "GSEA sidecar inputs_sha256 is missing 'de_results' key.\n\
         Got keys: {:?}",
        gsea_inputs.keys().collect::<Vec<_>>()
    );
    let expected_de_hash = sha256_of_file(&de_results_path);
    let gsea_de_hash = gsea_inputs["de_results"]
        .strip_prefix("sha256:")
        .unwrap_or(&gsea_inputs["de_results"]);
    assert_eq!(
        gsea_de_hash, expected_de_hash,
        "GSEA sidecar hash for de_results does not match actual de_results.tsv SHA-256.\n\
         sidecar={gsea_de_hash}\n  actual={expected_de_hash}"
    );

    // ── Step 5: Null calibration ─────────────────────────────────────────────
    let null_out = run(&[
        "null",
        "--test",
        "welch-t",
        "--input-dir",
        canonical_dir.to_str().unwrap(),
        "--output-dir",
        null_output_dir.to_str().unwrap(),
        "--groups",
        "mesenchymal-other",
        "--n",
        "200",
        "--seed",
        "42",
        "--min-pairs",
        "5",
    ]);
    assert_success(&null_out, "Step 5: null --test welch-t");

    let null_summary = null_output_dir.join("null_summary.tsv");
    let empirical_p = null_output_dir.join("empirical_p.tsv");
    assert!(null_summary.exists(), "null_summary.tsv missing");
    assert!(empirical_p.exists(), "empirical_p.tsv missing");
    let null_rows = read_tsv(&null_summary);
    assert!(!null_rows.is_empty(), "null_summary.tsv has no rows");
    let emp_rows = read_tsv(&empirical_p);
    assert!(!emp_rows.is_empty(), "empirical_p.tsv has no rows");

    // Sidecar for null (keyed on null_summary.tsv).
    let null_sidecar_path = {
        let mut s = null_summary.as_os_str().to_owned();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    assert!(null_sidecar_path.exists(), "null sidecar missing: {null_sidecar_path:?}");
    let null_sidecar_text = std::fs::read_to_string(&null_sidecar_path).unwrap();
    assert!(
        null_sidecar_text.contains("\"command\""),
        "null sidecar missing 'command' field"
    );

    // ── Load-bearing assertion C: null sidecar records canonical input hashes ──
    let null_inputs = parse_sidecar_inputs(&null_sidecar_path);
    assert!(
        null_inputs.contains_key("measurements.tsv"),
        "null sidecar inputs_sha256 is missing 'measurements.tsv' key.\n\
         Got keys: {:?}",
        null_inputs.keys().collect::<Vec<_>>()
    );
    let expected_meas_hash = sha256_of_file(&canonical_dir.join("measurements.tsv"));
    let null_meas_hash = null_inputs["measurements.tsv"]
        .strip_prefix("sha256:")
        .unwrap_or(&null_inputs["measurements.tsv"]);
    assert_eq!(
        null_meas_hash, expected_meas_hash,
        "null sidecar hash for measurements.tsv does not match actual file SHA-256.\n\
         sidecar={null_meas_hash}\n  actual={expected_meas_hash}"
    );

    // ── Summary diagnostic ───────────────────────────────────────────────────
    eprintln!(
        "\ne2e_demo_paper_chain PASSED:\n\
         NMF:  {} loadings rows, {} activations rows\n\
         DE:   {} gene rows (welch-t + adjust-for NMF programs)\n\
         GSEA: {} set rows\n\
         Null: {} summary rows, {} empirical rows",
        loadings_rows.len(),
        activations_rows.len(),
        de_rows.len(),
        gsea_rows.len(),
        null_rows.len(),
        emp_rows.len(),
    );
}
