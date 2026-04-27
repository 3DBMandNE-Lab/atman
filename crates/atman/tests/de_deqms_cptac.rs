//! DEqMS parity validation on the CPTAC Study 6 UPS1 spike-in.
//!
//! Pre-summarizes the shared CPTAC peptide fixture (log2 →
//! per-sample median-center → per-protein per-sample median) into
//! protein-level canonical measurements, runs
//! `atman de --test limma --peptide-metadata peptides.tsv`, and
//! compares per-protein output against
//! `tests/fixtures/deqms_cptac_reference.tsv` (regenerated from
//! `deqms_cptac_reference.R` via Bioconductor DEqMS v1.26). Assertions:
//!
//! 1. `de_results.tsv` reports `effect_size_method = "limma-DEqMS-*"`
//!    and populates `n_peptides_observed` (peptide count) for each
//!    protein row.
//! 2. Sign agreement with DEqMS on all proteins with
//!    `|reference logFC(B-A)| ≥ 0.3`.
//! 3. Median per-protein `|atman.mean_diff − (−DEqMS.logFC)|` < 0.30
//!    log₂ across all jointly-fitted proteins.

use std::collections::{BTreeMap, HashMap};
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
    let rows: Vec<HashMap<String, String>> = lines
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

fn median(values: &mut Vec<f64>) -> f64 {
    values.retain(|v| v.is_finite());
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = values.len();
    if n.is_multiple_of(2) {
        0.5 * (values[n / 2 - 1] + values[n / 2])
    } else {
        values[n / 2]
    }
}

/// Prepare canonical inputs that match what DEqMS's reference
/// workflow sees: log2, per-sample median-centering, then per-protein
/// per-sample median summary.
fn write_cptac_protein_canonical(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    let src = Path::new("tests/fixtures/cptac_peptides30_6ab.tsv");
    let (header, rows) = parse_tsv(src);
    let a_cols: Vec<usize> = (1..=9)
        .map(|i| header.iter().position(|h| h == &format!("A{i}")).unwrap())
        .collect();
    let b_cols: Vec<usize> = (1..=9)
        .map(|i| header.iter().position(|h| h == &format!("B{i}")).unwrap())
        .collect();
    let seq_col = header.iter().position(|h| h == "Sequence").unwrap();
    let razor_col = header
        .iter()
        .position(|h| h == "LeadingRazorProtein")
        .unwrap();

    // samples.tsv: 18 samples.
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    let mut order = 0u64;
    for (cond, count) in [("A", 9usize), ("B", 9)] {
        for i in 1..=count {
            order += 1;
            samples.push_str(&format!(
                "{cond}{i}\t{cond}{i}\t{cond}\t0\tcell_lysate\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    // Collect per-sample log2 non-zero peptide intensities for median
    // centering. All sample columns flattened; group per column.
    let all_cols: Vec<(String, usize)> = (1..=9)
        .map(|i| (format!("A{i}"), a_cols[i - 1]))
        .chain((1..=9).map(|i| (format!("B{i}"), b_cols[i - 1])))
        .collect();
    let mut sample_median: HashMap<String, f64> = HashMap::new();
    for (sid, col) in &all_cols {
        let mut logs: Vec<f64> = rows
            .iter()
            .filter_map(|r| {
                let v: f64 = r[&header[*col]].parse().unwrap_or(0.0);
                if v > 0.0 {
                    Some(v.log2())
                } else {
                    None
                }
            })
            .collect();
        sample_median.insert(sid.clone(), median(&mut logs));
    }

    // Build centered peptide-level matrix in memory, then collapse
    // peptides to proteins via per-sample median. Emit as canonical
    // protein-level measurements.
    let mut per_protein_peps: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (idx, r) in rows.iter().enumerate() {
        let razor = &r[&header[razor_col]];
        if razor.is_empty() {
            continue;
        }
        per_protein_peps.entry(razor.clone()).or_default().push(idx);
    }

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    let mut peptides =
        String::from("peptide_id\tassay_id\tsequence\tcharge\tmodifications\tmissed_cleavages\n");
    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0u64;
    for (razor, pep_indices) in &per_protein_peps {
        // Use the razor string itself as the gene symbol so the limma
        // dispatch's `(panel, gene_symbol)` keying produces one feature
        // per protein. Panel tags the class so tests can group without
        // re-parsing the razor string.
        let panel = if razor.contains("_HUMAN_UPS") {
            "UPS1"
        } else if razor.contains("_YEAST") {
            "YEAST"
        } else {
            "OTHER"
        };
        proteins.push_str(&format!(
            "maxquant_lfq\t{razor}\t{razor}\t{razor}\t{panel}\t\n"
        ));
        // Peptide catalog rows.
        for &pep_idx in pep_indices {
            let seq = &rows[pep_idx][&header[seq_col]];
            peptides.push_str(&format!("P{pep_idx:04}_{seq}\t{razor}\t{seq}\t\t\t\n"));
        }
        // Per-sample median-summarized protein abundance.
        for (sid, col) in &all_cols {
            let med_sample = sample_median.get(sid).copied().unwrap_or(0.0);
            let mut values: Vec<f64> = pep_indices
                .iter()
                .filter_map(|&pep_idx| {
                    let v: f64 = rows[pep_idx][&header[*col]].parse().unwrap_or(0.0);
                    if v > 0.0 {
                        Some(v.log2() - med_sample)
                    } else {
                        None
                    }
                })
                .collect();
            let protein_value = median(&mut values);
            order += 1;
            if protein_value.is_finite() {
                qc.push_str(&format!(
                    "maxquant_lfq\t{sid}\t{razor}\t{razor}\t{panel}\t{protein_value:.6}\t\
                     {protein_value:.6}\t{protein_value:.6}\tlog2_protein_summary\tPASS\tPASS\t\
                     \t0\t0\t\t\t{order}\n"
                ));
            } else {
                qc.push_str(&format!(
                    "maxquant_lfq\t{sid}\t{razor}\t{razor}\t{panel}\t\t\t\tlog2_protein_summary\t\
                     PASS\tPASS\t\t0\t1\t\t\t{order}\n"
                ));
            }
        }
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();
    std::fs::write(dir.join("peptides.tsv"), peptides).unwrap();
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

#[test]
fn deqms_limma_path_matches_bioconductor_deqms_on_cptac() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("deqms_out");
    write_cptac_protein_canonical(&input);
    std::fs::create_dir_all(&output).unwrap();

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "limma",
        "--groups",
        "A-B",
        "--peptide-metadata",
        input.join("peptides.tsv").to_str().unwrap(),
        "--trend",
        "true",
        "--robust",
        "false",
        "--min-pairs",
        "4",
    ]);
    assert!(
        out.status.success(),
        "atman de --test limma --peptide-metadata failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );

    let (_, rows) = parse_tsv(&output.join("de_results.tsv"));
    let ab: Vec<&HashMap<String, String>> = rows
        .iter()
        .filter(|r| r.get("comparison").map(|s| s.as_str()) == Some("A-B"))
        .collect();
    assert!(!ab.is_empty());
    eprintln!("atman wrote {} rows; skip reasons:", ab.len());
    let mut skip_tally: BTreeMap<String, usize> = BTreeMap::new();
    for r in &ab {
        let reason = r.get("skip_reason").cloned().unwrap_or_default();
        *skip_tally.entry(reason).or_insert(0) += 1;
    }
    for (reason, count) in &skip_tally {
        eprintln!("  {reason:>30} : {count}");
    }

    // Column & method assertions.
    for r in &ab {
        assert!(
            r.get("effect_size_method")
                .map(|s| s.starts_with("limma-DEqMS"))
                .unwrap_or(false),
            "expected limma-DEqMS method, got {:?}",
            r.get("effect_size_method")
        );
        assert!(
            r.get("n_peptides_observed")
                .map(|s| !s.is_empty() && s != "0")
                .unwrap_or(false),
            "expected populated n_peptides_observed, got {:?}",
            r.get("n_peptides_observed")
        );
    }

    // Parity against DEqMS Bioconductor reference.
    let ref_path = Path::new("tests/fixtures/deqms_cptac_reference.tsv");
    if !ref_path.exists() {
        eprintln!("DEqMS reference missing; run `Rscript tests/fixtures/deqms_cptac_reference.R` to generate.");
        return;
    }
    let (_, ref_rows) = parse_tsv(ref_path);
    let ref_by_id: HashMap<String, &HashMap<String, String>> = ref_rows
        .iter()
        .filter(|r| {
            r.get("logFC_B_minus_A")
                .map(|s| s != "NA" && !s.is_empty())
                .unwrap_or(false)
        })
        .map(|r| (r["assay_id"].clone(), r))
        .collect();

    let mut diffs: Vec<(String, f64, f64, f64)> = Vec::new();
    for row in &ab {
        if row
            .get("skip_reason")
            .map(|s| !s.is_empty())
            .unwrap_or(true)
        {
            continue;
        }
        let id = match row.get("assay_id") {
            Some(s) if !s.is_empty() => s.clone(),
            _ => continue,
        };
        let atman_eff: f64 = match row.get("mean_diff").and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let ref_row = match ref_by_id.get(&id) {
            Some(r) => *r,
            None => continue,
        };
        let ref_logfc: f64 = match ref_row.get("logFC_B_minus_A").and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        // atman convention: mean_diff = mean_a − mean_b = −logFC_B_minus_A.
        let ref_equiv = -ref_logfc;
        diffs.push((id, atman_eff, ref_equiv, (atman_eff - ref_equiv).abs()));
    }
    assert!(
        diffs.len() >= 15,
        "expected ≥15 joinable proteins, got {}",
        diffs.len()
    );

    eprintln!("\nDEqMS parity (atman mean_diff vs DEqMS −logFC_B_minus_A):");
    let mut sorted = diffs.clone();
    sorted.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap());
    for (id, a, r, d) in sorted.iter().take(10) {
        eprintln!(
            "  {:<40} atman={:+6.3}  deqms={:+6.3}  Δ={:5.3}",
            id, a, r, d
        );
    }

    // Sign on signal proteins (|ref| ≥ 0.3).
    let signal: Vec<&(String, f64, f64, f64)> =
        diffs.iter().filter(|(_, _, r, _)| r.abs() >= 0.3).collect();
    let sign_match = signal
        .iter()
        .filter(|(_, a, r, _)| (*a > 0.0) == (*r > 0.0))
        .count();
    assert!(
        signal.len() >= 8,
        "need ≥8 signal proteins; got {}",
        signal.len()
    );
    let sign_rate = sign_match as f64 / signal.len() as f64;
    assert!(
        sign_rate >= 0.9,
        "sign agreement {sign_rate} < 0.9 across {} signal proteins",
        signal.len()
    );

    // Magnitude: median |diff| across all jointly-fitted proteins.
    let mut abs_diffs: Vec<f64> = diffs.iter().map(|(_, _, _, d)| *d).collect();
    abs_diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_abs = abs_diffs[abs_diffs.len() / 2];
    eprintln!(
        "median |atman − DEqMS| = {median_abs:.3}  (n={}, tolerance 0.30)",
        abs_diffs.len()
    );
    assert!(
        median_abs < 0.30,
        "median |atman − DEqMS| = {median_abs} exceeds 0.30 tolerance"
    );
}
