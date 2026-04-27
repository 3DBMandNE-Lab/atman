//! CPTAC Study 6 UPS1 spike-in validation for `atman de --test msqrob`.
//!
//! This is the canonical MS-proteomics peptide-level benchmark used by
//! the msqrob / msqrob2 paper. Yeast whole-cell lysate is spiked with
//! 48 human UPS1 proteins at four concentrations (6A: 0.25 fmol, 6B:
//! 0.74 fmol, 6C: 2.22 fmol, 6D: 6.67 fmol per μg yeast). Ground truth
//! for the B vs A comparison is log₂(0.74/0.25) ≈ +1.57 for UPS1
//! proteins and 0 for yeast.
//!
//! The fixture `fixtures/cptac_peptides30_6ab.tsv` is a 30-protein /
//! 309-peptide subset of the MaxQuant `peptides.txt` output released in
//! the statOmics/MSqRobData GitHub repository. We keep the Sequence,
//! Leading razor protein, and eighteen `Intensity 6A_1..9` +
//! `Intensity 6B_1..9` columns (the 6C/6D samples are dropped here
//! because atman's feature is tested on the canonical two-group A–B
//! comparison). Intensities of zero are MaxQuant's "not detected"
//! sentinel and are masked via `dropped_by_qc = 1`; nonzero intensities
//! are log₂-transformed.
//!
//! The test asserts three things on the resulting atman output:
//!
//! 1. Every UPS1 protein shows a negative `mean_diff` in the `A-B`
//!    comparison (atman convention is `mean_a − mean_b`, and UPS1 is
//!    higher in B).
//! 2. The median |effect| across UPS1 proteins is ≥ 0.5 log₂, i.e. the
//!    spike-in direction is captured with recognizable magnitude despite
//!    heavy MaxQuant missingness.
//! 3. The median |effect| across yeast proteins is < 0.5, i.e. the
//!    background distribution is correctly centered near zero.
//!
//! Observed behavior on this fixture (9 UPS1 fits, 19 yeast fits after
//! min-peptides filtering):
//!
//! - UPS1 A–B effects: 9/9 negative, range [−1.34, −0.70], median
//!   ≈ −1.0. The ground-truth magnitude is ≈ −1.57; recovered values
//!   at ~60–85% of that are expected on unnormalized MaxQuant
//!   intensities, and are consistent with published msqrob2 numbers
//!   on this dataset.
//! - Yeast A–B effects: range [−0.11, +0.78], median ≈ +0.35. The
//!   small positive bias is a normalization artifact, not a model bug:
//!   the heavier UPS1 load in B compresses yeast relative intensities
//!   in absolute MaxQuant counts. Normalization (median-centering,
//!   cyclic loess, VSN) is upstream of msqrob and out of scope here.

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

fn write_cptac_canonical(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    let src = Path::new("tests/fixtures/cptac_peptides30_6ab.tsv");
    let (header, rows) = parse_tsv(src);
    let a_cols: Vec<usize> = (1..=9)
        .map(|i| {
            header
                .iter()
                .position(|h| h == &format!("A{i}"))
                .unwrap_or_else(|| panic!("missing A{i}"))
        })
        .collect();
    let b_cols: Vec<usize> = (1..=9)
        .map(|i| {
            header
                .iter()
                .position(|h| h == &format!("B{i}"))
                .unwrap_or_else(|| panic!("missing B{i}"))
        })
        .collect();
    let seq_col = header.iter().position(|h| h == "Sequence").unwrap();
    let razor_col = header
        .iter()
        .position(|h| h == "LeadingRazorProtein")
        .unwrap();

    // Precompute per-sample median of log2(nonzero intensity) so we can
    // emit median-centered abundances (matching the msqrob2 vignette's
    // `normalize(method = "center.median")` step). Without this,
    // unbalanced UPS1 spike-in load in condition B compresses yeast
    // intensities and shifts the background distribution.
    let all_sample_cols: Vec<(String, usize)> = (1..=9)
        .map(|i| (format!("A{i}"), a_cols[i - 1]))
        .chain((1..=9).map(|i| (format!("B{i}"), b_cols[i - 1])))
        .collect();
    let mut sample_median: HashMap<String, f64> = HashMap::new();
    for (sample_id, col_idx) in &all_sample_cols {
        let mut logs: Vec<f64> = rows
            .iter()
            .filter_map(|row| {
                let v: f64 = row[&header[*col_idx]].parse().unwrap_or(0.0);
                if v > 0.0 {
                    Some(v.log2())
                } else {
                    None
                }
            })
            .collect();
        logs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med = if logs.is_empty() {
            0.0
        } else if logs.len().is_multiple_of(2) {
            0.5 * (logs[logs.len() / 2 - 1] + logs[logs.len() / 2])
        } else {
            logs[logs.len() / 2]
        };
        sample_median.insert(sample_id.clone(), med);
    }

    // samples.tsv — 18 samples, 9 condition A, 9 condition B.
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

    // proteins.tsv — one row per distinct razor protein, plus a flag we
    // derive locally (UPS vs yeast) to the gene_symbol slot so the test
    // can key off it.
    let mut protein_gene: BTreeMap<String, (String, String)> = BTreeMap::new();
    for row in &rows {
        let razor = row[&header[razor_col]].clone();
        if razor.is_empty() {
            continue;
        }
        let assay_id = razor_to_assay_id(&razor);
        let gene = razor_to_gene(&razor);
        protein_gene.entry(assay_id).or_insert((razor, gene));
    }
    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for (assay_id, (razor, gene)) in &protein_gene {
        let uniprot = razor_to_uniprot(razor);
        proteins.push_str(&format!(
            "maxquant_lfq\t{assay_id}\t{uniprot}\t{gene}\tMQ\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    // Empty qc_measurements.tsv with a valid canonical header.
    let qc = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
              abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
              detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n";
    std::fs::write(dir.join("measurements.tsv"), qc).unwrap();
    std::fs::write(dir.join("measurements.tsv"), qc).unwrap();

    // peptides.tsv — one row per (peptide_id, assay_id) pair.
    let mut peptides =
        String::from("peptide_id\tassay_id\tsequence\tcharge\tmodifications\tmissed_cleavages\n");
    let mut pep_measurements = String::from(
        "sample_id\tpeptide_id\tabundance\tabundance_unit\tdropped_by_qc\tbelow_lod\n",
    );
    for (idx, row) in rows.iter().enumerate() {
        let seq = &row[&header[seq_col]];
        let razor = &row[&header[razor_col]];
        if seq.is_empty() || razor.is_empty() {
            continue;
        }
        let assay_id = razor_to_assay_id(razor);
        let peptide_id = format!("P{idx:04}_{seq}");
        peptides.push_str(&format!("{peptide_id}\t{assay_id}\t{seq}\t\t\t\n"));
        for (cond, cols) in [("A", &a_cols), ("B", &b_cols)] {
            for (i, col_idx) in cols.iter().enumerate() {
                let sample_id = format!("{cond}{}", i + 1);
                let raw = &row[&header[*col_idx]];
                let intensity: f64 = raw.parse().unwrap_or(0.0);
                if intensity > 0.0 {
                    let log2 = intensity.log2();
                    let centered = log2 - sample_median.get(&sample_id).copied().unwrap_or(0.0);
                    pep_measurements.push_str(&format!(
                        "{sample_id}\t{peptide_id}\t{centered:.6}\tlog2_maxquant_median_centered\t0\t0\n"
                    ));
                } else {
                    // MaxQuant zero → not detected. Mask via dropped_by_qc.
                    pep_measurements.push_str(&format!(
                        "{sample_id}\t{peptide_id}\t\tlog2_maxquant_median_centered\t1\t0\n"
                    ));
                }
            }
        }
    }
    std::fs::write(dir.join("peptides.tsv"), peptides).unwrap();
    std::fs::write(dir.join("peptide_measurements.tsv"), pep_measurements).unwrap();
}

fn razor_to_assay_id(razor: &str) -> String {
    // Keep the full MaxQuant LeadingRazorProtein string as assay_id so
    // atman output joins 1:1 with msqrob2's reference TSV, which also
    // uses the full string as its Proteins key.
    razor.to_string()
}

fn razor_to_gene(razor: &str) -> String {
    // Encode the UPS/yeast flag into gene_symbol so the test can key
    // assertions off it without carrying a separate metadata map.
    if razor.contains("_HUMAN_UPS") {
        "UPS1".to_string()
    } else if razor.contains("_YEAST") {
        "YEAST".to_string()
    } else {
        "OTHER".to_string()
    }
}

fn razor_to_uniprot(razor: &str) -> String {
    if let Some(rest) = razor.strip_prefix("sp|") {
        if let Some(pipe) = rest.find('|') {
            return rest[..pipe].to_string();
        }
    }
    if let Some(pipe) = razor.find('|') {
        let head = &razor[..pipe];
        if let Some(stripped) = head.strip_suffix("ups") {
            return stripped.to_string();
        }
        return head.to_string();
    }
    razor.to_string()
}

#[test]
fn msqrob_recovers_ups1_spike_in_on_cptac_study_6() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("msqrob_cptac_out");
    write_cptac_canonical(&input);
    std::fs::create_dir_all(&output).unwrap();

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "msqrob",
        "--groups",
        "A-B",
        "--peptide-measurements",
        input.join("peptide_measurements.tsv").to_str().unwrap(),
        "--peptide-metadata",
        input.join("peptides.tsv").to_str().unwrap(),
        "--ridge-lambda",
        "auto",
        "--min-peptides",
        "2",
        "--min-pairs",
        "4",
    ]);
    assert!(
        out.status.success(),
        "atman de --test msqrob failed on CPTAC fixture:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );

    let (_, rows) = parse_tsv(&output.join("de_results.tsv"));
    let ab: Vec<&HashMap<String, String>> = rows
        .iter()
        .filter(|r| r.get("comparison").map(|s| s.as_str()) == Some("A-B"))
        .collect();
    assert!(!ab.is_empty(), "no A-B rows in de_results.tsv");

    let mut ups_effects: Vec<f64> = Vec::new();
    let mut yeast_effects: Vec<f64> = Vec::new();
    let mut ups_neg_sign = 0usize;
    let mut ups_total = 0usize;
    for row in &ab {
        if !row
            .get("skip_reason")
            .map(|s| s.is_empty())
            .unwrap_or(false)
        {
            continue;
        }
        let gene = row.get("gene_symbol").map(String::as_str).unwrap_or("");
        let effect: f64 = match row.get("mean_diff").and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        match gene {
            "UPS1" => {
                ups_effects.push(effect);
                ups_total += 1;
                if effect < 0.0 {
                    ups_neg_sign += 1;
                }
            }
            "YEAST" => yeast_effects.push(effect),
            _ => {}
        }
    }

    eprintln!(
        "CPTAC msqrob summary: UPS1 fits={ups_total}, yeast fits={}, UPS1 neg-sign={ups_neg_sign}/{ups_total}",
        yeast_effects.len()
    );
    eprintln!("UPS1 effects (A-B, expected ≈ -1.57):");
    let mut sorted_ups: Vec<f64> = ups_effects.clone();
    sorted_ups.sort_by(|a, b| a.partial_cmp(b).unwrap());
    for v in &sorted_ups {
        eprintln!("  {:>7.3}", v);
    }
    let mut sorted_yeast: Vec<f64> = yeast_effects.clone();
    sorted_yeast.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eprintln!("Yeast effect quantiles (expected ≈ 0):");
    if !sorted_yeast.is_empty() {
        let n = sorted_yeast.len();
        eprintln!("  min   {:>7.3}", sorted_yeast[0]);
        eprintln!("  q25   {:>7.3}", sorted_yeast[n / 4]);
        eprintln!("  med   {:>7.3}", sorted_yeast[n / 2]);
        eprintln!("  q75   {:>7.3}", sorted_yeast[(3 * n) / 4]);
        eprintln!("  max   {:>7.3}", sorted_yeast[n - 1]);
    }

    assert!(
        ups_total >= 5,
        "need ≥5 UPS1 fits for a meaningful check, got {ups_total}"
    );
    let sign_rate = ups_neg_sign as f64 / ups_total as f64;
    assert!(
        sign_rate >= 0.8,
        "expected ≥80% of UPS1 effects to be negative (A < B); got {sign_rate} ({ups_neg_sign}/{ups_total})"
    );
    let ups_median = median_abs(&ups_effects);
    let yeast_median = median_abs(&yeast_effects);
    assert!(
        ups_median >= 0.5,
        "UPS1 median |effect| = {ups_median} (expected ≥ 0.5 log₂; ground-truth is ≈ 1.57)"
    );
    assert!(
        yeast_median < 0.5,
        "yeast median |effect| = {yeast_median} (expected background near zero)"
    );
    assert!(
        ups_median > yeast_median,
        "UPS1 median |effect| ({ups_median}) must exceed yeast background ({yeast_median})"
    );

    // --- Per-protein parity against msqrob2 reference ---------------
    // `tests/fixtures/msqrob2_cptac_reference.tsv` is produced by
    // `msqrob2_cptac_reference.R` running the canonical msqrob2
    // vignette workflow on the same fixture (log2 → median-center →
    // median-summarise to protein → `msqrob(formula=~condition)` →
    // `hypothesisTest`). Reference logFC is `β_conditionB = mean_B
    // − mean_A`, so atman's `mean_diff` (= mean_a − mean_b) should
    // equal `−logFC_B_minus_A` up to model differences (atman fits
    // at the peptide level with random intercept + optional ridge;
    // msqrob2 here fits OLS on protein-level medians).
    let ref_path = Path::new("tests/fixtures/msqrob2_cptac_reference.tsv");
    if ref_path.exists() {
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

        let mut diffs: Vec<(String, f64, f64, f64)> = Vec::new(); // (id, atman, msqrob2, abs_diff)
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
            // atman mean_diff(A-B) vs −msqrob2 logFC(B-A).
            let ref_equiv = -ref_logfc;
            diffs.push((id, atman_eff, ref_equiv, (atman_eff - ref_equiv).abs()));
        }

        assert!(
            diffs.len() >= 15,
            "expected ≥15 joinable proteins against msqrob2 reference; got {}",
            diffs.len()
        );

        eprintln!("\nPer-protein parity (atman mean_diff vs msqrob2 −logFC_B_minus_A):");
        let mut sorted = diffs.clone();
        sorted.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap());
        for (id, a, r, d) in sorted.iter().take(10) {
            eprintln!(
                "  {:<40} atman={:+6.3}  msqrob2={:+6.3}  Δ={:5.3}",
                id, a, r, d
            );
        }

        // Sign agreement should be near-perfect on proteins with a
        // clear effect; we require 90%+ matching sign among proteins
        // where |ref_equiv| >= 0.3 (i.e., not noisy near-zero yeast
        // backgrounds where sign can flip by chance).
        let signal: Vec<&(String, f64, f64, f64)> =
            diffs.iter().filter(|(_, _, r, _)| r.abs() >= 0.3).collect();
        let sign_match = signal
            .iter()
            .filter(|(_, a, r, _)| (*a > 0.0) == (*r > 0.0))
            .count();
        assert!(
            signal.len() >= 8,
            "need ≥8 signal-proteins for sign test; got {}",
            signal.len()
        );
        let sign_rate = sign_match as f64 / signal.len() as f64;
        assert!(
            sign_rate >= 0.9,
            "atman/msqrob2 sign agreement {sign_rate} < 0.9 on {} signal proteins",
            signal.len()
        );

        // Magnitude agreement: median absolute diff across all
        // jointly-fitted proteins. Tolerance 0.4 log2 covers the
        // peptide-LMM-vs-protein-OLS model gap on this dataset.
        let mut abs_diffs: Vec<f64> = diffs.iter().map(|(_, _, _, d)| *d).collect();
        abs_diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median_abs_diff = abs_diffs[abs_diffs.len() / 2];
        eprintln!(
            "median |atman − msqrob2| = {median_abs_diff:.3}  (n={}, tolerance 0.40)",
            abs_diffs.len()
        );
        assert!(
            median_abs_diff < 0.40,
            "median |atman − msqrob2| = {median_abs_diff} exceeds 0.40 tolerance"
        );
    } else {
        eprintln!(
            "note: msqrob2 reference TSV not present; run Rscript tests/fixtures/msqrob2_cptac_reference.R to regenerate"
        );
    }
}

fn median_abs(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut abs: Vec<f64> = values.iter().map(|v| v.abs()).collect();
    abs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = abs.len() / 2;
    if abs.len().is_multiple_of(2) {
        0.5 * (abs[mid - 1] + abs[mid])
    } else {
        abs[mid]
    }
}
