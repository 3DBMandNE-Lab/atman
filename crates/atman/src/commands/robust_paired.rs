//! Leave-one-pair-out robustness diagnostic for paired DE.
//!
//! For each protein and each requested A-vs-B comparison, runs a paired t
//! on the full set of subject pairs and on every n-1 jackknife replicate
//! (one pair dropped). Reports per-protein robustness metrics that the
//! parametric paired-t alone cannot show: sign stability, p-value
//! stability, range of t, and the most-influential pair (whose deletion
//! changes |t| most).
//!
//! Pairing is read from a column in `samples.tsv` (default `patient_id`,
//! e.g. emitted by `atman recover-plex --infer-pairs`). The biological
//! replicate is the unit; a "pair" is a single (A, B) observation linked
//! by the pairing column.

use anyhow::{bail, Context, Result};
use atman_core::de::{bh_fdr, paired_t, PairedTResult};
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use std::time::SystemTime;

use crate::commands::parse_comparisons;
use crate::io::{
    atomic_write, hash_canonical_inputs, need_col, read_measurements_long, sidecar_path_for,
    write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Comparison list in A-B form. The paired difference is A - B.
    #[arg(long)]
    groups: String,

    /// Column in samples.tsv whose value identifies a pair (one A sample
    /// and one B sample share the same value). Defaults to `patient_id`.
    #[arg(long, default_value = "patient_id")]
    paired_by: String,

    /// Minimum number of complete pairs required to compute the test.
    #[arg(long, default_value_t = 5)]
    min_pairs: usize,

    /// Output TSV path for per-protein robustness rows.
    #[arg(long)]
    output: PathBuf,
}

#[derive(Debug)]
struct OutRow {
    comparison: String,
    panel: String,
    assay_id: String,
    gene_symbol: String,
    n_pairs_full: usize,
    mean_diff_full: Option<f64>,
    t_full: Option<f64>,
    p_full: Option<f64>,
    q_full: Option<f64>,
    loso_sign_stability: Option<f64>,
    loso_p_lt_05_stability: Option<f64>,
    loso_min_t: Option<f64>,
    loso_max_t: Option<f64>,
    loso_max_abs_t_delta: Option<f64>,
    most_influential_pair_id: String,
    skip_reason: String,
}

#[derive(Debug)]
struct SummaryRow {
    comparison: String,
    n_proteins: usize,
    n_proteins_full_q05: usize,
    n_proteins_robust: usize,
    fraction_robust_of_significant: f64,
    median_sign_stability_of_significant: f64,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    let comparisons = parse_comparisons(&args.groups)?;
    let samples_path = args.input_dir.join("samples.tsv");
    let pairing = read_pairing(&samples_path, &args.paired_by)?;
    let measurements_path = args.input_dir.join("measurements.tsv");
    let measurements = read_measurements_long(&measurements_path)?;

    // Index abundance by (sample_id, assay_id). We use the canonical `abundance` field;
    // detectability already covers the missingness layer separately.
    let mut by_sample_assay: HashMap<(String, String, String), f64> = HashMap::new();
    let mut assay_meta: BTreeMap<(String, String), (String, String)> = BTreeMap::new(); // (platform, assay) -> (panel, gene)
    for m in &measurements {
        let v = m.abundance.as_f64();
        if !v.is_finite() {
            continue;
        }
        if m.dropped_by_qc {
            continue;
        }
        let platform = m.platform.as_str().to_string();
        by_sample_assay.insert(
            (m.sample_id.clone(), platform.clone(), m.assay_id.0.clone()),
            v,
        );
        assay_meta
            .entry((platform, m.assay_id.0.clone()))
            .or_insert_with(|| {
                (
                    m.panel.clone().unwrap_or_default(),
                    m.gene_symbol.clone().unwrap_or_default(),
                )
            });
    }

    let mut all_rows: Vec<OutRow> = Vec::new();
    let mut summaries: Vec<SummaryRow> = Vec::new();

    for (cond_a, cond_b) in &comparisons {
        let pair_units = collect_pairs(&pairing, cond_a, cond_b)?;
        if pair_units.len() < args.min_pairs {
            bail!(
                "comparison {}-{}: only {} pairs found (need >= {})",
                cond_a,
                cond_b,
                pair_units.len(),
                args.min_pairs
            );
        }
        let assay_ids: Vec<(String, String)> = assay_meta.keys().cloned().collect();
        let mut rows: Vec<OutRow> = Vec::with_capacity(assay_ids.len());

        for assay in &assay_ids {
            let (platform, assay_id) = assay;
            let (panel, gene) = assay_meta.get(assay).cloned().unwrap_or_default();
            let mut pair_values: Vec<(f64, f64, String)> = Vec::new();
            for unit in &pair_units {
                let a = by_sample_assay
                    .get(&(unit.sample_a.clone(), platform.clone(), assay_id.clone()))
                    .copied();
                let b = by_sample_assay
                    .get(&(unit.sample_b.clone(), platform.clone(), assay_id.clone()))
                    .copied();
                if let (Some(av), Some(bv)) = (a, b) {
                    pair_values.push((av, bv, unit.pair_id.clone()));
                }
            }
            let n_full = pair_values.len();

            if n_full < args.min_pairs {
                rows.push(OutRow {
                    comparison: format!("{}-{}", cond_a, cond_b),
                    panel: panel.clone(),
                    assay_id: assay_id.clone(),
                    gene_symbol: gene.clone(),
                    n_pairs_full: n_full,
                    mean_diff_full: None,
                    t_full: None,
                    p_full: None,
                    q_full: None,
                    loso_sign_stability: None,
                    loso_p_lt_05_stability: None,
                    loso_min_t: None,
                    loso_max_t: None,
                    loso_max_abs_t_delta: None,
                    most_influential_pair_id: String::new(),
                    skip_reason: format!("insufficient_pairs (n={})", n_full),
                });
                continue;
            }

            let pairs_only: Vec<(f64, f64)> =
                pair_values.iter().map(|(a, b, _)| (*a, *b)).collect();
            let full = paired_t(&pairs_only, args.min_pairs);

            let (mean_diff_full, t_full, p_full) = match full {
                PairedTResult::Computed {
                    mean_diff,
                    t,
                    p_value,
                    ..
                } => (Some(mean_diff), Some(t), Some(p_value)),
                PairedTResult::Skipped { .. } => (None, None, None),
            };

            // LOSO replicates.
            let mut sign_keep = 0usize;
            let mut p_lt_05 = 0usize;
            let mut t_min = f64::INFINITY;
            let mut t_max = f64::NEG_INFINITY;
            let mut max_dt = 0.0f64;
            let mut influential_pair_id = String::new();
            let mut counted_loso = 0usize;
            if let (Some(t_full_v), Some(mean_diff_v)) = (t_full, mean_diff_full) {
                // Sign of the full-sample effect we are testing LOSO stability
                // against. `f64::signum(0.0)` is +1.0 in Rust, so use an explicit
                // ternary that maps an exactly-zero effect to 0.0 (no direction).
                // A zero target can never match a LOSO replicate's signum (which
                // is only ever counted when it is +1.0 or -1.0; a zero replicate
                // is likewise treated as non-matching below), so a no-direction
                // full effect honestly yields sign_stability = 0.
                let target_sign = if mean_diff_v > 0.0 {
                    1.0
                } else if mean_diff_v < 0.0 {
                    -1.0
                } else {
                    0.0
                };
                for drop_idx in 0..n_full {
                    let mut sub: Vec<(f64, f64)> = Vec::with_capacity(n_full - 1);
                    for (i, (a, b)) in pairs_only.iter().enumerate() {
                        if i != drop_idx {
                            sub.push((*a, *b));
                        }
                    }
                    if let PairedTResult::Computed {
                        mean_diff,
                        t,
                        p_value,
                        ..
                    } = paired_t(&sub, args.min_pairs)
                    {
                        counted_loso += 1;
                        // An exactly-zero replicate has no direction, so it is
                        // NOT sign-stable with any target (including a positive
                        // one): `f64::signum(0.0)` would spuriously report +1.0,
                        // so compare an explicit direction instead. A zero
                        // replicate counts as non-matching.
                        let replicate_sign = if mean_diff > 0.0 {
                            1.0
                        } else if mean_diff < 0.0 {
                            -1.0
                        } else {
                            0.0
                        };
                        if replicate_sign != 0.0 && replicate_sign == target_sign {
                            sign_keep += 1;
                        }
                        if p_value < 0.05 {
                            p_lt_05 += 1;
                        }
                        if t < t_min {
                            t_min = t;
                        }
                        if t > t_max {
                            t_max = t;
                        }
                        let dt = (t - t_full_v).abs();
                        if dt > max_dt {
                            max_dt = dt;
                            influential_pair_id = pair_values[drop_idx].2.clone();
                        }
                    }
                }
            }
            let (sign_stab, p_stab, min_t_o, max_t_o, max_dt_o) = if counted_loso == 0 {
                (None, None, None, None, None)
            } else {
                (
                    Some(sign_keep as f64 / counted_loso as f64),
                    Some(p_lt_05 as f64 / counted_loso as f64),
                    Some(t_min),
                    Some(t_max),
                    Some(max_dt),
                )
            };

            rows.push(OutRow {
                comparison: format!("{}-{}", cond_a, cond_b),
                panel,
                assay_id: assay_id.clone(),
                gene_symbol: gene,
                n_pairs_full: n_full,
                mean_diff_full,
                t_full,
                p_full,
                q_full: None,
                loso_sign_stability: sign_stab,
                loso_p_lt_05_stability: p_stab,
                loso_min_t: min_t_o,
                loso_max_t: max_t_o,
                loso_max_abs_t_delta: max_dt_o,
                most_influential_pair_id: influential_pair_id,
                skip_reason: String::new(),
            });
        }

        // BH-FDR over the full p-values, per comparison.
        let p_vec: Vec<Option<f64>> = rows.iter().map(|r| r.p_full).collect();
        let q_vec = bh_fdr(&p_vec);
        for (r, q) in rows.iter_mut().zip(q_vec.iter()) {
            r.q_full = *q;
        }

        // Summary.
        let n_proteins = rows.len();
        let n_q05 = rows
            .iter()
            .filter(|r| r.q_full.map(|q| q < 0.05).unwrap_or(false))
            .count();
        let n_robust = rows
            .iter()
            .filter(|r| {
                r.q_full.map(|q| q < 0.05).unwrap_or(false)
                    && r.loso_sign_stability.map(|s| s >= 0.999).unwrap_or(false)
                    && r.loso_p_lt_05_stability.map(|s| s >= 0.95).unwrap_or(false)
            })
            .count();
        let frac_robust = if n_q05 == 0 {
            0.0
        } else {
            n_robust as f64 / n_q05 as f64
        };
        let median_sign = median_f(
            &rows
                .iter()
                .filter(|r| r.q_full.map(|q| q < 0.05).unwrap_or(false))
                .filter_map(|r| r.loso_sign_stability)
                .collect::<Vec<_>>(),
        );
        summaries.push(SummaryRow {
            comparison: format!("{}-{}", cond_a, cond_b),
            n_proteins,
            n_proteins_full_q05: n_q05,
            n_proteins_robust: n_robust,
            fraction_robust_of_significant: frac_robust,
            median_sign_stability_of_significant: median_sign,
        });

        all_rows.extend(rows);
    }

    write_rows(&args.output, &all_rows)?;
    let summary_path = summary_path_for(&args.output);
    write_summary(&summary_path, &summaries)?;

    eprintln!(
        "robust_paired: comparisons={} proteins={} output={} summary={}",
        comparisons.len(),
        all_rows.len(),
        args.output.display(),
        summary_path.display()
    );
    for s in &summaries {
        eprintln!(
            "robust_paired: {} significant_at_q05={} robust_of_significant={} ({:.1}%)",
            s.comparison,
            s.n_proteins_full_q05,
            s.n_proteins_robust,
            100.0 * s.fraction_robust_of_significant
        );
    }

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "robust-paired",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "groups": args.groups,
            "paired-by": args.paired_by,
            "min-pairs": args.min_pairs,
            "output": args.output.display().to_string(),
            "summary-output": summary_path.display().to_string(),
        }),
        &inputs_sha256,
        &[args.output.clone(), summary_path.clone()],
        started_at,
        finished_at,
        None,
    )?;

    Ok(())
}

#[derive(Debug, Clone)]
struct PairUnit {
    pair_id: String,
    sample_a: String,
    sample_b: String,
}

fn read_pairing(
    samples_path: &std::path::Path,
    paired_by: &str,
) -> Result<Vec<(String, String, String)>> {
    // Returns Vec<(sample_id, condition, pair_value)>.
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(samples_path)
        .with_context(|| format!("opening {:?}", samples_path))?;
    let headers = reader.headers()?.clone();
    let c_sample = need_col(&headers, "sample_id", samples_path)?;
    let c_condition = need_col(&headers, "condition", samples_path)?;
    let c_pair = need_col(&headers, paired_by, samples_path).with_context(|| {
        format!(
            "samples.tsv must include the pairing column {:?}; rerun atman recover-plex --infer-pairs to add it",
            paired_by
        )
    })?;

    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        let pair_value = row[c_pair].to_string();
        if pair_value.is_empty() {
            continue;
        }
        out.push((
            row[c_sample].to_string(),
            row[c_condition].to_string(),
            pair_value,
        ));
    }
    Ok(out)
}

fn collect_pairs(
    pairing: &[(String, String, String)],
    cond_a: &str,
    cond_b: &str,
) -> Result<Vec<PairUnit>> {
    let mut by_pair: BTreeMap<String, (Option<String>, Option<String>)> = BTreeMap::new();
    for (sid, cond, pid) in pairing {
        let entry = by_pair.entry(pid.clone()).or_insert((None, None));
        if cond == cond_a {
            if entry.0.is_some() {
                bail!(
                    "pair {:?} has more than one {} sample ({} and {})",
                    pid,
                    cond_a,
                    entry.0.as_ref().unwrap(),
                    sid
                );
            }
            entry.0 = Some(sid.clone());
        } else if cond == cond_b {
            if entry.1.is_some() {
                bail!(
                    "pair {:?} has more than one {} sample ({} and {})",
                    pid,
                    cond_b,
                    entry.1.as_ref().unwrap(),
                    sid
                );
            }
            entry.1 = Some(sid.clone());
        }
    }
    let mut units = Vec::new();
    for (pid, (a, b)) in by_pair {
        if let (Some(a), Some(b)) = (a, b) {
            units.push(PairUnit {
                pair_id: pid,
                sample_a: a,
                sample_b: b,
            });
        }
    }
    Ok(units)
}

fn write_rows(path: &std::path::Path, rows: &[OutRow]) -> Result<()> {
    let mut buf = String::from(
        "comparison\tpanel\tassay_id\tgene_symbol\tn_pairs_full\tmean_diff_full\tt_full\tp_full\tq_full\tloso_sign_stability\tloso_p_lt_05_stability\tloso_min_t\tloso_max_t\tloso_max_abs_t_delta\tmost_influential_pair_id\tskip_reason\n",
    );
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.comparison,
            r.panel,
            r.assay_id,
            r.gene_symbol,
            r.n_pairs_full,
            opt_f(r.mean_diff_full),
            opt_f(r.t_full),
            opt_f(r.p_full),
            opt_f(r.q_full),
            opt_f(r.loso_sign_stability),
            opt_f(r.loso_p_lt_05_stability),
            opt_f(r.loso_min_t),
            opt_f(r.loso_max_t),
            opt_f(r.loso_max_abs_t_delta),
            r.most_influential_pair_id,
            r.skip_reason,
        ));
    }
    atomic_write(path, buf.as_bytes()).with_context(|| format!("writing {:?}", path))
}

fn summary_path_for(out: &std::path::Path) -> PathBuf {
    let stem = out
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "robust_paired".to_string());
    let ext = out
        .extension()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "tsv".to_string());
    out.with_file_name(format!("{}_summary.{}", stem, ext))
}

fn write_summary(path: &std::path::Path, summaries: &[SummaryRow]) -> Result<()> {
    let mut buf = String::from(
        "comparison\tn_proteins\tn_proteins_full_q05\tn_proteins_robust\tfraction_robust_of_significant\tmedian_sign_stability_of_significant\n",
    );
    for s in summaries {
        // Route both floats through the same empty-cell formatter used for
        // per-row values: an undefined median (no significant proteins) becomes
        // an empty cell rather than the literal `NaN`.
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            s.comparison,
            s.n_proteins,
            s.n_proteins_full_q05,
            s.n_proteins_robust,
            opt_f(Some(s.fraction_robust_of_significant)),
            opt_f(Some(s.median_sign_stability_of_significant)),
        ));
    }
    atomic_write(path, buf.as_bytes()).with_context(|| format!("writing {:?}", path))
}

fn opt_f(x: Option<f64>) -> String {
    match x {
        None => String::new(),
        Some(v) if v.is_nan() => String::new(),
        Some(v) if v.is_infinite() => {
            if v > 0.0 {
                "inf".to_string()
            } else {
                "-inf".to_string()
            }
        }
        Some(v) => format!("{:.6}", v),
    }
}

fn median_f(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        0.5 * (s[n / 2 - 1] + s[n / 2])
    }
}
