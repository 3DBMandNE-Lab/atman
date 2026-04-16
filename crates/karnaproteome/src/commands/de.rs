//! Differential abundance command. Paired Student's t-test at subject level
//! with BH-FDR per comparison family. See spec §13.2 and the v0.2 analysis
//! design.

use anyhow::{Context, Result};
use clap::Args as ClapArgs;
use proteome_core::de::{bh_fdr, paired_t, PairedTResult, SkipReason};
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;

use super::parse_comparisons;
use crate::io::{
    read_measurements_long, read_proteins, read_samples, write_de_report, write_de_results,
    DeReportRow, DeResultRow,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing qc_measurements.tsv, samples.tsv, proteins.tsv.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for de_results.tsv + de_report.tsv.
    #[arg(long)]
    output_dir: PathBuf,

    /// Test kind: `paired-t` (default) or `moderated`.
    #[arg(long, default_value = "paired-t")]
    test: String,

    /// Biological-replicate key column on samples.tsv. v0.1-analysis expects
    /// `participant` (the Dube subject_id).
    #[arg(long, default_value = "participant")]
    paired_by: String,

    /// Comma-separated comparisons in `A-B` form. Each is a separate
    /// hypothesis family for FDR. Example: "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2".
    #[arg(long)]
    groups: String,

    /// Minimum number of paired subjects for a test to run. Below this, the
    /// protein is emitted as Skipped(InsufficientPairs) with NaN p/q.
    #[arg(long, default_value_t = 5)]
    min_pairs: usize,

    /// Prior degrees of freedom for `--test moderated` variance shrinkage.
    #[arg(long, default_value_t = 4.0)]
    moderation_prior_df: f64,
}

pub fn run(args: Args) -> Result<()> {
    if args.test != "paired-t" && args.test != "moderated" {
        anyhow::bail!("test {:?} not supported; use paired-t or moderated", args.test);
    }
    if args.paired_by != "participant" {
        anyhow::bail!(
            "paired-by {:?} not supported in v0.1-analysis (expected `participant`)",
            args.paired_by
        );
    }
    if args.min_pairs < 2 {
        anyhow::bail!("min-pairs must be >= 2 for paired t-test");
    }
    if args.test == "moderated"
        && (!args.moderation_prior_df.is_finite() || args.moderation_prior_df <= 0.0)
    {
        anyhow::bail!("moderation-prior-df must be > 0 for moderated test");
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    // Parse comparisons.
    let comparisons = parse_comparisons(&args.groups)?;

    // Read inputs.
    let measurements = read_measurements_long(&args.input_dir.join("qc_measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let proteins = read_proteins(&args.input_dir.join("proteins.tsv"))?;

    // Lookup: sample_id → (subject, condition, is_control)
    let sample_by_id: HashMap<&str, &proteome_core::Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    // Lookup: (panel, gene_symbol) → (first_assay_id, uniprot) for reporting.
    // Gene symbols are unique within a panel for Olink Explore; we pick the
    // first OlinkID encountered as the canonical AssayId for the row.
    let mut gene_meta: BTreeMap<(String, String), (String, Vec<String>)> = BTreeMap::new();
    for p in &proteins {
        if let Some(gene) = &p.gene_symbol {
            if let Some(panel) = &p.panel {
                gene_meta
                    .entry((panel.clone(), gene.clone()))
                    .or_insert_with(|| (p.assay_id.0.clone(), p.uniprot.clone()));
            }
        }
    }

    // Build per-(panel, gene_symbol, condition) → Vec<(subject, value)>.
    // This is the paired-t input shape. We drop QC-masked measurements,
    // control samples, and rows missing subject/condition metadata.
    type CellsByCondition = BTreeMap<String, Vec<(String, f64)>>;
    let mut cells: BTreeMap<(String, String), CellsByCondition> = BTreeMap::new();

    for m in &measurements {
        let abundance = match m.effective_abundance() {
            Some(a) => a,
            None => continue, // QC-masked
        };
        let s = match sample_by_id.get(m.sample_id.as_str()) {
            Some(s) => *s,
            None => continue,
        };
        if s.is_control {
            continue;
        }
        let subject = match &s.subject_id {
            Some(id) => id.clone(),
            None => continue,
        };
        let condition = match &s.condition {
            Some(c) => c.clone(),
            None => continue,
        };
        let panel = match &m.panel {
            Some(p) => p.clone(),
            None => continue,
        };
        let gene = match &m.gene_symbol {
            Some(g) => g.clone(),
            None => continue,
        };
        cells
            .entry((panel, gene))
            .or_default()
            .entry(condition)
            .or_default()
            .push((subject, abundance));
    }

    // For each comparison, run paired_t on every (panel, gene) and collect.
    // Results are accumulated in a flat Vec for the final TSV + a per-family
    // Vec for BH-FDR correction. A family = one comparison.
    let mut all_rows: Vec<DeResultRow> = Vec::new();
    // For the report sidecar.
    let mut report_rows: Vec<DeReportRow> = Vec::new();

    for (comp_a, comp_b) in &comparisons {
        let comparison_label = format!("{}-{}", comp_a, comp_b);
        // Per-family p-value vector aligned with `family_rows` order.
        let mut family_p: Vec<Option<f64>> = Vec::new();
        let mut family_rows: Vec<DeResultRow> = Vec::new();

        for ((panel, gene), by_condition) in &cells {
            let (assay_id, uniprot) = match gene_meta.get(&(panel.clone(), gene.clone())) {
                Some(meta) => meta.clone(),
                None => (String::new(), vec![]),
            };

            // Build paired vector: for each subject present in BOTH conditions.
            let empty: Vec<(String, f64)> = Vec::new();
            let va = by_condition.get(comp_a).unwrap_or(&empty);
            let vb = by_condition.get(comp_b).unwrap_or(&empty);

            // Subject → value, for fast paired join.
            let a_by_subj: HashMap<&str, f64> = va.iter().map(|(s, v)| (s.as_str(), *v)).collect();
            let b_by_subj: HashMap<&str, f64> = vb.iter().map(|(s, v)| (s.as_str(), *v)).collect();

            let pairs: Vec<(f64, f64)> = a_by_subj
                .iter()
                .filter_map(|(subj, a)| b_by_subj.get(subj).map(|b| (*a, *b)))
                .collect();

            let result = paired_t(&pairs, args.min_pairs);
            let row = match &result {
                PairedTResult::Computed {
                    n_pairs,
                    mean_a,
                    mean_b,
                    mean_diff,
                    t,
                    df,
                    p_value,
                } => DeResultRow {
                    panel: panel.clone(),
                    assay_id: assay_id.clone(),
                    gene_symbol: gene.clone(),
                    uniprot: uniprot.join(","),
                    comparison: comparison_label.clone(),
                    n_pairs: *n_pairs,
                    mean_a: Some(*mean_a),
                    mean_b: Some(*mean_b),
                    mean_diff: Some(*mean_diff),
                    t: Some(*t),
                    df: Some(*df),
                    p_value: Some(*p_value),
                    bh_q: None, // filled after BH sweep
                    skip_reason: String::new(),
                },
                PairedTResult::Skipped { reason, n_pairs } => DeResultRow {
                    panel: panel.clone(),
                    assay_id: assay_id.clone(),
                    gene_symbol: gene.clone(),
                    uniprot: uniprot.join(","),
                    comparison: comparison_label.clone(),
                    n_pairs: *n_pairs,
                    mean_a: None,
                    mean_b: None,
                    mean_diff: None,
                    t: None,
                    df: None,
                    p_value: None,
                    bh_q: None,
                    skip_reason: match reason {
                        SkipReason::InsufficientPairs => "insufficient_pairs".to_string(),
                        SkipReason::ZeroVariance => "zero_variance".to_string(),
                        SkipReason::NonFiniteInput => "non_finite_input".to_string(),
                    },
                },
            };
            family_p.push(row.p_value);
            family_rows.push(row);
        }

        if args.test == "moderated" {
            apply_moderated_shrinkage(&mut family_rows, args.moderation_prior_df)?;
            family_p = family_rows.iter().map(|r| r.p_value).collect();
        }

        // BH-FDR within this comparison family.
        let qs = bh_fdr(&family_p);
        for (row, q) in family_rows.iter_mut().zip(qs.into_iter()) {
            row.bh_q = q;
        }

        // Per-panel report summary.
        let mut per_panel: BTreeMap<String, ReportAccumulator> = BTreeMap::new();
        for row in &family_rows {
            let acc = per_panel.entry(row.panel.clone()).or_default();
            acc.n_tests += 1;
            if let Some(q) = row.bh_q {
                if q < 0.05 {
                    acc.n_q_lt_05 += 1;
                }
                if q < 0.10 {
                    acc.n_q_lt_10 += 1;
                }
                if acc.min_q.map(|m| q < m).unwrap_or(true) {
                    acc.min_q = Some(q);
                }
            } else {
                acc.n_skipped += 1;
            }
            if let Some(d) = row.mean_diff {
                let abs_d = d.abs();
                if acc.max_abs_effect.map(|m| abs_d > m).unwrap_or(true) {
                    acc.max_abs_effect = Some(abs_d);
                }
            }
        }
        for (panel, acc) in per_panel {
            report_rows.push(DeReportRow {
                comparison: comparison_label.clone(),
                panel,
                n_tests: acc.n_tests,
                n_skipped: acc.n_skipped,
                n_q_lt_05: acc.n_q_lt_05,
                n_q_lt_10: acc.n_q_lt_10,
                min_q: acc.min_q,
                max_abs_effect: acc.max_abs_effect,
            });
        }

        all_rows.extend(family_rows);
    }

    // Stable sort for inspection: (comparison, bh_q asc with None last, |mean_diff| desc).
    all_rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| match (a.bh_q, b.bh_q) {
                (Some(qa), Some(qb)) => qa.partial_cmp(&qb).unwrap_or(std::cmp::Ordering::Equal),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            })
            .then_with(|| match (a.mean_diff, b.mean_diff) {
                (Some(da), Some(db)) => db
                    .abs()
                    .partial_cmp(&da.abs())
                    .unwrap_or(std::cmp::Ordering::Equal),
                _ => std::cmp::Ordering::Equal,
            })
            .then_with(|| a.gene_symbol.cmp(&b.gene_symbol))
    });

    write_de_results(&args.output_dir.join("de_results.tsv"), &all_rows)?;
    write_de_report(&args.output_dir.join("de_report.tsv"), &report_rows)?;

    let computed = all_rows.iter().filter(|r| r.p_value.is_some()).count();
    let skipped = all_rows.len() - computed;
    eprintln!(
        "de: test={} comparisons={} rows={} computed={} skipped={} min_pairs={} prior_df={}",
        args.test,
        comparisons.len(),
        all_rows.len(),
        computed,
        skipped,
        args.min_pairs,
        args.moderation_prior_df,
    );
    Ok(())
}

#[derive(Default)]
struct ReportAccumulator {
    n_tests: usize,
    n_skipped: usize,
    n_q_lt_05: usize,
    n_q_lt_10: usize,
    min_q: Option<f64>,
    max_abs_effect: Option<f64>,
}

fn apply_moderated_shrinkage(rows: &mut [DeResultRow], prior_df: f64) -> Result<()> {
    let mut vars: Vec<f64> = rows
        .iter()
        .filter_map(paired_variance_from_row)
        .filter(|v| v.is_finite() && *v > 0.0)
        .collect();
    if vars.is_empty() {
        return Ok(());
    }
    vars.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let prior_var = if vars.len() % 2 == 1 {
        vars[vars.len() / 2]
    } else {
        (vars[vars.len() / 2 - 1] + vars[vars.len() / 2]) / 2.0
    };

    for row in rows.iter_mut() {
        let mean_diff = match row.mean_diff {
            Some(v) => v,
            None => continue,
        };
        let df_i = match row.df {
            Some(v) if v > 0.0 => v,
            _ => continue,
        };
        let n = row.n_pairs as f64;
        if n < 2.0 {
            continue;
        }
        let var_i = match paired_variance_from_row(row) {
            Some(v) if v.is_finite() && v >= 0.0 => v,
            _ => continue,
        };
        let post_var = (prior_df * prior_var + df_i * var_i) / (prior_df + df_i);
        if !(post_var.is_finite() && post_var > 0.0) {
            continue;
        }
        let t_mod = mean_diff / (post_var / n).sqrt();
        let df_mod = prior_df + df_i;
        let dist = StudentsT::new(0.0, 1.0, df_mod)
            .with_context(|| format!("building t distribution with df={}", df_mod))?;
        let p_mod = 2.0 * (1.0 - dist.cdf(t_mod.abs()));
        if p_mod.is_finite() {
            row.t = Some(t_mod);
            row.df = Some(df_mod);
            row.p_value = Some(p_mod);
        }
    }
    Ok(())
}

fn paired_variance_from_row(row: &DeResultRow) -> Option<f64> {
    let t = row.t?;
    let mean = row.mean_diff?;
    let n = row.n_pairs as f64;
    if n < 2.0 {
        return None;
    }
    if t == 0.0 {
        return Some(0.0);
    }
    let se = mean / t;
    Some((se * n.sqrt()).powi(2))
}
