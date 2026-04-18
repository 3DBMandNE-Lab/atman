//! Differential abundance command.

use anyhow::{anyhow, Context, Result};
use atman_core::de::{bh_fdr, ols, paired_t, welch_t, OlsOutcome, PairedTResult, SkipReason};
use atman_core::Sample;
use clap::Args as ClapArgs;
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

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

    /// Test kind: `paired-t` (default), `moderated` (paired variance
    /// shrinkage), `welch-t` (unpaired two-sample with Welch-Satterthwaite
    /// df), or `ols` (covariate-adjusted unpaired linear model).
    #[arg(long, default_value = "paired-t")]
    test: String,

    /// Biological-replicate key column on samples.tsv. Required for
    /// `paired-t` and `moderated`; ignored for `welch-t` and `ols`.
    #[arg(long, default_value = "participant")]
    paired_by: String,

    /// Comma-separated covariate column names from samples.tsv for
    /// `--test ols`. Numeric columns (parseable as f64 in every non-missing
    /// row) are used as continuous predictors; everything else is treated
    /// as categorical and one-hot encoded, dropping the alphabetically
    /// first level as reference. Samples with any missing covariate value
    /// are dropped before fitting. Ignored for other tests.
    #[arg(long)]
    covariates: Option<String>,

    /// Comma-separated comparisons in `A-B` form. Each is a separate
    /// hypothesis family for FDR. Example: "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2".
    #[arg(long)]
    groups: String,

    /// Minimum number of samples required per group for a test to run.
    /// For paired tests this is the number of matched pairs; for `welch-t`
    /// it is the minimum of `|group_a|` and `|group_b|`. Below this the
    /// protein is emitted as Skipped(InsufficientPairs) with NaN p/q.
    #[arg(long, default_value_t = 5)]
    min_pairs: usize,

    /// Prior degrees of freedom for `--test moderated` variance shrinkage.
    #[arg(long, default_value_t = 4.0)]
    moderation_prior_df: f64,
}

pub fn run(args: Args) -> Result<()> {
    if args.test != "paired-t"
        && args.test != "moderated"
        && args.test != "welch-t"
        && args.test != "ols"
    {
        anyhow::bail!(
            "test {:?} not supported; use paired-t, moderated, welch-t, or ols",
            args.test
        );
    }
    let is_unpaired = args.test == "welch-t" || args.test == "ols";
    if !is_unpaired && args.paired_by != "participant" {
        anyhow::bail!(
            "paired-by {:?} not supported for paired tests (expected `participant`)",
            args.paired_by
        );
    }
    if args.min_pairs < 2 {
        anyhow::bail!("min-pairs must be >= 2");
    }
    if args.test == "moderated"
        && (!args.moderation_prior_df.is_finite() || args.moderation_prior_df <= 0.0)
    {
        anyhow::bail!("moderation-prior-df must be > 0 for moderated test");
    }
    if args.test != "ols" && args.covariates.is_some() {
        anyhow::bail!("--covariates is only valid with --test ols");
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
    let sample_by_id: HashMap<&str, &atman_core::Sample> =
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

    // Build per-(panel, gene_symbol, condition) → Vec<(subject_or_sample, value)>.
    // For paired-t / moderated the first element is subject_id (for the paired
    // join); for welch-t it is unused. For `ols` we build a separate
    // sample-level structure below and use `cells_by_sample` instead.
    type CellsByCondition = BTreeMap<String, Vec<(String, f64)>>;
    let mut cells: BTreeMap<(String, String), CellsByCondition> = BTreeMap::new();
    // (panel, gene, sample_id) → abundance — only built for `--test ols` so
    // per-protein fits can look up abundance per sample and drop QC-masked
    // cells without needing to pair by subject.
    let mut cells_by_sample: HashMap<(String, String, String), f64> = HashMap::new();

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
        if args.test == "ols" {
            cells_by_sample.insert(
                (panel.clone(), gene.clone(), m.sample_id.clone()),
                abundance,
            );
        }
        cells
            .entry((panel, gene))
            .or_default()
            .entry(condition)
            .or_default()
            .push((subject, abundance));
    }

    // Covariate setup for `--test ols`. Parse the comma list, read
    // samples.tsv as raw records to pull arbitrary columns by header name,
    // then classify each as numeric or categorical. Categorical columns get
    // one-hot-encoded later, per comparison, after the sample subset is
    // known.
    let (cov_names, cov_raw) = if args.test == "ols" {
        let names: Vec<String> = args
            .covariates
            .as_deref()
            .unwrap_or("")
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();
        let raw = if names.is_empty() {
            HashMap::new()
        } else {
            read_covariate_columns(&args.input_dir.join("samples.tsv"), &names)?
        };
        (names, raw)
    } else {
        (Vec::new(), HashMap::new())
    };

    // For each comparison, run paired_t on every (panel, gene) and collect.
    // Results are accumulated in a flat Vec for the final TSV + a per-family
    // Vec for BH-FDR correction. A family = one comparison.
    let mut all_rows: Vec<DeResultRow> = Vec::new();
    // For the report sidecar.
    let mut report_rows: Vec<DeReportRow> = Vec::new();

    // Per-covariate estimates across all proteins × comparisons, emitted to
    // de_covariates.tsv when --test ols. One row per (comparison, panel,
    // gene, covariate) with beta, se, t, p for the non-intercept and
    // non-group columns.
    let mut covariate_rows: Vec<CovariateRow> = Vec::new();

    for (comp_a, comp_b) in &comparisons {
        let comparison_label = format!("{}-{}", comp_a, comp_b);
        // Per-family p-value vector aligned with `family_rows` order.
        let mut family_p: Vec<Option<f64>> = Vec::new();
        let mut family_rows: Vec<DeResultRow> = Vec::new();

        // OLS mode: precompute the encoded design matrix (excluding y) for
        // this comparison so every per-protein fit reuses it, only swapping
        // the y vector (with complete-case filtering on non-finite abundance).
        let ols_design: Option<OlsDesign> = if args.test == "ols" {
            Some(build_ols_design(
                &samples, comp_a, comp_b, &cov_names, &cov_raw,
            )?)
        } else {
            None
        };

        for ((panel, gene), by_condition) in &cells {
            let (assay_id, uniprot) = match gene_meta.get(&(panel.clone(), gene.clone())) {
                Some(meta) => meta.clone(),
                None => (String::new(), vec![]),
            };

            // Build per-condition cell vectors. In paired mode we intersect
            // by subject_id; in unpaired (welch-t) mode we treat each
            // condition's samples as an independent group and drop the
            // subject linkage.
            let empty: Vec<(String, f64)> = Vec::new();
            let va = by_condition.get(comp_a).unwrap_or(&empty);
            let vb = by_condition.get(comp_b).unwrap_or(&empty);

            let result = if args.test == "ols" {
                let design = ols_design
                    .as_ref()
                    .expect("ols_design populated when test=ols");
                // Look up per-sample abundance for this protein; drop
                // samples where the cell is QC-masked (missing from
                // cells_by_sample).
                let mut y: Vec<f64> = Vec::with_capacity(design.rows.len());
                let mut rows: Vec<Vec<f64>> = Vec::with_capacity(design.rows.len());
                for (sid, row) in design.rows.iter() {
                    if let Some(abundance) =
                        cells_by_sample.get(&(panel.clone(), gene.clone(), sid.clone()))
                    {
                        y.push(*abundance);
                        rows.push(row.clone());
                    }
                }
                let fit = ols(&rows, &y, args.min_pairs);
                // Compute raw group means for schema compatibility with
                // paired-t / welch-t.
                let (mean_a_raw, mean_b_raw) = raw_group_means(&y, &rows, design.group_col);
                ols_to_paired_t_result(
                    fit,
                    design,
                    mean_a_raw,
                    mean_b_raw,
                    panel,
                    gene,
                    &comparison_label,
                    &mut covariate_rows,
                )
            } else if is_unpaired {
                let a_vals: Vec<f64> = va.iter().map(|(_, v)| *v).collect();
                let b_vals: Vec<f64> = vb.iter().map(|(_, v)| *v).collect();
                welch_t(&a_vals, &b_vals, args.min_pairs)
            } else {
                // Subject → value, for fast paired join.
                let a_by_subj: HashMap<&str, f64> =
                    va.iter().map(|(s, v)| (s.as_str(), *v)).collect();
                let b_by_subj: HashMap<&str, f64> =
                    vb.iter().map(|(s, v)| (s.as_str(), *v)).collect();
                let pairs: Vec<(f64, f64)> = a_by_subj
                    .iter()
                    .filter_map(|(subj, a)| b_by_subj.get(subj).map(|b| (*a, *b)))
                    .collect();
                paired_t(&pairs, args.min_pairs)
            };
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
    if args.test == "ols" && !covariate_rows.is_empty() {
        write_covariate_rows(&args.output_dir.join("de_covariates.tsv"), &covariate_rows)?;
    }

    let computed = all_rows.iter().filter(|r| r.p_value.is_some()).count();
    let skipped = all_rows.len() - computed;
    eprintln!(
        "de: test={} comparisons={} rows={} computed={} skipped={} min_pairs={} prior_df={} covariates={}",
        args.test,
        comparisons.len(),
        all_rows.len(),
        computed,
        skipped,
        args.min_pairs,
        args.moderation_prior_df,
        cov_names.join(","),
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

// -----------------------------------------------------------------
// OLS (covariate-adjusted) helpers
// -----------------------------------------------------------------

/// Covariate kind determined at column-classification time. `Numeric`
/// columns are passed through as f64. `Categorical` columns are one-hot
/// encoded using the supplied levels; the first element of `levels` is the
/// reference and is dropped.
#[derive(Debug, Clone)]
enum CovKind {
    Numeric,
    Categorical { levels: Vec<String> },
}

impl CovKind {
    /// Number of design columns this covariate contributes.
    fn width(&self) -> usize {
        match self {
            CovKind::Numeric => 1,
            CovKind::Categorical { levels } => levels.len().saturating_sub(1),
        }
    }

    /// Human label for each design column this covariate contributes.
    /// Used for the header of de_covariates.tsv.
    fn col_labels(&self, base: &str) -> Vec<String> {
        match self {
            CovKind::Numeric => vec![base.to_string()],
            CovKind::Categorical { levels } => levels
                .iter()
                .skip(1)
                .map(|lvl| format!("{}={}", base, lvl))
                .collect(),
        }
    }
}

/// Encoded design matrix for one comparison. `rows` lists every sample
/// that (a) is in `comp_a` or `comp_b`, (b) is non-control, and (c) has a
/// non-empty value for every covariate; each entry is
/// `(sample_id, [intercept, group, cov_col_0, cov_col_1, …])`. `group_col`
/// is the index of the group column (always 1 under the current design)
/// so the caller knows which β to extract as the effect of interest.
/// `design_labels[i]` names the i-th design column for de_covariates.tsv.
struct OlsDesign {
    rows: Vec<(String, Vec<f64>)>,
    group_col: usize,
    design_labels: Vec<String>,
}

/// Row shape for de_covariates.tsv — per (comparison × protein × covariate
/// column) estimate, skipping the intercept and the group indicator.
struct CovariateRow {
    comparison: String,
    panel: String,
    gene_symbol: String,
    covariate: String,
    beta: f64,
    se: f64,
    t: f64,
    p_value: f64,
    df: f64,
    n: usize,
}

/// Raw-read `samples.tsv` and return `sample_id → [cov_value]` for the
/// requested covariate names. Values that are missing (empty cell) are
/// returned as `None`; the caller drops them before fitting.
fn read_covariate_columns(
    path: &Path,
    names: &[String],
) -> Result<HashMap<String, Vec<Option<String>>>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let sample_col = headers
        .iter()
        .position(|h| h == "sample_id")
        .ok_or_else(|| anyhow!("samples.tsv missing 'sample_id' column"))?;
    let cov_cols: Vec<usize> = names
        .iter()
        .map(|n| {
            headers
                .iter()
                .position(|h| h == n.as_str())
                .ok_or_else(|| anyhow!("covariate {:?} not found in samples.tsv", n))
        })
        .collect::<Result<_>>()?;
    let mut out: HashMap<String, Vec<Option<String>>> = HashMap::new();
    for row in reader.records() {
        let row = row?;
        let sid = row
            .get(sample_col)
            .ok_or_else(|| anyhow!("samples.tsv row missing sample_id"))?
            .to_string();
        let vals: Vec<Option<String>> = cov_cols
            .iter()
            .map(|&i| {
                let v = row.get(i).unwrap_or("");
                if v.is_empty() {
                    None
                } else {
                    Some(v.to_string())
                }
            })
            .collect();
        out.insert(sid, vals);
    }
    Ok(out)
}

/// Classify each covariate as Numeric or Categorical based on whether
/// every non-missing value parses as f64.
fn classify_covariates(
    names: &[String],
    values_per_sample: &[Vec<Option<String>>],
) -> Vec<CovKind> {
    let mut out = Vec::with_capacity(names.len());
    for (i, _) in names.iter().enumerate() {
        let mut all_numeric = true;
        let mut levels: BTreeSet<String> = BTreeSet::new();
        for row in values_per_sample {
            if let Some(v) = &row[i] {
                if v.parse::<f64>().is_err() {
                    all_numeric = false;
                }
                levels.insert(v.clone());
            }
        }
        if all_numeric {
            out.push(CovKind::Numeric);
        } else {
            out.push(CovKind::Categorical {
                levels: levels.into_iter().collect(),
            });
        }
    }
    out
}

/// Build the design matrix for one comparison. Samples not in `comp_a` or
/// `comp_b`, control samples, and samples with any missing covariate are
/// dropped before encoding.
fn build_ols_design(
    samples: &[Sample],
    comp_a: &str,
    comp_b: &str,
    cov_names: &[String],
    cov_raw: &HashMap<String, Vec<Option<String>>>,
) -> Result<OlsDesign> {
    // Pick samples in either group.
    let in_group: Vec<&Sample> = samples
        .iter()
        .filter(|s| {
            !s.is_control
                && match &s.condition {
                    Some(c) => c == comp_a || c == comp_b,
                    None => false,
                }
        })
        .collect();

    // Pull covariate values aligned with `in_group`; drop samples that
    // are missing any covariate.
    let mut kept: Vec<(&Sample, Vec<Option<String>>)> = Vec::new();
    for s in &in_group {
        let vals = if cov_names.is_empty() {
            Vec::new()
        } else {
            match cov_raw.get(&s.sample_id) {
                Some(v) => v.clone(),
                None => continue,
            }
        };
        if vals.iter().any(|v| v.is_none()) {
            continue;
        }
        kept.push((s, vals));
    }

    // Classify covariate kinds from the retained samples.
    let kinds = classify_covariates(
        cov_names,
        &kept.iter().map(|(_, v)| v.clone()).collect::<Vec<_>>(),
    );

    // Assemble encoded rows: [1 (intercept), group (0/1 for comp_b/comp_a), encoded covariate columns...].
    // We encode group = 1 when the sample's condition is comp_a and 0 when
    // comp_b so `beta_group` is the effect of being in `comp_a` relative to
    // `comp_b` — i.e., positive β means comp_a is higher, matching the
    // paired-t / welch-t convention where `mean_diff = mean_a − mean_b`.
    let mut rows: Vec<(String, Vec<f64>)> = Vec::with_capacity(kept.len());
    for (s, vals) in &kept {
        let mut row: Vec<f64> =
            Vec::with_capacity(2 + kinds.iter().map(|k| k.width()).sum::<usize>());
        row.push(1.0);
        row.push(if s.condition.as_deref() == Some(comp_a) {
            1.0
        } else {
            0.0
        });
        for (i, kind) in kinds.iter().enumerate() {
            let v = vals[i]
                .as_ref()
                .expect("missing covariate filtered out above");
            match kind {
                CovKind::Numeric => {
                    row.push(v.parse::<f64>().unwrap_or(f64::NAN));
                }
                CovKind::Categorical { levels } => {
                    // one-hot drop-first
                    for lvl in levels.iter().skip(1) {
                        row.push(if v == lvl { 1.0 } else { 0.0 });
                    }
                }
            }
        }
        rows.push((s.sample_id.clone(), row));
    }

    // Column labels for reporting.
    let mut design_labels = vec![
        "(Intercept)".to_string(),
        format!("{}_vs_{}", comp_a, comp_b),
    ];
    for (name, kind) in cov_names.iter().zip(kinds.iter()) {
        design_labels.extend(kind.col_labels(name));
    }

    Ok(OlsDesign {
        rows,
        group_col: 1,
        design_labels,
    })
}

/// Raw per-group means of the response within the fitted rows. Used to
/// populate `mean_a` / `mean_b` on the `DeResultRow` for schema
/// compatibility with paired-t / welch-t output.
fn raw_group_means(y: &[f64], rows: &[Vec<f64>], group_col: usize) -> (Option<f64>, Option<f64>) {
    let mut sum_a = 0.0;
    let mut n_a = 0;
    let mut sum_b = 0.0;
    let mut n_b = 0;
    for (yi, row) in y.iter().zip(rows.iter()) {
        if row[group_col] > 0.5 {
            sum_a += yi;
            n_a += 1;
        } else {
            sum_b += yi;
            n_b += 1;
        }
    }
    let ma = if n_a > 0 {
        Some(sum_a / n_a as f64)
    } else {
        None
    };
    let mb = if n_b > 0 {
        Some(sum_b / n_b as f64)
    } else {
        None
    };
    (ma, mb)
}

/// Convert an `OlsOutcome` into the `PairedTResult` enum used by the
/// downstream writer. The group coefficient becomes `mean_diff`; residual
/// df becomes `df`; raw group means (supplied) become `mean_a` / `mean_b`.
/// Non-intercept, non-group coefficients are recorded into
/// `covariate_rows` for emission to de_covariates.tsv.
#[allow(clippy::too_many_arguments)]
fn ols_to_paired_t_result(
    outcome: OlsOutcome,
    design: &OlsDesign,
    mean_a: Option<f64>,
    mean_b: Option<f64>,
    panel: &str,
    gene: &str,
    comparison: &str,
    covariate_rows: &mut Vec<CovariateRow>,
) -> PairedTResult {
    match outcome {
        OlsOutcome::Computed(fit) => {
            let g = design.group_col;
            // Record per-covariate estimates (everything except intercept and group).
            for (col, label) in design.design_labels.iter().enumerate() {
                if col == 0 || col == g {
                    continue;
                }
                covariate_rows.push(CovariateRow {
                    comparison: comparison.to_string(),
                    panel: panel.to_string(),
                    gene_symbol: gene.to_string(),
                    covariate: label.clone(),
                    beta: fit.beta[col],
                    se: fit.se[col],
                    t: fit.t[col],
                    p_value: fit.p_value[col],
                    df: fit.df,
                    n: fit.n,
                });
            }
            PairedTResult::Computed {
                n_pairs: fit.n,
                mean_a: mean_a.unwrap_or(f64::NAN),
                mean_b: mean_b.unwrap_or(f64::NAN),
                mean_diff: fit.beta[g],
                t: fit.t[g],
                df: fit.df,
                p_value: fit.p_value[g],
            }
        }
        OlsOutcome::Skipped { reason, n } => PairedTResult::Skipped { reason, n_pairs: n },
    }
}

/// Writer for `de_covariates.tsv`. Tab-separated, one row per
/// (comparison × panel × gene × covariate column).
fn write_covariate_rows(path: &Path, rows: &[CovariateRow]) -> Result<()> {
    let mut buf =
        String::from("comparison\tpanel\tgene_symbol\tcovariate\tbeta\tse\tt\tp_value\tdf\tn\n");
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.comparison,
            r.panel,
            r.gene_symbol,
            r.covariate,
            r.beta,
            r.se,
            r.t,
            r.p_value,
            r.df,
            r.n,
        ));
    }
    crate::io::atomic_write(path, buf.as_bytes())
}
