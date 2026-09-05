//! Design-matrix construction and OLS plumbing shared by the
//! covariate-adjusted (`--test ols`), mixed (`--test mixed`), and
//! post-hoc OLS pipelines: formula parsing, comparison resolution,
//! covariate classification/encoding, full-rank checks, and the
//! `OlsOutcome → PairedTResult` bridge.

use anyhow::{anyhow, Context, Result};
use atman_core::de::{OlsOutcome, PairedTResult};
use atman_core::Sample;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

use super::parse_comparisons;
use super::{CovariateRow, DesignReportRow};

// -----------------------------------------------------------------
// OLS (covariate-adjusted) helpers
// -----------------------------------------------------------------

#[derive(Debug, Clone)]
pub(super) enum DesignTerm {
    Condition,
    Covariate(String),
}

pub(super) struct OlsSetup {
    pub(super) terms: Vec<DesignTerm>,
    pub(super) cov_names: Vec<String>,
    pub(super) cov_raw: HashMap<String, Vec<Option<String>>>,
    pub(super) contrast: Option<String>,
    pub(super) label: String,
}

pub(super) fn build_ols_setup(
    samples_path: &Path,
    design: Option<&str>,
    covariates: Option<&str>,
    contrast: Option<&str>,
) -> Result<OlsSetup> {
    if let Some(formula) = design {
        let terms = parse_design_terms(formula)?;
        if !terms.iter().any(|t| matches!(t, DesignTerm::Condition)) {
            anyhow::bail!("--design must include `condition` for DE contrasts");
        }
        let cov_names = covariate_names_from_terms(&terms);
        let cov_raw = if cov_names.is_empty() {
            HashMap::new()
        } else {
            read_covariate_columns(samples_path, &cov_names)?
        };
        Ok(OlsSetup {
            terms,
            cov_names,
            cov_raw,
            contrast: contrast.map(str::to_string),
            label: formula.to_string(),
        })
    } else {
        let cov_names: Vec<String> = covariates
            .unwrap_or("")
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();
        let cov_raw = if cov_names.is_empty() {
            HashMap::new()
        } else {
            read_covariate_columns(samples_path, &cov_names)?
        };
        let mut terms = vec![DesignTerm::Condition];
        terms.extend(cov_names.iter().cloned().map(DesignTerm::Covariate));
        Ok(OlsSetup {
            terms,
            cov_names,
            cov_raw,
            contrast: contrast.map(str::to_string),
            label: if covariates.unwrap_or("").trim().is_empty() {
                "~ condition".to_string()
            } else {
                format!("~ condition + {}", covariates.unwrap_or("").trim())
            },
        })
    }
}

pub(super) fn parse_design_terms(formula: &str) -> Result<Vec<DesignTerm>> {
    let formula = formula.trim();
    let rhs = formula
        .strip_prefix('~')
        .ok_or_else(|| anyhow!("--design must start with `~`"))?
        .trim();
    if rhs.is_empty() {
        anyhow::bail!("--design must contain at least `condition`");
    }
    let mut terms = Vec::new();
    let mut seen = BTreeSet::new();
    for raw in rhs.split('+') {
        let term = raw.trim();
        if term.is_empty() {
            anyhow::bail!("empty term in --design");
        }
        if term == "1" {
            continue;
        }
        if term == "0" || term == "-1" {
            anyhow::bail!("intercept removal is not supported; Atman includes an intercept");
        }
        if !seen.insert(term.to_string()) {
            anyhow::bail!("duplicate term {:?} in --design", term);
        }
        if term == "condition" {
            terms.push(DesignTerm::Condition);
        } else {
            terms.push(DesignTerm::Covariate(term.to_string()));
        }
    }
    if terms.is_empty() {
        anyhow::bail!("--design must contain at least `condition`");
    }
    Ok(terms)
}

pub(super) fn covariate_names_from_terms(terms: &[DesignTerm]) -> Vec<String> {
    terms
        .iter()
        .filter_map(|t| match t {
            DesignTerm::Condition => None,
            DesignTerm::Covariate(name) => Some(name.clone()),
        })
        .collect()
}

pub(super) fn resolve_comparisons(
    groups: Option<&str>,
    samples: &[Sample],
    test: &str,
    design: Option<&str>,
    contrast: Option<&str>,
) -> Result<Vec<(String, String)>> {
    if let Some(groups) = groups {
        return parse_comparisons(groups);
    }
    if (test != "ols" && test != "mixed") || design.is_none() {
        anyhow::bail!("--groups is required unless using --test ols --design with --contrast");
    }
    let contrast =
        contrast.ok_or_else(|| anyhow!("--contrast is required when --groups is omitted"))?;
    let comp_a = contrast
        .strip_prefix("condition")
        .filter(|s| !s.is_empty())
        .ok_or_else(|| {
            anyhow!("can only infer --groups from condition contrasts like conditionCase")
        })?;
    let conditions: BTreeSet<String> = samples
        .iter()
        .filter(|s| !s.is_control)
        .filter_map(|s| s.condition.clone())
        .collect();
    if !conditions.contains(comp_a) {
        anyhow::bail!(
            "contrast condition {:?} is not present in samples.tsv",
            comp_a
        );
    }
    let others: Vec<String> = conditions.into_iter().filter(|c| c != comp_a).collect();
    if others.len() != 1 {
        anyhow::bail!(
            "cannot infer --groups for contrast {}; provide --groups A-B when more than two conditions are present",
            contrast
        );
    }
    Ok(vec![(comp_a.to_string(), others[0].clone())])
}

pub(super) fn validate_random_intercept_design(
    samples: &[Sample],
    comparisons: &[(String, String)],
    min_samples: usize,
) -> Result<()> {
    for (a, b) in comparisons {
        let mut by_subject: BTreeMap<String, usize> = BTreeMap::new();
        for sample in samples {
            if sample.is_control {
                continue;
            }
            let Some(condition) = &sample.condition else {
                continue;
            };
            if condition != a && condition != b {
                continue;
            }
            let subject = sample.subject_id.as_ref().unwrap_or(&sample.sample_id);
            *by_subject.entry(subject.clone()).or_default() += 1;
        }
        let repeated_subjects = by_subject.values().filter(|n| **n >= 2).count();
        let n_samples: usize = by_subject.values().sum();
        if n_samples < min_samples || repeated_subjects < 2 {
            anyhow::bail!(
                "mixed model for {a}-{b} requires at least {} samples and at least two repeated subjects for --random '1|subject_id'",
                min_samples
            );
        }
    }
    Ok(())
}

/// Covariate kind determined at column-classification time. `Numeric`
/// columns are passed through as f64. `Categorical` columns are one-hot
/// encoded using the supplied levels; the first element of `levels` is the
/// reference and is dropped.
#[derive(Debug, Clone)]
pub(super) enum CovKind {
    Numeric,
    Categorical { levels: Vec<String> },
}

impl CovKind {
    /// Human label for each design column this covariate contributes.
    /// Used for the header of de_covariates.tsv.
    pub(super) fn col_labels(&self, base: &str) -> Vec<String> {
        match self {
            CovKind::Numeric => vec![base.to_string()],
            CovKind::Categorical { levels } => levels
                .iter()
                .skip(1)
                .map(|lvl| format!("{}{}", base, lvl))
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
pub(super) struct OlsDesign {
    pub(super) rows: Vec<(String, Vec<f64>)>,
    pub(super) group_col: usize,
    pub(super) design_labels: Vec<String>,
    pub(super) report_rows: Vec<DesignReportRow>,
    /// True when `group_col` is the 0/1 condition indicator (so raw group
    /// means and Cohen d are meaningful); false for a continuous contrast.
    pub(super) group_is_indicator: bool,
}

/// Raw-read `samples.tsv` and return `sample_id → [cov_value]` for the
/// requested covariate names. Values that are missing (empty cell) are
/// returned as `None`; the caller drops them before fitting.
pub(super) fn read_covariate_columns(
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
pub(super) fn classify_covariates(
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
pub(super) fn build_ols_design(
    samples: &[Sample],
    comp_a: &str,
    comp_b: &str,
    setup: &OlsSetup,
) -> Result<OlsDesign> {
    let comparison = format!("{comp_a}-{comp_b}");
    let condition_label = format!("condition{comp_a}");

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

    // Pull covariate values aligned with `in_group`; drop samples missing
    // any formula covariate before classifying/encoding.
    let mut kept: Vec<(&Sample, Vec<Option<String>>)> = Vec::new();
    let mut report_rows = Vec::new();
    for s in &in_group {
        let vals = if setup.cov_names.is_empty() {
            Vec::new()
        } else {
            match setup.cov_raw.get(&s.sample_id) {
                Some(v) => v.clone(),
                None => {
                    report_rows.push(DesignReportRow {
                        comparison: comparison.clone(),
                        sample_id: s.sample_id.clone(),
                        condition: s.condition.clone().unwrap_or_default(),
                        included: false,
                        drop_reason: "missing_covariate_row".to_string(),
                        columns: String::new(),
                        values: String::new(),
                    });
                    continue;
                }
            }
        };
        if vals.iter().any(|v| v.is_none()) {
            report_rows.push(DesignReportRow {
                comparison: comparison.clone(),
                sample_id: s.sample_id.clone(),
                condition: s.condition.clone().unwrap_or_default(),
                included: false,
                drop_reason: "missing_covariate".to_string(),
                columns: String::new(),
                values: String::new(),
            });
            continue;
        }
        kept.push((s, vals));
    }

    // Classify covariate kinds from the retained samples.
    let cov_kinds = classify_covariates(
        &setup.cov_names,
        &kept.iter().map(|(_, v)| v.clone()).collect::<Vec<_>>(),
    );

    let mut design_labels = vec!["(Intercept)".to_string()];
    for term in &setup.terms {
        match term {
            DesignTerm::Condition => design_labels.push(condition_label.clone()),
            DesignTerm::Covariate(name) => {
                let idx = cov_index(&setup.cov_names, name)?;
                design_labels.extend(cov_kinds[idx].col_labels(name));
            }
        }
    }
    let contrast = setup
        .contrast
        .as_deref()
        .unwrap_or(condition_label.as_str());
    let group_col = design_labels
        .iter()
        .position(|label| label == contrast)
        .ok_or_else(|| {
            anyhow!(
                "contrast {:?} not found in design columns: {}",
                contrast,
                design_labels.join(",")
            )
        })?;

    let mut rows: Vec<(String, Vec<f64>)> = Vec::with_capacity(kept.len());
    for (s, vals) in kept {
        let mut row: Vec<f64> = Vec::with_capacity(design_labels.len());
        row.push(1.0);
        for term in &setup.terms {
            match term {
                DesignTerm::Condition => {
                    row.push(if s.condition.as_deref() == Some(comp_a) {
                        1.0
                    } else {
                        0.0
                    });
                }
                DesignTerm::Covariate(name) => {
                    let idx = cov_index(&setup.cov_names, name)?;
                    push_encoded_covariate(&mut row, name, &cov_kinds[idx], &vals[idx])?;
                }
            }
        }
        report_rows.push(DesignReportRow {
            comparison: comparison.clone(),
            sample_id: s.sample_id.clone(),
            condition: s.condition.clone().unwrap_or_default(),
            included: true,
            drop_reason: String::new(),
            columns: design_labels.join(","),
            values: row
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(","),
        });
        rows.push((s.sample_id.clone(), row));
    }
    ensure_full_rank(&rows, &design_labels, &comparison)?;

    Ok(OlsDesign {
        rows,
        group_col,
        design_labels,
        report_rows,
        group_is_indicator: true,
    })
}

fn cov_index(names: &[String], name: &str) -> Result<usize> {
    names
        .iter()
        .position(|n| n == name)
        .ok_or_else(|| anyhow!("internal error: missing covariate {:?}", name))
}

pub(super) fn push_encoded_covariate(
    row: &mut Vec<f64>,
    name: &str,
    kind: &CovKind,
    value: &Option<String>,
) -> Result<()> {
    let v = value
        .as_ref()
        .expect("missing covariate filtered out above");
    match kind {
        CovKind::Numeric => {
            let parsed = v
                .parse::<f64>()
                .with_context(|| format!("numeric covariate {:?} contains {:?}", name, v))?;
            row.push(parsed);
        }
        CovKind::Categorical { levels } => {
            for lvl in levels.iter().skip(1) {
                row.push(if v == lvl { 1.0 } else { 0.0 });
            }
        }
    }
    Ok(())
}

fn ensure_full_rank(
    rows: &[(String, Vec<f64>)],
    labels: &[String],
    comparison: &str,
) -> Result<()> {
    let rank = matrix_rank(&rows.iter().map(|(_, row)| row.clone()).collect::<Vec<_>>());
    if rank < labels.len() {
        anyhow::bail!(
            "singular design matrix for {comparison}: rank {} < {} columns ({})",
            rank,
            labels.len(),
            labels.join(",")
        );
    }
    Ok(())
}

fn matrix_rank(rows: &[Vec<f64>]) -> usize {
    if rows.is_empty() {
        return 0;
    }
    let mut a = rows.to_vec();
    let n_rows = a.len();
    let n_cols = a[0].len();
    let mut rank = 0;
    let eps = 1e-10;
    for col in 0..n_cols {
        let pivot = (rank..n_rows).max_by(|&i, &j| {
            a[i][col]
                .abs()
                .partial_cmp(&a[j][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let Some(pivot) = pivot else {
            continue;
        };
        if a[pivot][col].abs() <= eps {
            continue;
        }
        a.swap(rank, pivot);
        let pivot_value = a[rank][col];
        for value in &mut a[rank][col..n_cols] {
            *value /= pivot_value;
        }
        let pivot_tail = a[rank][col..n_cols].to_vec();
        for (r, row) in a.iter_mut().enumerate().take(n_rows) {
            if r == rank {
                continue;
            }
            let factor = row[col];
            for (value, pivot_value) in row[col..n_cols].iter_mut().zip(pivot_tail.iter()) {
                *value -= factor * pivot_value;
            }
        }
        rank += 1;
        if rank == n_rows {
            break;
        }
    }
    rank
}

/// Raw per-group means of the response within the fitted rows. Used to
/// populate `mean_a` / `mean_b` on the `DeResultRow` for schema
/// compatibility with paired-t / welch-t output.
pub(super) fn raw_group_means(
    y: &[f64],
    rows: &[Vec<f64>],
    group_col: usize,
) -> (Option<f64>, Option<f64>) {
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
pub(super) fn ols_to_paired_t_result(
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

/// Build the two-group design matrix for a limma `(a, b)` comparison.
///
/// Rows correspond to non-control samples whose `condition` is `a` or `b`,
/// in the order they appear in the input `samples` slice. Columns are
/// `[intercept, group_b_indicator]` — `group_b_indicator = 1.0` when the
/// sample's condition is `b`, `0.0` when it is `a`. The returned
/// `sample_ids` vector is in the same row order and is used to align
/// abundance lookups in the `y` matrix.
pub(super) fn build_two_group_design(
    samples: &[Sample],
    a: &str,
    b: &str,
) -> (Vec<Vec<f64>>, Vec<String>) {
    let mut design = Vec::new();
    let mut sample_ids = Vec::new();
    for s in samples {
        if s.is_control {
            continue;
        }
        let cond = match s.condition.as_deref() {
            Some(c) => c,
            None => continue,
        };
        let group_b = if cond == a {
            0.0
        } else if cond == b {
            1.0
        } else {
            continue;
        };
        design.push(vec![1.0, group_b]);
        sample_ids.push(s.sample_id.clone());
    }
    (design, sample_ids)
}

/// Augment an already-built `OlsDesign` with external covariate columns from
/// `--adjust-for`. Each covariate is appended as a numeric column to every
/// design row; samples whose `sample_id` is absent from `ext` are dropped
/// (with a drop-reason record). The `group_col` index is unchanged because
/// external columns are appended after the existing columns.
///
/// External covariates are always treated as numeric (the caller has already
/// parsed them as `f64` via `read_external_covariates`).
pub(super) fn augment_ols_design_with_external(
    design: &mut OlsDesign,
    ext: &std::collections::BTreeMap<String, std::collections::BTreeMap<String, f64>>,
    comp_a: &str,
    comp_b: &str,
) -> Result<()> {
    if ext.is_empty() {
        return Ok(());
    }

    // Collect sorted covariate names from the first sample that's in the design.
    // All samples must have the same covariate names (enforced by join_external_covariates).
    let ext_cov_names: Vec<String> = {
        let sample_with_covs = design
            .rows
            .iter()
            .find_map(|(sid, _)| ext.get(sid.as_str()));
        match sample_with_covs {
            Some(covs) => covs.keys().cloned().collect(),
            None => return Ok(()), // no samples in ext — nothing to augment
        }
    };

    if ext_cov_names.is_empty() {
        return Ok(());
    }

    // Drop rows whose sample_id is missing from ext, keep the rest.
    let comparison = format!("{comp_a}-{comp_b}");
    let mut new_rows: Vec<(String, Vec<f64>)> = Vec::new();
    for (sid, mut row) in design.rows.drain(..) {
        match ext.get(&sid) {
            Some(cov_map) => {
                for name in &ext_cov_names {
                    let val = cov_map.get(name).copied().unwrap_or(f64::NAN);
                    row.push(val);
                }
                new_rows.push((sid, row));
            }
            None => {
                // Sample missing from external covariate file: drop with report.
                design.report_rows.push(DesignReportRow {
                    comparison: comparison.clone(),
                    sample_id: sid,
                    condition: String::new(),
                    included: false,
                    drop_reason: "missing_external_covariate_row".to_string(),
                    columns: String::new(),
                    values: String::new(),
                });
            }
        }
    }
    design.rows = new_rows;

    // Append column labels for the external covariates.
    for name in &ext_cov_names {
        design.design_labels.push(name.clone());
    }

    Ok(())
}

/// Paired-t with covariates via paired-difference ANCOVA.
///
/// For each matched subject k: `d_k = abundance_b_k - abundance_a_k`.
/// For each external covariate c: `dc_k = cov_c(sample_b_k) - cov_c(sample_a_k)`.
/// Fits OLS: `d ~ 1 + dc_1 + dc_2 + ...` per protein.
/// The intercept is the adjusted mean paired difference (log_fc).
/// t and p are from the intercept coefficient with df = n_pairs - (1 + n_cov).
///
/// `va` and `vb` are `(subject_id, abundance)` pairs for conditions a and b.
/// `ext` maps sample_id → covariate_name → value.
/// `comp_a` and `comp_b` are the condition labels (used to look up sample_ids).
#[allow(clippy::too_many_arguments)]
pub(super) fn paired_t_covariate_adjusted(
    va: &[(String, f64)],
    vb: &[(String, f64)],
    ext: &std::collections::BTreeMap<String, std::collections::BTreeMap<String, f64>>,
    sample_by_id: &HashMap<&str, &atman_core::Sample>,
    comp_a: &str,
    comp_b: &str,
    min_pairs: usize,
) -> PairedTResult {
    use atman_core::de::{OlsOutcome, SkipReason};

    // Build subject → abundance maps.
    let a_by_subj: HashMap<&str, f64> = va.iter().map(|(s, v)| (s.as_str(), *v)).collect();
    let b_by_subj: HashMap<&str, f64> = vb.iter().map(|(s, v)| (s.as_str(), *v)).collect();

    // Build subject → sample_id lookup for each condition from sample_by_id.
    let mut subj_to_sid_a: HashMap<&str, &str> = HashMap::new();
    let mut subj_to_sid_b: HashMap<&str, &str> = HashMap::new();
    for (sid, s) in sample_by_id.iter() {
        if s.is_control {
            continue;
        }
        let subj = match s.subject_id.as_deref() {
            Some(id) => id,
            None => continue,
        };
        match s.condition.as_deref() {
            Some(c) if c == comp_a => {
                subj_to_sid_a.insert(subj, sid);
            }
            Some(c) if c == comp_b => {
                subj_to_sid_b.insert(subj, sid);
            }
            _ => {}
        }
    }

    // Collect sorted covariate names.
    let ext_cov_names: Vec<String> = {
        let any = ext.values().next();
        match any {
            Some(m) => m.keys().cloned().collect(),
            None => {
                return PairedTResult::Skipped {
                    reason: SkipReason::InsufficientPairs,
                    n_pairs: 0,
                }
            }
        }
    };
    let n_cov = ext_cov_names.len();

    // Build per-subject: (d_k, [dc_k_1, dc_k_2, ...]).
    let mut diffs: Vec<f64> = Vec::new();
    let mut cov_diffs: Vec<Vec<f64>> = Vec::new();

    let subjects: std::collections::BTreeSet<&str> = a_by_subj
        .keys()
        .copied()
        .filter(|s| b_by_subj.contains_key(s))
        .collect();

    for subj in &subjects {
        let abund_a = match a_by_subj.get(subj) {
            Some(v) => *v,
            None => continue,
        };
        let abund_b = match b_by_subj.get(subj) {
            Some(v) => *v,
            None => continue,
        };
        let sid_a = match subj_to_sid_a.get(subj) {
            Some(s) => *s,
            None => continue,
        };
        let sid_b = match subj_to_sid_b.get(subj) {
            Some(s) => *s,
            None => continue,
        };

        let cov_a = match ext.get(sid_a) {
            Some(m) => m,
            None => continue,
        };
        let cov_b = match ext.get(sid_b) {
            Some(m) => m,
            None => continue,
        };

        let mut cds: Vec<f64> = Vec::with_capacity(n_cov);
        let mut valid = true;
        for name in &ext_cov_names {
            let va_cov = match cov_a.get(name) {
                Some(v) => *v,
                None => {
                    valid = false;
                    break;
                }
            };
            let vb_cov = match cov_b.get(name) {
                Some(v) => *v,
                None => {
                    valid = false;
                    break;
                }
            };
            // cov_a is the comp_a sample (e.g. condition "b"), cov_b is comp_b (e.g. condition "a").
            // Covariate diff = comp_a_cov - comp_b_cov, consistent with abundance diff sign.
            let d = va_cov - vb_cov;
            if !d.is_finite() {
                valid = false;
                break;
            }
            cds.push(d);
        }
        if !valid {
            continue;
        }
        if !abund_a.is_finite() || !abund_b.is_finite() {
            continue;
        }

        // diff = comp_a_value - comp_b_value, consistent with regular paired_t
        // which computes diffs as a_i - b_i where a = va[subj] (comp_a condition).
        // With --groups b-a: comp_a = "b", so diff = condition_b - condition_a
        // (positive when b > a, matching atman sign convention).
        diffs.push(abund_a - abund_b);
        cov_diffs.push(cds);
    }

    let n_pairs = diffs.len();
    if n_pairs < min_pairs || n_pairs < 2 {
        return PairedTResult::Skipped {
            reason: SkipReason::InsufficientPairs,
            n_pairs,
        };
    }

    // Build design matrix for OLS: [1, dc_1, dc_2, ...]
    // The intercept (column 0) is the adjusted mean paired difference.
    let p = 1 + n_cov;
    let design: Vec<Vec<f64>> = diffs
        .iter()
        .zip(cov_diffs.iter())
        .map(|(_, cds)| {
            let mut row = Vec::with_capacity(p);
            row.push(1.0);
            row.extend_from_slice(cds);
            row
        })
        .collect();

    // Use atman_core::de::ols() to fit d ~ 1 + cov_diffs.
    // The intercept (column 0) is the adjusted mean diff.
    let fit = match atman_core::de::ols(&design, &diffs, min_pairs) {
        OlsOutcome::Computed(f) => f,
        OlsOutcome::Skipped { reason, n } => {
            return PairedTResult::Skipped { reason, n_pairs: n };
        }
    };

    // mean_a / mean_b are not readily available from the paired-diff design
    // (we only have diffs, not the paired subset's individual abundances).
    // Report NaN; the key outputs (mean_diff, t, p) are fully determined.

    PairedTResult::Computed {
        n_pairs,
        mean_a: f64::NAN,
        mean_b: f64::NAN,
        mean_diff: fit.beta[0], // intercept = adjusted mean diff
        t: fit.t[0],
        df: fit.df,
        p_value: fit.p_value[0],
    }
}
