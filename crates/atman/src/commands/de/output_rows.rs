//! Sidecar row shapes and their TSV writers (`de_omnibus.tsv`,
//! `de_covariates.tsv`, `de_proxy_summary.tsv`), the per-panel report
//! accumulator, float formatting helpers, and the `--paired-by`
//! subject-id override.

use anyhow::{anyhow, Context, Result};
use atman_core::Sample;
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

use super::median;
use crate::io::atomic_write;

/// Replace each sample's `condition` with the value of `column` in samples.tsv.
pub(super) fn override_condition(
    samples: &mut [Sample],
    samples_path: &Path,
    column: &str,
) -> Result<()> {
    let raw = read_raw_samples(samples_path)?;
    for s in samples.iter_mut() {
        let row = raw.get(&s.sample_id);
        if row.map(|m| !m.contains_key(column)).unwrap_or(true) {
            anyhow::bail!(
                "--condition-col {:?} not found in {:?}",
                column,
                samples_path
            );
        }
        s.condition = row
            .and_then(|m| m.get(column))
            .map(|v| v.trim().to_string())
            .filter(|v| !v.is_empty());
    }
    Ok(())
}

/// Retain samples whose raw samples.tsv row satisfies every predicate.
/// Returns the number of samples dropped.
pub(super) fn apply_subset(
    samples: &mut Vec<Sample>,
    samples_path: &Path,
    preds: &[crate::design::Predicate],
) -> Result<usize> {
    let raw = read_raw_samples(samples_path)?;
    let before = samples.len();
    samples.retain(|s| {
        let row = raw.get(&s.sample_id);
        crate::design::predicates_match(preds, &|col| {
            row.and_then(|m| m.get(col))
                .map(|v| v.trim().to_string())
                .filter(|v| !v.is_empty())
        })
    });
    Ok(before - samples.len())
}

/// Retain samples whose raw samples.tsv row has a non-empty value in every
/// listed column. Returns the number of samples dropped; a column absent
/// from the header is an error.
pub(super) fn require_columns(
    samples: &mut Vec<Sample>,
    samples_path: &Path,
    columns: &[String],
    numeric: bool,
) -> Result<usize> {
    let raw = read_raw_samples(samples_path)?;
    for col in columns {
        if raw
            .values()
            .next()
            .map(|m| !m.contains_key(col))
            .unwrap_or(false)
        {
            anyhow::bail!("--require-cols {:?} not found in {:?}", col, samples_path);
        }
    }
    let before = samples.len();
    let present = |v: &String| -> bool {
        if numeric {
            v.parse::<f64>().map(|x| x.is_finite()).unwrap_or(false)
        } else {
            !v.is_empty()
        }
    };
    samples.retain(|s| {
        raw.get(&s.sample_id)
            .map(|m| {
                columns
                    .iter()
                    .all(|c| m.get(c).map(present).unwrap_or(false))
            })
            .unwrap_or(false)
    });
    Ok(before - samples.len())
}

/// sample_id → {column → raw trimmed value} for every column of samples.tsv.
fn read_raw_samples(path: &Path) -> Result<HashMap<String, BTreeMap<String, String>>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers: Vec<String> = reader
        .headers()?
        .iter()
        .map(|h| h.trim().to_string())
        .collect();
    let sid_col = headers
        .iter()
        .position(|h| h == "sample_id")
        .ok_or_else(|| anyhow!("missing sample_id in {:?}", path))?;
    let mut out = HashMap::new();
    for row in reader.records() {
        let row = row?;
        let mut m = BTreeMap::new();
        for (i, h) in headers.iter().enumerate() {
            m.insert(h.clone(), row.get(i).unwrap_or("").trim().to_string());
        }
        out.insert(row.get(sid_col).unwrap_or("").trim().to_string(), m);
    }
    Ok(out)
}

pub(super) fn override_subject_id(
    samples: &mut [Sample],
    samples_path: &std::path::Path,
    column: &str,
) -> Result<()> {
    use crate::io::need_col;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(samples_path)
        .with_context(|| format!("opening {:?}", samples_path))?;
    let headers = reader.headers()?.clone();
    let c_sample = need_col(&headers, "sample_id", samples_path)?;
    let c_pair = need_col(&headers, column, samples_path)
        .with_context(|| format!("--paired-by column {column:?} not found in samples.tsv"))?;
    let mut by_sample: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    for row in reader.records() {
        let row = row?;
        by_sample.insert(row[c_sample].to_string(), row[c_pair].to_string());
    }
    for s in samples.iter_mut() {
        if let Some(v) = by_sample.get(&s.sample_id) {
            if v.is_empty() {
                s.subject_id = None;
            } else {
                s.subject_id = Some(v.clone());
            }
        }
    }
    Ok(())
}

#[derive(Default)]
pub(super) struct ReportAccumulator {
    pub(super) n_tests: usize,
    pub(super) n_skipped: usize,
    pub(super) n_q_strict: usize,
    pub(super) n_q_relaxed: usize,
    pub(super) min_q: Option<f64>,
    pub(super) max_abs_effect: Option<f64>,
}

/// Row shape for de_omnibus.tsv — per (comparison × protein) omnibus
/// F-test of the user-specified factor's columns.
pub(super) struct OmnibusRow {
    pub(super) panel: String,
    pub(super) assay_id: String,
    pub(super) gene_symbol: String,
    pub(super) factor: String,
    pub(super) comparison: String,
    pub(super) f_statistic: f64,
    pub(super) df_num: usize,
    pub(super) df_den: f64,
    pub(super) p_value: f64,
    pub(super) bh_q: Option<f64>,
}

pub(super) fn write_omnibus_rows(path: &Path, rows: &[OmnibusRow]) -> Result<()> {
    let mut buf = String::from(
        "panel\tassay_id\tgene_symbol\tfactor\tcomparison\tf_statistic\tdf_num\tdf_den\tp_value\tbh_q\n",
    );
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.panel,
            r.assay_id,
            r.gene_symbol,
            r.factor,
            r.comparison,
            format_opt(r.f_statistic),
            r.df_num,
            format_opt(r.df_den),
            format_opt(r.p_value),
            r.bh_q.map(format_opt).unwrap_or_default(),
        ));
    }
    atomic_write(path, buf.as_bytes())
}

pub(super) fn format_opt(v: f64) -> String {
    if !v.is_finite() {
        if v.is_nan() {
            "NaN".into()
        } else if v > 0.0 {
            "Inf".into()
        } else {
            "-Inf".into()
        }
    } else {
        format!("{v}")
    }
}

/// Row shape for de_covariates.tsv — per (comparison × protein × covariate
/// column) estimate, skipping the intercept and the group indicator.
pub(super) struct CovariateRow {
    pub(super) comparison: String,
    pub(super) panel: String,
    pub(super) gene_symbol: String,
    pub(super) covariate: String,
    pub(super) beta: f64,
    pub(super) se: f64,
    pub(super) t: f64,
    pub(super) p_value: f64,
    pub(super) df: f64,
    pub(super) n: usize,
}

#[derive(Clone)]
pub(super) struct DesignReportRow {
    pub(super) comparison: String,
    pub(super) sample_id: String,
    pub(super) condition: String,
    pub(super) included: bool,
    pub(super) drop_reason: String,
    pub(super) columns: String,
    pub(super) values: String,
}

/// Writer for `de_covariates.tsv`. Tab-separated, one row per
/// (comparison × panel × gene × covariate column).
pub(super) fn write_covariate_rows(path: &Path, rows: &[CovariateRow]) -> Result<()> {
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
    atomic_write(path, buf.as_bytes())
}

pub(super) fn write_proxy_summary(
    path: &Path,
    proxy: &str,
    rows: &[CovariateRow],
    p_threshold: f64,
) -> Result<()> {
    let mut by_comparison: BTreeMap<&str, Vec<&CovariateRow>> = BTreeMap::new();
    for row in rows {
        if row.covariate == proxy {
            by_comparison
                .entry(row.comparison.as_str())
                .or_default()
                .push(row);
        }
    }
    let mut buf =
        String::from("proxy\tcomparison\tn_tests\tn_p_strict\tmedian_abs_beta\tmedian_p_value\n");
    for (comparison, rows) in by_comparison {
        let n_tests = rows.len();
        let n_p_strict = rows.iter().filter(|row| row.p_value < p_threshold).count();
        let median_abs_beta = median(rows.iter().map(|row| row.beta.abs()).collect());
        let median_p_value = median(rows.iter().map(|row| row.p_value).collect());
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            proxy,
            comparison,
            n_tests,
            n_p_strict,
            fmt_opt(median_abs_beta),
            fmt_opt(median_p_value),
        ));
    }
    atomic_write(path, buf.as_bytes())
}

pub(super) fn fmt_opt(value: Option<f64>) -> String {
    value
        .filter(|v| v.is_finite())
        .map(|v| format!("{v:.6}"))
        .unwrap_or_else(|| "NA".to_string())
}
