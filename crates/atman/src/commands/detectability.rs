use anyhow::{anyhow, bail, Context, Result};
use atman_core::de::{bh_fdr, ols, OlsOutcome};
use atman_core::{MeasurementRecord, Sample};
use clap::Args as ClapArgs;
use csv::ReaderBuilder;
use serde_json::json;
use statrs::distribution::{ContinuousCDF, Normal};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::commands::parse_comparisons;
use crate::io::{
    atomic_write, hash_canonical_inputs, need_col, read_measurements_long, read_samples,
    sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Comparison list in A-B form. Detection odds are A relative to B.
    #[arg(long)]
    groups: String,

    /// Covariate formula appended to condition, e.g. "~ age + sex + batch".
    #[arg(long, default_value = "~ 1")]
    design: String,

    /// Output TSV path.
    #[arg(long)]
    output: PathBuf,

    /// Minimum complete samples per condition after covariate filtering.
    #[arg(long, default_value_t = 5)]
    min_samples: usize,

    /// Minimum detected observations required for interpretable abundance testing.
    #[arg(long, default_value_t = 3)]
    min_detected: usize,

    /// Minimum detected observations required in EACH condition before an abundance
    /// effect is reported. Guards against one-sided splits (e.g. 3-vs-0) where the
    /// condition coefficient carries no between-group contrast.
    #[arg(long, default_value_t = 2)]
    min_detected_per_condition: usize,

    /// Minimum non-detected observations required for logistic detection testing.
    #[arg(long, default_value_t = 3)]
    min_missing: usize,

    /// FDR threshold used for the classification column.
    #[arg(long, default_value_t = 0.05)]
    q_threshold: f64,

    /// Logistic IRLS iteration limit.
    #[arg(long, default_value_t = 50)]
    max_iter: usize,

    /// Logistic IRLS convergence tolerance on max coefficient change.
    #[arg(long, default_value_t = 1e-7)]
    tol: f64,
}

#[derive(Debug, Clone)]
struct SampleDesign {
    sample: Sample,
    row: Vec<f64>,
}

#[derive(Debug, Clone)]
enum CovariateSpec {
    Numeric { name: String },
    Categorical { name: String, levels: Vec<String> },
}

#[derive(Debug, Clone)]
struct Row {
    comparison: String,
    panel: String,
    assay_id: String,
    gene_symbol: String,
    n_a: usize,
    n_b: usize,
    n_detected_a: usize,
    n_detected_b: usize,
    n_missing_a: usize,
    n_missing_b: usize,
    n_qc_excluded: usize,
    detection_log_or: Option<f64>,
    detection_or: Option<f64>,
    detection_z: Option<f64>,
    detection_p: Option<f64>,
    detection_q: Option<f64>,
    abundance_effect_detected_only: Option<f64>,
    abundance_t: Option<f64>,
    abundance_p: Option<f64>,
    abundance_q: Option<f64>,
    classification: String,
    skip_reason: String,
}

#[derive(Debug, Clone)]
struct SummaryRow {
    comparison: String,
    n_proteins: usize,
    n_detection_tested: usize,
    n_detection_q_lt_threshold: usize,
    n_abundance_q_lt_threshold: usize,
    n_identical_detection_counts: usize,
    fraction_identical_detection_counts: f64,
    n_zero_detection_delta: usize,
    fraction_zero_detection_delta: f64,
    median_abs_detection_delta: f64,
    max_abs_detection_delta: f64,
    diagnostic: String,
}

pub fn run(args: Args) -> Result<()> {
    if !(args.q_threshold > 0.0 && args.q_threshold < 1.0) {
        bail!("--q-threshold must be between 0 and 1");
    }
    if args.max_iter == 0 {
        bail!("--max-iter must be > 0");
    }
    let started_at = SystemTime::now();
    let comparisons = parse_comparisons(&args.groups)?;
    let terms = parse_design_terms(&args.design)?;
    let samples_path = args.input_dir.join("samples.tsv");
    let samples = read_samples(&samples_path)?;
    let metadata = read_sample_metadata(&samples_path, &terms)?;
    let specs = infer_covariates(&terms, &metadata)?;
    let all_design = build_sample_design(samples, &metadata, &specs)?;
    if all_design.is_empty() {
        bail!(
            "no samples have complete design metadata for {}",
            args.design
        );
    }
    let measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;
    let mut rows = compute_rows(&args, &comparisons, &all_design, &measurements)?;
    apply_q_values(&mut rows, args.q_threshold);
    write_rows(&args.output, &rows)?;
    let summary = summarize_rows(&rows, args.q_threshold);
    let summary_path = summary_path_for(&args.output);
    write_summary(&summary_path, &summary)?;
    for row in &summary {
        if row.diagnostic != "detection_informative" {
            eprintln!(
                "detectability: WARNING {} diagnostic={} identical_counts={:.3} zero_delta={:.3}",
                row.comparison,
                row.diagnostic,
                row.fraction_identical_detection_counts,
                row.fraction_zero_detection_delta
            );
        }
    }

    eprintln!(
        "detectability: comparisons={} proteins={} output={} summary={}",
        comparisons.len(),
        rows.iter()
            .map(|r| (r.assay_id.as_str(), r.gene_symbol.as_str()))
            .collect::<BTreeSet<_>>()
            .len(),
        args.output.display(),
        summary_path.display()
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "detectability",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "groups": args.groups,
            "design": args.design,
            "output": args.output.display().to_string(),
            "summary-output": summary_path.display().to_string(),
            "min-samples": args.min_samples,
            "min-detected": args.min_detected,
            "min-detected-per-condition": args.min_detected_per_condition,
            "min-missing": args.min_missing,
            "q-threshold": args.q_threshold,
            "max-iter": args.max_iter,
            "tol": args.tol,
        }),
        &inputs_sha256,
        &[args.output.clone(), summary_path],
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("detectability: sidecar={}", sidecar.display());
    Ok(())
}

fn compute_rows(
    args: &Args,
    comparisons: &[(String, String)],
    all_design: &[SampleDesign],
    measurements: &[MeasurementRecord],
) -> Result<Vec<Row>> {
    let mut by_protein: BTreeMap<
        (String, String, String, String),
        HashMap<String, MeasurementRecord>,
    > = BTreeMap::new();
    for m in measurements {
        let gene = m
            .gene_symbol
            .clone()
            .unwrap_or_else(|| m.assay_id.0.clone());
        let panel = m.panel.clone().unwrap_or_default();
        by_protein
            .entry((
                m.platform.as_str().to_string(),
                m.assay_id.0.clone(),
                gene,
                panel,
            ))
            .or_default()
            .insert(m.sample_id.clone(), m.clone());
    }

    // Built once; reused for every protein's detection z-test.
    let standard_normal = Normal::new(0.0, 1.0).expect("standard normal");

    let mut out = Vec::new();
    for (a, b) in comparisons {
        let comparison = format!("{a}-{b}");
        let samples: Vec<&SampleDesign> = all_design
            .iter()
            .filter(|sd| {
                matches!(
                    sd.sample.condition.as_deref(),
                    Some(condition) if condition == a || condition == b
                )
            })
            .collect();
        let n_a_total = samples
            .iter()
            .filter(|sd| sd.sample.condition.as_deref() == Some(a.as_str()))
            .count();
        let n_b_total = samples
            .iter()
            .filter(|sd| sd.sample.condition.as_deref() == Some(b.as_str()))
            .count();
        if n_a_total < args.min_samples || n_b_total < args.min_samples {
            bail!(
                "comparison {} has insufficient complete samples after covariate filtering: {}={}, {}={}; min={}",
                comparison,
                a,
                n_a_total,
                b,
                n_b_total,
                args.min_samples
            );
        }

        for ((_platform, assay_id, gene_symbol, panel), per_sample) in &by_protein {
            let mut detection_design = Vec::new();
            let mut y_detect = Vec::new();
            let mut abundance_design = Vec::new();
            let mut y_abundance = Vec::new();
            let mut n_a = 0usize;
            let mut n_b = 0usize;
            let mut n_detected_a = 0usize;
            let mut n_detected_b = 0usize;
            let mut n_missing_a = 0usize;
            let mut n_missing_b = 0usize;
            let mut n_qc_excluded = 0usize;

            for sd in &samples {
                let is_a = sd.sample.condition.as_deref() == Some(a.as_str());
                let mut row = Vec::with_capacity(sd.row.len() + 1);
                row.push(1.0);
                row.push(if is_a { 1.0 } else { 0.0 });
                row.extend(sd.row.iter().skip(1).copied());
                let measurement = per_sample.get(sd.sample.sample_id.as_str());
                if matches!(measurement, Some(m) if m.dropped_by_qc) {
                    n_qc_excluded += 1;
                    continue;
                }
                if is_a {
                    n_a += 1;
                } else {
                    n_b += 1;
                }
                let abundance = measurement.and_then(detected_abundance);
                let detected = abundance.is_some();
                detection_design.push(row.clone());
                y_detect.push(if detected { 1.0 } else { 0.0 });
                match (is_a, detected) {
                    (true, true) => n_detected_a += 1,
                    (false, true) => n_detected_b += 1,
                    (true, false) => n_missing_a += 1,
                    (false, false) => n_missing_b += 1,
                }
                if let Some(value) = abundance {
                    abundance_design.push(row);
                    y_abundance.push(value);
                }
            }

            let mut skip = Vec::new();
            let (detection_log_or, detection_or, detection_z, detection_p) =
                if n_a < args.min_samples || n_b < args.min_samples {
                    skip.push("insufficient_samples");
                    (None, None, None, None)
                } else if n_detected_a + n_detected_b < args.min_detected {
                    skip.push("insufficient_detected");
                    (None, None, None, None)
                } else if n_missing_a + n_missing_b < args.min_missing {
                    skip.push("insufficient_missing");
                    (None, None, None, None)
                } else {
                    match logistic(&detection_design, &y_detect, args.max_iter, args.tol) {
                        Some(fit) if fit.beta.len() > 1 => {
                            let log_or = fit.beta[1];
                            let se = fit.se[1];
                            let z = log_or / se;
                            let p = normal_two_sided_p(&standard_normal, z);
                            (Some(log_or), Some(log_or.exp()), Some(z), Some(p))
                        }
                        _ => {
                            skip.push("detection_model_failed");
                            (None, None, None, None)
                        }
                    }
                };

            let (abundance_effect, abundance_t, abundance_p) =
                if y_abundance.len() < args.min_detected {
                    if !skip.contains(&"insufficient_detected") {
                        skip.push("insufficient_detected_abundance");
                    }
                    (None, None, None)
                } else if n_detected_a < args.min_detected_per_condition
                    || n_detected_b < args.min_detected_per_condition
                {
                    // One-sided detection (e.g. all detected samples in a single
                    // condition): the OLS condition coefficient has no real
                    // between-group contrast even if the rank guard passes.
                    skip.push("one_sided_detection");
                    (None, None, None)
                } else {
                    match ols(&abundance_design, &y_abundance, args.min_detected) {
                        OlsOutcome::Computed(fit) if fit.beta.len() > 1 => (
                            Some(fit.beta[1]),
                            Some(fit.t[1]),
                            finite_opt(fit.p_value[1]),
                        ),
                        _ => {
                            skip.push("abundance_model_failed");
                            (None, None, None)
                        }
                    }
                };

            out.push(Row {
                comparison: comparison.clone(),
                panel: panel.clone(),
                assay_id: assay_id.clone(),
                gene_symbol: gene_symbol.clone(),
                n_a,
                n_b,
                n_detected_a,
                n_detected_b,
                n_missing_a,
                n_missing_b,
                n_qc_excluded,
                detection_log_or,
                detection_or,
                detection_z,
                detection_p,
                detection_q: None,
                abundance_effect_detected_only: abundance_effect,
                abundance_t,
                abundance_p,
                abundance_q: None,
                classification: "pending".to_string(),
                skip_reason: if skip.is_empty() {
                    String::new()
                } else {
                    skip.join(";")
                },
            });
        }
    }
    out.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| a.panel.cmp(&b.panel))
            .then_with(|| a.gene_symbol.cmp(&b.gene_symbol))
            .then_with(|| a.assay_id.cmp(&b.assay_id))
    });
    Ok(out)
}

fn detected_abundance(m: &MeasurementRecord) -> Option<f64> {
    if m.dropped_by_qc || m.below_lod {
        return None;
    }
    let value = m.abundance.as_f64();
    value.is_finite().then_some(value)
}

#[derive(Debug, Clone)]
struct LogisticFit {
    beta: Vec<f64>,
    se: Vec<f64>,
}

fn logistic(design: &[Vec<f64>], y: &[f64], max_iter: usize, tol: f64) -> Option<LogisticFit> {
    let n = design.len();
    if n == 0 || n != y.len() || y.iter().any(|v| *v != 0.0 && *v != 1.0) {
        return None;
    }
    let p = design[0].len();
    if p == 0 || n <= p || design.iter().any(|row| row.len() != p) {
        return None;
    }
    let positives = y.iter().filter(|v| **v == 1.0).count();
    if positives == 0 || positives == n {
        return None;
    }
    let mut beta = vec![0.0; p];
    for _ in 0..max_iter {
        let mut xtwx = vec![vec![0.0; p]; p];
        let mut xtwz = vec![0.0; p];
        for (row, yi) in design.iter().zip(y.iter()) {
            let eta = dot(row, &beta).clamp(-30.0, 30.0);
            let mu = 1.0 / (1.0 + (-eta).exp());
            let w = (mu * (1.0 - mu)).max(1e-9);
            let z = eta + (yi - mu) / w;
            for j in 0..p {
                for k in 0..p {
                    xtwx[j][k] += w * row[j] * row[k];
                }
                xtwz[j] += w * row[j] * z;
            }
        }
        let l = cholesky_lower(&xtwx)?;
        let next = solve_cholesky(&l, &xtwz);
        let max_delta = beta
            .iter()
            .zip(next.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        if next.iter().any(|v| !v.is_finite() || v.abs() > 30.0) {
            return None;
        }
        beta = next;
        if max_delta < tol {
            let cov = inverse_from_cholesky(&l);
            let mut se = Vec::with_capacity(p);
            // `j` indexes the diagonal `cov[j][j]`, not a whole row, so enumerate() does not apply.
            #[allow(clippy::needless_range_loop)]
            for j in 0..p {
                let v = cov[j][j];
                if !v.is_finite() || v <= 0.0 {
                    return None;
                }
                se.push(v.sqrt());
            }
            return Some(LogisticFit { beta, se });
        }
    }
    None
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn cholesky_lower(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    if n == 0 || a.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut l = vec![vec![0.0_f64; n]; n];
    for j in 0..n {
        let diag_sum: f64 = l[j].iter().take(j).map(|v| v * v).sum();
        let diag = a[j][j] - diag_sum;
        if !diag.is_finite() || diag <= 0.0 {
            return None;
        }
        l[j][j] = diag.sqrt();
        for i in (j + 1)..n {
            let off_sum: f64 = l[i]
                .iter()
                .zip(l[j].iter())
                .take(j)
                .map(|(li, lj)| li * lj)
                .sum();
            l[i][j] = (a[i][j] - off_sum) / l[j][j];
        }
    }
    Some(l)
}

fn solve_cholesky(l: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let n = l.len();
    let mut z = vec![0.0; n];
    for i in 0..n {
        let mut sum = rhs[i];
        for k in 0..i {
            sum -= l[i][k] * z[k];
        }
        z[i] = sum / l[i][i];
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = z[i];
        for k in (i + 1)..n {
            sum -= l[k][i] * x[k];
        }
        x[i] = sum / l[i][i];
    }
    x
}

fn inverse_from_cholesky(l: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = l.len();
    let mut out = vec![vec![0.0; n]; n];
    for target in 0..n {
        let mut rhs = vec![0.0; n];
        rhs[target] = 1.0;
        let col = solve_cholesky(l, &rhs);
        for i in 0..n {
            out[i][target] = col[i];
        }
    }
    out
}

fn normal_two_sided_p(dist: &Normal, z: f64) -> f64 {
    if !z.is_finite() {
        return f64::NAN;
    }
    (2.0 * dist.sf(z.abs())).clamp(0.0, 1.0)
}

fn finite_opt(value: f64) -> Option<f64> {
    value.is_finite().then_some(value)
}

fn apply_q_values(rows: &mut [Row], q_threshold: f64) {
    let mut by_comparison: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (idx, row) in rows.iter().enumerate() {
        by_comparison
            .entry(row.comparison.clone())
            .or_default()
            .push(idx);
    }
    for indices in by_comparison.values() {
        let detection_p: Vec<Option<f64>> =
            indices.iter().map(|idx| rows[*idx].detection_p).collect();
        let detection_q = bh_fdr(&detection_p);
        for (idx, q) in indices.iter().zip(detection_q) {
            rows[*idx].detection_q = q;
        }
        let abundance_p: Vec<Option<f64>> =
            indices.iter().map(|idx| rows[*idx].abundance_p).collect();
        let abundance_q = bh_fdr(&abundance_p);
        for (idx, q) in indices.iter().zip(abundance_q) {
            rows[*idx].abundance_q = q;
        }
    }
    for row in rows {
        let detection_sig = row.detection_q.map(|q| q < q_threshold).unwrap_or(false);
        let abundance_sig = row.abundance_q.map(|q| q < q_threshold).unwrap_or(false);
        row.classification = match (detection_sig, abundance_sig) {
            (true, true) => "coupled_shift",
            (true, false) => "detection_shifted",
            (false, true) => "abundance_shifted",
            (false, false) if row.detection_p.is_some() || row.abundance_p.is_some() => {
                "not_significant"
            }
            _ => "uninformative_sparse",
        }
        .to_string();
    }
}

fn summarize_rows(rows: &[Row], q_threshold: f64) -> Vec<SummaryRow> {
    let mut by_comparison: BTreeMap<String, Vec<&Row>> = BTreeMap::new();
    for row in rows {
        by_comparison
            .entry(row.comparison.clone())
            .or_default()
            .push(row);
    }
    let mut out = Vec::new();
    for (comparison, group) in by_comparison {
        let n_proteins = group.len();
        let n_detection_tested = group.iter().filter(|r| r.detection_p.is_some()).count();
        let n_detection_q_lt_threshold = group
            .iter()
            .filter(|r| r.detection_q.map(|q| q < q_threshold).unwrap_or(false))
            .count();
        let n_abundance_q_lt_threshold = group
            .iter()
            .filter(|r| r.abundance_q.map(|q| q < q_threshold).unwrap_or(false))
            .count();
        let n_identical_detection_counts = group
            .iter()
            .filter(|r| r.n_detected_a == r.n_detected_b && r.n_missing_a == r.n_missing_b)
            .count();
        let mut abs_delta = Vec::new();
        for row in &group {
            let pa = if row.n_a == 0 {
                f64::NAN
            } else {
                row.n_detected_a as f64 / row.n_a as f64
            };
            let pb = if row.n_b == 0 {
                f64::NAN
            } else {
                row.n_detected_b as f64 / row.n_b as f64
            };
            if pa.is_finite() && pb.is_finite() {
                abs_delta.push((pa - pb).abs());
            }
        }
        abs_delta.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n_zero_detection_delta = abs_delta.iter().filter(|v| **v < 1e-12).count();
        let fraction_identical_detection_counts =
            fraction(n_identical_detection_counts, n_proteins);
        let fraction_zero_detection_delta = fraction(n_zero_detection_delta, abs_delta.len());
        let median_abs_detection_delta = median_sorted(&abs_delta);
        let max_abs_detection_delta = abs_delta.last().copied().unwrap_or(f64::NAN);
        let diagnostic = if n_detection_tested == 0 {
            "no_detection_tests"
        } else if fraction_identical_detection_counts >= 0.95
            || (fraction_zero_detection_delta >= 0.95 && n_detection_q_lt_threshold == 0)
        {
            "condition_invariant_detection"
        } else if n_detection_q_lt_threshold == 0 && median_abs_detection_delta < 0.02 {
            "weak_detection_contrast"
        } else {
            "detection_informative"
        }
        .to_string();
        out.push(SummaryRow {
            comparison,
            n_proteins,
            n_detection_tested,
            n_detection_q_lt_threshold,
            n_abundance_q_lt_threshold,
            n_identical_detection_counts,
            fraction_identical_detection_counts,
            n_zero_detection_delta,
            fraction_zero_detection_delta,
            median_abs_detection_delta,
            max_abs_detection_delta,
            diagnostic,
        });
    }
    out
}

fn fraction(n: usize, denom: usize) -> f64 {
    if denom == 0 {
        f64::NAN
    } else {
        n as f64 / denom as f64
    }
}

fn median_sorted(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        values[mid]
    } else {
        (values[mid - 1] + values[mid]) / 2.0
    }
}

fn summary_path_for(output: &Path) -> PathBuf {
    let mut path = output.to_path_buf();
    let file_name = output
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("detectability.tsv");
    let summary_name = if let Some(stem) = file_name.strip_suffix(".tsv") {
        format!("{stem}_summary.tsv")
    } else {
        format!("{file_name}.summary.tsv")
    };
    path.set_file_name(summary_name);
    path
}

fn write_summary(path: &Path, rows: &[SummaryRow]) -> Result<()> {
    let mut out = String::from(
        "comparison\tn_proteins\tn_detection_tested\tn_detection_q_lt_threshold\t\
         n_abundance_q_lt_threshold\tn_identical_detection_counts\t\
         fraction_identical_detection_counts\tn_zero_detection_delta\t\
         fraction_zero_detection_delta\tmedian_abs_detection_delta\t\
         max_abs_detection_delta\tdiagnostic\n",
    );
    for row in rows {
        out.push_str(&row.comparison);
        out.push('\t');
        out.push_str(&row.n_proteins.to_string());
        out.push('\t');
        out.push_str(&row.n_detection_tested.to_string());
        out.push('\t');
        out.push_str(&row.n_detection_q_lt_threshold.to_string());
        out.push('\t');
        out.push_str(&row.n_abundance_q_lt_threshold.to_string());
        out.push('\t');
        out.push_str(&row.n_identical_detection_counts.to_string());
        out.push('\t');
        out.push_str(&fmt_f(row.fraction_identical_detection_counts));
        out.push('\t');
        out.push_str(&row.n_zero_detection_delta.to_string());
        out.push('\t');
        out.push_str(&fmt_f(row.fraction_zero_detection_delta));
        out.push('\t');
        out.push_str(&fmt_f(row.median_abs_detection_delta));
        out.push('\t');
        out.push_str(&fmt_f(row.max_abs_detection_delta));
        out.push('\t');
        out.push_str(&row.diagnostic);
        out.push('\n');
    }
    atomic_write(path, out.as_bytes())
}

fn write_rows(path: &Path, rows: &[Row]) -> Result<()> {
    let mut out = String::from(
        "comparison\tpanel\tassay_id\tgene_symbol\tn_a\tn_b\t\
         n_detected_a\tn_detected_b\tn_missing_a\tn_missing_b\tn_qc_excluded\t\
         detection_log_or\tdetection_or\tdetection_z\tdetection_p\tdetection_q\t\
         abundance_effect_detected_only\tabundance_t\tabundance_p\tabundance_q\t\
         classification\tskip_reason\n",
    );
    for row in rows {
        out.push_str(&row.comparison);
        out.push('\t');
        out.push_str(&row.panel);
        out.push('\t');
        out.push_str(&row.assay_id);
        out.push('\t');
        out.push_str(&row.gene_symbol);
        out.push('\t');
        out.push_str(&row.n_a.to_string());
        out.push('\t');
        out.push_str(&row.n_b.to_string());
        out.push('\t');
        out.push_str(&row.n_detected_a.to_string());
        out.push('\t');
        out.push_str(&row.n_detected_b.to_string());
        out.push('\t');
        out.push_str(&row.n_missing_a.to_string());
        out.push('\t');
        out.push_str(&row.n_missing_b.to_string());
        out.push('\t');
        out.push_str(&row.n_qc_excluded.to_string());
        push_opt(&mut out, row.detection_log_or);
        push_opt(&mut out, row.detection_or);
        push_opt(&mut out, row.detection_z);
        push_opt(&mut out, row.detection_p);
        push_opt(&mut out, row.detection_q);
        push_opt(&mut out, row.abundance_effect_detected_only);
        push_opt(&mut out, row.abundance_t);
        push_opt(&mut out, row.abundance_p);
        push_opt(&mut out, row.abundance_q);
        out.push('\t');
        out.push_str(&row.classification);
        out.push('\t');
        out.push_str(&row.skip_reason);
        out.push('\n');
    }
    atomic_write(path, out.as_bytes())
}

fn fmt_f(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.12}")
    } else {
        String::new()
    }
}

fn push_opt(out: &mut String, value: Option<f64>) {
    out.push('\t');
    if let Some(v) = value {
        if v.is_finite() {
            out.push_str(&format!("{v:.12}"));
        }
    }
}

fn parse_design_terms(design: &str) -> Result<Vec<String>> {
    let rhs = design
        .trim()
        .strip_prefix('~')
        .ok_or_else(|| anyhow!("--design must start with `~`"))?
        .trim();
    if rhs.is_empty() || rhs == "1" {
        return Ok(Vec::new());
    }
    let mut terms = Vec::new();
    let mut seen = BTreeSet::new();
    for raw in rhs.split('+') {
        let term = raw.trim();
        if term.is_empty() {
            bail!("empty term in --design");
        }
        if term == "1" {
            continue;
        }
        if term == "0" || term == "-1" {
            bail!("intercept removal is not supported; Atman includes an intercept");
        }
        if term == "condition" {
            bail!("condition is included automatically from --groups; omit it from --design");
        }
        if !seen.insert(term.to_string()) {
            bail!("duplicate term {:?} in --design", term);
        }
        terms.push(term.to_string());
    }
    Ok(terms)
}

fn read_sample_metadata(
    path: &Path,
    terms: &[String],
) -> Result<BTreeMap<String, BTreeMap<String, String>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let sample_col = need_col(&headers, "sample_id", path)?;
    let mut term_cols = Vec::new();
    for term in terms {
        term_cols.push((
            term.clone(),
            need_col(&headers, term, path)
                .with_context(|| format!("required by detectability --design {}", term))?,
        ));
    }
    let mut out = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let mut values = BTreeMap::new();
        for (term, col) in &term_cols {
            values.insert(term.clone(), row[*col].trim().to_string());
        }
        out.insert(row[sample_col].to_string(), values);
    }
    Ok(out)
}

fn infer_covariates(
    terms: &[String],
    metadata: &BTreeMap<String, BTreeMap<String, String>>,
) -> Result<Vec<CovariateSpec>> {
    let mut specs = Vec::new();
    for term in terms {
        let values: Vec<&str> = metadata
            .values()
            .filter_map(|row| row.get(term).map(String::as_str))
            .filter(|value| !value.is_empty())
            .collect();
        if values.is_empty() {
            bail!("covariate {:?} has no non-empty values", term);
        }
        if values.iter().all(|value| value.parse::<f64>().is_ok()) {
            specs.push(CovariateSpec::Numeric { name: term.clone() });
        } else {
            let levels: Vec<String> = values
                .into_iter()
                .map(str::to_string)
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            if levels.len() < 2 {
                bail!("categorical covariate {:?} has fewer than two levels", term);
            }
            specs.push(CovariateSpec::Categorical {
                name: term.clone(),
                levels,
            });
        }
    }
    Ok(specs)
}

fn build_sample_design(
    samples: Vec<Sample>,
    metadata: &BTreeMap<String, BTreeMap<String, String>>,
    specs: &[CovariateSpec],
) -> Result<Vec<SampleDesign>> {
    let mut out = Vec::new();
    for sample in samples {
        let Some(meta) = metadata.get(&sample.sample_id) else {
            continue;
        };
        let mut row = vec![1.0];
        let mut complete = true;
        for spec in specs {
            match spec {
                CovariateSpec::Numeric { name } => {
                    let value = meta.get(name).map(String::as_str).unwrap_or("").trim();
                    if value.is_empty() {
                        complete = false;
                        break;
                    }
                    row.push(value.parse()?);
                }
                CovariateSpec::Categorical { name, levels } => {
                    let value = meta.get(name).map(String::as_str).unwrap_or("").trim();
                    if value.is_empty() {
                        complete = false;
                        break;
                    }
                    if !levels.iter().any(|level| level == value) {
                        bail!("unknown level {:?} for covariate {:?}", value, name);
                    }
                    for level in levels.iter().skip(1) {
                        row.push((value == level) as u8 as f64);
                    }
                }
            }
        }
        if complete {
            out.push(SampleDesign { sample, row });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn logistic_recovers_positive_detection_shift() {
        let mut design = Vec::new();
        let mut y = Vec::new();
        for _ in 0..20 {
            design.push(vec![1.0, 1.0]);
            y.push(1.0);
        }
        for _ in 0..5 {
            design.push(vec![1.0, 1.0]);
            y.push(0.0);
        }
        for _ in 0..5 {
            design.push(vec![1.0, 0.0]);
            y.push(1.0);
        }
        for _ in 0..20 {
            design.push(vec![1.0, 0.0]);
            y.push(0.0);
        }
        let fit = logistic(&design, &y, 50, 1e-8).expect("fit");
        assert!(fit.beta[1] > 2.0, "log OR should be positive: {:?}", fit);
        assert!(fit.se[1].is_finite());
    }

    #[test]
    fn design_rejects_condition_term() {
        let err = parse_design_terms("~ condition + age").expect_err("must fail");
        assert!(err
            .to_string()
            .contains("condition is included automatically"));
    }
}
