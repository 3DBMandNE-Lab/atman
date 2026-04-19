//! Null calibration command.

use anyhow::{bail, Context, Result};
use atman_core::de::{bh_fdr, paired_t, welch_t, PairedTResult, SkipReason};
use atman_core::stats::mean;
use atman_core::Sample;
use clap::Args as ClapArgs;
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use super::parse_comparisons;
use crate::io::{atomic_write, read_measurements_long, read_proteins, read_samples};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for null_summary.tsv and empirical_p.tsv.
    #[arg(long)]
    output_dir: PathBuf,

    /// Comma-separated comparisons in A-B form.
    #[arg(long)]
    groups: String,

    /// Null mode: welch-t uses unpaired label permutation; paired-t uses sign flips.
    #[arg(long, default_value = "welch-t")]
    test: String,

    /// Number of null resamples.
    #[arg(long, default_value_t = 1000)]
    n: usize,

    /// Minimum subjects per group, or matched subjects in paired mode.
    #[arg(long, default_value_t = 2)]
    min_pairs: usize,

    /// Deterministic RNG seed.
    #[arg(long, default_value_t = 1)]
    seed: u64,

    /// Comma-separated BH-q thresholds to summarize.
    #[arg(long, default_value = "0.05,0.10")]
    q_thresholds: String,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ProteinKey {
    panel: String,
    gene_symbol: String,
}

#[derive(Debug, Clone)]
struct ProteinMeta {
    assay_id: String,
    uniprot: String,
}

#[derive(Debug, Clone)]
struct ProteinCells {
    key: ProteinKey,
    meta: ProteinMeta,
    by_condition: BTreeMap<String, BTreeMap<String, Vec<f64>>>,
}

#[derive(Debug)]
struct EmpiricalRow {
    comparison: String,
    panel: String,
    assay_id: String,
    gene_symbol: String,
    uniprot: String,
    n_pairs: usize,
    mean_diff: Option<f64>,
    t: Option<f64>,
    p_value: Option<f64>,
    bh_q: Option<f64>,
    empirical_p: Option<f64>,
    skip_reason: String,
}

#[derive(Debug)]
struct SummaryRow {
    comparison: String,
    test: String,
    q_threshold: f64,
    observed_hits: usize,
    mean_null_hits: f64,
    p95_null_hits: usize,
    empirical_fdr: f64,
    n_permutations: usize,
    seed: u64,
}

pub fn run(args: Args) -> Result<()> {
    if args.test != "welch-t" && args.test != "paired-t" {
        bail!("null supports --test welch-t or paired-t");
    }
    if args.n == 0 {
        bail!("--n must be > 0");
    }
    if args.min_pairs < 2 {
        bail!("--min-pairs must be >= 2");
    }
    let q_thresholds = parse_q_thresholds(&args.q_thresholds)?;
    let comparisons = parse_comparisons(&args.groups)?;
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    validate_design(&samples, &comparisons, &args.test, args.min_pairs)?;
    let proteins = read_proteins(&args.input_dir.join("proteins.tsv"))?;
    let measurements_path = if args.input_dir.join("qc_measurements.tsv").exists() {
        args.input_dir.join("qc_measurements.tsv")
    } else {
        args.input_dir.join("measurements.tsv")
    };
    let measurements = read_measurements_long(&measurements_path)?;
    let protein_cells = build_cells(&samples, &proteins, &measurements);

    let mut rng = Rng64::new(args.seed);
    let mut summary_rows = Vec::new();
    let mut empirical_rows = Vec::new();

    for (a, b) in &comparisons {
        let comparison = format!("{a}-{b}");
        let mut observed_p = Vec::with_capacity(protein_cells.len());
        let mut observed_rows = Vec::with_capacity(protein_cells.len());
        for protein in &protein_cells {
            let result = observed_test(protein, a, b, &args.test, args.min_pairs);
            observed_p.push(result.p_value);
            observed_rows.push(result);
        }
        let observed_q = bh_fdr(&observed_p);
        for (row, q) in observed_rows.iter_mut().zip(observed_q.into_iter()) {
            row.bh_q = q;
        }

        let mut extreme_counts = vec![0usize; protein_cells.len()];
        let mut null_hits_by_threshold = vec![Vec::with_capacity(args.n); q_thresholds.len()];
        for _ in 0..args.n {
            let mut null_p = Vec::with_capacity(protein_cells.len());
            let mut null_abs_t = Vec::with_capacity(protein_cells.len());
            for protein in &protein_cells {
                let null = null_test(protein, a, b, &args.test, args.min_pairs, &mut rng);
                null_p.push(null.p_value);
                null_abs_t.push(null.t.map(f64::abs));
            }
            let null_q = bh_fdr(&null_p);
            for (idx, row) in observed_rows.iter().enumerate() {
                if let (Some(obs_t), Some(null_t)) = (row.t.map(f64::abs), null_abs_t[idx]) {
                    if null_t >= obs_t {
                        extreme_counts[idx] += 1;
                    }
                }
            }
            for (threshold_idx, threshold) in q_thresholds.iter().enumerate() {
                let hits = null_q
                    .iter()
                    .filter(|q| q.map(|v| v < *threshold).unwrap_or(false))
                    .count();
                null_hits_by_threshold[threshold_idx].push(hits);
            }
        }

        for (idx, row) in observed_rows.iter_mut().enumerate() {
            if row.t.is_some() {
                row.empirical_p = Some((extreme_counts[idx] + 1) as f64 / (args.n + 1) as f64);
            }
        }
        for threshold_idx in 0..q_thresholds.len() {
            let threshold = q_thresholds[threshold_idx];
            let observed_hits = observed_rows
                .iter()
                .filter(|r| r.bh_q.map(|q| q < threshold).unwrap_or(false))
                .count();
            let mut null_hits = null_hits_by_threshold[threshold_idx].clone();
            null_hits.sort_unstable();
            let mean_null_hits = null_hits.iter().sum::<usize>() as f64 / null_hits.len() as f64;
            let p95_null_hits = quantile_usize_sorted(&null_hits, 0.95);
            let empirical_fdr = if observed_hits == 0 {
                0.0
            } else {
                (mean_null_hits / observed_hits as f64).min(1.0)
            };
            summary_rows.push(SummaryRow {
                comparison: comparison.clone(),
                test: args.test.clone(),
                q_threshold: threshold,
                observed_hits,
                mean_null_hits,
                p95_null_hits,
                empirical_fdr,
                n_permutations: args.n,
                seed: args.seed,
            });
        }
        empirical_rows.extend(observed_rows);
    }

    empirical_rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| match (a.bh_q, b.bh_q) {
                (Some(qa), Some(qb)) => qa.partial_cmp(&qb).unwrap_or(std::cmp::Ordering::Equal),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            })
            .then_with(|| a.panel.cmp(&b.panel))
            .then_with(|| a.gene_symbol.cmp(&b.gene_symbol))
    });

    write_summary(&args.output_dir.join("null_summary.tsv"), &summary_rows)?;
    write_empirical(&args.output_dir.join("empirical_p.tsv"), &empirical_rows)?;
    eprintln!(
        "null: test={} comparisons={} proteins={} n={} seed={}",
        args.test,
        comparisons.len(),
        protein_cells.len(),
        args.n,
        args.seed
    );
    Ok(())
}

fn validate_design(
    samples: &[Sample],
    comparisons: &[(String, String)],
    test: &str,
    min_pairs: usize,
) -> Result<()> {
    for (a, b) in comparisons {
        let a_subjects = subjects_for_condition(samples, a);
        let b_subjects = subjects_for_condition(samples, b);
        if test == "welch-t" {
            if a_subjects.len() < min_pairs || b_subjects.len() < min_pairs {
                bail!(
                    "comparison {a}-{b} has insufficient subjects for welch-t: {} vs {}, min {}",
                    a_subjects.len(),
                    b_subjects.len(),
                    min_pairs
                );
            }
        } else {
            let matched = a_subjects.intersection(&b_subjects).count();
            if matched < min_pairs {
                bail!(
                    "comparison {a}-{b} has insufficient matched subjects for paired-t: {}, min {}",
                    matched,
                    min_pairs
                );
            }
        }
    }
    Ok(())
}

fn subjects_for_condition(samples: &[Sample], condition: &str) -> BTreeSet<String> {
    samples
        .iter()
        .filter(|s| !s.is_control && s.condition.as_deref() == Some(condition))
        .map(|s| s.subject_id.as_ref().unwrap_or(&s.sample_id).clone())
        .collect()
}

fn build_cells(
    samples: &[Sample],
    proteins: &[atman_core::ProteinIdentity],
    measurements: &[atman_core::MeasurementRecord],
) -> Vec<ProteinCells> {
    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();
    let mut meta: BTreeMap<ProteinKey, ProteinMeta> = BTreeMap::new();
    for p in proteins {
        if let (Some(panel), Some(gene)) = (&p.panel, &p.gene_symbol) {
            meta.entry(ProteinKey {
                panel: panel.clone(),
                gene_symbol: gene.clone(),
            })
            .or_insert_with(|| ProteinMeta {
                assay_id: p.assay_id.0.clone(),
                uniprot: p.uniprot.join(","),
            });
        }
    }

    let mut cells: BTreeMap<ProteinKey, BTreeMap<String, BTreeMap<String, Vec<f64>>>> =
        BTreeMap::new();
    for m in measurements {
        let Some(value) = m.effective_abundance() else {
            continue;
        };
        let Some(sample) = sample_by_id.get(m.sample_id.as_str()) else {
            continue;
        };
        if sample.is_control {
            continue;
        }
        let Some(condition) = &sample.condition else {
            continue;
        };
        let Some(panel) = &m.panel else {
            continue;
        };
        let Some(gene) = &m.gene_symbol else {
            continue;
        };
        let subject = sample.subject_id.as_ref().unwrap_or(&sample.sample_id);
        cells
            .entry(ProteinKey {
                panel: panel.clone(),
                gene_symbol: gene.clone(),
            })
            .or_default()
            .entry(condition.clone())
            .or_default()
            .entry(subject.clone())
            .or_default()
            .push(value);
    }

    cells
        .into_iter()
        .map(|(key, by_condition)| {
            let meta = meta.get(&key).cloned().unwrap_or_else(|| ProteinMeta {
                assay_id: String::new(),
                uniprot: String::new(),
            });
            ProteinCells {
                key,
                meta,
                by_condition,
            }
        })
        .collect()
}

fn observed_test(
    protein: &ProteinCells,
    a: &str,
    b: &str,
    test: &str,
    min_pairs: usize,
) -> EmpiricalRow {
    let result = if test == "paired-t" {
        let pairs = paired_values(protein, a, b);
        paired_t(&pairs, min_pairs)
    } else {
        let a_values = subject_mean_values(protein.by_condition.get(a));
        let b_values = subject_mean_values(protein.by_condition.get(b));
        welch_t(&a_values, &b_values, min_pairs)
    };
    result_to_empirical_row(protein, format!("{a}-{b}"), result)
}

fn null_test(
    protein: &ProteinCells,
    a: &str,
    b: &str,
    test: &str,
    min_pairs: usize,
    rng: &mut Rng64,
) -> TestStat {
    if test == "paired-t" {
        let pairs = paired_values(protein, a, b);
        if pairs.len() < min_pairs {
            return TestStat::skipped();
        }
        let diffs: Vec<f64> = pairs
            .iter()
            .map(|(av, bv)| {
                let diff = av - bv;
                if rng.gen_bool() {
                    diff
                } else {
                    -diff
                }
            })
            .collect();
        t_from_diffs(&diffs, min_pairs)
    } else {
        let a_values = subject_mean_values(protein.by_condition.get(a));
        let b_values = subject_mean_values(protein.by_condition.get(b));
        if a_values.len() < min_pairs || b_values.len() < min_pairs {
            return TestStat::skipped();
        }
        let na = a_values.len();
        let mut pooled = Vec::with_capacity(a_values.len() + b_values.len());
        pooled.extend(a_values);
        pooled.extend(b_values);
        shuffle(&mut pooled, rng);
        let (null_a, null_b) = pooled.split_at(na);
        result_to_stat(welch_t(null_a, null_b, min_pairs))
    }
}

fn paired_values(protein: &ProteinCells, a: &str, b: &str) -> Vec<(f64, f64)> {
    let a_values = subject_mean_map(protein.by_condition.get(a));
    let b_values = subject_mean_map(protein.by_condition.get(b));
    a_values
        .iter()
        .filter_map(|(subject, av)| b_values.get(subject).map(|bv| (*av, *bv)))
        .collect()
}

fn subject_mean_values(input: Option<&BTreeMap<String, Vec<f64>>>) -> Vec<f64> {
    input
        .into_iter()
        .flat_map(|m| m.values())
        .map(|values| mean(values))
        .collect()
}

fn subject_mean_map(input: Option<&BTreeMap<String, Vec<f64>>>) -> BTreeMap<String, f64> {
    input
        .into_iter()
        .flat_map(|m| m.iter())
        .map(|(subject, values)| (subject.clone(), mean(values)))
        .collect()
}

fn result_to_empirical_row(
    protein: &ProteinCells,
    comparison: String,
    result: PairedTResult,
) -> EmpiricalRow {
    let stat = result_to_stat(result.clone());
    let (n_pairs, mean_diff, skip_reason) = match result {
        PairedTResult::Computed {
            n_pairs, mean_diff, ..
        } => (n_pairs, Some(mean_diff), String::new()),
        PairedTResult::Skipped { reason, n_pairs } => {
            (n_pairs, None, skip_reason_to_string(reason))
        }
    };
    EmpiricalRow {
        comparison,
        panel: protein.key.panel.clone(),
        assay_id: protein.meta.assay_id.clone(),
        gene_symbol: protein.key.gene_symbol.clone(),
        uniprot: protein.meta.uniprot.clone(),
        n_pairs,
        mean_diff,
        t: stat.t,
        p_value: stat.p_value,
        bh_q: None,
        empirical_p: None,
        skip_reason,
    }
}

#[derive(Debug, Clone, Copy)]
struct TestStat {
    t: Option<f64>,
    p_value: Option<f64>,
}

impl TestStat {
    fn skipped() -> Self {
        Self {
            t: None,
            p_value: None,
        }
    }
}

fn result_to_stat(result: PairedTResult) -> TestStat {
    match result {
        PairedTResult::Computed { t, p_value, .. } => TestStat {
            t: Some(t),
            p_value: Some(p_value),
        },
        PairedTResult::Skipped { .. } => TestStat::skipped(),
    }
}

fn t_from_diffs(diffs: &[f64], min_pairs: usize) -> TestStat {
    if diffs.len() < min_pairs || diffs.len() < 2 {
        return TestStat::skipped();
    }
    if diffs.iter().any(|v| !v.is_finite()) {
        return TestStat::skipped();
    }
    let n = diffs.len() as f64;
    let mean_diff = mean(diffs);
    let var = diffs.iter().map(|d| (d - mean_diff).powi(2)).sum::<f64>() / (n - 1.0);
    if !var.is_finite() || var == 0.0 {
        return TestStat::skipped();
    }
    let t = mean_diff / (var.sqrt() / n.sqrt());
    let df = n - 1.0;
    if !t.is_finite() || df <= 0.0 {
        return TestStat::skipped();
    }
    let Ok(dist) = StudentsT::new(0.0, 1.0, df) else {
        return TestStat::skipped();
    };
    TestStat {
        t: Some(t),
        p_value: Some(2.0 * (1.0 - dist.cdf(t.abs()))),
    }
}

fn skip_reason_to_string(reason: SkipReason) -> String {
    match reason {
        SkipReason::InsufficientPairs => "insufficient_pairs".to_string(),
        SkipReason::ZeroVariance => "zero_variance".to_string(),
        SkipReason::NonFiniteInput => "non_finite_input".to_string(),
    }
}

fn parse_q_thresholds(input: &str) -> Result<Vec<f64>> {
    let mut thresholds = Vec::new();
    for raw in input.split(',') {
        let value: f64 = raw
            .trim()
            .parse()
            .with_context(|| format!("invalid q threshold {:?}", raw))?;
        if !value.is_finite() || value <= 0.0 || value > 1.0 {
            bail!("q thresholds must satisfy 0 < q <= 1");
        }
        thresholds.push(value);
    }
    if thresholds.is_empty() {
        bail!("at least one q threshold is required");
    }
    Ok(thresholds)
}

fn write_summary(path: &Path, rows: &[SummaryRow]) -> Result<()> {
    let mut buf = String::from(
        "comparison\ttest\tq_threshold\tobserved_hits\tmean_null_hits\tp95_null_hits\tempirical_fdr\tn_permutations\tseed\n",
    );
    for r in rows {
        buf.push_str(&r.comparison);
        buf.push('\t');
        buf.push_str(&r.test);
        buf.push('\t');
        buf.push_str(&r.q_threshold.to_string());
        buf.push('\t');
        buf.push_str(&r.observed_hits.to_string());
        buf.push('\t');
        buf.push_str(&r.mean_null_hits.to_string());
        buf.push('\t');
        buf.push_str(&r.p95_null_hits.to_string());
        buf.push('\t');
        buf.push_str(&r.empirical_fdr.to_string());
        buf.push('\t');
        buf.push_str(&r.n_permutations.to_string());
        buf.push('\t');
        buf.push_str(&r.seed.to_string());
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn write_empirical(path: &Path, rows: &[EmpiricalRow]) -> Result<()> {
    let mut buf = String::from(
        "comparison\tpanel\tassay_id\tgene_symbol\tuniprot\tn_pairs\tmean_diff\tt\tp_value\tbh_q\tempirical_p\tskip_reason\n",
    );
    for r in rows {
        buf.push_str(&r.comparison);
        buf.push('\t');
        buf.push_str(&r.panel);
        buf.push('\t');
        buf.push_str(&r.assay_id);
        buf.push('\t');
        buf.push_str(&r.gene_symbol);
        buf.push('\t');
        buf.push_str(&r.uniprot);
        buf.push('\t');
        buf.push_str(&r.n_pairs.to_string());
        buf.push('\t');
        push_opt(&mut buf, r.mean_diff);
        buf.push('\t');
        push_opt(&mut buf, r.t);
        buf.push('\t');
        push_opt(&mut buf, r.p_value);
        buf.push('\t');
        push_opt(&mut buf, r.bh_q);
        buf.push('\t');
        push_opt(&mut buf, r.empirical_p);
        buf.push('\t');
        buf.push_str(&r.skip_reason);
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn push_opt(buf: &mut String, value: Option<f64>) {
    if let Some(value) = value {
        buf.push_str(&value.to_string());
    }
}


fn quantile_usize_sorted(values: &[usize], q: f64) -> usize {
    let idx = ((values.len() - 1) as f64 * q).round() as usize;
    values[idx.min(values.len() - 1)]
}

fn shuffle<T>(values: &mut [T], rng: &mut Rng64) {
    for i in (1..values.len()).rev() {
        let j = rng.gen_range(i + 1);
        values.swap(i, j);
    }
}

struct Rng64 {
    state: u64,
}

impl Rng64 {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 0x9e3779b97f4a7c15 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545f4914f6cdd1d)
    }

    fn gen_bool(&mut self) -> bool {
        self.next_u64() & 1 == 1
    }

    fn gen_range(&mut self, upper: usize) -> usize {
        (self.next_u64() as usize) % upper
    }
}
