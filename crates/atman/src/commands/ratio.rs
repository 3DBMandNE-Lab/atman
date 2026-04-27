use anyhow::{bail, Result};
use atman_core::stats::mean;
use atman_core::Sample;
use clap::{Args as ClapArgs, ValueEnum};
use serde_json::json;
use statrs::distribution::{ContinuousCDF, Normal, StudentsT};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::SystemTime;

use crate::io::{
    atomic_write, hash_canonical_inputs, read_measurements_long, read_samples, sidecar_path_for,
    write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Numerator quantity, e.g. modules.intrathecal_IgV or intrathecal_IgV.
    #[arg(long)]
    numerator: String,

    /// Denominator quantity, e.g. modules.plasma_IgG or plasma_IgG.
    #[arg(long)]
    denominator: String,

    /// Group comparison encoded as GroupA-GroupB.
    #[arg(long)]
    groups: String,

    /// Ratio test.
    #[arg(long, value_enum, default_value_t = Test::Wilcoxon)]
    test: Test,

    /// Bootstrap replicates for median-difference CI.
    #[arg(long, default_value_t = 1000)]
    n_bootstrap: usize,

    /// RNG seed for bootstrap/permutation.
    #[arg(long, default_value_t = 20260418)]
    seed: u64,

    /// Permutations used when --test permutation.
    #[arg(long, default_value_t = 1000)]
    n_permutations: usize,

    /// Output ratio result TSV.
    #[arg(long)]
    output: PathBuf,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum Test {
    Wilcoxon,
    TTest,
    Permutation,
}

#[derive(Debug, Clone)]
struct SubjectRatio {
    condition: String,
    ratio: f64,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    if args.n_bootstrap == 0 {
        bail!("--n-bootstrap must be at least 1");
    }
    if args.test == Test::Permutation && args.n_permutations == 0 {
        bail!("--n-permutations must be at least 1");
    }
    let (group_a, group_b) = parse_groups(&args.groups)?;
    let numerator = normalize_quantity(&args.numerator);
    let denominator = normalize_quantity(&args.denominator);
    if numerator == denominator {
        bail!("--numerator and --denominator resolve to the same quantity");
    }

    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let measurements_path = args.input_dir.join("measurements.tsv");
    let measurements = read_measurements_long(&measurements_path)?;
    let ratios = subject_ratios(&samples, &measurements, &numerator, &denominator)?;
    let a: Vec<f64> = ratios
        .iter()
        .filter(|row| row.condition == group_a)
        .map(|row| row.ratio)
        .collect();
    let b: Vec<f64> = ratios
        .iter()
        .filter(|row| row.condition == group_b)
        .map(|row| row.ratio)
        .collect();
    if a.len() < 2 || b.len() < 2 {
        bail!(
            "ratio comparison {}-{} has n_a={} n_b={}; need at least 2 per group",
            group_a,
            group_b,
            a.len(),
            b.len()
        );
    }

    let median_a = median(a.clone()).expect("non-empty");
    let median_b = median(b.clone()).expect("non-empty");
    let point = median_a - median_b;
    let mut rng = Rng64::new(args.seed);
    let mut effects = bootstrap_median_diffs(&a, &b, args.n_bootstrap, &mut rng);
    effects.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let ci_low = percentile_sorted(&effects, 0.025);
    let ci_high = percentile_sorted(&effects, 0.975);
    let sign_stability = sign_stability(point, &effects);
    let p_value = match args.test {
        Test::Wilcoxon => wilcoxon_rank_sum_p(&a, &b),
        Test::TTest => welch_t_p(&a, &b),
        Test::Permutation => permutation_p(&a, &b, args.n_permutations, &mut rng),
    };

    let mut out = String::from(
        "comparison\tnumerator\tdenominator\ttest\tn_subjects\tn_a\tn_b\tmedian_a\tmedian_b\tmedian_diff\tci_low\tci_high\tp_value\tsign_stability\tn_bootstrap\n",
    );
    out.push_str(&format!(
        "{}-{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        group_a,
        group_b,
        numerator,
        denominator,
        test_name(args.test),
        ratios.len(),
        a.len(),
        b.len(),
        fmt(median_a),
        fmt(median_b),
        fmt(point),
        fmt(ci_low),
        fmt(ci_high),
        fmt(p_value),
        fmt(sign_stability),
        args.n_bootstrap
    ));
    atomic_write(&args.output, out.as_bytes())?;
    eprintln!(
        "ratio: comparison={}-{} numerator={} denominator={} n_subjects={} output={}",
        group_a,
        group_b,
        numerator,
        denominator,
        ratios.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &[
            "measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
    )?;
    let sidecar = sidecar_path_for(&args.output);
    let test_label = test_name(args.test);
    write_run_sidecar(
        &sidecar,
        "ratio",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "numerator": args.numerator,
            "denominator": args.denominator,
            "groups": args.groups,
            "test": test_label,
            "n-bootstrap": args.n_bootstrap,
            "seed": args.seed,
            "n-permutations": args.n_permutations,
            "output": args.output.display().to_string(),
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("ratio: sidecar={}", sidecar.display());
    Ok(())
}

fn parse_groups(groups: &str) -> Result<(String, String)> {
    let Some((a, b)) = groups.split_once('-') else {
        bail!("--groups must be encoded as GroupA-GroupB");
    };
    let a = a.trim();
    let b = b.trim();
    if a.is_empty() || b.is_empty() || a == b {
        bail!("--groups must contain two distinct non-empty group names");
    }
    Ok((a.to_string(), b.to_string()))
}

fn normalize_quantity(value: &str) -> String {
    value
        .strip_prefix("modules.")
        .unwrap_or(value)
        .trim()
        .to_string()
}

fn subject_ratios(
    samples: &[Sample],
    measurements: &[atman_core::MeasurementRecord],
    numerator: &str,
    denominator: &str,
) -> Result<Vec<SubjectRatio>> {
    let mut sample_meta: BTreeMap<&str, (&str, &str)> = BTreeMap::new();
    for sample in samples {
        let Some(subject) = sample.subject_id.as_deref() else {
            continue;
        };
        let Some(condition) = sample.condition.as_deref() else {
            continue;
        };
        sample_meta.insert(sample.sample_id.as_str(), (subject, condition));
    }

    let mut per_sample: BTreeMap<&str, (Vec<f64>, Vec<f64>)> = BTreeMap::new();
    for measurement in measurements {
        let Some(value) = measurement.effective_abundance() else {
            continue;
        };
        if !sample_meta.contains_key(measurement.sample_id.as_str()) {
            continue;
        }
        let is_num = measurement_matches(measurement, numerator);
        let is_den = measurement_matches(measurement, denominator);
        if !is_num && !is_den {
            continue;
        }
        let entry = per_sample
            .entry(measurement.sample_id.as_str())
            .or_default();
        if is_num {
            entry.0.push(value);
        }
        if is_den {
            entry.1.push(value);
        }
    }

    let mut per_subject: BTreeMap<&str, (&str, Vec<f64>)> = BTreeMap::new();
    for (sample_id, (num, den)) in per_sample {
        if num.is_empty() || den.is_empty() {
            continue;
        }
        let (subject, condition) = sample_meta
            .get(sample_id)
            .copied()
            .expect("sample filtered above");
        let ratio = mean(&num) - mean(&den);
        let entry = per_subject
            .entry(subject)
            .or_insert((condition, Vec::new()));
        if entry.0 != condition {
            bail!(
                "subject {:?} has conflicting conditions in samples.tsv",
                subject
            );
        }
        entry.1.push(ratio);
    }

    let rows: Vec<SubjectRatio> = per_subject
        .into_iter()
        .map(|(_subject, (condition, ratios))| SubjectRatio {
            condition: condition.to_string(),
            ratio: mean(&ratios),
        })
        .collect();
    if rows.is_empty() {
        bail!(
            "no complete subject ratios found for numerator {:?} and denominator {:?}",
            numerator,
            denominator
        );
    }
    Ok(rows)
}

fn measurement_matches(measurement: &atman_core::MeasurementRecord, target: &str) -> bool {
    measurement.assay_id.0 == target
        || measurement.gene_symbol.as_deref() == Some(target)
        || format!("modules.{}", measurement.assay_id.0) == target
        || measurement
            .gene_symbol
            .as_deref()
            .map(|gene| format!("modules.{gene}") == target)
            .unwrap_or(false)
}

fn bootstrap_median_diffs(a: &[f64], b: &[f64], n: usize, rng: &mut Rng64) -> Vec<f64> {
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let mut aa = Vec::with_capacity(a.len());
        let mut bb = Vec::with_capacity(b.len());
        for _ in 0..a.len() {
            aa.push(a[rng.gen_range(a.len())]);
        }
        for _ in 0..b.len() {
            bb.push(b[rng.gen_range(b.len())]);
        }
        out.push(median(aa).expect("resampled") - median(bb).expect("resampled"));
    }
    out
}

fn permutation_p(a: &[f64], b: &[f64], n: usize, rng: &mut Rng64) -> f64 {
    let observed = median(a.to_vec()).expect("a") - median(b.to_vec()).expect("b");
    let mut combined = a.to_vec();
    combined.extend_from_slice(b);
    let n_a = a.len();
    let mut extreme = 0usize;
    for _ in 0..n {
        shuffle(&mut combined, rng);
        let perm_a = median(combined[..n_a].to_vec()).expect("perm a");
        let perm_b = median(combined[n_a..].to_vec()).expect("perm b");
        if (perm_a - perm_b).abs() >= observed.abs() {
            extreme += 1;
        }
    }
    (extreme as f64 + 1.0) / (n as f64 + 1.0)
}

fn wilcoxon_rank_sum_p(a: &[f64], b: &[f64]) -> f64 {
    if a.len() < 2 || b.len() < 2 {
        return f64::NAN;
    }
    let mut pooled: Vec<(f64, usize)> = a
        .iter()
        .copied()
        .map(|v| (v, 0))
        .chain(b.iter().copied().map(|v| (v, 1)))
        .collect();
    pooled.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut rank_sum_a = 0.0;
    let mut tie_term = 0.0;
    let mut i = 0usize;
    while i < pooled.len() {
        let start = i;
        let value = pooled[i].0;
        while i < pooled.len() && pooled[i].0 == value {
            i += 1;
        }
        let end = i;
        let avg_rank = (start + 1 + end) as f64 / 2.0;
        for item in pooled.iter().take(end).skip(start) {
            if item.1 == 0 {
                rank_sum_a += avg_rank;
            }
        }
        let tie = (end - start) as f64;
        tie_term += tie.powi(3) - tie;
    }
    let n_a = a.len() as f64;
    let n_b = b.len() as f64;
    let n_total = n_a + n_b;
    let u_a = rank_sum_a - n_a * (n_a + 1.0) / 2.0;
    let mean_u = n_a * n_b / 2.0;
    let var_u = n_a * n_b / 12.0 * ((n_total + 1.0) - tie_term / (n_total * (n_total - 1.0)));
    if var_u <= 0.0 {
        return f64::NAN;
    }
    let z = (u_a - mean_u).abs() / var_u.sqrt();
    let normal = Normal::new(0.0, 1.0).expect("standard normal");
    2.0 * (1.0 - normal.cdf(z))
}

fn welch_t_p(a: &[f64], b: &[f64]) -> f64 {
    let mean_a = mean(a);
    let mean_b = mean(b);
    let var_a = sample_var(a, mean_a);
    let var_b = sample_var(b, mean_b);
    let n_a = a.len() as f64;
    let n_b = b.len() as f64;
    let se2 = var_a / n_a + var_b / n_b;
    if se2 <= 0.0 {
        return f64::NAN;
    }
    let t = (mean_a - mean_b) / se2.sqrt();
    let df_num = se2.powi(2);
    let df_den = (var_a / n_a).powi(2) / (n_a - 1.0) + (var_b / n_b).powi(2) / (n_b - 1.0);
    if df_den <= 0.0 {
        return f64::NAN;
    }
    let dist = StudentsT::new(0.0, 1.0, df_num / df_den).expect("student t");
    2.0 * (1.0 - dist.cdf(t.abs()))
}

fn median(mut values: Vec<f64>) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        Some(values[mid])
    } else {
        Some((values[mid - 1] + values[mid]) / 2.0)
    }
}

fn sample_var(values: &[f64], mean: f64) -> f64 {
    values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64
}

fn percentile_sorted(values: &[f64], q: f64) -> f64 {
    let idx = ((values.len() - 1) as f64 * q).round() as usize;
    values[idx.min(values.len() - 1)]
}

fn sign_stability(point: f64, effects: &[f64]) -> f64 {
    if point > 0.0 {
        effects.iter().filter(|v| **v > 0.0).count() as f64 / effects.len() as f64
    } else if point < 0.0 {
        effects.iter().filter(|v| **v < 0.0).count() as f64 / effects.len() as f64
    } else {
        0.0
    }
}

fn shuffle(values: &mut [f64], rng: &mut Rng64) {
    for i in (1..values.len()).rev() {
        let j = rng.gen_range(i + 1);
        values.swap(i, j);
    }
}

fn test_name(test: Test) -> &'static str {
    match test {
        Test::Wilcoxon => "wilcoxon",
        Test::TTest => "t-test",
        Test::Permutation => "permutation",
    }
}

fn fmt(value: f64) -> String {
    if value.is_nan() {
        "NA".to_string()
    } else if value == 0.0 {
        "0".to_string()
    } else {
        format!("{value:.6}")
    }
}

struct Rng64 {
    state: u64,
}

impl Rng64 {
    fn new(seed: u64) -> Self {
        Self {
            state: seed ^ 0x9E37_79B9_7F4A_7C15,
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn gen_range(&mut self, upper: usize) -> usize {
        (self.next_u64() as usize) % upper
    }
}
