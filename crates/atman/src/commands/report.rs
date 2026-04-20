use anyhow::{Context, Result};
use atman_core::{MeasurementRecord, ProteinIdentity, Sample};
use clap::{Args as ClapArgs, Subcommand};
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, hash_canonical_inputs, read_measurements_long, read_proteins, read_samples,
    sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Summarize QC, missingness, and effective sample support.
    Qc(QcArgs),
}

#[derive(ClapArgs, Debug)]
pub struct QcArgs {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Directory where report TSV files will be written.
    #[arg(long)]
    output_dir: PathBuf,

    /// Warn when a condition has fewer effective subjects than this threshold.
    #[arg(long, default_value_t = 2)]
    min_subjects: usize,

    /// Warn when a protein has an effective measurement fraction below this threshold.
    #[arg(long, default_value_t = 0.5)]
    sparse_threshold: f64,
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Qc(args) => run_qc(args),
    }
}

fn run_qc(args: QcArgs) -> Result<()> {
    let started_at = SystemTime::now();
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let measurements_path = if args.input_dir.join("qc_measurements.tsv").exists() {
        args.input_dir.join("qc_measurements.tsv")
    } else {
        args.input_dir.join("measurements.tsv")
    };
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let proteins = read_proteins(&args.input_dir.join("proteins.tsv"))?;
    let measurements = read_measurements_long(&measurements_path)?;

    let report = QcReport::build(&samples, &proteins, &measurements);
    report.write(&args.output_dir)?;

    let warnings = report.emit_warnings(args.min_subjects, args.sparse_threshold);
    eprintln!(
        "report qc: samples={} proteins={} measurements={} warnings={}",
        samples.len(),
        proteins.len(),
        measurements.len(),
        warnings
    );

    let finished_at = SystemTime::now();
    let qc_summary_path = args.output_dir.join("qc_summary.tsv");
    let sample_qc_path = args.output_dir.join("sample_qc.tsv");
    let protein_qc_path = args.output_dir.join("protein_qc.tsv");
    let condition_counts_path = args.output_dir.join("condition_counts.tsv");
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &[
            "qc_measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
    )?;
    let sidecar = sidecar_path_for(&qc_summary_path);
    write_run_sidecar(
        &sidecar,
        "report qc",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
            "min-subjects": args.min_subjects,
            "sparse-threshold": args.sparse_threshold,
        }),
        &inputs_sha256,
        &[
            qc_summary_path.clone(),
            sample_qc_path.clone(),
            protein_qc_path.clone(),
            condition_counts_path.clone(),
        ],
        started_at,
        finished_at,
        None,
)?;
    eprintln!("report qc: sidecar={}", sidecar.display());
    Ok(())
}

#[derive(Default)]
struct QcReport {
    summary: Vec<(String, String)>,
    samples: Vec<SampleQcRow>,
    proteins: Vec<ProteinQcRow>,
    conditions: Vec<ConditionCountRow>,
}

#[derive(Default)]
struct Count {
    total: usize,
    effective: usize,
    qc_masked: usize,
    below_lod: usize,
}

#[derive(Default)]
struct SampleQcRow {
    sample_id: String,
    subject_id: String,
    condition: String,
    is_control: bool,
    count: Count,
}

#[derive(Default)]
struct ProteinQcRow {
    platform: String,
    assay_id: String,
    gene_symbol: String,
    panel: String,
    count: Count,
    effective_conditions: BTreeSet<String>,
}

#[derive(Default)]
struct ConditionScratch {
    samples: BTreeSet<String>,
    subjects: BTreeSet<String>,
    effective_samples: BTreeSet<String>,
    effective_subjects: BTreeSet<String>,
    effective_measurements: usize,
}

struct ConditionCountRow {
    condition: String,
    n_samples: usize,
    n_subjects: usize,
    n_effective_samples: usize,
    n_effective_subjects: usize,
    n_effective_measurements: usize,
}

impl QcReport {
    fn build(
        samples: &[Sample],
        proteins: &[ProteinIdentity],
        measurements: &[MeasurementRecord],
    ) -> Self {
        let sample_by_id: HashMap<&str, &Sample> =
            samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

        let mut sample_rows: BTreeMap<String, SampleQcRow> = samples
            .iter()
            .map(|s| {
                (
                    s.sample_id.clone(),
                    SampleQcRow {
                        sample_id: s.sample_id.clone(),
                        subject_id: s.subject_id.clone().unwrap_or_default(),
                        condition: s.condition.clone().unwrap_or_default(),
                        is_control: s.is_control,
                        count: Count::default(),
                    },
                )
            })
            .collect();

        let protein_key = |platform: &str, assay_id: &str| format!("{platform}\t{assay_id}");
        let mut protein_rows: BTreeMap<String, ProteinQcRow> = proteins
            .iter()
            .map(|p| {
                let platform = p.platform.as_str().to_string();
                let assay_id = p.assay_id.0.clone();
                (
                    protein_key(&platform, &assay_id),
                    ProteinQcRow {
                        platform,
                        assay_id,
                        gene_symbol: p.gene_symbol.clone().unwrap_or_default(),
                        panel: p.panel.clone().unwrap_or_default(),
                        count: Count::default(),
                        effective_conditions: BTreeSet::new(),
                    },
                )
            })
            .collect();

        let mut conditions: BTreeMap<String, ConditionScratch> = BTreeMap::new();
        for s in samples.iter().filter(|s| !s.is_control) {
            if let Some(condition) = &s.condition {
                let id = s.subject_id.as_ref().unwrap_or(&s.sample_id);
                let row = conditions.entry(condition.clone()).or_default();
                row.samples.insert(s.sample_id.clone());
                row.subjects.insert(id.clone());
            }
        }

        for m in measurements {
            let effective = m.effective_abundance().is_some();

            if let Some(row) = sample_rows.get_mut(&m.sample_id) {
                add_count(&mut row.count, m, effective);
            }

            let platform = m.platform.as_str();
            let key = protein_key(platform, &m.assay_id.0);
            let protein_row = protein_rows.entry(key).or_insert_with(|| ProteinQcRow {
                platform: platform.to_string(),
                assay_id: m.assay_id.0.clone(),
                gene_symbol: m.gene_symbol.clone().unwrap_or_default(),
                panel: m.panel.clone().unwrap_or_default(),
                count: Count::default(),
                effective_conditions: BTreeSet::new(),
            });
            add_count(&mut protein_row.count, m, effective);

            let Some(sample) = sample_by_id.get(m.sample_id.as_str()) else {
                continue;
            };
            if sample.is_control {
                continue;
            }
            let Some(condition) = &sample.condition else {
                continue;
            };
            if effective {
                protein_row.effective_conditions.insert(condition.clone());
                let id = sample.subject_id.as_ref().unwrap_or(&sample.sample_id);
                let row = conditions.entry(condition.clone()).or_default();
                row.effective_samples.insert(sample.sample_id.clone());
                row.effective_subjects.insert(id.clone());
                row.effective_measurements += 1;
            }
        }

        let total_measurements = measurements.len();
        let effective_measurements = measurements
            .iter()
            .filter(|m| m.effective_abundance().is_some())
            .count();
        let qc_masked = measurements.iter().filter(|m| m.dropped_by_qc).count();
        let below_lod = measurements.iter().filter(|m| m.below_lod).count();
        let control_samples = samples.iter().filter(|s| s.is_control).count();
        let biological_samples = samples.len().saturating_sub(control_samples);
        let condition_count = conditions.len();

        let summary = vec![
            ("sample_count".to_string(), samples.len().to_string()),
            (
                "biological_sample_count".to_string(),
                biological_samples.to_string(),
            ),
            (
                "control_sample_count".to_string(),
                control_samples.to_string(),
            ),
            ("protein_count".to_string(), proteins.len().to_string()),
            ("condition_count".to_string(), condition_count.to_string()),
            (
                "measurement_count".to_string(),
                total_measurements.to_string(),
            ),
            (
                "effective_measurement_count".to_string(),
                effective_measurements.to_string(),
            ),
            (
                "missing_measurement_count".to_string(),
                total_measurements
                    .saturating_sub(effective_measurements)
                    .to_string(),
            ),
            (
                "missing_fraction".to_string(),
                fraction(
                    total_measurements.saturating_sub(effective_measurements),
                    total_measurements,
                ),
            ),
            ("qc_masked_count".to_string(), qc_masked.to_string()),
            (
                "qc_masked_fraction".to_string(),
                fraction(qc_masked, total_measurements),
            ),
            ("below_lod_count".to_string(), below_lod.to_string()),
            (
                "below_lod_fraction".to_string(),
                fraction(below_lod, total_measurements),
            ),
        ];

        let condition_rows = conditions
            .into_iter()
            .map(|(condition, c)| ConditionCountRow {
                condition,
                n_samples: c.samples.len(),
                n_subjects: c.subjects.len(),
                n_effective_samples: c.effective_samples.len(),
                n_effective_subjects: c.effective_subjects.len(),
                n_effective_measurements: c.effective_measurements,
            })
            .collect();

        Self {
            summary,
            samples: sample_rows.into_values().collect(),
            proteins: protein_rows.into_values().collect(),
            conditions: condition_rows,
        }
    }

    fn write(&self, output_dir: &Path) -> Result<()> {
        write_qc_summary(&output_dir.join("qc_summary.tsv"), &self.summary)?;
        write_sample_qc(&output_dir.join("sample_qc.tsv"), &self.samples)?;
        write_protein_qc(&output_dir.join("protein_qc.tsv"), &self.proteins)?;
        write_condition_counts(&output_dir.join("condition_counts.tsv"), &self.conditions)?;
        Ok(())
    }

    fn emit_warnings(&self, min_subjects: usize, sparse_threshold: f64) -> usize {
        let mut warnings = 0;
        for c in &self.conditions {
            if c.n_effective_subjects < min_subjects {
                warnings += 1;
                eprintln!(
                    "warning [small_condition] condition {:?} has {} effective subject(s), below {}",
                    c.condition, c.n_effective_subjects, min_subjects
                );
            }
        }
        for p in &self.proteins {
            if p.count.total == 0 {
                warnings += 1;
                eprintln!(
                    "warning [unmeasured_protein] protein {:?}/{:?} has no measurements",
                    p.platform, p.assay_id
                );
                continue;
            }
            let effective_fraction = p.count.effective as f64 / p.count.total as f64;
            if effective_fraction < sparse_threshold {
                warnings += 1;
                eprintln!(
                    "warning [sparse_protein] protein {:?}/{:?} effective fraction {} is below {}",
                    p.platform,
                    p.assay_id,
                    format_float(effective_fraction),
                    format_float(sparse_threshold)
                );
            }
        }
        warnings
    }
}

fn add_count(count: &mut Count, m: &MeasurementRecord, effective: bool) {
    count.total += 1;
    if effective {
        count.effective += 1;
    }
    if m.dropped_by_qc {
        count.qc_masked += 1;
    }
    if m.below_lod {
        count.below_lod += 1;
    }
}

fn write_qc_summary(path: &Path, rows: &[(String, String)]) -> Result<()> {
    let mut buf = String::from("metric\tvalue\n");
    for (metric, value) in rows {
        buf.push_str(metric);
        buf.push('\t');
        buf.push_str(value);
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn write_sample_qc(path: &Path, rows: &[SampleQcRow]) -> Result<()> {
    let mut buf = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tn_measurements\tn_effective\tn_missing\tmissing_fraction\tn_qc_masked\tn_below_lod\n",
    );
    for r in rows {
        push_qc_prefix(&mut buf, &[&r.sample_id, &r.subject_id, &r.condition]);
        buf.push_str(&(r.is_control as u8).to_string());
        buf.push('\t');
        push_count(&mut buf, &r.count);
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn write_protein_qc(path: &Path, rows: &[ProteinQcRow]) -> Result<()> {
    let mut buf = String::from(
        "platform\tassay_id\tgene_symbol\tpanel\tn_measurements\tn_effective\tn_missing\tmissing_fraction\tn_qc_masked\tn_below_lod\tn_effective_conditions\n",
    );
    for r in rows {
        push_qc_prefix(
            &mut buf,
            &[&r.platform, &r.assay_id, &r.gene_symbol, &r.panel],
        );
        push_count(&mut buf, &r.count);
        buf.push('\t');
        buf.push_str(&r.effective_conditions.len().to_string());
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn write_condition_counts(path: &Path, rows: &[ConditionCountRow]) -> Result<()> {
    let mut buf = String::from(
        "condition\tn_samples\tn_subjects\tn_effective_samples\tn_effective_subjects\tn_effective_measurements\n",
    );
    for r in rows {
        buf.push_str(&r.condition);
        buf.push('\t');
        buf.push_str(&r.n_samples.to_string());
        buf.push('\t');
        buf.push_str(&r.n_subjects.to_string());
        buf.push('\t');
        buf.push_str(&r.n_effective_samples.to_string());
        buf.push('\t');
        buf.push_str(&r.n_effective_subjects.to_string());
        buf.push('\t');
        buf.push_str(&r.n_effective_measurements.to_string());
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn push_qc_prefix(buf: &mut String, fields: &[&str]) {
    for field in fields {
        buf.push_str(field);
        buf.push('\t');
    }
}

fn push_count(buf: &mut String, count: &Count) {
    buf.push_str(&count.total.to_string());
    buf.push('\t');
    buf.push_str(&count.effective.to_string());
    buf.push('\t');
    buf.push_str(&count.total.saturating_sub(count.effective).to_string());
    buf.push('\t');
    buf.push_str(&fraction(
        count.total.saturating_sub(count.effective),
        count.total,
    ));
    buf.push('\t');
    buf.push_str(&count.qc_masked.to_string());
    buf.push('\t');
    buf.push_str(&count.below_lod.to_string());
}

fn fraction(num: usize, denom: usize) -> String {
    if denom == 0 {
        return String::new();
    }
    format_float(num as f64 / denom as f64)
}

fn format_float(v: f64) -> String {
    format!("{v}")
}
