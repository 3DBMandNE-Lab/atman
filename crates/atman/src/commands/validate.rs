use anyhow::{Context, Result};
use atman_core::Platform;
use clap::Args as ClapArgs;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use super::parse_comparisons;
use crate::io::{atomic_write, read_measurements_long, read_samples};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Optional output path for a TSV report. Defaults to no file output.
    #[arg(long)]
    report: Option<PathBuf>,

    /// Optional comma-separated comparisons in A-B form for sample-count checks.
    #[arg(long)]
    groups: Option<String>,

    /// Minimum effective samples per group or paired subjects per comparison.
    #[arg(long, default_value_t = 2)]
    min_pairs: usize,

    /// Treat warnings as validation failures.
    #[arg(long, default_value_t = false)]
    strict: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Severity {
    Error,
    Warning,
}

impl Severity {
    fn as_str(self) -> &'static str {
        match self {
            Self::Error => "error",
            Self::Warning => "warning",
        }
    }
}

#[derive(Debug)]
struct Finding {
    severity: Severity,
    code: &'static str,
    message: String,
}

#[derive(Default)]
struct Findings {
    rows: Vec<Finding>,
}

impl Findings {
    fn error(&mut self, code: &'static str, message: impl Into<String>) {
        self.rows.push(Finding {
            severity: Severity::Error,
            code,
            message: message.into(),
        });
    }

    fn warning(&mut self, code: &'static str, message: impl Into<String>) {
        self.rows.push(Finding {
            severity: Severity::Warning,
            code,
            message: message.into(),
        });
    }

    fn error_count(&self) -> usize {
        self.rows
            .iter()
            .filter(|f| f.severity == Severity::Error)
            .count()
    }

    fn warning_count(&self) -> usize {
        self.rows
            .iter()
            .filter(|f| f.severity == Severity::Warning)
            .count()
    }
}

#[derive(Debug)]
struct RawTable {
    headers: Vec<String>,
    rows: Vec<HashMap<String, String>>,
}

pub fn run(args: Args) -> Result<()> {
    let mut findings = Findings::default();
    validate_dir(&args, &mut findings)?;

    for f in &findings.rows {
        eprintln!("{} [{}] {}", f.severity.as_str(), f.code, f.message);
    }
    eprintln!(
        "validate: {} error(s), {} warning(s)",
        findings.error_count(),
        findings.warning_count()
    );

    if let Some(path) = &args.report {
        write_report(path, &findings)?;
    }

    let fails = findings.error_count() > 0 || (args.strict && findings.warning_count() > 0);
    if fails {
        anyhow::bail!("validation failed");
    }
    Ok(())
}

fn validate_dir(args: &Args, findings: &mut Findings) -> Result<()> {
    let samples_path = args.input_dir.join("samples.tsv");
    let proteins_path = args.input_dir.join("proteins.tsv");
    let measurements_path = if args.input_dir.join("qc_measurements.tsv").exists() {
        args.input_dir.join("qc_measurements.tsv")
    } else {
        args.input_dir.join("measurements.tsv")
    };

    for path in [&samples_path, &proteins_path, &measurements_path] {
        if !path.exists() {
            findings.error(
                "missing_file",
                format!("required file {:?} does not exist", path),
            );
        }
    }
    if findings.error_count() > 0 {
        return Ok(());
    }

    let samples = read_raw_tsv(&samples_path)?;
    let proteins = read_raw_tsv(&proteins_path)?;
    let measurements = read_raw_tsv(&measurements_path)?;

    require_columns(
        findings,
        &samples_path,
        &samples,
        &[
            "sample_id",
            "subject_id",
            "condition",
            "is_control",
            "sample_type",
            "ingest_order",
        ],
    );
    require_columns(
        findings,
        &proteins_path,
        &proteins,
        &[
            "platform",
            "assay_id",
            "uniprot",
            "gene_symbol",
            "panel",
            "panel_lot",
        ],
    );
    require_columns(
        findings,
        &measurements_path,
        &measurements,
        &[
            "sample_id",
            "assay_id",
            "gene_symbol",
            "panel",
            "npx_source_str",
            "abundance",
            "abundance_raw",
            "abundance_unit",
            "qc_sample",
            "qc_assay",
            "detection_limit",
            "below_lod",
            "dropped_by_qc",
            "plate_id",
            "panel_lot",
            "ingest_order",
        ],
    );
    if findings.error_count() > 0 {
        return Ok(());
    }

    validate_samples(findings, &samples);
    validate_proteins(findings, &proteins);
    validate_measurements(findings, &measurements);
    validate_references(findings, &samples, &proteins, &measurements);
    validate_typed_readers(findings, &samples_path, &proteins_path, &measurements_path);
    validate_groups(args, findings, &samples_path, &measurements_path);

    Ok(())
}

fn read_raw_tsv(path: &Path) -> Result<RawTable> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers: Vec<String> = reader.headers()?.iter().map(|s| s.to_string()).collect();
    let mut rows = Vec::new();
    for result in reader.records() {
        let record = result.with_context(|| format!("reading record from {:?}", path))?;
        let mut row = HashMap::new();
        for (header, value) in headers.iter().zip(record.iter()) {
            row.insert(header.clone(), value.to_string());
        }
        rows.push(row);
    }
    Ok(RawTable { headers, rows })
}

fn require_columns(findings: &mut Findings, path: &Path, table: &RawTable, required: &[&str]) {
    let have: BTreeSet<&str> = table.headers.iter().map(|s| s.as_str()).collect();
    for col in required {
        if !have.contains(col) {
            findings.error(
                "missing_column",
                format!("{:?} is missing required column {}", path, col),
            );
        }
    }
}

fn validate_samples(findings: &mut Findings, samples: &RawTable) {
    let mut ids = BTreeSet::new();
    for (i, row) in samples.rows.iter().enumerate() {
        let n = i + 2;
        let sample_id = value(row, "sample_id");
        if sample_id.is_empty() {
            findings.error(
                "empty_sample_id",
                format!("samples.tsv line {n}: empty sample_id"),
            );
        } else if !ids.insert(sample_id.to_string()) {
            findings.error(
                "duplicate_sample_id",
                format!("samples.tsv line {n}: duplicate sample_id {sample_id:?}"),
            );
        }
        if !matches!(value(row, "is_control"), "0" | "1") {
            findings.error(
                "invalid_is_control",
                format!("samples.tsv line {n}: is_control must be 0 or 1"),
            );
        }
        if value(row, "is_control") == "0" && value(row, "condition").is_empty() {
            findings.error(
                "missing_condition",
                format!("samples.tsv line {n}: biological sample has empty condition"),
            );
        }
    }
}

fn validate_proteins(findings: &mut Findings, proteins: &RawTable) {
    let mut keys = BTreeSet::new();
    for (i, row) in proteins.rows.iter().enumerate() {
        let n = i + 2;
        let platform = value(row, "platform");
        let assay = value(row, "assay_id");
        if platform.parse::<Platform>().is_err() {
            findings.error(
                "invalid_platform",
                format!("proteins.tsv line {n}: unsupported platform {platform:?}"),
            );
        }
        if assay.is_empty() {
            findings.error(
                "empty_assay_id",
                format!("proteins.tsv line {n}: empty assay_id"),
            );
        }
        if !platform.is_empty() && !assay.is_empty() && !keys.insert((platform, assay)) {
            findings.error(
                "duplicate_protein_key",
                format!("proteins.tsv line {n}: duplicate platform/assay_id key"),
            );
        }
    }
}

fn validate_measurements(findings: &mut Findings, measurements: &RawTable) {
    let mut keys = BTreeSet::new();
    for (i, row) in measurements.rows.iter().enumerate() {
        let n = i + 2;
        let platform = value(row, "platform");
        if !platform.is_empty() && platform.parse::<Platform>().is_err() {
            findings.error(
                "invalid_platform",
                format!("measurements line {n}: unsupported platform {platform:?}"),
            );
        }
        let assay = value(row, "assay_id");
        let sample = value(row, "sample_id");
        if sample.is_empty() || assay.is_empty() {
            findings.error(
                "empty_measurement_key",
                format!("measurements line {n}: sample_id and assay_id are required"),
            );
        }
        let key_platform = if platform.is_empty() {
            "olink_explore_ngs"
        } else {
            platform
        };
        if !sample.is_empty() && !assay.is_empty() && !keys.insert((key_platform, assay, sample)) {
            findings.error(
                "duplicate_measurement_key",
                format!("measurements line {n}: duplicate platform/assay_id/sample_id key"),
            );
        }
        for col in ["qc_sample", "qc_assay"] {
            if !matches!(value(row, col), "PASS" | "WARN" | "FAIL" | "") {
                findings.error(
                    "invalid_qc_flag",
                    format!("measurements line {n}: {col} must be PASS, WARN, or FAIL"),
                );
            }
        }
        for col in ["below_lod", "dropped_by_qc"] {
            if !matches!(value(row, col), "0" | "1" | "") {
                findings.error(
                    "invalid_boolean",
                    format!("measurements line {n}: {col} must be 0 or 1"),
                );
            }
        }
        let unit = value(row, "abundance_unit");
        if !is_known_abundance_unit(unit) {
            findings.warning(
                "unknown_abundance_unit",
                format!("measurements line {n}: abundance_unit {unit:?} will be read as raw"),
            );
        }
        for col in ["abundance", "abundance_raw", "detection_limit"] {
            let v = value(row, col);
            if !v.is_empty() && v.parse::<f64>().is_err() {
                findings.error(
                    "invalid_number",
                    format!("measurements line {n}: {col} is not numeric"),
                );
            }
        }
    }
}

fn validate_references(
    findings: &mut Findings,
    samples: &RawTable,
    proteins: &RawTable,
    measurements: &RawTable,
) {
    let sample_ids: BTreeSet<&str> = samples.rows.iter().map(|r| value(r, "sample_id")).collect();
    let protein_keys: BTreeSet<(&str, &str)> = proteins
        .rows
        .iter()
        .map(|r| (value(r, "platform"), value(r, "assay_id")))
        .collect();
    for (i, row) in measurements.rows.iter().enumerate() {
        let n = i + 2;
        let sample = value(row, "sample_id");
        if !sample_ids.contains(sample) {
            findings.error(
                "orphan_measurement_sample",
                format!("measurements line {n}: sample_id {sample:?} not found in samples.tsv"),
            );
        }
        let platform = if value(row, "platform").is_empty() {
            "olink_explore_ngs"
        } else {
            value(row, "platform")
        };
        let assay = value(row, "assay_id");
        if !protein_keys.contains(&(platform, assay)) {
            findings.error(
                "orphan_measurement_assay",
                format!(
                    "measurements line {n}: platform/assay_id ({platform:?}, {assay:?}) not found in proteins.tsv"
                ),
            );
        }
    }
}

fn validate_typed_readers(
    findings: &mut Findings,
    samples_path: &Path,
    proteins_path: &Path,
    measurements_path: &Path,
) {
    if let Err(e) = read_samples(samples_path) {
        findings.error(
            "sample_reader_failed",
            format!("read_samples failed: {e:#}"),
        );
    }
    if let Err(e) = crate::io::read_proteins(proteins_path) {
        findings.error(
            "protein_reader_failed",
            format!("read_proteins failed: {e:#}"),
        );
    }
    if let Err(e) = read_measurements_long(measurements_path) {
        findings.error(
            "measurement_reader_failed",
            format!("read_measurements_long failed: {e:#}"),
        );
    }
}

fn validate_groups(
    args: &Args,
    findings: &mut Findings,
    samples_path: &Path,
    measurements_path: &Path,
) {
    let Some(groups) = &args.groups else {
        return;
    };
    let comparisons = match parse_comparisons(groups) {
        Ok(v) => v,
        Err(e) => {
            findings.error("invalid_groups", e.to_string());
            return;
        }
    };
    let samples = match read_samples(samples_path) {
        Ok(v) => v,
        Err(_) => return,
    };
    let measurements = match read_measurements_long(measurements_path) {
        Ok(v) => v,
        Err(_) => return,
    };
    let sample_by_id: HashMap<&str, &atman_core::Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();
    let mut by_condition: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for m in &measurements {
        if m.effective_abundance().is_none() {
            continue;
        }
        let Some(sample) = sample_by_id.get(m.sample_id.as_str()) else {
            continue;
        };
        if sample.is_control {
            continue;
        }
        let Some(condition) = &sample.condition else {
            continue;
        };
        let id = sample.subject_id.as_ref().unwrap_or(&sample.sample_id);
        by_condition
            .entry(condition.clone())
            .or_default()
            .insert(id.clone());
    }
    for (a, b) in comparisons {
        let na = by_condition.get(&a).map(BTreeSet::len).unwrap_or(0);
        let nb = by_condition.get(&b).map(BTreeSet::len).unwrap_or(0);
        if na < args.min_pairs || nb < args.min_pairs {
            findings.warning(
                "insufficient_group_samples",
                format!(
                    "comparison {a}-{b}: effective sample counts are {na} and {nb}, below min-pairs {}",
                    args.min_pairs
                ),
            );
        }
    }
}

fn is_known_abundance_unit(unit: &str) -> bool {
    matches!(
        unit,
        "log2_npx"
            | "npx"
            | "log2_intensity"
            | "log2_lfq"
            | "log2_rfu"
            | "log2_pg_quantity"
            | "log2_diann_pg_quantity"
            | "ibaq"
            | "raw"
    )
}

fn value<'a>(row: &'a HashMap<String, String>, key: &str) -> &'a str {
    row.get(key).map(|s| s.as_str()).unwrap_or("")
}

fn write_report(path: &Path, findings: &Findings) -> Result<()> {
    let mut buf = String::from("severity\tcode\tmessage\n");
    for f in &findings.rows {
        buf.push_str(f.severity.as_str());
        buf.push('\t');
        buf.push_str(f.code);
        buf.push('\t');
        buf.push_str(&f.message.replace(['\t', '\n'], " "));
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}
