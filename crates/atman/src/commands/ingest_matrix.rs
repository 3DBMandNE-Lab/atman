use anyhow::{bail, Context, Result};
use atman_core::Platform;
use clap::{Args as ClapArgs, ValueEnum};
use serde_json::json;
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{atomic_write, hash_labeled_inputs, sidecar_path_for, write_run_sidecar};

const MEASUREMENT_HEADER: &[&str] = &[
    "platform",
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
];

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Orientation {
    /// Matrix rows are proteins/features and matrix columns are samples.
    ProteinsRows,
    /// Matrix rows are samples and matrix columns are proteins/features.
    SamplesRows,
}

#[derive(Clone, Copy, Debug, ValueEnum, PartialEq, Eq)]
enum Normalize {
    /// Pass abundances through unchanged. Appropriate for data that is
    /// already cross-sample normalized (e.g. Olink NPX).
    None,
    /// Per-sample median centering on log-scale: subtract each sample's
    /// median, add the mean of sample medians back. Preserves relative
    /// protein differences within samples.
    Median,
    /// Force identical empirical distributions across samples. Complete-case
    /// assays (observed in every sample) are quantile-normalized; remaining
    /// assays fall back to median centering so no data is silently dropped.
    Quantile,
}

impl Normalize {
    fn as_str(&self) -> &'static str {
        match self {
            Normalize::None => "none",
            Normalize::Median => "median",
            Normalize::Quantile => "quantile",
        }
    }
}

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Wide abundance matrix, CSV or TSV.
    #[arg(long)]
    matrix: PathBuf,

    /// Sample metadata table, CSV or TSV.
    #[arg(long)]
    samples: PathBuf,

    /// Protein metadata table. Required for --orientation samples-rows.
    #[arg(long)]
    proteins: Option<PathBuf>,

    /// Output directory for canonical Atman TSV files.
    #[arg(long)]
    output_dir: PathBuf,

    /// Layout of the matrix.
    #[arg(long, value_enum, default_value_t = Orientation::ProteinsRows)]
    orientation: Orientation,

    /// Proteomics platform identifier.
    #[arg(long)]
    platform: String,

    /// Unit tag for abundance values.
    #[arg(long)]
    abundance_unit: String,

    /// Matrix/protein metadata column holding assay IDs.
    #[arg(long)]
    assay_id_col: String,

    /// Matrix/protein metadata column holding gene symbols.
    #[arg(long)]
    gene_col: Option<String>,

    /// Matrix/protein metadata column holding UniProt IDs.
    #[arg(long)]
    uniprot_col: Option<String>,

    /// Constant panel name when --panel-col is not supplied.
    #[arg(long, default_value = "")]
    panel: String,

    /// Matrix/protein metadata column holding panel names.
    #[arg(long)]
    panel_col: Option<String>,

    /// Sample metadata column holding sample IDs.
    #[arg(long, default_value = "sample_id")]
    sample_id_col: String,

    /// Sample metadata column holding subject IDs.
    #[arg(long)]
    subject_id_col: Option<String>,

    /// Sample metadata column holding conditions.
    #[arg(long)]
    condition_col: String,

    /// Sample metadata column holding sample types.
    #[arg(long)]
    sample_type_col: Option<String>,

    /// Sample metadata column holding 0/1 control flags.
    #[arg(long)]
    is_control_col: Option<String>,

    /// Log2-transform positive linear abundances. Non-positive values are QC-masked.
    #[arg(long, default_value_t = false)]
    log2_transform: bool,

    /// Cross-sample normalization applied after any `--log2-transform`.
    /// Default `none` passes values through. `median` subtracts each sample's
    /// median and adds the grand mean of sample medians. `quantile` forces
    /// identical empirical distributions across samples for complete-case
    /// assays and falls back to median centering for partially-observed
    /// assays. Olink NPX is already normalized — use `none`. For SomaScan,
    /// MaxQuant/LFQ, DIA-NN, and Spectronaut, either normalize upstream or
    /// pass `--normalize median` (recommended default for log-scale MS).
    #[arg(long, value_enum, default_value_t = Normalize::None)]
    normalize: Normalize,

    /// Disable the cross-sample-median-spread warning that fires when
    /// `--normalize none` appears to be applied to un-normalized data.
    #[arg(long, default_value_t = false)]
    skip_normalization_check: bool,
}

#[derive(Debug)]
struct Table {
    headers: Vec<String>,
    index: HashMap<String, usize>,
    rows: Vec<Vec<String>>,
}

#[derive(Debug)]
struct ProteinRow {
    assay_id: String,
    uniprot: String,
    gene_symbol: String,
    panel: String,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    let platform: Platform = args.platform.parse().map_err(anyhow::Error::msg)?;
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let matrix = read_table(&args.matrix)?;
    let samples = read_table(&args.samples)?;
    let sample_ids = write_samples(&args, &samples)?;

    let (proteins, mut measurements) = match args.orientation {
        Orientation::ProteinsRows => {
            let proteins = protein_rows_from_table(&args, &matrix)?;
            let measurements =
                measurements_from_proteins_rows(&args, &platform, &matrix, &sample_ids)?;
            (proteins, measurements)
        }
        Orientation::SamplesRows => {
            let protein_path = args
                .proteins
                .as_ref()
                .context("--orientation samples-rows requires --proteins")?;
            let protein_meta = read_table(protein_path)?;
            let proteins = protein_rows_from_table(&args, &protein_meta)?;
            let measurements = measurements_from_samples_rows(&args, &platform, &matrix, &proteins)?;
            (proteins, measurements)
        }
    };

    if !args.skip_normalization_check && args.normalize == Normalize::None {
        warn_if_unnormalized(&measurements);
    }
    apply_normalization(&mut measurements, args.normalize)?;

    write_proteins(&args.output_dir.join("proteins.tsv"), &platform, &proteins)?;
    write_measurements(
        &args.output_dir.join("measurements.tsv"),
        &measurements,
        &args.abundance_unit,
    )?;

    eprintln!(
        "ingest-matrix: {} samples, {} proteins, {} measurements",
        sample_ids.len(),
        proteins.len(),
        measurements.len()
    );

    let finished_at = SystemTime::now();
    let mut input_entries: Vec<(&str, &Path)> = vec![
        ("matrix", args.matrix.as_path()),
        ("samples", args.samples.as_path()),
    ];
    if let Some(p) = args.proteins.as_ref() {
        input_entries.push(("proteins", p.as_path()));
    }
    let inputs_sha256 = hash_labeled_inputs(&input_entries)?;
    let measurements_out = args.output_dir.join("measurements.tsv");
    let proteins_out = args.output_dir.join("proteins.tsv");
    let samples_out = args.output_dir.join("samples.tsv");
    let outputs: Vec<PathBuf> = vec![
        measurements_out.clone(),
        proteins_out.clone(),
        samples_out.clone(),
    ];
    let orientation = match args.orientation {
        Orientation::ProteinsRows => "proteins-rows",
        Orientation::SamplesRows => "samples-rows",
    };
    let sidecar = sidecar_path_for(&measurements_out);
    write_run_sidecar(
        &sidecar,
        "ingest-matrix",
        json!({
            "matrix": args.matrix.display().to_string(),
            "samples": args.samples.display().to_string(),
            "proteins": args.proteins.as_ref().map(|p| p.display().to_string()),
            "output-dir": args.output_dir.display().to_string(),
            "orientation": orientation,
            "platform": args.platform,
            "abundance-unit": args.abundance_unit,
            "assay-id-col": args.assay_id_col,
            "gene-col": args.gene_col,
            "uniprot-col": args.uniprot_col,
            "panel": args.panel,
            "panel-col": args.panel_col,
            "sample-id-col": args.sample_id_col,
            "subject-id-col": args.subject_id_col,
            "condition-col": args.condition_col,
            "sample-type-col": args.sample_type_col,
            "is-control-col": args.is_control_col,
            "log2-transform": args.log2_transform,
            "normalize": args.normalize.as_str(),
            "skip-normalization-check": args.skip_normalization_check,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("ingest-matrix: sidecar={}", sidecar.display());
    Ok(())
}

fn read_table(path: &Path) -> Result<Table> {
    let delimiter = match path.extension().and_then(|s| s.to_str()) {
        Some("tsv" | "txt") => b'\t',
        _ => b',',
    };
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers: Vec<String> = reader.headers()?.iter().map(str::to_string).collect();
    let index = headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.clone(), i))
        .collect();
    let mut rows = Vec::new();
    for row in reader.records() {
        rows.push(row?.iter().map(str::to_string).collect());
    }
    Ok(Table {
        headers,
        index,
        rows,
    })
}

fn write_samples(args: &Args, table: &Table) -> Result<Vec<String>> {
    let sample_idx = table.need(&args.sample_id_col)?;
    let subject_idx = args
        .subject_id_col
        .as_ref()
        .map(|c| table.need(c))
        .transpose()?;
    let condition_idx = table.need(&args.condition_col)?;
    let sample_type_idx = args
        .sample_type_col
        .as_ref()
        .map(|c| table.need(c))
        .transpose()?;
    let is_control_idx = args
        .is_control_col
        .as_ref()
        .map(|c| table.need(c))
        .transpose()?;

    let consumed: BTreeSet<&str> = [
        Some(args.sample_id_col.as_str()),
        args.subject_id_col.as_deref(),
        Some(args.condition_col.as_str()),
        args.sample_type_col.as_deref(),
        args.is_control_col.as_deref(),
    ]
    .into_iter()
    .flatten()
    .collect();
    let extra_cols: Vec<&str> = table
        .headers
        .iter()
        .map(String::as_str)
        .filter(|h| !consumed.contains(h))
        .collect();

    let mut seen = BTreeSet::new();
    let mut sample_ids = Vec::new();
    let mut buf =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order");
    for col in &extra_cols {
        buf.push('\t');
        buf.push_str(col);
    }
    buf.push('\n');

    for (i, row) in table.rows.iter().enumerate() {
        let sample_id = cell(row, sample_idx).to_string();
        if sample_id.is_empty() {
            bail!("samples row {} has empty sample ID", i + 2);
        }
        if !seen.insert(sample_id.clone()) {
            bail!("duplicate sample ID {:?}", sample_id);
        }
        sample_ids.push(sample_id.clone());
        let subject_id = subject_idx
            .map(|idx| cell(row, idx).to_string())
            .unwrap_or_else(|| sample_id.clone());
        let condition = cell(row, condition_idx);
        let sample_type = sample_type_idx.map(|idx| cell(row, idx)).unwrap_or("");
        let is_control = is_control_idx
            .map(|idx| parse_is_control(cell(row, idx)))
            .unwrap_or(false);

        buf.push_str(&sample_id);
        buf.push('\t');
        buf.push_str(&subject_id);
        buf.push('\t');
        buf.push_str(condition);
        buf.push('\t');
        buf.push_str(&(is_control as u8).to_string());
        buf.push('\t');
        buf.push_str(sample_type);
        buf.push('\t');
        buf.push_str(&(i + 1).to_string());
        for col in &extra_cols {
            buf.push('\t');
            buf.push_str(cell(row, table.need(col)?));
        }
        buf.push('\n');
    }

    atomic_write(&args.output_dir.join("samples.tsv"), buf.as_bytes())?;
    Ok(sample_ids)
}

fn protein_rows_from_table(args: &Args, table: &Table) -> Result<Vec<ProteinRow>> {
    let assay_idx = table.need(&args.assay_id_col)?;
    let gene_idx = args.gene_col.as_ref().map(|c| table.need(c)).transpose()?;
    let uniprot_idx = args
        .uniprot_col
        .as_ref()
        .map(|c| table.need(c))
        .transpose()?;
    let panel_idx = args.panel_col.as_ref().map(|c| table.need(c)).transpose()?;

    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for (i, row) in table.rows.iter().enumerate() {
        let assay_id = cell(row, assay_idx).to_string();
        if assay_id.is_empty() {
            bail!("protein row {} has empty assay ID", i + 2);
        }
        if !seen.insert(assay_id.clone()) {
            bail!("duplicate assay ID {:?}", assay_id);
        }
        out.push(ProteinRow {
            assay_id,
            uniprot: uniprot_idx
                .map(|idx| cell(row, idx).to_string())
                .unwrap_or_default(),
            gene_symbol: gene_idx
                .map(|idx| cell(row, idx).to_string())
                .unwrap_or_default(),
            panel: panel_idx
                .map(|idx| cell(row, idx).to_string())
                .unwrap_or_else(|| args.panel.clone()),
        });
    }
    Ok(out)
}

fn measurements_from_proteins_rows(
    args: &Args,
    platform: &Platform,
    matrix: &Table,
    sample_ids: &[String],
) -> Result<Vec<MeasurementRow>> {
    let mut meta_cols = BTreeSet::new();
    meta_cols.insert(args.assay_id_col.as_str());
    if let Some(col) = &args.gene_col {
        meta_cols.insert(col.as_str());
    }
    if let Some(col) = &args.uniprot_col {
        meta_cols.insert(col.as_str());
    }
    if let Some(col) = &args.panel_col {
        meta_cols.insert(col.as_str());
    }
    let sample_set: BTreeSet<&str> = sample_ids.iter().map(String::as_str).collect();
    let value_cols: Vec<(usize, String)> = matrix
        .headers
        .iter()
        .enumerate()
        .filter(|(_, h)| !meta_cols.contains(h.as_str()) && sample_set.contains(h.as_str()))
        .map(|(i, h)| (i, h.clone()))
        .collect();
    if value_cols.is_empty() {
        bail!("no matrix columns matched sample IDs from samples table");
    }

    let proteins = protein_rows_from_table(args, matrix)?;
    let mut out = Vec::new();
    let mut order = 1_u64;
    for (protein, row) in proteins.iter().zip(&matrix.rows) {
        for (idx, sample_id) in &value_cols {
            out.push(measurement_row(
                platform,
                sample_id,
                protein,
                cell(row, *idx),
                args.log2_transform,
                order,
            )?);
            order += 1;
        }
    }
    Ok(out)
}

fn measurements_from_samples_rows(
    args: &Args,
    platform: &Platform,
    matrix: &Table,
    proteins: &[ProteinRow],
) -> Result<Vec<MeasurementRow>> {
    let sample_idx = matrix.need(&args.sample_id_col)?;
    let protein_by_assay: HashMap<&str, &ProteinRow> =
        proteins.iter().map(|p| (p.assay_id.as_str(), p)).collect();
    let value_cols: Vec<(usize, &ProteinRow)> = matrix
        .headers
        .iter()
        .enumerate()
        .filter_map(|(i, h)| protein_by_assay.get(h.as_str()).map(|p| (i, *p)))
        .collect();
    if value_cols.is_empty() {
        bail!("no matrix columns matched assay IDs from protein metadata");
    }

    let mut out = Vec::new();
    let mut order = 1_u64;
    for row in &matrix.rows {
        let sample_id = cell(row, sample_idx);
        for (idx, protein) in &value_cols {
            out.push(measurement_row(
                platform,
                sample_id,
                protein,
                cell(row, *idx),
                args.log2_transform,
                order,
            )?);
            order += 1;
        }
    }
    Ok(out)
}

#[derive(Debug)]
struct MeasurementRow {
    platform: Platform,
    sample_id: String,
    assay_id: String,
    gene_symbol: String,
    panel: String,
    source: String,
    abundance: Option<f64>,
    dropped_by_qc: bool,
    ingest_order: u64,
}

fn measurement_row(
    platform: &Platform,
    sample_id: &str,
    protein: &ProteinRow,
    source: &str,
    log2_transform: bool,
    ingest_order: u64,
) -> Result<MeasurementRow> {
    let abundance = parse_abundance(source, log2_transform)?;
    Ok(MeasurementRow {
        platform: platform.clone(),
        sample_id: sample_id.to_string(),
        assay_id: protein.assay_id.clone(),
        gene_symbol: protein.gene_symbol.clone(),
        panel: protein.panel.clone(),
        source: source.to_string(),
        dropped_by_qc: abundance.is_none(),
        abundance,
        ingest_order,
    })
}

fn parse_abundance(raw: &str, log2_transform: bool) -> Result<Option<f64>> {
    let raw = raw.trim();
    if raw.is_empty() || raw.eq_ignore_ascii_case("nan") || raw.eq_ignore_ascii_case("na") {
        return Ok(None);
    }
    let mut value: f64 = raw
        .parse()
        .with_context(|| format!("parsing abundance value {raw:?}"))?;
    if !value.is_finite() {
        return Ok(None);
    }
    if log2_transform {
        if value <= 0.0 {
            return Ok(None);
        }
        value = value.log2();
    }
    Ok(Some(value))
}

fn write_proteins(path: &Path, platform: &Platform, proteins: &[ProteinRow]) -> Result<()> {
    let mut buf = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for p in proteins {
        buf.push_str(platform.as_str());
        buf.push('\t');
        buf.push_str(&p.assay_id);
        buf.push('\t');
        buf.push_str(&p.uniprot);
        buf.push('\t');
        buf.push_str(&p.gene_symbol);
        buf.push('\t');
        buf.push_str(&p.panel);
        buf.push_str("\t\n");
    }
    atomic_write(path, buf.as_bytes())
}

fn write_measurements(path: &Path, rows: &[MeasurementRow], abundance_unit: &str) -> Result<()> {
    let mut buf = String::new();
    buf.push_str(&MEASUREMENT_HEADER.join("\t"));
    buf.push('\n');
    for r in rows {
        let abundance = r.abundance.map(format_float).unwrap_or_default();
        buf.push_str(r.platform.as_str());
        buf.push('\t');
        buf.push_str(&r.sample_id);
        buf.push('\t');
        buf.push_str(&r.assay_id);
        buf.push('\t');
        buf.push_str(&r.gene_symbol);
        buf.push('\t');
        buf.push_str(&r.panel);
        buf.push('\t');
        buf.push_str(&r.source);
        buf.push('\t');
        buf.push_str(&abundance);
        buf.push('\t');
        buf.push_str(&abundance);
        buf.push('\t');
        buf.push_str(abundance_unit);
        buf.push_str("\tPASS\tPASS\t\t0\t");
        buf.push_str(&(r.dropped_by_qc as u8).to_string());
        buf.push_str("\t\t\t");
        buf.push_str(&r.ingest_order.to_string());
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

impl Table {
    fn need(&self, column: &str) -> Result<usize> {
        self.index
            .get(column)
            .copied()
            .with_context(|| format!("missing column {column:?}"))
    }
}

fn cell(row: &[String], idx: usize) -> &str {
    row.get(idx).map(String::as_str).unwrap_or("")
}

fn parse_is_control(value: &str) -> bool {
    matches!(
        value.trim().to_ascii_lowercase().as_str(),
        "1" | "true" | "yes" | "y"
    )
}

fn format_float(value: f64) -> String {
    format!("{value}")
}

/// Emit a stderr warning when per-sample medians span more than 1.0 log2
/// unit (≈2× fold). Only meaningful for log-scale inputs; on linear data
/// one log2 unit of median spread is still a strong sign of loading bias.
fn warn_if_unnormalized(rows: &[MeasurementRow]) {
    let medians = per_sample_medians(rows);
    if medians.len() < 2 {
        return;
    }
    let min = medians.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = medians.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let spread = max - min;
    if spread.is_finite() && spread > 1.0 {
        eprintln!(
            "ingest-matrix: WARNING per-sample median abundance spans {spread:.2} log2 units \
across {} samples (max={max:.3}, min={min:.3}). This pattern usually reflects \
un-normalized loading or intensity differences across samples. Re-run with \
`--normalize median` (or `--normalize quantile`), normalize upstream, or pass \
`--skip-normalization-check` to suppress this warning.",
            medians.len()
        );
    }
}

fn per_sample_medians(rows: &[MeasurementRow]) -> Vec<f64> {
    let mut by_sample: HashMap<String, Vec<f64>> = HashMap::new();
    for r in rows.iter() {
        if let Some(v) = r.abundance {
            by_sample.entry(r.sample_id.clone()).or_default().push(v);
        }
    }
    by_sample
        .into_values()
        .map(|mut v| median_in_place(&mut v))
        .filter(|m| m.is_finite())
        .collect()
}

fn median_in_place(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    }
}

fn apply_normalization(rows: &mut [MeasurementRow], method: Normalize) -> Result<()> {
    match method {
        Normalize::None => Ok(()),
        Normalize::Median => {
            apply_median_normalization(rows);
            Ok(())
        }
        Normalize::Quantile => apply_quantile_normalization(rows),
    }
}

fn apply_median_normalization(rows: &mut [MeasurementRow]) {
    let mut by_sample: HashMap<String, Vec<f64>> = HashMap::new();
    for r in rows.iter() {
        if let Some(v) = r.abundance {
            by_sample.entry(r.sample_id.clone()).or_default().push(v);
        }
    }
    if by_sample.is_empty() {
        return;
    }
    let mut sample_medians: HashMap<String, f64> = HashMap::new();
    for (sid, mut vs) in by_sample.into_iter() {
        let m = median_in_place(&mut vs);
        if m.is_finite() {
            sample_medians.insert(sid, m);
        }
    }
    if sample_medians.is_empty() {
        return;
    }
    let grand: f64 = sample_medians.values().copied().sum::<f64>() / sample_medians.len() as f64;
    for r in rows.iter_mut() {
        if let Some(v) = r.abundance.as_mut() {
            if let Some(&m) = sample_medians.get(&r.sample_id) {
                *v = *v - m + grand;
            }
        }
    }
}

fn apply_quantile_normalization(rows: &mut [MeasurementRow]) -> Result<()> {
    // Build sample and assay orderings.
    let mut samples: BTreeSet<String> = BTreeSet::new();
    let mut assays: BTreeSet<String> = BTreeSet::new();
    for r in rows.iter() {
        samples.insert(r.sample_id.clone());
        assays.insert(r.assay_id.clone());
    }
    let sample_list: Vec<String> = samples.into_iter().collect();
    let assay_list: Vec<String> = assays.into_iter().collect();
    let sample_idx: HashMap<String, usize> = sample_list
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), i))
        .collect();
    let assay_idx: HashMap<String, usize> = assay_list
        .iter()
        .enumerate()
        .map(|(i, a)| (a.clone(), i))
        .collect();
    let n_s = sample_list.len();
    let n_a = assay_list.len();
    if n_s < 2 {
        return Ok(());
    }
    let mut grid: Vec<Vec<Option<f64>>> = vec![vec![None; n_s]; n_a];
    for r in rows.iter() {
        let s = sample_idx[&r.sample_id];
        let a = assay_idx[&r.assay_id];
        grid[a][s] = r.abundance;
    }
    // Complete-case: assays with a value in every sample.
    let complete_rows: Vec<usize> = (0..n_a)
        .filter(|&a| grid[a].iter().all(|v| v.is_some()))
        .collect();
    let n_c = complete_rows.len();
    if n_c < 2 {
        eprintln!(
            "ingest-matrix: quantile normalization fell back to median centering for every assay \
(only {n_c} complete-case assay(s) across {n_s} samples)."
        );
        apply_median_normalization(rows);
        return Ok(());
    }
    // Per-sample sorted vectors and argsort (which complete-row idx sits at each rank).
    let mut sorted: Vec<Vec<f64>> = vec![Vec::with_capacity(n_c); n_s];
    let mut argsort: Vec<Vec<usize>> = vec![Vec::with_capacity(n_c); n_s];
    for s in 0..n_s {
        let mut pairs: Vec<(usize, f64)> = complete_rows
            .iter()
            .map(|&a| (a, grid[a][s].unwrap()))
            .collect();
        pairs.sort_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal));
        for (a_idx, v) in pairs {
            sorted[s].push(v);
            argsort[s].push(a_idx);
        }
    }
    // Target per rank = mean across samples of the rank-th sorted value.
    let target: Vec<f64> = (0..n_c)
        .map(|r| sorted.iter().map(|row| row[r]).sum::<f64>() / n_s as f64)
        .collect();
    // Write target back into the grid for complete-case assays.
    let mut new_values: HashMap<(usize, usize), f64> = HashMap::new();
    for (s, ranks) in argsort.iter().enumerate() {
        for (r, &a_idx) in ranks.iter().enumerate() {
            new_values.insert((a_idx, s), target[r]);
        }
    }
    // Median-fallback for partial-case assays: compute per-sample medians
    // using ALL non-dropped values in grid (stable reference).
    let sample_medians: Vec<f64> = (0..n_s)
        .map(|s| {
            let mut vs: Vec<f64> = grid.iter().filter_map(|row| row[s]).collect();
            median_in_place(&mut vs)
        })
        .collect();
    let finite_med: Vec<f64> = sample_medians
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    let grand: f64 = if finite_med.is_empty() {
        0.0
    } else {
        finite_med.iter().sum::<f64>() / finite_med.len() as f64
    };
    let n_partial = n_a - n_c;
    if n_partial > 0 {
        eprintln!(
            "ingest-matrix: quantile-normalized {n_c} complete-case assay(s); \
median-centered {n_partial} partial-case assay(s)."
        );
    }
    for r in rows.iter_mut() {
        if let Some(v) = r.abundance.as_mut() {
            let s = sample_idx[&r.sample_id];
            let a = assay_idx[&r.assay_id];
            if let Some(nv) = new_values.get(&(a, s)) {
                *v = *nv;
            } else {
                let m = sample_medians[s];
                if m.is_finite() {
                    *v = *v - m + grand;
                }
            }
        }
    }
    Ok(())
}
