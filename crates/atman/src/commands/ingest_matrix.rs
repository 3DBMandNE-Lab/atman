use anyhow::{bail, Context, Result};
use atman_core::Platform;
use clap::{Args as ClapArgs, ValueEnum};
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use crate::io::atomic_write;

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

    /// Do not copy measurements.tsv to qc_measurements.tsv.
    #[arg(long, default_value_t = false)]
    no_copy_measurements_to_qc: bool,
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
    let platform: Platform = args.platform.parse().map_err(anyhow::Error::msg)?;
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let matrix = read_table(&args.matrix)?;
    let samples = read_table(&args.samples)?;
    let sample_ids = write_samples(&args, &samples)?;

    let (proteins, measurements) = match args.orientation {
        Orientation::ProteinsRows => {
            let proteins = protein_rows_from_table(&args, &matrix)?;
            let measurements =
                measurements_from_proteins_rows(&args, platform, &matrix, &sample_ids)?;
            (proteins, measurements)
        }
        Orientation::SamplesRows => {
            let protein_path = args
                .proteins
                .as_ref()
                .context("--orientation samples-rows requires --proteins")?;
            let protein_meta = read_table(protein_path)?;
            let proteins = protein_rows_from_table(&args, &protein_meta)?;
            let measurements = measurements_from_samples_rows(&args, platform, &matrix, &proteins)?;
            (proteins, measurements)
        }
    };

    write_proteins(&args.output_dir.join("proteins.tsv"), platform, &proteins)?;
    write_measurements(
        &args.output_dir.join("measurements.tsv"),
        &measurements,
        &args.abundance_unit,
    )?;
    if !args.no_copy_measurements_to_qc {
        write_measurements(
            &args.output_dir.join("qc_measurements.tsv"),
            &measurements,
            &args.abundance_unit,
        )?;
    }

    eprintln!(
        "ingest-matrix: {} samples, {} proteins, {} measurements",
        sample_ids.len(),
        proteins.len(),
        measurements.len()
    );
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
    platform: Platform,
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
    platform: Platform,
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
    platform: Platform,
    sample_id: &str,
    protein: &ProteinRow,
    source: &str,
    log2_transform: bool,
    ingest_order: u64,
) -> Result<MeasurementRow> {
    let abundance = parse_abundance(source, log2_transform)?;
    Ok(MeasurementRow {
        platform,
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

fn write_proteins(path: &Path, platform: Platform, proteins: &[ProteinRow]) -> Result<()> {
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
