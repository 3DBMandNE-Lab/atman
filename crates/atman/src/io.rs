//! File I/O helpers: atomic writes, long-TSV read/write, CSV write.
//!
//! The long TSV is header-indexed on read (not positional) so adding columns
//! in the future does not silently mis-read older files. The writer emits a
//! stable column order but the reader tolerates either ordering.

use anyhow::{anyhow, bail, Context, Result};
use atman_core::{
    fold_change::FoldChangePanel, matrix::WidePanel, Abundance, AssayId, Batch, DetectionLimit,
    MeasurementRecord, PeptideIdentity, PeptideMeasurementRecord, Platform, ProteinIdentity,
    QcFlag, Sample,
};
use csv::StringRecord;
use sha2::{Digest, Sha256};
use std::{
    collections::HashMap,
    io::Write,
    path::{Path, PathBuf},
};

pub use crate::run_sidecar::{
    format_iso8601_utc, hash_canonical_inputs, hash_labeled_inputs, sidecar_path_for,
    write_run_sidecar,
};

/// Locate a required column by header name; bail with a readable error when missing.
pub fn need_col(headers: &StringRecord, name: &str, path: &Path) -> Result<usize> {
    headers
        .iter()
        .position(|h| h == name)
        .ok_or_else(|| anyhow!("missing column {:?} in {:?}", name, path))
}

/// Locate the first matching column from a list of candidate names.
pub fn find_col(headers: &StringRecord, names: &[&str]) -> Option<usize> {
    names
        .iter()
        .find_map(|name| headers.iter().position(|h| h == *name))
}

/// Fetch an optional cell value, treating empty strings as absent.
pub fn optional_cell(row: &StringRecord, col: Option<usize>) -> Option<&str> {
    col.and_then(|idx| row.get(idx))
        .filter(|value| !value.is_empty())
}

/// Fixed-6-decimal f64 formatter for human-facing report/QC columns
/// (e.g. report.rs fraction columns, decompose `majority_sign`). Handles
/// `NaN`/`Inf` explicitly and collapses zero to `"0"`. This rounding loses
/// precision and is NOT round-trip-safe — use [`format_f64`] for value columns
/// that must reparse to the exact same `f64`.
/// Format a float for TSV output at SIX DECIMALS.
///
/// This is the narrower of atman's two float writers, and the pair is a
/// footgun because the names look interchangeable and are not:
///
/// | Writer | Precision | Used by |
/// |---|---|---|
/// | `format_float` | `{:.6}` | `decompose`, `align`, `harmonize`, `bootstrap`, most analysis commands |
/// | [`format_f64`] | `{}`, shortest round-trip | the `axes` family, `de/output_rows`, `enrich`, `concordance`, `residuals`, `scale`, `score_weighted`, `report` |
///
/// At six decimals a value read back from the file is not the value that
/// was computed: `7.5/8.5` writes as `0.882353` and parses back 5.9e-8
/// away from the original. That is fine for a report and wrong for a
/// round-trip, so a test that parses this output must not assert a
/// tolerance tighter than about 1e-6 — and tighter still if it sums
/// several values, since the rounding accumulates.
///
/// It also means a byte-comparison of two files written by THIS function
/// is an agreement test at 1e-6, not a bit-identity test. Files written
/// by [`format_f64`] compare at full precision. When citing a
/// byte-comparison as evidence, say which writer produced the file.
///
/// Zero is written as `0` rather than `0.000000`, and non-finite values
/// as `NaN` / `Inf` / `-Inf`.
pub fn format_float(value: f64) -> String {
    if !value.is_finite() {
        return if value.is_nan() {
            "NaN".to_string()
        } else if value > 0.0 {
            "Inf".to_string()
        } else {
            "-Inf".to_string()
        };
    }
    if value == 0.0 {
        "0".to_string()
    } else {
        format!("{value:.6}")
    }
}

/// Replace tab / newline / carriage-return characters with spaces so the value
/// fits safely into a TSV cell.
pub fn escape_tsv(value: &str) -> String {
    value.replace(['\t', '\n', '\r'], " ")
}

/// SHA-256 of a byte slice as lower-case hex.
pub fn sha256_hex(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    let digest = hasher.finalize();
    let mut out = String::with_capacity(digest.len() * 2);
    for b in digest {
        out.push_str(&format!("{b:02x}"));
    }
    out
}

/// Atomic file write: write to a sibling tempfile, then rename over target.
pub fn atomic_write(path: &Path, bytes: &[u8]) -> Result<()> {
    let parent = path.parent().unwrap_or(Path::new("."));
    let mut tmp = tempfile::NamedTempFile::new_in(parent)
        .with_context(|| format!("creating tempfile in {:?}", parent))?;
    tmp.write_all(bytes)
        .with_context(|| format!("writing tempfile for {:?}", path))?;
    tmp.persist(path)
        .with_context(|| format!("atomic rename to {:?}", path))?;
    Ok(())
}

/// Streaming variant of [`atomic_write`]: `write` receives a buffered
/// writer over a tempfile in the target directory; on success the
/// tempfile is renamed into place, on error nothing is left at `path`.
/// Use it when the content is too large to assemble in memory first.
pub fn atomic_write_with(
    path: &Path,
    write: impl FnOnce(&mut dyn Write) -> Result<()>,
) -> Result<()> {
    let parent = path.parent().unwrap_or(Path::new("."));
    let tmp = tempfile::NamedTempFile::new_in(parent)
        .with_context(|| format!("creating tempfile in {:?}", parent))?;
    {
        let mut buf = std::io::BufWriter::with_capacity(1 << 20, tmp.as_file());
        write(&mut buf)?;
        buf.flush()
            .with_context(|| format!("flushing tempfile for {:?}", path))?;
    }
    tmp.persist(path)
        .with_context(|| format!("atomic rename to {:?}", path))?;
    Ok(())
}

const LONG_TSV_HEADER: &[&str] = &[
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

/// Write the canonical long-format measurements TSV.
pub fn write_measurements_long(path: &Path, records: &[MeasurementRecord]) -> Result<()> {
    let mut buf = String::new();
    buf.push_str(&LONG_TSV_HEADER.join("\t"));
    buf.push('\n');
    for r in records {
        let abundance = if r.dropped_by_qc {
            String::new()
        } else {
            format_f64(r.abundance.as_f64())
        };
        let abundance_raw = format_f64(r.abundance_raw.as_f64());
        let lod = match r.detection_limit.0 {
            Some(v) => format_f64(v),
            None => String::new(),
        };
        let plate = r.batch.plate.clone().unwrap_or_default();
        let panel = r.panel.clone().unwrap_or_default();
        let lot = r.batch.lot.clone().unwrap_or_default();
        let gene = r.gene_symbol.clone().unwrap_or_default();
        buf.push_str(r.platform.as_str());
        buf.push('\t');
        buf.push_str(&r.sample_id);
        buf.push('\t');
        buf.push_str(&r.assay_id.0);
        buf.push('\t');
        buf.push_str(&gene);
        buf.push('\t');
        buf.push_str(&panel);
        buf.push('\t');
        buf.push_str(&r.npx_source_str);
        buf.push('\t');
        buf.push_str(&abundance);
        buf.push('\t');
        buf.push_str(&abundance_raw);
        buf.push('\t');
        buf.push_str(abundance_unit(&r.abundance_raw));
        buf.push('\t');
        buf.push_str(r.qc_sample.as_str());
        buf.push('\t');
        buf.push_str(r.qc_assay.as_str());
        buf.push('\t');
        buf.push_str(&lod);
        buf.push('\t');
        buf.push_str(&(r.below_lod as u8).to_string());
        buf.push('\t');
        buf.push_str(&(r.dropped_by_qc as u8).to_string());
        buf.push('\t');
        buf.push_str(&plate);
        buf.push('\t');
        buf.push_str(&lot);
        buf.push('\t');
        buf.push_str(&r.ingest_order.to_string());
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

/// Read the long-format measurements TSV, indexing columns by header name.
pub fn read_measurements_long(path: &Path) -> Result<Vec<MeasurementRecord>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;

    let headers = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .clone();
    let col: HashMap<String, usize> = headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.to_string(), i))
        .collect();
    let need = |name: &str| -> Result<usize> {
        col.get(name)
            .copied()
            .ok_or_else(|| anyhow!("missing column {:?} in {:?}", name, path))
    };
    let c_sample = need("sample_id")?;
    let c_assay = need("assay_id")?;
    let c_gene = need("gene_symbol")?;
    let c_panel = need("panel")?;
    let c_src = need("npx_source_str")?;
    let c_abund = need("abundance")?;
    let c_abund_raw = need("abundance_raw")?;
    let c_unit = need("abundance_unit")?;
    let c_qcs = need("qc_sample")?;
    let c_qca = need("qc_assay")?;
    let c_lod = need("detection_limit")?;
    let c_below = need("below_lod")?;
    let c_drop = need("dropped_by_qc")?;
    let c_plate = need("plate_id")?;
    let c_lot = need("panel_lot")?;
    let c_order = need("ingest_order")?;
    let c_platform = col.get("platform").copied();

    let mut out = Vec::new();
    for result in reader.records() {
        let row = result.with_context(|| format!("reading record from {:?}", path))?;
        let unit = &row[c_unit];
        let dropped_by_qc: bool = row[c_drop]
            .parse::<u8>()
            .map(|v| v != 0)
            .with_context(|| format!("parsing dropped_by_qc as 0/1 in {:?}", path))?;
        let abund_str = &row[c_abund];
        let abundance_raw_value: f64 = if row[c_abund_raw].is_empty() {
            if dropped_by_qc {
                f64::NAN
            } else {
                anyhow::bail!("abundance_raw parse: empty value on non-dropped row");
            }
        } else {
            row[c_abund_raw].parse().context("abundance_raw parse")?
        };
        // NOTE: a non-finite (`NaN`) value in the `abundance` / `abundance_raw`
        // columns is INTENTIONAL — it is the canonical MNAR/below-LOD missingness
        // marker that the missingness-aware paths (e.g. missingness-ICA) read and
        // model. Finite-assuming consumers must go through
        // `MeasurementRecord::effective_abundance()`, which returns `None` for
        // non-finite values; the reader therefore does NOT reject them.
        let abundance = if abund_str.is_empty() {
            abundance_from_unit(unit, abundance_raw_value)
        } else {
            abundance_from_unit(unit, abund_str.parse().context("abundance parse")?)
        };
        let abundance_raw = abundance_from_unit(unit, abundance_raw_value);
        let detection_limit = if row[c_lod].is_empty() {
            DetectionLimit(None)
        } else {
            DetectionLimit(Some(row[c_lod].parse().context("detection_limit parse")?))
        };
        let below_lod: bool = row[c_below]
            .parse::<u8>()
            .map(|v| v != 0)
            .with_context(|| format!("parsing below_lod as 0/1 in {:?}", path))?;
        let platform = match c_platform
            .and_then(|idx| row.get(idx))
            .filter(|s| !s.is_empty())
        {
            Some(value) => value.parse::<Platform>().map_err(anyhow::Error::msg)?,
            None => Platform::OlinkExploreNgs,
        };
        out.push(MeasurementRecord {
            platform,
            assay_id: AssayId(row[c_assay].to_string()),
            gene_symbol: empty_to_none(&row[c_gene]),
            sample_id: row[c_sample].to_string(),
            abundance,
            abundance_raw,
            npx_source_str: row[c_src].to_string(),
            qc_sample: parse_qc(&row[c_qcs])
                .with_context(|| format!("parsing qc_sample in {:?}", path))?,
            qc_assay: parse_qc(&row[c_qca])
                .with_context(|| format!("parsing qc_assay in {:?}", path))?,
            detection_limit,
            below_lod,
            batch: Batch {
                plate: empty_to_none(&row[c_plate]),
                lot: empty_to_none(&row[c_lot]),
                run: None,
            },
            dropped_by_qc,
            ingest_order: row[c_order]
                .parse()
                .with_context(|| format!("parsing ingest_order as integer in {:?}", path))?,
            panel: empty_to_none(&row[c_panel]),
        });
    }
    Ok(out)
}

/// samples.tsv writer.
pub fn write_samples(path: &Path, samples: &[Sample]) -> Result<()> {
    let mut buf =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for s in samples {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            s.sample_id,
            s.subject_id.clone().unwrap_or_default(),
            s.condition.clone().unwrap_or_default(),
            s.is_control as u8,
            s.sample_type.clone().unwrap_or_default(),
            s.ingest_order,
        ));
    }
    atomic_write(path, buf.as_bytes())
}

pub fn read_samples(path: &Path) -> Result<Vec<Sample>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .clone();
    let col: HashMap<String, usize> = headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.to_string(), i))
        .collect();
    let need = |name: &str| -> Result<usize> {
        col.get(name)
            .copied()
            .ok_or_else(|| anyhow!("missing column {:?} in {:?}", name, path))
    };
    let c_sample = need("sample_id")?;
    let c_subject = need("subject_id")?;
    let c_condition = need("condition")?;
    let c_control = need("is_control")?;
    let c_type = need("sample_type")?;
    let c_order = need("ingest_order")?;
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row.with_context(|| format!("reading record from {:?}", path))?;
        out.push(Sample {
            sample_id: row[c_sample].to_string(),
            subject_id: empty_to_none(&row[c_subject]),
            condition: empty_to_none(&row[c_condition]),
            is_control: row[c_control]
                .parse::<u8>()
                .with_context(|| format!("is_control parse in {:?}", path))?
                != 0,
            sample_type: empty_to_none(&row[c_type]),
            ingest_order: row[c_order]
                .parse()
                .with_context(|| format!("parsing ingest_order as integer in {:?}", path))?,
        });
    }
    Ok(out)
}

/// proteins.tsv reader, returning `ProteinIdentity` rows. Header-indexed.
pub fn read_proteins(path: &Path) -> Result<Vec<ProteinIdentity>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .clone();
    let col: HashMap<String, usize> = headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.to_string(), i))
        .collect();
    let need = |name: &str| -> Result<usize> {
        col.get(name)
            .copied()
            .ok_or_else(|| anyhow!("missing column {:?} in {:?}", name, path))
    };
    let c_platform = need("platform")?;
    let c_assay = need("assay_id")?;
    let c_uniprot = need("uniprot")?;
    let c_gene = need("gene_symbol")?;
    let c_panel = need("panel")?;
    let c_lot = need("panel_lot")?;

    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        let uniprot: Vec<String> = if row[c_uniprot].is_empty() {
            vec![]
        } else {
            row[c_uniprot].split(',').map(|s| s.to_string()).collect()
        };
        out.push(ProteinIdentity {
            platform: row[c_platform]
                .parse::<Platform>()
                .map_err(anyhow::Error::msg)?,
            assay_id: AssayId(row[c_assay].to_string()),
            uniprot,
            gene_symbol: empty_to_none(&row[c_gene]),
            panel: empty_to_none(&row[c_panel]),
            panel_lot: empty_to_none(&row[c_lot]),
        });
    }
    Ok(out)
}

/// peptides.tsv reader: peptide catalog keyed by `peptide_id`, with the
/// parent protein's `assay_id` required and sequence/charge/modification
/// columns optional. Header-indexed.
pub fn read_peptides(path: &Path) -> Result<Vec<PeptideIdentity>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .clone();
    let col: HashMap<String, usize> = headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.to_string(), i))
        .collect();
    let need = |name: &str| -> Result<usize> {
        col.get(name)
            .copied()
            .ok_or_else(|| anyhow!("missing column {:?} in {:?}", name, path))
    };
    let c_pep = need("peptide_id")?;
    let c_assay = need("assay_id")?;
    let c_seq = col.get("sequence").copied();
    let c_charge = col.get("charge").copied();
    let c_mods = col.get("modifications").copied();
    let c_mc = col.get("missed_cleavages").copied();
    let mut out = Vec::new();
    for result in reader.records() {
        let row = result.with_context(|| format!("reading peptide row from {:?}", path))?;
        let charge = match c_charge.and_then(|i| row.get(i)) {
            Some(s) if !s.is_empty() => Some(
                s.parse::<i32>()
                    .with_context(|| format!("charge parse in {:?}", path))?,
            ),
            _ => None,
        };
        let missed = match c_mc.and_then(|i| row.get(i)) {
            Some(s) if !s.is_empty() => Some(
                s.parse::<u32>()
                    .with_context(|| format!("missed_cleavages parse in {:?}", path))?,
            ),
            _ => None,
        };
        out.push(PeptideIdentity {
            peptide_id: row[c_pep].to_string(),
            assay_id: AssayId(row[c_assay].to_string()),
            sequence: c_seq.and_then(|i| row.get(i)).and_then(empty_to_none),
            charge,
            modifications: c_mods.and_then(|i| row.get(i)).and_then(empty_to_none),
            missed_cleavages: missed,
        });
    }
    Ok(out)
}

/// peptides.tsv writer.
pub fn write_peptides(path: &Path, peptides: &[PeptideIdentity]) -> Result<()> {
    let mut buf =
        String::from("peptide_id\tassay_id\tsequence\tcharge\tmodifications\tmissed_cleavages\n");
    for p in peptides {
        buf.push_str(&p.peptide_id);
        buf.push('\t');
        buf.push_str(&p.assay_id.0);
        buf.push('\t');
        buf.push_str(&p.sequence.clone().unwrap_or_default());
        buf.push('\t');
        if let Some(c) = p.charge {
            buf.push_str(&c.to_string());
        }
        buf.push('\t');
        buf.push_str(&p.modifications.clone().unwrap_or_default());
        buf.push('\t');
        if let Some(m) = p.missed_cleavages {
            buf.push_str(&m.to_string());
        }
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

/// peptide_measurements.tsv reader. One row per (sample, peptide) pair.
/// Header-indexed.
pub fn read_peptide_measurements(path: &Path) -> Result<Vec<PeptideMeasurementRecord>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .clone();
    let col: HashMap<String, usize> = headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.to_string(), i))
        .collect();
    let need = |name: &str| -> Result<usize> {
        col.get(name)
            .copied()
            .ok_or_else(|| anyhow!("missing column {:?} in {:?}", name, path))
    };
    let c_sample = need("sample_id")?;
    let c_pep = need("peptide_id")?;
    let c_abund = need("abundance")?;
    let c_unit = need("abundance_unit")?;
    let c_drop = need("dropped_by_qc")?;
    let c_below = need("below_lod")?;
    let mut out = Vec::new();
    for result in reader.records() {
        let row = result.with_context(|| format!("reading peptide measurement from {:?}", path))?;
        let dropped: bool = row[c_drop]
            .parse::<u8>()
            .map(|v| v != 0)
            .with_context(|| format!("parsing dropped_by_qc as 0/1 in {:?}", path))?;
        let abund_str = &row[c_abund];
        let abundance: f64 = if abund_str.is_empty() {
            f64::NAN
        } else {
            abund_str
                .parse()
                .with_context(|| format!("peptide abundance parse in {:?}", path))?
        };
        let below_lod: bool = row[c_below]
            .parse::<u8>()
            .map(|v| v != 0)
            .with_context(|| format!("parsing below_lod as 0/1 in {:?}", path))?;
        out.push(PeptideMeasurementRecord {
            sample_id: row[c_sample].to_string(),
            peptide_id: row[c_pep].to_string(),
            abundance,
            abundance_unit: row[c_unit].to_string(),
            dropped_by_qc: dropped,
            below_lod,
        });
    }
    Ok(out)
}

/// peptide_measurements.tsv writer.
pub fn write_peptide_measurements(path: &Path, records: &[PeptideMeasurementRecord]) -> Result<()> {
    let mut buf = String::from(
        "sample_id\tpeptide_id\tabundance\tabundance_unit\tdropped_by_qc\tbelow_lod\n",
    );
    for r in records {
        buf.push_str(&r.sample_id);
        buf.push('\t');
        buf.push_str(&r.peptide_id);
        buf.push('\t');
        if r.dropped_by_qc || !r.abundance.is_finite() {
            // masked — leave cell empty so reader sees NaN
        } else {
            buf.push_str(&format!("{}", r.abundance));
        }
        buf.push('\t');
        buf.push_str(&r.abundance_unit);
        buf.push('\t');
        buf.push_str(&(r.dropped_by_qc as u8).to_string());
        buf.push('\t');
        buf.push_str(&(r.below_lod as u8).to_string());
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

/// proteins.tsv writer.
pub fn write_proteins(path: &Path, proteins: &[ProteinIdentity]) -> Result<()> {
    let mut buf = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for p in proteins {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            p.platform.as_str(),
            p.assay_id.0,
            p.uniprot.join(","),
            p.gene_symbol.clone().unwrap_or_default(),
            p.panel.clone().unwrap_or_default(),
            p.panel_lot.clone().unwrap_or_default(),
        ));
    }
    atomic_write(path, buf.as_bytes())
}

/// Wide-panel CSV writer. Starts from the source NPX string (no f64
/// round-trip) and applies the trailing-zero trim documented on
/// [`trim_npx_string`]: `0.0980` → `0.098`, `1.4200` → `1.42`,
/// `0.4000` → `0.4`. This matches the formatting rule observed
/// empirically in the published Olink Explore NPX files used as the
/// byte-exact reproduction reference.
/// Sanitize a *data-derived* string (e.g. a `panel` value read from an
/// untrusted `measurements.tsv`) for safe use as a single output-filename
/// component. Without this, a crafted panel like `../../evil` would make
/// `out_dir.join(...)` resolve outside `out_dir`. Strict-failure: reject any
/// value carrying a path separator or directory-traversal rather than silently
/// rewriting it. A legitimate panel name (e.g. `Cardiometabolic`) is unaffected.
fn safe_filename_component(raw: &str) -> Result<String> {
    let s = raw.trim();
    if s.is_empty() {
        bail!("empty panel name cannot form an output filename");
    }
    if s == "." || s == ".." || s.contains('/') || s.contains('\\') || s.contains('\0') {
        bail!(
            "panel name {raw:?} is not a safe output-filename component \
             (contains a path separator or traversal)"
        );
    }
    Ok(s.to_ascii_lowercase())
}

pub fn write_wide_panel(out_dir: &Path, panel: &WidePanel) -> Result<PathBuf> {
    let filename = format!("{}_npx.csv", safe_filename_component(&panel.panel)?);
    let path = out_dir.join(filename);
    let mut buf = String::new();
    buf.push_str("Participant,Exposure,SampleID");
    for a in &panel.assays {
        buf.push(',');
        buf.push_str(a);
    }
    buf.push('\n');
    for row in &panel.rows {
        buf.push_str(&row.participant);
        buf.push(',');
        buf.push_str(&row.exposure);
        buf.push(',');
        buf.push_str(&row.sample_id);
        for v in &row.values {
            buf.push(',');
            if let Some(s) = v {
                buf.push_str(&trim_npx_string(s));
            }
        }
        buf.push('\n');
    }
    atomic_write(&path, buf.as_bytes())?;
    Ok(path)
}

/// Dube's trailing-zero trim: strips trailing `0` characters from the
/// fractional part, but always leaves **at least one digit after the decimal
/// point**. Examples (verified against Dube's published filtered NPX):
///   `0.0980` → `0.098`
///   `0.4000` → `0.4`
///   `1.4200` → `1.42`
///   `0.0000` → `0.0`   (NOT `0`)
///   `-0.0000` → `-0.0`
fn trim_npx_string(s: &str) -> String {
    let dot_pos = match s.find('.') {
        Some(p) => p,
        None => return s.to_string(),
    };
    let mut end = s.len();
    while end > dot_pos + 2 && s.as_bytes()[end - 1] == b'0' {
        end -= 1;
    }
    s[..end].to_string()
}

#[cfg(test)]
mod trim_tests {
    use super::trim_npx_string;
    #[test]
    fn trims_trailing_zero_keeping_nonzero_fractional() {
        assert_eq!(trim_npx_string("0.0980"), "0.098");
        assert_eq!(trim_npx_string("0.4000"), "0.4");
        assert_eq!(trim_npx_string("1.4200"), "1.42");
    }
    #[test]
    fn preserves_one_fractional_digit_for_all_zeros() {
        assert_eq!(trim_npx_string("0.0000"), "0.0");
        assert_eq!(trim_npx_string("-0.0000"), "-0.0");
        assert_eq!(trim_npx_string("7.0000"), "7.0");
    }
    #[test]
    fn leaves_values_without_trailing_zeros_unchanged() {
        assert_eq!(trim_npx_string("0.5777"), "0.5777");
        assert_eq!(trim_npx_string("-0.1219"), "-0.1219");
        assert_eq!(trim_npx_string("7.4289"), "7.4289");
    }
    #[test]
    fn integer_strings_pass_through() {
        assert_eq!(trim_npx_string("42"), "42");
        assert_eq!(trim_npx_string("-0"), "-0");
    }
}

/// Fold-change panel CSV writer.
pub fn write_fold_change_panel(
    out_dir: &Path,
    panel: &FoldChangePanel,
    comparisons: &[atman_core::fold_change::Comparison],
) -> Result<PathBuf> {
    let filename = format!("{}_log2_fc.csv", safe_filename_component(&panel.panel)?);
    let path = out_dir.join(filename);
    let mut buf = String::from("Assay");
    for c in comparisons {
        buf.push(',');
        buf.push_str(&format!("{}-{}", c.a, c.b));
    }
    buf.push('\n');
    for (i, assay) in panel.assays.iter().enumerate() {
        buf.push_str(assay);
        for v in &panel.values[i] {
            buf.push(',');
            if let Some(x) = v {
                buf.push_str(&format_f64(*x));
            }
        }
        buf.push('\n');
    }
    atomic_write(&path, buf.as_bytes())?;
    Ok(path)
}

/// Parse a QC-flag cell. Strict: only the canonical spellings emitted by
/// [`QcFlag::as_str`] (`PASS`/`WARN`/`FAIL`) are accepted. Any other value —
/// a lowercase, truncated, or garbled cell — is an error rather than a silent
/// fall-through to `Pass`, which would let a failed sample masquerade as passed.
fn parse_qc(s: &str) -> Result<QcFlag> {
    match s {
        // Empty is the canonical "no flag specified" form and means Pass; it
        // is accepted by validate_measurements too. Strict-failure applies to
        // *unrecognized non-empty* values only (typos, lowercase, garbage).
        "PASS" | "" => Ok(QcFlag::Pass),
        "WARN" => Ok(QcFlag::Warn(String::new())),
        "FAIL" => Ok(QcFlag::Fail(String::new())),
        other => bail!("unrecognized QC flag {other:?} (expected PASS, WARN, FAIL, or empty)"),
    }
}

fn abundance_unit(abundance: &Abundance) -> &'static str {
    match abundance {
        Abundance::Log2Npx(_) => "log2_npx",
        Abundance::Log2Intensity(_) => "log2_intensity",
        Abundance::Ibaq(_) => "ibaq",
        Abundance::Raw(_) => "raw",
    }
}

fn abundance_from_unit(unit: &str, value: f64) -> Abundance {
    match unit.trim().to_ascii_lowercase().as_str() {
        "log2_npx" | "npx" => Abundance::Log2Npx(value),
        "log2_intensity"
        | "log2_lfq"
        | "log2_rfu"
        | "log2_pg_quantity"
        | "log2_diann_pg_quantity" => Abundance::Log2Intensity(value),
        "ibaq" => Abundance::Ibaq(value),
        _ => Abundance::Raw(value),
    }
}

fn empty_to_none(s: &str) -> Option<String> {
    if s.is_empty() {
        None
    } else {
        Some(s.to_string())
    }
}

/// Full-precision, round-trip f64 formatter for value columns (long-TSV
/// `abundance`/`abundance_raw`, fold-change panel cells). Uses Rust's default
/// `{}` formatting, which emits the shortest decimal string that parses back to
/// the exact same `f64`. Contrast with [`format_float`], whose fixed-6-decimal
/// rounding is for human-facing report columns and is NOT round-trip-safe.
/// Format a float for TSV output at FULL round-trip precision.
///
/// The wider of atman's two float writers. `{}` on an `f64` emits the
/// shortest decimal string that parses back to the same bits, so a value
/// written by this function survives a write/read cycle exactly.
///
/// See [`format_float`] for the six-decimal writer and for why the two
/// are easy to confuse. Do not swap one for the other to make a test
/// pass: the choice changes the precision of every column in the file,
/// and downstream byte-comparisons inherit it.
pub fn format_f64(v: f64) -> String {
    format!("{}", v)
}

/// One row of `de_results.tsv`. Every field that can be missing for a
/// skipped test is `Option<f64>`.
pub struct DeResultRow {
    pub panel: String,
    pub assay_id: String,
    pub gene_symbol: String,
    pub uniprot: String,
    pub comparison: String,
    pub n_pairs: usize,
    pub mean_a: Option<f64>,
    pub mean_b: Option<f64>,
    pub mean_diff: Option<f64>,
    pub t: Option<f64>,
    pub df: Option<f64>,
    pub p_value: Option<f64>,
    pub bh_q: Option<f64>,
    pub effect_size: Option<f64>,
    pub effect_size_method: String,
    pub ci_low: Option<f64>,
    pub ci_high: Option<f64>,
    pub wilcoxon_p: Option<f64>,
    pub wilcoxon_method: String,
    pub median_diff: Option<f64>,
    pub trimmed_mean_diff: Option<f64>,
    /// Non-empty when the test was skipped.
    pub skip_reason: String,
    // Appended by the limma eBayes + F-tests feature. Each Option<f64>
    // renders as an empty TSV cell when `None`.
    pub s2_trend: Option<f64>,
    pub s2_prior: Option<f64>,
    pub s2_posterior: Option<f64>,
    pub df_prior: Option<f64>,
    pub df_total: Option<f64>,
    pub f_statistic: Option<f64>,
    pub f_p_value: Option<f64>,
    pub f_bh_q: Option<f64>,
    pub lfc_threshold: Option<f64>,
    // Appended by the msqrob peptide-level ridge mixed model. `None`
    // for non-msqrob tests; empty TSV cells on write.
    pub n_peptides_observed: Option<usize>,
    pub peptide_variance_ratio: Option<f64>,
    pub ridge_lambda: Option<f64>,
    /// Emitting dispatch name (e.g. "paired-t", "welch-t", "ols",
    /// "mixed", "limma", "msqrob"). Populated for every row so the
    /// `--test ensemble` aggregator can bucket rows by method without
    /// parsing the free-form `effect_size_method` string.
    pub method: String,
    /// Post-hoc adjustment method name (e.g. "sidak", "tukey",
    /// "dunnett") for rows emitted by the `--post-hoc` path. Empty
    /// for every other row.
    pub posthoc_method: String,
    /// Raw per-contrast p-value before family-wise adjustment.
    /// Populated only for post-hoc rows.
    pub posthoc_p: Option<f64>,
    /// Family-wise adjusted p-value within the contrast family (e.g.
    /// Sidak `1 - (1-p)^m`). Populated only for post-hoc rows.
    pub posthoc_adj_p: Option<f64>,
    /// Observed subjects in group A (first side of the comparison) that
    /// entered the test. Filled for `welch-t` and the OLS path; empty for
    /// paired tests.
    pub n_a: Option<usize>,
    /// Observed subjects in group B.
    pub n_b: Option<usize>,
}

/// Parse `de_results.tsv` back into `DeResultRow`s. Used by
/// `atman de --test ensemble` to pull each sub-method's output out of
/// its tempdir for aggregation. Header-indexed so adding new columns
/// doesn't break older files.
pub fn read_de_results(path: &Path) -> Result<Vec<DeResultRow>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .clone();
    let col: HashMap<String, usize> = headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.to_string(), i))
        .collect();
    let get = |row: &csv::StringRecord, name: &str| -> String {
        col.get(name)
            .and_then(|&i| row.get(i))
            .unwrap_or("")
            .to_string()
    };
    let parse_opt_f64 = |row: &csv::StringRecord, name: &str| -> Result<Option<f64>> {
        let s = get(row, name);
        if s.is_empty() {
            Ok(None)
        } else {
            s.parse()
                .map(Some)
                .with_context(|| format!("parsing {name} as f64 in {:?}", path))
        }
    };
    let parse_opt_usize = |row: &csv::StringRecord, name: &str| -> Result<Option<usize>> {
        let s = get(row, name);
        if s.is_empty() {
            Ok(None)
        } else {
            s.parse()
                .map(Some)
                .with_context(|| format!("parsing {name} as integer in {:?}", path))
        }
    };
    let mut out = Vec::new();
    for result in reader.records() {
        let row = result.with_context(|| format!("reading de_results row from {:?}", path))?;
        out.push(DeResultRow {
            panel: get(&row, "panel"),
            assay_id: get(&row, "assay_id"),
            gene_symbol: get(&row, "gene_symbol"),
            uniprot: get(&row, "uniprot"),
            comparison: get(&row, "comparison"),
            n_pairs: get(&row, "n_pairs")
                .parse()
                .with_context(|| format!("parsing n_pairs as integer in {:?}", path))?,
            mean_a: parse_opt_f64(&row, "mean_a")?,
            mean_b: parse_opt_f64(&row, "mean_b")?,
            mean_diff: parse_opt_f64(&row, "mean_diff")?,
            t: parse_opt_f64(&row, "t")?,
            df: parse_opt_f64(&row, "df")?,
            p_value: parse_opt_f64(&row, "p_value")?,
            bh_q: parse_opt_f64(&row, "bh_q")?,
            effect_size: parse_opt_f64(&row, "effect_size")?,
            effect_size_method: get(&row, "effect_size_method"),
            ci_low: parse_opt_f64(&row, "ci_low")?,
            ci_high: parse_opt_f64(&row, "ci_high")?,
            wilcoxon_p: parse_opt_f64(&row, "wilcoxon_p")?,
            wilcoxon_method: get(&row, "wilcoxon_method"),
            median_diff: parse_opt_f64(&row, "median_diff")?,
            trimmed_mean_diff: parse_opt_f64(&row, "trimmed_mean_diff")?,
            skip_reason: get(&row, "skip_reason"),
            s2_trend: parse_opt_f64(&row, "s2_trend")?,
            s2_prior: parse_opt_f64(&row, "s2_prior")?,
            s2_posterior: parse_opt_f64(&row, "s2_posterior")?,
            df_prior: parse_opt_f64(&row, "df_prior")?,
            df_total: parse_opt_f64(&row, "df_total")?,
            f_statistic: parse_opt_f64(&row, "f_statistic")?,
            f_p_value: parse_opt_f64(&row, "f_p_value")?,
            f_bh_q: parse_opt_f64(&row, "f_bh_q")?,
            lfc_threshold: parse_opt_f64(&row, "lfc_threshold")?,
            n_peptides_observed: parse_opt_usize(&row, "n_peptides_observed")?,
            peptide_variance_ratio: parse_opt_f64(&row, "peptide_variance_ratio")?,
            ridge_lambda: parse_opt_f64(&row, "ridge_lambda")?,
            method: get(&row, "method"),
            posthoc_method: get(&row, "posthoc_method"),
            posthoc_p: parse_opt_f64(&row, "posthoc_p")?,
            posthoc_adj_p: parse_opt_f64(&row, "posthoc_adj_p")?,
            n_a: parse_opt_usize(&row, "n_a")?,
            n_b: parse_opt_usize(&row, "n_b")?,
        });
    }
    Ok(out)
}

pub fn write_de_results(path: &Path, rows: &[DeResultRow]) -> Result<()> {
    let mut buf = String::from(
        "panel\tassay_id\tgene_symbol\tuniprot\tcomparison\tn_pairs\t\
         mean_a\tmean_b\tmean_diff\tt\tdf\tp_value\tbh_q\tskip_reason\t\
         effect_size\teffect_size_method\tci_low\tci_high\twilcoxon_p\twilcoxon_method\t\
         median_diff\ttrimmed_mean_diff\t\
         s2_trend\ts2_prior\ts2_posterior\tdf_prior\tdf_total\t\
         f_statistic\tf_p_value\tf_bh_q\tlfc_threshold\t\
         n_peptides_observed\tpeptide_variance_ratio\tridge_lambda\tmethod\t\
         posthoc_method\tposthoc_p\tposthoc_adj_p\tn_a\tn_b\n",
    );
    for r in rows {
        buf.push_str(&r.panel);
        buf.push('\t');
        buf.push_str(&r.assay_id);
        buf.push('\t');
        buf.push_str(&r.gene_symbol);
        buf.push('\t');
        buf.push_str(&r.uniprot);
        buf.push('\t');
        buf.push_str(&r.comparison);
        buf.push('\t');
        buf.push_str(&r.n_pairs.to_string());
        buf.push('\t');
        push_opt_f64(&mut buf, r.mean_a);
        buf.push('\t');
        push_opt_f64(&mut buf, r.mean_b);
        buf.push('\t');
        push_opt_f64(&mut buf, r.mean_diff);
        buf.push('\t');
        push_opt_f64(&mut buf, r.t);
        buf.push('\t');
        push_opt_f64(&mut buf, r.df);
        buf.push('\t');
        push_opt_f64(&mut buf, r.p_value);
        buf.push('\t');
        push_opt_f64(&mut buf, r.bh_q);
        buf.push('\t');
        buf.push_str(&r.skip_reason);
        buf.push('\t');
        push_opt_f64(&mut buf, r.effect_size);
        buf.push('\t');
        buf.push_str(&r.effect_size_method);
        buf.push('\t');
        push_opt_f64(&mut buf, r.ci_low);
        buf.push('\t');
        push_opt_f64(&mut buf, r.ci_high);
        buf.push('\t');
        push_opt_f64(&mut buf, r.wilcoxon_p);
        buf.push('\t');
        buf.push_str(&r.wilcoxon_method);
        buf.push('\t');
        push_opt_f64(&mut buf, r.median_diff);
        buf.push('\t');
        push_opt_f64(&mut buf, r.trimmed_mean_diff);
        buf.push('\t');
        push_opt_f64(&mut buf, r.s2_trend);
        buf.push('\t');
        push_opt_f64(&mut buf, r.s2_prior);
        buf.push('\t');
        push_opt_f64(&mut buf, r.s2_posterior);
        buf.push('\t');
        push_opt_f64(&mut buf, r.df_prior);
        buf.push('\t');
        push_opt_f64(&mut buf, r.df_total);
        buf.push('\t');
        push_opt_f64(&mut buf, r.f_statistic);
        buf.push('\t');
        push_opt_f64(&mut buf, r.f_p_value);
        buf.push('\t');
        push_opt_f64(&mut buf, r.f_bh_q);
        buf.push('\t');
        push_opt_f64(&mut buf, r.lfc_threshold);
        buf.push('\t');
        if let Some(n) = r.n_peptides_observed {
            buf.push_str(&n.to_string());
        }
        buf.push('\t');
        push_opt_f64(&mut buf, r.peptide_variance_ratio);
        buf.push('\t');
        push_opt_f64(&mut buf, r.ridge_lambda);
        buf.push('\t');
        buf.push_str(&r.method);
        buf.push('\t');
        buf.push_str(&r.posthoc_method);
        buf.push('\t');
        push_opt_f64(&mut buf, r.posthoc_p);
        buf.push('\t');
        push_opt_f64(&mut buf, r.posthoc_adj_p);
        buf.push('\t');
        if let Some(n) = r.n_a {
            buf.push_str(&n.to_string());
        }
        buf.push('\t');
        if let Some(n) = r.n_b {
            buf.push_str(&n.to_string());
        }
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

/// One row of `de_report.tsv`. Per-panel, per-comparison summary.
///
/// `n_q_strict` and `n_q_relaxed` count BH-q values below the strict
/// and relaxed thresholds resolved at the CLI (`--report-q-strict` and
/// `--report-q-relaxed` on `atman de`); the resolved numeric values
/// are stamped in the run sidecar's `report_thresholds` block.
pub struct DeReportRow {
    pub comparison: String,
    pub panel: String,
    pub n_tests: usize,
    pub n_skipped: usize,
    pub n_q_strict: usize,
    pub n_q_relaxed: usize,
    pub min_q: Option<f64>,
    pub max_abs_effect: Option<f64>,
    // Appended by the limma eBayes + F-tests feature.
    pub limma_trend_fallback_used: Option<bool>,
}

pub fn write_de_report(path: &Path, rows: &[DeReportRow]) -> Result<()> {
    let mut buf = String::from(
        "comparison\tpanel\tn_tests\tn_skipped\tn_q_strict\tn_q_relaxed\tmin_q\tmax_abs_effect\t\
         limma_trend_fallback_used\n",
    );
    for r in rows {
        buf.push_str(&r.comparison);
        buf.push('\t');
        buf.push_str(&r.panel);
        buf.push('\t');
        buf.push_str(&r.n_tests.to_string());
        buf.push('\t');
        buf.push_str(&r.n_skipped.to_string());
        buf.push('\t');
        buf.push_str(&r.n_q_strict.to_string());
        buf.push('\t');
        buf.push_str(&r.n_q_relaxed.to_string());
        buf.push('\t');
        push_opt_f64(&mut buf, r.min_q);
        buf.push('\t');
        push_opt_f64(&mut buf, r.max_abs_effect);
        buf.push('\t');
        match r.limma_trend_fallback_used {
            Some(true) => buf.push('1'),
            Some(false) => buf.push('0'),
            None => {}
        }
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

/// One row of `de_ensemble.tsv` — per (comparison, protein)
/// cross-method summary. Produced only by `atman de --test ensemble`.
pub struct EnsembleRow {
    pub comparison: String,
    pub panel: String,
    pub assay_id: String,
    pub gene_symbol: String,
    pub uniprot: String,
    pub n_applied: usize,
    pub n_significant: usize,
    pub n_sign_consistent: usize,
    pub majority_sign: f64,
    pub ensemble_p: Option<f64>,
    pub ensemble_q: Option<f64>,
    pub grade: String,
    pub methods_applied: String,
    pub methods_skipped: String,
}

pub fn write_de_ensemble(path: &Path, rows: &[EnsembleRow]) -> Result<()> {
    let mut buf = String::from(
        "comparison\tpanel\tassay_id\tgene_symbol\tuniprot\t\
         n_applied\tn_significant\tn_sign_consistent\tmajority_sign\t\
         ensemble_p\tensemble_q\tgrade\tmethods_applied\tmethods_skipped\n",
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
        buf.push_str(&r.n_applied.to_string());
        buf.push('\t');
        buf.push_str(&r.n_significant.to_string());
        buf.push('\t');
        buf.push_str(&r.n_sign_consistent.to_string());
        buf.push('\t');
        buf.push_str(&format_float(r.majority_sign));
        buf.push('\t');
        push_opt_f64(&mut buf, r.ensemble_p);
        buf.push('\t');
        push_opt_f64(&mut buf, r.ensemble_q);
        buf.push('\t');
        buf.push_str(&r.grade);
        buf.push('\t');
        buf.push_str(&r.methods_applied);
        buf.push('\t');
        buf.push_str(&r.methods_skipped);
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn push_opt_f64(buf: &mut String, v: Option<f64>) {
    if let Some(x) = v {
        buf.push_str(&format!("{}", x));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn safe_filename_component_blocks_traversal_and_separators() {
        // Legitimate panel names pass (lowercased).
        assert_eq!(
            safe_filename_component("Cardiometabolic").unwrap(),
            "cardiometabolic"
        );
        assert_eq!(safe_filename_component("Panel 2").unwrap(), "panel 2");
        // Data-derived traversal / separators are rejected (strict-failure),
        // so a crafted `panel` column cannot escape the output directory.
        for bad in ["..", ".", "../evil", "a/b", "a\\b", "../../tmp/x", ""] {
            assert!(
                safe_filename_component(bad).is_err(),
                "expected rejection for {bad:?}"
            );
        }
    }

    #[test]
    fn atomic_write_creates_file_with_content() {
        let d = TempDir::new().unwrap();
        let p = d.path().join("out.tsv");
        atomic_write(&p, b"hello\n").unwrap();
        assert_eq!(std::fs::read_to_string(&p).unwrap(), "hello\n");
    }

    #[test]
    fn atomic_write_overwrites_existing() {
        let d = TempDir::new().unwrap();
        let p = d.path().join("out.tsv");
        atomic_write(&p, b"old").unwrap();
        atomic_write(&p, b"new").unwrap();
        assert_eq!(std::fs::read(&p).unwrap(), b"new");
    }

    #[test]
    fn parse_qc_accepts_canonical_spellings() {
        assert_eq!(parse_qc("PASS").unwrap(), QcFlag::Pass);
        assert_eq!(parse_qc("WARN").unwrap(), QcFlag::Warn(String::new()));
        assert_eq!(parse_qc("FAIL").unwrap(), QcFlag::Fail(String::new()));
        // Empty is the "no flag" form and means Pass (consistent with
        // validate_measurements, which also accepts "").
        assert_eq!(parse_qc("").unwrap(), QcFlag::Pass);
    }

    #[test]
    fn parse_qc_rejects_unrecognized_values() {
        // Lowercase, truncated, and garbage cells must error rather than silently
        // mapping to Pass — a garbled "fail" cell must not masquerade as passed.
        // (Empty is a valid "no flag" form and is tested as accepted above.)
        for bad in ["pass", "fail", "FA", "OK", "PASS "] {
            assert!(parse_qc(bad).is_err(), "expected error for {bad:?}");
        }
    }

    #[test]
    fn de_results_round_trip_preserves_limma_columns() {
        let d = TempDir::new().unwrap();
        let p = d.path().join("de_results.tsv");
        let row = DeResultRow {
            panel: "P1".into(),
            assay_id: "A001".into(),
            gene_symbol: "GENE1".into(),
            uniprot: "Q00001".into(),
            comparison: "A-B".into(),
            n_pairs: 5,
            mean_a: Some(1.0),
            mean_b: Some(0.0),
            mean_diff: Some(1.0),
            t: Some(2.5),
            df: Some(10.0),
            p_value: Some(0.01),
            bh_q: Some(0.05),
            effect_size: Some(1.2),
            effect_size_method: "limma-eBayes-robust-trend".into(),
            ci_low: Some(0.3),
            ci_high: Some(1.7),
            wilcoxon_p: None,
            wilcoxon_method: "".into(),
            median_diff: None,
            trimmed_mean_diff: None,
            skip_reason: "".into(),
            // New limma columns:
            s2_trend: Some(0.5),
            s2_prior: Some(0.4),
            s2_posterior: Some(0.45),
            df_prior: Some(8.0),
            df_total: Some(18.0),
            f_statistic: Some(6.25),
            f_p_value: Some(0.05),
            f_bh_q: Some(0.1),
            lfc_threshold: Some(0.0),
            n_peptides_observed: None,
            peptide_variance_ratio: None,
            ridge_lambda: None,
            method: "limma".into(),
            posthoc_method: String::new(),
            posthoc_p: None,
            posthoc_adj_p: None,
            n_a: None,
            n_b: None,
        };
        write_de_results(&p, &[row]).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let header = text.lines().next().unwrap();
        for col in [
            "s2_trend",
            "s2_prior",
            "s2_posterior",
            "df_prior",
            "df_total",
            "f_statistic",
            "f_p_value",
            "f_bh_q",
            "lfc_threshold",
        ] {
            assert!(header.contains(col), "header missing {col}: {header}");
        }
        let body = text.lines().nth(1).unwrap();
        assert!(body.contains("0.5")); // s2_trend
        assert!(body.contains("limma-eBayes-robust-trend"));
    }

    #[test]
    fn de_ensemble_tsv_roundtrips_all_columns() {
        let d = TempDir::new().unwrap();
        let p = d.path().join("de_ensemble.tsv");
        let row = EnsembleRow {
            comparison: "PT2-PR2".into(),
            panel: "Inflammation".into(),
            assay_id: "OID12345".into(),
            gene_symbol: "HSPA1A".into(),
            uniprot: "P08107".into(),
            n_applied: 4,
            n_significant: 4,
            n_sign_consistent: 4,
            majority_sign: 1.0,
            ensemble_p: Some(1e-9),
            ensemble_q: Some(1e-8),
            grade: "VALIDATED".into(),
            methods_applied: "welch-t,ols,limma,msqrob".into(),
            methods_skipped: "paired-t".into(),
        };
        write_de_ensemble(&p, &[row]).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let header = text.lines().next().unwrap();
        for col in [
            "comparison",
            "panel",
            "assay_id",
            "gene_symbol",
            "uniprot",
            "n_applied",
            "n_significant",
            "n_sign_consistent",
            "majority_sign",
            "ensemble_p",
            "ensemble_q",
            "grade",
            "methods_applied",
            "methods_skipped",
        ] {
            assert!(header.contains(col), "missing column {col}: {header}");
        }
        let body = text.lines().nth(1).unwrap();
        assert!(body.contains("VALIDATED"));
        assert!(body.contains("HSPA1A"));
    }

    #[test]
    fn de_results_method_column_roundtrips() {
        let d = TempDir::new().unwrap();
        let p = d.path().join("de_results.tsv");
        let row = DeResultRow {
            panel: "P1".into(),
            assay_id: "A001".into(),
            gene_symbol: "GENE1".into(),
            uniprot: "Q00001".into(),
            comparison: "A-B".into(),
            n_pairs: 5,
            mean_a: Some(1.0),
            mean_b: Some(0.0),
            mean_diff: Some(1.0),
            t: Some(2.5),
            df: Some(10.0),
            p_value: Some(0.01),
            bh_q: Some(0.05),
            effect_size: None,
            effect_size_method: "welch-t".into(),
            ci_low: None,
            ci_high: None,
            wilcoxon_p: None,
            wilcoxon_method: "".into(),
            median_diff: None,
            trimmed_mean_diff: None,
            skip_reason: "".into(),
            s2_trend: None,
            s2_prior: None,
            s2_posterior: None,
            df_prior: None,
            df_total: None,
            f_statistic: None,
            f_p_value: None,
            f_bh_q: None,
            lfc_threshold: None,
            n_peptides_observed: None,
            peptide_variance_ratio: None,
            ridge_lambda: None,
            method: "welch-t".into(),
            posthoc_method: String::new(),
            posthoc_p: None,
            posthoc_adj_p: None,
            n_a: None,
            n_b: None,
        };
        write_de_results(&p, &[row]).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let header = text.lines().next().unwrap();
        // `method` is no longer the final column after DEBT-1
        // introduced posthoc_method/posthoc_p/posthoc_adj_p columns
        // behind it. Assert the column exists instead.
        assert!(
            header.split('\t').any(|c| c == "method"),
            "method column missing from header: {header}"
        );
        let body = text.lines().nth(1).unwrap();
        assert!(
            body.split('\t').any(|c| c == "welch-t"),
            "method value missing from body: {body}"
        );
    }
}
