//! File I/O helpers: atomic writes, long-TSV read/write, CSV write.
//!
//! The long TSV is header-indexed on read (not positional) so adding columns
//! in the future does not silently mis-read older files. The writer emits a
//! stable column order but the reader tolerates either ordering.

use anyhow::{anyhow, Context, Result};
use proteome_core::{
    fold_change::FoldChangePanel,
    matrix::DubeWidePanel,
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, ProteinIdentity,
    QcFlag, Sample,
};
use std::{
    collections::HashMap,
    io::Write,
    path::{Path, PathBuf},
};

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

const LONG_TSV_HEADER: &[&str] = &[
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
        buf.push_str(&r.sample_id);             buf.push('\t');
        buf.push_str(&r.assay_id.0);            buf.push('\t');
        buf.push_str(&gene);                    buf.push('\t');
        buf.push_str(&panel);                   buf.push('\t');
        buf.push_str(&r.npx_source_str);        buf.push('\t');
        buf.push_str(&abundance);               buf.push('\t');
        buf.push_str(&abundance_raw);           buf.push('\t');
        buf.push_str("log2_npx");               buf.push('\t');
        buf.push_str(r.qc_sample.as_str());     buf.push('\t');
        buf.push_str(r.qc_assay.as_str());      buf.push('\t');
        buf.push_str(&lod);                     buf.push('\t');
        buf.push_str(&(r.below_lod as u8).to_string());          buf.push('\t');
        buf.push_str(&(r.dropped_by_qc as u8).to_string());      buf.push('\t');
        buf.push_str(&plate);                   buf.push('\t');
        buf.push_str(&lot);                     buf.push('\t');
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
    let c_assay  = need("assay_id")?;
    let c_gene   = need("gene_symbol")?;
    let c_panel  = need("panel")?;
    let c_src    = need("npx_source_str")?;
    let c_abund  = need("abundance")?;
    let c_abund_raw = need("abundance_raw")?;
    let c_qcs    = need("qc_sample")?;
    let c_qca    = need("qc_assay")?;
    let c_lod    = need("detection_limit")?;
    let c_below  = need("below_lod")?;
    let c_drop   = need("dropped_by_qc")?;
    let c_plate  = need("plate_id")?;
    let c_lot    = need("panel_lot")?;
    let c_order  = need("ingest_order")?;

    let mut out = Vec::new();
    for result in reader.records() {
        let row = result.with_context(|| format!("reading record from {:?}", path))?;
        let abundance_raw: f64 = row[c_abund_raw]
            .parse()
            .context("abundance_raw parse")?;
        let dropped_by_qc: bool = row[c_drop].parse::<u8>().map(|v| v != 0).unwrap_or(false);
        let abund_str = &row[c_abund];
        let abundance = if abund_str.is_empty() {
            Abundance::Log2Npx(abundance_raw)
        } else {
            Abundance::Log2Npx(abund_str.parse().context("abundance parse")?)
        };
        let detection_limit = if row[c_lod].is_empty() {
            DetectionLimit(None)
        } else {
            DetectionLimit(Some(row[c_lod].parse().context("detection_limit parse")?))
        };
        let below_lod: bool = row[c_below].parse::<u8>().map(|v| v != 0).unwrap_or(false);
        out.push(MeasurementRecord {
            platform: Platform::OlinkExploreNgs,
            assay_id: AssayId(row[c_assay].to_string()),
            gene_symbol: empty_to_none(&row[c_gene]),
            sample_id: row[c_sample].to_string(),
            abundance,
            abundance_raw: Abundance::Log2Npx(abundance_raw),
            npx_source_str: row[c_src].to_string(),
            qc_sample: parse_qc(&row[c_qcs]),
            qc_assay: parse_qc(&row[c_qca]),
            detection_limit,
            below_lod,
            batch: Batch {
                plate: empty_to_none(&row[c_plate]),
                lot: empty_to_none(&row[c_lot]),
                run: None,
            },
            dropped_by_qc,
            ingest_order: row[c_order].parse().unwrap_or(0),
            panel: empty_to_none(&row[c_panel]),
        });
    }
    Ok(out)
}

/// samples.tsv writer.
pub fn write_samples(path: &Path, samples: &[Sample]) -> Result<()> {
    let mut buf = String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
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
        .from_path(path)?;
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        out.push(Sample {
            sample_id: row[0].to_string(),
            subject_id: empty_to_none(&row[1]),
            condition: empty_to_none(&row[2]),
            is_control: row[3].parse::<u8>()? != 0,
            sample_type: empty_to_none(&row[4]),
            ingest_order: row[5].parse().unwrap_or(0),
        });
    }
    Ok(out)
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

/// Dube-wide panel CSV writer. Writes source strings verbatim — no f64
/// round-trip — so reproduction against Dube's published files is byte-exact.
pub fn write_dube_wide_panel(
    out_dir: &Path,
    panel: &DubeWidePanel,
) -> Result<PathBuf> {
    let filename = format!("{}_npx.csv", panel.panel.to_ascii_lowercase());
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
                buf.push_str(s);
            }
        }
        buf.push('\n');
    }
    atomic_write(&path, buf.as_bytes())?;
    Ok(path)
}

/// Fold-change panel CSV writer.
pub fn write_fold_change_panel(
    out_dir: &Path,
    panel: &FoldChangePanel,
    comparisons: &[proteome_core::fold_change::Comparison],
) -> Result<PathBuf> {
    let filename = format!("{}_log2_fc.csv", panel.panel.to_ascii_lowercase());
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
                buf.push_str(&format_f64_fc(*x));
            }
        }
        buf.push('\n');
    }
    atomic_write(&path, buf.as_bytes())?;
    Ok(path)
}

fn parse_qc(s: &str) -> QcFlag {
    match s {
        "PASS" => QcFlag::Pass,
        "WARN" => QcFlag::Warn(String::new()),
        "FAIL" => QcFlag::Fail(String::new()),
        _ => QcFlag::Pass,
    }
}

fn empty_to_none(s: &str) -> Option<String> {
    if s.is_empty() { None } else { Some(s.to_string()) }
}

fn format_f64(v: f64) -> String {
    format!("{}", v)
}

fn format_f64_fc(v: f64) -> String {
    format!("{}", v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

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
}
