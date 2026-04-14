//! Olink Explore NGS long-CSV ingest adapter. See spec §5.1.

use crate::{
    errors::IngestError,
    ingest::{IngestOutput, ProteomeIngest},
    sample_id::{ParsedSampleId, SampleIdParser},
    types::{
        Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, ProteinIdentity,
        QcFlag, Sample,
    },
};
use serde::Deserialize;
use std::{
    collections::{HashMap, HashSet},
    path::{Path, PathBuf},
};

pub struct OlinkExploreLongCsv;

const EXPECTED_HEADERS: &[&str] = &[
    "SampleID",
    "Index",
    "OlinkID",
    "UniProt",
    "Assay",
    "MissingFreq",
    "Panel",
    "Panel_Lot_Nr",
    "PlateID",
    "QC_Warning",
    "LOD",
    "NPX",
    "Normalization",
    "Assay_Warning",
];

#[derive(Debug, Deserialize)]
struct RawRow {
    #[serde(rename = "SampleID")]
    sample_id: String,
    #[serde(rename = "Index")]
    _index: String,
    #[serde(rename = "OlinkID")]
    olink_id: String,
    #[serde(rename = "UniProt")]
    uniprot: String,
    #[serde(rename = "Assay")]
    assay: String,
    #[serde(rename = "MissingFreq")]
    _missing_freq: String,
    #[serde(rename = "Panel")]
    panel: String,
    #[serde(rename = "Panel_Lot_Nr")]
    panel_lot: String,
    #[serde(rename = "PlateID")]
    plate_id: String,
    #[serde(rename = "QC_Warning")]
    qc_warning: String,
    #[serde(rename = "LOD")]
    lod: String,
    #[serde(rename = "NPX")]
    npx: String,
    #[serde(rename = "Normalization")]
    normalization: String,
    #[serde(rename = "Assay_Warning")]
    assay_warning: String,
}

impl ProteomeIngest for OlinkExploreLongCsv {
    fn platform(&self) -> Platform {
        Platform::OlinkExploreNgs
    }

    fn read(
        &self,
        inputs: &[PathBuf],
        sample_id_parser: &dyn SampleIdParser,
    ) -> Result<IngestOutput, IngestError> {
        let mut out = IngestOutput::default();
        let mut seen_pk: HashSet<(String, String)> = HashSet::new();
        let mut seen_sample: HashMap<String, u64> = HashMap::new();
        let mut protein_catalog: HashMap<String, ProteinIdentity> = HashMap::new();
        let mut ingest_counter: u64 = 0;

        for input in inputs {
            read_one_file(
                input,
                sample_id_parser,
                &mut out,
                &mut seen_pk,
                &mut seen_sample,
                &mut protein_catalog,
                &mut ingest_counter,
            )?;
        }

        // Finalize protein catalog — stable sort by assay_id.
        let mut proteins: Vec<ProteinIdentity> = protein_catalog.into_values().collect();
        proteins.sort_by(|a, b| a.assay_id.0.cmp(&b.assay_id.0));
        out.proteins = proteins;

        Ok(out)
    }
}

fn read_one_file(
    path: &Path,
    sample_id_parser: &dyn SampleIdParser,
    out: &mut IngestOutput,
    seen_pk: &mut HashSet<(String, String)>,
    seen_sample: &mut HashMap<String, u64>,
    protein_catalog: &mut HashMap<String, ProteinIdentity>,
    ingest_counter: &mut u64,
) -> Result<(), IngestError> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b';')
        .has_headers(true)
        .from_path(path)
        .map_err(|e| IngestError::CsvError {
            file: path.to_path_buf(),
            source: e,
        })?;

    let headers: Vec<String> = reader
        .headers()
        .map_err(|e| IngestError::CsvError {
            file: path.to_path_buf(),
            source: e,
        })?
        .iter()
        .map(|s| s.to_string())
        .collect();
    validate_headers(path, &headers)?;

    for (line_idx, result) in reader.deserialize::<RawRow>().enumerate() {
        let raw: RawRow = result.map_err(|e| IngestError::CsvError {
            file: path.to_path_buf(),
            source: e,
        })?;
        let line_number = line_idx + 2;

        if raw.normalization != "Plate control" {
            return Err(IngestError::UnexpectedNormalization {
                file: path.to_path_buf(),
                line: line_number,
                value: raw.normalization,
            });
        }

        let npx: f64 = raw.npx.parse().map_err(|_| IngestError::InvalidAbundance {
            file: path.to_path_buf(),
            line: line_number,
            raw: raw.npx.clone(),
        })?;
        let lod: f64 = raw.lod.parse().map_err(|_| IngestError::InvalidAbundance {
            file: path.to_path_buf(),
            line: line_number,
            raw: raw.lod.clone(),
        })?;

        protein_catalog
            .entry(raw.olink_id.clone())
            .or_insert_with(|| ProteinIdentity {
                platform: Platform::OlinkExploreNgs,
                assay_id: AssayId(raw.olink_id.clone()),
                uniprot: vec![raw.uniprot.clone()],
                gene_symbol: Some(raw.assay.clone()),
                panel: Some(raw.panel.clone()),
                panel_lot: Some(raw.panel_lot.clone()),
            });

        let pk = (raw.olink_id.clone(), raw.sample_id.clone());
        if !seen_pk.insert(pk) {
            return Err(IngestError::DuplicatePrimaryKey {
                platform: Platform::OlinkExploreNgs.as_str().to_string(),
                assay_id: raw.olink_id.clone(),
                sample_id: raw.sample_id.clone(),
            });
        }

        if !seen_sample.contains_key(&raw.sample_id) {
            let parsed = sample_id_parser.parse(&raw.sample_id)?;
            let (subject_id, condition, is_control) = match parsed {
                ParsedSampleId::Biological { subject, condition } => {
                    (Some(subject), Some(condition), false)
                }
                ParsedSampleId::Control => (None, None, true),
            };
            let order = *ingest_counter;
            seen_sample.insert(raw.sample_id.clone(), order);
            *ingest_counter += 1;
            out.samples.push(Sample {
                sample_id: raw.sample_id.clone(),
                subject_id,
                condition,
                is_control,
                sample_type: None,
                ingest_order: order,
            });
        }

        let qc_sample = parse_qc_flag(&raw.qc_warning);
        let qc_assay = parse_qc_flag(&raw.assay_warning);

        let npx_source_str = raw.npx.clone();
        let gene_symbol = Some(raw.assay.clone());
        let panel = Some(raw.panel.clone());

        out.measurements.push(MeasurementRecord {
            platform: Platform::OlinkExploreNgs,
            assay_id: AssayId(raw.olink_id),
            gene_symbol,
            sample_id: raw.sample_id,
            abundance: Abundance::Log2Npx(npx),
            abundance_raw: Abundance::Log2Npx(npx),
            npx_source_str,
            qc_sample,
            qc_assay,
            detection_limit: DetectionLimit(Some(lod)),
            below_lod: npx < lod,
            batch: Batch {
                plate: Some(raw.plate_id),
                lot: Some(raw.panel_lot),
                run: None,
            },
            dropped_by_qc: false,
            ingest_order: *ingest_counter,
            panel,
        });
        *ingest_counter += 1;
    }

    Ok(())
}

fn validate_headers(path: &Path, headers: &[String]) -> Result<(), IngestError> {
    let got: HashSet<&str> = headers.iter().map(|s| s.as_str()).collect();
    let want: HashSet<&str> = EXPECTED_HEADERS.iter().copied().collect();
    let missing: Vec<String> = want.difference(&got).map(|s| (*s).to_string()).collect();
    let extra: Vec<String> = got.difference(&want).map(|s| (*s).to_string()).collect();
    if !missing.is_empty() || !extra.is_empty() {
        return Err(IngestError::SchemaMismatch {
            file: path.to_path_buf(),
            missing,
            extra,
        });
    }
    Ok(())
}

fn parse_qc_flag(s: &str) -> QcFlag {
    match s {
        "PASS" => QcFlag::Pass,
        "WARN" => QcFlag::Warn(String::new()),
        "FAIL" => QcFlag::Fail(String::new()),
        other => QcFlag::Warn(other.to_string()),
    }
}
