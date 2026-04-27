//! Long → wide-per-panel pivot.
//!
//! Produces one `WidePanel` per raw `Panel` value, keyed by gene symbol
//! (the `gene_symbol` column on `MeasurementRecord`), with string cell
//! values copied verbatim from `MeasurementRecord.npx_source_str`. This
//! is a display pivot: no numeric arithmetic, no reformatting, so the
//! output reproduces source byte-for-byte and is suitable for diffing
//! against published reference NPX tables.

use crate::types::{MeasurementRecord, Sample};
use std::collections::{BTreeMap, BTreeSet, HashMap};

#[derive(Debug, Clone, PartialEq)]
pub struct WidePanelRow {
    /// Empty string for control samples; subject_id for biological samples.
    pub participant: String,
    /// Empty string for control samples; condition for biological samples.
    pub exposure: String,
    pub sample_id: String,
    /// Aligned with `WidePanel::assays`. `None` = missing / masked.
    pub values: Vec<Option<String>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct WidePanel {
    /// Raw panel name from the source data, e.g. `Inflammation_II`.
    pub panel: String,
    /// Gene symbols, sorted by Unicode codepoint order.
    pub assays: Vec<String>,
    pub rows: Vec<WidePanelRow>,
}

pub fn pivot_wide_panels(measurements: &[MeasurementRecord], samples: &[Sample]) -> Vec<WidePanel> {
    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    let mut by_panel: BTreeMap<String, PanelBuilder> = BTreeMap::new();

    for m in measurements {
        let panel = match &m.panel {
            Some(p) => p.clone(),
            None => continue,
        };
        let gene = match &m.gene_symbol {
            Some(g) => g.clone(),
            None => continue,
        };
        let entry = by_panel.entry(panel).or_default();
        entry.note_assay(gene.clone());
        entry.set_value(&m.sample_id, &gene, m);
    }

    by_panel
        .into_iter()
        .map(|(panel_name, builder)| builder.finalize(panel_name, &sample_by_id))
        .collect()
}

#[derive(Default)]
struct PanelBuilder {
    assays: BTreeSet<String>,
    /// sample_id → gene_symbol → Option<npx_source_str>. Each entry is a
    /// sample-level measurement; the biological replicate (subject_id) is
    /// the statistical unit and lives on the `Sample` struct, never here.
    measurements: HashMap<String, HashMap<String, Option<String>>>,
}

impl PanelBuilder {
    fn note_assay(&mut self, gene: String) {
        self.assays.insert(gene);
    }

    fn set_value(&mut self, sample_id: &str, gene: &str, m: &MeasurementRecord) {
        let value = if m.dropped_by_qc {
            None
        } else {
            Some(m.npx_source_str.clone())
        };
        self.measurements
            .entry(sample_id.to_string())
            .or_default()
            .insert(gene.to_string(), value);
    }

    fn finalize(self, panel: String, sample_by_id: &HashMap<&str, &Sample>) -> WidePanel {
        let assays: Vec<String> = self.assays.into_iter().collect();

        let mut bio: Vec<WidePanelRow> = Vec::new();
        let mut ctl: Vec<(u64, WidePanelRow)> = Vec::new();

        for (sample_id, assay_values) in self.measurements {
            let s = match sample_by_id.get(sample_id.as_str()) {
                Some(s) => *s,
                None => continue,
            };
            let values: Vec<Option<String>> = assays
                .iter()
                .map(|gene| assay_values.get(gene).cloned().unwrap_or(None))
                .collect();

            if s.is_control {
                let row = WidePanelRow {
                    participant: String::new(),
                    exposure: String::new(),
                    sample_id: s.sample_id.clone(),
                    values,
                };
                ctl.push((s.ingest_order, row));
            } else {
                let row = WidePanelRow {
                    participant: s.subject_id.clone().unwrap_or_default(),
                    exposure: s.condition.clone().unwrap_or_default(),
                    sample_id: s.sample_id.clone(),
                    values,
                };
                bio.push(row);
            }
        }

        bio.sort_by(|a, b| {
            a.participant
                .cmp(&b.participant)
                .then_with(|| a.exposure.cmp(&b.exposure))
                .then_with(|| a.sample_id.cmp(&b.sample_id))
        });
        ctl.sort_by_key(|(order, _)| *order);

        let mut rows = bio;
        rows.extend(ctl.into_iter().map(|(_, r)| r));

        WidePanel {
            panel,
            assays,
            rows,
        }
    }
}
