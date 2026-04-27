//! QC rules. The `mask-warn-fail` rule marks `dropped_by_qc=true`
//! when either `qc_sample` or `qc_assay` is not `Pass`. `abundance` and
//! `abundance_raw` are never mutated in memory; downstream code uses
//! `MeasurementRecord::effective_abundance()` which honors the flag.

use crate::types::MeasurementRecord;

/// Apply the `mask-warn-fail` rule to a single record in-place. Idempotent.
pub fn apply_mask_warn_fail(record: &mut MeasurementRecord) {
    let masked = !record.qc_sample.is_pass() || !record.qc_assay.is_pass();
    if masked {
        record.dropped_by_qc = true;
    }
}

/// Convenience: apply rule to all records in a slice.
pub fn apply_mask_warn_fail_all(records: &mut [MeasurementRecord]) {
    for r in records.iter_mut() {
        apply_mask_warn_fail(r);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, QcFlag};

    #[test]
    fn rule_sets_dropped_by_qc_on_warn() {
        let mut r = MeasurementRecord {
            platform: Platform::OlinkExploreNgs,
            assay_id: AssayId("X".into()),
            gene_symbol: Some("X".into()),
            sample_id: "S".into(),
            abundance: Abundance::Log2Npx(1.0),
            abundance_raw: Abundance::Log2Npx(1.0),
            npx_source_str: "1.0".into(),
            qc_sample: QcFlag::Warn("".into()),
            qc_assay: QcFlag::Pass,
            detection_limit: DetectionLimit(None),
            below_lod: false,
            batch: Batch {
                plate: None,
                lot: None,
                run: None,
            },
            dropped_by_qc: false,
            ingest_order: 0,
            panel: None,
        };
        apply_mask_warn_fail(&mut r);
        assert!(r.dropped_by_qc);
    }
}
