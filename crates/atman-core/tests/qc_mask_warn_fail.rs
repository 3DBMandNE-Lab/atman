use atman_core::qc::apply_mask_warn_fail;
use atman_core::{Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, QcFlag};

fn rec(qc_s: QcFlag, qc_a: QcFlag) -> MeasurementRecord {
    MeasurementRecord {
        platform: Platform::OlinkExploreNgs,
        assay_id: AssayId("OID1".into()),
        gene_symbol: Some("X".into()),
        sample_id: "SSNA-001B-PR1".into(),
        abundance: Abundance::Log2Npx(0.5),
        abundance_raw: Abundance::Log2Npx(0.5),
        npx_source_str: "0.5".into(),
        qc_sample: qc_s,
        qc_assay: qc_a,
        detection_limit: DetectionLimit(Some(0.2)),
        below_lod: false,
        batch: Batch {
            plate: None,
            lot: None,
            run: None,
        },
        dropped_by_qc: false,
        ingest_order: 0,
        panel: Some("Neurology".into()),
    }
}

#[test]
fn pass_pass_unchanged() {
    let mut r = rec(QcFlag::Pass, QcFlag::Pass);
    apply_mask_warn_fail(&mut r);
    assert!(!r.dropped_by_qc);
    assert!(matches!(r.abundance, Abundance::Log2Npx(v) if (v - 0.5).abs() < 1e-12));
}

#[test]
fn warn_sample_masks_abundance() {
    let mut r = rec(QcFlag::Warn("".into()), QcFlag::Pass);
    apply_mask_warn_fail(&mut r);
    assert!(r.dropped_by_qc);
    // abundance_raw must remain intact
    assert!(matches!(r.abundance_raw, Abundance::Log2Npx(v) if (v - 0.5).abs() < 1e-12));
}

#[test]
fn warn_assay_masks_abundance() {
    let mut r = rec(QcFlag::Pass, QcFlag::Warn("".into()));
    apply_mask_warn_fail(&mut r);
    assert!(r.dropped_by_qc);
}

#[test]
fn warn_both_masks() {
    let mut r = rec(QcFlag::Warn("".into()), QcFlag::Warn("".into()));
    apply_mask_warn_fail(&mut r);
    assert!(r.dropped_by_qc);
}

#[test]
fn idempotent() {
    let mut r = rec(QcFlag::Warn("".into()), QcFlag::Pass);
    apply_mask_warn_fail(&mut r);
    let snapshot = r.clone();
    apply_mask_warn_fail(&mut r);
    assert_eq!(r, snapshot);
}
