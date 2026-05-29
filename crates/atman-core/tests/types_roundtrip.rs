use atman_core::{
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, ProteinIdentity,
    QcFlag, Sample,
};

fn make_record() -> MeasurementRecord {
    MeasurementRecord {
        platform: Platform::OlinkExploreNgs,
        assay_id: AssayId("OID20838".to_string()),
        gene_symbol: Some("GGA1".into()),
        sample_id: "SSNA-001B-PR1".to_string(),
        abundance: Abundance::Log2Npx(0.1234),
        abundance_raw: Abundance::Log2Npx(0.1234),
        npx_source_str: "0.1234".to_string(),
        qc_sample: QcFlag::Pass,
        qc_assay: QcFlag::Pass,
        detection_limit: DetectionLimit(Some(0.2178)),
        below_lod: true,
        batch: Batch {
            plate: Some("plate1".into()),
            lot: Some("B04414".into()),
            run: None,
        },
        dropped_by_qc: false,
        ingest_order: 0,
        panel: Some("Neurology".into()),
    }
}

#[test]
fn construct_measurement_record() {
    let rec = make_record();
    assert_eq!(rec.assay_id.0, "OID20838");
    assert_eq!(rec.gene_symbol.as_deref(), Some("GGA1"));
    assert_eq!(rec.npx_source_str, "0.1234");
    assert!(matches!(rec.abundance, Abundance::Log2Npx(v) if (v - 0.1234).abs() < 1e-12));
    assert!(matches!(rec.qc_sample, QcFlag::Pass));
}

#[test]
fn effective_abundance_accessor_honors_qc() {
    let mut rec = make_record();
    assert_eq!(rec.effective_abundance(), Some(0.1234));
    rec.dropped_by_qc = true;
    assert_eq!(rec.effective_abundance(), None);
}

#[test]
fn effective_abundance_rejects_non_finite_on_non_dropped_row() {
    // `f64::parse` accepts the literal strings "NaN"/"inf", so a non-dropped
    // row can carry a non-finite abundance. The accessor must return `None`
    // rather than leaking NaN/inf into downstream models, matching
    // `PeptideMeasurementRecord::effective_abundance`.
    let mut rec = make_record();
    assert!(!rec.dropped_by_qc);

    rec.abundance = Abundance::Log2Npx(f64::NAN);
    assert_eq!(rec.effective_abundance(), None);

    rec.abundance = Abundance::Log2Npx(f64::INFINITY);
    assert_eq!(rec.effective_abundance(), None);

    rec.abundance = Abundance::Log2Npx(f64::NEG_INFINITY);
    assert_eq!(rec.effective_abundance(), None);

    // A finite value on the same non-dropped row still passes through.
    rec.abundance = Abundance::Log2Npx(2.5);
    assert_eq!(rec.effective_abundance(), Some(2.5));
}

#[test]
fn qc_flag_warn_carries_reason() {
    let f = QcFlag::Warn("plate deviation".to_string());
    match f {
        QcFlag::Warn(r) => assert_eq!(r, "plate deviation"),
        _ => panic!("expected Warn"),
    }
}

#[test]
fn protein_identity_allows_multi_uniprot() {
    let p = ProteinIdentity {
        platform: Platform::OlinkExploreNgs,
        assay_id: AssayId("OID30000".to_string()),
        uniprot: vec!["P01375".into(), "P01376".into()],
        gene_symbol: Some("MICB_MICA".into()),
        panel: Some("Inflammation".into()),
        panel_lot: Some("B04414".into()),
    };
    assert_eq!(p.uniprot.len(), 2);
    assert_eq!(p.gene_symbol.as_deref(), Some("MICB_MICA"));
}

#[test]
fn sample_struct_supports_control() {
    let s = Sample {
        sample_id: "CONTROL_SAMPLE_US_CS_AS_2-1".into(),
        subject_id: None,
        condition: None,
        is_control: true,
        sample_type: Some("plasma".into()),
        ingest_order: 40,
    };
    assert!(s.is_control);
}
