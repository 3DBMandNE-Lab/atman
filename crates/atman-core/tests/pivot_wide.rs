use atman_core::matrix::{pivot_wide_panels, WidePanel};
use atman_core::{
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, QcFlag, Sample,
};

#[allow(clippy::too_many_arguments)]
fn m(
    sample: &str,
    assay_id: &str,
    gene: &str,
    panel: &str,
    source: &str,
    val: f64,
    dropped: bool,
    order: u64,
) -> MeasurementRecord {
    MeasurementRecord {
        platform: Platform::OlinkExploreNgs,
        assay_id: AssayId(assay_id.into()),
        gene_symbol: Some(gene.into()),
        sample_id: sample.into(),
        abundance: Abundance::Log2Npx(val),
        abundance_raw: Abundance::Log2Npx(val),
        npx_source_str: source.into(),
        qc_sample: QcFlag::Pass,
        qc_assay: QcFlag::Pass,
        detection_limit: DetectionLimit(None),
        below_lod: false,
        batch: Batch {
            plate: None,
            lot: None,
            run: None,
        },
        dropped_by_qc: dropped,
        ingest_order: order,
        panel: Some(panel.into()),
    }
}

fn sample_bio(id: &str, subj: &str, cond: &str, order: u64) -> Sample {
    Sample {
        sample_id: id.into(),
        subject_id: Some(subj.into()),
        condition: Some(cond.into()),
        is_control: false,
        sample_type: None,
        ingest_order: order,
    }
}

fn sample_ctl(id: &str, order: u64) -> Sample {
    Sample {
        sample_id: id.into(),
        subject_id: None,
        condition: None,
        is_control: true,
        sample_type: None,
        ingest_order: order,
    }
}

#[test]
fn pivot_two_panels_bio_then_control_by_gene_symbol() {
    let measurements = vec![
        m("SSNA-001B-PR1", "OID2", "ZZZ", "P1", "1.0", 1.0, false, 0),
        m("SSNA-001B-PR1", "OID1", "AAA", "P1", "2.0", 2.0, false, 1),
        m("SSNA-001B-PR2", "OID1", "AAA", "P1", "3.0", 3.0, false, 2),
        m("SSNA-001B-PR2", "OID2", "ZZZ", "P1", "4.0", 4.0, false, 3),
        m(
            "CONTROL_SAMPLE_X",
            "OID1",
            "AAA",
            "P1",
            "9.0",
            9.0,
            false,
            4,
        ),
        m(
            "CONTROL_SAMPLE_X",
            "OID2",
            "ZZZ",
            "P1",
            "8.0",
            8.0,
            false,
            5,
        ),
        m("SSNA-001B-PR1", "OID3", "BBB", "P2", "5.0", 5.0, false, 6),
    ];
    let samples = vec![
        sample_bio("SSNA-001B-PR1", "001B", "PR1", 0),
        sample_bio("SSNA-001B-PR2", "001B", "PR2", 1),
        sample_ctl("CONTROL_SAMPLE_X", 2),
    ];

    let panels = pivot_wide_panels(&measurements, &samples);
    assert_eq!(panels.len(), 2);

    let p1: &WidePanel = panels.iter().find(|p| p.panel == "P1").unwrap();
    assert_eq!(p1.assays, vec!["AAA".to_string(), "ZZZ".to_string()]);
    assert_eq!(p1.rows.len(), 3);
    assert_eq!(p1.rows[0].participant, "001B");
    assert_eq!(p1.rows[0].exposure, "PR1");
    assert_eq!(p1.rows[0].sample_id, "SSNA-001B-PR1");
    assert_eq!(p1.rows[0].values[0], Some("2.0".to_string())); // AAA
    assert_eq!(p1.rows[0].values[1], Some("1.0".to_string())); // ZZZ

    assert_eq!(p1.rows[1].sample_id, "SSNA-001B-PR2");
    assert_eq!(p1.rows[2].participant, "");
    assert_eq!(p1.rows[2].exposure, "");
    assert_eq!(p1.rows[2].sample_id, "CONTROL_SAMPLE_X");
}

#[test]
fn dropped_by_qc_becomes_missing() {
    let m1 = m("SSNA-001B-PR1", "OID1", "AAA", "P1", "1.0", 1.0, true, 0);
    let m2 = m("SSNA-001B-PR1", "OID2", "BBB", "P1", "2.0", 2.0, false, 1);
    let s = vec![sample_bio("SSNA-001B-PR1", "001B", "PR1", 0)];
    let panels = pivot_wide_panels(&[m1, m2], &s);
    assert_eq!(panels.len(), 1);
    let row = &panels[0].rows[0];
    assert_eq!(row.values[0], None); // AAA masked
    assert_eq!(row.values[1], Some("2.0".to_string())); // BBB kept
}

#[test]
fn control_rows_preserve_ingest_order() {
    let measurements = vec![
        m("CONTROL_B", "OID1", "AAA", "P1", "1.0", 1.0, false, 0),
        m("CONTROL_A", "OID1", "AAA", "P1", "2.0", 2.0, false, 1),
    ];
    let samples = vec![sample_ctl("CONTROL_B", 0), sample_ctl("CONTROL_A", 1)];
    let panels = pivot_wide_panels(&measurements, &samples);
    assert_eq!(panels[0].rows[0].sample_id, "CONTROL_B"); // ingest order, not alpha
    assert_eq!(panels[0].rows[1].sample_id, "CONTROL_A");
}
