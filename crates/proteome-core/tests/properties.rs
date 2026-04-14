use proteome_core::fold_change::{compute_log2_fc, Comparison, FoldChangeInput};
use proteome_core::matrix::dube_wide_pivot;
use proteome_core::{
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, QcFlag, Sample,
};
use proptest::prelude::*;

fn arb_value() -> impl Strategy<Value = f64> {
    (-5.0_f64..5.0_f64).prop_filter("finite", |v| v.is_finite())
}

// Property: fold_change(x, x) == 0 when both sides draw from the same set of
// non-missing participants (self-vs-self with the same participant universe).
proptest! {
    #[test]
    fn self_vs_self_is_zero_when_symmetric(values in prop::collection::vec(arb_value(), 1..6)) {
        // Build owned participant names so the input outlives from_cells.
        let participants: Vec<String> = (0..values.len())
            .map(|i| format!("P{}", i))
            .collect();
        let cells: Vec<(String, String, String, String, Option<f64>)> = participants
            .iter()
            .zip(values.iter())
            .map(|(p, v)| ("P".to_string(), "A".to_string(), p.clone(), "PR1".to_string(), Some(*v)))
            .collect();
        let input = FoldChangeInput::from_cells(cells);
        let comps = vec![Comparison { a: "PR1".into(), b: "PR1".into() }];
        let out = compute_log2_fc(&input, &comps);
        let fc = out.panels[0].values[0][0].unwrap();
        prop_assert!(fc.abs() < 1e-12, "self vs self was {}", fc);
    }
}

// Property: the pivot is deterministic regardless of input row order.
// With unique (participant, exposure) per biological sample and a single
// ingest-order for each control, reversing the measurement input must
// produce the same output.
#[test]
fn pivot_is_deterministic_across_input_order() {
    let s = vec![
        Sample {
            sample_id: "SSNA-1-PR1".into(),
            subject_id: Some("1".into()),
            condition: Some("PR1".into()),
            is_control: false,
            sample_type: None,
            ingest_order: 0,
        },
        Sample {
            sample_id: "SSNA-2-PR1".into(),
            subject_id: Some("2".into()),
            condition: Some("PR1".into()),
            is_control: false,
            sample_type: None,
            ingest_order: 1,
        },
    ];
    let base = vec![
        mk("SSNA-1-PR1", "Z", "P", 1.0, 0),
        mk("SSNA-1-PR1", "A", "P", 2.0, 1),
        mk("SSNA-2-PR1", "A", "P", 3.0, 2),
        mk("SSNA-2-PR1", "Z", "P", 4.0, 3),
    ];
    let mut reversed = base.clone();
    reversed.reverse();
    let a = dube_wide_pivot(&base, &s);
    let b = dube_wide_pivot(&reversed, &s);
    assert_eq!(a, b, "pivot output must not depend on input ordering");
}

fn mk(sample: &str, assay: &str, panel: &str, val: f64, order: u64) -> MeasurementRecord {
    MeasurementRecord {
        platform: Platform::OlinkExploreNgs,
        assay_id: AssayId(assay.into()),
        gene_symbol: Some(assay.into()),
        sample_id: sample.into(),
        abundance: Abundance::Log2Npx(val),
        abundance_raw: Abundance::Log2Npx(val),
        npx_source_str: format!("{}", val),
        qc_sample: QcFlag::Pass,
        qc_assay: QcFlag::Pass,
        detection_limit: DetectionLimit(None),
        below_lod: false,
        batch: Batch { plate: None, lot: None, run: None },
        dropped_by_qc: false,
        ingest_order: order,
        panel: Some(panel.into()),
    }
}
