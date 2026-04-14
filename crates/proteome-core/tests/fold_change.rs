use proteome_core::fold_change::{compute_log2_fc, Comparison, FoldChangeInput, FoldChangeOutput};

#[test]
fn hand_computed_three_participants_two_exposures() {
    // PR1 values: 0.0, 1.0, 2.0  -> 2^v = 1, 2, 4 -> mean = 7/3
    // PT1 values: 1.0, 2.0, 3.0  -> 2^v = 2, 4, 8 -> mean = 14/3
    // FC(PT1-PR1) = log2((14/3) / (7/3)) = log2(2) = 1.0
    let input = FoldChangeInput::from_cells(vec![
        ("P", "A", "P1", "PR1", Some(0.0)),
        ("P", "A", "P2", "PR1", Some(1.0)),
        ("P", "A", "P3", "PR1", Some(2.0)),
        ("P", "A", "P1", "PT1", Some(1.0)),
        ("P", "A", "P2", "PT1", Some(2.0)),
        ("P", "A", "P3", "PT1", Some(3.0)),
    ]);
    let comps = vec![Comparison { a: "PT1".into(), b: "PR1".into() }];
    let out: FoldChangeOutput = compute_log2_fc(&input, &comps);
    assert_eq!(out.panels.len(), 1);
    let panel = &out.panels[0];
    assert_eq!(panel.panel, "P");
    assert_eq!(panel.assays, vec!["A".to_string()]);
    let fc = panel.values[0][0];
    assert!((fc.unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn missing_in_one_group_drops_participant_not_panel() {
    let input = FoldChangeInput::from_cells(vec![
        ("P", "A", "P1", "PR1", Some(0.0)), // 2^0 = 1
        ("P", "A", "P2", "PR1", Some(2.0)), // 2^2 = 4; mean = 2.5
        ("P", "A", "P1", "PT1", Some(1.0)), // 2^1 = 2; mean = 2
    ]);
    let comps = vec![Comparison { a: "PT1".into(), b: "PR1".into() }];
    let out = compute_log2_fc(&input, &comps);
    let fc = out.panels[0].values[0][0].unwrap();
    let expected = (2.0_f64 / 2.5_f64).log2();
    assert!((fc - expected).abs() < 1e-12, "fc={}, expected={}", fc, expected);
}

#[test]
fn empty_group_produces_missing_fc() {
    let input = FoldChangeInput::from_cells(vec![
        ("P", "A", "P1", "PR1", Some(1.0)),
    ]);
    let comps = vec![Comparison { a: "PT1".into(), b: "PR1".into() }];
    let out = compute_log2_fc(&input, &comps);
    assert!(out.panels[0].values[0][0].is_none());
}

#[test]
fn self_vs_self_is_zero() {
    let input = FoldChangeInput::from_cells(vec![
        ("P", "A", "P1", "PR1", Some(2.0)),
        ("P", "A", "P2", "PR1", Some(3.0)),
    ]);
    let comps = vec![Comparison { a: "PR1".into(), b: "PR1".into() }];
    let out = compute_log2_fc(&input, &comps);
    assert_eq!(out.panels[0].values[0][0], Some(0.0));
}
