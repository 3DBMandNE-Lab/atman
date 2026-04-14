use proteome_core::{DubeSampleIdParser, ParsedSampleId, SampleIdParser};

#[test]
fn parses_dube_biological_sample_ids() {
    let p = DubeSampleIdParser;
    let r = p.parse("SSNA-001B-PR1").unwrap();
    assert_eq!(
        r,
        ParsedSampleId::Biological {
            subject: "001B".into(),
            condition: "PR1".into()
        }
    );

    let r = p.parse("SSNA-019-PT2").unwrap();
    assert_eq!(
        r,
        ParsedSampleId::Biological {
            subject: "019".into(),
            condition: "PT2".into()
        }
    );
}

#[test]
fn classifies_control_sample_as_control() {
    let p = DubeSampleIdParser;
    assert_eq!(
        p.parse("CONTROL_SAMPLE_US_CS_AS_2-1").unwrap(),
        ParsedSampleId::Control
    );
    assert_eq!(
        p.parse("CONTROL_SAMPLE_US_CS_AS_2-4").unwrap(),
        ParsedSampleId::Control
    );
}

#[test]
fn rejects_malformed_non_control() {
    let p = DubeSampleIdParser;
    assert!(p.parse("RANDOM-001B-PR1").is_err());
}

#[test]
fn rejects_ssna_with_unknown_timepoint() {
    let p = DubeSampleIdParser;
    assert!(p.parse("SSNA-001B-PX9").is_err());
}

#[test]
fn rejects_ssna_without_subject() {
    let p = DubeSampleIdParser;
    assert!(p.parse("SSNA--PR1").is_err());
}
