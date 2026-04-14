//! Sample ID parsing. The Dube parser matches `^SSNA-<subject>-<PR1|PR2|PT1|PT2>$`;
//! unparseable IDs with a `CONTROL_SAMPLE_` prefix are classified as controls;
//! anything else is an error.

use crate::errors::IngestError;
use regex::Regex;
use std::sync::OnceLock;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParsedSampleId {
    Biological { subject: String, condition: String },
    Control,
}

pub trait SampleIdParser {
    fn parse(&self, sample_id: &str) -> Result<ParsedSampleId, IngestError>;
}

pub struct DubeSampleIdParser;

impl SampleIdParser for DubeSampleIdParser {
    fn parse(&self, sample_id: &str) -> Result<ParsedSampleId, IngestError> {
        static RE: OnceLock<Regex> = OnceLock::new();
        let re = RE.get_or_init(|| {
            Regex::new(r"^SSNA-(?P<subject>[^-]+)-(?P<condition>PR1|PR2|PT1|PT2)$").unwrap()
        });

        if let Some(caps) = re.captures(sample_id) {
            return Ok(ParsedSampleId::Biological {
                subject: caps["subject"].to_string(),
                condition: caps["condition"].to_string(),
            });
        }

        if sample_id.starts_with("CONTROL_SAMPLE_") {
            return Ok(ParsedSampleId::Control);
        }

        Err(IngestError::UnparseableSampleId {
            sample_id: sample_id.to_string(),
        })
    }
}
