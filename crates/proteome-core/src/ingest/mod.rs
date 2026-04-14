//! Ingest boundary. Adapters implement `ProteomeIngest` to produce `IngestOutput`
//! from raw files. v0.1 ships exactly one: `olink_explore::OlinkExploreLongCsv`.

pub mod olink_explore;

use crate::{
    errors::IngestError,
    sample_id::SampleIdParser,
    types::{MeasurementRecord, Platform, ProteinIdentity, Sample},
};
use std::path::PathBuf;

/// Canonical output of any ingest adapter.
#[derive(Debug, Default)]
pub struct IngestOutput {
    pub measurements: Vec<MeasurementRecord>,
    pub proteins: Vec<ProteinIdentity>,
    pub samples: Vec<Sample>,
}

pub trait ProteomeIngest {
    fn platform(&self) -> Platform;

    /// Read one or more input files. `sample_id_parser` is consulted for every
    /// unique `sample_id`; biological samples populate `Sample::subject_id` and
    /// `condition`, controls get `is_control=true` with `subject_id=None`.
    fn read(
        &self,
        inputs: &[PathBuf],
        sample_id_parser: &dyn SampleIdParser,
    ) -> Result<IngestOutput, IngestError>;
}
