//! atman-core: platform-agnostic primitives and algorithms for proteomics data.
//!
//! The current package implements an Olink Explore NGS long-CSV adapter and
//! the transformations required to reproduce the Dube et al. Scientific Data
//! 2023 published filtered NPX + log2 fold-change files.

pub mod align;
pub mod de;
pub mod deqms;
pub mod errors;
pub mod fold_change;
pub mod ica;
pub mod ingest;
pub mod limma;
pub mod matrix;
pub mod msqrob;
pub mod qc;
pub mod sample_id;
pub mod stats;
pub mod types;

pub use de::{bh_fdr, paired_t, PairedTResult, SkipReason};
pub use deqms::{deqms_shrink, tricube_moving_average, DeqmsShrinkage};
pub use msqrob::{fit_msqrob, squeeze_variance, MsqrobFit, MsqrobOutcome};
pub use errors::IngestError;
pub use fold_change::{
    compute_log2_fc, Comparison, FoldChangeInput, FoldChangeOutput, FoldChangePanel,
};
pub use ingest::{IngestOutput, ProteomeIngest};
pub use matrix::{DubeWidePanel, DubeWideRow};
pub use sample_id::{DubeSampleIdParser, ParsedSampleId, SampleIdParser};
pub use types::{
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, PeptideIdentity,
    PeptideMeasurementRecord, Platform, ProteinIdentity, QcFlag, Sample,
};
