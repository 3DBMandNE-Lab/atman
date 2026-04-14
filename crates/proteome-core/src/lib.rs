//! proteome-core: platform-agnostic primitives and algorithms for proteomics data.
//!
//! v0.1 implements exactly one ingest adapter (Olink Explore NGS long CSV) and
//! the transformations required to reproduce the Dube et al. Scientific Data 2023
//! published filtered NPX + log2 fold-change files. See
//! `docs/superpowers/specs/2026-04-14-karnaproteome-olink-reproduction-design.md`
//! for the full design.

pub mod errors;
pub mod fold_change;
pub mod ingest;
pub mod matrix;
pub mod qc;
pub mod sample_id;
pub mod types;

pub use errors::IngestError;
pub use fold_change::{
    compute_log2_fc, Comparison, FoldChangeInput, FoldChangeOutput, FoldChangePanel,
};
pub use ingest::{IngestOutput, ProteomeIngest};
pub use matrix::{DubeWidePanel, DubeWideRow};
pub use sample_id::{DubeSampleIdParser, ParsedSampleId, SampleIdParser};
pub use types::{
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, ProteinIdentity,
    QcFlag, Sample,
};
