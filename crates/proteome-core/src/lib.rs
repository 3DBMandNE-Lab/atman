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
