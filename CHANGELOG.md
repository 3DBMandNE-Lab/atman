# Changelog

All notable changes to Atman are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] — 2026-04-XX

First public release of Atman as a standalone proteomics command-line tool.
Covers the built-in Olink Explore NPX reproduction path, the canonical TSV
adapter contract for other proteomics matrices, and the core analysis commands
exposed by the Rust binary.

### Added

- **Reproduction base (`ingest`, `qc`, `matrix`, `fold-change`).** Byte-exact
  reproduction of Dube's 8 filtered NPX files and 8 log2 fold-change files
  (max delta 1.05e-15, IEEE 754 last-bit drift).
- **Canonical TSV analysis contract.** Downstream commands operate on
  `samples.tsv`, `proteins.tsv`, `measurements.tsv`, and `qc_measurements.tsv`,
  allowing thin adapters for LFQ/intensity matrices, Spectronaut exports,
  DIA-NN protein-group matrices, SomaScan-style log2 abundance matrices, and
  already-normalized protein tables.
- **Canonical TSV validation (`validate`).** Checks schema, duplicate keys,
  sample/protein references, platform parsing, QC fields, abundance units, and
  effective sample counts before downstream analysis.
- **QC reporting (`report qc`).** Writes dataset, sample, protein, and
  condition-level QC/missingness TSVs for canonical Atman directories.
- **Wide matrix ingest (`ingest-matrix`).** Converts common protein matrix
  layouts plus metadata into canonical Atman TSVs, with optional log2
  transformation and sample covariate preservation.
- **Differential abundance (`de`).** Paired Student's t-test with BH-FDR,
  plus a moderated variance-shrinkage mode (`--test moderated
  --moderation-prior-df 4`).
- **Asymmetry (`asymmetry`).** Provocation-contrast asymmetry metrics for
  paired contrast sets (challenge-vs-rest style).
- **Robustness (`robustness`).** Leave-one-out rank stability summaries from
  baseline + LOO DE runs, including threshold sensitivity tables.
- **Module trajectory (`module-trajectory`).** User-supplied module scoring
  over per-subject gene deltas.
- **Module-level differential abundance (`module-de`).** Aggregate proteins
  into user-supplied modules and test at module level.

### Fixed

- Neumaier-compensated summation in `mean_linear` to eliminate
  order-dependent float drift (`cfa76c5`).
- Project package naming finalized as `atman` and `atman-core`.

### Notes

- 84 tests in `cargo test --workspace --release`, including an integration
  test that diffs every output against the published Dube reference files.
- No Python runtime is required for the Atman CLI.

[1.0.0]: https://github.com/kevinjoseph/atman/releases/tag/v1.0.0
