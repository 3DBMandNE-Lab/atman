# Changelog

All notable changes to Atman are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Multi-seed FastICA decomposition (`decompose ica`).** Pure-Rust FastICA
  (log-cosh contrast, symmetric decorrelation) with deterministic
  Xoshiro256++ seeding. Runs `--n-seeds` decompositions per cohort and reports
  top-N Jaccard seed stability against the reference seed, plus canonical
  loadings and activations TSVs. Supports `--max-missing-fraction` +
  `--impute mean` for sparse canonical matrices.
- **g:Profiler REST wrapper (`enrich gprofiler`).** Live POST to the
  g:Profiler `gost/profile` endpoint with on-disk cache keyed by a SHA-256
  of the canonicalized request (genes, background, organism, sources,
  threshold method, user threshold, pinned ontology version). `--offline`
  mode fails on cache miss so pre-populated caches produce byte-identical
  runs without a network round-trip. Accepts either a DE-results TSV
  (single query of significant genes) or a `query\tgene_symbol` TSV
  (per-program queries for ICA program annotation).
- **Cross-cohort program alignment (`align programs`).** Single-config and
  `--sweep` mode for Jaccard (top-N), cosine, and Spearman similarity with
  optional reciprocal-best and category-agreement filters. Archetypes are
  connected components of the filtered bipartite graph across all cohort
  pairs. Sweep mode emits an `(metric, top_n, tau, category_constraint)`
  sensitivity matrix so the annotation-constrained vs unconstrained
  comparison is citable rather than ad-hoc.

## [1.0.0] — 2026-04-19

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
- **Adapter examples.** Tiny synthetic Spectronaut, DIA-NN, MaxQuant/LFQ, and
  SomaScan-style fixtures demonstrate `ingest-matrix` mappings and validate
  cleanly.
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
- **Protein bootstrap (`bootstrap protein`).** Subject-level bootstrap
  intervals and sign stability for paired or unpaired protein effects.
- **Module bootstrap (`bootstrap module`).** Subject-level bootstrap intervals,
  sign stability, and gene coverage for user-defined module scores.
- **Null calibration (`null`).** Label permutation for unpaired Welch tests and
  paired sign flips for matched designs, with contrast-level null summaries and
  per-protein empirical p-values.
- **Formula OLS (`de --test ols --design`).** Formula-style covariate-adjusted
  linear models with explicit contrasts, design-matrix reporting, and
  compatibility with the existing `--covariates` shortcut.
- **Mixed-effects DE (`de --test mixed`).** Initial repeated-measures model
  with fixed-effect formulas and a REML-profiled random intercept for
  `subject_id`.
- **Robust DE sidecars.** `de_results.tsv` now appends effect sizes,
  confidence intervals, Wilcoxon p-values, median differences, and trimmed
  mean differences while preserving the original leading columns.
- **Module scoring (`score modules`).** Per-sample module scores with mean,
  median, z-score, and PC1 methods, coverage reporting, and optional canonical
  outputs for downstream DE.
- **ORA enrichment (`enrich ora`).** Built-in over-representation analysis for
  DE hits with measured-universe defaults, optional explicit universes, Fisher
  upper-tail p-values, odds ratios, and BH q-values.
- **Meta-analysis (`meta`).** Combine DE results across cohorts with
  fixed-effect, random-effects, Stouffer, sign-consistency, and heterogeneity
  summaries.

### Fixed

- Neumaier-compensated summation in `mean_linear` to eliminate
  order-dependent float drift (`cfa76c5`).
- Project package naming finalized as `atman` and `atman-core`.

### Notes

- 104 tests in `cargo test --workspace --release`, including an integration
  test that diffs every output against the published Dube reference files.
- No Python runtime is required for the Atman CLI.

[1.0.0]: https://github.com/kevinj24fr/atman/releases/tag/v1.0.0
