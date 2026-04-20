# Changelog

All notable changes to Atman are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Feature-covariance network influence (`atman network influence`).**
  Builds a feature × feature similarity graph over subjects
  (`--method {pearson,spearman,covariance}`), applies a hard
  threshold or WGCNA-style soft-power adjacency
  (`--threshold <f>` or `--soft-power <β>`), and scores each
  feature by its role as a hub via eigenvector centrality ×
  betweenness centrality (Burberry-Pillai 2026 construction).
  Optional `--stratify <column>` builds one graph per stratum.
  Output `network_influence.tsv` carries `feature_id, stratum,
  eigenvector_centrality, betweenness_centrality, influence_score,
  degree, n_subjects_used`. Pure math in
  `atman-core::network` (8 unit tests covering star-topology hub
  recovery for both centralities, clique symmetry, adjacency
  policies, and determinism). Eigenvector iteration uses an
  `A + αI` shift (α = max row sum) so bipartite graphs don't
  oscillate. Smoke-tested on Dube (2938 features, ~600 degree per
  hub).
- **Multi-level OLS omnibus F-test (`atman de --test ols
  --omnibus-factor <name>`).** For categorical factors with ≥ 3
  levels (e.g. `stage` ∈ {CN, MCI, AD} or `diagnosis` across six
  disease groups), emits a per-protein F-test of the joint
  hypothesis that every level of the factor has zero coefficient.
  Output `de_omnibus.tsv` carries `panel, assay_id, gene_symbol,
  factor, comparison, f_statistic, df_num, df_den, p_value, bh_q`
  (BH-adjusted within each `(comparison, panel)` family). The
  factor's design columns are identified by prefix match on
  `design_labels` from the existing OLS design builder, so no
  additional formula machinery is needed. Atomic integration test
  shows F≈21 / p<10⁻⁵ on a planted-responder protein and F≈0.2 /
  p=0.83 on a stage-independent protein across a 24-sample 3-stage
  fixture. Partially closes P9 of the methods track.
  *Deferred to a follow-on:* Tukey HSD, Dunnett, and Sidak post-hoc
  pairwise contrasts. Tukey/Dunnett require distributions not in
  statrs (studentized range, multivariate-t); Sidak on a
  user-supplied `--contrast-list` is tracked as the next increment.
- **Archetype variance decomposition (`atman decompose variance`).**
  Partitions each archetype's subject-level activation variance into
  per-fixed-factor contributions, a random-intercept component, and
  residual variance via a linear mixed model. Reuses
  `atman de --test mixed`'s REML engine; adds only the partitioning
  math on top. Formula syntax mirrors R:
  `--factors "cohort + condition + (1|subject_id)"`. Output
  `archetype_variance.tsv` has one row per archetype with
  `total_var`, `var_residual`, `var_random_<group>`,
  `icc_random_<group>`, and per-factor `var_<factor>`,
  `max_abs_t_<factor>`, `min_p_<factor>`. v1 partition is the
  Type-I projection variance (`var_f = Var(X_f · β_f)`);
  correlation-adjusted Type II/III partitions are a documented
  follow-on. Per-factor omnibus F-tests are deferred — v1 reports
  per-coefficient Wald summaries. Closes P5 of the decomposition
  methods track.
- **Subject-level bootstrap of cross-cohort archetype alignment
  (`atman align bootstrap`).** Takes `--cohorts dir_a,dir_b[,...]`
  (canonical Atman directories), runs a point-estimate FastICA per
  cohort + cosine-metric alignment to identify archetypes, then for
  each of `--n-boot` iterations resamples subjects within every
  cohort with replacement, re-runs ICA per cohort, re-aligns, and
  matches every point-estimate archetype to its best-cosine
  bootstrap counterpart. Emits
  `align_bootstrap_summary.tsv` with `archetype_id,
  observed_n_cohorts, observed_cohorts, bootstrap_mean_n_cohorts,
  bootstrap_prob_universal, bootstrap_prob_multi,
  ci_lower_n_cohorts, ci_upper_n_cohorts, bootstrap_match_rate`.
  Deterministic under `--seed` via SplitMix64`(seed, iter)` sub-seed
  derivation. v1 is cosine-only, percentile-CI (no BCa), and no
  `alignment_entropy` metric — noted as follow-ons. Partially
  closes P3 of the decomposition methods track.
- **Compositional transforms in `atman decompose ica`
  (`--transform {none,clr,alr,ilr,ratio-anchor}`).** Bulk proteomics
  is closed-sum (total intensity is a plate/panel artifact), so
  Euclidean-geometry ICA on raw log2 abundance mixes biology with
  scaling. The transform is now a first-class CLI option with a
  per-run `transform_applied.json` audit file recording the
  geometry the archetypes live in. `clr` is sample-wise
  mean-centering on log data; `alr` / `ratio-anchor` subtract a
  reference gene's column (use `--alr-reference <gene_symbol>`);
  `ilr` projects through a Helmert orthonormal basis (output
  loadings in `ilr_coord_*` coordinates rather than raw proteins,
  since ILR shrinks dimensionality by 1). Run sidecar captures
  `transform` and `alr-reference`. Zero-handling is deferred
  until linear-scale input routes exist — atman's canonical data
  arrives finite on a log scale after QC, so no zero mapping is
  needed on the common paths.
- **Archetype null calibration (`atman decompose null`).** Permutation
  null for FastICA archetype stability: runs the multi-seed Jaccard
  stability metric on the real matrix, then runs it again on
  `--n-perm` null matrices generated via `--null-mode` (one of
  `sample-shuffle`, `protein-shuffle`, `gaussian-matched`). Each null
  iteration contributes its max-across-programs stability to a null
  distribution; per-program p-value uses the `(1 + count)/(1 + n_perm)`
  permutation correction, BH-adjusted across programs within the run
  to produce `null_q`. Output `archetype_null.tsv` columns:
  `program, observed_stability, null_stability_mean,
  null_stability_p95, null_p, null_q, decision` (decision is
  `signal` when `null_q < --q-threshold` (default 0.05), else
  `noise`). Fully deterministic under `--seed`; each null iteration
  derives its sub-seed from SplitMix64`(seed, iter)`. Closes P1 of
  the decomposition methods track ("are these archetypes just
  ICA-finds-whatever?"). Directly addresses a reviewer defense that
  `atman decompose ica` alone could not provide.
- **Cross-method consensus DE (`atman de --test ensemble`).** Runs
  every applicable DE method on the same canonical inputs — paired-t,
  welch-t, ols, mixed, limma (with or without DEqMS), msqrob — and
  writes a per-(comparison, protein) agreement summary to
  `de_ensemble.tsv` alongside the tagged per-method rows in
  `de_results.tsv` (new `method` column). Per-method p-values are
  Stouffer-combined into an `ensemble_p`, BH-FDR'd to `ensemble_q`
  within each comparison; grade is **VALIDATED** when
  `ensemble_q < --ensemble-q-threshold` (default 0.05) AND every
  applicable method agrees on sign; **PROVISIONAL** when
  ensemble_q is significant but at least `--ensemble-provisional-fraction`
  (default 0.50) of methods agree on sign; **INSUFFICIENT**
  otherwise. Method applicability is auto-detected: methods whose
  required inputs are missing (e.g. msqrob without
  `--peptide-measurements`) are listed in `methods_skipped` rather
  than aborting the run. Thresholds are overridable
  (`--ensemble-q-threshold`, `--ensemble-validated-fraction`,
  `--ensemble-provisional-fraction`, `--ensemble-sign-fraction`).
  Validated on the bundled Dube heat-acclimation cohort — all three
  canonical heat-shock proteins present in the Olink panels (HSPA1A,
  HSPB1, DNAJB1) are graded VALIDATED in PT2-PR2 with
  `ensemble_q < 0.025` and positive majority sign.
- **DEqMS peptide-count-weighted variance on `--test limma`
  (`atman de --test limma --peptide-metadata peptides.tsv`).**
  Swaps limma's parametric mean-variance trend covariate for
  `log(peptide_count + 1)` smoothed by a tricube kernel
  (`atman_core::deqms::tricube_moving_average`, span 0.5), matching
  `DEqMS::spectraCounteBayes` (Zhu et al. 2020, MCP). Peptide counts
  come from the ingested `peptides.tsv`; proteins absent from the
  catalog receive count 0. `effect_size_method` in `de_results.tsv`
  becomes `limma-DEqMS-trend` (or `limma-DEqMS-robust-trend` under
  `--robust`) and `n_peptides_observed` is populated per row. Same
  run-sidecar surface — the `peptide-metadata` path is recorded in
  the JSON. Validated against Bioconductor DEqMS v1.26 on the CPTAC
  Study 6 UPS1 spike-in fixture
  (`crates/atman/tests/fixtures/deqms_cptac_reference.R`):
  **median per-protein |atman − DEqMS| = 0.000 log₂** across 29
  jointly-fitted proteins (three-decimal bit-for-bit match on every
  displayed protein).
- **Per-feature complete-case filtering inside `limma_fit`.**
  Features with `NaN` samples used to be rejected wholesale by the
  per-feature OLS call; now each feature drops its own `NaN`
  samples and fits OLS on the remainder (matching limma's
  `lm.series(ndups=1)` behavior). The residual-df used for the
  eBayes prior uses the maximum fitted df across features. Unlocks
  the DEqMS validation on the sparse MaxQuant fixture and
  generally improves coverage on datasets with per-protein
  missingness.
- **Peptide-level ridge mixed model (`atman de --test msqrob`).**
  Per-protein linear mixed model over peptide-level measurements with a
  random intercept per peptide, REML 1D profile over the variance
  ratio `τ = σ²_peptide / σ²_res`, and an L2 ridge penalty on
  non-intercept fixed-effect coefficients. Collapses protein rollup,
  imputation, and per-protein variance estimation into one joint fit in
  the msqrob2 tradition. Reads two new canonical TSVs —
  `peptide_measurements.tsv` (sample × peptide abundances) and
  `peptides.tsv` (peptide catalog with parent `assay_id`) — alongside
  `samples.tsv` and `proteins.tsv`. Output extends `de_results.tsv` with
  three additive columns: `n_peptides_observed`,
  `peptide_variance_ratio`, `ridge_lambda`. Sign convention matches the
  rest of `atman de` (`mean_a − mean_b`). `--ridge-lambda auto` is
  currently equivalent to `0.0` (data-driven selection is a follow-on).
  `--min-peptides` defaults to 2; proteins below the threshold emit
  `skip_reason = "insufficient_peptides"`. Per-comparison
  empirical-Bayes variance shrinkage via limma's `fit_f_dist` across
  all fitted proteins (`atman_core::msqrob::squeeze_variance`) matches
  msqrob2's `squeezeVarRob` step and populates the existing
  `s2_prior`, `s2_posterior`, `df_prior`, `df_total` columns in
  `de_results.tsv`. Validated on the CPTAC Study 6 UPS1 spike-in
  (30-protein MaxQuant peptides.txt subset from statOmics/MSqRobData):
  9/9 UPS1 proteins recover the expected negative A–B direction; 100%
  sign agreement with msqrob2 on all 11 signal proteins (|effect|
  ≥ 0.3 log₂); **median per-protein |atman − msqrob2| = 0.118 log₂**
  across 26 jointly-fitted proteins (test tolerance 0.40), using
  msqrob2 v1.16 run through the canonical vignette workflow
  (`log2 → center.median → median-summarise → msqrob(~condition)`) on
  the identical fixture. Reference script
  `crates/atman/tests/fixtures/msqrob2_cptac_reference.R` regenerates
  the `msqrob2_cptac_reference.tsv` baseline.
- **limma-grade eBayes + F-tests (`atman de --test limma`).** Pure-Rust
  port of limma 3.x's core `eBayes` path: moderated t-statistic with
  empirical-Bayes variance shrinkage, parametric-quadratic mean–variance
  trend (`--trend`), one-pass Winsorized robust prior fit (`--robust`),
  TREAT minimum-effect testing (`--lfc-threshold`). Numerically matches
  `limma::eBayes(trend=FALSE, robust=FALSE)` within `1e-4` on the
  100-feature × 20-sample regression fixture (F001 spot-check diff
  2e-8). Output drops into `de_results.tsv` via nine additive columns
  (`s2_trend`, `s2_prior`, `s2_posterior`, `df_prior`, `df_total`,
  `f_statistic`, `f_p_value`, `f_bh_q`, `lfc_threshold`); `de_report.tsv`
  gains `limma_trend_fallback_used`. Sign convention matches the rest
  of `atman de` (`mean_a − mean_b` for `A-B`). Trend-mode and `robust=true`
  numerical parity against R are deferred follow-ups — see the design
  spec's out-of-scope section.
- **Per-invocation run-sidecar (`<primary_output>.run.json`).** Shared
  helper `atman::io::write_run_sidecar` emits a JSON sidecar next to the
  primary output of each subcommand, capturing the atman version and git SHA,
  the fully-resolved argument dict (including defaults), a SHA-256 hash of
  the canonical input TSVs, SHA-256 hashes of every output file, ISO-8601
  UTC start/finish timestamps, and the build target triple. Moves "what
  parameters produced this file?" from "read the driver script" to "cat one
  JSON." Wired across: `decompose ica`, `align programs`, `coupling`,
  `null`, `enrich gprofiler`, `de`, `bootstrap protein`, `bootstrap module`,
  `bootstrap program`, `meta`, `ratio`, `validate` (when `--report` is
  set), `report qc`, and `ingest-matrix`.
- **Multi-seed FastICA decomposition (`decompose ica`).** Pure-Rust FastICA
  (log-cosh contrast, symmetric decorrelation) with deterministic
  Xoshiro256++ seeding. Runs `--n-seeds` decompositions per cohort and reports
  top-N Jaccard seed stability against the reference seed, plus canonical
  loadings and activations TSVs. Supports `--max-missing-fraction` +
  `--impute mean` for sparse canonical matrices.
- **Pre-registered analysis runner (`run`).** Executes a declarative plan
  file (YAML or JSON) stage by stage via `sh -c`, captures runtime, exit
  code, atman version, OS/arch, and SHA-256 hashes of each stage's
  declared inputs and outputs into `plan_manifest.tsv`. Detects plan
  content drift via `plan_hash` comparison against the previous manifest
  and refuses to overwrite unless the plan's `plan_commit` is bumped (or
  `--allow-drift` is passed). Converts "trust my git log" into a
  hash-verifiable provenance artifact.
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
