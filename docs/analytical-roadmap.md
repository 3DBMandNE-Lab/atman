# Atman Analytical Roadmap

This roadmap turns Atman from a reproducible proteomics CLI into a broader
standalone analysis tool. The order is intentional: validate data first,
improve input coverage second, then add stronger inference.

## Guiding Principles

- Keep Atman local and file-based.
- Keep the canonical TSV contract stable.
- Prefer explicit, auditable transformations over automatic correction.
- Make every new command testable on tiny fixtures.
- Treat biological replicate units, missingness, and batch structure as
  first-class analysis concerns.

## Phase 1: Trust The Input

### 1. `atman validate` (implemented)

Schema and consistency validator for canonical input directories.

Command:

```bash
atman validate --input-dir out
```

Checks:

- Required files exist: `samples.tsv`, `proteins.tsv`, and either
  `measurements.tsv` or `qc_measurements.tsv`.
- Required columns exist in each file.
- `sample_id` values match between samples and measurements.
- `assay_id` values match between proteins and measurements.
- `(platform, assay_id, sample_id)` measurement keys are unique.
- `platform` values parse.
- `abundance_unit` values are known or reported as raw.
- `condition` is present for non-control samples.
- `is_control` parses as `0` or `1`.
- `dropped_by_qc`, `below_lod`, and QC flags parse.
- Each requested comparison has enough effective samples.

Outputs:

- terminal summary suitable for CI
- optional `validate_report.tsv`

Acceptance:

- Fails with a clear message for missing columns.
- Fails on duplicate measurement keys.
- Warns, but does not fail, on unknown abundance units.
- Passes Dube fixture outputs.
- Passes a tiny non-Olink canonical TSV fixture.

Tests:

- Integration tests for valid, missing-column, duplicate-key, warning/strict,
  and non-Olink directories.

### 2. `atman report qc` (implemented)

Compact QC/missingness report.

Command:

```bash
atman report qc --input-dir out --output-dir out/report
```

Metrics:

- sample count, protein count, measurement count
- missing/effective abundance by sample
- missing/effective abundance by protein
- QC-masked fraction
- below-LOD fraction
- condition-level sample counts
- effective sample count per condition and protein
- top sample and protein outliers by missingness

Outputs:

- `qc_summary.tsv`
- `sample_qc.tsv`
- `protein_qc.tsv`
- `condition_counts.tsv`

Acceptance:

- Works on Olink and non-Olink canonical TSVs.
- No plotting dependency.
- Clear warnings for small groups and sparse proteins.
- Integration test covers summary, sample, protein, and condition tables.

## Phase 2: Make Input Broad And Boring

### 3. `atman ingest-matrix` (implemented)

Supported CLI command for common wide protein matrices.

Command:

```bash
atman ingest-matrix \
  --matrix matrix.tsv \
  --samples sample_metadata.tsv \
  --proteins protein_metadata.tsv \
  --orientation proteins-rows \
  --platform diann_report \
  --abundance-unit log2_diann_pg_quantity \
  --condition-col diagnosis \
  --sample-id-col sample_id \
  --assay-id-col Protein.Group \
  --gene-col Genes \
  --log2-transform \
  --output-dir out_ms
```

Supported layouts:

- proteins as rows, samples as columns
- samples as rows, proteins as columns

Supported sources by schema, not by brand-specific parser:

- log2 LFQ/intensity matrices
- Spectronaut protein-group quantity tables
- DIA-NN protein-group matrices
- SomaScan-style log2 abundance matrices
- already-normalized protein tables

Acceptance:

- Replaces the generic Python adapter for common cases.
- Writes `samples.tsv`, `proteins.tsv`, `measurements.tsv`, and
  `qc_measurements.tsv`.
- Preserves extra sample covariates after the canonical sample columns.
- Can log2-transform positive linear intensities.
- QC-masks non-positive values when log-transforming.

Tests:

- proteins-rows fixture
- samples-rows fixture
- log2-transform fixture
- covariate preservation fixture
- failure on unmatched sample IDs

### 4. Adapter Examples (implemented)

Keep `adapters/` as examples and templates for source-specific parsing that is
too irregular for `ingest-matrix`.

Add examples for:

- Spectronaut protein groups
- DIA-NN protein groups
- MaxQuant/LFQ matrix
- SomaScan-style abundance matrix

Acceptance:

- No absolute paths.
- No private datasets.
- Each example has a tiny synthetic fixture.
- Each example writes files that pass `atman validate`.

## Phase 3: Quantify Uncertainty

### 5. `atman bootstrap protein` (implemented)

Subject-level bootstrap uncertainty for protein effects.

Command:

```bash
atman bootstrap protein \
  --input-dir out \
  --groups Case-Control \
  --test welch-t \
  --n 2000 \
  --output out/protein_bootstrap.tsv
```

Outputs:

- point effect
- bootstrap mean effect
- confidence interval
- sign stability
- effective sample counts
- skip reason

Acceptance:

- Resamples subjects, not measurement rows.
- Paired mode resamples paired subjects.
- Unpaired mode resamples within condition.
- Seeded runs are deterministic.
- Integration tests cover seeded unpaired output and paired resampling.

### 6. `atman bootstrap module` (implemented)

Same as protein bootstrap, but over module scores.

Inputs:

- canonical TSV directory
- `modules.tsv` with `module, gene_symbol`

Outputs:

- module effect interval
- sign stability
- effective sample counts
- skip reason
- number of module genes declared
- number of module genes observed

Acceptance:

- Shares bootstrap engine with protein mode where possible.
- Handles missing module members.
- Seeded runs are deterministic.
- Integration tests cover seeded unpaired output, paired resampling, and missing module members.

### 7. `atman null` (implemented)

Add permutation and sign-flip calibration.

Commands:

```bash
atman null --input-dir out --output-dir out/null --groups Case-Control --test welch-t --n 1000
atman null --input-dir out --output-dir out/null --groups PT1-PR1 --test paired-t --n 1000
```

Null modes:

- unpaired label permutation
- paired sign flip
- optional covariate-preserving permutation later

Outputs:

- null hit counts at q thresholds
- empirical FDR summary
- per-protein empirical p-values where feasible

Acceptance:

- Seeded deterministic runs.
- Refuses invalid paired/unpaired designs.
- Reports null calibration by contrast.
- Integration tests cover seeded unpaired permutation, paired sign flips, and invalid paired designs.

## Phase 4: Better Models

### 8. Formula-Based OLS (implemented)

Replace the current covariate list with a formula interface while preserving
the existing `--covariates` path as a compatibility shortcut.

Command:

```bash
atman de --test ols --design "~ condition + age + sex + batch" \
  --contrast conditionCase
```

Capabilities:

- intercept
- numeric covariates
- categorical covariates with reference level
- missing covariate complete-case handling
- clear design-matrix report
- `--covariates` compatibility shortcut expands to formula-style OLS

Acceptance:

- Reproduces current OLS output for equivalent designs.
- Emits `de_covariates.tsv` with stable coefficient names.
- Fails clearly on singular designs.
- Emits `de_design.tsv`.
- Integration tests cover formula/shortcut equivalence, inferred two-condition
  contrasts, and singular design refusal.

### 9. Mixed-Effects Models (implemented)

Add repeated-measures modeling after formula OLS is stable.

Command:

```bash
atman de --test mixed --fixed "condition + age + sex" --random "1|subject_id"
```

Initial scope:

- random intercept only
- one grouping factor: `subject_id`
- restricted maximum likelihood profiled over the random-intercept variance
  ratio, followed by GLS fixed-effect inference

Acceptance:

- Tiny fixture matches a reference implementation within tolerance.
- Refuses underpowered or singular random-effects designs.
- Integration tests cover repeated-subject fitting and unrepeated-subject
  refusal.

### 10. Robust Tests And Effect Sizes (implemented)

Add effect-size-first outputs and robust alternatives.

Features:

- Cohen's dz for paired designs
- Hedges' g for unpaired designs
- log2 fold-change confidence intervals
- Wilcoxon signed-rank and rank-sum
- median difference
- trimmed mean difference

Acceptance:

- Existing `de_results.tsv` remains backward compatible, with extra columns
  appended or a versioned sidecar.
- Integration tests cover paired Cohen's dz/signed-rank output and unpaired
  Hedges' g/rank-sum output.

## Phase 5: Biological Summaries

### 11. Module And Pathway Scoring (implemented)

Add per-sample module scores.

Command:

```bash
atman score modules --input-dir out --modules-tsv modules.tsv --output module_scores.tsv
```

Methods:

- mean abundance
- mean z-score
- median abundance
- first principal component

Acceptance:

- Writes per-sample scores usable by `de`.
- Handles missing module members and reports coverage.
- Optional `--canonical-output-dir` writes module-score canonical TSVs for
  downstream DE.
- Integration tests cover all scoring methods, coverage reporting, and DE over
  canonical module scores.

### 12. Enrichment (implemented)

Add simple built-in over-representation analysis first.

Command:

```bash
atman enrich ora --de-results out/de_results.tsv --gene-sets gene_sets.tsv --output ora.tsv
```

Inputs:

- `gene_sets.tsv`: `set_name, gene_symbol`
- optional measured-universe file

Acceptance:

- Uses the measured protein universe by default.
- Reports set size, overlap size, odds ratio, p-value, and BH q-value.
- Optional explicit universe file.
- Integration tests cover measured-universe defaults, comparison filtering,
  overlap reporting, and explicit universe files.

### 13. Multi-Cohort Meta-Analysis (implemented)

Combine results across independently processed cohorts.

Command:

```bash
atman meta --inputs cohort1/de_results.tsv,cohort2/de_results.tsv \
  --method fixed-effect --output meta.tsv
```

Methods:

- fixed-effect inverse variance
- random-effects summary
- Stouffer p-value combination
- sign consistency table

Acceptance:

- Handles missing proteins across cohorts.
- Reports heterogeneity metrics.
- Reports fixed-effect, random-effects, Stouffer, and sign-consistency fields.
- Integration tests cover missing proteins, sign inconsistency, and
  heterogeneity metrics.

## Phase 6: Peptide-Level Inference

### 14. Peptide-level ridge mixed model (`atman de --test msqrob`) (implemented)

Per-protein ridge-regularized linear mixed model fit over peptide-level
measurements. Random intercept per peptide, REML 1D profile over the
variance ratio, L2 ridge on non-intercept fixed-effect coefficients.

Command:

```bash
atman de \
  --input-dir out \
  --output-dir out_msqrob \
  --test msqrob \
  --peptide-measurements out/peptide_measurements.tsv \
  --peptide-metadata out/peptides.tsv \
  --groups Case-Control \
  --ridge-lambda 0.5 \
  --min-peptides 2 \
  --min-pairs 5
```

New canonical files (additive; protein-level files unchanged):

- `peptide_measurements.tsv`: `sample_id, peptide_id, abundance,
  abundance_unit, dropped_by_qc, below_lod`
- `peptides.tsv`: `peptide_id, assay_id, sequence, charge, modifications,
  missed_cleavages`

Acceptance:

- `de_results.tsv` gains `n_peptides_observed`, `peptide_variance_ratio`,
  `ridge_lambda` columns for msqrob rows (empty for other test types).
- `effect_size_method = "msqrob-ridge"` on every fitted row.
- Sign convention matches the rest of `atman de`: reported effect is
  `mean_a − mean_b` for an `A-B` comparison.
- `--ridge-lambda auto` parses but currently equals `0.0`; data-driven
  selection is a follow-on.
- `--min-peptides` refuses proteins with too few observed peptides,
  emitting `skip_reason = "insufficient_peptides"`.
- Integration tests cover direction recovery on a 3-protein synthetic
  fixture, column presence, `auto` default, CLI rejection of peptide
  flags under non-msqrob tests, and failure on missing peptide inputs.

Follow-ons this unlocks:

- **DEqMS peptide-count-weighted variance** on top of `--test limma`
  **(implemented, 2026-04-20):** `atman de --test limma
  --peptide-metadata peptides.tsv` swaps the mean-variance trend
  covariate for `log(peptide_count + 1)` smoothed by a tricube
  kernel, mirroring `DEqMS::spectraCounteBayes` (Zhu et al. 2020).
  Validated to three-decimal parity with Bioconductor DEqMS v1.26 on
  the CPTAC Study 6 UPS1 spike-in fixture — median per-protein
  |atman − DEqMS| = 0.000 log₂ across 29 jointly-fitted proteins.
- **proDA-style probabilistic missingness** (joint abundance + detection
  likelihood) using the same peptide schema.

## Phase 7: Cross-method Consensus

### 16. Ensemble DE dispatcher (`atman de --test ensemble`) (implemented)

Runs every applicable DE method on the same canonical inputs,
Stouffer-combines per-method p-values within each (comparison,
protein), BH-FDRs ensemble p-values across proteins, and assigns a
VALIDATED / PROVISIONAL / INSUFFICIENT grade from the combined
evidence plus cross-method sign agreement. Directly answers the
"is this hit a model-choice artifact?" reviewer question that every
proteomics paper faces — no other bulk-proteomics tool ships this.

Command:

```bash
atman de --input-dir out --output-dir out_ensemble \
  --test ensemble \
  --ensemble-methods "paired-t,welch-t,limma,msqrob" \
  --peptide-measurements out/peptide_measurements.tsv \
  --peptide-metadata out/peptides.tsv \
  --groups "PT2-PR2" --paired-by participant --min-pairs 5
```

Grading logic (post-BH):

- **VALIDATED:** `ensemble_q < q_threshold` AND 100% of methods
  agree on mean_diff sign.
- **PROVISIONAL:** `ensemble_q < q_threshold` AND
  `≥ provisional_sign_fraction` (default 0.50) agree on sign.
- **INSUFFICIENT:** otherwise.

Acceptance:

- `de_ensemble.tsv` emits one row per (comparison, protein) with
  `n_applied`, `n_significant`, `n_sign_consistent`, `majority_sign`,
  `ensemble_p` (Stouffer), `ensemble_q` (BH within comparison),
  `grade`, `methods_applied`, `methods_skipped`.
- `de_results.tsv` gains a `method` column for bucketing per-method
  rows.
- Methods whose required inputs are missing are auto-skipped and
  reported in `methods_skipped`; the run does not abort.
- Implementation uses per-method sub-dispatch into tempdirs via
  `Args: Clone`, re-reading each method's `de_results.tsv` before
  aggregating. No refactor of the monolithic per-method paths was
  required.
- Integration tests cover (i) grade assignment on a synthetic 3-protein
  fixture (UP/DN/ST known ground truth); (ii) auto-skip behaviour
  when msqrob inputs are absent; (iii) real-data grading on the
  bundled Dube heat-acclimation cohort — HSPA1A, HSPB1, DNAJB1 all
  graded VALIDATED in PT2-PR2.

## Implemented Build Order

1. `atman validate`
2. `atman report qc`
3. `atman ingest-matrix`
4. Adapter examples with tiny fixtures
5. `atman bootstrap protein`
6. `atman bootstrap module`
7. `atman null`
8. Formula-based OLS
9. Robust effect-size outputs
10. Mixed-effects models
11. Module scoring
12. ORA enrichment
13. Meta-analysis
14. Peptide-level ridge mixed model
15. DEqMS peptide-count-weighted variance
16. Cross-method ensemble DE dispatcher

## Release Readiness Status

Atman now covers the standalone release surface targeted by this roadmap:

- Canonical validation and QC reporting are implemented.
- Olink and matrix ingest paths are covered by integration tests and adapter
  fixtures.
- Differential abundance, uncertainty, robustness, module, enrichment, and
  meta-analysis commands are implemented.
- The core CLI does not require Python, notebooks, services, or external
  workflow infrastructure.
