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

- Required files exist: `samples.tsv`, `proteins.tsv`, and `measurements.tsv`.
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
- Writes `samples.tsv`, `proteins.tsv`, and `measurements.tsv`.
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
- `--ridge-lambda auto` is **refused**, not silently zero. Data-driven
  penalty selection remains a follow-on.

  **History, 2026-09-06.** `auto` used to be this flag's default and
  silently resolved to `0.0`, so it read as working while applying no
  shrinkage — the dangerous state for a provenance tool. Fixed: the
  default is now literally `0.0` (numerically unchanged), `auto` is
  refused under `--test msqrob` with a message naming the numeric
  alternative, and it is reported as ignored under any other test. The
  sidecar gained `ridge-lambda-resolved`, the penalty that actually
  applied, which is null when the shrinkage path did not execute.

  Two things to keep in mind when reading **old** sidecars, which the
  fix deliberately does not rewrite:

  - `ridge-lambda: auto` appears in pre-fix sidecars regardless of test,
    because the sidecar records the resolved CLI argument set rather
    than the executed path. Ten such sidecars exist in the CSF
    CrossDisease project, all `test=ols`, none of which ran shrinkage.
    A pre-fix sidecar naming `auto` is not evidence that msqrob ran.
  - Pre-fix sidecars have no `ridge-lambda-resolved` field at all. Its
    absence means "written before 2026-09-06", not "no penalty".
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

Runs every applicable DE method on the same canonical inputs and, for
each (comparison, protein), counts how many methods individually clear
their own per-method BH-q (`n_significant`) and how many agree on the
`mean_diff` sign (`n_sign_consistent`). The VALIDATED / PROVISIONAL /
INSUFFICIENT grade is assigned from that method agreement — NOT from a
combined p-value. A Stouffer `ensemble_p`/`ensemble_q` is still emitted
as a non-calibrated ranking convenience (Stouffer assumes independent
p-values, but the methods share one abundance matrix and are positively
correlated, so the combined p is anti-conservative and must not gate a
grade). Directly answers the "is this hit a model-choice artifact?"
reviewer question that every proteomics paper faces — no other
bulk-proteomics tool ships this.

Command:

```bash
atman de --input-dir out --output-dir out_ensemble \
  --test ensemble \
  --ensemble-methods "paired-t,welch-t,limma,msqrob" \
  --peptide-measurements out/peptide_measurements.tsv \
  --peptide-metadata out/peptides.tsv \
  --groups "PT2-PR2" --min-pairs 5
```

Grading logic (anchored to per-method significance, not the combined p):

- **VALIDATED:** every applied method is individually significant
  (`n_significant / n_applied ≥ validated_sign_fraction`, default 1.0)
  AND the same fraction agree on `mean_diff` sign.
- **PROVISIONAL:** a majority of methods are individually significant
  (`≥ provisional_sign_fraction`, default 0.50) AND the same fraction
  agree on sign.
- **INSUFFICIENT:** otherwise (including when no method clears its own
  significance, regardless of how small the Stouffer `ensemble_q` is).

Acceptance:

- `de_ensemble.tsv` emits one row per (comparison, protein) with
  `n_applied`, `n_significant`, `n_sign_consistent`, `majority_sign`,
  `ensemble_p` (Stouffer — non-calibrated heuristic, not a grade
  input), `ensemble_q` (BH within comparison, same caveat),
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
  bundled Dube heat-acclimation cohort — the canonical HSPs (HSPA1A,
  HSPB1, DNAJB1) show consistent positive sign, but at n≈9–20 paired
  subjects no single method clears its own BH-q, so honest grading does
  NOT mark them VALIDATED. The test asserts sign consistency and that
  the grade never exceeds per-method significance support (the prior
  VALIDATED label came only from the anti-conservative combined p).

## Phase 8: Decomposition Rigor

### 17. Archetype null calibration (`atman decompose null`) (implemented)

Permutation-null calibration for FastICA archetype stability. Runs
multi-seed ICASSO stability on the real matrix and on `--n-perm`
null matrices (sample-shuffle / protein-shuffle / gaussian-matched
modes), records each null run's max-across-programs stability, and
reports a per-program permutation p-value with `(1+count)/(1+n_perm)`
correction plus BH-q across programs.

Command:

```bash
atman decompose null \
  --input-dir out --output out/archetype_null.tsv \
  --k 20 --n-perm 200 --n-seeds 3 --seed 20260418 \
  --null-mode protein-shuffle --top-n 20 \
  --max-missing-fraction 0.5 --impute mean
```

Null modes:

- `sample-shuffle`: permute the row order of the sample × protein
  matrix. Weakest null — preserves inter-protein covariance, so ICA
  finds the same archetypes in the null. Useful only for testing
  activation-sample coherence.
- `protein-shuffle`: independently permute samples within each
  protein column. Destroys inter-protein covariance. Recommended
  default null for archetype-existence questions.
- `gaussian-matched`: per-protein `N(μ, σ²)` with matched first two
  moments. Strongest null — no structure beyond mean/variance.

Acceptance:

- Deterministic under `--seed`; every null iteration derives its
  sub-seed from SplitMix64`(seed, iter)`.
- Output includes `decision` column (`signal` when `null_q <
  --q-threshold`, else `noise`); threshold configurable via
  `--q-threshold` (default 0.05).
- Run sidecar captures `null-mode`, `n-perm`, `seed`, `top-n`,
  `q-threshold`, and SHA-256 of the output TSV.
- Unit tests verify null-matrix generators preserve expected
  marginal invariants (row permutation preserves column multisets,
  column shuffle ditto per-column, Gaussian-matched recovers mean
  and variance within tolerance at n=2000).
- Integration tests cover synthetic planted-signal schema +
  determinism + k-size refusal, plus an end-to-end smoke test on
  the bundled Dube cohort.

Why it matters: no standalone proteomics ICA tool ships archetype-
level permutation null calibration. MOFA, consICA, fastICA+ICASSO
rely on user-side scripting. Atman now provides it as a first-class
CLI surface with sidecar provenance.

### 18. Compositional transforms in `atman decompose ica` (implemented)

First-class `--transform` flag exposing CLR, ALR, ILR, and
ratio-anchor projections for closed-sum proteomics data. Writes a
`transform_applied.json` audit file next to `loadings.tsv` so the
coordinate system the archetypes live in is recorded in the run's
provenance trail.

Command:

```bash
atman decompose ica \
  --input-dir out --k 20 \
  --transform clr \
  --seed 20260418 \
  --output-loadings out/decomp_clr/loadings.tsv \
  --output-activations out/decomp_clr/activations.tsv \
  --output-stability out/decomp_clr/stability.tsv
```

Transforms (atman's canonical abundance is already log-scale, so
these are linear operations on log values):

- `none` — pass-through.
- `clr` — per-sample mean centering of the log-abundance row.
- `alr` — subtract a reference protein's log-abundance column
  (`--alr-reference <gene_symbol>`, e.g. `ALB` for CSF).
- `ilr` — Helmert orthonormal projection; outputs are in
  `ilr_coord_*` coordinates (p − 1 dimensions), not raw proteins.
- `ratio-anchor` — algebraically identical to ALR on log input;
  preserved as a distinct tag for domain convention.

Acceptance:

- Pure transforms live in `atman-core::compositional`; 12 unit tests
  cover invariants (CLR rows sum to 0, ALR reference column is 0,
  ILR preserves CLR norm, constant-row maps to 0).
- Integration tests cover CLR, ALR (valid + missing reference +
  unknown gene), ILR dimensionality, and determinism across repeated
  runs.
- Run sidecar records `transform` and `alr-reference` alongside the
  existing ICA knobs.

Why it matters: closes the loop on a CSF-manuscript finding —
albumin normalization is not a preprocessing hack, it is a
compositional transform, and the decomposition result depends on
which transform you chose. No standalone proteomics ICA tool
currently exposes this as a first-class CLI option with provenance.

### 19. Bootstrap alignment uncertainty (`atman align bootstrap`) (implemented)

Subject-level bootstrap of cross-cohort archetype alignment.
Every iteration resamples subjects within each cohort with
replacement, re-runs FastICA per cohort, re-aligns under the
chosen `--metric`, and matches every point-estimate archetype to
its best-similarity bootstrap counterpart. Output summarises
per-archetype `bootstrap_prob_universal` /
`bootstrap_prob_multi`, both a percentile-CI band and a BCa CI on
the cohort-count distribution (alpha set by `--ci-alpha`,
default 0.05 ⇒ 95%), plus Shannon entropy (`alignment_entropy`)
over the same histogram.

Command:

```bash
atman align bootstrap \
  --cohorts out_cohort_a,out_cohort_b,out_cohort_c \
  --labels A,B,C \
  --k 20 --n-boot 200 --seed 20260418 \
  --metric cosine \
  --cosine-tau 0.30 --match-tau 0.50 \
  --ci-alpha 0.05 \
  --min-subjects 20 \
  --max-missing-fraction 0.5 --impute mean \
  --output out/align_bootstrap_summary.tsv
```

Scope:

- `--metric` accepts `cosine` (default), `jaccard` (top-N loading
  set overlap), or `spearman` (absolute rank correlation).
- BCa CI at the level set by `--ci-alpha`. Acceleration is
  estimated by pooled subject-level jackknife (one leave-one-out
  pass per subject across all cohorts). When the BCa denominator
  goes non-monotone the CI falls back to the percentile bound and
  the `bca_fallback_to_percentile` column flips to 1.
- Shannon entropy (`alignment_entropy`, bits) over the empirical
  histogram of `n_cohorts` across the bootstrap, counting misses
  as `n_cohorts = 0`. Low values ⇒ the archetype's cohort coverage
  is stable under resampling.
- Archetype matching is by representative-program best-cosine with
  a user-tunable `--match-tau` floor (default 0.50).

Acceptance:

- Pure bootstrap lives in `atman-core::align_bootstrap`; 11 unit
  tests cover row resampling, single-cohort refusal, min-subjects
  refusal, mismatched-universe refusal, finite-stats on toy data,
  determinism, Shannon-entropy edge cases, and BCa CI
  (symmetric-data percentile reduction, zero-spread jackknife
  fallback to bias correction, empty-bootstrap fallback).
- Integration tests cover CLI schema, single-cohort refusal,
  min-subjects refusal, determinism across repeated runs, and a
  **magnitude assertion**: on a strong planted fixture (40 subjects
  × 15 proteins, universal signal on 6 shared proteins) at least
  one archetype must exhibit `bootstrap_prob_multi ≥ 0.7`,
  `alignment_entropy < 1.0`, and BCa upper bound ≥ 2. This is the
  bootstrap-has-power contract — schema alone is not enough.
- Each cohort's canonical directory is intersected to a shared
  protein universe before alignment; `--impute mean` fills residual
  missingness within the `--max-missing-fraction` envelope.
- Run sidecar captures every knob plus per-cohort input hashes.

Why it matters: moves alignment from "point estimate at some
`--cosine-tau`" to "posterior probability the archetype is
cross-cohort universal." Directly answers the reviewer question
"would this survive a different random split of your data?"

### 20. Archetype variance decomposition (`atman decompose variance`) (implemented, v1)

Mixed-model variance partition per archetype. Takes a
`decompose ica` activations TSV plus a samples metadata file, fits
`y_archetype ~ fixed_terms + (1|group)` per archetype via the
existing REML engine, and reports per-fixed-factor projection
variance, random-intercept variance, residual variance, and the
intraclass correlation for the random factor.

Command:

```bash
atman decompose variance \
  --activations out/activations.tsv \
  --samples out/samples.tsv \
  --factors "cohort + condition + (1|subject_id)" \
  --output out/archetype_variance.tsv
```

v1 scope:

- Type-I projection variance (`var_f = Var(X_f · β_f)` across
  samples). Correlation-adjusted Type II / III partitions deferred.
- Per-coefficient Wald summaries (`max_abs_t_<factor>`,
  `min_p_<factor>`) rather than factor-level omnibus F-tests.
- At most one random-intercept term `(1|col)` per run.

Acceptance:

- Pure partition math lives in
  `atman-core::variance_decomposition`; 6 unit tests cover sample
  variance math, factor projection math, recovery of fixed-only
  structure, recovery of random-intercept variance on grouped data,
  skip-on-singular-design propagation, and determinism.
- Integration tests exercise the full CLI path on a synthetic
  2-cohort × 2-condition fixture: `cohort_confounded` archetype
  gets `var_cohort > var_condition > var_residual`;
  `condition_specific` gets the reverse; `cohort_shared` noise gets
  `var_residual` dominant. Random-intercept grouping on `cohort`
  gives `cohort_confounded` ICC > 0.3.
- Single-level-factor refusal with a clear error.

Why it matters: gives a quantitative answer to "what fraction of
this archetype is shared biology vs. between-cohort batch?" — the
hardest-to-defend question in any cross-cohort factor paper.

### 21. Multi-level omnibus F-test (`atman de --test ols --omnibus-factor`) (implemented, v1)

Per-protein joint hypothesis test on a categorical factor's
coefficients. Answers "does this protein respond to the factor at
all?" — the omnibus question that precedes any specific pairwise
contrast. Directly addresses the common workflow where users have
≥ 3 levels of a categorical (`stage`, `diagnosis`, `dose`) and
currently stitch pairwise DEs together by hand, losing omnibus
information and FWER control.

Command:

```bash
atman de --test ols \
  --design "~ condition + stage" \
  --contrast conditionCase \
  --omnibus-factor stage \
  --groups "Case-Control" \
  --output-dir out_omnibus
```

Output `de_omnibus.tsv`:

- `panel, assay_id, gene_symbol, factor, comparison,
  f_statistic, df_num, df_den, p_value, bh_q`
- `bh_q` is BH-adjusted within each `(comparison, panel)` family.
- `df_num` = number of factor columns (levels − 1);
  `df_den` = residual degrees of freedom from the OLS fit.

Acceptance:

- Unit tests verify `contrast_inference` matches OLS's own β / SE /
  t / p on a single-coefficient contrast to 1e-10 and correctly
  estimates `β_level_i − β_level_j` for a difference contrast.
- Unit test verifies `omnibus_f_test` produces F > 10 / p < 0.001
  on a planted 3-level signal and p > 0.05 on deterministic-noise
  null data.
- Integration test on a 3-stage 2-condition synthetic cohort shows
  RESPONDER F≈21, NONRESPONDER F≈0.2, BH-q separation at 0.05 to
  the expected side in both cases.
- Refuses when combined with `--test welch-t`, `--test paired-t`,
  or any test other than `ols`.
- Refuses when no `--design` is supplied.

v1 scope trim (**superseded — see below**):

- No automatic Tukey HSD or Dunnett pairwise contrasts (requires
  studentized-range and multivariate-t distributions not in
  statrs). User-supplied `--contrast-list` with Sidak adjustment
  is the next increment on the post-hoc track.

**Status correction, 2026-09-06.** That trim no longer describes the
code. All three post-hoc methods ship and are dispatched today:
`--post-hoc sidak`, `--post-hoc tukey` (studentized-range) and
`--post-hoc dunnett` (Dunnett–Hsu, with the unbalanced-design
correlation path). The distributions the trim called missing were
added in-tree as `atman_core::studentized_range` and
`atman_core::multivariate_t`; see the CHANGELOG entries and
`atman de --help`. Check capability claims against `--help` and the
CHANGELOG rather than against this section.

### 22. Feature-covariance network influence (`atman network influence`) (implemented)

New top-level command family (`atman network …`) for
feature-graph analysis. `influence` builds the subject-level
covariance graph over proteins (or assays via `--feature-col
assay_id`) under pearson / spearman / covariance similarity,
applies a hard-threshold or WGCNA-style soft-power adjacency
policy, and scores every feature by the Burberry-Pillai
construction `influence = eigenvector × betweenness`.

Command:

```bash
atman network influence \
  --input out/measurements.tsv \
  --samples out/samples.tsv \
  --method spearman \
  --threshold 0.3 \
  --stratify condition \
  --output out/network_influence.tsv
```

Pure math (`atman-core::network`):

- Similarity matrices from the shared `stats` primitives.
- `HardThreshold { threshold }` or `SoftPower { power }` adjacency
  construction.
- Eigenvector centrality via shifted power iteration
  (`A + αI`, α = max row sum) so bipartite graphs don't oscillate.
- Brandes' weighted betweenness with Dijkstra SSSP on reciprocal
  similarities (high-similarity = short path).
- `influence_score = eigenvector × betweenness`.

Outputs:

- `network_influence.tsv`: `feature_id, stratum,
  eigenvector_centrality, betweenness_centrality, influence_score,
  degree, n_subjects_used`.
- `network_edges.tsv` (when `--emit-edges`): every retained edge.
- Run sidecar captures method, adjacency policy, stratify column,
  eigen iteration params.

Acceptance:

- 8 unit tests cover star-topology hub recovery for both
  eigenvector and betweenness, clique symmetry, hard-threshold
  drop, soft-power scaling, zero-adjacency produces zero
  centralities, and determinism.
- 5 integration tests cover output schema, empty-graph refusal,
  stratify disjointness, repeated-run byte equality, and
  conflicting-adjacency-flag refusal.
- Smoke-tested on Dube (2938 features, ~600 degree per hub).

Why it matters: complements `align programs` by answering the
orthogonal question of which individual features scaffold the
covariance structure archetypes are built on. No standalone
proteomics tool currently ships feature-hub scoring as a
first-class CLI with provenance.

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
17. Archetype null calibration
18. Compositional transforms in decompose ica
19. Bootstrap alignment uncertainty
20. Archetype variance decomposition
21. Multi-level omnibus F-test
22. Feature-covariance network influence

## Release Readiness Status

Atman now covers the standalone release surface targeted by this roadmap:

- Canonical validation and QC reporting are implemented.
- Olink and matrix ingest paths are covered by integration tests and adapter
  fixtures.
- Differential abundance, uncertainty, robustness, module, enrichment, and
  meta-analysis commands are implemented.
- The core CLI does not require Python, notebooks, services, or external
  workflow infrastructure.

## Open question: does a bulk proteome have module structure at all?

Raised 2026-09-06 while fixing `modules discover`'s soft-power fallback.
On a 93-sample, 9,376-feature CPTAC GBM contrast, **no** soft power in
1..=20 satisfied the scale-free criterion: R² fell monotonically from
0.505 at β = 1 to 0.133 at β = 20. The fix changed what atman does when
that happens (WGCNA's default power by sample count, and a refusal to
emit an all-grey discovery), but it does not answer the underlying
question:

> Is the scale-free criterion simply the wrong topology model for a bulk
> proteome correlation matrix, or does that data genuinely lack module
> structure?

These are very different conclusions and they call for different tools.
The discriminating experiment is cheap: run the same data at a
conventional β (6 for this sample count) and count non-grey modules,
**with their size distribution** — a handful of tiny modules beside a
large grey remainder means something different from a genuine partition.

The atman methods-paper session was running exactly that when both
sessions wrapped on 2026-09-06; the result is expected in
`atman-paper/runs/determinism/`, alongside their soft-power sweep
evidence. Pick it up from there before designing anything new here. If
real modules do emerge at β = 6, the scale-free sweep is the wrong
criterion for this data type and that is worth stating in the docs. If
they do not, `modules discover` is the wrong instrument for bulk
proteomes and no amount of parameter tuning fixes that.

## Unscheduled: alternative decompositions

Method-transfer ideas, parked here rather than in `feature_requests.md`
because nobody has requested them. They were proposed on 2026-04-20 while
adjacent decomposition work was underway, and the 2026-09-06 survey of
every manuscript session found no project that wants one. They are not
queued; they are a menu to reach for if a future project needs an
alternative to the FastICA + ICASSO path the current work commits to.

Anything built from here starts by re-asking whether a project needs it,
and then goes through the normal spec route, not straight to code.

### PMF (positive matrix factorization) with per-cell uncertainty weighting

Port from atmospheric source apportionment (Paatero & Tapper 1994).
Sibling to `atman decompose unmix`: where VCA+FCLS finds geometric
endmembers, PMF finds non-negative factors while natively weighting each
cell by its reciprocal analytical uncertainty. Would suit DIA-MS inputs
that ship per-protein-per-sample CV or detection-limit flags. Plausible
surface: `atman decompose pmf --uncertainty measurements_cv.tsv`.

### MCR-ALS (multivariate curve resolution — alternating least squares)

From chemometrics. Sibling to `atman decompose unmix`: NMF with
chemistry-aware constraints (non-negativity, unimodality, closure) on the
protein-by-sample matrix. Plausible surface: `atman decompose mcr`.

### Longitudinal tensor decomposition

Samples × proteins × timepoints, where random-subject variance becomes
the signal rather than the nuisance. Would sit on top of the existing
`decompose variance` surface. Note the practical obstacle: the only
longitudinal arm in any current project is 12 subjects × 3 timepoints,
which the existing one-way ICC already covers.

## Dependency notes

### `statrs` (0.17) — symbol-usage audit

Audited 2026-05-29 (TASK-014). `statrs` is the workspace's distribution-CDF
provider; every use is a closed-form CDF/quantile or special function, never
matrix linear algebra and never random sampling.

Symbols actually imported across `atman` + `atman-core`:

- `distribution::ContinuousCDF` (trait) — `.cdf()` on the below distributions.
- `distribution::DiscreteCDF` (trait) — `.cdf()` for the discrete case.
- `distribution::Normal` — z-tests, normal-approx p-values
  (`align_bootstrap`, `ensemble`, `network_differential`, `studentized_range`,
  `multivariate_t`, `coupling`, `meta`, `ratio`, `detectability`, `robust_stats`,
  and the normal fallback path in `de`).
- `distribution::StudentsT` — moderated/empirical-Bayes t-tests
  (`limma`, `msqrob`, `de`, `moderated`, `robust_stats`, `ratio`, `null`).
- `distribution::FisherSnedecor` — F-test tail probabilities (`limma`, `de`).
- `distribution::Hypergeometric` — ORA / over-representation exact tail
  (`decompose_unmix`).
- `function::gamma::ln_gamma` — log-gamma for the studentized-range and
  multivariate-t densities (`studentized_range`, `multivariate_t`).

Transitive footprint and whether it is justified:

`statrs` pulls `nalgebra` (+ `nalgebra-macros`), `rand 0.8.5`, `approx`, and
`num-traits`. The workspace uses **none** of the features those transitive crates
exist to serve: no `MultivariateNormal`/`Dirichlet` (the only `statrs`
distributions backed by `nalgebra`), and no sampling (`rand` is `statrs`'s
`Distribution::sample` path). So the `nalgebra` + `rand` weight is paid for
features we do not use — it is *not* justified by current usage, only tolerated.

Decision: keep `statrs` as-is for this task (audit only). The scalar CDFs and
`ln_gamma` we depend on are correct, well-tested, and match published reference
values, and re-implementing them deterministically would re-litigate numerical
accuracy for no provenance gain. The transitive bloat is a candidate for a
future task — either a Cargo feature gate that drops `nalgebra`/`rand` if
`statrs` exposes one, or replacing the handful of CDFs with small in-tree
implementations (the `de.rs` normal-fallback comment already hints the team has
considered avoiding `StudentsT` allocation in hot paths). Do not act on that here.
