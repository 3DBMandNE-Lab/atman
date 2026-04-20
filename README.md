# Atman

[![CI](https://github.com/kevinj24fr/atman/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/kevinj24fr/atman/actions/workflows/ci.yml)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](./LICENSE-MIT)

Atman is a standalone Rust command-line tool for reproducible, deterministic
proteomics workflows across affinity and mass-spectrometry platforms. It ships
wide-matrix and platform-native ingest paths, a canonical TSV interchange
format, and analysis commands for differential abundance, decomposition,
module testing, enrichment, asymmetry, and robustness summaries.

Supported platforms: Olink Explore NGS, SomaScan, MaxQuant/LFQ, DIA-NN,
Spectronaut, and any wide abundance matrix via `ingest-matrix` or a
thin user-written adapter emitting the canonical TSV schema.

Atman does not require a service, database, notebook runtime, or workflow
system. The repository includes a Dube et al. 2023 Olink Explore fixture
dataset so the package can test its raw-to-result reproduction path end to end,
and `adapters/examples/` ships tiny synthetic fixtures for DIA-NN, MaxQuant,
and SomaScan-style matrices.

## Install

From source with stable Rust 1.94 or newer:

```bash
git clone https://github.com/kevinj24fr/atman.git
cd atman
cargo install --path crates/atman
```

Or build a container:

```bash
docker build -t atman:1.0.0 .
docker run --rm atman:1.0.0 --help
```

## Validation

Atman's statistical tests are validated to numerical parity against
canonical R references on published fixtures. Every reference file below
ships in the repository; every reference script regenerates the reference
TSV deterministically under a pinned R + Bioconductor version.

| Test path | R reference | Fixture | Agreement | atman wall-clock |
|---|---|---|---|---|
| `de --test limma --peptide-metadata` (DEqMS trend) | Bioconductor `DEqMS::spectraCounteBayes` v1.26 | CPTAC Study 6 UPS1 spike-in (29-protein subset) | median per-protein `|atman − DEqMS|` = **0.000 log₂** (bit-for-bit match at 3 d.p. on every displayed protein) | 0.60 s |
| `de --test msqrob` | `msqrob2::msqrob(~condition)` v1.16 via the canonical vignette workflow | CPTAC Study 6 UPS1 (30-protein subset from `statOmics/MSqRobData`) | median per-protein `|atman − msqrob2|` = **0.118 log₂** across 26 jointly-fitted proteins; **100% sign agreement** on all 11 signal proteins `\|effect\| ≥ 0.3 log₂`; 9/9 UPS1 proteins recover the expected negative direction | 0.01 s |
| `de --test limma` (parametric eBayes, `trend=false`, `robust=false`) | `limma::eBayes` 3.x | 100-feature × 20-sample regression fixture | `\|atman − limma\|` < **1e-4** on `t`, `p_value`, `df_total`, `s2_post` | 0.06 s |
| `de --post-hoc sidak --test ols` | `R lm()` + `pairwise.t.test(p.adjust = "bonferroni")` with Sidak re-derivation | Planted 3-level stage × 2-covariate fixture | `max abs diff` < **1e-6** on `estimate`, `posthoc_p`, `posthoc_adj_p` | 0.01 s |
| `decompose variance --omnibus-factor` (Type III F) | `car::Anova(type = 3)` (Wald form) | Planted variance-attribution fixture | `max abs p diff` < **1e-6** | < 0.01 s |
| Olink Explore NPX reproduction (`ingest`, `qc`, `matrix`, `fold-change`) | Dube et al. 2023 Scientific Data published tables | Dube heat-stress cohort (2 panels, 8 filtered NPX files, 8 log2-FC files) | filtered NPX files: **byte-exact** cell match; log2-FC files: match to **floating-point precision** (max delta 1.05e-15) | — |

Wall-clock numbers are end-to-end test times on the release build
(Apple M-series, single thread), including fixture I/O and binary
invocation. They are upper bounds on the analytical compute cost per
fixture.

Reference scripts live alongside the fixtures in
`crates/atman/tests/fixtures/`:

- `deqms_cptac_reference.R`
- `msqrob2_cptac_reference.R`
- `gen_limma_reference.R`
- `posthoc_sidak_reference.R`
- `variance_type3_reference.R`

The parity assertions above run on every `cargo test --workspace --release`
invocation; CI fails if any drifts.

## Input Support

Atman has three input modes, all landing on the same canonical TSV schema:

- **Matrix ingest:** `atman ingest-matrix` reads common wide protein matrices
  (DIA-NN, Spectronaut, MaxQuant/LFQ, SomaScan, any log2 intensity/abundance
  export) plus sample/protein metadata and writes Atman's canonical TSVs.
  Supports `--log2-transform` for linear-scale inputs and both
  `proteins-rows` / `samples-rows` orientations.
- **Platform-native ingest:** `atman ingest --platform olink-explore-ngs`
  reads raw long-format Olink Explore NPX CSV files. Additional native
  parsers can be added per platform as needed.
- **Canonical TSV adapters:** any script or converter can write
  `samples.tsv`, `proteins.tsv`, `measurements.tsv`, and optionally
  `qc_measurements.tsv`. Once those files exist, Atman's downstream commands
  (`de`, `module-de`, `robustness`, and related summaries) run the same way
  regardless of the original assay source.

Adapters are intentionally thin: normalize source metadata, map samples and
proteins, log-transform linear intensities when needed, and emit Atman's TSV
schema. See `adapters/examples/` for tiny synthetic DIA-NN, MaxQuant, and
SomaScan-style fixtures end-to-end.

## Platform Quirks: Adapter Responsibilities

Atman's canonical TSV absorbs scale (`abundance_unit` + `--log2-transform`),
per-sample/per-assay QC flags, LOD (`below_lod`, `detection_limit`), panel/plate
provenance, and peptide-level inputs for `msqrob`/`limma-DEqMS`. Analysis
commands are agnostic to those.

Analysis commands are **not** automatic on three fronts. The adapter (or your
upstream pipeline) is responsible:

1. **Cross-sample normalization.** `atman ingest-matrix --normalize
   median|quantile|none` applies per-sample median centering or quantile
   normalization on log-scale values. Default is `none`: Olink NPX arrives
   pre-normalized, so pass-through is correct there. When `--normalize none`
   is used and per-sample median abundance spans more than one log2 unit,
   Atman emits a stderr warning recommending `median`. For SomaScan RFU,
   MaxQuant/LFQ, DIA-NN, and Spectronaut exports, pass `--normalize median`
   (or normalize upstream). VSN, plate bridging, and ComBat-style correction
   are out of scope — do those upstream.
2. **Missingness semantics (MNAR vs MCAR).** Atman drops missing cells per
   protein/test. `below_lod` is recorded but the default tests still treat
   it as MCAR. `atman de` refuses to run when the fraction of `below_lod=1`
   rows among non-dropped measurements exceeds `--max-below-lod-fraction`
   (default 0.5). Either impute upstream (e.g. Perseus-style min-shifted),
   filter proteins to an acceptable observation rate, or pass
   `--allow-censored` after confirming the chosen test is appropriate.
   `decompose ica` exposes `--impute mean|none`; the other analysis
   commands do not.
3. **Protein-group ambiguity.** `proteins.tsv` is one `assay_id` per row.
   MaxQuant `P1;P2;P3` protein-group rows and similar ambiguous identifiers
   must be picked/split/deduped by the adapter. Atman will not disambiguate.

For plate or batch effects, pass `batch` as a covariate to `atman de --test ols
--design "~ condition + batch"`. There is no ComBat-style in-place correction.

## Commands

```text
atman ingest-matrix      wide protein matrix + metadata -> canonical TSVs
                         (DIA-NN, Spectronaut, MaxQuant/LFQ, SomaScan, ...)
atman ingest             platform-native parser (Olink Explore NPX)
                         -> canonical TSVs
atman align programs     cross-cohort program alignment + sensitivity sweep
atman align bootstrap    subject-level bootstrap of cross-cohort archetype alignment
atman decompose ica      multi-seed FastICA with seed-stability reporting
                         (supports --transform clr/alr/ilr/ratio-anchor)
atman decompose null     permutation-null calibration of archetype stability
atman decompose variance mixed-model variance partition per archetype
atman validate           check canonical TSV schema, keys, and sample support
atman qc                 apply QC masking rules
atman report qc          summarize QC, missingness, and condition support
atman matrix             canonical long TSV -> per-panel wide NPX CSVs
atman fold-change        compute per-panel log2 fold-change tables
atman de                 paired, Welch, OLS, mixed, limma, msqrob, or cross-method ensemble DE
atman bootstrap protein  subject-level bootstrap intervals for protein effects
atman bootstrap module   subject-level bootstrap intervals for module effects
atman null               permutation/sign-flip null calibration for DE effects
atman asymmetry          compare matched contrast pairs
atman robustness         summarize rerun/LOO rank and sign stability
atman module-trajectory  score user-defined modules from per-subject deltas
atman module-de          aggregate proteins into modules and test at module level
atman score modules      score modules per sample from canonical measurements
atman enrich ora         over-representation analysis from DE hits
atman enrich gprofiler   live g:Profiler REST wrapper with cached responses
atman meta               combine DE results across cohorts
atman network influence  feature-covariance hub scoring (eigenvector × betweenness)
atman run                execute a plan YAML/JSON and emit a hash manifest
```

## Quick Start

Pick one ingest path. The rest of the pipeline is platform-agnostic.

```bash
# Option A — wide-matrix ingest (DIA-NN, Spectronaut, MaxQuant/LFQ, SomaScan, ...)
atman ingest-matrix \
    --matrix matrix.tsv --samples sample_metadata.tsv \
    --orientation proteins-rows \
    --platform diann_report \
    --abundance-unit log2_diann_pg_quantity \
    --assay-id-col Protein.Group --gene-col Genes \
    --condition-col diagnosis \
    --output-dir out

# Option B — platform-native Olink Explore NPX ingest (bundled fixture)
atman ingest \
    --platform olink-explore-ngs --parser dube \
    --output-dir out \
    example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
    example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv

atman qc --input-dir out --output-dir out --rule dube

atman validate \
    --input-dir out \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" \
    --min-pairs 5

atman report qc \
    --input-dir out \
    --output-dir out/report

atman matrix \
    --input-dir out --output-dir out \
    --format dube-wide --split-by panel

atman fold-change \
    --input-dir out --output-dir out \
    --groups "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1"

atman de \
    --input-dir out --output-dir out \
    --test paired-t --paired-by participant \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" \
    --min-pairs 5

# Covariate-adjusted OLS with formula-style design.
atman de \
    --input-dir out --output-dir out_ols \
    --test ols \
    --groups "Case-Control" \
    --design "~ condition + age + sex + batch" \
    --contrast conditionCase \
    --min-pairs 5

# Multi-level omnibus F-test: report a per-protein F for whether a
# categorical factor (here: stage with 3 levels CN/MCI/AD) has ANY
# effect, alongside the standard pairwise contrast. Emits
# de_omnibus.tsv with BH-adjusted q within each (comparison, panel).
atman de \
    --input-dir out --output-dir out_ols_omnibus \
    --test ols \
    --groups "Case-Control" \
    --design "~ condition + stage" \
    --contrast conditionCase \
    --omnibus-factor stage \
    --min-pairs 5

# Repeated-measures random-intercept model.
atman de \
    --input-dir out --output-dir out_mixed \
    --test mixed \
    --groups "PT2-PT1" \
    --fixed "condition + age + sex" \
    --random "1|subject_id" \
    --min-pairs 5

# limma-grade eBayes with parametric mean-variance trend and robust prior.
atman de \
    --input-dir out --output-dir out_limma \
    --test limma \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" \
    --paired-by participant \
    --trend true --robust true \
    --min-pairs 5

# DEqMS peptide-count-weighted eBayes. Swaps the mean-variance trend
# covariate for log(peptide_count + 1) via a tricube kernel smoother
# (Zhu et al. 2020). Requires peptides.tsv for the peptide→protein
# mapping; any protein missing from peptides.tsv gets count 0.
atman de \
    --input-dir out --output-dir out_deqms \
    --test limma \
    --peptide-metadata out/peptides.tsv \
    --groups "Case-Control" \
    --trend true \
    --min-pairs 5

# Cross-method consensus DE. Runs every listed method on the same
# canonical inputs, Stouffer-combines per-method p-values within each
# (comparison, protein), BH-FDRs across proteins, and assigns a
# VALIDATED / PROVISIONAL / INSUFFICIENT grade based on ensemble_q +
# sign consistency. Methods with missing inputs (msqrob without
# --peptide-measurements, etc.) are auto-skipped, not an error.
#
# INTERPRETIVE CAVEAT. Ensemble is a within-dataset robustness check,
# not independent-study meta-analysis. paired-t, welch-t, ols, mixed,
# limma, and msqrob all share most of their signal on the same
# measurement matrix, so Stouffer p-values are not the meta-analytic
# combination of independent studies. Read VALIDATED as "the finding
# survives method swap on this dataset," not "independently replicated."
atman de \
    --input-dir out --output-dir out_ensemble \
    --test ensemble \
    --ensemble-methods "paired-t,welch-t,limma,msqrob" \
    --peptide-measurements out/peptide_measurements.tsv \
    --peptide-metadata out/peptides.tsv \
    --groups "PT2-PR2" \
    --paired-by participant \
    --min-pairs 5

# TREAT (minimum-effect test) at log2-FC threshold 0.5.
atman de \
    --input-dir out --output-dir out_limma_treat \
    --test limma \
    --groups "PT1-PR1" \
    --lfc-threshold 0.5 \
    --min-pairs 5

# msqrob peptide-level ridge mixed model. One model per protein over its
# peptides, random intercept per peptide, ridge penalty on non-intercept
# fixed effects. Requires peptide_measurements.tsv + peptides.tsv.
atman de \
    --input-dir out --output-dir out_msqrob \
    --test msqrob \
    --peptide-measurements out/peptide_measurements.tsv \
    --peptide-metadata out/peptides.tsv \
    --groups "Case-Control" \
    --ridge-lambda 0.5 \
    --min-peptides 2 \
    --min-pairs 5

atman bootstrap protein \
    --input-dir out \
    --groups "PT2-PT1" \
    --test paired-t \
    --n 2000 \
    --seed 1 \
    --output out/protein_bootstrap.tsv

atman bootstrap module \
    --input-dir out \
    --modules-tsv modules.tsv \
    --groups "PT2-PT1" \
    --test paired-t \
    --n 2000 \
    --seed 1 \
    --output out/module_bootstrap.tsv

atman null \
    --input-dir out \
    --output-dir out/null \
    --groups "PT2-PT1" \
    --test paired-t \
    --n 1000 \
    --seed 1

atman score modules \
    --input-dir out \
    --modules-tsv modules.tsv \
    --method mean \
    --output out/module_scores.tsv \
    --canonical-output-dir out/module_score_canonical

atman enrich ora \
    --de-results out/de_results.tsv \
    --gene-sets gene_sets.tsv \
    --comparison "PT2-PT1" \
    --output out/ora.tsv

atman meta \
    --inputs cohort1/de_results.tsv,cohort2/de_results.tsv \
    --output meta.tsv
```

Common outputs:

- `measurements.tsv`, `qc_measurements.tsv`, `samples.tsv`, `proteins.tsv`
- `peptide_measurements.tsv`, `peptides.tsv` (msqrob peptide-level input)
- optional `validate_report.tsv`
- `report/qc_summary.tsv`, `report/sample_qc.tsv`, `report/protein_qc.tsv`,
  `report/condition_counts.tsv`
- `<panel>_npx.csv`
- `<panel>_log2_fc.csv`
- `de_results.tsv`, `de_report.tsv`; OLS also writes `de_covariates.tsv` and
  `de_design.tsv`
- `de_results.tsv` appends effect-size, confidence-interval, Wilcoxon, median,
  and trimmed-mean columns after the original DE fields
- `protein_bootstrap.tsv`, `module_bootstrap.tsv`
- `null/null_summary.tsv`, `null/empirical_p.tsv`
- `module_scores.tsv`; optional module-score canonical TSVs for downstream DE
- `ora.tsv`
- `meta.tsv`
- `de_ensemble.tsv` (ensemble mode): per-(comparison, protein)
  cross-method agreement with `n_applied`, `n_significant`,
  `n_sign_consistent`, `majority_sign`, `ensemble_p` (Stouffer),
  `ensemble_q` (BH within comparison), `grade`
  (VALIDATED / PROVISIONAL / INSUFFICIENT), `methods_applied`,
  `methods_skipped`. `de_results.tsv` in ensemble mode additionally
  tags every row with the emitting `method`.

## Reproducibility Check

Run the package test suite:

```bash
cargo test --workspace --release
```

The integration tests execute `ingest -> qc -> matrix -> fold-change` against
the bundled Dube fixture data and compare outputs with the published reference
tables:

- filtered NPX files: byte-exact cell match
- log2 fold-change files: numerical match to floating-point precision

The DE sanity test also verifies that canonical heat-shock proteins increase
in both acute heat comparisons.

### Run sidecars

Every analysis subcommand that writes output files also writes a
`<primary_output>.run.json` sidecar next to its main artifact (e.g.
`loadings.tsv.run.json` alongside `loadings.tsv`). The sidecar captures:

- the atman binary version and git SHA,
- the fully-resolved argument dict (including defaults),
- a SHA-256 hash of the canonical input TSVs the command actually read,
- SHA-256 hashes of every output file from this invocation,
- ISO-8601 UTC start/finish timestamps and build target triple.

Reviewers can reproduce a run from the sidecar alone — no need to chase
the driver script. Commands with sidecar output: `ingest-matrix`,
`validate` (when `--report` is set), `report qc`, `decompose ica`,
`align programs`, `coupling`, `null`, `enrich gprofiler`, `de`, `bootstrap
protein`, `bootstrap module`, `bootstrap program`, `meta`, and `ratio`.
The plan-level `atman run` manifest already covers its own provenance.

## Data Model

Atman writes a small canonical TSV dataset between commands. This schema is
the adapter target for any proteomics source:

- `samples.tsv`: sample IDs, subject IDs, conditions, controls, ingest order
- `proteins.tsv`: assay IDs, UniProt IDs, gene symbols, panel metadata
- `measurements.tsv`: one row per sample-assay abundance measurement
- `qc_measurements.tsv`: same schema after QC masking
- `peptides.tsv`: peptide catalog mapping `peptide_id → assay_id` (parent
  protein) plus optional sequence, charge, modifications,
  missed_cleavages. Input to `atman de --test msqrob`.
- `peptide_measurements.tsv`: one row per sample × peptide abundance.
  Schema: `sample_id, peptide_id, abundance, abundance_unit,
  dropped_by_qc, below_lod`.

Downstream commands operate on these files, so each stage can be inspected,
rerun, or replaced independently.

## Layout

```text
crates/atman-core/       core data model and algorithms
crates/atman/            CLI, command orchestration, and file IO
adapters/                canonical TSV adapter helpers and templates
adapters/examples/       tiny synthetic matrix examples (DIA-NN, MaxQuant, SomaScan)
CITATION.cff             citation metadata for release archives
docs/tutorial.md         package tutorial using the bundled Dube fixture
docs/release-checklist.md
                         standalone release checklist
docs/analytical-roadmap.md
                         implemented analytical capability roadmap
example_data/dube_heat_2023/
                         Dube et al. 2023 Olink Explore fixture data
```

## Reference Dataset

Gagnon D, Barry H, Barhdadi A, Oussaid E, Mongrain I, Lemieux Perreault LP,
Dubé MP. *A dataset of proteomic changes during human heat stress and heat
acclimation.* Scientific Data (2023).
[10.1038/s41597-023-02809-5](https://doi.org/10.1038/s41597-023-02809-5).

## License

Atman is licensed under either [MIT](./LICENSE-MIT) or
[Apache-2.0](./LICENSE-APACHE), at your option.
