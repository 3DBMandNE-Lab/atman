# Atman

[![CI](https://github.com/kevinj24fr/atman/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/kevinj24fr/atman/actions/workflows/ci.yml)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](./LICENSE-MIT)

Atman is a Rust command-line tool for proteomics workflows across
affinity and mass-spectrometry platforms. It ships wide-matrix and
platform-native ingest paths, a canonical TSV interchange format, and
analysis commands for differential abundance, decomposition, module
testing, enrichment, asymmetry, and robustness summaries.

Atman builds on the reference statistics (limma, DEqMS, msqrob2,
`car::Anova`, and base `lm + pairwise.t.test`) and validates each
path to numerical parity against those references on published
fixtures — see [Validation](#validation). Every invocation writes a
SHA-256 hash manifest (`*.run.json`) next to its primary output so a
run's resolved args, input hashes, output hashes, atman version, and
git SHA stay attached to the artifact.

**Data flow**

```
Input                  Atman stage                 Output
-------------------------------------------------------------------
raw NPX / matrix  ──▶  ingest / ingest-matrix  ──▶ canonical TSV
canonical TSV     ──▶  qc / validate / report  ──▶ QC'd TSV + report
QC'd TSV          ──▶  de / decompose / bench  ──▶ results.tsv
any stage         ──▶  (every command)         ──▶ <output>.run.json
                                                   (SHA-256 manifest)
```

Supported platforms: Olink Explore NGS, SomaScan, MaxQuant/LFQ, DIA-NN,
Spectronaut, and any wide abundance matrix via `ingest-matrix` or a
thin user-written adapter emitting the canonical TSV schema.

The repository ships the Dube et al. 2023 Olink Explore fixture so the
raw-to-result reproduction path is tested end-to-end;
`adapters/examples/` ships tiny synthetic fixtures for DIA-NN,
MaxQuant, and SomaScan-style matrices.

## Install

Install from source or container. Atman is not yet on crates.io.

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

Atman's tests are validated against R reference implementations on
published fixtures. Every reference file ships in the repository and
every reference script regenerates the reference TSV deterministically.

In the table, `Δ` means the element-wise absolute difference
`abs(atman − R)` between atman's output column and the corresponding R
column on the same protein/contrast.

| Test path | R reference | Fixture | Agreement | Wall-clock |
|---|---|---|---|---|
| `de --test limma --peptide-metadata` (DEqMS trend) | `DEqMS::spectraCounteBayes` | CPTAC Study 6 UPS1 spike-in (29-protein subset) | median per-protein `Δ` = **0.000 log₂** (three-decimal match on every displayed protein) | 0.60 s |
| `de --test msqrob` | `msqrob2::msqrob(~condition)` via the QFeatures vignette workflow (log2 → median centering → median peptide summary → protein-level fit) | CPTAC Study 6 UPS1 (30-protein subset from `statOmics/MSqRobData`) | median per-protein `Δ` = **0.118 log₂** across 26 jointly-fitted proteins; **100% sign agreement** on all 11 signal proteins (`abs(effect) ≥ 0.3 log₂`); 9/9 UPS1 proteins recover the expected negative direction. See note below — atman and msqrob2 are different estimators on this comparison. | 0.01 s |
| `de --test limma` (parametric eBayes, `trend=false`, `robust=false`) | `limma::eBayes` | 100-feature × 20-sample regression fixture | `Δ` < **1e-4** on `t`, `p_value`, `df_total`, `s2_post` | 0.06 s |
| `de --post-hoc sidak --test ols` | `lm()` + `pairwise.t.test(p.adjust = "bonferroni")` with Sidak re-derivation | Planted 3-level `stage` × 2-covariate fixture | `max Δ` < **1e-6** on `estimate`, `posthoc_p`, `posthoc_adj_p` | 0.01 s |
| `decompose variance --omnibus-factor` (Type III F) | `car::Anova(type = 3)` Wald form | Planted variance-attribution fixture | `max Δ p` < **1e-6** | < 0.01 s |
| Olink Explore NPX reproduction (`ingest`, `qc`, `matrix`, `fold-change`) | Dube et al. 2023 Scientific Data published tables | Dube heat-stress cohort (2 panels, 8 filtered NPX files, 8 log2-FC files) | filtered NPX files: **byte-exact** cell match; log2-FC files: **floating-point precision** (max `Δ` = 1.05e-15) | — |

Wall-clock numbers are end-to-end integration-test times on the
release build, single-threaded, Apple M-series. They include fixture
I/O and binary invocation and are upper bounds on the analytical
compute cost per fixture.

### Note on msqrob2 parity

Atman's `--test msqrob` and msqrob2 in the canonical QFeatures vignette
workflow are not the same estimator; the 0.118 log₂ median difference
reflects the method choice rather than a numerical gap:

- **msqrob2 here (aggregate-then-fit).** Peptides are summarised to one
  value per (protein, sample) via `matrixStats::colMedians`, then
  `msqrob(~condition)` fits per-protein OLS with empirical-Bayes
  variance shrinkage.
- **Atman (joint fit).** One peptide-level linear mixed model per
  protein, random intercept per peptide, ridge penalty on non-intercept
  fixed effects, REML 1D profile over `τ = σ²_peptide / σ²_res`, with
  empirical-Bayes variance shrinkage across all fitted proteins.

The substantive claims are sign agreement and spike-in recovery; the
median-Δ is reported for transparency. The aggregate-then-fit path can
be reproduced inside atman by first summarising peptides to proteins
with `atman score modules --method median` and running
`atman de --test ols`, at which point the two paths converge.

### Reference environment

Reference scripts live alongside the fixtures in
`crates/atman/tests/fixtures/` and declare their package requirements
in the script header:

| Script | Packages |
|---|---|
| `deqms_cptac_reference.R` | limma (CRAN), DEqMS (Bioconductor ≥ 1.26), matrixStats |
| `msqrob2_cptac_reference.R` | msqrob2 (Bioconductor ≥ 1.16), QFeatures, SummarizedExperiment |
| `gen_limma_reference.R` | limma (CRAN) 3.x |
| `posthoc_sidak_reference.R` | base R (`lm`, `pairwise.t.test`) |
| `variance_type3_reference.R` | car (CRAN) |

Reference TSVs committed in the repo were generated under
R 4.3.x + Bioconductor 3.18. Pinning to an exact environment via
`renv.lock` or a Dockerfile is a later pass; in the meantime the
committed reference TSVs are what the parity tests diff against, so
drift on the R side cannot affect CI.

The parity assertions run on every `cargo test --workspace --release`;
CI fails if any drifts.

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

Minimal end-to-end: ingest → validate → QC → paired-t DE. For the
bundled Dube fixture:

```bash
atman ingest \
    --platform olink-explore-ngs --parser dube \
    --output-dir out \
    example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
    example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv

atman validate \
    --input-dir out \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" \
    --min-pairs 5

atman qc --input-dir out --output-dir out --rule dube

atman de \
    --input-dir out --output-dir out \
    --test paired-t --paired-by participant \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" \
    --min-pairs 5

# Provenance manifest written automatically alongside de_results.tsv:
cat out/de_results.tsv.run.json | head -40
```

Primary outputs after the `de` step: `de_results.tsv`, `de_report.tsv`,
`de_results.tsv.run.json` (the SHA-256 sidecar). A full per-command
output reference lives in the [Data Model](#data-model) section.

### Other ingest paths

Wide matrices (DIA-NN, Spectronaut, MaxQuant/LFQ, SomaScan, or any
log2-intensity export) land through `atman ingest-matrix` instead:

```bash
atman ingest-matrix \
    --matrix matrix.tsv --samples sample_metadata.tsv \
    --orientation proteins-rows \
    --platform diann_report \
    --abundance-unit log2_diann_pg_quantity \
    --assay-id-col Protein.Group --gene-col Genes \
    --condition-col diagnosis \
    --normalize median \
    --output-dir out
```

See [`docs/tutorial.md`](docs/tutorial.md) for the end-to-end Dube
walkthrough and [`docs/recipes.md`](docs/recipes.md) for the DE test
cookbook — OLS formulas, mixed models, limma + DEqMS, msqrob, TREAT,
post-hoc Sidak / Tukey / Dunnett, cross-method ensemble, bootstrap
intervals, null calibration, module scoring, ORA, meta-analysis, and
the `atman run` plan-manifest format.

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
`loadings.tsv.run.json` alongside `loadings.tsv`). The sidecar has the
following fields:

| Field | Type | Meaning |
|---|---|---|
| `command` | string | The atman subcommand that was run (e.g. `"decompose ica"`, `"de"`). |
| `args` | object | Fully-resolved argument dict, defaults included, kebab-case keys matching the CLI flag names. |
| `atman_version` | string | `CARGO_PKG_VERSION` at build time. |
| `atman_git_sha` | string | Git SHA at build time (set via `build.rs`; empty when built outside a git checkout). |
| `os_arch` | string | Build-time Rust target triple (e.g. `aarch64-apple-darwin`). |
| `input_dir_sha256` | string | SHA-256 over a sorted `basename\tsha256(contents)` roll-up of the canonical input TSVs the command actually read. Missing files are skipped. See note below on the planned per-file replacement. |
| `output_files` | object | `{<output_path>: "sha256:<hex>"}` for every artifact written, excluding the sidecar itself. |
| `started_at` | string | ISO-8601 UTC timestamp (`YYYY-MM-DDTHH:MM:SSZ`) at the start of the invocation. |
| `finished_at` | string | Same, at end of invocation. |
| `exit_code` | integer | `0` on success. Sidecar is only written on success, so the field is uniformly `0`; retained so batch auditors can filter this field without branching on its presence. |

Commands that write a sidecar: `ingest-matrix`, `validate` (when
`--report` is set), `report qc`, `decompose ica`, `decompose unmix`,
`align programs`, `align project`, `align bootstrap`, `coupling`,
`null`, `enrich gprofiler`, `de`, `bootstrap protein`,
`bootstrap module`, `bootstrap program`, `meta`, and `ratio`. The
plan-level `atman run` manifest covers its own provenance independently.

#### Reproducing a single seed

For multi-seed commands, per-seed PRNG state derives deterministically
from the master `--seed` recorded in the sidecar:

- **Multi-seed FastICA** (`decompose ica --n-seeds N`) uses
  `seed_i = seed.wrapping_add(i)` for `i ∈ [0, n_seeds)`. To reproduce
  the `i`-th seed's result in isolation, run a single-seed
  `decompose ica` with `--seed (sidecar.seed + i)` and `--n-seeds 1`.
- **Bootstrap iterations** (`align bootstrap`, `bootstrap protein`,
  `bootstrap module`, `decompose unmix --n-boot`) derive per-iteration
  sub-seeds from `(seed, iter)` via SplitMix64 (see
  `atman_core::decompose_unmix::bootstrap_sub_seed` / the identical
  function in `align_bootstrap`). Same scheme across all bootstrap
  modules, so the same top-level `--seed` produces the same iteration
  stream everywhere.

#### Resolved vs requested parameters

Some arguments carry a selection rule rather than a concrete value
(e.g. `decompose ica --k cumulative-variance=0.80` picks `k` at
runtime from the spectrum; `decompose unmix --k auto` sweeps and
picks the elbow). The sidecar today records the **rule** in `args`
but not the **resolved value**. If a reviewer needs the integer `k`
that was actually used, read it from the primary output TSV
(`loadings.tsv` column count) or the companion diagnostics TSV
(`k_selection.tsv` for `decompose unmix --k auto`). Recording the
resolved value alongside the rule is a planned sidecar schema
change; see below.

#### Planned schema evolution

Known asymmetries against `output_files` and build-environment gaps
will be addressed in a follow-on commit that bumps the sidecar to
`schema_version: 1`:

- `inputs_sha256: {relative_path: sha256}` replaces the roll-up
  `input_dir_sha256`, matching the per-file shape already used by
  `output_files`.
- `schema_version: 1` for forward-compatible parsers.
- `build_env` block — `rustc_version`, `cargo_lock_sha256`,
  `profile` — so reproductions can pin the FP-reduction-order-
  sensitive build environment behind the byte-exact Dube and
  `Δ = 1.05e-15` log2-fold-change parity claims.
- `reinvoke` — reconstructed command-line string for paste-and-run.
- `run_uuid` — UUID per invocation so reviewers can cross-reference
  independent reproductions.
- `cwd_at_start` — explicit so relative paths in the sidecar are
  anchored.
- Resolved-value fields alongside rule-based args (e.g.
  `k_resolved` when `--k` carried a rule).

Old sidecars stay readable; the schema is additive.

### Plan manifests (`atman run`)

For end-to-end reproducibility beyond the single-command sidecar,
`atman run` executes a declarative plan file (YAML or JSON) stage by
stage and writes `plan_manifest.tsv` with SHA-256 hashes of every
declared input and output per stage, atman version, OS/arch, exit
code, and per-stage wall-clock. If the plan content changes, atman
refuses to overwrite the manifest unless the plan's `plan_commit` tag
is bumped (or `--allow-drift` is passed). Full schema and example in
[`docs/recipes.md`](docs/recipes.md).

### Network-dependent commands

One command reaches external infrastructure: `atman enrich gprofiler`
calls the live g:Profiler REST endpoint. Responses are cached to disk
keyed by a SHA-256 of the canonicalized request (genes, background,
organism, sources, threshold method, user threshold, pinned ontology
version), and `--offline` fails on cache miss so a pre-populated cache
reproduces byte-identical runs without a network round-trip. For a
fully-offline pipeline, populate the cache once and then run with
`--offline`.

All other commands read only local files.

## Non-goals

Atman deliberately does not cover:

- **Single-cell proteomics.** No cell-level quantification, no
  SCP-specific normalization; the data model is sample × protein.
- **Bayesian DE / posteriors.** All variance shrinkage is
  empirical-Bayes (limma-style `fit_f_dist`), not hierarchical
  posterior sampling. No MCMC.
- **Network inference.** `atman network influence` scores hub
  centrality on a given adjacency; it does not learn the adjacency.
  Graphical-lasso / causal-discovery are out of scope.
- **MS raw file handling.** Ingest starts from peptide-level or
  protein-level matrices (MaxQuant `peptides.txt`, DIA-NN report,
  Spectronaut report, SomaScan RFU, Olink NPX). Upstream feature
  extraction (MaxQuant / FragPipe / DIA-NN) is assumed.
- **VSN, ComBat-style batch correction, plate bridging.** Explicit
  §Adapter Responsibilities above; belongs in the adapter or
  upstream pipeline.
- **Survival analysis, Cox regression, Kaplan-Meier.** Outside the
  DE-focused scope.

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
