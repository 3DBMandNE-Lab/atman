# Atman Reference

Detailed reference for validation, input formats, data model,
provenance, and scope. For getting started, see the
[README](../README.md). For workflows, see
[tutorial.md](tutorial.md) and [recipes.md](recipes.md).

## Validation Details

### msqrob2 parity note

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

### Reproducibility check

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

## Input Support

Atman has three input modes, all landing on the same canonical TSV schema:

- **Matrix ingest:** `atman ingest-matrix` reads common wide protein matrices
  (DIA-NN, Spectronaut, MaxQuant/LFQ, SomaScan, any log2 intensity/abundance
  export) plus sample/protein metadata and writes Atman's canonical TSVs.
  Supports `--log2-transform` for linear-scale inputs and both
  `proteins-rows` / `samples-rows` orientations.
- **Canonical TSV adapters:** Atman has no built-in platform-native
  ingest. Vendor long-formats land on the canonical schema via Python
  adapters in `adapters/`. The bundled `adapters/generic/olink_explore_to_atman.py`
  handles Olink Explore NGS NPX CSV exports. Other formats follow the
  same pattern: read upstream, emit `samples.tsv`, `proteins.tsv`,
  `measurements.tsv`. Once those exist, Atman's downstream commands
  (`de`, `module-de`, `robustness`, etc.) run the same way regardless of
  source.

Adapters are intentionally thin: normalize source metadata, map samples and
proteins, log-transform linear intensities when needed, and emit Atman's TSV
schema. See `adapters/examples/` for tiny synthetic DIA-NN, MaxQuant, and
SomaScan-style fixtures end-to-end.

### Adapter responsibilities

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

Canonical TSV readers are header-indexed rather than position-indexed.
Required columns must be present, but column order may vary and extra
columns may be appended without changing downstream behavior.

## Run Sidecars

Every analysis subcommand that writes output files also writes a
`<primary_output>.run.json` sidecar next to its main artifact (e.g.
`loadings.tsv.run.json` alongside `loadings.tsv`). The sidecar has the
following fields:

| Field | Type | Meaning |
|---|---|---|
| `schema_version` | integer | Sidecar schema version. Current value is `1`. |
| `run_uuid` | string | Per-invocation UUID. |
| `command` | string | The atman subcommand that was run (e.g. `"decompose ica"`, `"de"`). |
| `reinvoke` | string | Best-effort paste-and-run reconstruction of the CLI invocation. |
| `args` | object | Fully-resolved argument dict, defaults included, kebab-case keys matching the CLI flag names. |
| `atman_version` | string | `CARGO_PKG_VERSION` at build time. |
| `atman_git_sha` | string | Git SHA at build time (set via `build.rs`; empty when built outside a git checkout). |
| `build_env` | object | Build-time `rustc` version, `Cargo.lock` hash, profile, and target triple. |
| `cwd_at_start` | string | Working directory at invocation start. |
| `os_arch` | string | Build-time Rust target triple (e.g. `aarch64-apple-darwin`). |
| `inputs_sha256` | object | Per-input SHA-256 dict: `{<label>: "sha256:<hex>"}`. Missing optional inputs are omitted, not hashed as sentinel values. |
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

`reinvoke` is a convenience field, not the source of truth. The canonical
replay data is the `args` object. Standard presence flags are rendered as
bare flags such as `--offline`, while explicit bool-valued options that take
a value are rendered with `true` or `false`, for example
`de --trend false --robust true`.

### Reproducing a single seed

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

### Resolved vs requested parameters

Some arguments carry a selection rule rather than a concrete value
(e.g. `decompose ica --k cumulative-variance=0.80` picks `k` at
runtime from the spectrum; `decompose unmix --k auto` sweeps and
picks the elbow). The sidecar records both the **rule** in `args`
and the **resolved value** (e.g. `k_resolved`) so reviewers can
see the runtime decision without re-reading output files.

### Plan manifests (`atman run`)

For end-to-end reproducibility beyond the single-command sidecar,
`atman run` executes a declarative plan file (YAML or JSON) stage by
stage and writes `plan_manifest.tsv` with SHA-256 hashes of every
declared input and output per stage, atman version, OS/arch, exit
code, and per-stage wall-clock. If the plan content changes, atman
refuses to overwrite the manifest unless the plan's `plan_commit` tag
is bumped (or `--allow-drift` is passed). Full schema and example in
[recipes.md](recipes.md).

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
  adapter responsibilities above; belongs in the adapter or
  upstream pipeline.
- **Survival analysis, Cox regression, Kaplan-Meier.** Outside the
  DE-focused scope.
