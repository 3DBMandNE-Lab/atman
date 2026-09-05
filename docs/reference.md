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

### fgsea parity note

Atman's `enrich gsea` and `fgsea::fgseaSimple` use the same weighted
Kolmogorov-Smirnov enrichment statistic (Subramanian et al. 2005) and
the same gene-permutation null formulation (Korotkevich et al. 2019).
The agreement claim is therefore on the deterministic part of the
algorithm — the enrichment score (ES) — which depends only on the
ranked gene list and the set membership mask.

- **Enrichment score (ES).** Matches `fgseaSimple` within `1e-12` on
  the planted fixture in `gsea_reference.R` for all three sets
  (`top_loaded`, `bottom_loaded`, `scattered`).
- **Permutation p-value and NES.** Both quantities depend on the
  permutation null draws. Atman uses the deterministic `Xoshiro256pp`
  PRNG seeded via `--seed`; fgsea uses R's Mersenne Twister. The
  per-set p-value and NES are therefore reproducible across atman
  re-runs at the same seed but are not byte-equal to fgsea. For paper
  reporting, cite atman's seed and `n-permutations` so reviewers can
  re-derive the exact p-values from the published artifact.
- **NES sign and significance bucket.** ES sign matches fgsea exactly
  (deterministic), so the direction of enrichment is always
  consistent. Significance buckets agree on the planted fixture; on
  noisier real-world inputs the rank-order of marginal sets may
  differ at the third decimal due to permutation noise — increase
  `--n-permutations` if a tight rank-order needs to be reported.

### singscore parity note

Atman's `score signatures --method singscore` and `singscore::simpleScore`
use the same per-sample rank construction (average-rank tie handling,
matching R's `rank(..., ties.method = "average")`) and the same centered,
scale-normalized TotalScore formula:

```text
    score = (2 * mean_rank - n - 1) / (2 * (n - m))
```

where `n` is the count of measured proteins in the sample and `m` is the
count of signature genes observed in that sample. The score lies in
`[-0.5, 0.5]` and matches `singscore::simpleScore` with `knownDirection = TRUE`
and `centerScore = TRUE` (its default). The algorithm is fully
deterministic — no permutations, no randomness — so byte-stability is
expected across re-runs and across implementations.

Validated on the planted fixture in `singscore_reference.R` (4 samples × 30
genes, three signatures: up-loaded, down-loaded, scattered). All 12
(set × sample) reference rows match `singscore::simpleScore` within
`1e-12`. See `singscore_reference.tsv` for the reference values.

## Commands

### atman decompose nmf

Non-negative matrix factorization via multiplicative updates (Brunet et al.
2004). Recovers *k* non-negative latent protein programs whose product
reconstructs the input matrix. Multi-seed stability framework assesses
robustness across random initializations.

**Flags:**
- `--k <INT>` — fixed rank; required when `--k-selection fixed` (the default)
- `--k-selection <RULE>` — `fixed` (default, uses `--k`), `cophenetic-knee`, or `rss-knee`
- `--k-min <INT>` — minimum k for the automatic sweep (≥ 2); required when `--k-selection` is not `fixed` (no default)
- `--k-max <INT>` — maximum k for the automatic sweep (≤ 50, > `--k-min`); required when `--k-selection` is not `fixed` (no default)
- `--beta-loss <STR>` — loss function: `frobenius` (default, squared error) or `kullback-leibler` (alias `kl`)
- `--init <STR>` — initialization: `nndsvda` (default, deterministic) or `random`
- `--solver <STR>` — update algorithm: `mu` (default; the only supported value — multiplicative updates)
- `--max-iter <INT>` — iteration limit (default 400)
- `--tol <FLOAT>` — convergence tolerance on the per-iteration change in Frobenius error (default 1e-6)
- `--seed <INT>` — PRNG seed, used only when `--init random` (default 42)
- `--n-seeds <INT>` — number of independent runs (default 1); `> 1` enables multi-seed stability filtering
- `--seed-base <INT>` — base seed for multi-seed runs, `seed_i = seed_base + i` (default 42)
- `--min-stable-seed-fraction <FLOAT>` — fraction of seeds a program must survive in to pass multi-seed filtering (default 0.9)
- `--stability-metric <STR>` — reproducibility metric: `jaccard-top20` only (Jaccard overlap of the top-`--stability-top-n` loadings)
- `--stability-top-n <INT>` — top *n* loadings used by the stability metric (default 20)
- `--max-missing-fraction <FLOAT>` — drop assays whose missing-sample fraction exceeds this value (default 0.0, strict complete-case; mirrors `decompose ica`); the sidecar records the resolved value plus `n_assays_retained`/`n_assays_dropped_missingness`
- `--transform <STR>` — pre-decomposition transform: `none` (default; NMF requires non-negative input and rejects negative values loudly), `exp2-clip` (`2^clamp(x, -c, +c)` — restores a non-negative ratio scale from log2-ratio input while winsorizing extreme tails), or `shift-min` (`x - min(X)` over the whole matrix — a sensitivity alternative to `exp2-clip`)
- `--transform-clamp <FLOAT>` — clamp radius `c` for `--transform exp2-clip` (default 6.0 when omitted); only valid together with `--transform exp2-clip` — a hard error otherwise
- `--output-loadings <PATH>` — protein weights per component (required)
- `--output-activations <PATH>` — per-sample component activations (optional)
- `--output-stability <PATH>` — cross-seed reproducibility scores (written only when `--n-seeds > 1`)
- `--output-k-sweep <PATH>` — k-selection diagnostic metrics (written only when `--k-selection` is automatic)

**Output formats:**
- `loadings.tsv`: `program\tassay_id\tgene_symbol\tloading` (long format, one row per assay per program)
- `activations.tsv`: `sample_id\tprogram\tactivation`; long format compatible with `--adjust-for`
- `stability.tsv`: `program\tstable_seed_fraction\tn_seeds_present`
- `k_sweep.tsv`: `k\tcophenetic\tmean_rss\tmean_kl\tselected`

See `docs/recipes.md` for the canonical Fig 3 admixture-adjusted DE chain.

### atman decompose ica --missingness-model abundance-conditional

FastICA with MNAR (missing-not-at-random) awareness via abundance-conditional
detection-curve model. On fully-observed data, output is numerically identical
to standard FastICA (MAR-collapse property). When data contain structured
missingness (below-LOD entries), the method iterates jointly over component
estimation and detection-curve fit.

**New flag:**
- `--max-joint-iter <INT>` — iteration limit for joint MNAR fit (default 50)

Existing `--k-selection`, `--n-seeds`, `--seed`, `--output-loadings`,
`--output-activations`, `--output-stability` flags work as in standard ICA.

### atman de --adjust-for <PATH>

Admixture-adjusted differential abundance testing. Accepts per-sample
covariates from an external TSV and composes them into the design matrix
alongside condition and `samples.tsv` columns.

**Format detection:** `--adjust-for` auto-detects two formats:
- **Wide:** sample_id column plus covariate columns; one row per sample
- **Long:** sample_id, program, activation columns; multiple rows per sample (auto-pivoted to wide)

NMF activation outputs use long format and are consumed directly. Both formats
can be combined: if both `--adjust-for <EXTERNAL_TSV>` and `--covariates
age,sex,batch` are passed, the design matrix is `~ condition +
[external_covariates] + age + sex + batch`.

**Implementation:** Wired through `--test` paths:
- `limma` — design matrix composed and passed to `lmFit` + eBayes
- `msqrob` — external covariates added to protein-level model
- `ols` — design formula extended with external covariate columns
- `mixed` — fixed-effects formula extended
- `welch-t` — routes to plain OLS internally, with the external covariates in
  the design (no HC3 robust SE — that would require a new dependency)
- `paired-t` — routes to a paired-difference ANCOVA: per-subject
  `d = abundance_b − abundance_a` regressed via OLS on the per-subject
  covariate differences (`d ~ 1 + Δcov_1 + ...`); the intercept is the
  adjusted mean paired difference

(see tests K6/K7 in `crates/atman/tests/de_adjust_for.rs`, which parity-check
both routings against R references)

**Sidecar:** `de_results.tsv.run.json` records `inputs_sha256` entry for the
`--adjust-for` path and its input hash, enabling full reproducibility.

### atman bench decompose --tools atman.<method>

Native dispatch of decomposition methods under `atman bench decompose`. Allowed
methods: `ica`, `nmf`, `missingness-ica`. The bare `atman` method (no prefix)
remains an alias for `atman.ica`.

**New flag:**
- `--fixture-nmf <PATH>` — NMF-appropriate benchmark fixture (defaults to
  `bench/planted_archetypes_nmf_v1/`)

**Fixture requirements:** NMF fixtures must include a ground-truth component
matrix (`loadings_truth.tsv`) and expected reconstruction error targets for
correctness validation.

### atman align programs / bootstrap / project — decomposition method tracking

Sidecars for all `align` subcommands now record:

```json
{
  "args": {
    "decomposition_method": "ica | nmf | mixed | unknown"
  }
}
```

This enables downstream auditing: readers can trace whether a set of aligned
programs came from ICA, NMF, or a mixed ensemble, and reproduce the exact
alignment parameters from the sidecar.

**`align bootstrap --decomposition`:** per-resample decomposition method,
applied identically to the point estimate, every bootstrap resample, and
every jackknife replicate.
- `--decomposition <STR>` — `ica` (default, FastICA, byte-identical to prior
  releases) or `nmf` (single-seed multiplicative-updates NMF per resample)
- `--beta-loss <STR>` — NMF loss: `frobenius` (default) or `kullback-leibler`
  (alias `kl`). Ignored for `--decomposition ica`.
- `--init <STR>` — NMF initialization: `random` (default, seeded from the
  same SplitMix64 derivation the ICA branch uses) or `nndsvda`. Ignored for
  `--decomposition ica`.
- `--nmf-max-iter <INT>` — NMF max multiplicative-update iterations per seed
  per cohort (default 500). Only used with `--decomposition nmf`.
- `--nmf-tol <FLOAT>` — NMF convergence tolerance (default 1e-5). Only used
  with `--decomposition nmf`.
- `--transform <STR>` — pre-decomposition transform for `--decomposition
  nmf`: `none` (default; requires non-negative input, rejected loudly
  otherwise), `exp2-clip`, or `shift-min`. Re-applied fresh to every
  per-resample matrix (point estimate, each bootstrap resample, each
  jackknife replicate) rather than cached once. Ignored for `ica`.
- `--transform-clamp <FLOAT>` — clamp radius `c` for `--transform exp2-clip`
  (default 6.0 when omitted); only valid together with `--transform
  exp2-clip` — a hard error otherwise.

Non-finite (`NaN`/±∞) loadings from either decomposition are rejected loudly,
naming the cohort and call site (point estimate / bootstrap iteration /
jackknife replicate), rather than flowing silently into cosine similarity,
archetype grouping, and the bootstrap/BCa accumulators downstream. The
sidecar's `decomposition_method` echoes the `--decomposition` flag value, and
`transform_clamp` records the *resolved* clamp (e.g. 6.0 when `--transform
exp2-clip` is passed without an explicit `--transform-clamp`), not the raw,
possibly-null CLI value.

**`align project --transform`:** gains the NMF input transforms alongside the
existing compositional ones. Valid values: `none`/`clr`/`alr`/`ratio-anchor`
(compositional transforms; `ilr` is parsed but rejected in `align project` —
it reorders coordinates so atlas labels would no longer line up with the
transformed cohort columns) and `exp2-clip`/`shift-min` (the same NMF input
transforms as `decompose nmf`, shared via `atman_core::nmf::apply_transform`).
`--transform-clamp` follows the same rule as `decompose nmf` (only valid with
`exp2-clip`, default 6.0 when omitted). `shift-min` always recomputes its
shift on the cohort being projected — it does not reuse the atlas's
training-time shift. The sidecar's `TransformRecord` fields (`transform_clamp`,
`transform_shift`) record what was actually applied.

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
| `gsea_reference.R` | fgsea (Bioconductor ≥ 1.30) |
| `singscore_reference.R` | singscore (Bioconductor ≥ 1.24) |
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
`bootstrap module`, `bootstrap program`, `meta`, `ratio`,
`network influence`, `modules discover`, and `programs filter`. The
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

### atman axes

Subject-level score tables. Every `axes` subcommand reads a wide score TSV
(`--scores`: `sample_id`, `cohort`, `condition`, `is_control`, numeric score
columns), optionally joins `--cohort-dirs cohort=dir,...` (`samples.tsv`
covariates) and `--covariates-tsv PATH` files on `sample_id`, and writes a
TSV plus `<output>.run.json`.

Covariate expressions (used by `--designs`, `--anchors`, `--median-cols`):
identifiers, numbers, `+ - * /`, parentheses, `log10()`, `log2()`, `ln()`,
and `z()` (standardize over the rows entering the fit; outermost only).
Columns whose trimmed non-empty values all parse as numbers are numeric;
other columns are categorical, one-hot with the alphabetically first level as
reference (`sex` F/M ⇒ one `sexM` column).

#### atman axes build

`--activations` (one wide table, or `label=path,...` of `align project`
outputs joined on `sample_id` with `label_` prefixes), `--cohort-dirs`,
`--representatives axis=column,...`, `--orthogonalize axisA,axisB --against
axisR --within cohort`, `--output`. Orthogonalization is the OLS residual
(with intercept) of the representative on the reference, fitted separately
per `--within` group. Output: input columns, then per axis `<axis>_raw`,
`<axis>_z` (global z, ddof = 1) and, for orthogonalized axes,
`<axis>_unorth_raw`, `<axis>_unorth_z`.

#### atman axes contrast

`--scores`, `--cohort-dirs`, `--covariates-tsv` (repeatable), `--manifest`
(`label, cohort, case, control, family[, condition_col, subset]`;
`case`/`control` accept `a|b` lists, `control=*` = all other levels within the
subset), `--score-cols`, `--designs "~ case; ~ case + z(age) + sex"`,
`--n-bootstrap`, `--seed`, `--ci` (default 0.95), `--output`,
`--output-covariates`, `--output-bootstrap`. One row per contrast × score ×
design. Unadjusted columns (`n_case … auc`) use every subject with a score;
`n_fit` onward use the complete-case rows of the design. `q` is BH over rows
sharing (family, design); `welch_q` likewise over `welch_p`. Bootstrap:
subjects resampled within case and control separately, design rebuilt per
resample (so `z()` is re-standardized), percentile interval at `--ci`; one
SplitMix64 stream per contrast seeded by `derive_sub_seed(seed, contrast_index)`.
Multi-level categorical covariates are one-hot encoded; `--omnibus-factor
COL` (repeatable) adds `omnibus_f` rows to `--output-covariates` with the
joint F-test of that factor's columns (`f, df_num, df_den, p`).

#### atman axes groups

`--cohort`, `--group-by COL`, `--groups a,b,...` (order; default sorted
observed values), `--score-cols`, `--median-cols "age,QIgG/QAlb"`,
`--reference GROUP`, `--ci`, `--output` (`group, n, <expr>_median…,
<score>_mean, <score>_ci_lo, <score>_ci_hi`; t-based CI, n ≥ 3),
`--output-tests` (`kind, score, reference, group, n_ref, n_group, statistic,
p`: `kruskal` rows carry H and its chi-square p; `reference_vs_group` rows
carry pooled-SD Cohen d of reference minus group and the Welch p).
