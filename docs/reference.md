# Atman Reference

Detailed reference for validation, input formats, data model,
provenance, and scope. For getting started, see the
[README](../README.md). For workflows, see
[tutorial.md](tutorial.md) and [recipes.md](recipes.md).

## Validation Details

### msqrob2 parity note

Atman's `--test msqrob` and msqrob2 in the canonical QFeatures vignette
workflow are not the same estimator. The 0.118 log₂ median difference
reflects the method choice rather than a numerical gap:

- **msqrob2 here (aggregate-then-fit).** Peptides are summarised to one
  value per (protein, sample) via `matrixStats::colMedians`, then
  `msqrob(~condition)` fits per-protein OLS with empirical-Bayes
  variance shrinkage.
- **Atman (joint fit).** One peptide-level linear mixed model per
  protein, random intercept per peptide, ridge penalty on non-intercept
  fixed effects, REML 1D profile over `τ = σ²_peptide / σ²_res`, with
  empirical-Bayes variance shrinkage across all fitted proteins.

The substantive claims are sign agreement and spike-in recovery. The
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
  PRNG seeded via `--seed`. The fgsea package uses R's Mersenne
  Twister. The per-set p-value and NES are therefore reproducible
  across atman re-runs at the same seed, but are not byte-equal to
  fgsea. For paper
  reporting, cite atman's seed and `n-permutations` so reviewers can
  re-derive the exact p-values from the published artifact.
- **NES sign and significance bucket.** ES sign matches fgsea exactly
  (deterministic), so the direction of enrichment is always
  consistent. Significance buckets agree on the planted fixture. On
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
- `--k <INT>` — fixed rank. Required when `--k-selection fixed` (the default)
- `--k-selection <RULE>` — `fixed` (default, uses `--k`), `cophenetic-knee`, or `rss-knee`
- `--k-min <INT>` — minimum k for the automatic sweep (≥ 2). Required when `--k-selection` is not `fixed` (no default)
- `--k-max <INT>` — maximum k for the automatic sweep (≤ 50, > `--k-min`). Required when `--k-selection` is not `fixed` (no default)
- `--beta-loss <STR>` — loss function: `frobenius` (default, squared error) or `kullback-leibler` (alias `kl`)
- `--init <STR>` — initialization: `nndsvda` (default, deterministic) or `random`
- `--solver <STR>` — update algorithm: `mu`. This is the default and the only supported value (multiplicative updates)
- `--max-iter <INT>` — iteration limit (default 400)
- `--tol <FLOAT>` — convergence tolerance on the per-iteration change in Frobenius error (default 1e-6)
- `--seed <INT>` — PRNG seed, used only when `--init random` (default 42)
- `--n-seeds <INT>` — number of independent runs (default 1). `> 1` enables multi-seed stability filtering
- `--seed-base <INT>` — base seed for multi-seed runs, `seed_i = seed_base + i` (default 42)
- `--min-stable-seed-fraction <FLOAT>` — fraction of seeds a program must survive in to pass multi-seed filtering (default 0.9)
- `--stability-metric <STR>` — reproducibility metric: `jaccard-top20` only (Jaccard overlap of the top-`--stability-top-n` loadings)
- `--stability-top-n <INT>` — top *n* loadings used by the stability metric (default 20)
- `--max-missing-fraction <FLOAT>` — drop assays whose missing-sample fraction exceeds this value (default 0.0, strict complete-case, mirroring `decompose ica`). The sidecar records the resolved value plus `n_assays_retained`/`n_assays_dropped_missingness`
- `--transform <STR>` — pre-decomposition transform: `none` (the default: NMF requires non-negative input and rejects negative values loudly), `exp2-clip` (`2^clamp(x, -c, +c)` — restores a non-negative ratio scale from log2-ratio input while winsorizing extreme tails), or `shift-min` (`x - min(X)` over the whole matrix — a sensitivity alternative to `exp2-clip`)
- `--transform-clamp <FLOAT>` — clamp radius `c` for `--transform exp2-clip` (default 6.0 when omitted). Only valid together with `--transform exp2-clip` — a hard error otherwise
- `--output-loadings <PATH>` — protein weights per component (required)
- `--output-activations <PATH>` — per-sample component activations (optional)
- `--output-stability <PATH>` — cross-seed reproducibility scores (written only when `--n-seeds > 1`)
- `--output-k-sweep <PATH>` — k-selection diagnostic metrics (written only when `--k-selection` is automatic)

**Output formats:**
- `loadings.tsv`: `program\tassay_id\tgene_symbol\tloading` (long format, one row per assay per program)
- `activations.tsv`: `sample_id\tprogram\tactivation`. Long format, compatible with `--adjust-for`
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
- **Wide:** sample_id column plus covariate columns. One row per sample
- **Long:** sample_id, program, activation columns. Multiple rows per
  sample, auto-pivoted to wide

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
  covariate differences (`d ~ 1 + Δcov_1 + ...`). The intercept is the
  adjusted mean paired difference

(see tests K6/K7 in `crates/atman/tests/de_adjust_for.rs`, which parity-check
both routings against R references)

**Sidecar:** `de_results.tsv.run.json` records `inputs_sha256` entry for the
`--adjust-for` path and its input hash, enabling full reproducibility.

### atman de — sample selection and effect-size flags

- `--include-controls` — samples flagged `is_control=1` join a group when
  their condition is named in `--groups`. Off by default. The sidecar
  records the flag.
- `--condition-col COL` — use another samples.tsv column as the condition
  label (e.g. `diagnosis_group`).
- `--subset 'col!=value'` (repeatable, `;`-separated) — keep only samples
  whose raw samples.tsv values satisfy every predicate. A missing value fails
  `==` and passes `!=`.
- `--collapse-others LABEL` — every condition not named in `--groups` becomes
  `LABEL` (which must be one side of a comparison), so
  `--condition-col diagnosis_group --groups Headache-otherRef --collapse-others otherRef`
  contrasts one level against all others.
- `--max-missing-fraction F` — drop a protein from a comparison when more
  than `F` of that comparison's samples lack a value (default 1.0 = keep
  all). Counts land in the sidecar's `missingness_filter`.
- `--test ols` now reports `effect_size` = pooled-SD Cohen d over the fitted
  subjects (`effect_size_method = cohen_d`) with its large-sample CI in
  `ci_low`/`ci_high`, and `de_results.tsv` carries `n_a`/`n_b` (observed
  subjects per group, also filled for `welch-t`).
- `--require-cols col[,col]` (repeatable) — keep only samples whose
  samples.tsv row has a non-empty value in every listed column.
  `--require-numeric col[,col]` additionally requires a finite number (so
  placeholders such as `not measured` drop out), which makes a
  relative-scale run fit the same subjects as a `scale absolute` run
  (`n_require_cols_dropped` / `n_require_numeric_dropped` in the sidecar).
- `--collapse-genes none|mean|max-observed` — how protein groups (assays)
  that share a gene symbol become one value per sample. `none` (default)
  keeps the lexically first assay id. `mean` averages the assays observed in
  that sample. `max-observed` keeps the assay observed in the most samples.
  The missingness filter runs on the collapsed gene. The sidecar's
  `gene_symbol_collapse` block records `rule`, `n_assays`,
  `n_genes_with_multiple_assays`, `n_extra_assays`, and
  `n_genes_after_collapse` (before the per-comparison missingness filter,
  whose retained counts are in `missingness_filter`). `residuals` adds
  `n_rows_after_collapse`, and `score weighted` adds
  `n_genes_shared_after_filter`. Applies to paired-t, welch-t, ols, and
  mixed. (`align programs` is last-wins for a
  duplicated label within a program: the final loadings row for that label
  overwrites earlier ones.)
- `--design` accepts covariate expressions: `~ condition + z(age) + sex +
  log10(QAlb) + log10(leukocyte_count + 1)`; `z()` standardizes over the
  fitted samples. The condition coefficient and p are invariant to it, and the
  covariate row in `de_covariates.tsv` is then per SD.
- Continuous contrasts: `--design '~ log10(QAlb) + z(age) + sex' --contrast
  'log10(QAlb)'` without `--groups` fits every sample (complete case per
  protein) and reports the named coefficient: `comparison` = the term,
  `mean_diff` = beta, `t`, `p_value`, `bh_q` over the proteins. `mean_a`,
  `mean_b`, and `effect_size` are empty. The `--contrast` text must match the
  `--design` term exactly.

### atman residuals — variance and canonical outputs

`--covariates-tsv PATH` (repeatable) joins extra `sample_id`-keyed columns
(e.g. axis scores) so `--design '~ axis1_raw'` or `'~ axis1_raw + axis2_raw'`
can regress them out. `--design` accepts covariate expressions.
`--max-missing-fraction F` drops proteins observed in fewer than `1−F` of the
design samples. `--output-variance` writes per-protein `n, ss_model,
ss_total, r2`. `--output-variance-summary` writes `design, n_samples,
n_proteins, frac_variance (Σss_model/Σss_total), median_r2, q75_r2,
frac_r2_gt_0_25`. `--output-canonical-dir DIR` writes a canonical directory
whose `measurements.tsv` carries the residual as abundance (observed cells
only, same unit label as the input) with `samples.tsv`/`proteins.tsv`
copied, so `atman de` or `atman score weighted` run directly on the
residual matrix. `--collapse-genes none|mean|max-observed` reduces assays
that share a gene symbol as in `de`. `none` keeps every assay as its own
row. The other rules emit one row per gene whose `assay_id` is the
representative assay, also in the canonical output.

### atman concordance

`--manifest` (`label, path, effect_col[, feature_col, q_col, stage]`, where
`feature_col` defaults to `gene_symbol` and is resolved per table, so a
loading vector keyed by `protein` joins a DE table keyed by `gene_symbol`),
`--pairs a:b,...` (default all label pairs within each stage), `--top-n`,
`--q-threshold`, `--n-bootstrap`, `--seed`, `--ci`, `--output`,
`--output-delta`. Rows: `a, b, stage, n_features, rho, p, ci_lo, ci_hi,
n_hits_both, sign_concordance, jaccard_top_n`. A pair present in several
stages is evaluated on the features shared by every stage, and
`--output-delta` reports `delta_rho` (later stage minus earlier) with a CI
from the same feature resamples. Features are not independent units, so the
intervals are descriptive. `--output-tree-linkage` (with optional
`--output-tree-support`, `--output-tree-newick`, and `--tree-stage` when
the manifest has several stages) builds the average-linkage tree of the
stage's tables on `1 − Spearman` over pairwise-shared features, with
support from feature-bootstrap replicates (resampling the union of features
and re-evaluating each pair on the drawn features both tables carry). The
file conventions match `axes tree`.

### atman score weighted

`--input-dir` (canonical cohort, or a `residuals --output-canonical-dir`),
`--weights PATH --feature-col gene_symbol --weight-col cohen_d`,
`--signature LABEL`, `--max-missing-fraction`, `--min-shared` (default 10),
`--output` (`sample_id, subject_id, condition, is_control, signature,
n_shared, n_used, score`), and optionally `--groups A-B --output-summary
PATH` (`signature, comparison, n_shared_proteins, n_case, n_control,
cohen_d, d_ci_lo, d_ci_hi, welch_p, auc`; groups are matched on the
condition label regardless of `is_control`). Proteins are z-scored within
the scored cohort. `score = Σ w·z / Σ|w|` over the proteins the subject has.
`--collapse-genes none|mean|max-observed` reduces assays sharing a gene
symbol as in `de` (default `none` = lexically first assay).
`--summary-only` (with `--groups` and `--output-summary`) skips the
per-sample table. The sidecar then attaches to the summary file.

### atman scale absolute

`--input-dir <canonical> --total-col total_protein --output-canonical-dir <dir>`.
For every sample with a positive total-protein value, each log2 abundance
becomes `a − log2(Σ_assays 2^a) + log2(total_protein)`, the sum running over
that sample's measured assays, so the sample's linear sum equals its total
protein. Samples without a total-protein value are dropped from the output
measurements, and the sidecar counts them. `samples.tsv` and `proteins.tsv`
are copied verbatim, and the unit label is unchanged. Used for Reiber-style
per-protein exponents: `atman de --design "~ log2(QAlb) + z(age) + sex"
--contrast "log2(QAlb)"` on the rescaled directory gives the absolute-scale
slope. The same call on the input directory gives the relative-scale slope.

### atman enrich ora --query-tsv

`enrich ora` tests either the BH hits of a `--de-results` table or, with
`--query-tsv PATH` (a `gene_symbol` column plus an optional `query` column,
each distinct value its own family), a hand-made gene list. `--query-tsv`
requires `--universe` (the measured background: a `gene_symbol` column or
one gene per line). The output gains a trailing `query` column (`all`, the
`--comparison`, or the query label). The sidecar hashes the query, set, and
universe files, so an exported HPA/GO set is pinned by content.

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
- `--threads <INT>` — worker threads for the bootstrap iterations and the
  jackknife replicates (default 0 = one per available core, `1` is serial).
  Output bytes do not depend on this value: every iteration draws from its
  own `SplitMix64(seed, iter)` sub-seed and the accumulators are folded in
  iteration order. Recorded in the sidecar.
- `--transform <STR>` — pre-decomposition transform for `--decomposition
  nmf`: `none` (default; requires non-negative input, rejected loudly
  otherwise), `exp2-clip`, or `shift-min`. Re-applied fresh to every
  per-resample matrix (point estimate, each bootstrap resample, each
  jackknife replicate) rather than cached once. Ignored for `ica`.
- `--transform-clamp <FLOAT>` — clamp radius `c` for `--transform exp2-clip`
  (default 6.0 when omitted). Only valid together with `--transform
  exp2-clip`. Anything else is a hard error.

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
(compositional transforms, except that `ilr` is parsed and then rejected
in `align project`, because it reorders coordinates and atlas labels
would no longer line up with the transformed cohort columns) and
`exp2-clip`/`shift-min` (the same NMF input transforms as
`decompose nmf`, shared via `atman_core::nmf::apply_transform`).
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
`renv.lock` or a Dockerfile is a later pass. In the meantime the
committed reference TSVs are what the parity tests diff against, so
drift on the R side cannot affect CI.

The parity assertions run on every `cargo test --workspace --release`.
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

### atman harmonize

`harmonize` compares cross-cohort harmonisation methods by held-out
transfer. It has two subcommands and a strict contract between them.

`harmonize fit` reads the training cohorts. It writes one model file.

`harmonize apply` reads that model file and one cohort. It reads nothing
else. A held-out evaluation therefore cannot leak training data. The
command refuses to apply a model to a cohort that the model was fit on.

The model file records four things: the shared feature axis, the disease
direction in that method's representation, the labels of the training
cohorts, and a SHA-256 hash of their inputs. The hash matters because a
label can be reused over different data. A replay can check the hash and
prove that the held-out cohort was absent.

**Methods** (`--method`):

- `zscore` — per-cohort standardisation. This is the null
  harmonisation. It learns nothing from the training cohorts, because a
  new cohort is standardised by its own statistics.
- `rank` — within-sample rank, rescaled to `[0, 1]`. Also stateless.
- `quantile` — quantile normalisation against a reference profile
  learned on the training cohorts. Each training subject contributes its
  own observed values, interpolated onto the common grid. A subject with
  two or more observed features counts.
- `reference-protein` — per-sample division by the geometric mean of a
  reference set. `--reference-k` sets the size. The command selects the
  set on the training cohorts only. `apply` receives the chosen set in
  the model and cannot reselect.

Read `method_is_fitted` in the model file to tell the two stateless
methods from the two fitted ones. A stateless method that matches a
fitted one is a result. It means the fitted state did not carry the
signal.

**The negative control** (`--permute-labels`):

Every method produces a score on a held-out cohort. A benchmark whose
entries all beat zero measures nothing. `--permute-labels` shuffles the
case and control labels inside each training cohort before it learns the
direction. The permuted model travels the same `apply` path.

One permuted fit is one draw from the null. It is not the null. Run
several `--seed` values and compare the observed effect against that
distribution. A single draw can be large in either direction. On a
60-feature synthetic fixture one draw reached −4.45 against a real
effect of +11.90.

**Direction centring** (`--center-direction` on `apply`):

Every per-sample normalisation leaves a per-subject term in the
harmonised value. A direction with a non-zero mean projects onto that
term. For `reference-protein` the score carries
`−anchor × mean(direction)` exactly. If the anchor differs between the
arms of the held-out cohort, any direction separates them, including a
direction learned from shuffled labels.

`--center-direction` subtracts the direction's mean over each subject's
own used features. It is per-subject and not global, because missing
data makes the usable feature set differ between subjects.

The flag is off by default, so that runs made before it existed stay
reproducible. `apply` warns whenever the flag is off and the method is a
per-sample normalisation. Turn the flag on unless you have a reason not
to.

Centring has a cost. A disease effect that raises every protein by the
same amount is not separable from an anchor shift under a per-sample
normalisation. Centring does not discard recoverable signal. It declines
to claim unrecoverable signal.

### atman decompose ica — missingness and convergence flags

- `--weighted-whitening` — use the detection model to weight the
  whitening. An observed cell counts fully. An undetected cell counts by
  `1 − P(detected)` at its imputed value. The contrast function stays
  unweighted. Off by default.
- `--degeneracy-floor <FLOAT>` — stop the joint loop when the imputed
  cells come within this fraction of their first-iteration distance to
  their own reconstruction. Default 0.1. Set 0 to disable the guard.
- `--joint-tol <FLOAT>` — convergence tolerance for the joint loop.
  Default 1e-6.

The joint loop is iterated reconstruction-imputation. Its fixed point is
the state where the imputed cells equal their own reconstruction. At
that point the missing cells carry no independent information, and the
fit explains them by construction. Running the loop to convergence
therefore produces a worse model than stopping earlier.

Measured on the committed ground-truth fixture over six seeds, mean
absolute error against the planted sources:

| Configuration | MAE |
|---|---|
| mean-imputation, then ICA | 0.0768 |
| joint loop, detection curve discarded | 0.0759 |
| weighted whitening inside the joint loop | 0.0754 |
| `--weighted-whitening --max-joint-iter 1` | **0.0737** |

The last row is best on every seed. Prefer it.

`mnar_joint_stop_reason` in the sidecar records why the loop stopped. It
takes three values: `converged`, `degeneracy-floor`, and
`max-iterations`. A stop on the degeneracy floor is not a convergence,
and `ica_converged` stays false for it.

`mnar_joint_trace.tsv` records `beta0`, `beta1`, the step size, the
change in the imputed cells, and the reconstruction gap at every joint
iteration. Read `imputation_delta` to see whether the state is settling.
Read `reconstruction_gap` to see whether it is settling on the
degenerate fixed point.

### atman align — metric selectivity

`--metric cosine-centered` mean-centres each loading vector before the
cosine. Use it whenever you compare non-negative loadings against signed
ones.

Plain cosine has a positivity floor on non-negative vectors. Two NMF
programs cannot score below zero, and in practice they sit high because
they share a baseline. One threshold is then selective for signed
loadings and inert for non-negative ones. Measured across six CPTAC
cohorts on a shared 4,375-gene universe: the median cross-cohort cosine
was −0.004 for ICA and +0.813 for NMF, so `--tau 0.30` admitted 7.7% of
ICA pairs and 99.9% of NMF pairs.

`align programs` and `align bootstrap` report what fraction of
cross-cohort pairs the supplied tau admits. They warn above 90%. Above
that fraction the threshold is not thresholding, and reciprocal-best
matching alone decides the archetypes. A sweep over tau then looks
robust because it is saturated.

`align bootstrap` applies the metric twice. `--cosine-tau` groups
bootstrap programs into archetypes. `--match-tau` matches those
archetypes back to the point estimate. Both gates report their admitted
fraction.

### atman network differential — bounded edge-pairwise output

`--mode edge-pairwise` streams the leading `--top-rows` rows through a
bounded heap. A capped run holds `--top-rows` candidates whatever the
input size.

An uncapped run (`--top-rows 0`) counts the candidate rows first. It
refuses above 100 million rows and names the count. Pass `--top-rows N`
to keep the N largest-`|z_diff|` rows in bounded memory, or use
`--mode edge-summary`.

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
   `decompose ica` exposes `--impute mean|none`. The other analysis
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
| `atman_git_sha` | string | Git SHA at build time, set via `build.rs`. Empty when built outside a git checkout. |
| `build_env` | object | Build-time `rustc` version, `Cargo.lock` hash, profile, and target triple. |
| `cwd_at_start` | string | Working directory at invocation start. |
| `os_arch` | string | Build-time Rust target triple (e.g. `aarch64-apple-darwin`). |
| `inputs_sha256` | object | Per-input SHA-256 dict: `{<label>: "sha256:<hex>"}`. Missing optional inputs are omitted, not hashed as sentinel values. |
| `output_files` | object | `{<output_path>: "sha256:<hex>"}` for every artifact written, excluding the sidecar itself. |
| `started_at` | string | ISO-8601 UTC timestamp (`YYYY-MM-DDTHH:MM:SSZ`) at the start of the invocation. |
| `finished_at` | string | Same, at end of invocation. |
| `exit_code` | integer | `0` on success. Sidecar is only written on success, so the field is uniformly `0`. It is retained so batch auditors can filter this field without branching on its presence. |

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
runtime from the spectrum, and `decompose unmix --k auto` sweeps and
picks the elbow). The sidecar records both the **rule** in `args`
and the **resolved value** (e.g. `k_resolved`) so reviewers can
see the runtime decision without re-reading output files.

The same convention covers rules that fail rather than resolve. A
recorded value cannot show that a rule was never satisfied, so each of
these fields records what happened next to what was asked for:

| Field | Command | Meaning |
|---|---|---|
| `k_selection_bound_hit` | `decompose ica`, `decompose nmf`, `decompose unmix` | `k_max`, `k_min`, or null. A rule that returns its own search bound has not selected anything. It ran out of room, and its criterion may never have been met. An explicit `--k` records null, because that is not a selection. |
| `scale-free-fit-achieved` | `modules discover` | Whether any soft power met the scale-free criterion. `scale-free-best-r-squared` and `soft-power-fallback-rule` record how far short the data fell and which rule supplied the power. |
| `ica_converged`, `ica_n_iterations`, `ica_final_tol` | `decompose ica` | Whether the reference fit converged, and the evidence. `ica_n_seeds_not_converged` covers the alternative seeds behind the stability ranking. |
| `abundance_n_not_converged` | `decompose unmix` | Per-subject abundance solves that hit `--fcls-max-iter`. Always 0 for `--abundance ucls`, which is closed-form. |
| `ridge-lambda-resolved` | `de` | The penalty that applied. Null when the shrinkage path did not run, so a sidecar naming a ridge value is not evidence that msqrob ran. |
| `stability_alt_seed_method` | `decompose ica` | Which method the alternative seeds ran. The seed-stability column is only comparable when it matches the reference method. |
| `n_empty_outputs` | `atman run` | Declared outputs that exist and carry no data row. Existence and non-emptiness are different checks.

### Plan manifests (`atman run`)

For end-to-end reproducibility beyond the single-command sidecar,
`atman run` executes a declarative plan file (YAML or JSON) stage by
stage and writes `plan_manifest.tsv` with SHA-256 hashes of every
declared input and output per stage, atman version, OS/arch, exit
code, and per-stage wall-clock. If the plan content changes, atman
refuses to overwrite the manifest unless the plan's `plan_commit` tag
is bumped (or `--allow-drift` is passed). Full schema and example in
[recipes.md](recipes.md).

Plans may declare `vars:` (name → string), substituted as `${name}` in
every stage's `command`, `inputs`, and `outputs`. An undefined name is an
error, and a bare `$` is left to the shell. `--dry-run` validates the plan
and lists the resolved stages without executing: every input must exist
under `--input-dir` or be declared as an output of an earlier stage.
After each stage, any `<output>.run.json` sidecar found next to a declared
output is hashed into the manifest's trailing `sidecar_hash` column, and
`--strict-outputs` (default `true`) records `exit_code 2` and aborts when a
declared output is missing afterwards (`--strict-outputs false` restores the
permissive `MISSING` behaviour).

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
  SCP-specific normalization. The data model is sample × protein.
- **Bayesian DE / posteriors.** All variance shrinkage is
  empirical-Bayes (limma-style `fit_f_dist`), not hierarchical
  posterior sampling. No MCMC.
- **Network inference.** `atman network influence` scores hub
  centrality on a given adjacency. It does not learn the adjacency.
  Graphical-lasso / causal-discovery are out of scope.
- **MS raw file handling.** Ingest starts from peptide-level or
  protein-level matrices (MaxQuant `peptides.txt`, DIA-NN report,
  Spectronaut report, SomaScan RFU, Olink NPX). Upstream feature
  extraction (MaxQuant / FragPipe / DIA-NN) is assumed.
- **VSN, ComBat-style batch correction, plate bridging.** Explicit
  adapter responsibilities above. This belongs in the adapter or the
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
and `z()`, which standardizes over the rows entering the fit and is
accepted at the outermost level only. Columns whose trimmed non-empty
values all parse as numbers are numeric. Other columns are categorical,
one-hot with the alphabetically first level as reference (`sex` F/M ⇒
one `sexM` column).

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
(`label, cohort, case, control, family[, condition_col, subset]`, where
`case`/`control` accept `a|b` lists, `control=*` = all other levels within the
subset), `--score-cols`, `--designs "~ case; ~ case + z(age) + sex"`,
`--n-bootstrap`, `--seed`, `--ci` (default 0.95), `--output`,
`--output-covariates`, `--output-bootstrap`. One row per contrast × score ×
design. Unadjusted columns (`n_case … auc`) use every subject with a score.
`n_fit` onward use the complete-case rows of the design. `q` is BH over rows
sharing (family, design), and `welch_q` likewise over `welch_p`. Bootstrap:
subjects resampled within case and control separately, design rebuilt per
resample (so `z()` is re-standardized), percentile interval at `--ci`. There
is one SplitMix64 stream per contrast, seeded by
`derive_sub_seed(seed, contrast_index)`.
Multi-level categorical covariates are one-hot encoded. `--omnibus-factor
COL` (repeatable) adds `omnibus_f` rows to `--output-covariates` with the
joint F-test of that factor's columns (`f, df_num, df_den, p`). A bootstrap
replicate whose resample leaves a categorical term with one level (or whose
fit is singular) is skipped and counted in `boot_n_skipped`. `boot_n` is
the number of replicates that contributed.

#### atman axes groups

`--cohort`, `--group-by COL`, `--groups a,b,...` (order, defaulting to
sorted observed values), `--score-cols`, `--median-cols "age,QIgG/QAlb"`,
`--reference GROUP`, `--ci`, `--output` (`group, n, <expr>_median…,
<score>_mean, <score>_ci_lo, <score>_ci_hi`; t-based CI, n ≥ 3),
`--output-tests` (`kind, score, reference, group, n_ref, n_group, statistic,
p`). `kruskal` rows carry H and its chi-square p. `reference_vs_group` rows
carry pooled-SD Cohen d of reference minus group, and the Welch p.

#### atman axes anchor

`--cohorts`, `--score-cols`, `--anchors "log10(QAlb),log10(QIgG/QAlb),age"`,
`--scope all|controls|cases|condition=X` (repeatable), `--partial EXPR`,
`--n-bootstrap`, `--seed`, `--ci`, `--min-n` (default 10), `--output`
(`cohort, scope, given, score, anchor, n, rho, ci_lo, ci_hi, p`).
Spearman runs on pairwise-complete subjects. `p` comes from the t
approximation. The CI comes from a subject bootstrap, one stream per
output row, `derive_sub_seed(seed, row_index)`.
`--partial` adds `<scope>_partial` rows: Pearson correlation of the rank
residuals of score and anchor on the ranks of the given expression, no CI.

#### atman axes displacement

`--manifest` (contrast manifest), `--score-sets "4axes=axis1_z,...;13arch=primary_A0001,..."`,
`--standardize global|none` (default `global`: each column z-scored over
every subject in the table before centroids), `--n-bootstrap`, `--seed`,
`--ci`, `--output` (`a, b, space, cosine, ci_lo, ci_hi, norm_a, norm_b`),
`--output-vectors` (`label, space, score, displacement, n_case, n_control`).
Displacement = case centroid − control centroid within the contrast's cohort,
column-wise over subjects with a value. The bootstrap resamples case and
control subjects of both contrasts independently, one stream per pair.

#### atman axes loco

`--score-cols axis1_raw,...` (raw columns, with z recomputed), `--groups
"Controls=healthy|NoLeak|control|CU|nonMS;AD=AD;..."`, `--output` (`score,
group, dropped_cohort, n_remaining, full_z_mean, loo_z_mean, delta`). Cohorts
are dropped in order of first appearance. A group absent from the remaining
subjects yields no row.

#### atman axes icc

`--scores`, `--cohort-dirs` (supplies `subject_id`) or `--subject-col COL`,
`--cohorts`, `--score-cols`, `--output` (`cohort, score, n_subjects,
n_samples, k_mean, k0, n_singletons_excluded, icc1, between_sd, within_sd`).
One-way random-effects ICC(1): `icc1 = (MSB − MSW) / (MSB + (k0 − 1) MSW)`
with the unbalanced `k0 = (N − Σk_i²/N)/(n − 1)`. Then `between_sd =
sqrt(max(0, (MSB − MSW)/k0))` and `within_sd = sqrt(MSW)`. Subjects with a
single scored sample are excluded and counted.

#### atman axes tree

`--scores`, `--manifest` (each contrast a leaf), `--score-cols`, `--leaf
displacement|centroid` (case − control with `1 − cosine`, or the case-group
centroid with Euclidean distance, unless `--distance` overrides),
`--standardize global|none`, `--n-bootstrap` (default 1000), `--seed`, `--quote-labels`,
`--output-linkage` (`left, right, distance, n`, scipy convention: leaves
`0..n`, merge `i` is node `n+i`, height = merge distance, `left < right`),
`--output-support` (`node` = sorted leaf labels joined by ` | `, `n_leaves`,
`support`, `node_id`, `height`, `n_boot`), `--output-newick` (branch length
= parent height − child height, six decimals, and internal labels are
integer-percent support). Support is the fraction of subject-bootstrap
trees (resampling within each contrast's case and control groups) that
contain the reference clade.
