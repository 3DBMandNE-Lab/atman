# Atman feature requests

Open batch for the cross-cohort decomposition / alignment methods track.
The DE surface is mature (paired-t, welch, OLS-formula, mixed, limma
+ DEqMS, msqrob with Huber robust, ensemble consensus dispatcher) and
every post-hoc family for multi-level OLS (Sidak, Tukey HSD via
from-scratch studentized range, Dunnett via equicorrelated
multivariate-t). The next frontier for atman as a standalone methods
tool splits into two tracks: (i) rigor around the existing
`decompose ica` + `align programs` pair — projection onto a trained
atlas, benchmark harness, data-driven module discovery; and (ii) a
cross-domain port of geometric compartmental unmixing (VCA + FCLS)
from hyperspectral remote sensing, giving interpretable per-subject
compartment fractions where ICA gives abstract axes.

All items below keep commandment 8 (Rust, SQLite, local, deterministic,
no cloud except PubMed API) and the canonical TSV contract.

Open priorities — all shipped 2026-04-20:

- ~~Priority 4~~ — `atman align project` ✅ (commit `0337f15`)
- ~~Priority 6~~ — `atman bench decompose` ✅ (commit `e46bbed`)
- ~~Priority 8~~ — `atman modules discover` ✅ (commit `92961c4`)
- ~~Priority 10~~ — `atman decompose unmix` ✅ (commit `d61706c`)

Residual deferrals — **all cleared 2026-04-20**:

- ~~Jaccard/Spearman metrics in `align bootstrap`~~ ✅ (commit
  `cf271c5`) — BootstrapParams + CLI now dispatch on metric.
- ~~Unbalanced Dunnett–Hsu~~ ✅ (commit `968a25c`) — MC
  multivariate-t with per-pair correlation matrix, auto-detected
  from observed n_i when `max/min > 1.25`.
- ~~NFINDR endmember extraction~~ ✅ (commit `93a87f4`) —
  iterative simplex-volume max over VCA init.
- ~~`--k auto` HySime-style selection~~ ✅ (commit `b2533ba`) —
  residual-elbow sweep + `k_selection.tsv`.
- ~~`--n-boot` bootstrap CI~~ ✅ (commit `5ed5bdb`) — subject-
  level resample, matched+signed loadings, percentile CI on both
  loadings and abundances.
- ~~`--annotate-markers` ORA~~ ✅ (commit `aea7196`) —
  hypergeometric upper-tail on top-N proteins.

Only remaining deferral: the CSF plasma-endmember *real-data*
sanity check on `sih/qc_measurements.tsv`. That one requires the
live `atman_inputs_albnorm/` dataset and lives downstream of
atman's test surface.

Also shipped 2026-04-20 outside the original request batch:
**Cross-method consensus DE dispatcher (`atman de --test ensemble`)** —
runs every applicable DE method on the same inputs and emits a
per-protein VALIDATED / PROVISIONAL / INSUFFICIENT grade with Stouffer
ensemble p/q.

---

## ~~Priority 1: Archetype null calibration (`atman decompose null`)~~ *(shipped 2026-04-20)*

Permutation null for FastICA archetype stability with three null modes
(sample-shuffle, protein-shuffle, gaussian-matched), BH-adjusted `null_q`,
and a `decision` column. Integration tests cover schema, determinism,
k-size refusal, and an end-to-end smoke test against the Dube cohort. See
`docs/analytical-roadmap.md` §8 for details.

---

## ~~Priority 2: Compositional transforms in `atman decompose ica`~~ *(shipped 2026-04-20)*

Exposed `--transform {none,clr,alr,ilr,ratio-anchor}` with
`--alr-reference <gene>` and a `transform_applied.json` audit file
next to `loadings.tsv`. Pure transforms live in
`atman-core::compositional`. Zero-handling deferred until linear-scale
input routes exist; atman canonical data is already finite on a log
scale after QC. See `docs/analytical-roadmap.md` §8 for details.

---

## ~~Priority 3: Bootstrap alignment uncertainty (`atman align bootstrap`)~~ *(shipped 2026-04-20)*

Subject-level bootstrap that resamples subjects within every cohort
with replacement, re-runs FastICA per cohort, re-aligns with cosine
similarity, and emits per-PE-archetype `bootstrap_prob_universal`,
`bootstrap_prob_multi`, percentile CI band, Shannon entropy over
the n_cohorts histogram, and a BCa CI whose acceleration is
estimated by pooled subject-level jackknife. Deterministic
SplitMix64`(seed, iter)` sub-seeds. See
`docs/analytical-roadmap.md` §8 for details.

**Retired in DEBT-4 (commit `e418c07`, 2026-04-20):** BCa CI,
`alignment_entropy`, and planted-fixture magnitude assertion on
`bootstrap_prob_multi ≥ 0.7 ∧ entropy < 1.0 ∧ bca_upper ≥ 2` for
the universal archetype.

**Still deferred:**

- Jaccard and Spearman metrics (v1 is cosine-only).

---

## Priority 4: Cohort projection (`atman align project`)

Given a trained atlas (archetypes from cohorts A/B/C), project a new
cohort D into that archetype space without re-running the full
alignment. Essential for atlas / consortium workflows: an external site
contributes a new dataset, you tell them how their subjects score on the
atlas without rebuilding the alignment from scratch.

Command:

```bash
atman align project \
  --atlas-dir out_align/atlas \
  --cohort-dir out_new_cohort \
  --transform clr \
  --output-dir out_new_cohort/projected
```

`--atlas-dir` is expected to contain the outputs of a prior `atman align
programs` run: the canonical archetype loading matrix (one row per
archetype, one column per protein), the chosen `--transform`, and the
intersection protein universe.

Behavior:

1. Load the atlas archetype loadings and protein universe.
2. Load the new cohort's protein matrix, restrict to the atlas universe,
   apply the same transform stored in the atlas sidecar.
3. Solve for per-subject activations via a non-negative or ridge
   least-squares projection onto the atlas loading matrix (configurable
   via `--projection ls|nnls|ridge`, default `ridge` with
   `--ridge-lambda 0.01`).

Outputs:

- `projected_activations.tsv`: one row per subject, one column per
  atlas archetype
- `projection_qc.tsv`: per-subject residual norm, number of atlas
  proteins present, number missing, per-archetype activation CI from
  bootstrap (optional with `--n-boot`)
- run sidecar capturing atlas SHA-256, transform, projection method,
  ridge-lambda, etc.

Acceptance:

- Projecting a cohort's own training data back onto the atlas it
  contributed to must recover activations within numerical tolerance
  of the original `decompose ica` activations (sanity check).
- A cohort with only partial overlap in the protein universe must emit
  a clear warning in `projection_qc.tsv` and set per-archetype
  `coverage_fraction`.
- Integration test under `crates/atman/tests/align_project.rs` covering
  round-trip sanity, partial-coverage warnings, and mismatched-transform
  refusal.

Why it matters for the methods paper: the difference between a one-shot
alignment experiment and a reusable resource. Lets the CSF atlas be
queried by new datasets without re-running everything. Combined with
Priority 3, it also gives projected-archetype activations proper
uncertainty intervals.

---

## ~~Priority 5: Archetype variance decomposition (`atman decompose variance`)~~ *(shipped 2026-04-20)*

Mixed-model variance partition per archetype reusing the
`atman de --test mixed` REML engine. Accepts formulas like
`"cohort + condition + (1|subject_id)"`, outputs per-factor
Type III sums of squares, F-statistic, df, p-value, random-
intercept variance, ICC, and per-coefficient Wald summaries.
See `docs/analytical-roadmap.md` §8 for details.

**Retired in DEBT-2 (commit `c81e5cc`, 2026-04-20):** Type III
SS (Wald-form) and per-factor omnibus F-tests, with parity to
R's `car::Anova(type = 3)` at ≥ 6 decimals on a 48-sample 3-factor
balanced fixture (max |Δss| = 1.9e-6, max |ΔF| = 1.5e-5,
max |Δp| = 3.7e-7).

**Still deferred:** none.

---

## Priority 6: Benchmark harness (`atman bench decompose`)

Head-to-head recovery benchmark of `atman decompose ica` against named
reference tools, on a shared synthetic fixture with planted archetypes.
Needed for a methods-paper figure that compares atman to MOFA, consICA,
scICA, or fastICA+ICASSO reference implementations.

Command:

```bash
atman bench decompose \
  --fixture bench/planted_archetypes_v1 \
  --tools atman,fastica-icasso,mofa,consica \
  --metric recovery-jaccard,runtime,memory \
  --seed 20260418 \
  --output bench_results.tsv
```

Scope:

- Fixture is a pure TSV package shipped under `bench/` with a
  ground-truth archetype loading matrix, a synthetic sample-protein
  matrix generated from it plus calibrated noise, and a planted
  cross-cohort structure (two or three synthetic "cohorts" with a
  shared universal archetype + cohort-specific archetypes).
- Atman runs its own `decompose ica` + `align programs` natively.
- Other tools are invoked via thin adapter shells in `bench/adapters/`
  (one shell script per tool, wrapping the reference implementation
  from R/Python with a TSV-in/TSV-out contract). No Python on the
  atman analytical path; the adapters are outside atman's own runtime
  surface.
- `atman bench` only orchestrates and scores; it does not require the
  reference tools at test time. Missing tools emit a clear
  `tool_not_available` row rather than failing the run.

Metrics:

- `recovery_jaccard`: Jaccard overlap between recovered and planted
  top-N loading sets per archetype (N=50)
- `archetype_correlation`: Pearson between recovered and planted
  loading vectors (best-matching pair via reciprocal-best cosine)
- `runtime_seconds`, `peak_memory_mb`
- `determinism_score`: byte-equality of outputs under repeated
  identical seed/input

Acceptance:

- `atman` rows always populate (the tool can always score itself).
- Fixture committed under `bench/planted_archetypes_v1/` with its own
  README explaining the generation process and seed.
- Unit test covering the scoring math on a hand-constructed tiny
  recovered-vs-planted pair.

Why it matters: a methods paper without a benchmark figure is a tool
paper. With it, atman can be cited on "equal or better recovery at
fraction-of-the-runtime, with determinism reference tools can't match."

---

## ~~Priority 7: Network influence (`atman network influence`)~~ *(shipped 2026-04-20)*

Feature-covariance hub scoring via eigenvector × betweenness
centrality. New top-level command family (`atman network …`);
subject-level similarity with pearson/spearman/covariance,
hard-threshold or WGCNA-style soft-power adjacency, optional
stratification. See `docs/analytical-roadmap.md` §8 for details.

Spec below preserved for reference:

---

## Priority 7 (archived spec): Network influence (`atman network influence`)

Score each feature by its role as a hub in the subject-level covariance
graph. `atman align programs` tells you which archetypes recur across
cohorts; `atman network influence` answers the orthogonal question of
which individual proteins are the scaffolding the archetypes rest on.
Motivated by Burberry et al. 2026 (bioRxiv
doi:10.64898/2026.04.02.716122), who use the same construction on
CyTOF features to identify signaling hubs that organize the
cross-compartment immune network in AD.

Command:

```bash
atman network influence \
  --input out/qc_measurements.tsv \
  --samples out/samples.tsv \
  --method spearman \
  --threshold 0.3 \
  --output out/network_influence.tsv
```

Graph construction:

- `--method pearson|spearman|covariance`: feature × feature similarity
  computed over subjects. Reuse `atman-core::stats` primitives.
- `--threshold <f64>`: hard threshold on `|similarity|`. Below, edge is
  dropped.
- `--soft-power <u32>`: alternative to hard threshold. WGCNA-style
  `|similarity|^β` adjacency; preserves weighted edges. Mutually
  exclusive with `--threshold`.
- `--directed/--undirected`: default undirected.
- `--stratify <column>`: optional sample-metadata column (e.g.
  `cohort`, `condition`, `sex`). Builds one graph per stratum and emits
  one output row per (feature, stratum). Lets users ask "is this hub
  universal or stratum-specific?"

Node scoring:

- `eigenvector_centrality`: power-iteration solve on the weighted
  adjacency; converged to a fixed tolerance with a deterministic
  starting vector.
- `betweenness_centrality`: Brandes' algorithm on the thresholded
  graph; edge weights interpreted as reciprocal distances under
  `--betweenness-weighting reciprocal` (default).
- `influence_score = eigenvector_centrality × betweenness_centrality`
  (Burberry-Pillai construction).
- `degree`: number of retained edges.

Outputs:

- `network_influence.tsv`: `feature_id, stratum,
  eigenvector_centrality, betweenness_centrality, influence_score,
  degree, n_subjects_used`.
- `network_edges.tsv`: `feature_i, feature_j, edge_weight, stratum`
  (one row per retained edge). Emitted only when `--emit-edges` is
  set — can be large.
- run sidecar capturing method, threshold/soft-power, stratification,
  SHA-256 of input measurements, and any randomized-tiebreaker seed.

Acceptance:

- Deterministic under identical input + flags. Any tiebreaking in
  centrality iteration uses a fixed deterministic rule; repeated runs
  byte-equal.
- On a synthetic star-topology fixture (1 hub + 5 leaves, no other
  edges), the hub's `eigenvector_centrality` must be within numerical
  tolerance of 1.0 and its `betweenness_centrality` must strictly
  exceed every leaf's.
- On a two-block fixture (10 correlated features + 10 uncorrelated
  features), no uncorrelated feature may appear in the top-5 by
  `influence_score`.
- Stratified mode produces disjoint per-stratum outputs; the union
  must match non-stratified output when `--stratify` is dropped.
- Integration test under `crates/atman/tests/network_influence.rs`
  covering star-topology recovery, two-block separation, stratified
  mode, determinism, and clear failure when the graph is empty
  (threshold too high / no surviving edges).

Why it matters for the methods paper: gives atman a first-class
network-hub analysis that no standalone proteomics tool currently
exposes. Directly answers "which proteins are organizing the covariance
structure the archetypes are built on?" Complements `align programs`
rather than duplicating it.

---

## Priority 8: Data-driven module discovery (`atman modules discover`)

Atman currently *scores* user-provided modules (`score modules`,
`bootstrap module`, `module-de`). It cannot *find* them. This is the
main structural gap vs WGCNA / MOFA. Close it with a soft-thresholded
adjacency + topological overlap clustering command that emits a
canonical `modules.tsv` consumable by the existing scoring surface.
Motivated by Burberry et al. 2026 (bioRxiv
doi:10.64898/2026.04.02.716122), who use WGCNA-inspired module
construction to summarize correlated immune features into biologically
coherent modules per compartment.

Command:

```bash
atman modules discover \
  --input out/qc_measurements.tsv \
  --samples out/samples.tsv \
  --method wgcna-soft \
  --soft-power auto \
  --min-module-size 5 \
  --seed 20260418 \
  --output out/modules_discovered.tsv
```

Methods:

- `wgcna-soft`: soft-threshold `|correlation|` at power β. If
  `--soft-power auto`, sweep β ∈ {1, …, 20} and pick the smallest
  integer satisfying the scale-free topology criterion (R² ≥ 0.8,
  slope < 0). Compute topological overlap matrix (TOM); cluster
  `1 - TOM` with average linkage; dynamic tree cut for module
  assignment. Unassigned features get module `grey` (WGCNA
  convention).
- `hard-threshold`: binary adjacency at `|r| ≥ threshold`; connected
  components become modules. Simpler, faster, less robust to noise.
- `consensus`: resample subjects with replacement, re-cluster each
  resample, aggregate the co-assignment matrix, hierarchical-cluster
  on `1 - co_assignment_freq`. Slower but calibrates module
  stability. Requires `--n-resamples` (default 100).

Outputs:

- `modules_discovered.tsv`: schema matches the existing user-supplied
  `modules.tsv` contract (`module, gene_symbol`) so it is immediately
  consumable by `score modules`, `bootstrap module`, and `module-de`.
- `module_discovery_report.tsv`: per-module summary — `module, size,
  mean_within_|r|, hub_feature, eigenprotein_pc1_variance_explained`.
- `soft_power_diagnostics.tsv` (wgcna-soft only): β vs scale-free R²
  sweep, final β, sft_index.
- `consensus_stability.tsv` (consensus only): per-feature mean
  co-assignment frequency with its final module.
- run sidecar capturing method, power, min-size, seed.

Acceptance:

- Deterministic under `--seed`. Consensus-clustering resamples draw
  from the in-tree Xoshiro256++ via a fixed iteration-to-seed
  derivation rule.
- On a synthetic two-block fixture (two planted blocks of 20
  correlated features + 40 uncorrelated background features),
  `wgcna-soft` must recover two non-grey modules of size ≥ 18 each and
  zero spurious non-grey modules above `--min-module-size`.
- Round-trip sanity: `atman score modules` on the discovery output,
  scoring the training data, must produce per-subject module scores
  that correlate > 0.9 with the PC1 eigenprotein of the same module
  features (the WGCNA-canonical "module eigengene" definition).
- Refuses when no module meets `--min-module-size` with a clear error.
- Integration test under `crates/atman/tests/modules_discover.rs`
  covering two-block recovery, eigenprotein round-trip, soft-power
  auto-selection, consensus-stability emission, and
  insufficient-module-size refusal.

Why it matters for the methods paper: converts atman from "module
scoring tool that assumes you already have modules" to "module
discovery + scoring pipeline." Matches the WGCNA surface every
co-expression methods paper compares against, with atman's determinism
and canonical-TSV contract on top.

---

## ~~Priority 9: Post-hoc contrasts for multi-level OLS (`atman de --post-hoc`)~~ *(shipped 2026-04-20)*

**Shipped:**

- Per-protein omnibus F-test via `atman de --test ols
  --omnibus-factor <name>`. Output `de_omnibus.tsv` carries
  `f_statistic, df_num, df_den, p_value, bh_q` with BH-adjustment
  within each (comparison, panel) family.
- Sidak on a user-supplied `--contrast-list` — DEBT-1 (`2f9f61a`),
  parity to R `lm` + linear contrasts + Sidak at 1e-6 on a
  36-sample 3-level fixture.
- Tukey HSD via from-scratch studentized range distribution —
  DEBT-5 (`1ac2e06`), nested 128×48 Gauss–Legendre; parity to R's
  `ptukey` at ≥ 3 decimals on α = 0.05 critical values for
  `k ∈ {3, 4, 5}`.
- Dunnett via equicorrelated multivariate-t — DEBT-6 (`1c58d78`),
  Dunnett–Curnow one-variate-plus-idiosyncratic reduction to a 2D
  quadrature at `ρ = 0.5`; parity to R `qmvt`-based critical values
  at ≤ 1e-2.

See `docs/analytical-roadmap.md` §8 for details.

**Still deferred:** unbalanced Dunnett–Hsu (heterogeneous per-pair
correlations via Genz–Bretz on an arbitrary correlation matrix).

Spec below preserved for future reference:

---

`atman de --test ols` fits the linear model and emits a single contrast
specified via `--contrast`. For designs with ≥3 levels of a categorical
factor (e.g. CN vs CN-Aβ+ vs MCI vs AD; placebo vs low vs high dose),
users currently stitch together multiple pairwise DEs — which
double-adjusts across the protein axis, loses family-wise error control
across the contrast axis, and drops the omnibus F-test. Extend OLS /
mixed to emit all pairwise contrasts with proper FWER correction.

Command:

```bash
atman de \
  --input-dir out \
  --output-dir out_de_posthoc \
  --test ols \
  --design "~ stage + age + sex" \
  --factor stage \
  --post-hoc tukey \
  --alpha 0.05 \
  --min-pairs 5
```

Post-hoc methods:

- `tukey`: all-pairs Tukey HSD on `--factor`; studentized-range
  distribution; controls FWER at `--alpha`.
- `dunnett`: each non-reference level vs the reference level of
  `--factor`; one-step multivariate-t; controls FWER. Reference level
  configurable via `--reference-level`, default = alphabetically-first
  or `is_control == 1`.
- `sidak`: user-supplied contrast list via `--contrast-list
  "MCI-CN,AD-CN,AD-MCI"`; Sidak-adjusted across the specified list
  only.
- FDR across proteins (`--fdr bh`, already implemented) is orthogonal
  and still applies.

Behavior:

- Fit the OLS once per protein (as today).
- For `tukey|dunnett`, derive all defined contrasts on the fitted
  model automatically from the factor's levels.
- For `sidak`, parse `--contrast-list` and refuse unknown levels.
- Emit one `de_results` row per (protein × contrast).
- Also emit the omnibus F-test per protein.

Outputs:

- `de_results.tsv` gains columns: `contrast` (e.g. `"MCI-CN"`),
  `posthoc_method`, `posthoc_p`, `posthoc_adj_p` (family-wise). Empty
  when `--post-hoc` is not requested.
- `de_omnibus.tsv`: one row per protein with the `--factor` omnibus
  F-statistic, df, and p-value.
- `de_report.tsv` gains per-contrast counts at `posthoc_adj_p < alpha`
  and omnibus hit counts at `omnibus_q < 0.05`.

Acceptance:

- Tukey HSD results match R's `TukeyHSD(aov(...))` or
  `emmeans::pairs(lm_fit)` to 4 decimal places on a balanced
  three-level fixture.
- Dunnett results match `emmeans::contrast(method="dunnett")` to 4
  decimal places.
- Refuses with a clear error when `--post-hoc` is requested with
  `--factor` having fewer than 3 levels.
- Refuses when `--post-hoc` is combined with `--test welch-t` or
  `--test paired-t` (OLS / mixed only).
- Integration tests under `crates/atman/tests/de_posthoc.rs`:
  Tukey / Dunnett / Sidak parity to a reference-R golden file,
  2-level refusal, non-OLS-test refusal, omnibus F emission, and a
  covariate design preserving the factor's post-hoc contrasts after
  adjustment.

Why it matters for the methods paper: removes one of the most common
reasons users fall off the canonical `atman de` path for multi-level
designs. Burberry et al. 2026 (bioRxiv doi:10.64898/2026.04.02.716122)
use exactly this combination (two-way ANOVA + Tukey + Dunnett + Sidak)
as the workhorse for their four-stage AD continuum design; the CSF
cross-disease manuscript has the same shape (CN vs AD vs iNPH vs MS vs
SIH vs PD) and currently hand-codes it in Python.

---

## ~~Priority 10: Compartmental unmixing (`atman decompose unmix`)~~ *(shipped 2026-04-20, commit `d61706c`)*

**Shipped:** VCA endmember extraction + FCLS / UCLS abundance
estimation in `atman-core::decompose_unmix`; CLI at
`atman decompose unmix` with `--transform {none,log,clr,alr,
ratio-anchor}` (CLR/ALR + FCLS gated on
`--allow-unconstrained-simplex`; ILR refused). Integration tests
verify planted-3-endmember recovery at cosine ≥ 0.85, FCLS row-sum
= 1 ± 1e-5 with non-negative entries, determinism, and the `k > n/2`
+ CLR+FCLS refusal paths.

**Still deferred:** NFINDR endmember extraction, `--k auto`
(HySime-style virtual dimensionality), `--n-boot` subject-level
bootstrap CI for loadings + abundances, `--annotate-markers` with
ORA against supplied marker sets, and the CSF plasma-endmember
real-data sanity assertion on `sih/qc_measurements.tsv` (requires
the live atman_inputs_albnorm dataset).

Spec below preserved for future reference:

---

## Priority 10: Compartmental unmixing (`atman decompose unmix`)

Atman's `decompose ica` returns statistically-independent axes; their
biological interpretation is a post-hoc curation step. Hyperspectral
remote sensing solved the same "each sample is a weighted mixture of a
few pure endmember signatures" problem geometrically fifty years ago.
Port Vertex Component Analysis (Nascimento & Bioucas-Dias 2005) for
endmember extraction and Fully Constrained Least Squares (Heinz & Chang
2001) for abundance estimation to give atman a command that returns
interpretable *compartment fractions* — each subject lands as
`(α_endmember_1, α_endmember_2, …)` summing to 1 under non-negativity,
directly readable as compartmental composition without naming the
compartments a priori. Complements `decompose ica` rather than
replacing it: ICA gives statistically-independent abstract axes,
unmixing gives geometrically-identified pure endmembers.

Command:

```bash
atman decompose unmix \
  --input out/qc_measurements.tsv \
  --samples out/samples.tsv \
  --k auto \
  --method vca \
  --abundance fcls \
  --transform log \
  --seed 20260418 \
  --output-dir out_unmix
```

Methods:

- `--method vca|nfindr`: geometric endmember extraction. `vca`
  (default) iteratively finds extreme points by projecting onto the
  orthogonal complement of already-found endmembers. `nfindr` maximizes
  simplex volume spanned by a candidate endmember set.
- `--abundance fcls|ucls`: per-sample abundance estimation. `fcls`
  (default) solves `min ||x - Eα||²` subject to `α ≥ 0 ∧ Σα = 1`
  (non-negative + sum-to-one). `ucls` drops the simplex constraint for
  cases where sum-to-one is not meaningful (e.g. CLR-transformed
  input); returns unconstrained least-squares abundances.
- `--k <u32>|auto`: number of endmembers. `auto` sweeps
  `k ∈ {2, …, 10}` and picks the smallest k where reconstruction
  residual norm levels off (HySime-style virtual-dimensionality
  criterion adapted for proteomics) or where the permutation-null
  residual distribution no longer distinguishes the fit from noise
  (reuses the Priority 1 null machinery).
- `--transform log|clr|none`: pre-transform. `log` (default) for MS
  intensity, `clr` for NPX / rank-quantile inputs, `none` for inputs
  already on a log scale.
- `--n-boot <u32>`: optional subject-level bootstrap for per-endmember
  loading CI and per-abundance CI. Reuses the same
  `Xoshiro256++(seed, iter)` sub-seed scheme as `align bootstrap`.
- `--annotate-top <N>`: auto-annotate each endmember with top-N
  proteins by `|loading|`. If `--annotate-markers <markers.tsv>` is
  supplied, also runs ORA against the supplied marker sets and emits
  enrichment p-values.

Outputs:

- `endmembers.tsv`: `endmember_id, protein_id, gene_symbol, loading,
  loading_ci_lower, loading_ci_upper, rank_in_endmember`. One row per
  (endmember × protein); CI columns filled only when `--n-boot > 0`.
- `abundances.tsv`: `sample_id, subject_id, endmember_id, abundance,
  abundance_ci_lower, abundance_ci_upper`. Row-sums over `endmember_id`
  equal 1.0 within 1e-6 under `--abundance fcls`.
- `unmix_diagnostics.tsv`: per-sample `reconstruction_residual_norm,
  n_active_constraints, coverage_fraction` (fraction of endmember
  universe present in the sample).
- `endmember_annotations.tsv`: per-endmember top-N proteins plus ORA
  hits when `--annotate-markers` is supplied.
- `k_selection.tsv`: emitted when `--k auto`, reporting per-candidate-k
  residual norm and chosen k.
- run sidecar capturing `seed, k, method, abundance, transform,
  n_boot, input_sha256`.

Acceptance:

- Deterministic under identical `--seed`. VCA's initial projection
  direction and any tiebreaking use a fixed rule derived from the seed;
  repeated runs byte-equal under byte-equal inputs.
- Synthetic recovery fixture (3 planted pure endmembers × 50 proteins,
  100 synthetic subjects drawn from Dirichlet-mixed abundances plus
  Gaussian noise): `atman decompose unmix --k 3 --method vca` must
  recover endmembers at pairwise cosine ≥ 0.9 with the planted
  endmembers (after best-matching permutation) and reconstruct
  abundances at Pearson correlation ≥ 0.95 with the planted Dirichlet
  draws.
- Compositional constraints under `--abundance fcls`: every abundance
  must satisfy `abundance ≥ -1e-9`, every per-sample row-sum must equal
  `1.0 ± 1e-6`, and no sample may have more than `k - 1` abundances
  simultaneously at the zero constraint active-set (guards against
  pathological rank collapse).
- k auto-selection: on the planted-3 fixture, `--k auto` must select 3;
  on a planted-5 fixture, `--k auto` must select 5.
- CSF plasma-endmember sanity: on
  `atman_inputs_albnorm/sih/qc_measurements.tsv`, at least one
  recovered endmember must dominate-load a supplied plasma-protein
  marker set (e.g. `{ALB, IGHG1, TF, HP, APOA1}`) at cosine ≥ 0.5 with
  the marker-set indicator vector.
- Refuses with a clear error when `k > n_samples / 2`
  (under-determined).
- Refuses with a clear error when `--abundance fcls` is combined with a
  transform that does not admit a simplex interpretation (`clr`, `ilr`)
  unless `--allow-unconstrained-simplex` is set.
- Integration test under `crates/atman/tests/decompose_unmix.rs`
  covering: synthetic recovery at k = 3, compositional-constraint
  verification under FCLS, k-auto selection on planted-3 and planted-5
  fixtures, CSF plasma-endmember sanity on SIH, determinism across
  repeat runs, under-determined refusal, and
  transform-incompatibility refusal.

Why it matters for the methods paper: gives atman a first-class
interpretable compartmental-decomposition command that no standalone
proteomics tool currently exposes. ICA finds statistically independent
abstract axes; unmixing finds *geometrically-identified pure
endmembers* with per-subject fractional abundances that sum to 1 — the
natural output when the generating model is a mixture of a few pure
sources, which is exactly the model behind any convoluted biological
fluid (CSF, plasma, urine, bile, pleural, ascitic, synovial, amniotic).
Imported directly from a fifty-year-old solved problem in hyperspectral
remote sensing (Nascimento & Bioucas-Dias 2005 for VCA; Heinz & Chang
2001 for FCLS). Paper angle: *"Geometric compartmental deconvolution of
biological fluids via hyperspectral-style endmember extraction,"* with
atman's determinism, canonical-TSV contract, and per-abundance
bootstrap CI as differentiators against any R/Python reference
implementation ported straight out of remote-sensing textbooks.

---

## Future directions (not yet spec'd)

The following are deliberately held back until Priorities 1-10 ship,
because each extends one of them and the priority-ordered design is
cleaner once the foundations are tested:

- **Missingness-aware ICA** (priority 11 candidate). Extend the
  Priority 2 compositional transforms with a joint abundance +
  detection likelihood at decomposition time, same philosophy as
  `atman de --test msqrob` or proDA but at the component level. Lets
  archetypes be defined partly by "which proteins were detectable in
  this sample" — which is real biology for MNAR-heavy proteomics.
- **Hierarchical / nested alignment** (priority 12 candidate). Site →
  cohort → cross-cohort alignment for consortium data where a single
  cohort has multiple acquisition sites. Requires Priority 3 bootstrap
  uncertainty to be in place so each level of the hierarchy carries
  its own confidence interval.
- **PMF (positive matrix factorization) with per-cell uncertainty
  weighting**, ported from atmospheric source apportionment (Paatero &
  Tapper 1994). Sibling to Priority 10: where VCA+FCLS finds geometric
  endmembers, PMF finds non-negative factors while natively weighting
  each cell by its reciprocal analytical uncertainty. Useful for
  DIA-MS inputs that ship per-protein-per-sample CV or detection-limit
  flags. Implementable as `atman decompose pmf --uncertainty
  measurements_cv.tsv`.
- **MCR-ALS (multivariate curve resolution — alternating least
  squares)** from chemometrics. Sibling to Priority 10: NMF with
  chemistry-aware constraints (non-negativity, unimodality, closure)
  applied to the protein-by-sample matrix. Implementable as
  `atman decompose mcr`.
- **Pooled-QC drift correction**, ported from untargeted metabolomics.
  Standard practice there; near-absent in proteomics. Fit a per-protein
  LOESS on injection order using interleaved pooled QC samples, correct
  analytical drift before DE or decomposition. Implementable as
  `atman qc drift --injection-order injection_order.tsv --pooled-qc
  qc_samples.tsv`.
- **Longitudinal tensor decomposition**. Samples × proteins ×
  timepoints. Relevant for MS PILOT, iNPH repeat measures, any
  consortium with longitudinal arms. Depends on Priority 5 variance
  decomposition since random-subject variance becomes the signal, not
  the nuisance.
- **Counterfactual archetype simulation**. `atman decompose
  counterfactual --set-archetype A0004 --to 0 --input activations.tsv`
  reconstructs a predicted abundance matrix with one archetype zeroed
  out. Interpretive aid; low-priority until Priorities 1-4 land.
- **Compositional effect size for DE**. A `--effect-size-scale clr`
  flag on `atman de` that reports CLR-space coefficients alongside the
  standard log-FC. Natural follow-on to Priority 2.
