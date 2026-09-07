# Changelog

All notable changes to Atman are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`decompose null` writes an `n_perm` column.** It is the denominator
  of `null_stability_mean`, `null_stability_p95` and `null_p`, and it
  travelled only in the run sidecar. `null_stability_p95` is the reason:
  its index is `ceil(0.95 * n_perm) - 1`, so at `--n-perm 20` it lands on
  the 19th of 20 sorted draws and one or two decide it, with nothing in
  the row saying how many there were. Exposure was limited rather than
  absent — the default is 200, where the p95 rests on roughly the top
  ten, and `decision` keys on `null_q` from the properly floored
  `null_p`, so nothing rides on the p95 — but a reader of the table does
  not open the sidecar. `n_perm` is NOT the denominator of `null_q`,
  which is BH across programs and bounded by the row count, and the
  reference says so. The column is appended last so a positional reader
  of the original seven keeps working. Found by applying the meningioma
  session's percentile-of-few-draws shape to atman's own estimators.

- **`enrich gsea` reports `n_same_sign_perms`**, the number of
  permutation draws the NES denominator was averaged over, as a new
  trailing column and on `GseaResult`. NES is `es / mean(|es_perm|)`
  over same-sign draws only, and nothing recorded how many that was.
  The relationship is inverted, which is what makes it dangerous: fewer
  same-sign draws give a LARGER NES and a less determined one. Measured
  at `es = 0.80` over 1000 draws, varying only the same-sign count: 1
  draw gives NES 40.00 at p 0.50, 3 gives 38.10 at p 0.25, 10 gives
  32.65, 50 gives 17.98, 500 gives 2.97 at p 0.002. The biggest scores
  in a run are the least determined, and the p-value moves the opposite
  way because its denominator is `n_same_sign_perms + 1`. The command
  now warns below 19 draws, which is the point under which the p-value
  floor of `1 / (n + 1)` prevents `p <= 0.05` however extreme the ES.
  The normalisation is NOT changed and stays faithful to the reference
  method — the fgsea parity test passes unaltered. This reports how well
  determined a result is, it does not alter the result. The column is
  appended last so a positional reader of the existing seven columns
  keeps working. Extreme sign asymmetry turns out to be rare with real
  rankings (a perfectly top-loaded fixture still had 353 same-sign draws
  of 500), so the common cause is a low `--n-permutations`, which caps
  the denominator directly. Found by applying the CSF session's null
  calibration rule to atman's own estimators.

- **`harmonize apply` now reads the model file strictly**, and the
  CLI layer has tests. Every field the fit/apply contract rests on was
  read with `unwrap_or_default`, so the two safety claims FAILED OPEN on
  a truncated or hand-edited file.

  An absent `fit_cohorts` became an empty list. The held-out refusal is
  a membership test against that list, so an empty one matched nothing
  and permitted every cohort, including the ones the model was fit on.
  **The command then exits 0 and reports a held-out evaluation that is
  not one.** The failure preserves every surface property a user
  checks — exit code, output shape, subject count, score distribution
  are all normal. Only the guarantee is gone, and nothing in the output
  records its absence.

  An absent `permuted_labels` became `false`, so a negative-control
  model reported itself as a real-labels one in the log line and the run
  sidecar. That is the one arm whose entire purpose is to be
  distinguishable from the real one.

  Both failures were silent.
  `apply` now names the missing field and exits non-zero, and also
  refuses an empty feature axis, a non-numeric or non-finite direction
  entry, a method whose fitted state is absent or empty, and an unknown
  `schema_version`. It warns when the file records no
  `fit_inputs_sha256`. A model written by `harmonize fit` is unaffected.
  Ten CLI tests were added, against a layer that previously had none.
  Six of them fail against the old reader, which was verified rather
  than assumed.

- **`align bootstrap` now names the population behind every summary
  column**, on stderr and in the run sidecar. The percentile CI and
  `bootstrap_mean_n_cohorts` are computed over matched replicates only,
  while `bootstrap_prob_universal`, `bootstrap_prob_multi`,
  `alignment_entropy` and the BCa bounds run over all `n_boot` with each
  miss entered as zero. Both are correct and neither name said which.
  Found by the GBM manuscript session: an NMF archetype at
  `prob_universal` 0.205 and an ICA archetype at 0.805 reported the
  IDENTICAL interval `[6, 6]` over six cohorts, a four-fold difference
  in reproducibility the interval could not show. The percentile CI is
  conditional on recovery, so an archetype matched twice in 200
  resamples, spanning six cohorts both times, reports `[6, 6]` — which
  makes a pre-registered `ci_lower >= k` criterion close to untestable.
  The command now reports the median match rate and warns below 0.50,
  and the run sidecar carries a `column_populations` block so a figure
  script reading the TSV months later can assert its sort key is
  unconditional. Column names and the TSV bytes are unchanged.

- **Two more defects in the same file, found by sweeping the match rate
  from 0.005 to 1.000.** `alignment_entropy` is NOT monotone in
  reproducibility: when an archetype spans the same cohorts whenever
  matched, it reduces to the binary entropy of the match rate, peaks at
  0.5 and falls to zero at BOTH ends, so a barely recovered archetype
  scores as the most stable. Measured at `n_boot` 200, a match rate of
  0.025 gives 0.169 against 0.732 at 0.205. The docstring had claimed
  low entropy implies stability, which inverts any ranking built on it.
  Separately, the BCa endpoints are degenerate at both tails: the lower
  bound reads the conditional value below a match rate of about 0.03,
  and the interval collapses to `[0, 0]` at 0.99.
  `bca_fallback_to_percentile` does not fire, because it tests the
  acceleration denominator. The BCa numerics are deliberately NOT
  changed — the reporting session declined the fix so its completed
  figure surface would not move — and are pinned in a test instead, so
  any future change produces a visible diff.

- **`harmonize apply --center-direction`**, and a warning when it is off
  for a per-sample method. Every per-sample normalisation leaves a
  per-subject term in the harmonised value that a direction with a
  non-zero mean projects onto. For `reference-protein` the algebra is
  exact — the score carries `−anchor × mean(direction)` — so if the
  anchor differs between arms in the held-out cohort, ANY direction
  separates them, including one learned from shuffled labels.
  Found by the CSF session's permuted arm on real data: holding out one
  cohort, reference-protein gave a real effect of +0.611 at p = 1.5e-05
  and AUC 0.70, while its permuted null ran +0.542 to +0.698 with mean
  AUC 0.73 — the shuffled arm separated cases *better* than the real
  one, and standardised against its own null the result was z = −0.66.
  That z is weakly determined and should be read as a direction rather
  than a magnitude. The reporting session has since measured its null
  spread: the reference-protein arm on that fold has a null standard
  deviation of 0.048 across eight permuted seeds, an order narrower than
  the widest cell in the same sweep at 0.399. A denominator estimated
  from eight draws that barely vary makes the ratio itself uncertain,
  and a narrow null inflates any standardised score. The finding the
  entry describes does not depend on the z: the permuted arm's raw
  separation exceeded the real arm's, which needs no denominator.
  Their measured null means order exactly as the mechanism predicts:
  reference-protein +0.267, quantile +0.193, rank +0.103, and the one
  per-protein method, zscore, at −0.103.
  Centring subtracts the direction's mean over each subject's own used
  features, which under MNAR dropout differs per subject, so a globally
  centred direction would still leave a residual. Off by default so
  earlier runs stay reproducible; the warning names the risk. What it
  costs is honest rather than hidden: a disease effect raising every
  protein uniformly is not separable from an anchor shift under a
  per-sample normalisation, so centring declines to claim unrecoverable
  signal rather than discarding recoverable signal.
- **`decompose ica --weighted-whitening`**: the detection model now
  reaches the estimator. Per-cell reliability scales each cell's
  contribution to the whitening moments — an observed cell counts fully,
  an undetected one counts by `1 − P(detected)` at its imputed value —
  while the contrast function stays unweighted, because a per-cell
  weight has a clean meaning in a covariance and none in a contrast.
  `atman_core::ica::pca_whiten_weighted` uses the unbiased
  reliability-weight denominator `Σw − Σw²/Σw`, so uniform weights
  reproduce `pca_whiten` exactly rather than differing by `n/(n−1)`.
  Off by default.
  Measured on the committed ground-truth fixture over six seeds, mean
  absolute error against the planted sources: mean-imputation 0.0768,
  the existing joint loop 0.0759, weighted whitening in that loop
  0.0754, and **weighted whitening with `--max-joint-iter 1` 0.0737**,
  best on every seed. So per-cell detection weighting helps, and the
  joint loop costs recovery: it drives imputed cells toward their own
  reconstruction, and running it undoes part of what the weighting
  gains. The recommended configuration is `--weighted-whitening
  --max-joint-iter 1`. All four orderings are asserted, not just
  reported.
- **`atman harmonize fit` / `atman harmonize apply`**: cross-cohort
  harmonisation with a fit-then-apply contract, for benchmarking
  harmonisation methods by held-out transfer. `fit` reads the training
  cohorts and writes a model; `apply` reads that model and one cohort
  and nothing else, so a held-out evaluation cannot leak by
  construction rather than by discipline, and applying a model to one of
  its own training cohorts is refused. The model carries the shared
  feature axis, the disease direction learned in that method's
  representation, the training cohort labels, and a SHA-256 of their
  inputs — labels alone would let a replay reuse a label over different
  data.
  Four methods: `zscore` (the null harmonisation), `rank`, `quantile`
  (reference profile learned on training cohorts), and
  `reference-protein` (per-sample division by the geometric mean of a
  low-variance reference set selected on training cohorts only, with
  `--reference-k`). `HarmonizeMethod::is_fitted` distinguishes the two
  stateless methods, because a stateless comparator matching a fitted
  one means the fitted state was not what carried the signal.
  `--permute-labels` shuffles case/control within each training cohort
  before learning the direction, giving a negative-control arm that
  travels the identical apply path. Every method scores something on a
  held-out cohort, so a benchmark whose entries all beat zero measures
  nothing. Note that ONE permuted fit is one draw from the null, not the
  null, and the flag's help says so: on a 60-feature fixture a single
  draw produced −4.45 against a real effect of +11.90. Use several
  seeds.
- **`decompose ica --degeneracy-floor`** (default 0.1) stops the MNAR
  joint loop when the imputed cells come within that fraction of their
  first-iteration RMS distance to their own rank-`k` reconstruction. The
  loop is iterated reconstruction-imputation, whose fixed point is the
  state where imputed cells equal their reconstruction exactly, carrying
  no independent information while the fit explains them by
  construction; measured on a 110-subject cohort the gap decayed
  geometrically at ratio ~0.992 per iteration with no floor, so running
  longer converged steadily toward that state and no criterion noticed.
  Stopping here is reported as `degeneracy-floor`, never as convergence:
  `MnarIcaResult` gains a three-valued `stop_reason` (`converged`,
  `degeneracy-floor`, `max-iterations`), recorded in the sidecar as
  `mnar_joint_stop_reason`, because `joint_converged` alone cannot
  distinguish outcomes that mean different things about the fit.
  `--degeneracy-floor 0` restores the previous behaviour.
- **`align programs --metric cosine-centered`, and a tau-selectivity
  report.** Plain cosine has a positivity floor on non-negative loadings:
  two NMF programs cannot score below zero and in practice sit high just
  from sharing that baseline. Measured across six CPTAC cohorts on a
  shared 4,375-gene universe, the median cross-cohort cosine was +0.813
  for NMF against −0.004 for ICA, so `--tau 0.30` admitted 99.9% of NMF
  pairs and 7.7% of ICA pairs. The same threshold was selective for one
  method and inert for the other, which makes a method comparison at a
  shared tau meaningless, and it also makes a tau sweep look robust when
  it is merely saturated. `cosine-centered` mean-centres each loading
  vector first, removing the floor rather than compensating for it
  (centred NMF median −0.002, 9.4% above 0.30, against ICA's 7.7%).
  `align programs` now also reports what fraction of cross-cohort pairs
  the supplied tau admits, and warns above 90% that the threshold is not
  thresholding and the structure is being decided by reciprocal-best
  matching alone. The choice of metric is recorded in the sidecar, since
  a reader cannot otherwise tell which similarity was computed.
  `align bootstrap` gets the same treatment, and needed it more: it
  applies the metric twice internally (`--cosine-tau` to group bootstrap
  programs into archetypes, `--match-tau` to match them back to the
  point estimate), both on the same loading vectors, so on non-negative
  loadings the bootstrap was deciding how often an archetype recurs
  through gates that admitted essentially everything. Both gates now
  report their admitted fraction and warn above 90%, and
  `cosine-centered` is documented in `--metric` there — it was already
  accepted by the parser but absent from the help text, so it could not
  be found.
  On what centring does to the result, now measured atman-to-atman with
  the metric as the only variable: almost nothing. NMF goes from 31
  archetypes to 29 with the same 3 multi-cohort and same single
  six-cohort archetype; ICA goes from 214 to 213 with 15 multi-cohort
  becoming 13. So the threshold was inert and the comparison was
  specified wrongly, and correcting it does not overturn the finding —
  the structure was never threshold-determined, reciprocal-best matching
  was doing the work. Both halves are worth stating: the comparison is
  now at comparable selectivity, and the result is insensitive to the
  choice.
- **`align project` now reports loading-weighted coverage per archetype**
  (`projection_archetype_coverage.tsv`, plus a stderr warning below 50%
  mass or fewer than 10 of the top 20 defining proteins). The existing
  per-sample `coverage_fraction` counts atlas proteins and weights them
  equally, so it stays high when hundreds of low-loading proteins are
  present and every defining one is absent — and absent proteins are
  filled with zero before the solve, so a missing defining feature
  actively contributes to the fit rather than merely going unused. Found
  in the wild: an archetype characterised as blood contamination was
  projected into a cohort where no haemoglobin is measured at all, and
  nothing in the output distinguished that from a faithful transfer.
  The new columns are weighted coverage, how many of the top 20
  loadings survive, and the largest absent loading as a share of the
  archetype's maximum.
- **MNAR joint trace now records the state, not only the readout.**
  `mnar_joint_trace.tsv` gains `imputation_delta` (largest change in an
  imputed cell against the previous round) and `reconstruction_gap`
  (RMS distance between imputed cells and their own ICA
  reconstruction). The existing `delta` column tracks the detection
  curve, which is refit every round but never re-enters the update —
  the inner ICA runs with uniform weights — so the iterated state is the
  imputed matrix and `delta` was reporting an observable rather than the
  thing that moves. `reconstruction_gap` exists because iterated
  reconstruction-imputation has a degenerate attractor: missing cells
  converge to exactly their rank-`k` reconstruction, at which point they
  carry no independent information and the fit explains them by
  construction. That gap approaching zero is a reason to stop, not a
  sign of health, and no previous output distinguished it from
  convergence.
- **`atman decompose ica --joint-tol`, and a per-iteration MNAR joint
  trace.** The outer joint loop's tolerance was hard-coded at 1e-6 and
  only its final verdict was reported, so a run that stopped at
  `--max-joint-iter` could not be described as closing, oscillating, or
  diverging — which is the difference between a reportable sentence and
  a bare failure. `--joint-tol` is now a flag, non-convergence warns with
  the achieved step size, and `mnar_joint_trace.tsv` is written beside
  the loadings with the detection-curve coefficients and step size at
  every joint iteration. The trace is a recorded output, so its hash is
  in the sidecar.
- **`decompose unmix --method spa`, now the default: deterministic
  endmember extraction.** The Successive Projection Algorithm (Araújo et
  al. 2001; robustness analysed by Gillis & Vavasis 2014) takes the
  sample of largest residual norm and projects it out, repeatedly. It
  uses no random direction, so the endmember set is a function of the
  data alone and `--seed` cannot reach it — verified on real data, where
  four different seeds give identical vertices on both matrices of a CSF
  cohort. `vca` and `nfindr` remain available under `--method`; both are
  seed-dependent on noisy data, and the sidecar now records
  `seed-affects-endmembers` so a reader does not need to know which is
  which. SPA also refuses loudly when `k` exceeds the dimensionality the
  data supports, instead of returning a vertex made of numerical noise:
  the relative residual falls by seven orders of magnitude at the true
  rank, and the threshold sits in that gap.
- **Real-data sanity test for `decompose unmix`** on the CSF CrossDisease
  SIH cohort (`crates/atman/tests/unmix_csf_plasma_endmember.rs`):
  albumin-ratio normalization should collapse the bipolar
  plasma↔neuronal axis, so a plasma endmember must appear on the primary
  matrix and must not appear on the albumin-normalized one. `#[ignore]`
  and gated on `ATMAN_CSF_CANONICAL_DIR`, since the inputs are
  gitignored and live outside the repo. Scores three disjoint marker
  panels on per-protein z-scores rather than a cosine against the
  published Axis 1 consensus vector: that vector's negative pole carries
  the immunoglobulins (IGHA2 −0.230, IGHG3 −0.151, IGHA1 −0.096)
  alongside the plasma proteins, so a whole-vector cosine reaches only
  −0.35 for the genuine plasma endmember and −0.66 for the
  immunoglobulin one — it cannot separate the two poles.
- **`atman align bootstrap --threads N`** — the bootstrap iterations and
  the pooled-subject jackknife replicates now run as a rayon parallel map
  (default: one thread per core; `--threads 1` restores serial execution).
  Every iteration already drew from its own `derive_sub_seed(seed, iter)`
  stream, and the accumulators are now folded in iteration order, so the
  output is byte-identical for any thread count (verified against the
  pre-change binary on ICA and NMF fixtures; pinned by a core
  thread-invariance test and CLI byte-identity tests for both
  decompositions). The first error is reported by lowest iteration index,
  as before. The flag is recorded in the sidecar. `rayon` becomes a
  workspace dependency, used only here.
- **`atman axes`** command group for subject-level score tables:
  `build` (representative archetype columns → axes, within-cohort
  orthogonalization, global z), `contrast` (Cohen d with CI, Welch, AUC,
  covariate-adjusted OLS with nested designs, BH within family, subject
  bootstrap, `--omnibus-factor` joint F-tests), `groups` (per-group means
  with t CI, medians, Kruskal–Wallis, reference-vs-group effects), `anchor`
  (Spearman with bootstrap CI and partial rows), `displacement`
  (case-minus-control vectors, pairwise cosine with bootstrap CI), `loco`
  (leave-one-cohort-out centroid stability). Covariate expressions (`z()`,
  `log10()`, `log2()`, `ln()`, arithmetic) in every design.
- **`atman concordance`**: Spearman/sign/Jaccard agreement of per-feature
  effect tables over shared features, with feature-resampling CIs and
  between-stage delta rows.
- **`atman score weighted`**: signed-weight signature transfer with an
  optional two-group summary (Cohen d, Welch p, AUC).
- **`atman de`**: `--include-controls`, `--condition-col`, `--subset`,
  `--collapse-others`, `--max-missing-fraction`, covariate expressions in
  `--design`, continuous `--contrast` terms; the OLS path reports pooled-SD
  Cohen d (`effect_size_method = cohen_d`) with CI; `de_results.tsv` gains
  `n_a`/`n_b`.
- **`atman residuals`**: `--covariates-tsv`, `--max-missing-fraction`,
  expressions in `--design`, `--output-variance`,
  `--output-variance-summary`, `--output-canonical-dir`.
- `atman-core`: `contrast` (Cohen d, AUC, Kruskal–Wallis, Spearman p, t CI,
  percentile, z-score) and `expr` (covariate expressions) modules.
- **`--collapse-genes none|mean|max-observed`** on `de`, `residuals`, and
  `score weighted`: protein groups sharing a gene symbol are reduced to one
  value per sample by a stated rule (lexically first assay, per-sample mean,
  or most-observed assay); affected-gene counts land in the sidecar.
- **`axes contrast`** skips bootstrap replicates whose resample leaves a
  categorical term with a single level (new `boot_n_skipped` column)
  instead of aborting.
- **`atman run`**: `vars:` with `${name}` substitution, `--dry-run`
  validation (input provenance across stages), a trailing `sidecar_hash`
  manifest column hashing each output's `.run.json`, and
  `--strict-outputs` (default on) failing a stage whose declared output is
  missing.
- **`atman scale absolute`**: per-sample total-protein rescaling
  (`a − log2(Σ 2^a) + log2(total_protein)`) for Reiber-style exponents.
- **`atman axes icc`**: one-way random-effects ICC(1) with unbalanced `k0`.
- **`atman axes tree`** and **`atman concordance --tree-*`**: average-linkage
  trees (cosine or Euclidean over displacement/centroid leaves; `1 −
  Spearman` over effect tables) with bootstrap clade support, scipy-style
  linkage tables, and Newick output.
- **`atman enrich ora --query-tsv`**: hand-made gene lists (one family per
  `query` value) against a `--universe`; output gains a trailing `query`
  column.
- **`atman de`**: `--require-cols` / `--require-numeric` sample filters.
- **`atman score weighted --summary-only`**.

### Changed

- `de_results.tsv` has two new trailing columns `n_a`, `n_b`; readers that
  index by header are unaffected.
- `plan_manifest.tsv` (`atman run`) has a new trailing column
  `sidecar_hash`, and a stage whose declared output is missing now records
  `exit_code 2` and aborts unless `--strict-outputs false`.
- `enrich ora` output has a new trailing column `query`.

- **`atman decompose nmf --max-missing-fraction`** — drop assays whose
  missing-sample fraction exceeds a configurable threshold (default 0.0,
  strict complete-case), mirroring the existing `decompose ica` assay
  filter. The sidecar records the resolved threshold plus
  `n_assays_retained` / `n_assays_dropped_missingness` so a run's effective
  protein universe is auditable without re-deriving it from the input.

- **`atman decompose nmf` / `atman align project` — `exp2-clip` and
  `shift-min` input transforms.** NMF requires non-negative input; these
  two pre-decomposition transforms let log2-ratio proteomics data feed NMF
  without a separate upstream step. `exp2-clip`: `2^clamp(x, -c, +c)`
  (`--transform-clamp`, default 6.0) restores a non-negative ratio scale
  while winsorizing extreme tails. `shift-min`: `x - min(X)` over the whole
  matrix, a sensitivity alternative. `align project --transform` gains the
  same two values for projecting new cohorts onto an NMF-trained atlas;
  `shift-min` always recomputes its shift on the cohort being projected
  rather than reusing the atlas's training-time shift. Both commands'
  sidecars record the resolved `transform_clamp` / `transform_shift` via a
  shared `TransformRecord`.

- **`atman align bootstrap --decomposition ica|nmf`** (default `ica`,
  byte-identical to prior releases). `nmf` runs single-seed
  multiplicative-updates NMF per resample — applied identically to the
  point estimate, every bootstrap resample, and every jackknife
  replicate — with dedicated `--beta-loss`, `--init`, `--nmf-max-iter`,
  `--nmf-tol`, `--transform`, and `--transform-clamp` flags, and the same
  SplitMix64 per-resample seed derivation the ICA path uses. Lets
  cross-cohort archetype stability (CI width, sign stability) be assessed
  under NMF exactly as it already is under ICA.

### Fixed

- **`harmonize fit --method quantile` no longer requires a complete
  subject.** It built its reference profile only from training subjects
  observed on every shared feature, which at real proteomic missingness
  is close to unsatisfiable — reported failing outright on the fold the
  comparator was most needed for. Each subject now contributes its own
  observed values interpolated onto the common quantile grid, so any
  subject with two or more observations counts. On a complete matrix the
  interpolation is the identity and the profile is unchanged, which is
  asserted rather than assumed.
- **`atman run --strict-outputs` checked that declared outputs exist, not
  that they contain anything.** A stage emitting a header-only table was
  recorded with `exit_code 0` and its empty output flowed downstream,
  where almost every aggregation reads an empty unit as a run of
  negative observations rather than as a missing input. Existence,
  non-emptiness and non-degeneracy are three different audits, and only
  the first was being run. Empty declared outputs are now detected
  (zero bytes, or a `.tsv` holding only its tab-separated header), warned
  about by name, counted in the new `n_empty_outputs` manifest column,
  and — under `--strict-outputs`, which is on by default — fail the
  stage. `--strict-outputs false` keeps the old permissive behaviour.
  The header rule is deliberately narrow: a one-line file without a tab
  is treated as a legitimate single-value output, since erring toward
  not flagging is the right direction for a check that can fail a
  pipeline.
  **Manifest schema change**: `plan_manifest.tsv` gains
  `n_empty_outputs` as its final column. Readers keying on column name
  are unaffected; anything parsing positionally should be checked.
  Prompted by the karna manuscript session, whose own coverage audit
  reported a layer at "100%" by counting files rather than rows, while
  279 of 886 units in it were empty — the check intended to catch this
  class reproducing the class one level up.
- **Bootstrap replicates could score an incomputable distance as a
  maximal one.** In `concordance --output-tree-linkage` and `axes tree`,
  a pair whose distance could not be computed fell back to
  `1 - 0.0 = 1.0` — maximally dissimilar. That turns an absence of
  evidence into a definite claim. The point estimate was already safe
  (`concordance` refuses table pairs sharing fewer than three features
  up front, and a zero-length displacement vector is rare), but the
  bootstrap was not: a feature resample can strip a pair that passes the
  full-data check below the floor, and a subject resample can collapse a
  leaf to a zero vector. Those replicates were silently pushing pairs
  apart inside clade support. Both distance functions now return an
  error instead of a fallback; the point estimate refuses with the pair
  named, and a replicate that cannot be computed is skipped, warned
  about, and excluded from the support denominator, which previously
  used the requested replicate count rather than the usable one.
  The shape was contributed by the karna manuscript session, which found
  it in an unrelated pipeline: an empty unit is silently recoded as a
  negative observation by almost every aggregation, and no parameter is
  involved.
- **MNAR-ICA seed stability was structurally unable to report
  stability.** With `--missingness-model abundance-conditional`, the
  alternative seeds behind the stability ranking ran plain FastICA while
  the reference ran the missingness-aware fit, so `n_stable_runs`
  compared two different methods at a Jaccard threshold. Reported from a
  real tree as exactly 0.000 for all 120 programs across three cohorts —
  a column named for a measurement that could not produce one, and which
  reads as "nothing was stable" rather than "this was never computed".
  Alternative seeds now run the reference's method, so `--n-seeds` means
  the same thing under both missingness models and the numbers are
  comparable across methods. The sidecar records
  `stability_alt_seed_method`.
- **Iterative solvers that ran out of iterations said nothing.** Found
  by auditing for the pattern behind the three defects above, rather
  than from a failing run. `decompose ica` computed `n_iterations` and
  `final_tol` on every fit and no caller ever read them, so a FastICA
  run that stopped at `--max-iter` produced loadings indistinguishable
  from a converged decomposition — in stderr, in the sidecar, and in the
  output file. It now warns naming the final tolerance and the requested
  one, and records `ica_converged`, `ica_n_iterations`, `ica_final_tol`
  and `ica_n_seeds_not_converged` (the last covering the alternative
  seeds behind the stability ranking, which were equally silent).
  `decompose unmix` had no convergence tracking at all in its FCLS
  abundance solver; it now counts per-subject solves that hit
  `--fcls-max-iter`, warns, and records `abundance_n_not_converged`.
  `decompose nmf` and the missingness-aware ICA path already reported
  convergence and are unchanged.
- **`modules discover --soft-power auto` fell back to β = 1, producing a
  degenerate result that looked well-formed.** When no power satisfies
  the scale-free criterion the sweep returned the β with the *largest*
  R², which sounds reasonable and is a trap: on a bulk proteome the fit
  degrades monotonically, so the largest R² sits at β = 1 — no
  soft-thresholding at all. Reported from a real run: on a 93-sample,
  9,376-feature contrast, R² fell from 0.505 at β = 1 to 0.133 at β = 20,
  the sweep chose β = 1, every feature landed in the grey catch-all, and
  `module-de` then emitted a single row testing the mean of all 9,376
  proteins between subtypes at q = 0.04. Exit code 0, well-formed output,
  meaningless number.
  Three changes. The fallback is now WGCNA's documented default by
  sample count for unsigned networks (9/8/7/6 as n crosses 20/30/40),
  clamped to the swept range, instead of the best-fitting β. The sidecar
  records `scale-free-fit-achieved`, `scale-free-best-r-squared`,
  `scale-free-best-r-squared-beta` and `soft-power-fallback-rule`, so β
  alone no longer has to carry the distinction between a selection and a
  fallback, and stderr warns when the criterion failed. And a discovery
  in which every feature is grey is now refused before anything is
  written, with `module-de` independently refusing a module set whose
  only member is grey — a single all-features "module" is not a module,
  and every statistic computed on it is meaningless while looking
  perfectly well-formed.
- **Automatic `k`-selection that lands on its own search bound now says
  so.** A rule returning `--k-max` has not selected anything: it ran out
  of room, and the criterion it was asked to satisfy may never have been
  met. Found in the wild — a six-cohort `decompose ica` run with
  `--k-selection cumulative-variance=0.80 --k-max 30` resolved to
  exactly 30 in five of six cohorts, which needed 28 to 56, so those
  five captured 70-79% variance while the configuration claimed 80%. The
  sidecar recorded `k_resolved: 30` beside `k-max: 30`, indistinguishable
  from a genuine selection landing on 30. `decompose ica`, `decompose
  nmf` and `decompose unmix --k auto` now warn on stderr naming the rule
  and the bound, and record `k_selection_bound_hit` (`k_max`, `k_min`, or
  null) in the sidecar. An explicit `--k` records null: it is not a
  selection.
- **`decompose unmix` (VCA) picked endmembers almost arbitrarily.** Two
  defects compounded. First, the vertex search ran in the full
  `p`-dimensional feature space instead of the `k`-dimensional signal
  subspace — the dimensionality-reduction stage of Nascimento &
  Bioucas-Dias (2005), which a code comment had dismissed as a
  constant-factor speed step and removed. In `p` dimensions a random
  unit direction is almost orthogonal to a `k`-dimensional simplex, so
  `argmax |u·y|` was decided by noise. Second, the Gram-Schmidt step made
  a single pass over a *non-orthogonal* basis, so the search direction
  was never actually orthogonal to the vertices already chosen. Together
  they made the selected vertex set a function of the RNG seed rather
  than of the data: on a 110-subject CPTAC GBM cohort with no imputation
  at all, seeds 42/43/44 produced three *disjoint* vertex sets.
  The subspace projection is now computed from the `n × n` Gram matrix
  (cheap for `p ≫ n`, same subspace), with eigenvector signs
  canonicalised for reproducibility, followed by the projective
  transform onto the simplex hyperplane; the mean direction seeds the
  orthogonal basis for the first draw only, as in the published
  formulation; and the basis is kept orthonormal. On a planted simplex
  the vertex set is now identical across seeds and recovers the planted
  samples exactly. **This changes `decompose unmix` output.** No result
  in either manuscript project depends on it (both confirmed they have
  never run it for a reported number), but a prior unmix run is not
  reproducible against this version and should be re-run. Real data
  remains harder than a planted simplex: on a 26-subject cohort eight
  seeds still gave six distinct vertex sets, so `docs/recipes.md` now
  tells users to sweep seeds and report selection frequency rather than
  a single winner.
- **`atman de --ridge-lambda auto` no longer silently means "no
  shrinkage".** `auto` was the flag's default and resolved to `0.0`
  without comment, so a user who selected it — or who simply accepted the
  default under `--test msqrob` — got an unpenalised fit while the CLI
  help and the run sidecar both said `auto`, which reads as data-driven
  selection. The default is now literally `0.0` (numerically identical,
  so no result changes), `auto` is refused under `--test msqrob` with a
  message naming the numeric alternative, and under any other test it is
  reported on stderr as ignored. The sidecar gains
  `ridge-lambda-resolved`: the penalty that actually applied, null when
  the shrinkage path never executed — so a sidecar mentioning a ridge
  value can no longer be misread as evidence that msqrob ran. Old
  sidecars are left as they are; they carry `ridge-lambda: auto`
  regardless of test and no resolved field, and the absence of the
  resolved field means "written before this change", not "no penalty".
  Verified against the pre-change binary: `de_results.tsv` is
  byte-identical for `--test ols` and `--test welch-t`, output hashes
  match, and the sidecar differs only by the added
  `ridge-lambda-resolved` key and the `ridge-lambda` default string. The
  in-repo CPTAC msqrob fixture was passing `auto` and silently getting no
  shrinkage; it now states `0.0` explicitly, with identical numerics.
- **`atman align bootstrap` — non-finite loadings now rejected loudly.**
  A degenerate resample (duplicate rows, zero-variance column) could in
  principle drive either decomposition to diverge; `NaN`/±∞ loadings
  previously flowed silently into cosine similarity, archetype grouping,
  and the bootstrap/BCa accumulators. A shared finiteness gate now rejects
  them immediately, naming the cohort and call site (point estimate /
  bootstrap iteration / jackknife replicate) in the error. The sidecar also
  now records the *resolved* `--transform-clamp` (6.0 when `--transform
  exp2-clip` is used without an explicit clamp) instead of the raw,
  possibly-null CLI value.
- **`atman network differential --mode edge-pairwise` no longer dies
  silently on large inputs.** The full per-edge × per-cohort-pair row
  vector (four owned labels per row, doubling on growth, then a second
  full copy as the output buffer) was allocated *before* `--top-rows`
  truncation; on the CPTAC pan-cancer workload (22.4M shared edges × 15
  cohort pairs ≈ 336M rows) both the uncapped and the `--top-rows 100000`
  runs exited with no output and no error on a 128 GB machine. The
  enumeration now streams through a bounded binary heap of compact
  index-based candidates (`atman_core::network_differential::
  edge_pairwise_top_k` / `edge_pairwise_for_each_top_k`), so a capped run
  holds `--top-rows` candidates regardless of input size, and rows are
  written straight to the tempfile instead of assembled in memory. Output
  is byte-identical to the previous sort-then-truncate path (verified
  against the 1.1.0+ binary on a three-cohort fixture at every cap; the
  canonical order is now `edge_pairwise_order` in core and pinned by a
  randomized equivalence test including exact `|z_diff|` ties). The
  uncapped path (`--top-rows 0`) counts candidates first and refuses
  above 100M rows with the count and a pointer to `--top-rows`, replacing
  the silent exit with a stated limit. An ignored integration test
  reproduces the CPTAC envelope (6 cohorts × 6,700 features).

## [1.1.0] — 2026-05-30

### Added

- **`atman detectability` — two-layer detection + abundance differential
  screen.** Separates "is this protein detected differently between
  conditions?" (a logistic regression on the present-vs-missing
  pattern, with log-odds-ratio effect size) from "is it abundant
  differently where detected?" (OLS on the detected-only subset). Treats
  below-LOD / missing as genuinely missing (MNAR), never imputed. The
  abundance layer requires `--min-detected-per-condition` (default 2)
  detected samples in BOTH conditions, otherwise it is suppressed with
  `skip_reason = one_sided_detection` (a 3-vs-0 split has no real
  between-group contrast). Both layers report effect sizes alongside
  BH-adjusted q-values.

- **`atman absence-topology` — protein-side co-absence clustering.** Groups
  proteins that go missing together across samples via union-find over the
  Jaccard distance of their missingness patterns, surfacing structured
  dropout (shared depletion, panel/plate effects) as continuous,
  auditable clusters rather than per-protein noise. Deterministic
  (cluster ids assigned by sorted membership). The protein-side dual of
  `recover-plex`.

- **`atman recover-plex` — sample-side plex / pair recovery.** Clusters
  samples by co-detection structure (union-find over Jaccard distance) to
  recover acquisition plexes or paired designs from the data when the
  metadata is missing or suspect. `--infer-pairs` proposes a pairing on a
  detection-coherent plex and **fails loudly** if inference cannot
  complete on a coherent plex (rather than silently emitting no pairs);
  detection-incoherent plexes are reported as a skip.

- **`atman robust-paired` — leave-one-subject-out sign-stability for paired
  DE.** Jackknifes each paired comparison, recomputing the effect with one
  pair dropped, and reports per-protein sign-stability (the fraction of
  LOSO replicates preserving the full-data effect direction) plus the
  most-influential pair. An exactly-zero replicate effect is treated as
  non-directional (not spuriously sign-stable). BH-FDR at the biological-
  replicate (pair) level.

- **`atman decompose nmf` — Brunet 2004 multiplicative-updates NMF.**
  Frobenius or KL divergence loss (selectable via `--loss`) with two key features:
  (1) **Multi-seed stability framework** (`--n-seeds`, `--seed-stability-threshold`)
  matching the ICA stability interface — reports top-N Jaccard seed stability for
  each recovered factor; (2) **k-selection** (`--k-selection`) with two methods:
  cophenetic correlation (samples permuted clustering patterns to measure
  robustness of the factor rank) and RSS elbow (mean reconstruction residual
  per k, picks smallest k where marginal improvement drops below threshold).
  Outputs `nmf_loadings.tsv`, `nmf_activations.tsv`, and sidecar JSON.
  Parity-validated against `sklearn.decomposition.NMF` (Frobenius and KL):
  max |Δ loadings| = 0.0 at machine epsilon scale (sklearn and atman reach
  identical convergence thresholds on the test fixture). Closes a long-standing
  gap in atman's factor-search toolkit: ICA finds independent axes, but
  semi-supervised / archetypal / non-negative applications need NMF.

- **`atman decompose ica --missingness-model abundance-conditional`
  — joint-iteration MNAR-aware ICA.** Fits a closed-form logistic detection-curve
  model on the observed-vs-missing pattern per sample, treating detection
  probability as sample-and-protein-specific. Runs weighted FastICA where each
  cell's contribution scales by its estimated detection probability, iterating
  the detection-curve fit and the ICA until convergence. MAR-collapse (no
  missingness) is bit-identical to plain FastICA (max diff = 0.0); MNAR recovery
  beats column-mean impute-then-decompose on synthetic ground truth (~1.2% MAE
  improvement at 2.2% missingness). Note: uses joint iteration with
  reconstruction-based imputation rather than detection-probability-weighted
  FastICA — empirical finding that probability-weighting biases recovery on
  heavy-tailed sources. Closes the MNAR gap for high-missingness proteomics
  (e.g. DIA-MS, targeted assays where detection is protein-×-sample-specific).

- **`atman de --adjust-for <covariates-tsv>`
  — admixture-adjusted differential abundance.** Composes external covariates
  from TSV into the design matrix (e.g., `nmf_activations.tsv` as a covariate
  for batch/confounding adjustment). Supports wide format (`sample_id` column
  followed by numeric covariate columns) or long format (auto-detected). Wired
  through all DE test paths: limma (parity to R limma ≤ 1e-6), msqrob2
  (log_fc parity 1e-14, t-stat divergence max 0.40 due to pooled-df eBayes
  vs msqrob2 per-protein Satterthwaite), OLS (1e-13 parity), mixed-model
  (golden-section REML; max Δt 9.3e-5 vs lme4), Welch-t (routes to plain OLS
  internally; HC3 robust SE would require new Cargo dep), and paired-t (1e-14
  parity). External covariate-TSV reader auto-detects wide vs long format.
  Sidecar records each `--adjust-for` path AND its input SHA-256.

- **`atman bench decompose --tools atman.<method>` prefix dispatch
  and `--fixture-nmf`.** Benchmark harness now recognizes native tool prefixes
  (`atman.ica`, `atman.nmf`, `atman.missingness-ica`) for explicit method
  routing. New `--fixture-nmf` flag directs NMF (which requires non-negative
  input) to the planted non-negative fixture at `bench/planted_archetypes_nmf_v1/`
  instead of the regular fixture. NMF achieves ≥0.99 archetype correlation on
  the appropriate fixture.

- **Alignment and program sidecars now record `decomposition_method`.**
  `atman align programs` and `atman align project` output sidecars now include
  a `decomposition_method` field (one of `ica`, `nmf`, `mixed`, `unknown`)
  documenting which decomposition method(s) were used to generate the inputs.
  Supports auditability chains: "which programs are these, and what algorithm
  found them?"

### Fixed

- **Numerically stable two-sided p-value tails across DE / limma / msqrob
  / ensemble / network_differential / the permutation null.** Replaced
  `2 * (1 - cdf(|t|))` with the survival function `sf()` to avoid
  rounding-to-zero on extreme tails (e.g. t = 39, df = 326). Affects
  p-values reported in `de_results.tsv`, `ensemble_de.tsv`, and
  `network_differential.tsv` (mode edge-pairwise) when effect sizes are
  large or residual degrees of freedom are high. The `atman null`
  permutation path previously used the unstable form in its per-iteration
  t-statistic while the observed path used `sf()`, putting the two arms of
  the empirical-FDR calibration on different numerical scales in the tail;
  both now share `atman_core::de::two_sided_t_p_value`.

- **`atman_core::de::fit_f_dist` homogeneous-variance case.** Now returns
  `(Inf, arithmetic_mean(s²))` when `excess <= 0` (residual variance equal
  across features), matching R limma's `fitFDist` fallback. Previously fell
  through to no shrinkage (df_prior = Inf implicitly, but s2_prior was not
  set to the pooled mean). Affects all atman DE tests using eBayes shrinkage
  when residual variance is near-homogeneous (e.g., high-replicate assays,
  processed matrices with post-QC low variance). `de_results.tsv` columns
  `df_prior` and `s2_prior` now populate correctly in this case.

- **`ingest-matrix` `abundance_raw` column now holds the pre-normalization
  value.** It previously received the same post-`--log2-transform`,
  post-normalization value as the `abundance` column, defeating any
  downstream raw-vs-effective comparison. `MeasurementRow` now carries the
  pre-normalization value separately and emits it in `abundance_raw`.
  `ingest-matrix` also now warns on matrix value-columns matching no
  sample and hard-fails on a protein/row length mismatch instead of
  silently truncating.

- **`MeasurementRecord::effective_abundance()` rejects non-finite values.**
  It returned `Some(NaN)`/`Some(inf)` for a non-dropped row whose abundance
  parsed to a non-finite literal, leaking into `module-de`, `bootstrap`,
  `ratio`, `residuals`, and `fold-change`. It now returns `None` for
  non-finite abundance, matching the peptide-record accessor.

- **`meta` sign-test trial count.** The binomial sign test used the total
  cohort count as `n` while excluding zero-effect cohorts from `k`,
  deflating the tail. It now uses `n = n_positive + n_negative`.

- **`de --test limma` F-test q-values.** The per-feature `f_p_value` column
  was emitted with `f_bh_q` always empty — the only p-column without
  multiplicity control. `f_bh_q` is now BH-adjusted per (comparison, panel)
  consistently with the moderated-t `bh_q`.

- **Confidence intervals use the correct critical value (msqrob, Tukey,
  Dunnett, Šidák post-hoc).** These emitted `effect ± 1.96·se` z-intervals
  while their p-values were t-based or family-wise-adjusted, so a CI could
  disagree with its own p-value at small df. Each CI now uses the matching
  critical value (Student-t at the fit df; studentized-range for Tukey;
  Dunnett critical value for Dunnett). p-values are unchanged; only
  `ci_low`/`ci_high` move.

- **Deterministic tie-breaking in `robustness` top-k selection.** The
  top-k changed-feature set was sorted by `|mean_diff|` from a `HashMap`
  with no secondary key, so features tied at the k-boundary were selected
  by hash iteration order (non-reproducible Jaccard/overlap output, and
  non-finite effects could hijack the set). It now uses a total-order
  comparator with an ascending feature-id tie-break and excludes
  non-finite effects; two further hash-ordered output paths in the same
  command were made deterministic.

- **`align project` keys activations by the matrix's own sample order.**
  It re-derived sample IDs from `samples.tsv` in file order, independent of
  the order `load_cohort_matrix` built the columns from; on a count match
  with a different order this attached activations to the wrong samples.
  The matrix's column order is now the single source of truth.

- **Provenance input-hash lists de-duplicated** (`validate`, `de`,
  `de --test ensemble`, `null`) — `measurements.tsv` was listed twice in
  the canonical-input set recorded in the run sidecar.

### Security

- **Updated `rustls-webpki` to 0.103.13** — resolves `RUSTSEC-2026-0104`, a
  reachable panic in certificate-revocation-list parsing reachable via the
  `enrich gprofiler` HTTPS path.
- **Hardened data-derived output paths.** `atman matrix` / `atman fold-change`
  built per-panel CSV filenames from the untrusted `panel` column; a crafted
  value (e.g. `../../x`) could write outside `--output-dir`. The data-derived
  filename component is now sanitized and rejected on a path separator or `..`.
- Added a `cargo audit` CI gate (`.cargo/audit.toml` scopes two
  transitive-and-unreachable advisories with justification), a CycloneDX SBOM
  (`docs/sbom-1.1.0.cdx.json`), a `SECURITY.md` threat model, and
  adversarial-input robustness tests + a `fuzz/` scaffold.

### Documented divergences (atman vs the named reference)

- **msqrob2 path with `--adjust-for`:** log_fc is bit-identical (Δ ≤ 2.89e-14);
  t-stat diverges (max Δ 0.40) due to atman's pooled-df eBayes shrinkage
  across all proteins vs msqrob2's per-protein Satterthwaite degree-of-freedom
  from lme4. Effect-size parity (log_fc) is the load-bearing claim; the
  t-statistic rank-order is preserved (no sign flips).

- **mixed-model REML with `--adjust-for`:** atman uses golden-section line search
  optimizer vs lme4's L-BFGS-B; max Δt = 9.3e-5 on the test fixture. REML
  profile log-likelihood is identical (to 1e-8) and the fixed-effect point
  estimates are bit-identical.

- **Welch-t with `--adjust-for`:** routes through plain OLS internally.
  HC3 robust standard errors (which are the "Welch heteroskedasticity-aware"
  property) would require a new Cargo dependency (currently forbidden).
  The Welch-like unequal-variance property is not preserved when
  `--adjust-for` is used. Consider this when interpreting
  `--test welch-t --adjust-for` output for high-heteroskedasticity data.

- **Missingness-aware ICA:** uses joint iteration with reconstruction-based
  imputation and cell-weight scaling rather than detection-probability-weighted
  FastICA — empirical finding that probability-weighting biases recovery on
  heavy-tailed sources.

- **`atman network differential` — cross-cohort differential
  coexpression.** Given N cohorts on a shared feature universe (input
  as `--inputs LABEL1=path1,LABEL2=path2,...`), computes one signed
  correlation matrix per cohort and summarizes how each protein-pair
  edge varies across cohorts in one of three runtime-selectable modes:
  - `--mode edge-pairwise`: per-edge × per-cohort-pair Fisher-z test
    of `H_0: corr_A = corr_B`. Output shape `n_edges × n_cohort_pairs`.
  - `--mode edge-summary`: per-edge cross-cohort summary (mean / sd /
    range / sign-flip count, plus conservation_score = min |r| and
    divergence_score = sd / (|mean| + ε)). One row per edge; suited to
    pan-cancer "conserved vs. cohort-specific" maps.
  - `--mode module`: per-module within-module connectivity per cohort
    plus a cross-cohort rewiring score (sd of connectivity). Requires
    `--gene-sets` in the same `set_name` / `gene_symbol` schema used
    by `enrich ora` / `enrich gsea` / `score signatures`.
  Signed Pearson (default) or Spearman via `--method`. `--min-overlap`
  drops edges with sparse subject support in any cohort. `--top-rows`
  caps emitted rows for the edge modes (sorted by `|z_diff|` and
  `divergence_score` respectively, with stable tie-breaks). Emits a
  `*.run.json` SHA-256 sidecar that records every cohort input.

- **CPTAC TMT proteome adapter
  (`adapters/cptac/cptac_tmt_proteome_to_atman.py`).** Converts CPTAC's
  protein-level TMT proteome wide TSVs (the `<TUMOR>_proteome.tsv`
  files distributed by the CPTAC Data Coordination Center) to the
  canonical samples / proteins / measurements TSV triple. Skips the
  leading `Mean` / `Median` / `StdDev` summary rows; discriminates
  `Log Ratio` (kept) from `Unshared Log Ratio` (dropped); filters
  non-biological columns (`TumorOnlyIR`, `NormalOnlyIR`, `QC*`,
  `Pool*`, `Reference`, `RefMix*`); promotes NCBIGeneID into
  `assay_id`. For cohorts whose sample IDs encode condition inline
  (e.g. CPTAC HCC `T<N>` / `P<N>`), pass `--sample-id-regex` plus
  `--condition-key-map 'T=tumor,P=paired_non_tumor'` to recover the
  condition labels at ingest. Validated end-to-end on a real 4-gene ×
  4-sample HCC subset by `crates/atman/tests/cptac_tmt_adapter.rs`,
  including a downstream `atman validate` for the recovered
  `tumor` / `paired_non_tumor` contrast.

- **`atman score signatures` — per-sample gene-set signature scoring
  (singscore).** Sibling to `atman score modules` with the same
  canonical input directory (`samples.tsv` + `measurements.tsv`) and
  the `gene_sets.tsv` schema (`set_name`, `gene_symbol`) shared with
  `enrich ora` / `enrich gsea`. v1 implements singscore (Foroutan et
  al. 2018): per-sample protein ranking with average-rank tie
  handling, normalized mean-rank score per signature, centered to
  `[-0.5, 0.5]` matching `singscore::simpleScore` with
  `centerScore = TRUE`. Deterministic — no randomness, no
  permutations. Missing values handled per-sample by exclusion.
  Validated against `singscore::simpleScore` on a planted fixture:
  TotalScore agreement within `1e-12` on all 12 (set × sample) rows.
  `--method` CLI flag reserved for future ssGSEA / GSVA
  implementations. Emits a `*.run.json` SHA-256 sidecar.

- **`atman enrich gsea` — pre-ranked gene-set enrichment analysis.**
  Sibling to `atman enrich ora` with the same `gene_sets.tsv`
  schema (`set_name`, `gene_symbol`) and `de_results.tsv` input.
  Implements the weighted Kolmogorov-Smirnov-style enrichment
  statistic of Subramanian et al. (2005) using the `fgseaSimple`
  position-mask permutation formulation (Korotkevich et al. 2019).
  Reports per-set `es`, `nes`, `p_value`, `bh_q`, and
  `leading_edge` genes. Ranking statistic configurable via
  `--rank-by` (default `signed_log10_p`; alternatives `t`,
  `log2_fc`); permutation null seeded via `--seed` for full
  determinism; set-size filtering via `--min-set-size` /
  `--max-set-size`. Validated against `fgsea::fgseaSimple` on a
  planted fixture: ES agreement within `1e-12` on all three
  fixtures (top-loaded, bottom-loaded, scattered). NES depends on
  the permutation null and is reproducible across atman re-runs at
  the same seed but is not byte-equal to fgsea (different PRNGs).
  Emits a `*.run.json` SHA-256 sidecar with all CLI parameters.

### Changed

- **Unified pseudo-random number generation across the workspace.** Six
  independently hand-rolled splitmix/LCG generators (in `null`, `bootstrap`,
  `ratio`, and core `ica` / `multivariate_t` / `decompose_unmix` /
  `align_bootstrap` / `gsea`) — which had already drifted apart and all used
  a modulo-biased bounded draw — are replaced by one reviewed
  `atman_core::rng` (`SplitMix64` / `Xoshiro256pp`) with an unbiased Lemire
  bounded draw. The raw xoshiro stream is preserved byte-for-byte, so
  decomposition point estimates (ICA / NMF / VCA) are unchanged. The
  unbiased bounded draw **does change the absolute outputs** of the
  permutation/bootstrap commands (`null`, `bootstrap`, `ratio`,
  `decompose null`, `decompose unmix --n-boot`, GSEA permutation,
  `align bootstrap`, Dunnett-Hsu Monte-Carlo) — same-seed reruns remain
  byte-identical, but values differ from any pre-1.0 internal runs of those
  commands, which should be regenerated.

- **Strict QC-flag parsing.** Reading `measurements.tsv` now errors on an
  unrecognized `qc_sample` / `qc_assay` value instead of silently treating
  it as `PASS`. The accepted set is `PASS` / `WARN` / `FAIL` / empty (empty
  = no flag = pass), matching `validate`. A garbled or lowercase QC cell is
  now a loud failure rather than a silent "everything passed".

- **Crash-safe (atomic) output writes.** `asymmetry`, `robustness`,
  `module-de`, and `module-trajectory` now write their outputs via the
  shared `atomic_write` (temp file + rename) like the rest of the toolset,
  so an interrupted run cannot leave a truncated TSV.

- **Honest build provenance.** The baked-in `ATMAN_GIT_SHA` now carries a
  `-dirty` suffix when the working tree has uncommitted changes at build
  time (`build.rs` is also worktree-aware, resolving the real git dir for
  linked worktrees). A clean checkout records a bare SHA; CI/release builds
  are unaffected.

- **Minimum supported / pinned Rust is 1.94.** Added a committed
  `rust-toolchain.toml` (channel 1.94) so local, Docker, and CI builds
  cannot drift; CI now also pins the supported toolchain and exports
  `BLAS_NUM_THREADS`/`OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS=1` so the
  byte-identity determinism tests run under their documented conditions.
  Dropped the unused `thiserror` workspace dependency.

- **Internal: large command modules split (behavior-preserving).**
  `commands/decompose.rs` (≈3,350 LOC) became a `decompose/` submodule tree
  (one file per subcommand) and `commands/de/mod.rs` (≈2,500 LOC) was
  reduced to dispatch plus `args` / `ols_design` / `output_rows` submodules;
  the duplicated union-find / Jaccard / median helpers shared by
  `absence-topology` and `recover-plex` were extracted to one
  `jaccard_cluster` module. Output is byte-identical (guarded by the
  determinism + golden suites). A new `determinism_rng_commands` test suite
  asserts cross-run byte-identity for every RNG-bearing command
  (`null`, `bootstrap`, `ratio`, `align bootstrap`, `de --test ensemble`,
  `decompose ica`, `robust-paired`).

- **`Platform` accepts custom identifiers without recompile.** The
  canonical `platform` field in `samples.tsv` / `proteins.tsv` /
  `measurements.tsv` previously required one of a closed list of
  well-known values (`olink_explore_ngs`, `maxquant_lfq`,
  `diann_report`, ...); any other identifier was rejected at validate-
  time. New `Platform::Custom(String)` variant accepts any non-empty
  identifier and round-trips byte-stably through the canonical schema.
  Adapters for new platforms (e.g. `cptac_tmt_proteome`,
  `bruker_timstof_diapasef`) can now declare their own platform string
  without modifying atman-core, restoring the adapter-as-extensibility-
  seam intent. Well-known variants are unchanged and still receive
  platform-specific code paths where they exist (Olink NPX QC, etc.).
  Empty platform strings remain a schema error.

- **Audit silent `0.0` fallbacks in similarity computations
  (`network.rs`, `modules_discover.rs`, `align.rs`).** Previously a
  pair whose metric was undefined (e.g. zero-variance Pearson input)
  or had fewer than 2 finite observations after the complete-case
  filter was silently treated as zero similarity, indistinguishable
  from a real uncorrelated pair. `pairwise_similarity` and
  `pairwise_abs_similarity` now return a `SimilarityAudit` alongside
  the matrix, recording `n_pairs_total`,
  `n_pairs_insufficient_overlap` (network only), and
  `n_pairs_undefined_metric`. The `network influence` and
  `modules discover` commands log the audit to their run sidecars
  (per-stratum where applicable) and warn on stderr when any pair
  fell back. `align::similarity` now returns `NaN` instead of `0.0`
  on undefined metrics; the downstream `pairwise_edges` filter
  already discards non-finite similarities, so behavior is unchanged
  while a failed compute is distinguishable from a true zero at the
  value level.
- **Expose CI alpha as a CLI flag on `atman align bootstrap`
  (`--ci-alpha`).** The percentile and BCa CIs over the bootstrap
  `n_cohorts` distribution were both hardcoded to a two-sided 95%
  band (α=0.05). Now configurable; default 0.05 preserves prior
  behavior. `BootstrapParams::ci_alpha` carries the value through
  `align_bootstrap` and into `to_row`, replacing the literal
  `0.025`/`0.975`/`0.05` constants. Validation rejects values
  outside `(0, 1)` at both the CLI and core entry points; the
  sidecar records the chosen `ci-alpha`. Doc comments on
  `BootstrapRow.ci_lower_n_cohorts` / `bca_lower_n_cohorts` no
  longer claim a fixed 95% band.
- **Expose program-stability flag fraction as a CLI flag on
  `atman decompose ica` (`--min-stable-seed-fraction`).** Previously a
  hardcoded constant of 0.9 inside an indirection
  (`threshold_fraction(_threshold) -> 0.9`) that ignored its only
  argument; the comment admitted it was kept fixed "for now". A
  program is now flagged unstable when the fraction of alternative
  seeds that recover it (best-Jaccard ≥ `--seed-stability-threshold`)
  falls below `--min-stable-seed-fraction` (default 0.9, preserving
  prior behavior). The sidecar records the value so re-runs are
  reproducible. Validation rejects values outside `[0, 1]`.
- **Generalize keratin contamination filter on `atman programs filter`
  to a configurable regex (`--contamination-pattern`,
  `--max-contamination-fraction`).** Replaces the hardcoded
  `KRT*`/`KERATIN` gene-name match. The default pattern
  (`(?i)^KRT|KERATIN`) preserves the previous behavior for
  plasma/serum/CSF cohorts where keratin is a skin-shedding
  contaminant. Override for tissue contexts where keratin is biology
  (`--contamination-pattern '^_NEVER_'` to disable), or to flag other
  contaminants — `'(?i)^HB[AB]'` for hemolysis, `'(?i)^IG[HKL]'` for
  immunoglobulin carryover, `'(?i)^MT-'` for mitochondrial. The
  regex is compiled up front (invalid pattern → loud error before
  any work). The renamed flag, schema columns
  (`contamination_fraction_top_n`, `contamination_pass`), and
  fail-reason (`contamination_signature`) are breaking changes; the
  sidecar now records both `max-contamination-fraction` and the
  effective `contamination-pattern`.

### Added

- **Below-LOD safety gate on `atman de`
  (`--max-below-lod-fraction`, `--allow-censored`).** Refuses to run
  when the fraction of `below_lod=1` rows among non-dropped
  measurements exceeds `--max-below-lod-fraction` (default 0.5).
  Atman's default DE tests treat missing as MCAR; left-censored data
  at higher rates biases `mean_diff` toward zero and distorts
  t-statistics. The gate reads only `dropped_by_qc` and `below_lod`
  from `qc_measurements.tsv` up front so the full record parse is
  deferred until the chosen test path loads it. Escape hatches:
  raise the threshold after confirming the test is appropriate, or
  pass `--allow-censored` (emits a stderr warning instead of
  aborting). Commit `58595a8`.
- **`atman ingest-matrix --normalize {none,median,quantile}` +
  unnormalized-loading warning.** Adds cross-sample normalization
  after any `--log2-transform`. Median centering subtracts each
  sample's median and adds the grand mean of sample medians back,
  preserving relative protein differences within samples. Quantile
  forces identical empirical distributions for complete-case assays
  and falls back to median centering for partially-observed assays
  so no data is silently dropped. Default is `none` (correct for
  Olink NPX, which arrives pre-normalized). When `--normalize none`
  is used and per-sample median abundance spans more than one log2
  unit, a stderr warning recommends `median` (suppressible via
  `--skip-normalization-check`). For SomaScan, MaxQuant/LFQ, DIA-NN,
  and Spectronaut, `--normalize median` is the recommended default.
  VSN, plate bridging, and ComBat-style correction remain out of
  scope — those belong upstream. Commit `7733a36`.
- **Geometric compartmental unmixing (`atman decompose unmix`).**
  Ports Vertex Component Analysis (Nascimento & Bioucas-Dias 2005) +
  Fully Constrained Least Squares (Heinz & Chang 2001) from
  hyperspectral remote sensing into atman's canonical TSV surface.
  Where `decompose ica` returns statistically-independent abstract
  axes, unmix returns *geometrically-identified pure endmembers*
  with per-subject fractional abundances that sum to 1 under
  non-negativity. VCA iteratively picks the most extreme sample
  projection onto the orthogonal complement of already-found
  endmembers (deterministic under `--seed` via SplitMix64); FCLS
  solves `argmin_α ||x − Eα||² s.t. α ≥ 0 ∧ Σα = 1` via projected
  gradient with Euclidean simplex projection (Duchi et al. 2008).
  `--abundance ucls` drops the simplex constraint for compositional
  transforms where `Σα = 1` has no meaning. Supports
  `--transform {none,log,clr,alr,ratio-anchor}`; ILR is refused
  because it changes feature ordering, and CLR/ALR combined with
  FCLS is refused without `--allow-unconstrained-simplex` because
  log-ratio coordinates don't admit a simplex interpretation.
  Outputs `endmembers.tsv` (rank-ordered per endmember),
  `abundances.tsv`, `unmix_diagnostics.tsv`, plus the standard
  sidecar. Initial v1 deferred `nfindr`, `--k auto` (HySime),
  `--n-boot` CI, and `--annotate-markers` to follow-ons; all four
  landed as separate entries below. Integration tests verify
  planted-3-endmember recovery
  at cosine ≥ 0.85, FCLS row-sum = 1 ± 1e-5 with non-negative
  entries, determinism under fixed seed, refusal on `k > n/2`, and
  refusal on CLR+FCLS without the escape hatch. Closes Priority 10.
- **`atman decompose unmix --annotate-markers` (ORA on endmember
  loadings).** Hypergeometric upper-tail test on the top-N proteins
  per endmember against user-supplied marker sets. Reads a marker-set
  TSV (`set_name\tgene_symbol`) and emits `unmix_enrichment.tsv` with
  columns `endmember_id`, `set_name`, `set_size`, `overlap`,
  `top_n`, `universe_size`, `p_value`, `bh_q`. Integration test
  confirms planted markers hit at `p < 1e-6` and decoy sets stay
  `p > 0.5`. Commit `aea7196` (residual 10d).
- **`atman decompose unmix --n-boot` (subject-level bootstrap CI).**
  Adds percentile CIs on both endmember loadings and per-subject
  abundances by resampling subjects with replacement, re-running
  VCA+FCLS on each resample, matching each bootstrap endmember to
  its best-`|cosine|` point-estimate counterpart, sign-correcting
  the match, then re-solving abundances on the original subjects
  under the matched+signed boot endmember matrix. Sub-seeds derive
  deterministically from the master `--seed` via SplitMix64
  (matching `align bootstrap`). Outputs `loading_ci.tsv` and
  `abundance_ci.tsv` alongside the point-estimate artifacts.
  Commit `5ed5bdb` (residual 10c).
- **`atman decompose unmix --k auto` (HySime-style elbow
  selection).** Sweeps `k` in `[--k-min, --k-max]`, computes the
  mean reconstruction residual norm per `k`, and picks the smallest
  `k` whose marginal improvement over `k − 1` drops below
  `--k-elbow-threshold × peak_improvement` (default 0.10). Emits
  `k_selection.tsv` recording every `k` in the sweep with
  `mean_residual_norm` and `marginal_improvement`. Falls back to
  `k_max` when every marginal improvement stays above the
  threshold. Commit `b2533ba` (residual 10b).
- **`atman decompose unmix --endmember nfindr` (N-FINDR
  extraction).** Iterative simplex-volume maximization (Winter 1999)
  initialized from VCA's picks. Swaps each endmember with the
  candidate sample that most increases the Gram-matrix determinant
  of the simplex; repeats until a full pass produces no swap or
  `--nfindr-max-passes` is exhausted. Supplements VCA when the PE
  simplex is sub-optimal under pure projection search. Commit
  `93a87f4` (residual 10a).
- **Unbalanced Dunnett–Hsu via deterministic MC on `atman de
  --post-hoc dunnett`.** When observed per-group sample sizes vary
  by more than `max/min > 1.25`, the test auto-switches from the
  equicorrelated multivariate-t (which assumes balanced design) to
  Dunnett–Hsu's per-pair correlation matrix
  `ρ_{ij} = √(n_i n_j / ((n_0 + n_i)(n_0 + n_j)))`, evaluated via
  deterministic Monte Carlo on a per-pair correlation structure.
  Added `atman-core::multivariate_t::{dunnett_hsu_correlation_matrix,
  pdunnett_hsu}`. Balanced designs still dispatch to the faster
  equicorrelated path. Commit `968a25c` (residual 9).
- **Jaccard + Spearman metrics in `atman align bootstrap`.**
  `--metric {cosine,jaccard,spearman}` extends the v1 cosine-only
  subject-level bootstrap. Jaccard uses top-N overlap on signed
  loading ranks; Spearman uses rank-correlation on the full
  loading vectors. Both flow through the existing matching,
  percentile CI, entropy, and BCa pipeline unchanged. Commit
  `cf271c5` (residual 3a).
- **Cross-tool benchmark harness (`atman bench decompose`).** New
  top-level `atman bench` command family. Scores `atman decompose
  ica` against named external tools on a shared planted-archetype
  fixture. Scoring math in `atman-core::bench_decompose`: recovered
  vs planted loading vectors matched via reciprocal-best absolute
  cosine, then scored via `archetype_correlation`
  (`|pearson|`), `recovery_jaccard` (top-N overlap), `runtime_seconds`,
  and `determinism_score` (1.0 iff two repeated runs yield byte-equal
  recovered loadings). Atman always runs natively via `fast_ica` +
  `canonicalize_ica`. Other tools are invoked via thin shell
  adapters at `bench/adapters/<tool>.sh` with a strict TSV-in /
  TSV-out contract (`recovered_loadings.tsv` in the same 3-column
  schema as `planted_loadings.tsv`); missing or failing adapters
  produce a single `tool_not_available = 1` row rather than
  aborting, so `--tools atman,fastica-icasso` still works on
  machines where fastica-icasso isn't installed. Ships the v1
  fixture at `bench/planted_archetypes_v1/` (2 archetypes × 20
  proteins × 20 subjects, deterministic cubed-Gaussian activations);
  atman recovers both planted archetypes at `|pearson| > 0.999`,
  Jaccard 1.0, determinism 1.0. Closes Priority 6.
- **Data-driven module discovery (`atman modules discover`).** New
  top-level command family. WGCNA-style pipeline in pure Rust:
  subject-level pairwise |pearson| / |spearman| correlation → soft-
  power adjacency `|r|^β` (β explicit or auto-selected to the
  smallest integer with scale-free topology `R² ≥ 0.8` and slope
  < 0) → topological overlap matrix (TOM) → UPGMA hierarchical
  clustering on `1 − TOM` → fixed-height tree cut with
  `--min-module-size` floor; features below threshold land in the
  `grey` catch-all. Also supports `--method hard-threshold` (binary
  adjacency + connected components). Outputs
  `modules_discovered.tsv` (schema matches the existing
  `modules.tsv` contract so it drops straight into
  `atman score modules`/`module-de`), `module_discovery_report.tsv`
  (per-module size, mean within-|r|, hub feature, module
  eigengene PC1 variance fraction), and
  `soft_power_diagnostics.tsv` (β sweep). Integration test recovers
  two planted correlation blocks on a 40-subject × 80-protein
  fixture: both 20-feature blocks concentrate into distinct
  non-grey modules of size ≥ 15; refuses cleanly when retained
  features fall below `--min-module-size`. Consensus resampling
  and WGCNA's full dynamic tree cut are documented follow-ons.
  Closes Priority 8.
- **Cohort projection onto a trained atlas (`atman align project`).**
  New subcommand alongside `align programs` / `align bootstrap`. Reads
  the archetype TSV emitted by a prior `align programs` plus the
  per-cohort loading TSVs it was built from, assembles an atlas loading
  matrix by averaging member programs per multi-member archetype on
  the intersection protein universe, loads a new cohort's canonical
  abundance matrix (remapping `assay_id`→label via the cohort's
  `proteins.tsv`), applies a compositional transform, and solves for
  per-subject activations via Cholesky on the normal equations under
  either `--projection ls` (plain least-squares) or `--projection ridge
  --ridge-lambda <λ>` (default). Outputs `projected_activations.tsv`
  (subject × archetype), `projection_qc.tsv` (per-subject
  `residual_norm`, `coverage_fraction`, `n_present`, `n_missing`),
  and a run sidecar recording the atlas inputs, transform, projection
  method, and the list of atlas proteins absent from the cohort.
  Integration tests verify planted-coefficient recovery to 1e-3 on
  full coverage, partial-coverage warning when the cohort drops an
  atlas protein, and clean handling of zero-intersection cohorts.
  Closes Priority 4.
- **Dunnett post-hoc (`atman de --post-hoc dunnett`).** Compares
  every non-reference level of the factor to the (alphabetically
  first) reference level under the joint multivariate-t distribution
  of the `m = levels − 1` correlated t-statistics. Uses the
  Dunnett–Curnow one-variate-plus-idiosyncratic representation to
  reduce the equicorrelated multivariate-t CDF to a 2D
  Gauss–Legendre quadrature (128-node outer × 64-node inner), which
  matches `emmeans(..., adjust = "dunnett")` at ≤ 1e-2 on canonical
  critical values for `m ∈ {2, 3, 5}` and `ν ∈ {10, 30, ∞}`. v1 uses
  the balanced-design correlation `ρ = 0.5`; unbalanced Dunnett–Hsu
  (heterogeneous correlations via Genz–Bretz) is a documented
  follow-on. Integration tests verify self-consistency of the
  dispatch against `atman_core::pdunnett` to floating-point precision
  and magnitude on the planted stage fixture (strong responder
  `adj_p < 0.05` on every Dunnett contrast; null protein stays
  `adj_p > 0.10`). Closes DEBT-6.
- **Tukey HSD post-hoc (`atman de --post-hoc tukey`).** Implements the
  studentized range distribution from scratch in
  `atman-core::studentized_range` via nested Gauss–Legendre
  quadrature (128-node inner × 48-node outer) with Legendre zeros
  computed from Bonnet's recursion — no external data tables. P-value
  adjustment uses `1 − ptukey(|estimate|·√2 / SE, nmeans=k,
  df=residual)`, which matches `emmeans(..., adjust = "tukey")` under
  arbitrary covariate adjustment and unbalanced `n_i`. CLI auto-emits
  all ordered pairs of the factor's observed levels unless
  `--contrast-list` restricts the family. Parity to R's `ptukey` is
  ≥ 3 decimals at canonical α=0.05 critical values for
  `k ∈ {3, 4, 5}`; integration tests verify self-consistency of the
  dispatch against `atman_core::ptukey` to floating-point precision
  and magnitude on a 36-sample planted fixture (planted stage effect
  recovers `adj_p < 0.05` on every pair; null protein stays
  `adj_p > 0.10`). Closes DEBT-5.
- **BCa CI + Shannon entropy for `atman align bootstrap`.** Retires the
  v1 deferral on the subject-level bootstrap. The summary now carries
  four new columns: `alignment_entropy` (bits of Shannon entropy over
  the bootstrap `n_cohorts` histogram — low = stable cohort coverage),
  `bca_lower_n_cohorts` / `bca_upper_n_cohorts` (bias-corrected
  accelerated 95% CI, replacing the percentile-only v1 band), and
  `bca_fallback_to_percentile` (1 when the BCa denominator goes
  non-monotone and the CI silently falls back to the percentile bound).
  Acceleration is estimated by pooled subject-level jackknife —
  dropping each subject across all cohorts, re-running the full
  ICA-per-cohort + alignment pipeline, and recording the
  matched-archetype cohort count. Closes DEBT-4.
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
  `de_results.tsv` (new `method` column). The grade is assigned from
  **per-method agreement, not a combined p-value**: for each
  (comparison, protein) it counts how many applied methods individually
  clear their own per-method BH-q < `--ensemble-q-threshold` (default
  0.05) — `n_significant` — and how many agree on the `mean_diff` sign.
  Grade is **VALIDATED** when both the significant fraction and the
  sign-consistent fraction reach `--ensemble-sign-fraction` (default
  1.00 — every applied method); **PROVISIONAL** when both reach
  `--ensemble-provisional-fraction` (default 0.50 — a majority);
  **INSUFFICIENT** otherwise (including when no method clears its own
  significance, regardless of how small the combined p is). A Stouffer
  `ensemble_p` / BH `ensemble_q` is still emitted but is a
  **non-calibrated ranking heuristic only — never an input to the
  grade**: the methods share one abundance matrix, so their p-values are
  positively correlated and the combined p is anti-conservative. Method
  applicability is auto-detected: methods whose required inputs are
  missing (e.g. msqrob without `--peptide-measurements`) are listed in
  `methods_skipped` rather than aborting the run. Thresholds are
  overridable (`--ensemble-q-threshold`, `--ensemble-provisional-fraction`,
  `--ensemble-sign-fraction`). On the bundled Dube heat-acclimation
  cohort the canonical heat-shock proteins (HSPA1A, HSPB1, DNAJB1) show
  consistent positive sign in PT2-PR2, but at n≈9–20 paired subjects no
  single method clears its own BH-q, so honest grading does **not** mark
  them VALIDATED — the prior VALIDATED label came only from the
  anti-conservative combined `ensemble_q`.
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

[1.1.0]: https://github.com/kevinj24fr/atman/releases/tag/v1.1.0
[1.0.0]: https://github.com/kevinj24fr/atman/releases/tag/v1.0.0
