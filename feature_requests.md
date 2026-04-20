# Atman feature requests

Open batch for the cross-cohort decomposition / alignment methods track.
The DE surface is mature (paired-t, welch, OLS-formula, mixed, limma,
DEqMS, msqrob). The next frontier for atman as a standalone methods tool
is rigor around the `decompose ica` + `align programs` pair: null
calibration, compositional handling, alignment uncertainty, projection,
and a benchmark harness that lets the tool be cited on its own merits.

All items below keep commandment 8 (Rust, SQLite, local, deterministic,
no cloud except PubMed API) and the canonical TSV contract.

---

## ~~Priority 1: Archetype null calibration (`atman decompose null`)~~ *(shipped 2026-04-20)*

Permutation null for FastICA archetype stability with three null modes
(sample-shuffle, protein-shuffle, gaussian-matched), BH-adjusted `null_q`,
and a `decision` column. Integration tests cover schema, determinism,
k-size refusal, and an end-to-end smoke test against the Dube cohort. See
`docs/analytical-roadmap.md` §8 for details.

---

## Priority 2: Compositional transforms in `atman decompose ica`

Bulk proteomics is closed-sum (total intensity is a plate/panel-level
artifact), so euclidean-geometry ICA on raw log2 abundance mixes real
biology with scaling. Downstream analysts currently hack this at the
`ingest-matrix` layer (per-sample albumin division, CLR, etc.). Atman
should expose the transform as a first-class `decompose ica` option so
the run sidecar records which geometry the archetypes live in.

Command:

```bash
atman decompose ica \
  --input-dir out \
  --k 20 \
  --transform clr \
  --seed 20260418 \
  --output-dir out/decomp_clr
```

Transforms:

- `none` (default; current behavior; sample-centered log2)
- `clr`: centered log-ratio. Per-sample: `x_i = log2(y_i) − mean_j
  log2(y_j)` over proteins with `y > 0`. Zero handling:
  `--zero-handling minprob` (replace zeros by per-protein 1st-percentile
  non-zero value before log) or `--zero-handling keep-na` (drop from the
  sample's log-ratio mean; preserves NA).
- `alr`: additive log-ratio against a reference protein.
  `--alr-reference ALB` or `--alr-reference-col reference_protein`. For
  CSF: albumin is the domain-standard reference.
- `ilr`: isometric log-ratio against an orthonormal Helmert basis.
  Numerically better-conditioned than CLR when K is large.
- `ratio-anchor`: per-sample divide by a chosen protein, then log2. No
  transform on the others. Equivalent to the CSF-manuscript "albumin
  normalization" pipeline.

All transforms must be pure, idempotent, and reversible (where
applicable; `ilr` has a closed-form inverse given the basis).

Outputs:

- Existing `loadings.tsv`, `activations.tsv`, `archetypes.tsv`, plus
- `transform_applied.json` capturing the transform name, parameters,
  reference choice, and zero-handling.

Acceptance:

- On a synthetic fixture with a planted compositional archetype (one
  gene up, all others down proportionally to conserve the sample sum),
  `--transform clr` recovers it cleanly while `--transform none` splits
  it across multiple components.
- CLR of a dataset with one protein set to its per-protein mean in
  every sample must yield zero loading on that protein for every
  archetype (invariance check).
- Run sidecar records transform identity. Two runs with the same
  `--transform` on the same inputs are byte-equal.
- Integration tests under `crates/atman/tests/decompose_ica_transforms.rs`
  covering the four non-default transforms + zero-handling branches +
  the compositional recovery fixture.

Why it matters for the methods paper: lets atman claim "geometry-aware
decomposition for closed-sum data." No standalone proteomics ICA tool
currently does this end-to-end. It also closes the loop on one of the
biggest CSF-manuscript findings — albumin normalization isn't a
preprocessing hack; it is a compositional transform, and the
decomposition result depends on which transform you chose.

---

## Priority 3: Bootstrap alignment uncertainty (`atman align bootstrap`)

`atman align programs` currently returns a point estimate at
`--cosine-tau 0.30` (or whatever the user picks). This is fragile: reviewers
will ask how confident we are that a given archetype is cross-cohort
universal. Add a subject-level bootstrap that resamples within cohort,
re-runs `decompose ica` + `align programs`, and reports a per-archetype
posterior probability of universality.

Command:

```bash
atman align bootstrap \
  --cohorts out_cohort1,out_cohort2,out_cohort3 \
  --k 20 \
  --cosine-tau 0.30 \
  --n-boot 200 \
  --seed 20260418 \
  --output-dir out_align_boot
```

For each bootstrap iteration:

1. Resample subjects within each cohort with replacement.
2. Re-run `decompose ica` per cohort.
3. Re-run `align programs` at the chosen `--cosine-tau`.
4. Record per-archetype membership (which cohorts recovered it).

Aggregate outputs (`align_bootstrap_summary.tsv`):

- `archetype_id` (consensus ID from the point-estimate run)
- `observed_n_cohorts`: cohort count at the point estimate
- `bootstrap_mean_n_cohorts`: mean over `--n-boot` iterations
- `bootstrap_prob_universal`: fraction of bootstraps where the
  archetype was recovered in all cohorts
- `bootstrap_prob_multi`: fraction recovered in ≥2 cohorts
- `ci_lower_n_cohorts`, `ci_upper_n_cohorts`: BCa or percentile interval
- `alignment_entropy`: Shannon entropy over recovered cohort-membership
  patterns (captures "this archetype sometimes aligns to X, sometimes
  to Y")

Acceptance:

- Deterministic under `--seed`; each bootstrap iteration draws its own
  deterministic sub-seed from a fixed derivation rule (e.g.
  `hash(seed, iter)` via blake3 or similar pure-Rust primitive).
- On a synthetic fixture with one "truly universal" archetype and one
  "cohort-specific" archetype, the universal one must have
  `bootstrap_prob_universal > 0.9`, the cohort-specific one `< 0.2`.
- Integration test under `crates/atman/tests/align_bootstrap.rs`
  covering the universal/cohort-specific recovery fixture, determinism
  under repeated seed, and refusal when cohorts have fewer than a
  minimum subject count (configurable, default 20).
- Run sidecar captures `n_boot`, `cosine_tau`, and per-iteration seed
  derivation rule for full reproducibility.

Why it matters for the methods paper: moves atman from "point-estimate
alignment" to "calibrated posterior probability." Directly addresses the
single most common cross-cohort concern: "you found this at τ=0.30, would
it survive a different random split of the same data?"

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

## Priority 5: Archetype variance decomposition (`atman decompose variance`)

Partition each archetype's subject-level activation variance into
cohort, condition, subject (if repeated measures), and residual. Answers
"how cohort-specific is this archetype?" quantitatively.

Command:

```bash
atman decompose variance \
  --activations out/activations.tsv \
  --samples out/samples.tsv \
  --factors "cohort + condition + (1|subject_id)" \
  --output out/archetype_variance.tsv
```

Model: one linear mixed model per archetype over the activation vector,
with user-specified fixed and random effects. Reuses the existing mixed
engine from `atman de --test mixed`.

Outputs (`archetype_variance.tsv`):

- `archetype_id`
- one column per fixed factor: `var_cohort`, `var_condition`, ...
- one column per random factor: `var_subject`, ...
- `var_residual`
- `icc_cohort`, `icc_condition`: intraclass correlations
- `f_test_<factor>_p`: Wald or LRT p-value per fixed factor

Acceptance:

- Reuses `atman de --test mixed` numerical engine; no new REML
  implementation.
- Refuses rank-deficient designs with the same error reporting as
  `atman de --test mixed`.
- Integration test covering a two-cohort fixture where one archetype is
  cohort-shared and another is cohort-confounded.

Why it matters: gives a principled quantitative answer to "what fraction
of this archetype is shared biology vs between-cohort batch?" — one of
the hardest-to-defend questions in any cross-cohort factor paper.

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

## Future directions (not yet spec'd)

The following are deliberately held back until Priorities 1-6 ship,
because each extends one of them and the priority-ordered design is
cleaner once the foundations are tested:

- **Missingness-aware ICA** (priority 7 candidate). Extend the Priority 2
  compositional transforms with a joint abundance + detection
  likelihood at decomposition time, same philosophy as `atman de --test
  msqrob` or proDA but at the component level. Lets archetypes be
  defined partly by "which proteins were detectable in this sample" —
  which is real biology for MNAR-heavy proteomics.
- **Hierarchical / nested alignment** (priority 8 candidate). Site →
  cohort → cross-cohort alignment for consortium data where a single
  cohort has multiple acquisition sites. Requires Priority 3 bootstrap
  uncertainty to be in place so each level of the hierarchy carries
  its own confidence interval.
- **Longitudinal tensor decomposition**. Samples × proteins ×
  timepoints. Relevant for MS PILOT, iNPH repeat measures, any
  consortium with longitudinal arms. Depends on Priority 5 variance
  decomposition since random-subject variance becomes the signal, not
  the nuisance.
- **Counterfactual archetype simulation**. `atman decompose
  counterfactual --set-archetype A0004 --to 0 --input activations.tsv`
  reconstructs a predicted abundance matrix with one archetype zeroed
  out. Interpretive aid; low-priority until Priorities 1-4 land.
- **Compositional effect size for DE**. A `--effect-size-scale clr` flag
  on `atman de` that reports CLR-space coefficients alongside the
  standard log-FC. Natural follow-on to Priority 2.
