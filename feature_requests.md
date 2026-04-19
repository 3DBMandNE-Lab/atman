# Atman feature requests

Filed 2026-04-19. Source: cross-cohort CSF proteomics manuscript
(`/Users/kevinjoseph/Cursor/CSF_CrossDisease`, four cohorts, ICA +
annotation-constrained alignment + subject-level program coupling). Each
request maps to one or more reviewer-exposed weaknesses in that manuscript
and is scoped to a concrete atman CLI surface.

Priority ranking reflects how much the manuscript depends on the
capability. "Wait-worthy" = I would pause manuscript revision until the
feature lands. "Nice" = implementable in Python in a day; atman version is
preferred for reproducibility posture but not blocking.

---

## P0 — Wait-worthy

### 1. `atman decompose ica` with multi-seed stability

**Command sketch.**

```bash
atman decompose ica \
  --input-dir out \
  --k-selection cumulative-variance=0.80 \
  --k-min 3 --k-max 30 \
  --n-seeds 50 \
  --seed-stability-threshold 0.9 \
  --stability-metric jaccard-top20 \
  --output-loadings out/ica_loadings.tsv \
  --output-activations out/ica_activations.tsv \
  --output-stability out/ica_stability.tsv
```

**What it does.** Runs FastICA (or equivalently robust ICA) `n_seeds`
times per cohort with different initializations. For each program in the
reference-seed decomposition, computes recovery stability across the
other `n_seeds − 1` runs (top-20 loading Jaccard ≥ threshold against the
reference). Emits a `program, seed_stability_fraction, n_stable_runs`
table alongside the canonical loadings/activations. Programs below
`seed-stability-threshold` are flagged but not silently dropped (user can
gate downstream).

**Why we need it.** Single-seed FastICA is initialization-dependent; the
CSF manuscript currently uses `random_state=20260418` on one run and
reviewer raised this as the single most exposed technical item in the
paper. SIH at n=24/K=11 is the highest-risk case — if the decomposition
is unstable across seeds, the entire cross-cohort archetype alignment
needs re-derivation.

**Matches reviewer items.** #5 (single-seed ICA), partially #4 (forced
matched-K — same flag structure).

**Scope hints.** Rust FastICA implementations exist (e.g. via `ndarray`
and a symmetric decorrelation loop). Alternative: robust ICA variants
(FastICA with bootstrap-style restart clustering — ICASSO) provide the
stability metric natively. Either path is acceptable; the CLI contract
matters more than the internal algorithm.

---

## P1 — High-leverage, implementable in Python if atman doesn't ship soon

### 3. `atman align programs` with sensitivity sweep

**Command sketch.**

```bash
atman align programs \
  --loadings cohort_a/ica_loadings.tsv,cohort_b/ica_loadings.tsv,... \
  --annotations cohort_a/program_annotations.tsv,... \
  --metric jaccard \
  --top-n 40 \
  --tau 0.15 \
  --reciprocal-best \
  --category-constraint \
  --output out/archetypes.tsv

atman align programs sweep \
  --loadings ... \
  --metrics jaccard:top_n=[20,40,60]:tau=[0.10,0.15,0.20,0.25],cosine:tau=[0.20,0.30,0.40],spearman:tau=[0.20,0.30,0.40] \
  --compare-constrained-vs-unconstrained \
  --output-matrix out/alignment_sensitivity.tsv
```

**What it does.** Cross-cohort program-to-program matching by top-N
Jaccard / full-vector Spearman / cosine, with optional reciprocal-best
constraint and optional category-agreement constraint (from an
`annotations.tsv` that maps program → category). Sweep mode runs all
metric × threshold combinations and returns a matrix of
`(metric, top_n, threshold, n_multi_cohort_archetypes,
n_universal_archetypes, category_recovery_per_archetype)`.

**Why we need it.** The CSF manuscript currently has three bespoke Python
scripts for this (`ica_program_alignment.py`,
`ica_annotation_constrained_alignment.py`,
`ica_alignment_sensitivity.py`). Reviewer raised the circularity concern
— the annotation-constrained criterion is load-bearing — and demanded
the unconstrained vs constrained sweep be shown in main text. An atman
command makes the whole apparatus a single citation rather than three.

**Matches reviewer items.** #2 (annotation-constrained alignment
circularity).

**Scope hints.** Annotations input is cohort-agnostic: one TSV of
`program, category, top_annotation, top_annotation_p_value` per cohort.
Category-agreement constraint is boolean. Sweep is a nested loop over
metric/threshold combinations.

---

### 4. `atman enrich gprofiler` — live g:Profiler REST wrapper

**Command sketch.**

```bash
atman enrich gprofiler \
  --de-results out/de_results.tsv \
  --background out/measured_universe.tsv \
  --organism hsapiens \
  --sources GO:BP,GO:MF,GO:CC,KEGG,REAC,WP \
  --threshold-method fdr \
  --user-threshold 0.05 \
  --cache-dir .gprofiler_cache \
  --ontology-version pinned \
  --output out/gprofiler_enrichment.tsv
```

**What it does.** Wraps the g:Profiler REST POST endpoint currently
called from Python (`scripts/ica_program_enrichment.py`). Canonical TSV
output with `program, source, native, name, p_value, intersection_size,
query_size, term_size, effective_domain_size`. A `--cache-dir` + pinned
ontology version makes re-runs deterministic.

**Why we need it.** The CSF manuscript's M.9 currently calls g:Profiler
directly from Python with a `urllib.request` POST. This is cited as a
"live dependency" in M.14 with the acknowledgment that re-runs may vary
slightly. Atman wrapping the call with a cache + pinned ontology version
closes that caveat, and lets us cite a single atman command for program
annotation instead of ad-hoc Python.

**Matches reviewer items.** No direct review item, but strengthens
Methods reproducibility and removes a live-dependency footnote in M.14.

**Scope hints.** HTTP POST via `reqwest` or similar. JSON body matches
the g:Profiler REST spec. Rate-limit handling (0.15 s sleep between
queries is what the Python script does). Cache key = hash of (query
gene set, background, organism, sources, threshold, ontology version).

---

### 5. `atman run --plan analysis_plan.yaml` — pre-registered analysis runner

**Command sketch.**

```bash
atman run \
  --plan analysis_plan.yaml \
  --input-dir out \
  --output-dir out/provenance \
  --emit-manifest
```

with `analysis_plan.yaml`:

```yaml
name: csf_crossdisease_v3
plan_commit: <git-sha-of-this-file>
stages:
  - id: de_primary
    command: atman de --test ols --design "~ condition + age + sex"
  - id: bootstrap
    command: atman bootstrap protein --n 1000 --seed 20260418
  - id: null_permutation
    command: atman null --n 1000 --seed 20260418
  - id: ica_decompose
    command: atman decompose ica --n-seeds 50
  - id: program_coupling
    command: atman coupling --pairs pairs.tsv --method spearman
```

**What it does.** Executes each stage in order, captures `(stage_id,
command, input_hash, output_hash, runtime, exit_code, atman_version,
system_info)` into a manifest TSV. Re-running against the same plan
file checks that input hashes match and refuses to overwrite if the
plan has drifted without a new SHA.

**Why we need it.** Pre-specification provenance is currently a
git-log/trust claim in the CSF manuscript. Reviewer raised this as
essential for the "built-in positive control" framing of §5.3. An atman
plan manifest with hashes converts "trust my git log" into a citable,
hash-verifiable artifact.

**Matches reviewer items.** #9 (pre-specification provenance).

**Scope hints.** Plan YAML is declarative; atman shells out or
dispatches internally. Hash the plan file itself, hash each input
artifact, hash each output artifact. Manifest is one TSV, append-only
per run.

---

### 6. `atman residuals` — explicit covariate-adjusted residual matrix export

**Command sketch.**

```bash
atman residuals \
  --input-dir out \
  --design "~ age + sex" \
  --output out/residuals_long.tsv \
  --output-wide out/residuals_wide.tsv
```

**What it does.** Per-protein OLS fit with the requested nuisance
covariates, returns per-subject covariate-adjusted residuals. Long
format (`subject_id, gene_symbol, residual`) + wide format
(`gene_symbol` rows × `subject_id` columns, zero-filled at detection
threshold).

**Why we need it.** `scripts/build_module_deltas.py` does this in Python
and feeds it into ICA. Making residuals a first-class atman output
means `atman decompose ica` can consume it directly and the Methods
section cites one atman command instead of custom Python.

**Matches reviewer items.** No direct review item; cleaner Methods.

**Scope hints.** `atman de --test ols` already fits the model; exposing
the residuals is a small addition.

---

### 7. `atman ratio` — log-ratio testing between protein classes

**Command sketch.**

```bash
atman ratio \
  --input-dir out \
  --numerator modules.intrathecal_IgV \
  --denominator modules.plasma_IgG \
  --test wilcoxon \
  --groups MS-nonMS \
  --output out/ratio_results.tsv
```

**What it does.** Computes per-subject log-ratio of two module-scored
quantities, then tests disease-vs-control with Wilcoxon / t-test /
permutation. Effect size = median log-ratio difference + its
bootstrap CI.

**Why we need it.** The CSF manuscript's §5.3 plasmablast vs
plasma-derived humoral finding is currently framed as a Spearman
coupling, which is vulnerable to reviewer #7's "same Ig proteins load
both sides, sign reversal is tautological" critique. A log-ratio test
(intrathecal IgV ÷ plasma-derived IgG, subject-level) matches the
clinical MS literature's QIgG / IgG-index framing and is immune to
the sub-component bookkeeping issue. Also relevant for the A04
lysosomal vs A01 plasma cascade depression framing in §4.4.

**Matches reviewer items.** #7 (A02 bookkeeping) — ratio framing
sidesteps the subset-definition ambiguity.

**Scope hints.** Log-ratio computation is per-subject; tests are
standard non-parametric.

---

### 9. `atman bootstrap program` — signed-loading program bootstrap

**Command sketch.**

```bash
atman bootstrap program \
  --input-dir out \
  --loadings out/ica_loadings.tsv \
  --groups Case-Control \
  --n 1000 \
  --seed 20260418 \
  --output out/program_bootstrap.tsv
```

**What it does.** Like `bootstrap module`, but takes signed loadings
instead of binary membership. Per-subject program score is the
loading-weighted mean of per-subject protein abundances. Subject-level
bootstrap returns the sign-stability score that manuscript §3 already
uses.

**Why we need it.** Current `scripts/ica_program_de.py` does this by
averaging ICA activations across subjects and bootstrapping. Moving it
into atman with signed-loading weighting improves the statistic
(activation averaging loses loading magnitude) and replaces one more
Python script with a tested Rust run.

**Matches reviewer items.** None directly, cleaner Methods.

**Scope hints.** Generalization of `bootstrap module` with a loading
weight vector instead of an indicator vector.

---

## P2 — Nice to have

### 10. `atman programs filter` — interpretability filter

```bash
atman programs filter \
  --loadings out/ica_loadings.tsv \
  --annotations out/program_annotations.tsv \
  --min-annotation-pvalue 0.05 \
  --min-top-loading 0.5 \
  --max-keratin-fraction 0.2 \
  --output out/programs_interpretable.tsv
```

Flags programs that fail annotation significance, have diffuse loadings,
or show keratin-contamination. Lets Methods report "N of K programs
carry interpretable biology" (reviewer #12).

### 11. `atman within-cohort-rank` — explicit platform-heterogeneity preprocessing

```bash
atman within-cohort-rank --input-dir out --output-dir out/ranked
```

Replaces the DIA-vs-DDA intensity-scale heterogeneity across cohorts
with explicit within-cohort rank-transformation as an auditable pipeline
step (reviewer #11).

### 12. `atman de --per-subject-proxy QAlb` — barrier/clearance proxy regression

Auto-detects known physiological proxies in `samples.tsv` (QAlb, Evans
index, ventricular volume, `Q_{IgG}`) and includes them as continuous
regressors with a compartmental-physics summary section in the DE
output. Makes reviewer #8 a one-liner rather than a manual regression.

### 13. Archive deposit + `CITATION.cff` + Zenodo DOI

Not a CLI feature — but atman v1.0.0 with a Zenodo DOI and a
`CITATION.cff` at repo root lets Methods cite a single versioned
reference for the pipeline. Currently the CSF manuscript says
"implementation in other languages is straightforward"; a DOI-tagged
atman release makes that concrete.

---

## Feature-to-review-item matrix

| Review item | Addressed by | Priority |
|---|---|---|
| #2 annotation-constrained alignment circularity | Feature 3 | P1 |
| #4 forced matched-K sensitivity | Feature 1 (same flag structure) | P0 |
| #5 multi-seed FastICA stability | Feature 1 | P0 |
| #7 A02 bookkeeping | Feature 7 (ratio framing sidesteps) | P1 |
| #8 r×τ quantitative fit | Feature 12 | P2 |
| #9 pre-specification provenance | Feature 5 | P1 |
| #10 sign-test framing | Covered by completed module-level meta support | P1 |
| #11 platform heterogeneity | Feature 11 | P2 |
| #12 program interpretability count | Feature 10 | P2 |
| #13 category continuum rather than dichotomy | Covered by completed module-level meta support | P1 |
| Methods reproducibility posture | Features 4, 5, 6, 13 | P1 + P2 |

---

## Implementation sequencing that would maximally help the CSF manuscript

1. **Ship remaining P0 first (feature 1).** This converts the strongest
   reviewer liability into a one-line atman call. Without it, I have to
   implement multi-seed ICA stability in Python and hope reviewers accept
   the Python implementation.

2. **Ship P1 features 3 and 5 in a follow-up release.** These close the
   annotation-alignment circularity, the pre-specification provenance
   artifact. Both are implementable in Python as a fallback but reviewers respond
   better to "cited atman command" than "custom Python script in
   supplement."

3. **Features 4, 6, 7, 9 are methods-tightening.** They improve the
   Methods section's reproducibility posture but the paper can ship
   without them.

4. **P2 features are long-term polish.** If released in atman v1.1 or
   v1.2, they strengthen follow-up work.

---

## Notes for atman maintainers

- The CSF_CrossDisease manuscript is the first external consumer of atman
  v1.0.0 beyond the Dube NPX reproduction path. Any feature work here
  benefits other proteomics-analysis consumers of the canonical TSV
  contract.
- Contact: see `/Users/kevinjoseph/Cursor/CSF_CrossDisease/manuscript/`
  for the downstream analytic context each feature targets.
- Seed: all CSF manuscript analyses use `20260418` for direct
  compatibility with existing Python output; keeping that as the default
  seed in bootstrap/ICA/null commands simplifies cross-validation.
