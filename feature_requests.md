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

*(No open P0 items — Feature 1 `atman decompose ica` shipped 2026-04-19.)*

---

## P1 — High-leverage, implementable in Python if atman doesn't ship soon

### ~~3. `atman align programs` with sensitivity sweep~~ *(shipped 2026-04-19)*

Single-config and `--sweep` mode with Jaccard/cosine/Spearman similarity,
reciprocal-best and category-agreement filters, union-find archetype
assembly, and a sensitivity matrix over metric × top_n × tau ×
`--compare-constrained-vs-unconstrained`. Replaces
`ica_program_alignment.py`, `ica_annotation_constrained_alignment.py`,
`ica_alignment_sensitivity.py` with one atman command.

---

### ~~4. `atman enrich gprofiler` — live g:Profiler REST wrapper~~ *(shipped 2026-04-19)*

Live POST via `ureq` with on-disk cache keyed by a SHA-256 of the
canonicalized request (genes, background, organism, sources, threshold
method, user threshold, pinned ontology version). `--offline` mode
produces byte-identical runs from a pre-populated cache.

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

## P2 — Nice to have

## Feature-to-review-item matrix

| Review item | Addressed by | Priority |
|---|---|---|
| #2 annotation-constrained alignment circularity | Feature 3 (shipped) | P1 (done) |
| #4 forced matched-K sensitivity | Feature 1 (shipped) | P0 (done) |
| #5 multi-seed FastICA stability | Feature 1 (shipped) | P0 (done) |
| #7 A02 bookkeeping | Covered by completed ratio command | P1 |
| #8 r×τ quantitative fit | Covered by completed per-subject proxy regression | P2 |
| #9 pre-specification provenance | Feature 5 | P1 |
| #10 sign-test framing | Covered by completed module-level meta support | P1 |
| #11 platform heterogeneity | Covered by completed within-cohort rank transform | P2 |
| #12 program interpretability count | Covered by completed programs filter | P2 |
| #13 category continuum rather than dichotomy | Covered by completed module-level meta support | P1 |
| Methods reproducibility posture | Features 4, 5 plus completed residuals | P1 + P2 |

---

## Implementation sequencing that would maximally help the CSF manuscript

1. ~~**Ship remaining P0 first (feature 1).**~~ **Shipped 2026-04-19.**
   `atman decompose ica` lands the reviewer-exposed multi-seed FastICA
   stability requirement as a one-line CLI call.

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
