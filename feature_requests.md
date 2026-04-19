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

### ~~5. `atman run --plan analysis_plan.yaml` — pre-registered analysis runner~~ *(shipped 2026-04-19)*

Reads YAML/JSON plan, dispatches each stage via `sh -c`, SHA-256 hashes
declared inputs and outputs per stage, and writes `plan_manifest.tsv`
with `plan_hash, stage_id, command, input_hash, output_hash, runtime_s,
exit_code, atman_version, system, started_at_unix_s`. Plan content drift
vs a previous manifest is detected and refused unless `plan_commit`
changes (or `--allow-drift` is passed).

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
| #9 pre-specification provenance | Feature 5 (shipped) | P1 (done) |
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

2. ~~**Ship P1 features 3 and 5 in a follow-up release.**~~ **Shipped
   2026-04-19.** `atman align programs` (feature 3) closes the
   annotation-alignment circularity, and `atman run --plan` (feature 5)
   lands pre-specification provenance as a hash-verifiable manifest.

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
