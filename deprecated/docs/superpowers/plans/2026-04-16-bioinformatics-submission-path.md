# Atman Bioinformatics Submission — Execution Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Drive Atman v1.0.0 from current state to a ScholarOne submission package for *Bioinformatics* (OUP), with *Bioinformatics Advances* as a declared fallback.

**Architecture:** Five sequenced stages with explicit dependencies. **Stage 0 (Zenodo release)** blocks every later stage because DOI back-fill is a text edit across manuscript/README/supp. **Stages 1 (benchmark), 2 (manuscript restructure), 3 (second-dataset validation)** run in parallel after Stage 0. **Stage 4 (supplementary + cover letter finalisation)** consumes outputs of 1/2/3. **Stage 5 (PDF + portal submission)** is final. Fallback to *Bioinformatics Advances* requires no reformat.

**Comparator set (revised):** OlinkAnalyze 5.0.0 (CRAN) + limma 3.66.0. **DEP dropped** (CRAN-archived, not maintained; redundant with limma for moderated-t). **Python comparators dropped** (Olink has no first-party Python DE package; `pyOpenMS`/`sanbomics` do not target NPX). Confirmed 2026-04-16 via CRAN + Olink R-universe; no new Olink-specific DE tool released since Jan 2026 (Olink's only package remains OlinkAnalyze; its 2026 version is 5.0.0).

**Tech stack:** Rust 1.75+ (atman + atman-core), Python 3.11 (adapters + figures), R 4.3+ (benchmark comparators only), Docker, GitHub Actions CI, Zenodo-GitHub integration, ScholarOne portal.

---

## Stage 0 — Zenodo release (blocks everything)

**Why this is the blocking stage:** The manuscript's Data/Code-availability block, `README.md`, and `docs/findings/cover_letter.md` all cite a Zenodo DOI that does not yet exist. Nothing ships without it. Rename blocks DOI because Zenodo mints the DOI against the current repo name.

**Files touched:**
- Verify: repo is named `atman` in GitHub settings
- Tag: `v1.0.0` (git)
- Create: Zenodo release record (web UI — user action)
- Modify after DOI issues: `README.md`, `docs/findings/2026-04-16-bioinformatics-manuscript.md`, `docs/findings/cover_letter.md`, `docs/findings/submission_checklist.md`, `CHANGELOG.md`

### Tasks

- [ ] **0.1 Pre-rename green-tests gate**

  Run: `cargo test --workspace --release`
  Expected: 56+ tests pass (summary notes 67 passing; verify). Abort stage if any test fails — don't rename over a red build.

- [ ] **0.2 Rename repo on GitHub (user action)**

  GitHub → Settings → Repository name → `atman` → Rename. GitHub auto-redirects the old URL.

  Post-rename local hygiene:
  ```bash
  git remote set-url origin git@github.com:kevinjoseph/atman.git
  git remote -v   # expect origin → atman.git
  ```

- [ ] **0.3 Grep for old name, update local references**

  Run Grep for `Atman` across repo; update any that are not intentional (e.g. archive paths, historical CHANGELOG entries should stay; live references in `README.md`, `Cargo.toml`, `docs/findings/2026-04-16-bioinformatics-manuscript.md`, `scripts/reproduce_paper.sh` must update).

  Expected residual: only archive and git-history references.

- [ ] **0.4 Bump version + CHANGELOG entry**

  Confirm `Cargo.toml:6` is `version = "1.0.0"`. Confirm `CHANGELOG.md` has `## [1.0.0] — 2026-04-XX` section with the benchmark, multiscale, stability-score, and platform-adapter line-items. If missing a final dated line, add it.

- [ ] **0.5 Commit release-prep state**

  ```bash
  git add -A
  git commit -m "release: Atman v1.0.0"
  ```

- [ ] **0.6 Tag v1.0.0 and push**

  ```bash
  git tag -a v1.0.0 -m "Atman v1.0.0 — paper release"
  git push origin main
  git push origin v1.0.0
  ```

- [ ] **0.7 Enable Zenodo-GitHub integration (user action)**

  https://zenodo.org/account/settings/github/ → flip toggle for `kevinjoseph/atman` to ON. **Must be toggled ON before the release in 0.8 — Zenodo only archives releases created after the toggle flips.**

- [ ] **0.8 Create GitHub release against v1.0.0 tag (user action)**

  GitHub → Releases → Draft a new release → tag `v1.0.0` → title "Atman v1.0.0 — paper release" → use CHANGELOG text for the body → Publish release. Zenodo ingests within ~1 minute.

- [ ] **0.9 Retrieve Zenodo DOI (user action)**

  https://zenodo.org/deposit (or the repo's Zenodo badge). Copy the version-specific DOI (ends in `.v1`) and the concept DOI. Record both below in the DOI back-fill task.

- [ ] **0.10 Back-fill DOI across all artefacts**

  Record:
  - Version DOI: `10.5281/zenodo.XXXXXXX`
  - Concept DOI: `10.5281/zenodo.YYYYYYY`

  Edit with Edit tool (don't rewrite files):
  - `README.md` — replace placeholder Zenodo badge + "Zenodo DOI pending" text
  - `docs/findings/2026-04-16-bioinformatics-manuscript.md` — replace both occurrences ("Zenodo DOI pending" in Abstract/Availability and in Data/Code Availability block)
  - `docs/findings/cover_letter.md` — replace "(DOI pending deposition, to be backfilled at revision)" with the concrete DOI
  - `docs/findings/submission_checklist.md` — mark the four Zenodo `- [ ]` items as `- [x]`

- [ ] **0.11 Commit DOI back-fill**

  ```bash
  git add README.md docs/findings/ CHANGELOG.md
  git commit -m "docs: back-fill Zenodo DOI (10.5281/zenodo.XXXXXXX) into manuscript + README"
  git push origin main
  ```

  **Note:** this commit is post-tag and deliberately NOT retagged — v1.0.0 archived on Zenodo does not contain its own DOI (chicken-and-egg). That is expected and standard.

**Stage 0 exit gate:** DOI present in README, manuscript, cover letter. Repo at `github.com/kevinjoseph/atman`. v1.0.0 tag pushed. CI green on main.

---

## Stage 1 — Benchmark completion (parallel with Stage 2 + 3)

**Why separate from manuscript restructure:** Stage 1 regenerates numerical outputs consumed by Stage 2's text. Both can be drafted in parallel, but Stage 2 swaps in final numbers only after Stage 1 lands its TSVs.

**Scope:** Complete Floor 3 (moderated-t as Atman default + null-FPR three-way comparison). Comparator cleanup — confirm DEP and Python comparators are absent from `benchmarks/` and manuscript text. The Floor 3 reviewer plan: re-run Dube benchmark with moderated-t as Atman's default, report containment both directions, add a permutation-null FPR comparison across Atman (paired-t, moderated, Welch-t) vs OlinkAnalyze `olink_ttest` vs limma `eBayes` at nominal q cutoffs.

**Files touched:**
- Create: `benchmarks/null_fpr/run_null_fpr.py` (orchestrator), `benchmarks/null_fpr/README.md`
- Create: `benchmarks/out/null_fpr_results.tsv` (one row per (tool, mode, n_perm, nominal_q) cell)
- Create: `docs/findings/supp/benchmark_null_fpr.md` (supplementary writeup)
- Modify: `benchmarks/atman/run.sh` (add moderated-mode invocation as the default block), `benchmarks/run_all.sh`
- Modify: `benchmarks/out/table1.tsv`, `benchmarks/out/table1_counts.tsv`, `benchmarks/out/table1_overlap.tsv` (regenerate with moderated-t as Atman default)
- Verify: no `benchmarks/dep/` or `benchmarks/python/` directories; no DEP/Python references in capability_matrix.tsv, table1.tsv, or manuscript.

### Tasks

- [ ] **1.1 Confirm comparator-set cleanup**

  Run Grep for `DEP|pyOpenMS|sanbomics|benchmarks/dep|benchmarks/python` across repo. Expected: zero hits in live paths (archive/git-history acceptable). If any live references remain in `benchmarks/`, `docs/`, or `README.md`, remove them.

  Capability matrix check:
  ```
  Read benchmarks/capability_matrix.tsv
  ```
  Expected columns: `capability`, `atman`, `OlinkAnalyze`, `limma`. No `DEP`, no Python columns.

- [ ] **1.2 Switch Atman benchmark default to moderated-t**

  Edit `benchmarks/atman/run.sh`:
  - Primary invocation: `atman de --test moderated --moderation-prior-df 4` (now the Atman default for the benchmark)
  - Secondary invocation retained as `--test paired-t` with output written to `benchmarks/out/atman/paired_t/` for the within-Atman mode-comparison figure

  No code change to `atman-core` is required — `--test moderated` already exists and is covered by the 75-test suite.

- [ ] **1.3 Re-run Dube benchmark with moderated-t default**

  ```bash
  bash benchmarks/run_all.sh
  ```

  Expected: `benchmarks/out/atman_de.tsv` now carries moderated-t q-values as primary; `benchmarks/out/table1.tsv` regenerated. Atman paired-t output preserved at `benchmarks/out/atman/paired_t/*.tsv` for supplementary comparison.

- [ ] **1.4 Refresh containment + overlap tables**

  `benchmarks/aggregate.py` should compute bidirectional containment: Atman ⊆ limma AND limma ⊆ Atman at q<0.05, q<0.10 per contrast. If the script doesn't currently print both directions, extend it — it's ~15 lines of pandas.

  Expected numerical behaviour: moderated-t Atman default should narrow the set-size gap vs limma and tighten bidirectional containment substantially (both tools are now empirical-Bayes-style moderated). OlinkAnalyze gap largely persists (Satterthwaite-paired vs moderated).

- [ ] **1.5 Implement null-FPR permutation driver**

  Create `benchmarks/null_fpr/run_null_fpr.py`. Reuse the sign-flip machinery from `scripts/calibrate_module_de.py` (sign-flip per subject with probability 0.5, preserves within-subject correlation). For each of B=1,000 permutations, write sign-flipped long-format TSV to a tempdir and invoke each tool:

  - **Atman (paired-t, moderated, welch-t)** — shell out to `atman de --test {mode} ...` over the permuted input TSV, read `de.tsv`, count rows with `bh_q < {0.01, 0.05, 0.10}`
  - **OlinkAnalyze** — shell out to a tiny R driver (mirror `benchmarks/olinkanalyze/run.R`) over the permuted NPX shape, count `BH-adjusted` rows below cutoff
  - **limma** — shell out to a tiny R driver (mirror `benchmarks/limma/run.R`), count `adj.P.Val` rows below cutoff

  Output (`benchmarks/out/null_fpr_results.tsv`): one row per `(tool, mode, contrast, perm_id, q_cutoff, n_rejected)`. Aggregation: mean rejection count and empirical FPR per (tool, mode, q_cutoff) across permutations and contrasts.

  Reuse opportunity: the permutation loop for OlinkAnalyze/limma is expensive (~10s each per perm × 1000 perms × 2 tools = ~6 hours). Parallelise with `joblib` or pre-stage all B=1000 permuted TSVs once and sweep tools × perms. Budget: scale B down to 500 permutations per tool if wall-clock exceeds 8 hours; document the power floor (empirical FPR resolution 1/501 ≈ 0.002).

- [ ] **1.6 Run null-FPR benchmark**

  ```bash
  python3 benchmarks/null_fpr/run_null_fpr.py \
    --measurements out/qc_measurements.tsv \
    --samples out/samples.tsv \
    --contrasts PT1-PR1 PT2-PR2 PT2-PT1 PR2-PR1 \
    --n-perm 1000 \
    --output benchmarks/out/null_fpr_results.tsv \
    --summary benchmarks/out/null_fpr_summary.tsv
  ```

  Expected result (hypothesis to either confirm or honestly report): paired-t at n=10 has somewhat inflated type-I at q=0.05 (expected ~0.05 × 2938 = 147 rejections per contrast; if observed > 200 → inflation); moderated-t and limma are close to or below nominal; OlinkAnalyze close to paired-t.

  **Honest-reporting commitment:** if paired-t is inflated, document it plainly and position paired-t as "available for n ≥ 20; moderated-t is Atman's small-n default". This is already the direction Stage 1.2 moves in.

- [ ] **1.7 Write Floor 3 supplementary**

  Create `docs/findings/supp/benchmark_null_fpr.md`. Tables: (a) empirical FPR per (tool, mode, q_cutoff) across 1,000 perms; (b) mean false-positive count per contrast per tool; (c) paired-t vs moderated-t within-Atman FPR. Reuse the prose template from `docs/findings/supp/null_calibration.md`.

- [ ] **1.8 Update Table 1 capability row + benchmark prose inputs**

  If paired-t shows inflation, add a row to `benchmarks/capability_matrix.tsv` (or a note to Table 1) explicitly stating default mode per tool: Atman = moderated, OlinkAnalyze = paired-t (Satterthwaite), limma = moderated. Drives Stage 2.5 prose.

- [ ] **1.9 Commit Stage 1 output**

  ```bash
  git add benchmarks/ docs/findings/supp/benchmark_null_fpr.md
  git commit -m "benchmark: Floor 3 — moderated-t default + null FPR comparison"
  ```

**Stage 1 exit gate:** `benchmarks/out/null_fpr_summary.tsv` exists. `benchmarks/out/table1.tsv` regenerated with moderated-t Atman default. `docs/findings/supp/benchmark_null_fpr.md` drafted. Capability matrix and `run_all.sh` reflect comparator set (OlinkAnalyze + limma only).

---

## Stage 2 — Manuscript restructure (parallel with Stage 1 + 3)

**Why parallel with Stage 1:** prose edits for Floors 1 (null calibration) and 2 (S-score simulation) are independent of Stage 1's numerical refresh. Swap in Stage 1's final numbers in 2.5 once they land.

**Scope:**
- Integrate Floor 1 null-calibration reference into Methods + Results (already partially in manuscript; confirm complete).
- Integrate Floor 2 S-score simulation honest reframe into Results + Discussion (already done; final sweep).
- Consume Floor 3 outputs (Stage 1): update benchmark prose to reflect moderated-t-default + null-FPR table.
- Demolition: cut Wei Alzheimer case-study full section to a short "Methods capability note" (preserves Welch-t + scipy concordance line but drops the Results prose and figure).
- Framing: retarget to *Bioinformatics* Original Paper (currently the target); keep *Bioinformatics Advances* as prepared fallback. Title + abstract should not change across venues.
- Word-count trim to ≤5,000 countable body (checklist notes 4,168 as of last count; Stage 1 + Floor integration likely push this up — needs audit).

**Files touched:**
- Modify: `docs/findings/2026-04-16-bioinformatics-manuscript.md`
- Modify: `docs/findings/submission_checklist.md` (tick completed items, add Floor-3 rows)
- Delete or demote: any standalone Wei section content (Fig S-wei if present; existing manuscript has Wei as one of three case studies in the abstract and methods)

### Tasks

- [ ] **2.1 Manuscript audit: counted-body word count**

  Run: `python3 scripts/count_manuscript_words.py` if it exists, otherwise use pandoc or `wc -w` on a pre-stripped version. Target: ≤5,000 words excluding abstract + references + figure/table captions. Record baseline.

  If the script does not exist, create `scripts/count_manuscript_words.py` in under 30 lines (strip YAML/markdown; split abstract/body/refs by section header; print each count).

- [ ] **2.2 Verify Floor 1 (null calibration) integration**

  Search manuscript for the Floor 1 sentence pattern: "sign-flip permutation null that re-learns modules on permuted data for B = 1,000 permutations per contrast". Expected: present in Methods §Multiscale and Results §Module-DE. If absent, add a single sentence in each (one-line insertion, no structural change).

  Confirm supplementary pointer `docs/findings/supp/null_calibration.md` is cited exactly once in Methods.

- [ ] **2.3 Verify Floor 2 (S-score simulation) integration**

  Search manuscript for the honest S-score phrasing: "$S$ is not proposed as a discovery-optimal ranking" and "stability-prioritised alternative". Expected: present in Results §Stability and Discussion. Confirm supplementary pointer `docs/findings/supp/s_score_simulation.md` is cited exactly once in Methods.

- [ ] **2.4 Demolition: collapse Wei Alzheimer case study**

  Currently the manuscript references Wei 2025 in:
  - Abstract ("Wei et al. 2025 Alzheimer plasma (Olink Target 96, n=10/group, unpaired)")
  - Introduction (one mention, reference [18])
  - Results (its own cross-platform generalizability subsection — the one to demolish)
  - Methods (`--test welch-t` validation)
  - Figures (Fig 4C likely references Wei if present)

  Replacement policy: **keep** the Abstract mention (needed for cross-platform claim), **keep** reference [18], **keep** Methods sentence on Welch-t + scipy concordance. **Cut** the Wei Results subsection to one sentence: "The pipeline's Welch two-sample mode reproduces `scipy.stats.ttest_ind` to IEEE-754 machine precision; applied to an unpaired Alzheimer plasma cohort (Wei et al. 2025, Olink Target 96, n=10/group) [18] the Atman Welch-t output matches the scipy reference exactly. A fuller case study in a disease cohort with this design remains future work." **Cut** Fig 4C panel if Wei-specific; renumber subsequent panels.

  Estimated word reduction: 350–500 words. Feed into 2.1 trim budget.

- [ ] **2.5 Swap in Stage 1 numbers**

  After Stage 1 lands `null_fpr_summary.tsv`, update:
  - §Benchmark — mention "under moderated-t as Atman's default at the cohort size tested (n=10)"
  - §Benchmark — add one sentence with the null-FPR headline (e.g. "Under 1,000 sign-flip permutations of condition labels, empirical false-positive rates at nominal q<0.05 were X (Atman moderated), Y (limma), Z (OlinkAnalyze paired-t); Atman paired-t was W (reported for completeness; moderated-t is the small-n default)")
  - §Methods — one sentence pointing to `docs/findings/supp/benchmark_null_fpr.md` for full tables
  - Figure 2 caption — if UpSet or containment figure needs refresh, note the moderated-t change

- [ ] **2.6 Framing pass: tighten abstract + intro + discussion for Bioinformatics**

  *Bioinformatics* Original Paper reader expects: novel method, validation against established tools, demonstrated utility. Current framing is already aligned. Light-touch pass:
  - Abstract's Motivation should name the gap in one sentence (pooled significance loses subject-level structure; small-cohort q-ranking is unstable) rather than two.
  - Introduction's named-gap paragraph should end with "We address these two gaps with an integrated deterministic implementation" rather than restating three contributions.
  - Discussion: lead with the single most defensible claim (byte-exact reproduction + containment + honest S-score characterisation + null-calibrated module-DE), not a contribution laundry list.

  This is last-mile framing — do not rewrite. Edit in place.

  *Bioinformatics Advances* fallback: no textual changes required. Both venues accept the same format; portal re-upload only.

- [ ] **2.7 Final word-count audit**

  Re-run 2.1. Target: ≤5,000 countable body words, ≥200-word headroom. If over cap after 2.4 + 2.5 inserts, prioritise trimming the exploratory Dube Fig S3 prose (already labelled "hypothesis-generating only") and the detailed LOO numbers which are duplicated in supplementary.

- [ ] **2.8 Citation integrity sweep**

  Grep references [1]–[19] each appear in body at least once. Confirm references list has no stray entries after Wei demolition (ref [18] should stay since Wei is still cited once). If any reference is orphaned after demolition, remove it and renumber (but prefer not to renumber — keep [18] live via the Methods-note Wei citation).

- [ ] **2.9 Commit Stage 2 output**

  ```bash
  git add docs/findings/2026-04-16-bioinformatics-manuscript.md docs/findings/submission_checklist.md
  git commit -m "manuscript: Floor integration, Wei demolition to methods note, framing pass"
  ```

**Stage 2 exit gate:** word count ≤5,000 countable body. Floors 1/2/3 integrated. Wei demoted to one-sentence methods note. Submission checklist shows green on all manuscript-body rows except user-fill stubs.

---

## Stage 3 — Second-dataset validation (parallel with Stage 1)

**Why separate:** Gisby 2021 (Olink Target 96) and Gisby 2022 (SomaScan v4.1) are already ingested and analysed. This stage is a reproducibility re-run post-Stage-0-rename plus numerical sanity check. It's **not** a new-dataset ingestion task; that work landed earlier in the session.

**Files touched:**
- Re-run only: `scripts/reproduce_gisby.sh`, `scripts/reproduce_gisby_somascan.sh`
- Verify: `docs/findings/supp/gisby_generalizability.md`, `docs/findings/supp/gisby_somascan_generalizability.md`

### Tasks

- [ ] **3.1 Re-run Gisby 2021 Olink pipeline**

  ```bash
  bash scripts/reproduce_gisby.sh
  ```

  Expected: `out_gisby/de/late-early.tsv` matches previously-recorded hit count (120 at q<0.05). `out_gisby/robustness/loo_sign_stability.tsv` matches LOO sign match 0.972, top-20 Jaccard 0.911. If any numerical drift, treat as a regression: bisect from v1.0.0 tag back to the commit that produced the original numbers.

- [ ] **3.2 Re-run Gisby 2022 SomaScan pipeline**

  ```bash
  bash scripts/reproduce_gisby_somascan.sh
  ```

  Expected: `module_14` interferon/chemokine module at q ≈ 3.5 × 10⁻⁴; CXCL10, TNFSF13B, GRN with q<0.01. Numerical drift protocol identical to 3.1.

- [ ] **3.3 Cross-platform claim sanity check**

  Verify abstract sentence "Core inflammation biology (CXCL10, TNFSF13B, GRN) emerges with q<0.01 on both Olink and SomaLogic platforms on the same patient cohort." Confirm the three proteins meet q<0.01 on Gisby-2021-Olink AND Gisby-2022-SomaScan outputs.

- [ ] **3.4 Commit validation log (if any numerical updates)**

  Only if a number drifted: update the manuscript sentence, update the supplementary doc, commit with `docs: refresh Gisby cohort numerical claims post-v1.0.0 rerun`.

**Stage 3 exit gate:** Both Gisby runs reproduce original numbers bit-for-bit or manuscript is updated to the new numbers. Cross-platform claim verified.

---

## Stage 4 — Supplementary + cover letter finalisation (blocked by 1, 2, 3)

**Scope:** back-fill final numbers into supplementary docs, finalise cover letter, tick submission checklist, request user-fill stubs.

**Files touched:**
- Modify: `docs/findings/supp/benchmark_full.md` (moderated-t default + null-FPR section)
- Modify: `docs/findings/cover_letter.md`
- Modify: `docs/findings/submission_checklist.md` (tick)
- No new files.

### Tasks

- [ ] **4.1 Update benchmark full supplementary**

  Ensure `docs/findings/supp/benchmark_full.md` reflects:
  - Atman's benchmark default is now moderated-t
  - Paired-t remains available, reported as secondary
  - Reference to `docs/findings/supp/benchmark_null_fpr.md` for the null-FPR sweep

- [ ] **4.2 Final cover letter pass**

  Edit `docs/findings/cover_letter.md`:
  - Replace "(DOI pending deposition, to be backfilled at revision)" with the concrete Zenodo DOI from 0.10
  - Confirm the three-contribution paragraph still aligns with the post-demolition manuscript (Wei is no longer one of the three case studies; Gisby 2021 + Gisby 2022 are)
  - Concordance sentence should now cite moderated-t Atman vs OlinkAnalyze + limma

- [ ] **4.3 User-action stubs (explicit request block)**

  Generate a single markdown block for the user listing exactly what they must complete (no inference on my part):
  - Corresponding-author ORCID on ScholarOne portal
  - Affiliation statement for title page
  - Funding statement
  - Competing-interests declaration
  - CRediT author-contributions taxonomy
  - 2–3 suggested reviewer names (Olink/proteomics tools community, no direct collaborators)
  - Optional opposed reviewer names

  Save the block into `docs/findings/cover_letter.md` under the existing "Suggested reviewers" / "Opposed reviewers" sections; don't duplicate elsewhere.

- [ ] **4.4 Submission checklist sweep**

  Walk `docs/findings/submission_checklist.md` top to bottom. Tick every row whose evidence file now exists. Leave user-action rows untouched.

- [ ] **4.5 Commit Stage 4 output**

  ```bash
  git add docs/findings/cover_letter.md docs/findings/submission_checklist.md docs/findings/supp/
  git commit -m "submission: finalise cover letter + supplementary + checklist"
  git push origin main
  ```

**Stage 4 exit gate:** Cover letter references concrete DOI. Submission checklist: every non-user-action row green.

---

## Stage 5 — PDF build + portal submission (after Stage 4)

**Scope:** user-action heavy. I render the PDF locally; user uploads.

**Files touched:**
- Create: `docs/findings/2026-04-16-bioinformatics-manuscript.pdf` (generated artefact — do NOT commit)
- Verify: all figure PDFs in `docs/findings/figures/` and `docs/findings/heterogeneity/figure_pack/`

### Tasks

- [ ] **5.1 Render manuscript PDF**

  pandoc with a Bioinformatics-leaning style: double-spaced, 12 pt, numbered lines, standard margins. Example (adjust if a template exists):
  ```bash
  pandoc docs/findings/2026-04-16-bioinformatics-manuscript.md \
    -o docs/findings/2026-04-16-bioinformatics-manuscript.pdf \
    --pdf-engine=xelatex \
    -V mainfont="Times New Roman" -V fontsize=12pt -V linestretch=2 \
    -V geometry:margin=1in -V linkcolor=blue \
    -V header-includes="\usepackage{lineno}\linenumbers"
  ```

- [ ] **5.2 Figure 300 dpi sweep**

  Confirm every Fig 1–6 panel is available as vector PDF or 300 dpi PNG. Checklist item already notes Fig 1 complete; audit Figs 2–6 against `docs/findings/heterogeneity/figure_pack/`. If any are below 300 dpi PNG, regenerate from the per-figure Python script with `dpi=300` explicit.

- [ ] **5.3 Pre-submission self-review**

  Read the rendered PDF end-to-end once. Confirm: no "Zenodo DOI pending" string, no "USER-FILL" stub strings, no `Atman`, no stray Fig 7/Fig 8 references, all in-text citations resolve.

- [ ] **5.4 Upload to ScholarOne (user action)**

  https://mc.manuscriptcentral.com/bioinformatics → Original Paper submission. Upload:
  - Manuscript PDF
  - Cover letter (from `docs/findings/cover_letter.md`, converted to PDF or pasted)
  - Supplementary materials bundle (zip of `docs/findings/supp/`)
  - Figure files individually if portal requires
  - Suggested/opposed reviewers in portal fields

- [ ] **5.5 Confirmation + archive**

  Capture the ScholarOne submission ID. Update `docs/findings/submission_checklist.md` bottom with submission date and ID.

**Stage 5 exit gate:** ScholarOne submission confirmation email received.

---

## Fallback: *Bioinformatics Advances*

If *Bioinformatics* returns a desk rejection (or a verdict the user does not want to address):

- [ ] **F.1 Retarget cover letter address line**

  One edit in `docs/findings/cover_letter.md`: "The Editors, *Bioinformatics*" → "The Editors, *Bioinformatics Advances*". Remove the "declared fallback" sentence at the end of the venue rationale paragraph (it's redundant once submitted to BiA).

- [ ] **F.2 Upload to *Bioinformatics Advances* portal**

  https://mc.manuscriptcentral.com/bioinformaticsadv → same format, same manuscript PDF, same supplementary. No reformatting required.

---

## Dependency graph (summary)

```
Stage 0 (Zenodo release)
      │
      ├──▶ Stage 1 (benchmark: Floor 3) ─────────┐
      │                                           │
      ├──▶ Stage 2 (manuscript restructure) ─────┤ (Stage 2 consumes Stage 1 numbers at 2.5)
      │                                           │
      └──▶ Stage 3 (second-dataset validation) ──┤
                                                  │
                                                  ▼
                                      Stage 4 (supp + cover letter finalisation)
                                                  │
                                                  ▼
                                      Stage 5 (PDF + portal submission)
                                                  │
                                                  ▼
                                      (Fallback: Bioinformatics Advances)
```

## Critical-path estimate

- Stage 0: ~2 hours wall-clock (most of it waits on user web-UI actions)
- Stage 1: ~8 hours wall-clock (driven by 1,000-perm × 2-tool shell-out cost)
- Stage 2: ~3 hours wall-clock (prose editing)
- Stage 3: ~0.5 hours wall-clock (two reproduce-scripts)
- Stage 4: ~1 hour
- Stage 5: ~2 hours (PDF render + upload)

Parallelising 1/2/3 after Stage 0 puts the critical path at Stage 0 + Stage 1 + Stage 4 + Stage 5 ≈ 13 hours elapsed across sessions, assuming the user executes 0.7, 0.8, 0.9, 5.4 promptly.

## Outstanding user-action list (consolidated)

From Stages 0–5 above:
1. Verify GitHub repo is named `atman` (Stage 0.2)
2. Enable Zenodo-GitHub integration (Stage 0.7)
3. Publish GitHub Release against v1.0.0 tag (Stage 0.8)
4. Retrieve Zenodo DOI (Stage 0.9)
5. Fill ORCID, affiliation, funding, COI, CRediT on portal (Stage 4.3)
6. Name suggested + opposed reviewers (Stage 4.3)
7. ScholarOne portal upload (Stage 5.4)
