# Bioinformatics Submission — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship a `Bioinformatics` Original Paper submission for Atman — frozen Zenodo release, two-tool benchmark (OlinkAnalyze + limma), restructured method-first manuscript, cover letter, and supplementary package.

**Architecture:** Six parallel workstreams. Workstream A (release + DOI) is the unblocker — must land first. B (benchmark), C (manuscript), D (second dataset) run in parallel after A. E (supplementary) assembles after B and C. F (submission package) closes last. Bioinformatics Advances is the pre-declared fallback — zero reformat cost.

**Tech stack:** Rust workspace (`atman` CLI + `atman-core` lib), R 4.x for comparators (`OlinkAnalyze` 5.0.0, `limma` 3.66), Python drivers for pathway + phenotype passes, Docker for reproducible build, Zenodo-GitHub integration for DOI.

---

## Status snapshot (2026-04-17)

| Workstream | State |
|---|---|
| Task 0 (open questions) | **Done.** All three resolved on 2026-04-16. |
| A — release + DOI | **Files prepared, no commit/tag.** Cargo.toml v1.0.0, LICENSE-{MIT,APACHE}, CHANGELOG, README expansion, CI workflow, Dockerfile, `scripts/reproduce_paper.sh`. A7 (tag + push + Zenodo) awaits user. |
| B — benchmark | **Complete with real numbers.** OlinkAnalyze 5.0.0 + limma 3.66 executed; Table 1 + overlap + UpSet at `benchmarks/out/`. Capability matrix written. |
| C — manuscript | **Complete.** Retitle, abstract restructure, named-gap intro (cites [15,16]), command surface + Fig 1 DAG, benchmark subsection with containment framing (31/33, 44/44, 145/150), 6-figure consolidation, countable body 4370w. |
| D — second dataset | **D1 candidates written; D2/D3 deferred to post-review per PI decision on 2026-04-17.** Gisby 2021 (eLife) recommended if a revision request arrives. |
| E — supplementary | **Complete.** Tutorial (`docs/tutorial.md`), benchmark full methodology (`docs/findings/supp/benchmark_full.md` + `benchmark_raw/` TSVs), module membership (`docs/findings/supp/modules_full.{md,tsv}`). Module gap **closed 2026-04-17** via re-curation — modules.tsv written, atman module-trajectory re-run, latent axes + physiology regenerated, manuscript numbers updated (PC1 22.66%, PC2 20.66%). |
| F — submission | **F1 + F3 complete.** Cover letter at `docs/findings/cover_letter.md`; submission checklist at `docs/findings/submission_checklist.md`. F2 (ORCID, affiliation, funding, COI, CRediT — placeholders in manuscript) awaits user. F4 (portal submit) is user action. |

### Blocking items before submission

1. **v1.0.0 tag + Zenodo DOI** — user commit + push + tag + enable Zenodo + backfill DOI.
2. **GitHub repo rename** completed to `atman` — verify repo settings and badges.
3. **Author metadata stubs** — affiliation, ORCID, funding, COI in manuscript title page + Contact/Funding/Competing-interests sections (placeholders marked `USER-FILL`).
4. **Suggested reviewers** in cover letter (currently placeholder).

---

## Open questions to resolve before executing (Task 0)

Three factual uncertainties must be resolved before workstream C can proceed. Do these first.

### Task 0: Resolve open questions

**Files:** Conversation / repo inspection only. No code changes.

- [x] **0.1** Identify the authoritative manuscript draft. **RESOLVED 2026-04-16.**
  - Base: `docs/findings/2026-04-16-bioinformatics-manuscript.md` (renamed from the 3876w elife-labeled draft, structurally venue-agnostic).
  - Archived: the two shorter drafts moved to `docs/findings/archive/`.
- [x] **0.2** Verify test count. **RESOLVED 2026-04-16.**
  - Actual: 56 tests (verified via `cargo test --workspace --release`). README updated to cite 56. CI badge carries truth automatically.
- [x] **0.3** Confirm word-budget pressure. **RESOLVED 2026-04-16.**
  - Base is 3876w (77% of 5000 ceiling). ~1100w headroom exists for the command-surface subsection + benchmark subsection + named-gap intro paragraph. The spec's anticipated 500–1000w trim is not required unless C4+C5 overshoot.

**Commit:** No commit for Task 0 — these are decisions, not changes.

---

## Workstream A — Release + Zenodo DOI (UNBLOCKER)

Everything downstream cites the DOI. Land A before starting B/C/D.

### Task A1: Version bump

**Files:**
- Modify: `Cargo.toml:7`

- [ ] **A1.1** Bump `version = "0.1.0"` → `version = "1.0.0"`. Rationale per spec: cleanest signal of stable public release.
- [ ] **A1.2** Run `cargo build --workspace --release`. Expected: clean build.
- [ ] **A1.3** Run `cargo test --workspace --release`. Expected: all pass, count unchanged.
- [ ] **A1.4** Commit.
  ```bash
  git add Cargo.toml Cargo.lock
  git commit -m "chore(release): bump to v1.0.0 for paper release"
  ```

### Task A2: LICENSE file

**Files:**
- Create: `LICENSE-MIT`
- Create: `LICENSE-APACHE`
- Modify (if missing): top-level `README.md` license section

- [ ] **A2.1** Check whether `LICENSE`, `LICENSE-MIT`, or `LICENSE-APACHE` already exist at repo root. Cargo.toml declares dual MIT/Apache-2.0 but spec requires the files on disk.
- [ ] **A2.2** If missing, fetch upstream license texts (standard templates) and write both files.
- [ ] **A2.3** Commit.
  ```bash
  git add LICENSE-MIT LICENSE-APACHE
  git commit -m "chore: add LICENSE-MIT and LICENSE-APACHE"
  ```

### Task A3: CHANGELOG entry

**Files:**
- Create or modify: `CHANGELOG.md`

- [ ] **A3.1** Create `CHANGELOG.md` with Keep-a-Changelog format if absent.
- [ ] **A3.2** Add `## [1.0.0] — 2026-04-XX` section summarizing: reproduction base, DE (paired-t + moderated), asymmetry, robustness, module-trajectory, pathway ORA (universe-corrected), phenotype regression (null result).
- [ ] **A3.3** Commit.
  ```bash
  git add CHANGELOG.md
  git commit -m "docs: add CHANGELOG entry for v1.0.0"
  ```

### Task A4: README expansion

**Files:**
- Modify: `README.md`

- [ ] **A4.1** Add "Install" section: `git clone`, `cargo install --path crates/atman`, required R packages for benchmark.
- [ ] **A4.2** Add "Quickstart — reproduce the Dube results" — a single bash block running the full pipeline end-to-end on the example data.
- [ ] **A4.3** Add "Reproduce the paper" one-liner — wraps the quickstart plus enrichment, phenotype, heterogeneity, robustness passes.
- [ ] **A4.4** Add "Full command reference" section listing every `atman <subcommand>` with one-line description.
- [ ] **A4.5** Add CI status badge at top (GitHub Actions). Badge URL format: `https://github.com/kevinjoseph/atman/actions/workflows/<ci>.yml/badge.svg`. Confirm workflow name exists — if no CI workflow exists yet, this becomes a blocker for A5.
- [ ] **A4.6** Add placeholder Zenodo badge line (DOI filled in after A6).
- [ ] **A4.7** Commit.
  ```bash
  git add README.md
  git commit -m "docs: expand README with install, quickstart, command reference, badges"
  ```

### Task A5: CI workflow check

**Files:**
- Inspect: `.github/workflows/*.yml`
- Create if absent: `.github/workflows/ci.yml`

- [ ] **A5.1** Check whether a CI workflow exists. If not, create `.github/workflows/ci.yml` running `cargo fmt --check`, `cargo clippy --workspace -- -D warnings`, `cargo test --workspace --release`.
- [ ] **A5.2** Push and confirm green on main before tagging.
- [ ] **A5.3** Commit if created.
  ```bash
  git add .github/workflows/ci.yml
  git commit -m "ci: add workflow running fmt, clippy, tests"
  ```

### Task A6: Dockerfile

**Files:**
- Create: `Dockerfile`
- Create: `.dockerignore`

- [ ] **A6.1** Write two-stage `Dockerfile`: builder stage on `rust:1.75-slim` compiles the workspace; runtime stage on `debian:slim` includes the compiled binary, R + OlinkAnalyze + limma for benchmark reproduction, and Python + pandas + statsmodels for the scripts.
- [ ] **A6.2** Write `.dockerignore` excluding `target/`, `.git/`, `docs/findings/`.
- [ ] **A6.3** Build locally: `docker build -t atman:1.0.0 .`. Expected: builds clean.
- [ ] **A6.4** Run the Dube pipeline inside the container end-to-end:
  ```bash
  docker run --rm -v $(pwd)/example_data:/data -v $(pwd)/out_docker:/out atman:1.0.0 \
    atman ingest --platform olink-explore-ngs --parser dube --output-dir /out /data/dube_heat_2023/*.csv
  ```
  Expected: produces filtered NPX files byte-identical to the existing `out/` directory.
- [ ] **A6.5** Commit.
  ```bash
  git add Dockerfile .dockerignore
  git commit -m "build: add reproducible Dockerfile for paper release"
  ```

### Task A7: Frozen tag + Zenodo DOI

**Files:**
- Git tag only; Zenodo is external.

- [ ] **A7.1** **STOP. User confirmation required before tagging.** Confirm: repo state clean, all tests green on CI, README Zenodo placeholder present, v1.0.0 agreed.
- [ ] **A7.2** Tag and push.
  ```bash
  git tag -a v1.0.0 -m "Atman v1.0.0 — paper release"
  git push origin v1.0.0
  ```
- [ ] **A7.3** Confirm Zenodo-GitHub integration is enabled (`Settings → Integrations → Zenodo` on GitHub repo). If not, user enables it once, then deletes & re-pushes the tag to trigger archival.
- [ ] **A7.4** Retrieve DOI from Zenodo dashboard. Record in `docs/findings/zenodo_doi.txt`.
- [ ] **A7.5** Backfill DOI into README badge (Task A4.6).
- [ ] **A7.6** Commit the DOI backfill.
  ```bash
  git add README.md docs/findings/zenodo_doi.txt
  git commit -m "docs: add Zenodo DOI badge for v1.0.0"
  ```

**A complete. B/C/D can now start in parallel.**

---

## Workstream B — Benchmark (parallel after A)

Build the Table 1 comparison. Two comparators confirmed by literature search 2026-04-16:
- **OlinkAnalyze v4.5.0** — `olink_ttest`, `olink_lmer`, `olink_pathway_enrichment`
- **limma** — on the NPX matrix directly
- **DEP and "Olink Python" DROPPED** (off-label for NPX / does not exist)

### Task B1: Benchmark harness skeleton

**Files:**
- Create: `benchmarks/README.md`
- Create: `benchmarks/run_all.sh`
- Create: `benchmarks/out/` (gitignored or gitkept)

- [ ] **B1.1** Create `benchmarks/` directory with a README stating the comparator set, the Dube contrasts used (PT1-PR1, PR2-PR1, PT2-PT1, PT2-PR2), and the metrics (DE counts at q<0.05 and q<0.10, pairwise hit overlap, wall-clock, peak RSS).
- [ ] **B1.2** Write `run_all.sh` stub that invokes each comparator's script and writes per-tool TSVs to `benchmarks/out/`.
- [ ] **B1.3** Commit.
  ```bash
  git add benchmarks/README.md benchmarks/run_all.sh
  git commit -m "bench: scaffold benchmark harness"
  ```

### Task B2: OlinkAnalyze comparator

**Files:**
- Create: `benchmarks/olinkanalyze/run.R`
- Create: `benchmarks/olinkanalyze/README.md`

- [ ] **B2.1** Write `run.R` that loads raw Dube NPX (via `OlinkAnalyze::read_NPX`), runs `olink_ttest` for each of the four contrasts with BH-FDR, and writes results to `benchmarks/out/olinkanalyze_de.tsv` with columns `contrast, gene, logfc, pvalue, qvalue, tool, runtime_s`.
- [ ] **B2.2** Separately run `olink_lmer` for repeated-measures equivalence and write to `benchmarks/out/olinkanalyze_lmer.tsv`.
- [ ] **B2.3** Wrap with `/usr/bin/time -v` or `Sys.time()` + `pryr::mem_used()` for runtime and peak RSS. Record in output TSV.
- [ ] **B2.4** Sanity check: canonical HSPs (HSPB1, HSPA1A, DNAJB1, HSPG2) should come up as up-regulated in PT1-PR1 and PT2-PR2 at q<0.10.
- [ ] **B2.5** Commit.
  ```bash
  git add benchmarks/olinkanalyze/
  git commit -m "bench: OlinkAnalyze comparator (olink_ttest + olink_lmer)"
  ```

### Task B3: limma comparator

**Files:**
- Create: `benchmarks/limma/run.R`

- [ ] **B3.1** Write `run.R` that builds a wide NPX matrix from Dube output, fits `lmFit` + `eBayes` with a paired design (`~ subject + timepoint`), extracts `topTable` for each contrast, writes to `benchmarks/out/limma_de.tsv` in same schema as B2.
- [ ] **B3.2** Wrap with timing + memory instrumentation.
- [ ] **B3.3** Sanity check same as B2.4.
- [ ] **B3.4** Commit.
  ```bash
  git add benchmarks/limma/
  git commit -m "bench: limma comparator on NPX matrix"
  ```

### Task B4: Atman baseline harness

**Files:**
- Create: `benchmarks/atman/run.sh`

- [ ] **B4.1** Write `run.sh` that runs `atman de --test paired-t` and `--test moderated` on the four contrasts, wraps with `/usr/bin/time -v`, writes to `benchmarks/out/atman_de.tsv` in the same schema.
- [ ] **B4.2** Sanity check same as B2.4.
- [ ] **B4.3** Commit.
  ```bash
  git add benchmarks/atman/
  git commit -m "bench: atman baseline harness"
  ```

### Task B5: Aggregation + Table 1

**Files:**
- Create: `benchmarks/aggregate.py`
- Create: `benchmarks/table1.tsv`
- Create: `benchmarks/upset_de_overlap.pdf`

- [ ] **B5.1** Write `aggregate.py` that reads the three per-tool TSVs, computes:
  - DE count at q<0.05 and q<0.10 per tool per contrast
  - Pairwise Jaccard overlap of DE hit sets per contrast
  - Wall-clock and peak RSS per tool
  - Lines-of-code-to-reproduce per tool (count lines of the invocation scripts)
  - Writes `table1.tsv`.
- [ ] **B5.2** Write an UpSet plot (matplotlib-upsetplot) of DE hit overlap across the three tools for the acute-heat contrast (PT1-PR1). Save `upset_de_overlap.pdf` at 300dpi.
- [ ] **B5.3** Commit.
  ```bash
  git add benchmarks/aggregate.py benchmarks/table1.tsv benchmarks/upset_de_overlap.pdf
  git commit -m "bench: Table 1 + UpSet figure for manuscript"
  ```

### Task B6: Capability gap table

**Files:**
- Create: `benchmarks/capability_matrix.tsv`

- [ ] **B6.1** Build a rows-are-capabilities, cols-are-tools matrix covering: paired-t, moderated-t, LMM, ORA with universe correction, asymmetry, LOO rank stability, module trajectory, subject response burden, fingerprint, permutation null, deterministic single binary.
- [ ] **B6.2** Mark each cell ✓ / ✗ / partial. This becomes the differentiation column of Table 1 or a standalone supplementary table.
- [ ] **B6.3** Commit.
  ```bash
  git add benchmarks/capability_matrix.tsv
  git commit -m "bench: capability matrix (Atman vs OlinkAnalyze vs limma)"
  ```

---

## Workstream C — Manuscript restructure (parallel after A)

Depends on Task 0.1 (authoritative draft resolved) and ideally B5 (Table 1 data ready for C5). C1–C3, C6, C8 can proceed without B.

### Task C1: Retitle method-first

**Files:**
- Modify: the authoritative draft (TBD in Task 0.1). Call it `2026-04-16-bioinformatics-manuscript.md` below.

- [ ] **C1.1** Replace current title with a method-first line. Candidate: "Atman: a deterministic Rust pipeline for reproducible and heterogeneity-aware Olink proteomics reanalysis."
- [ ] **C1.2** Commit.
  ```bash
  git add docs/findings/2026-04-16-bioinformatics-manuscript.md
  git commit -m "manuscript: retitle method-first for Bioinformatics"
  ```

### Task C2: Restructure abstract

**Files:**
- Modify: `2026-04-16-bioinformatics-manuscript.md` (abstract section)

- [ ] **C2.1** Rewrite abstract to Bioinformatics-native shape with four labeled sections: Motivation, Results, Availability and Implementation, Contact.
- [ ] **C2.2** Results section retains the two-layer finding (reproduction + heterogeneity tooling); framing is tool capability, not heat physiology.
- [ ] **C2.3** Availability section cites Zenodo DOI (from A7) and GitHub repo.
- [ ] **C2.4** Check abstract word count against OUP limit (typically 150w for Bioinformatics Original Paper).
- [ ] **C2.5** Commit.
  ```bash
  git add docs/findings/2026-04-16-bioinformatics-manuscript.md
  git commit -m "manuscript: restructure abstract to motivation/results/availability/contact"
  ```

### Task C3: Named-gap introduction

**Files:**
- Modify: `2026-04-16-bioinformatics-manuscript.md` (introduction section)

- [ ] **C3.1** Add an explicit paragraph naming the comparators (OlinkAnalyze, limma) and the gap: "pooled DE is well-served by existing tools; subject-level heterogeneity-aware analysis with deterministic reproducibility is not."
- [ ] **C3.2** Add citations for OlinkAnalyze (CRAN) and limma (Ritchie et al. 2015) if not already present.
- [ ] **C3.3** Commit.
  ```bash
  git add docs/findings/2026-04-16-bioinformatics-manuscript.md
  git commit -m "manuscript: name comparators and state the capability gap in intro"
  ```

### Task C4: Command-surface Results subsection + Fig 1

**Files:**
- Modify: `2026-04-16-bioinformatics-manuscript.md` (Results, new subsection at top)
- Create: `docs/findings/figures/fig1_command_dag.pdf` (and source .tex, .py, or .dot)

- [ ] **C4.1** Write a short paragraph describing the command surface: `ingest`, `qc`, `matrix`, `fold-change`, `de` (with `paired-t` and `moderated` modes), `asymmetry`, `robustness`, `module-trajectory`, plus the Python driver scripts for ORA and phenotype regression.
- [ ] **C4.2** Build the command DAG schematic (raw NPX → ingest → qc → matrix → {fold-change, de} → {asymmetry, robustness, module-trajectory, enrichment, phenotype}). Render at 300dpi vector PDF. Tools: graphviz `dot` or matplotlib.
- [ ] **C4.3** Insert figure ref into manuscript as Fig 1.
- [ ] **C4.4** Commit.
  ```bash
  git add docs/findings/2026-04-16-bioinformatics-manuscript.md docs/findings/figures/fig1_command_dag.*
  git commit -m "manuscript: add command-surface subsection and Fig 1 DAG"
  ```

### Task C5: Benchmark Results subsection (depends on B5, B6)

**Files:**
- Modify: `2026-04-16-bioinformatics-manuscript.md` (Results)
- Copy into `docs/findings/figures/`: UpSet figure from B5.2

- [ ] **C5.1** Insert benchmark subsection with two paragraphs: (1) concordance — per-tool DE counts and overlap on the Dube contrasts, stating Atman agrees with OlinkAnalyze and limma on the canonical heat-shock signal; (2) differentiation — reference Table 1 and the capability matrix, stating which outputs the comparators cannot produce.
- [ ] **C5.2** Insert Table 1 (from `benchmarks/table1.tsv`) and the UpSet figure.
- [ ] **C5.3** Commit.
  ```bash
  git add docs/findings/2026-04-16-bioinformatics-manuscript.md docs/findings/figures/upset_de_overlap.pdf
  git commit -m "manuscript: add benchmark subsection, Table 1, UpSet figure"
  ```

### Task C6: Figure consolidation 8 → 6

**Files:**
- Modify: `2026-04-16-bioinformatics-manuscript.md`
- Modify: `docs/findings/figures/*.pdf`

- [ ] **C6.1** Inventory the 8 current main figures from the authoritative draft. Decide consolidation:
  - Fig 1 = command DAG (new, C4)
  - Fig 2 = subject architecture (merge of current Fig 1 + 2)
  - Fig 3 = DE + canonical heat-shock markers
  - Fig 4 = robustness (LOO + rank stability + threshold sensitivity)
  - Fig 5 = variable proteins + fingerprints (merge of current Fig 5 + 8)
  - Fig 6 = latent/physiology
  - Benchmark/UpSet → separate figure or part of Fig 3.
- [ ] **C6.2** Move displaced panels to supplementary.
- [ ] **C6.3** Renumber figure refs throughout manuscript.
- [ ] **C6.4** Commit.
  ```bash
  git add docs/findings/2026-04-16-bioinformatics-manuscript.md docs/findings/figures/
  git commit -m "manuscript: consolidate to 6 main figures, move overflow to supp"
  ```

### Task C7: Trim to word budget

**Files:**
- Modify: `2026-04-16-bioinformatics-manuscript.md`

- [ ] **C7.1** Run word count. Compute headroom = 5000 - current. If Task 0.3 confirmed the 3876w base, headroom is ~1124w; the added sections (C4 + C5, roughly 400–600w) should fit without trimming.
- [ ] **C7.2** If over budget, trim per spec suggestions:
  - Condense fingerprint biology paragraph
  - Shorten physiology-mapping discussion
  - Move second half of Limitations to supplementary
- [ ] **C7.3** Final word count under 5000. Verify.
- [ ] **C7.4** Commit.
  ```bash
  git add docs/findings/2026-04-16-bioinformatics-manuscript.md
  git commit -m "manuscript: trim to Bioinformatics word budget"
  ```

---

## Workstream D — Second dataset (optional, parallel after A)

Pre-empts the "only one dataset" reviewer comment. Abbreviated treatment: one paragraph + one supp figure. Skip if time-constrained.

### Task D1: Identify candidate dataset

**Files:** None yet.

- [ ] **D1.1** Search GEO and ArrayExpress for public Olink/NPX datasets with repeated-measures structure, n ≥ 8 within-subject. Criteria: CC-BY or equivalent, raw NPX available, at least two timepoints per subject. Candidates to evaluate first: any UK Biobank Olink sub-study release, published Olink longitudinal cohorts from 2023–2025.
- [ ] **D1.2** Pick one. Record choice in `docs/findings/second_dataset_choice.md` with citation and access URL.

### Task D2: Ingest + run pipeline

**Files:**
- Create: `example_data/<dataset>/` (gitignored if large; checked in if small)
- Create: `scripts/run_second_dataset.sh`

- [ ] **D2.1** If the dataset's NPX format matches Dube, reuse the parser. If not, add a minimal parser to `crates/atman-core/src/ingest/`. New parser gets full TDD coverage.
- [ ] **D2.2** Run ingest → qc → matrix → fold-change → de → asymmetry → robustness.
- [ ] **D2.3** Verify heterogeneity commands produce meaningful output (non-trivial asymmetry, stable LOO ranks at top-k ≥ 10).
- [ ] **D2.4** Commit pipeline scripts.
  ```bash
  git add scripts/run_second_dataset.sh crates/atman-core/src/ingest/
  git commit -m "feat(ingest): parser for <dataset>; run heterogeneity pipeline"
  ```

### Task D3: One paragraph + one supp figure

**Files:**
- Modify: `2026-04-16-bioinformatics-manuscript.md` (Results, short addition)
- Create: `docs/findings/figures/figS_second_dataset.pdf`

- [ ] **D3.1** Write one paragraph summarizing: dataset, n, contrast, concordance check (e.g., does Atman recover the dataset's published DE hits?), and one novel heterogeneity observation.
- [ ] **D3.2** Build one supplementary figure — typically the heterogeneity archetype plot or LOO stability.
- [ ] **D3.3** Commit.

---

## Workstream E — Supplementary materials (after B, C)

### Task E1: Tutorial walkthrough

**Files:**
- Create: `docs/tutorial.md`

- [ ] **E1.1** Write command-by-command walkthrough of the Dube pipeline with expected output excerpts at each step. Link from manuscript Availability section.
- [ ] **E1.2** Commit.
  ```bash
  git add docs/tutorial.md
  git commit -m "docs: add tutorial walkthrough"
  ```

### Task E2: Full benchmark tables

**Files:**
- Create: `docs/findings/supp/benchmark_full.md`
- Copy: `benchmarks/out/*.tsv` → `docs/findings/supp/benchmark_raw/`

- [ ] **E2.1** Write the full benchmark methodology (versions, seeds, machine spec) and the per-tool output tables from B2/B3/B4.
- [ ] **E2.2** Commit.

### Task E3: Module membership tables

**Files:**
- Create: `docs/findings/supp/modules_full.tsv`

- [ ] **E3.1** Export the full module definitions used by `module-trajectory` into a single TSV. Cite source of each module.
- [ ] **E3.2** Commit.

### Task E4: Extended biology

**Files:**
- Create: `docs/findings/supp/extended_biology.md`

- [ ] **E4.1** Move any fingerprint biology / physiology trimmed in C7.2 into this supplementary document, referenced from the main manuscript.
- [ ] **E4.2** Commit.

---

## Workstream F — Submission package (last)

All of A, B, C, E complete before starting F.

### Task F1: Cover letter

**Files:**
- Create: `docs/findings/cover_letter.md`

- [ ] **F1.1** Write 2–3 paragraph cover letter: (1) what the software does, (2) why Bioinformatics is the right venue, (3) why the benchmark and deterministic reanalysis justify the page budget, (4) code + data deposited at Zenodo DOI.
- [ ] **F1.2** Add 2–3 suggested reviewers from the Olink/proteomics tools community. Exclude anyone with a conflict.
- [ ] **F1.3** Commit.
  ```bash
  git add docs/findings/cover_letter.md
  git commit -m "manuscript: cover letter + reviewer suggestions"
  ```

### Task F2: ORCID + formatting

**Files:**
- Modify: `2026-04-16-bioinformatics-manuscript.md` title page

- [ ] **F2.1** Confirm corresponding author ORCID on file. Add to title page.
- [ ] **F2.2** When rendering the submission PDF: double-spaced, 12pt, numbered lines.
- [ ] **F2.3** Confirm all figures render at 300dpi or vector.
- [ ] **F2.4** Commit any manuscript changes.
  ```bash
  git add docs/findings/2026-04-16-bioinformatics-manuscript.md
  git commit -m "manuscript: ORCID and submission formatting"
  ```

### Task F3: Pre-submission checklist

**Files:**
- Create: `docs/findings/submission_checklist.md`

- [ ] **F3.1** Build a checklist derived from current OUP Bioinformatics Instructions to Authors: word count ≤ 5000 (+20% tolerance), abstract in required shape, cover letter present, ORCID, 300dpi figs, data/code availability statement, competing interests, funding statement.
- [ ] **F3.2** Walk through every item and tick. Flag any miss for the user.
- [ ] **F3.3** Commit.

### Task F4: Submit

- [ ] **F4.1** **STOP. User submits via OUP portal.** Claude does not submit.

---

## Fallback — Bioinformatics Advances

If Bioinformatics desk-rejects or returns with verdict the user doesn't want to address:

- [ ] **FB.1** Confirm near-zero reformat cost. Advances accepts the same format and scope is more permissive.
- [ ] **FB.2** Reuse cover letter with target-journal line adjusted.
- [ ] **FB.3** Submit.

---

## Dependency graph (summary)

```
Task 0 (decisions)
  └─> A (release + Zenodo DOI) — BLOCKING UNBLOCKER
         ├─> B (benchmark)          ─┐
         ├─> C (manuscript)          ├─> E (supplementary) ─> F (submission)
         └─> D (second dataset, opt) ─┘
```

## Open questions deferred to execution

- Test count reconciliation (Task 0.2) — if 658 is aspirational, spec's CI-badge language needs adjusting in Task A4.5 to cite the real count.
- Second dataset choice (D1) — user may know of a preferred public Olink cohort; if so, skip the GEO search.
- Reviewer suggestions (F1.2) — requires user input from the Olink community.
