# Pre-submission checklist — OUP Bioinformatics Original Paper

Work through every item before submitting via the ScholarOne portal. Each
item cites the responsible artefact in the repository.

## Manuscript body

- [x] Method-first title — *Atman: a multiscale, stability-aware framework for small-cohort quantitative proteomics* (`docs/findings/2026-04-16-bioinformatics-manuscript.md:1`)
- [x] Abstract in Motivation / Results / Availability and Implementation / Contact structure (~175 words, within typical Bioinformatics cap) (`docs/findings/2026-04-16-bioinformatics-manuscript.md` abstract section)
- [x] Named-gap paragraph in Introduction citing OlinkAnalyze and limma
- [x] Results subsection: Command surface (Fig 1 command DAG)
- [x] Results subsection: Benchmark (concordance + capability gap, real Table 1 numbers)
- [x] 6 main figures (command DAG, subject architecture, DE, robustness, variable proteins + fingerprints, latent + physiology)
- [x] Main-body word count under the 5000 cap (countable body = 4168 words; 832-word headroom at time of writing)
- [x] Abstract, references excluded from the 5000-word cap
- [x] Ethics statement (secondary analysis, no new consent required)
- [x] Data availability (Figshare + Zenodo DOI placeholder)
- [x] Code availability (GitHub URL + Zenodo DOI placeholder + Dockerfile reference)
- [x] Figure availability
- [x] References (15 entries, standard numbered format)

## Release & reproducibility

- [x] Cargo.toml version bumped to `1.0.0` (`Cargo.toml:6`)
- [x] `LICENSE-MIT` and `LICENSE-APACHE` present at repo root
- [x] `CHANGELOG.md` with v1.0.0 entry
- [x] `README.md` with install, quickstart, reproduce-paper one-liner, command reference, CI badge, Zenodo placeholder
- [x] `.github/workflows/ci.yml` runs fmt + clippy + test on push/PR
- [x] `Dockerfile` + `.dockerignore` for reproducible rebuilds
- [x] `scripts/reproduce_paper.sh` end-to-end one-liner
- [x] 75 tests pass on `cargo test --workspace --release`
- [ ] GitHub repo is named `atman` **(verify in GitHub settings)**
- [ ] `git tag -a v1.0.0 -m "Atman v1.0.0 — paper release"` + `git push origin v1.0.0` **(user action, after commit)**
- [ ] Zenodo-GitHub integration enabled in repo settings **(user action)**
- [ ] Zenodo DOI retrieved and backfilled into `README.md` + manuscript Data/Code availability **(user action)**
- [ ] CI workflow green on `main` after push **(user action — verify after push)**

## Benchmark & supplementary

- [x] Table 1 real numbers (`benchmarks/out/table1.tsv`)
- [x] Pairwise overlap table (`benchmarks/out/table1_overlap.tsv`)
- [x] UpSet figure rendered (`benchmarks/out/upset_de_overlap.pdf`)
- [x] Capability matrix (`benchmarks/capability_matrix.tsv`)
- [x] Benchmark supplementary doc (`docs/findings/supp/benchmark_full.md`) with full methodology + raw per-tool TSVs in `docs/findings/supp/benchmark_raw/`
- [x] Within-Atman paired-t vs moderated comparison (`benchmarks/out/atman_pairedt_vs_moderated.tsv`) showing perfect sign + rank concordance
- [x] Tutorial walkthrough (`docs/tutorial.md`)
- [x] Module membership supplementary — **resolved 2026-04-17 via path 2 (re-curation)**. `docs/findings/heterogeneity/modules.tsv` + `docs/findings/supp/modules_full.{tsv,md}` + regenerated `module_trajectory_scores.tsv`, `latent_axes.tsv`, `physiology_mapping.tsv`. Manuscript numerical claims updated (PC1 22.66%, PC2 20.66%, PC2~acc_dTbody r=−0.854, PC1~d1_dSSNA r=−0.725). v1 archived under `*_v1_archived.tsv`.
- [x] Second-dataset generalizability executed (2026-04-17). Gisby 2021 Olink Target 96 COVID dialysis cohort ingested via `scripts/ingest_gisby.py` and analysed through the same DAG without atman changes. 23 patients with paired early/late samples; 120 DE hits at q<0.05; LOO sign match 0.972, top-20 Jaccard 0.911, 82% hit retention. Documented in `docs/findings/supp/gisby_generalizability.md`. Manuscript now includes a dedicated Results subsection on generalizability with biology and scaling claims grounded in real numbers.

## Submission package

- [x] Cover letter (`docs/findings/cover_letter.md`) — 2–3 paragraphs, states what the software does, venue rationale, deposition statement
- [ ] Suggested reviewers (2–3 names from Olink/proteomics tools community) **(user action)**
- [ ] Opposed reviewers — optional **(user action)**
- [ ] Corresponding-author ORCID on file and entered in portal **(user action)**
- [ ] Affiliation statement in title page **(user action — add to manuscript title page)**
- [ ] Funding statement filled in **(user action — currently stub)**
- [ ] Competing-interests declaration filled in **(user action — currently stub)**
- [ ] CRediT author-contributions taxonomy filled in **(user action — currently stub)**

## PDF formatting for portal upload

- [ ] Manuscript PDF rendered double-spaced, 12pt, numbered lines **(user action at PDF build time)**
- [ ] All figures at 300 dpi or vector PDF
  - Fig 1 (pipeline schematic, matplotlib swim-lane): PDF vector, PNG 300 dpi — `docs/findings/figures/fig1_pipeline.{pdf,png}` ✓
  - Fig 2a–c, Fig 3a–c, Fig 4a–b, Fig 5a–b — each panel exists as standalone PDF+PNG in `docs/findings/figures/` (300 dpi PNG, vector PDF); Fig S4 (demoted robustness geometry) in same directory.
- [ ] Table 1 rendered legibly (either typeset or as a figure/table image at 300 dpi)

## Manuscript claim integrity — final sweep

- [x] Reference [16] cleanup (clusterProfiler mention stripped; references list is clean at 15 items)
- [x] All Fig references use the 5-figure main scheme (Fig 1, 2a–c, 3a–c, 4a–b, 5a–b) plus Fig S1–S4 supplementary; robustness geometry demoted Fig 3 → Fig S4
- [x] Table 1 numbers in manuscript body match `benchmarks/out/table1_counts.tsv`
- [x] OlinkAnalyze version (5.0.0) consistent across manuscript, benchmark README, and supplementary
- [x] Test count (56) consistent in README and manuscript-adjacent claims

## Submit

Once every `[ ]` in this checklist is ticked, submit via the ScholarOne
portal for *Bioinformatics*. If desk-rejected or returned with a verdict
the user does not want to address, the manuscript is ready to redirect to
*Bioinformatics Advances* with no reformat required.
