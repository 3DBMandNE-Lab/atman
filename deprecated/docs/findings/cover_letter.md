# Cover letter — Atman Bioinformatics submission

**To:** The Editors, *Bioinformatics* (Oxford University Press)

**From:** Kevin Joseph (<kevin.joseph08@gmail.com>), ORCID [pending]

**Date:** 2026-04-XX

**Re:** Original Paper submission — *Atman: a multiscale, stability-aware framework for small-cohort quantitative proteomics*

---

Dear Editors,

Please find enclosed an Original Paper submission introducing **Atman**, a
local-first Rust engine for Olink NPX proteomics. The software closes a
specific gap in the current toolchain: established packages — notably
OlinkAnalyze (CRAN) and limma — provide well-validated pooled differential
abundance on NPX data, but no existing tool emits subject-level heterogeneity
and robustness artefacts (asymmetry, leave-one-subject-out rank stability,
response burden, fingerprint, module trajectory) from a single deterministic
engine. For small-cohort human proteomics, where pooled significance lists
compress biologically meaningful person-level structure into a group-average
narrative, this gap materially changes interpretation.

The manuscript makes three contributions appropriate for *Bioinformatics*.
First, a head-to-head benchmark on the Dube et al. (2023, *Scientific Data*)
heat-acclimation dataset against OlinkAnalyze 5.0.0 (paired `olink_ttest` and
`olink_lmer`) and limma 3.66.0 (paired `lmFit` / `eBayes` on the NPX matrix),
showing concordant acute-dominant structure across all three tools and
substantial hit-set overlap at q<0.10. Second, a capability matrix making
the differentiation claim concrete: eight heterogeneity / robustness /
reproduction outputs that neither comparator emits. Third, a reanalysis of
the Dube dataset itself, demonstrating that the heterogeneity layer — all
four acute/adaptation responder classes populated, recurrent fingerprint
proteins (DAND5, GH1, ITGAL, AGBL2), and a permutation-supported excess
dispersion null at p<0.002 — changes the biological narrative beyond what a
pooled significance list can carry.

Reproducibility claims are backed by code and data depositions. Source code
is archived at Zenodo under the v1.0.0 release tag (DOI pending deposition,
to be backfilled at revision); a Dockerfile ships with the release for
hermetic rebuilds; a single `scripts/reproduce_paper.sh` invocation
regenerates every main-text artefact from the Dube raw inputs in roughly one
minute. Atman reproduces Dube's published filtered NPX files cell-exactly
and log2 fold-change tables to 1.05 × 10⁻¹⁵ — the machine-epsilon floor for
IEEE 754 arithmetic. All benchmark driver scripts, per-tool output tables,
capability matrix, and module/tutorial supplementary material are archived
in the same repository.

We believe *Bioinformatics* is the right venue because the work is primarily
a methods/software contribution with a fully reproducible case-study
reanalysis: the Original Paper format allows the space required to present
both the benchmark and the biological findings that motivate the
heterogeneity-aware approach. We have prepared the submission with
*Bioinformatics Advances* as a declared fallback target requiring no
reformatting.

Suggested reviewers from the Olink/proteomics tools and heat-acclimation
communities are listed in the submission portal. We have no competing
interests to declare. The paper has not been submitted elsewhere and is not
under consideration at another journal.

Thank you for considering our submission.

Yours sincerely,

Kevin Joseph

---

## Suggested reviewers (portal entry)

[User to complete — 2–3 names from the Olink/proteomics tools community; exclude anyone with a direct collaborative conflict.]

## Opposed reviewers

[User to complete if applicable — otherwise leave blank.]
