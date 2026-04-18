# Candidate second datasets for the Bioinformatics submission

**Status 2026-04-17:** Deferred from the initial submission by PI decision.
The single most common reviewer comment on Bioinformatics tool papers is
"please demonstrate on an independent dataset." This document pre-identifies
candidates so that if the comment comes back at review, a one-paragraph
addition can close it in days rather than weeks. Pick a candidate here and
stage the parser work in advance of the first revision request.

Selection criterion (from plan D1): public Olink/NPX dataset with
repeated-measures within-subject structure, n ≥ 8, open license (CC BY or
equivalent), raw NPX available.

## Shortlist — require user pick before ingesting

### Candidate 1 — Gisby et al. 2021 (eLife)

- Title: *Longitudinal proteomic profiling of dialysis patients with COVID-19
  reveals markers of severity and predictors of death.*
- DOI: 10.7554/eLife.64827
- Code/data: https://github.com/jackgisby/longitudinal_olink_proteomics
- Platform: Olink Target 96 (**not** Olink Explore NGS like Dube; same NPX
  format, different parser quirks possible).
- Design: repeated within-subject sampling over the course of COVID-19 admission
  (up to 20+ timepoints per subject for some). **Truly longitudinal**, not just
  pre/post.
- n: ~55 patients (see paper for exact).
- Strengths: extensively reused in the literature, peer-reviewed, open.
- Risks: Target 96 parser not yet implemented — would need a new ingest
  adapter. Non-trivial engineering but legitimate extension.

### Candidate 2 — Long COVID IRIS Study (Zenodo)

- DOI: 10.5281/zenodo.13363117
- Platform: Olink Target 96 Inflammation + Immune Response panels.
- Design: reported as longitudinal; exact sampling cadence to confirm by
  inspecting the Zenodo manifest.
- Strengths: self-contained Zenodo deposition, explicit CC license.
- Risks: Target 96 again; need to verify repeated-measures structure meets the
  n ≥ 8 within-subject criterion.

### Candidate 3 — UK Biobank Olink (UKB-PPP)

- Platform: Olink Explore 3072 (**matches Atman's current target**).
- n: >50,000, but primarily cross-sectional baseline samples; follow-up
  Olink re-measures exist in a small subset.
- Strengths: same platform as Dube; enormous cohort.
- Risks: requires UKB Tier-2 access, not freely reusable. **Disqualified for
  a same-day supplementary analysis** — the access process alone is weeks.

## Recommendation

**Candidate 1 (Gisby 2021)** gives the strongest replication story — peer-reviewed,
extensively longitudinal, different clinical domain (COVID vs heat stress)
which strengthens the "method generalizes" claim. Cost: one new Target 96
ingest adapter (~0.5 day of work).

**Fallback to Candidate 2** if the user wants to avoid writing a new parser;
the Zenodo structure may already be close enough to Dube NPX CSV to reuse
the existing Dube-adjacent parser with minor edits.

User action: pick one and I will ingest + produce one paragraph + one
supplementary figure.
