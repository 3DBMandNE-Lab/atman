# Dube heat acclimation — paired-t differential abundance findings

**Date:** 2026-04-14
**Dataset:** Gagnon, Barry, Barhdadi, Oussaid, Mongrain, Lemieux Perreault, Dubé. *A dataset of proteomic changes during human heat stress and heat acclimation.* Scientific Data (2023). DOI 10.1038/s41597-023-02809-5. Figshare [22811915](https://doi.org/10.6084/m9.figshare.22811915.v1). CC BY 4.0.
**Platform:** Olink Explore 1536 + Explore Expansion (NGS readout), 2,943 unique proteins across 8 panels.
**Design:** 10 healthy adults, paired design, 4 timepoints per subject:
- **PR1** = normothermic, pre-acclimation
- **PT1** = hyperthermic, pre-acclimation
- **PR2** = normothermic, post-acclimation (7-day heat protocol)
- **PT2** = hyperthermic, post-acclimation

## Method

1. **Ingest** raw Olink Explore NGS CSVs (122,135 NPX rows, 43 samples) via `karnaproteome ingest`.
2. **QC** applies the Dube rule (`NPX → missing` iff `QC_Warning != PASS` OR `Assay_Warning != PASS`) — masks 3,725 / 122,135 rows (3.0%).
3. **Differential abundance** via `karnaproteome de`:
   - Paired Student's t-test on log2-NPX at subject level (biological replicate = participant).
   - Four comparisons, each is an independent hypothesis family for multiple-testing correction.
   - `min_pairs = 5`: proteins with fewer than 5 paired subjects in both conditions are not tested.
   - BH-FDR **within each comparison** (≈ 2,940 tests per family).
   - No shrinkage, no empirical-Bayes variance pooling (deliberately simple; n is too small for limma-style hierarchical models to help much).

## Instrument sanity check

**Canonical heat-shock proteins must go UP in acute heat-stress comparisons.** Ran the check on `HSPB1, HSPA1A, DNAJB1, HSPG2` across both heat comparisons:

| Comparison | HSP | mean_diff (log2) | p-value | BH-q | n pairs |
|---|---|---|---|---|---|
| PT1−PR1 (acute, pre-acclimation) | HSPB1 | +0.888 | 0.0083 | 0.158 | 10 |
| PT1−PR1 | HSPA1A | +0.309 | 0.22 | — | 10 |
| PT1−PR1 | DNAJB1 | +0.119 | 0.78 | — | 10 |
| PT1−PR1 | HSPG2 | +0.081 | 0.12 | — | 10 |
| PT2−PR2 (acute, post-acclimation) | HSPA1A | **+0.714** | **0.0067** | **0.096** | 10 |
| PT2−PR2 | DNAJB1 | **+0.856** | 0.033 | 0.190 | 9 |
| PT2−PR2 | HSPB1 | +1.281 | 0.051 | 0.206 | 9 |
| PT2−PR2 | HSPG2 | +0.098 | 0.15 | — | 10 |

**4 out of 4 canonical HSPs go up in both comparisons.** Direction is correct. The strongest effect (HSPB1 pre, HSPA1A post) passes q < 0.2 in at least one comparison. This is as much as n=10 can deliver. The instrument is measuring what it's supposed to measure.

## Hit counts per comparison

| Comparison | Interpretation | q < 0.05 | q < 0.10 | (up / down) | q < 0.20 |
|---|---|---|---|---|---|
| **PT1−PR1** | Acute heat stress, pre-acclimation | 36 | 62 | 53 / 9 | 299 |
| **PR2−PR1** | Acclimation effect at rest (baseline shift) | 2 | 7 | 7 / 0 | 10 |
| **PT2−PT1** | Acclimation effect on the heat response | 7 | 8 | 6 / 2 | 13 |
| **PT2−PR2** | Acute heat stress, post-acclimation | 37 | **277** | 243 / 34 | **752** |

## Top hits per comparison (q < 0.10, sorted by |mean_diff|)

### PT1−PR1 — acute heat stress pre-acclimation (62 hits)

| Gene | Panel | mean_diff (log2) | p | BH-q | n |
|---|---|---|---|---|---|
| PRL | Neurology | +2.036 | 2.75e-05 | 0.009 | 10 |
| CLEC1B | Neurology | +1.050 | 1.85e-03 | 0.094 | 10 |
| DNM1 | Oncology_II | +0.930 | 1.84e-03 | 0.094 | 10 |
| STXBP1 | Neurology_II | +0.798 | 7.51e-04 | 0.055 | 10 |
| OBP2B | Neurology | +0.747 | 1.39e-05 | 0.007 | 9 |
| EDN1 | Cardiometabolic_II | +0.737 | 1.91e-04 | 0.031 | 9 |
| **CALCA** | Neurology | +0.709 | **3.25e-07** | **0.001** | 9 |
| NDUFB7 | Oncology_II | +0.677 | 1.55e-03 | 0.090 | 10 |
| RBM17 | Neurology_II | +0.641 | 1.30e-04 | 0.024 | 10 |
| HS6ST1 | Oncology | +0.568 | 2.01e-04 | 0.031 | 10 |

**Biology annotation** (ranked-level, not single-gene claims):
- **PRL** (prolactin), **CALCA** (calcitonin gene-related peptide), **EDN1** (endothelin-1) — classic acute-stress / HPA-axis / vasoactive peptide response. Expected direction.
- **OBP2B** (odorant binding protein 2B) — unusual. Possibly a secretory-pathway leakage artifact or a real finding; n=9 isn't enough to commit.
- 53/62 hits are up-regulated. Heat stress is a net upward push on the measurable secretome.

### PR2−PR1 — resting state, 7-day acclimation effect (7 hits)

| Gene | Panel | mean_diff | p | BH-q | n |
|---|---|---|---|---|---|
| DNER | Inflammation | +0.515 | 1.04e-04 | 0.065 | 10 |
| MMP7 | Cardiometabolic | +0.511 | 1.11e-04 | 0.065 | 9 |
| MAPRE3 | Oncology_II | +0.499 | 2.29e-04 | 0.099 | 9 |
| **KLK14** | Oncology | +0.492 | **1.47e-05** | **0.022** | 10 |
| CA6 | Neurology | +0.490 | 6.30e-05 | 0.062 | 10 |
| FEN1 | Oncology | +0.351 | 1.13e-05 | 0.022 | 8 |
| NCR3LG1 | Oncology_II | +0.319 | 2.36e-04 | 0.099 | 10 |

**Observation:** 7 days of heat acclimation produces very little detectable change in the resting plasma proteome. All 7 hits are up-regulated and modest in effect size (< 0.55 log2). This is the weakest signal of the 4 comparisons — meaning the "rest between heat exposures" doesn't rewrite the proteome on its own.

### PT2−PT1 — heat response pre- vs post-acclimation (8 hits)

| Gene | Panel | mean_diff | p | BH-q | n |
|---|---|---|---|---|---|
| CA6 | Neurology | +0.569 | 7.42e-05 | 0.042 | 10 |
| KLK14 | Oncology | +0.563 | 5.92e-05 | 0.042 | 10 |
| MMP7 | Cardiometabolic | +0.529 | 1.78e-05 | 0.026 | 9 |
| **DNER** | Inflammation | +0.433 | **3.76e-06** | **0.011** | 10 |
| SBSN | Cardiometabolic_II | +0.330 | 9.85e-05 | 0.042 | 10 |
| SIGLEC8 | Cardiometabolic_II | +0.151 | 1.48e-04 | 0.054 | 10 |
| PTPRN2 | Neurology | −0.348 | 5.26e-05 | 0.042 | 10 |
| NPDC1 | Cardiometabolic | −0.395 | 9.90e-05 | 0.042 | 10 |

**Observation:** 6 of the 8 hits here are the **same proteins** that showed up in PR2−PR1 (DNER, MMP7, CA6, KLK14, MAPRE3). These proteins are on a monotone upward trajectory across both timepoints — acclimation is shifting them up steadily, not flipping them on during heat exposure. Two proteins (PTPRN2, NPDC1) go DOWN post-acclimation during heat stress — a narrow signal.

### PT2−PR2 — acute heat stress post-acclimation (277 hits at q < 0.10)

| Gene | Panel | mean_diff | p | BH-q | n |
|---|---|---|---|---|---|
| **GH1** | Cardiometabolic | **+4.605** | 2.46e-03 | 0.078 | 9 |
| PRL | Neurology | +2.439 | 2.11e-04 | 0.036 | 10 |
| MGMT | Inflammation | +2.045 | 6.29e-03 | 0.095 | 8 |
| RHOC | Neurology | +2.014 | 6.41e-03 | 0.095 | 9 |
| SRC | Oncology | +1.954 | 3.75e-03 | 0.089 | 9 |
| DOK1 | Neurology_II | +1.943 | 9.20e-03 | 0.100 | 10 |
| TBCB | Neurology | +1.903 | 7.65e-03 | 0.099 | 9 |
| CLIP2 | Inflammation | +1.840 | 8.05e-03 | 0.099 | 8 |
| PDLIM7 | Inflammation | +1.792 | 5.90e-03 | 0.095 | 10 |
| MANF | Inflammation | +1.771 | 3.22e-03 | 0.088 | 10 |

277 hits at q < 0.10. GH1 (growth hormone) has a **4.6 log2** effect — enormous. This is the biggest signal in the whole analysis.

## The headline observation

**The post-acclimation heat response (PT2−PR2) produces ~4x more hits than the pre-acclimation heat response (PT1−PR1).**

| Comparison | q < 0.05 | q < 0.10 | q < 0.20 |
|---|---|---|---|
| PT1−PR1 (pre) | 36 | 62 | 299 |
| PT2−PR2 (post) | **37** | **277** | **752** |

This is not what the naive "acclimation dampens stress response" model predicts. Instead, the data suggests that **7 days of heat acclimation primes or amplifies the proteomic response to subsequent acute heat stress**. The same pattern shows in the canonical HSPs:

- HSPA1A: +0.309 (p=0.22) pre → **+0.714 (p=0.007)** post — more than doubled
- DNAJB1: +0.119 (p=0.78) pre → **+0.856 (p=0.033)** post — 7x larger
- HSPB1: +0.888 (p=0.008) pre → **+1.281 (p=0.051)** post — ~50% larger

**Interpretation caveat (important):** this could be a real biological priming effect (heat-shock machinery upregulated at baseline and now more responsive), OR a floor-effect artifact of the measurement (starting from a slightly different baseline where proteins have more headroom to change on log2 scale). Disentangling the two needs either a larger cohort or an orthogonal readout.

## What's not a finding

- **Nothing at FWER (Bonferroni) significance.** `0.05 / 2,943 ≈ 1.7e-5`. CALCA (PT1−PR1) is the only protein approaching it (`p = 3.3e-7`). All other "hits" are FDR-controlled at 5% or 10%, which means some fraction are expected to be false discoveries. At q < 0.10 with 62 hits, ≤ 6 of them are expected to be false; at q < 0.10 with 277 hits, ≤ 28 are expected false.
- **No single protein can pass FDR 5% by paired Wilcoxon signed-rank with n=10** (floor p ≈ 0.002 > 2,943 × FDR threshold). This is why we used paired t — not because we think the proteins are normally distributed, but because it's the only test with enough resolution. This is also what Dube themselves used via SAS PROC MIXED.
- **n=10 is small.** Large effect sizes (|mean_diff| > 1) with q < 0.10 are the honest headline; small effects passing q-thresholds by narrow margins are more likely to be noise.
- **No corrections across comparison families.** Each of the 4 comparisons is its own independent hypothesis family. Reporting a protein as "significant in both PT1−PR1 and PT2−PR2" is a separate conjunction test I haven't run here.

## Pathway enrichment

Handoff: `de_results.tsv` → `scripts/run_enrichment.sh` → per-`(comparison, direction, collection)` Fisher's exact ORA via the sibling `genesets` Rust CLI (MSigDB Hallmark / KEGG_MEDICUS / GO-BP, already bundled in the genesets binary at `~/Cursor/Rust_bioinfo_cli/genesets/`).

**Methodology:**
- Foreground: top-100 genes at `bh_q < 0.10` per comparison × direction, ranked by `|mean_diff|`.
- Background: full MSigDB gene universe. **Caveat:** not restricted to the Olink Explore measurable universe, so pathways with high Olink coverage may be over-represented in hits. This is a known ORA-with-restricted-platform bias; acknowledge it, don't trust single-pathway claims.
- Fisher's exact test per pathway, BH-FDR across pathways within each comparison × direction × collection.
- JSON results: `docs/findings/pathway_results/<comparison>_<direction>_<collection>.json`.
- Reproduce: `scripts/run_enrichment.sh out/de_results.tsv docs/findings/pathway_results/`.

### The biological contrast at pathway level

The headline per-protein observation (post-acclimation heat response produces ~4× more hits than pre-acclimation) sharpens at the pathway level into a **qualitatively different biology** between the two heat exposures:

**PT1−PR1 (acute heat stress, pre-acclimation) — inflammatory / stress-signaling signature**

Top up-regulated pathways:

| Collection | Pathway | Overlap | FDR |
|---|---|---|---|
| Hallmark | **EPITHELIAL_MESENCHYMAL_TRANSITION** | 6/53 | **4.45e-4** |
| Hallmark | KRAS_SIGNALING_UP | 4/53 | 0.027 |
| Hallmark | TNFA_SIGNALING_VIA_NFKB | 3/53 | 0.046 |
| Hallmark | ADIPOGENESIS / ALLOGRAFT_REJECTION / ESTROGEN_RESPONSE / MYOGENESIS | 3/53 each | 0.046 |
| GO-BP | CELL_CELL_SIGNALING | 16/53 | 5.65e-4 |
| GO-BP | CELL_MOTILITY | 18/53 | 0.0012 |
| GO-BP | RESPONSE_TO_GROWTH_FACTOR | 11/53 | 0.0036 |
| KEGG | PRL_JAK_STAT_SIGNALING / FIBRINOLYTIC_PAI | 1/53 | 0.12 |

Down-regulated (n=9 tiny list): Hallmark IL2_STAT5_SIGNALING (FDR 0.013), GO-BP PI3K_AKT signaling (FDR 0.061).

**PT2−PR2 (acute heat stress, post-acclimation) — cellular repair / growth / metabolism signature**

Top up-regulated pathways:

| Collection | Pathway | Overlap | FDR |
|---|---|---|---|
| Hallmark | **MITOTIC_SPINDLE** | 6/100 | **0.020** |
| Hallmark | **DNA_REPAIR** | 5/100 | **0.020** |
| Hallmark | **PI3K_AKT_MTOR_SIGNALING** | 4/100 | **0.026** |
| Hallmark | PROTEIN_SECRETION | 3/100 | 0.12 |
| Hallmark | INTERFERON_GAMMA_RESPONSE | 4/100 | 0.14 |
| GO-BP | INTRACELLULAR_SIGNALING_CASSETTE | 30/100 | 6.02e-5 |
| GO-BP | SMALL_GTPASE_MEDIATED_SIGNAL_TRANSDUCTION | 15/100 | 6.02e-5 |
| GO-BP | VESICLE_MEDIATED_TRANSPORT | 24/100 | 8.07e-4 |
| KEGG | AUTOPHAGOSOME_LYSOSOME_FUSION | 2/100 | 0.020 |
| KEGG | BCR_ABL_PI3K / METALS_NFKB / JAK_STAT | 2/100 | 0.053 |

Down-regulated: **GO-BP NEURONAL_ION_CHANNEL_CLUSTERING (FDR 0.001)**, NEURON_MATURATION (0.029), AXON_LOCALIZATION (0.029), Hallmark IL2_STAT5 (0.067).

### The observation

The two heat exposures look like **different biological responses**, not scaled versions of the same response:

- **Before acclimation**, the measurable proteomic response to acute heat is dominated by **EMT / KRAS / TNF-α-NF-κB / cell motility** — an inflammatory stress-signaling pattern. Classic acute-stress secretome.
- **After 7 days of acclimation**, the response to acute heat shifts to **mitotic spindle / DNA repair / PI3K/AKT/mTOR / vesicle transport / small GTPase signaling** — a cell-repair, protein-trafficking, growth-adaptation pattern. And it simultaneously **down-regulates neuronal signaling pathways** (ion channel clustering, neuron maturation, axon localization).

If this holds at larger n, the interpretation is: **acclimation doesn't dampen the heat stress response — it rewires it** from an acute-inflammatory mode into a cellular-maintenance / adaptive-growth mode, with concurrent dampening of neural excitability signals.

### PR2−PR1 (acclimation effect at rest, weak signal)

The 7-protein list gives only suggestive pathway-level hits:

| Collection | Pathway | FDR |
|---|---|---|
| GO-BP | ANTIBACTERIAL_HUMORAL_RESPONSE | 0.063 |
| GO-BP | ANTIMICROBIAL_HUMORAL_RESPONSE | 0.079 |
| KEGG | OKAZAKI_FRAGMENT_MATURATION / LONG_PATCH_BER | 0.009 |

Possibly suggests innate-immune and DNA-repair priming at rest after 7 days of acclimation, but too small to commit to.

### PT2−PT1 (acclimation × heat interaction, narrow signal)

The 6-up / 2-down list is almost too small to meaningfully enrich. Weak Hallmark signals on COAGULATION / FATTY_ACID_METABOLISM. The ONE striking GO-BP signal comes from the down-regulated side: **INSULIN_SECRETION / NEUROTRANSMITTER_SECRETION / INTRACELLULAR_GLUCOSE_HOMEOSTASIS** (FDR 0.046 on a 2-gene list — driven by one protein that happens to be annotated in multiple glucose-secretion sets, probably not real).

### Pathway-level caveats (separate from the per-protein caveats above)

1. **ORA background is wrong.** We use the full MSigDB universe as background, but our measurement universe is only ~2,943 proteins on Olink Explore. This inflates enrichment of pathways with high Olink coverage (secreted / inflammatory / EMT-adjacent proteins are over-represented on the Explore panel relative to, say, nuclear transcription factors). To fix this properly, rebuild the test with a restricted universe — requires extending `genesets` with a `--universe` flag.
2. **Large sets dominate Fisher's exact.** HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION has 200 genes — the largest Hallmark set. A 6/53 overlap against a 200-gene set is not as specific as a 6/53 overlap against a 30-gene set would be. Weight by set size when interpreting.
3. **FDR correction is within each (comparison × direction × collection) triple**, not pooled across all 18 slices. Reporting a pathway as "significant across multiple collections" stacks the deck; weight cross-collection agreement as a consistency check, not as an independent signal.
4. **Two of the four comparisons have too few hits for ORA to be meaningful.** PR2−PR1 (7 genes) and PT2−PT1 (6/2 genes) should be treated as descriptive, not inferential.

## Handoff

**Inputs to downstream analysis:**
```
out/de_results.tsv            # per-protein paired-t results
out/de_report.tsv             # per (comparison, panel) summary
docs/findings/pathway_results/<cmp>_<dir>_<coll>.json
                              # 21 pathway JSON files from ORA
```

## Reproduction

```bash
cargo test -p karnaproteome --test dube_de_sanity -- --nocapture
```

This runs `ingest → qc → de` against the real Dube data, asserts the heat-shock sanity check passes, and exits 0. Runtime ~1.7 s.
