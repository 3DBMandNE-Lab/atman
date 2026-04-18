# Supplementary — cross-platform generalizability: Gisby 2022 SomaLogic SomaScan v4.1

## Dataset

Gisby JS, Buang NB, Papadaki A, Clarke CL, Malik TH, Medjeral-Thomas N, et al. *Multi-omics identify LRRC15 as a COVID-19 severity predictor and persistent pro-thrombotic signals in convalescence.* Nat Commun. 2022;13:7775. Zenodo DOI [10.5281/zenodo.6497251](https://doi.org/10.5281/zenodo.6497251), CC BY 4.0.

- Platform: **SomaLogic SomaScan v4.1** (7,240 SOMAmer reagents → 6,360 unique protein gene symbols)
- Design: **Paired longitudinal**, same dialysis cohort as Gisby 2021 [17]
- Patients retained: **n = 28** with both an early (≤7 d from first swab) and a late (≥14 d) sample
- Units: raw SomaScan RFU → log₂-transformed by the Atman ingest adapter (log₂ of median ≈ 9.3, a typical normalized-RFU scale)
- Same biological context (COVID dialysis) as the Gisby 2021 Olink demonstration [17], so cross-platform concordance can be assessed on the same patient biology.

## Adapter

`scripts/ingest_gisby_somascan.py` (~130 lines of Python) reads `soma_abundance.csv` + `sample_technical_meta.csv` + `feature_meta.csv`, log₂-transforms the RFU matrix, maps `seq.XXXX.YY` SOMAmer identifiers to UniProt and `entrez_gene_symbol` via the feature metadata, constructs a paired late − early contrast on the 28 patients with both timepoints, and writes Atman's long-format TSVs. **No Rust engine changes are required**; the same `atman de`, `atman robustness`, `atman module-de`, and `atman stability-ranked` commands run unchanged.

## Pipeline

```bash
python3 scripts/ingest_gisby_somascan.py example_data/gisby_somascan_2022 out_gisby_soma
atman de --input-dir out_gisby_soma --output-dir out_gisby_soma \
    --test paired-t --paired-by participant \
    --groups "late-early" --min-pairs 10
python3 scripts/learn_modules.py out_gisby_soma/per_subject_gene_deltas.tsv \
    out_gisby_soma/modules_learned.tsv --k-modules 25
atman module-de --input-dir out_gisby_soma --output-dir out_gisby_soma \
    --modules-tsv out_gisby_soma/modules_learned.tsv \
    --test paired-t --paired-by participant \
    --groups "late-early" --min-pairs 10
```

or one-shot: `bash scripts/reproduce_gisby_somascan.sh`.

## Per-protein DE results

- 6,360 SOMAmers tested (paired Student's $t$, min_pairs=10)
- **226 hits at q<0.05**; **407 at q<0.10**

### Top down-regulated (late vs early) — acute-inflammation resolution

| SOMAmer (gene) | mean_diff (log₂-RFU) | q |
|---|---|---|
| C1QC          | −2.87 | 2.1×10⁻⁵ |
| CCL7          | −0.81 | 4.3×10⁻⁴ |
| IFNA1         | −1.03 | 1.2×10⁻³ |
| TNFSF13B      | −0.49 | 1.2×10⁻³ |
| LARGE1        | −0.47 | 1.2×10⁻³ |
| BTC           | −0.46 | 1.2×10⁻³ |
| CD209 (DC-SIGN) | −0.40 | 1.2×10⁻³ |
| TLR5          | −0.62 | 1.2×10⁻³ |
| TNFAIP6       | −0.54 | 1.2×10⁻³ |
| HLA-C         | −0.45 | 1.2×10⁻³ |
| PSAP          | −0.46 | 1.3×10⁻³ |
| **CXCL10**    | **−1.28** | **1.6×10⁻³** |
| CTSO          | −0.66 | 1.7×10⁻³ |
| GRN           | −0.41 | 2.1×10⁻³ |
| IL34          | −0.18 | 2.1×10⁻³ |

## Cross-platform concordance with Gisby 2021 Olink

Core inflammation markers reproduced on both platforms from the **same patient cohort** (different patients overlap; the SomaScan cohort is a superset/adjacent sample set to the Olink 2021 cohort):

Olink numbers sourced from `out_gisby/de_results.tsv` (filter `comparison="late-early"`, `skip_reason` empty). SomaScan numbers sourced from `out_gisby_soma/de_results.tsv` (same filter).

| gene | Olink Target 96 (n=23, paired) | SomaScan v4.1 (n=28, paired) |
|---|---|---|
| CXCL10   | mean_diff −1.88, q = 5.5×10⁻⁶        | mean_diff −1.28, q = 1.6×10⁻³ |
| TNFSF13B | mean_diff −0.88, q = 5.3×10⁻⁴        | mean_diff −0.49, q = 1.2×10⁻³ |
| GRN      | mean_diff −0.61, q = 1.2×10⁻⁵        | mean_diff −0.41, q = 2.1×10⁻³ |
| CCL16    | mean_diff +1.26, q = 1.6×10⁻⁴        | mean_diff +0.45, q = 0.079     |
| CCL17    | mean_diff +0.73, q = 0.023           | not significant                |
| CCL24    | mean_diff +1.60, q = 2.2×10⁻⁵        | not significant                |

All three core inflammation hits (CXCL10, TNFSF13B, GRN) reproduce on SomaScan at q<0.01 with matching direction. Effect-size magnitudes differ as expected because NPX and RFU use different internal normalisations; the sign and the rank ordering within shared proteins are concordant.

## SomaScan-only biology surfaced

SOMAmer coverage extends beyond Olink's Inflammation panel and the adapter pulls those additional proteins into atman's output unchanged. Several appear at the top of the SomaScan DE and are not tested on Olink:

- **C1QC** (q = 2.1×10⁻⁵) and the complement module below: the complement arm of acute-phase resolution.
- **IFNA1** (q = 1.2×10⁻³): type-I interferon, absent from Olink Target 96 Inflammation.
- **CD209, TLR5, HLA-C, CTSO, PSAP**: innate-immune pattern-recognition, MHC class I, lysosomal hydrolases.

This is not a biological discovery but a direct consequence of platform-agnostic ingest: Atman tests whatever proteins are present in the input matrix.

## Data-driven module-DE on SomaScan

K=25 modules learned by complete-linkage clustering on per-patient late−early delta correlations. **2 modules reject at q<0.05, 3 at q<0.10**.

| module | n_genes | mean_diff | q | notable members |
|---|---|---|---|---|
| `module_21` | 82  | −0.130 | **0.046** | **Complement + innate immune**: C1RL, C2, C4BPA, C7, CD14, CD48, CD72, CXCL12, ENTPD1, FSTL1, TNFAIP6, … |
| `module_12` | 165 | −0.100 | 0.046 | Broad acute-phase dispersion |
| `module_10` | 179 | −0.101 | 0.075 | Additional acute-response proteins |

`module_21` is the standout: 82 SOMAmer reagents, anchored by a complement-pathway core (C1RL, C2, C4BPA, C7) and early-immune components (CD14, CD48, CD72). The module emerges entirely unsupervised from the correlation structure of per-patient late−early deltas and collapses to a single testable unit at q=0.046 — the complement-arm equivalent of the Olink `module_14` interferon regulon, on the same patient biology, recovered by the same pipeline on a different platform.

## Interpretation

The SomaScan application is not a biological discovery paper. Its purpose is methodological: demonstrate that Atman's statistical machinery (paired-t, LOO, S-score, module-DE) runs unchanged on a non-Olink platform through a thin adapter, and that the small-cohort failure modes the framework addresses (independence assumption, q-only ranking, LOO fragility) are platform-general.

## Deposited artefacts

```
example_data/gisby_somascan_2022/          ← raw CSVs mirrored from Zenodo (CC BY 4.0)
scripts/ingest_gisby_somascan.py           ← ~130-line Python adapter
scripts/reproduce_gisby_somascan.sh        ← one-command end-to-end rerun
out_gisby_soma/                             ← full pipeline outputs
  measurements.tsv, qc_measurements.tsv, proteins.tsv, samples.tsv
  de_results.tsv, de_report.tsv
  per_subject_gene_deltas.tsv
  modules_learned.tsv, modules_learned_summary.tsv
  module_de_results.tsv
```
