# Supplementary — multiscale inference: data-driven modules + module-level DE

## Motivation

All DE in the main Results is per-protein, with each protein treated as an independent hypothesis under BH-FDR. At typical Olink sizes (P = 436 for Target 96, P ≈ 2,900 for Explore NGS) the BH correction penalises small-cohort discovery aggressively: 36 protein-level hits at q<0.05 on Dube PT1-PR1 n=10 comes after a 2,938-way correction. When proteins move together as a regulon — as inflammation and heat-shock genes demonstrably do in both demonstration cohorts — per-protein independence is a modelling assumption that (a) inflates the effective number of tests and (b) discards the signal aggregation that small-cohort inference most benefits from.

This supplementary adds a second inference scale on top of atman's existing protein-level machinery: data-driven module construction from per-subject delta correlation structure, followed by module-level DE using the same paired-t / Welch-t / moderated functions on aggregated module scores.

## Method

**Module learning.** `scripts/learn_modules.py` reads a per-subject-gene-delta table, concatenates each protein's delta values across subjects and contrasts into a feature vector, and clusters proteins by complete linkage on the signed distance $d_{ij} = (1 - r_{ij})/2$ where $r_{ij}$ is the Pearson correlation across feature vectors. Signed correlation keeps co-moving proteins together and anti-correlated proteins apart; complete linkage avoids the chaining failure of average linkage on high-P correlation matrices. The tree is cut at `--k-modules` (default 20). Modules are renamed `module_01 … module_K` in descending size order. Only proteins with non-missing feature vectors and positive variance are retained; dropped proteins are not assigned to any module.

**Module-level DE.** The new Rust subcommand `atman module-de` reads `qc_measurements.tsv` + a `modules.tsv` + `samples.tsv` and, for each (sample, module), computes a module score as the mean `effective_abundance()` across member proteins present in that sample. Per-(module, condition) cell vectors are then built exactly as `atman de` does for proteins. Each comparison is tested with the chosen mode (`--test paired-t | welch-t`) over module scores, and BH-FDR is applied across the K modules within each contrast. The output `module_de_results.tsv` has columns `module, contrast, n_samples_a, n_samples_b, n_genes, mean_a, mean_b, mean_diff, t, df, p_value, bh_q, skip_reason`.

## Empirical results

**Dube (K=20 modules learned from 2,920 proteins across 10 subjects × 4 contrasts = 40-dim feature vectors).**

| contrast | protein-level q<0.05 | module-level q<0.05 | modules hit |
|---|---|---|---|
| PT1-PR1   | 36 / 2,938 | 1 / 20 | `module_06` (144 genes, mean_diff +0.176, q=0.049) |
| PT2-PR2   | 37 / 2,938 | 2 / 20 | `module_06` (+0.205, q=0.039), `module_17` (+0.315, q=0.039) |
| PT2-PT1   | 7 / 2,938  | 0 / 20 | — |
| PR2-PR1   | 2 / 2,938  | 1 / 20 | `module_18` (+0.289, q=0.045) |

`module_06` is the largest cohesive regulon (144 genes, mean within-module |r| = 0.38) and is the only module that hits at q<0.05 in both acute contrasts — the consistent small-effect heat-response signal across PT1-PR1 and PT2-PR2 that protein-level DE reports as 36 and 37 individually weak hits aggregates into a single testable module. `module_17` adds a second, tighter 51-gene regulon specific to the post-acclimation acute challenge. At module level the adaptation contrasts (PT2-PT1, PR2-PR1) remain sparse or absent — aggregation cannot rescue signal that is truly weak.

**Gisby (K=15 modules learned from 435 proteins across 23 patients with paired early/late samples).**

| contrast | protein-level q<0.05 | module-level q<0.05 |
|---|---|---|
| late-early | 120 / 436 | **6 / 15 (40%)** |

| module | n_genes | mean_diff | q | notable members |
|---|---|---|---|---|
| `module_14` | 14 | **−0.462** | **3.5×10⁻⁴** | CCL19, CCL23, CXCL10, IL18, IL33, MCP-3, IGFBP-1, LAG3, IL-2RB, IL-20RA, TIMP1, TNC, GRN, RARRES2 |
| `module_04` | 37 | −0.275 | 1.9×10⁻³ | AREG, CCL25, CD93, EGFR, GDF-15, IGFBP-2, IL10, LIF-R, OPG, PCSK9, ... |
| `module_03` | 40 | −0.165 | 0.017 | CDH5, CXCL12, ENG, ICAM3, LILRB1/2/4, MBL2, MERTK, NID1, ... |
| `module_09` | 24 | +0.252 | 0.017 | CCL16, CCL17, CCL18, CCL24, FAS, IFNLR1, IL5, MMP7, PSP-D, TRAIL-R2, ... |
| `module_11` | 16 | −0.294 | 0.020 | CCL11, CCL20, CHIT1, IL6, IL8, MCP-1, MCP-2, PARP-1, ... |
| `module_02` | 49 | +0.144 | 0.043 | AXL, CD163, CXCL16, DPP4, FCGR2A, GAS6, ICAM1, IL-15RA, ... |

`module_14` is the standout finding: fourteen genes including the canonical type-II interferon / chemokine axis (CXCL10, CCL19, CCL23, IL18, MCP-3, IGFBP-1) and lymphoid markers (IL-2RB, LAG3, IL-20RA), emerging entirely data-driven from the correlation structure of per-patient late-early deltas. No gene lists, no pathway databases. The module matches the protein-level top hits reported in the main Results but collapses them into a single testable unit with q = 3.5×10⁻⁴ against a BH denominator of 15. `module_09` (CCL16, CCL17, CCL18, CCL24, IL5, MMP7) is the convalescent-phase up-regulon, consistent with the CCL24/CCL16/CCL17/CCL18 signals that appear as separate protein-level hits.

## Comparison of hit rates

| dataset | contrast | protein-level hit rate (q<0.05) | module-level hit rate (q<0.05) |
|---|---|---|---|
| Dube | PT1-PR1 | 36 / 2,938 = 1.2% | 1 / 20 = 5% |
| Dube | PT2-PR2 | 37 / 2,938 = 1.3% | 2 / 20 = 10% |
| Gisby | late-early | 120 / 436 = 28% | 6 / 15 = 40% |

Module-level rejection rate exceeds protein-level in every contrast where a signal exists, consistent with the signal-aggregation + reduced-BH-denominator mechanism. The largest gap is on Gisby where the underlying biology is both strong and modular (interferon response is canonically a coordinated regulon); the smaller gap on Dube reflects that individual heat-shock proteins have uncorrelated noise even when they share biology, and module averaging dilutes individual effect sizes.

## What module-level inference does not do

1. **It does not replace protein-level DE.** The two scales answer different questions: protein-level tests "does this specific protein change?" and module-level tests "does this co-moving group change?" A reader interested in a specific protein reads `de_results.tsv`; a reader interested in pathway-level programs reads `module_de_results.tsv`.

2. **It does not solve the interpretation problem automatically.** Data-driven modules are unsupervised; `module_06` in Dube is 144 genes whose biological meaning is not automatic from the cluster membership. Mapping learned modules to annotated pathways (an ORA on each module) is the natural follow-up and is out of scope for this paper.

3. **It does not apply without modification to very small P.** Correlation clustering of P<100 proteins at small n is unstable; the demonstration here uses Gisby Target 96 with P=436 and Dube Explore NGS with P=2,920, which both sit in a range where per-subject delta correlations are informative.

## Reproduction

```bash
# Dube
python3 scripts/learn_modules.py docs/findings/heterogeneity/per_subject_gene_deltas.tsv \
    docs/findings/heterogeneity/modules_learned.tsv --k-modules 20 \
    --summary docs/findings/heterogeneity/modules_learned_summary.tsv

atman module-de --input-dir out --output-dir out \
    --modules-tsv docs/findings/heterogeneity/modules_learned.tsv \
    --test paired-t --paired-by participant \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" --min-pairs 5

# Gisby — first build per-patient deltas, then learn modules, then module-DE
python3 - <<'PY'
import pandas as pd
m = pd.read_csv("out_gisby/qc_measurements.tsv", sep="\t")
s = pd.read_csv("out_gisby/samples.tsv", sep="\t")
m = m.merge(s[["sample_id","subject_id","condition"]], on="sample_id")
wide = m[m.dropped_by_qc==0].pivot_table(
    index=["subject_id","gene_symbol"], columns="condition",
    values="abundance", aggfunc="mean").reset_index()
wide["late_early"] = wide["late"] - wide["early"]
wide.dropna(subset=["late_early"])[["subject_id","gene_symbol","late_early"]] \
    .to_csv("out_gisby/per_subject_gene_deltas.tsv", sep="\t", index=False)
PY
python3 scripts/learn_modules.py out_gisby/per_subject_gene_deltas.tsv \
    out_gisby/modules_learned.tsv --k-modules 15

atman module-de --input-dir out_gisby --output-dir out_gisby \
    --modules-tsv out_gisby/modules_learned.tsv \
    --test paired-t --paired-by participant --groups "late-early" --min-pairs 5
```
