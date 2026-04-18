# Supplementary — stability-weighted ranking (`S`-score)

All numbers in this document are recomputed from the primary TSVs
`out/de_results.tsv`, `out/robustness/loo_sign_stability.tsv`,
`out_gisby/de_results.tsv`, `out_gisby/robustness/loo_sign_stability.tsv`,
and the per-LOO rerun directories `out_dube_loo/` and `out_gisby_loo/`.
The verification script is `scripts/figures/…` — no wait, that path is
wrong; the recomputation used `/tmp/atman_truth/compute_jaccard_by_rule.py`
which reads the TSVs listed above.

## Definition

For each protein `i` in contrast `c`, Atman's `robustness` stage emits

```
S_{i,c} = |mean_diff_i| · SSR_{i,c} · (1 − bh_q_{i,c})
```

where `SSR_{i,c}` is the per-protein LOO sign-match rate from
`out/robustness/loo_sign_stability.tsv` (column `sign_match_rate`). Under
the multiplicative interpretation:

- `|mean_diff_i|` is the raw effect magnitude.
- `SSR_{i,c} ∈ [0, 1]` is the empirical probability the sign of the effect
  is preserved under independent subject resampling.
- `(1 − bh_q_{i,c}) ∈ [0, 1]` is a monotonic proxy for significance
  reliability at protein `i`'s BH rank.

Each term zeroes the score at an edge condition (sign-unstable, or
q=1), and the product lives in effect-size units. The ranking column
`rank_by_stability` in `out/robustness/stability_ranked.tsv` is proteins
sorted by `S` within contrast.

## Top-K Jaccard stability by ranking rule

Baseline top-K sets (K = count of baseline q<0.05 hits) were ranked by
three rules and compared to the top-K sets produced by each LOO rerun
ranked under the same rule. `S` uses baseline SSR. Values are mean ± SD
of Jaccard across LOO reruns.

| dataset          | contrast    | K   | rank by `q`     | rank by `|d̄|`  | rank by `S`     |
|------------------|-------------|-----|-----------------|-----------------|-----------------|
| Dube             | PT1-PR1     |  36 | 0.737 ± 0.110   | 0.721 ± 0.113   | 0.691 ± 0.096   |
| Dube             | PT2-PR2     |  37 | 0.646 ± 0.161   | **0.841 ± 0.092** | 0.808 ± 0.117 |
| Dube             | PT2-PT1     |   7 | **0.761 ± 0.149** | 0.736 ± 0.123 | 0.519 ± 0.155   |
| Dube             | PR2-PR1     |   2 | **0.733 ± 0.344** | 0.667 ± 0.351 | 0.300 ± 0.105   |
| Gisby            | late-early  | 120 | 0.911 ± 0.029   | 0.904 ± 0.026   | 0.911 ± 0.021   |

Bold = largest mean Jaccard in that row.

### What the numbers say

1. **Mean Jaccard.** The mean top-K Jaccard is higher for `q`-ranking in
   three of the five contrast/cohort rows (PT1-PR1, PT2-PT1, PR2-PR1), higher
   for `|d̄|`-ranking in one (PT2-PR2), and effectively tied between `q`
   and `S` on Gisby (0.9111 vs 0.9109). `S`-ranking is never the strict
   mean winner. The manuscript's headline claim about `S` is therefore
   not "`S` improves mean Jaccard"; it is about variance — below.

2. **SD (variance).** `S`-ranking has **lower SD than `q`-ranking in 3 of
   4 Dube contrasts and in Gisby** (PT1-PR1 0.096 vs 0.110; PT2-PR2 0.117
   vs 0.161; Gisby 0.021 vs 0.029). The single exception is PT2-PT1 where
   `S` has slightly higher SD (0.155 vs 0.149) at K=7, a small-K corner
   where any effect-size-driven selection is volatile.

3. **Small-K corner.** Dube PT2-PT1 (K=7) and PR2-PR1 (K=2) are the
   regimes where `S` breaks down most clearly — at K=2, the baseline top-2
   vs LOO top-2 Jaccard is 0.300 under `S` versus 0.733 under `q`. Below
   ~10 baseline hits, any ranking that couples `|d̄|` and stability is
   dominated by single-protein volatility.

4. **Dominance by effect size within the standard hit window.** Every
   q<0.20 hit in both Dube acute contrasts and Gisby has SSR=1.00, so
   within that window the SSR factor does not bite and `S` reduces to
   `|d̄| · (1 − q)`. The practical divergence between `rank_by_q` and
   `rank_by_stability` on hit lists is therefore driven by `|d̄| · (1 − q)`
   relative to `q` alone, which is small when q-values span many orders of
   magnitude (Gisby n=23) and larger when q-values cluster near the BH
   threshold (Dube acute contrasts at n=10).

### Bottom line

`S` is not a discovery-optimal ranking. The simulation grid (see
`s_score_simulation.md`) shows `|d̄|` is F1-optimal in 258/270 cells and
`S` is systematically between `q` and `|d̄|`. `S`'s empirical virtue is
**reduced ranking variance across LOO reruns in 3 of 4 Dube contrasts and
on Gisby**, at the cost of mean Jaccard in the same contrasts. Under the
bootstrap paradigm introduced in the main text, the directional
probability `P(sign stable)` supersedes both `SSR` (its LOO analogue) and
`S` (the combined ranking) as the primary per-feature directional
quantity; `S` is retained for deployments without bootstrap compute.

## Reproducing the numbers

```
python3 /tmp/atman_truth/compute_jaccard_by_rule.py
```

reads `out/de_results.tsv`, `out_gisby/de_results.tsv`, and the LOO rerun
`de_results.tsv` files under `out_dube_loo/` and `out_gisby_loo/` and
writes `jaccard_by_rule.tsv` with the numbers above. The script is
deterministic: running it twice on the same TSVs produces the same output.

## Source files

| artefact                                         | path                                                         |
|--------------------------------------------------|--------------------------------------------------------------|
| Dube baseline DE                                 | `out/de_results.tsv`                                         |
| Dube LOO reruns (one subdir per dropped subject) | `out_dube_loo/loo_{001B,007,008,009,011,014,015,016,018,019}/de_results.tsv` |
| Dube LOO sign-match                              | `out/robustness/loo_sign_stability.tsv`                      |
| Dube S-ranked table                              | `out/robustness/stability_ranked.tsv`                        |
| Gisby baseline DE                                | `out_gisby/de_results.tsv`                                   |
| Gisby LOO reruns                                 | `out_gisby_loo/loo_{…23 patient IDs…}/de_results.tsv`         |
| Gisby LOO sign-match                             | `out_gisby/robustness/loo_sign_stability.tsv`                |
| Gisby S-ranked table                             | `out_gisby/robustness/stability_ranked.tsv`                  |
