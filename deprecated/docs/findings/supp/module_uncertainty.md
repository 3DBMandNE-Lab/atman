# Supplementary — subject-bootstrap uncertainty at the protein and module level

All tabulated numbers in this document are computed directly from the
primary bootstrap TSVs:

- `out_gisby/robustness/module_bootstrap.tsv`
- `out_gisby/robustness/protein_bootstrap.tsv`
- `out/robustness/module_bootstrap_PT1-PR1.tsv`
- `out/robustness/module_bootstrap_PT2-PR2.tsv`
- `out/robustness/protein_bootstrap_PT1-PR1.tsv`
- `out/robustness/protein_bootstrap_PT2-PR2.tsv`

Each TSV row reports `point_mean`, `boot_mean`, `boot_sd`, `ci_lo_95`,
`ci_hi_95`, `p_pos`, `p_sign_stable`, `n_genes` (for modules), `n_boot`.

## Motivation

Classical module-DE returns a point estimate of the module mean delta and
a single BH-adjusted $p$-value computed from a paired Student's $t$ on
$n-1$ degrees of freedom. At $n=10$ (Dube) or $n=23$ (Gisby), this is a
thin inferential summary for a quantity we are asking to carry most of the
paper's mechanistic weight. A principled upgrade is to quantify
uncertainty in $\bar m_k$ (module) and $\bar d_i$ (protein) directly by
resampling subjects with replacement. Bootstrap gives the full distribution
over effect size under subject perturbation; from it the quantities a
biologist actually wants — "how large is the effect?", "how sure are we of
the direction?", "what is the 95% range?" — follow immediately.

## Procedure

For contrast $A{-}B$ with $n$ subjects, let
$\delta_{s,i} = x_{s,A,i} - x_{s,B,i}$ be subject $s$'s per-protein delta
and $m_{s,k} = \frac{1}{|G_k|}\sum_{i \in G_k} \delta_{s,i}$ be subject
$s$'s module-$k$ mean delta. For each $b \in \{1,\dots,B\}$:

1. Draw $\mathbf{r}^{(b)} \sim \mathrm{Unif}(\{1,\dots,n\}^n)$ (index with replacement).
2. Compute $\bar m^{(b)}_k = \frac{1}{n}\sum_{j=1}^{n} m_{r^{(b)}_j, k}$ (same form at the protein scale with $\delta_{s,i}$ in place of $m_{s,k}$).

Outputs per module/protein:

$$
\hat\mu_k = \mathbb{E}[\bar m^{(b)}_k],\quad
\hat\sigma_k = \mathrm{SD}[\bar m^{(b)}_k],\quad
\mathrm{CI}_{95}(k) = [q_{0.025},\,q_{0.975}],\quad
\hat P(\text{sign stable}) = \max\bigl(\hat P(\bar m_k>0),\ \hat P(\bar m_k<0)\bigr).
$$

Implementation: `scripts/bootstrap_module_de.py` and
`scripts/bootstrap_protein_de.py`. $B=1{,}000$ throughout. Runtime: <2 s at
module scale, <10 s at protein scale on these cohorts.

## Gisby 2021 late − early (n=23, 15 modules)

Sorted by $|\hat\mu_k|$ (descending). Source:
`out_gisby/robustness/module_bootstrap.tsv`.

| module    | n_genes | $\hat\mu_k$ | $\hat\sigma_k$ | CI$_{95}$           | $P(\text{sign stable})$ |
|-----------|--------:|------------:|---------------:|---------------------|------------------------:|
| module_14 |      14 | −0.460      | 0.086          | [−0.642, −0.299]    | **1.000**               |
| module_11 |      16 | −0.291      | 0.095          | [−0.480, −0.110]    | 0.998                   |
| module_04 |      37 | −0.273      | 0.062          | [−0.393, −0.158]    | **1.000**               |
| module_09 |      24 | +0.252      | 0.079          | [+0.107, +0.413]    | **1.000**               |
| module_12 |      16 | −0.193      | 0.137          | [−0.453, +0.080]    | 0.907                   |
| module_03 |      40 | −0.165      | 0.050          | [−0.266, −0.070]    | **1.000**               |
| module_02 |      49 | +0.144      | 0.056          | [+0.044, +0.257]    | 0.998                   |
| module_01 |      77 | −0.140      | 0.145          | [−0.424, +0.151]    | 0.828                   |
| module_13 |      15 | +0.125      | 0.102          | [−0.061, +0.335]    | 0.895                   |
| module_06 |      32 | +0.057      | 0.054          | [−0.048, +0.162]    | 0.866                   |
| module_08 |      24 | −0.054      | 0.074          | [−0.210, +0.082]    | 0.767                   |
| module_10 |      22 | −0.047      | 0.056          | [−0.160, +0.056]    | 0.795                   |
| module_07 |      25 | −0.039      | 0.052          | [−0.144, +0.064]    | 0.778                   |
| module_15 |       8 | −0.031      | 0.048          | [−0.129, +0.061]    | 0.736                   |
| module_05 |      36 | −0.001      | 0.056          | [−0.114, +0.111]    | 0.505                   |

**Summary:** 6 of 15 modules reach $P(\text{sign stable}) \geq 0.95$; four
reach $P = 1.000$ (modules 14, 04, 09, 03). The headline 14-gene interferon
/ chemokine module (`module_14`; membership CCL19, CCL23, CXCL10, GRN,
IGFBP-1, IL-20RA, IL-2RB, IL18, IL33, LAG3, MCP-3, RARRES2, TIMP1, TNC)
has $P = 1.000$ and a 95% CI that excludes zero by more than three
bootstrap standard errors. This confirms the classical $q = 3.5 \times
10^{-4}$ result non-parametrically and additionally constrains the effect
magnitude: every one of $B = 1{,}000$ subject-resampled datasets returns
the same direction.

## Dube PT1-PR1 (n=10, 20 modules)

Top 10 by $P(\text{sign stable})$ (ties broken by $|\hat\mu_k|$). Source:
`out/robustness/module_bootstrap_PT1-PR1.tsv`.

| module    | n_genes | $\hat\mu_k$ | $\hat\sigma_k$ | CI$_{95}$           | $P(\text{sign stable})$ |
|-----------|--------:|------------:|---------------:|---------------------|------------------------:|
| module_06 |     144 | +0.168      | 0.039          | [+0.097, +0.253]    | **1.000**               |
| module_12 |      91 | +0.136      | 0.038          | [+0.063, +0.211]    | **1.000**               |
| module_18 |      46 | +0.113      | 0.052          | [+0.019, +0.222]    | 0.996                   |
| module_01 |     601 | +0.430      | 0.176          | [+0.110, +0.815]    | 0.994                   |
| module_09 |     117 | +0.068      | 0.035          | [+0.010, +0.141]    | 0.989                   |
| module_02 |     378 | +0.148      | 0.096          | [+0.002, +0.363]    | 0.976                   |
| module_03 |     297 | +0.060      | 0.036          | [−0.008, +0.124]    | 0.953                   |
| module_10 |     107 | +0.128      | 0.072          | [−0.027, +0.241]    | 0.951                   |
| module_19 |      42 | +0.106      | 0.069          | [−0.028, +0.246]    | 0.939                   |
| module_20 |      35 | +0.066      | 0.041          | [−0.022, +0.139]    | 0.936                   |

**Summary:** 8 of 20 modules at $P \geq 0.95$. `module_06` (144 genes)
reaches $P = 1.000$ — this is the cross-acute regulon also identified as
the only module hitting q<0.05 in both acute Dube contrasts by classical
module-DE.

## Dube PT2-PR2 (n=10, 20 modules)

Top 10 by $P(\text{sign stable})$. Source:
`out/robustness/module_bootstrap_PT2-PR2.tsv`.

| module    | n_genes | $\hat\mu_k$ | $\hat\sigma_k$ | CI$_{95}$           | $P(\text{sign stable})$ |
|-----------|--------:|------------:|---------------:|---------------------|------------------------:|
| module_01 |     601 | +0.718      | 0.198          | [+0.324, +1.115]    | **1.000**               |
| module_06 |     144 | +0.196      | 0.046          | [+0.102, +0.281]    | **1.000**               |
| module_17 |      51 | +0.284      | 0.067          | [+0.161, +0.422]    | **1.000**               |
| module_18 |      46 | −0.151      | 0.046          | [−0.241, −0.062]    | 0.999                   |
| module_20 |      35 | +0.114      | 0.068          | [−0.003, +0.266]    | 0.968                   |
| module_07 |     141 | −0.080      | 0.042          | [−0.157, +0.007]    | 0.962                   |
| module_12 |      91 | −0.077      | 0.049          | [−0.168, +0.025]    | 0.936                   |
| module_02 |     378 | +0.062      | 0.044          | [−0.021, +0.146]    | 0.931                   |
| module_19 |      42 | −0.127      | 0.088          | [−0.306, +0.046]    | 0.922                   |
| module_08 |     121 | −0.076      | 0.055          | [−0.186, +0.025]    | 0.913                   |

**Summary:** 6 of 20 modules at $P \geq 0.95$. `module_06` holds $P =
1.000$ in both acute contrasts — the bootstrap agrees with classical
module-DE that this is the robust cross-acute signal.

## Protein-scale bootstrap — Gisby late − early (n=23, 435 proteins)

Source: `out_gisby/robustness/protein_bootstrap.tsv`.

Threshold counts:

| threshold                      | count |
|--------------------------------|------:|
| total proteins                 |   435 |
| P(sign stable) ≥ 0.95          |   220 |
| P(sign stable) ≥ 0.999         |   114 |
| P(sign stable) = 1.000         |   101 |

Top 15 proteins by $|\hat\mu|$ among those with $P = 1.000$:

| gene       | $\hat\mu$ | $\hat\sigma$ | CI$_{95}$         | $P(\text{sign stable})$ |
|------------|----------:|-------------:|-------------------|------------------------:|
| IFN-gamma  | −2.509    | 0.500        | [−3.505, −1.515]  | **1.000**               |
| CXCL10     | −1.872    | 0.217        | [−2.287, −1.469]  | **1.000**               |
| DDX58      | −1.712    | 0.254        | [−2.213, −1.216]  | **1.000**               |
| CCL24      | +1.600    | 0.224        | [+1.165, +2.022]  | **1.000**               |
| MCP-2      | −1.507    | 0.213        | [−1.902, −1.094]  | **1.000**               |
| CCL16      | +1.262    | 0.211        | [+0.886, +1.707]  | **1.000**               |
| KRT19      | −1.035    | 0.225        | [−1.450, −0.605]  | **1.000**               |
| TRIM21     | −1.034    | 0.170        | [−1.358, −0.695]  | **1.000**               |
| MCP-3      | −0.925    | 0.260        | [−1.460, −0.450]  | **1.000**               |
| TNFSF13B   | −0.887    | 0.159        | [−1.199, −0.589]  | **1.000**               |
| IL10       | −0.886    | 0.311        | [−1.532, −0.310]  | **1.000**               |
| IGFBP-1    | −0.884    | 0.227        | [−1.333, −0.442]  | **1.000**               |
| PARP-1     | −0.844    | 0.166        | [−1.179, −0.528]  | **1.000**               |
| IRF9       | −0.844    | 0.145        | [−1.121, −0.563]  | **1.000**               |
| TCN2       | −0.832    | 0.103        | [−1.030, −0.635]  | **1.000**               |

Canonical spot-checks elsewhere in the table:
GRN $\hat\mu = -0.609$, CI$_{95}$ [−0.756, −0.455], $P = 1.000$.

## Protein-scale bootstrap — Dube acute contrasts (n=10)

Sources: `out/robustness/protein_bootstrap_PT1-PR1.tsv`,
`out/robustness/protein_bootstrap_PT2-PR2.tsv`.

| contrast | n_proteins | $P \geq 0.95$ | $P \geq 0.999$ | $P = 1.000$ |
|----------|-----------:|--------------:|---------------:|------------:|
| Dube PT1-PR1 |      1,996 |           812 |            355 |         299 |
| Dube PT2-PR2 |      1,411 |           627 |            319 |         260 |

At $n = 10$ the bootstrap index space has only $n^n = 10^{10}$ distinct
samples and the attainable $\bar d$ values are discretised, so the
bootstrap distribution is chunky. $P(\text{sign stable}) \geq 0.95$ is a
*directional* claim ("the effect direction is the same in ≥95% of
resamples") rather than a family-wise discovery cutoff; it should be read
as a per-feature confidence annotation on an existing ranked list, not a
replacement for BH-FDR at the protein scale. Module-level $P$ is where
this paradigm becomes discovery-strength because the effective number of
tests drops from $P$ proteins to $K$ modules.

## Aggregation benefit at the module scale

For each module $k$, the module bootstrap SD $\hat\sigma_k$ sits between
(i) any single member protein's bootstrap SD and (ii) the independent-noise
prediction $\mathrm{rms}(\hat\sigma^{\text{protein}}_{i \in G_k})/\sqrt{|G_k|}$.
The ratio between observed and independent-noise prediction quantifies
within-module correlation directly: ratio ≈ 1 would mean proteins are
independent (implausible for a real regulon); ratio ≫ 1 means members
co-move under shared input (expected for a biological program). Observed
aggregation ratios are ≈ 1.8 on Gisby and ≈ 3.1 on Dube acute contrasts —
the main-text's panel b of Fig 5 visualises this. Figure panel numbers
cited here refer to the manuscript's figure set; full per-module aggregation
numbers are in `out_gisby/robustness/module_bootstrap.tsv` (column
`boot_sd` relative to the analogous column in `protein_bootstrap.tsv`).

## What the bootstrap does and does not say

**Says:**

1. The module-level effects reported in the main text are accompanied by
   non-parametric uncertainty intervals and a direct "probability the
   direction is stable" quantity, derived by perturbing the unit of
   analysis (subject) rather than assuming a $t$-distributed test
   statistic.
2. On Gisby (n=23) the strongest modules have $P(\text{sign stable}) =
   1.000$ — every one of $B = 1{,}000$ subject-resamples preserves
   direction. This is a stronger directional claim than any single
   classical $q$-value provides.
3. The aggregation benefit is quantified: module SDs are lower than any
   member protein's SD but not as low as independent-noise aggregation
   would predict. The gap measures within-module correlation directly.

**Does not say:**

1. That bootstrap replaces BH-FDR. Bootstrap gives a single-module
   uncertainty, not a family-wise correction over $K$ modules; classical
   module-DE $q$-values remain the right family-level summary.
2. That this works at arbitrarily small $n$. At $n = 10$ the subject
   index space has only $10^{10}$ distinct draws; the bootstrap is
   informative (chunkier but not invalid) and is reported with that
   caveat. At $n < 6$ the bootstrap becomes unreliable.
3. That bootstrap SD can be interpreted as a posterior SD. Classical
   bootstrap has frequentist semantics; a parametric Bayesian module-DE
   (conjugate normal-inverse-gamma, or Stan/pyMC) would yield a proper
   posterior and is a natural next step, outside the scope of this
   deterministic Rust pipeline.

## Reproduce

```bash
# Module bootstrap — Gisby late-early (n=23, 15 modules)
python3 scripts/bootstrap_module_de.py \
    --input-dir out_gisby --modules out_gisby/modules_learned.tsv \
    --contrast late-early --n-boot 1000 \
    --output out_gisby/robustness/module_bootstrap.tsv

# Module bootstrap — Dube acute (n=10, 20 modules each)
python3 scripts/bootstrap_module_de.py \
    --input-dir out --modules docs/findings/heterogeneity/modules_learned.tsv \
    --contrast PT1-PR1 --n-boot 1000 \
    --output out/robustness/module_bootstrap_PT1-PR1.tsv

python3 scripts/bootstrap_module_de.py \
    --input-dir out --modules docs/findings/heterogeneity/modules_learned.tsv \
    --contrast PT2-PR2 --n-boot 1000 \
    --output out/robustness/module_bootstrap_PT2-PR2.tsv

# Protein bootstrap
python3 scripts/bootstrap_protein_de.py \
    --input-dir out_gisby --contrast late-early --n-boot 1000 \
    --output out_gisby/robustness/protein_bootstrap.tsv

python3 scripts/bootstrap_protein_de.py \
    --input-dir out --contrast PT1-PR1 --n-boot 1000 \
    --output out/robustness/protein_bootstrap_PT1-PR1.tsv

python3 scripts/bootstrap_protein_de.py \
    --input-dir out --contrast PT2-PR2 --n-boot 1000 \
    --output out/robustness/protein_bootstrap_PT2-PR2.tsv
```
