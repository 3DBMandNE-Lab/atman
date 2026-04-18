# Supplementary — full benchmark methodology and per-tool outputs

Complete methodology, software versions, and per-tool output tables for the
Table 1 comparison in the main text.

## Software versions

| Tool          | Version        | Source                                      |
|---------------|----------------|---------------------------------------------|
| Atman         | 1.0.0          | this repository, `crates/atman/`            |
| OlinkAnalyze  | 5.0.0          | CRAN, installed 2026-04-16                  |
| limma         | 3.66.0         | Bioconductor                                |
| R             | 4.5.0          | system                                      |
| Python        | 3.13           | miniconda                                   |
| pandas        | current        | for aggregator                              |
| upsetplot     | 0.9.0          | for UpSet figure                            |
| matplotlib    | 3.10.5         | figure rendering                            |

## Hardware

MacBook-class hardware, Darwin 25.2.0, Apple Silicon. All tools run on a
single CPU with no explicit parallelization. Runtime numbers are single-run
wall-clock; no warmup discarded, no median over repeats. Peak RSS was not
collected in this pass (instrumentation relied on `Sys.time()` and
`time.time()` rather than `/usr/bin/time -v` which is Linux-only); a future
pass on a Linux host will add peak RSS for the final Table 1 column.

## Dataset

Dube et al. 2023 *Scientific Data* — Olink Explore NGS heat-stress/acclimation
data. Raw NPX:
- `example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv`
- `example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv`

10 participants × 4 conditions (PR1, PT1, PR2, PT2) + 3 technical controls.

Contrasts tested: `PT1-PR1`, `PT2-PR2`, `PT2-PT1`, `PR2-PR1` with
`min_pairs = 5`. Each comparison family tests 2,938 proteins.

## Per-tool protocols

### Atman

- Parser: `atman ingest --platform olink-explore-ngs --parser dube`
- QC: `atman qc --rule dube`
- Differential abundance: `atman de --test paired-t --paired-by participant`
- Multiple testing: Benjamini-Hochberg within each contrast family
- Timing: wrapper in `benchmarks/atman/run.sh`, `time.time()`

### OlinkAnalyze 5.0.0

- Parser: `OlinkAnalyze::read_NPX` (handles Olink's native NPX export)
- QC: same Dube rule applied via `filter(Assay_Warning == "PASS", QC_Warning == "PASS")`
- Differential abundance: `olink_ttest(df, variable="condition", pair_id="SubjectID")`
- Multiple testing: OlinkAnalyze applies its own BH (`Adjusted_pval`)
- Secondary model: `olink_lmer(df, variable="condition", random="SubjectID")` across all four conditions, for repeated-measures sanity check. Full `olinkanalyze_lmer.tsv` attached.
- Timing: `Sys.time()` around the four-contrast loop
- Script: `benchmarks/olinkanalyze/run.R`

### limma 3.66.0

- Parser: same `OlinkAnalyze::read_NPX` to ingest Dube CSVs into a long table
- Wide matrix keyed by **OlinkID** (not gene symbol) because a single gene can appear on multiple panels (e.g., `cardiometabolic` + `cardiometabolic_ii`); gene symbols are mapped back at write time
- Per-row NA policy: rows with at least one non-NA NPX kept; limma handles per-row NAs internally. Partial-NA warnings for 71 probes (PT1-PR1) and 105 probes (PT2-PR2) are expected at this n.
- Design: `~ participant + condition` with `participant` as blocking factor — standard paired design
- Fit: `lmFit` → `eBayes` → `topTable(adjust.method="BH")`
- Timing: `Sys.time()` around the four-contrast loop
- Script: `benchmarks/limma/run.R`

## Results — hit counts

Per-tool differential abundance hits per contrast:

| tool                 | contrast | n_tested | n_q_lt_005 | n_q_lt_010 | runtime_s |
|---------------------|----------|----------|------------|------------|-----------|
| atman               | PT1-PR1  | 2938     | 36         | 62         | 0.104     |
| atman               | PT2-PR2  | 2938     | 37         | 277        | 0.104     |
| atman               | PT2-PT1  | 2938     | 7          | 8          | 0.104     |
| atman               | PR2-PR1  | 2938     | 2          | 7          | 0.104     |
| OlinkAnalyze_ttest  | PT1-PR1  | 2938     | 14         | 33         | 10.808    |
| OlinkAnalyze_ttest  | PT2-PR2  | 2938     | 6          | 44         | 10.808    |
| OlinkAnalyze_ttest  | PT2-PT1  | 2938     | 0          | 3          | 10.808    |
| OlinkAnalyze_ttest  | PR2-PR1  | 2938     | 1          | 2          | 10.808    |
| limma               | PT1-PR1  | 2938     | 11         | 15         | 0.324     |
| limma               | PT2-PR2  | 2938     | 3          | 150        | 0.324     |
| limma               | PT2-PT1  | 2938     | 0          | 5          | 0.324     |
| limma               | PR2-PR1  | 2938     | 0          | 0          | 0.324     |

(Also in `benchmark_raw/table1_counts.tsv`.)

## Results — pairwise hit-set overlap at q < 0.10

| contrast | tool_a              | tool_b              | n_a | n_b | n_shared | jaccard |
|----------|--------------------|--------------------|-----|-----|----------|---------|
| PT1-PR1  | atman              | OlinkAnalyze_ttest | 62  | 33  | 31       | 0.484   |
| PT1-PR1  | atman              | limma              | 62  | 15  | 15       | 0.242   |
| PT1-PR1  | OlinkAnalyze_ttest | limma              | 33  | 15  | 15       | 0.455   |
| PT2-PR2  | atman              | OlinkAnalyze_ttest | 277 | 44  | 44       | 0.159   |
| PT2-PR2  | atman              | limma              | 277 | 150 | 145      | 0.514   |
| PT2-PR2  | OlinkAnalyze_ttest | limma              | 44  | 150 | 38       | 0.244   |
| PT2-PT1  | atman              | OlinkAnalyze_ttest | 8   | 3   | 2        | 0.222   |
| PT2-PT1  | atman              | limma              | 8   | 5   | 4        | 0.444   |
| PT2-PT1  | OlinkAnalyze_ttest | limma              | 3   | 5   | 2        | 0.333   |
| PR2-PR1  | atman              | OlinkAnalyze_ttest | 7   | 2   | 2        | 0.286   |
| PR2-PR1  | atman              | limma              | 7   | 0   | 0        | 0.000   |
| PR2-PR1  | OlinkAnalyze_ttest | limma              | 2   | 0   | 0        | 0.000   |

(Also in `benchmark_raw/table1_overlap.tsv`.)

### Interpretation

- All three tools agree on the acute-dominant structure (PT1-PR1 and PT2-PR2 carry almost all signal).
- **Containment is the informative concordance signal.** At q<0.10 on acute contrasts, Atman's hit set contains 31/33 (94%) of OlinkAnalyze's hits on PT1-PR1, 44/44 (100%) of OlinkAnalyze's on PT2-PR2, and 145/150 (97%) of limma's on PT2-PR2. Atman additionally calls 31 hits on PT1-PR1 and 233 on PT2-PR2 beyond OlinkAnalyze's set.
- **Jaccard is mechanically depressed by set-size imbalance.** The PT2-PR2 Atman–OlinkAnalyze Jaccard of 0.16 reflects the 277-vs-44 hit-count gap (max-possible Jaccard bounded by 44/277 ≈ 0.16) rather than discordance — the 44 OlinkAnalyze hits are fully contained in Atman's 277. Read containment, not Jaccard, as the concordance check.
- **Count differences track statistical-model differences, not divergent biology.** Atman uses paired Student's t with exact df; OlinkAnalyze 5.0.0 applies a Satterthwaite correction; limma uses empirical-Bayes moderation. Under different models at a fixed q threshold on identical input, count differences are expected.

## Capability matrix (reproduced here for completeness)

See `benchmarks/capability_matrix.tsv`. Rows are capabilities; columns mark
yes / no / partial per tool with a short note.

## Raw per-tool output tables

All in `benchmark_raw/`:

| file                                | contents                                                    |
|-------------------------------------|-------------------------------------------------------------|
| `atman_de.tsv`                      | per-contrast per-protein Atman DE paired-t (11,752 rows)    |
| `olinkanalyze_de.tsv`               | per-contrast per-protein OlinkAnalyze `olink_ttest` (11,752)|
| `olinkanalyze_lmer.tsv`             | OlinkAnalyze `olink_lmer` per-protein LMM (2,938)           |
| `limma_de.tsv`                      | per-contrast per-protein limma DE (11,752)                  |
| `table1_counts.tsv`                 | wide-form per-tool per-contrast count + runtime table       |
| `table1_overlap.tsv`                | pairwise Jaccard / overlap table                            |
| `table1.tsv`                        | consolidated Table 1 (per-tool wide format)                 |
| `atman_pairedt_vs_moderated.tsv`    | within-Atman comparison of paired-t vs moderated modes on the four Dube contrasts: sign concordance (1.00 across all), Spearman rank correlation of effect sizes (ρ=1.00 across all), q<0.05 hit counts, hit overlap, and containment in both directions |

All shared schema columns: `contrast, gene, mean_diff, pvalue, qvalue, tool,
runtime_s, peak_rss_mb`.

## Reproduction

```bash
bash benchmarks/run_all.sh
```

Requires `cargo`, `Rscript` with `OlinkAnalyze` + `limma` + `dplyr` + `readr`
+ `tidyr`, `python3` with `pandas` + `matplotlib` + `upsetplot`.
