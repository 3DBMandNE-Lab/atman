# Atman comparator benchmark

Head-to-head against established Olink/proteomics pipelines on the Dube
heat-acclimation dataset. Produces Table 1 and the overlap figure used in the
Bioinformatics submission.

## Comparator set

| Tool           | Version    | Role                                              |
|----------------|------------|---------------------------------------------------|
| Atman  | 1.0.0      | subject under test (`de --test paired-t`, `moderated`) |
| OlinkAnalyze   | CRAN 5.0.0 | Olink's official R package; `olink_ttest` + `olink_lmer` + `olink_pathway_enrichment` |
| limma          | BioC 3.66  | moderated-t workhorse, fit against the NPX matrix |

`DEP` and any Olink Python package were evaluated and dropped — DEP is
off-label for NPX (built for imputation-first MS intensity workflows) and
no official Olink Python DE package exists.

## Contrasts

All four Dube comparison families:
`PT1-PR1`, `PT2-PR2`, `PT2-PT1`, `PR2-PR1` with `min_pairs = 5`.

## Metrics

1. DE hit count at `q < 0.05` and `q < 0.10` per tool per contrast.
2. Pairwise Jaccard overlap of DE hit sets per contrast.
3. Wall-clock runtime (single run, no warmup discarded).
4. Peak resident-set memory (`/usr/bin/time -v` on Linux; R `pryr::mem_used()` fallback).
5. Lines-of-code to reproduce (number of non-blank lines in the tool's driver script).

## Layout

```
benchmarks/
  README.md                this file
  run_all.sh               top-level driver — runs all three tools and aggregator
  atman/run.sh     Atman baseline harness
  olinkanalyze/run.R       OlinkAnalyze comparator (olink_ttest + olink_lmer)
  limma/run.R              limma on NPX matrix comparator
  aggregate.py             per-tool outputs → Table 1 + UpSet figure
  capability_matrix.tsv    Atman vs others — capability gap table
  out/                     per-tool TSVs + Table 1 + figure (gitkept)
```

## Per-tool output schema

Each tool writes `benchmarks/out/<tool>_de.tsv` with columns:

```
contrast  gene  mean_diff  pvalue  qvalue  tool  runtime_s  peak_rss_mb
```

Aggregator reads these and emits `benchmarks/out/table1.tsv` and
`benchmarks/out/upset_de_overlap.pdf` for insertion into the manuscript.

## Reproduce

Prerequisites: `cargo` (Rust 1.75+), `Rscript`, OlinkAnalyze + limma installed,
`python3` with `pandas`, `matplotlib`, `upsetplot`. The Dube example data at
`example_data/dube_heat_2023/` must be present, and the core Dube pipeline
under `out/` must have been run (e.g. via `scripts/reproduce_paper.sh`).

```bash
bash benchmarks/run_all.sh
```
