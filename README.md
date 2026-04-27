# Atman

[![CI](https://github.com/kevinj24fr/atman/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/kevinj24fr/atman/actions/workflows/ci.yml)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](./LICENSE-MIT)

Cross-platform proteomics CLI. One tool for Olink Explore, SomaScan,
MaxQuant/LFQ, DIA-NN, and Spectronaut — same commands, same output
schema, same reproducibility guarantees regardless of assay source.
Ships differential abundance (paired, Welch, OLS, mixed, limma,
DEqMS, msqrob, ensemble), decomposition (ICA, VCA+FCLS unmixing,
variance partition), module testing, enrichment, bootstrap intervals,
null calibration, and meta-analysis. Every invocation writes a
SHA-256 provenance sidecar so the run is auditable end-to-end.

Each statistical path is validated against R reference implementations
on published fixtures:

| Test path | R reference | Fixture | Agreement |
|---|---|---|---|
| `de --test limma --peptide-metadata` (DEqMS trend) | `DEqMS::spectraCounteBayes` | CPTAC Study 6 UPS1 spike-in | median Δ = **0.000 log₂** (three-decimal match) |
| `de --test msqrob` | `msqrob2::msqrob(~condition)` | CPTAC Study 6 UPS1 | **100% sign agreement** on all 11 signal proteins; 9/9 UPS1 spike-ins recover expected direction ([details](docs/reference.md#msqrob2-parity-note)) |
| `de --test limma` (parametric eBayes) | `limma::eBayes` | 100-feature × 20-sample regression | Δ < **1e-4** on `t`, `p_value`, `df_total`, `s2_post` |
| `de --post-hoc sidak --test ols` | `lm()` + `pairwise.t.test` | Planted 3-level fixture | max Δ < **1e-6** on `estimate`, `posthoc_p`, `posthoc_adj_p` |
| `decompose variance --omnibus-factor` (Type III F) | `car::Anova(type = 3)` | Planted variance fixture | max Δ p < **1e-6** |
| Olink Explore NPX reproduction | Dube et al. 2023 published tables | Dube heat-stress cohort | filtered NPX: **byte-exact**; log2-FC: Δ = 1.05e-15 |

Parity assertions run on every `cargo test --workspace --release`; CI
fails if any drifts. See [docs/reference.md](docs/reference.md) for
validation details, reference scripts, and the msqrob2 estimator note.

## Supported Platforms

All platforms land on the same canonical TSV schema. Use `ingest-matrix`
for stock wide-format exports, or run a Python adapter from `adapters/`
for vendor-specific long formats:

| Platform | Ingest path | Notes |
|---|---|---|
| Olink Explore NGS | `python3 adapters/generic/olink_explore_to_atman.py` | Long-format NPX CSV adapter |
| DIA-NN | `ingest-matrix --platform diann_report` | Protein-group matrix |
| MaxQuant/LFQ | `ingest-matrix --platform maxquant_lfq` | proteinGroups.txt or LFQ matrix |
| Spectronaut | `ingest-matrix --platform spectronaut_report` | Protein-group quantity table |
| SomaScan | `ingest-matrix --platform somascan` | RFU or log2-RFU matrix |
| Any wide matrix | `ingest-matrix` | Any log2 protein×sample table with metadata |

## Quick Start: MS Matrix (DIA-NN, MaxQuant, SomaScan, ...)

Three commands from a protein matrix to differential abundance results:

```bash
# 1. Ingest: wide matrix + sample metadata → canonical TSVs
atman ingest-matrix \
    --matrix protein_matrix.tsv \
    --samples sample_metadata.tsv \
    --orientation proteins-rows \
    --platform diann_report \
    --abundance-unit log2_diann_pg_quantity \
    --assay-id-col Protein.Group --gene-col Genes \
    --condition-col diagnosis \
    --normalize median \
    --output-dir out

# 2. Validate: check schema, keys, sample support
atman validate --input-dir out --groups "Case-Control" --min-pairs 5

# 3. DE: covariate-adjusted OLS with limma eBayes shrinkage
atman de \
    --input-dir out --output-dir out \
    --test limma --groups "Case-Control" \
    --design "~ condition + age + sex" \
    --contrast conditionCase \
    --min-pairs 5
```

Outputs: `de_results.tsv`, `de_report.tsv`, `de_results.tsv.run.json`
(SHA-256 provenance sidecar). See [docs/recipes.md](docs/recipes.md)
for OLS formulas, mixed models, DEqMS, msqrob, TREAT, post-hoc
corrections, ensemble consensus, bootstrap, null calibration, and more.

## Quick Start: Olink Explore

For the bundled Dube et al. 2023 fixture:

```bash
python3 adapters/generic/olink_explore_to_atman.py \
    --output-dir out \
    example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
    example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv

atman validate --input-dir out --groups "PT1-PR1,PR2-PR1" --min-pairs 5
atman qc --input-dir out --output-dir out

atman de \
    --input-dir out --output-dir out \
    --test paired-t --paired-by participant \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" \
    --min-pairs 5
```

See [docs/tutorial.md](docs/tutorial.md) for the full Dube walkthrough.

## Install

From source with stable Rust 1.94 or newer:

```bash
git clone https://github.com/kevinj24fr/atman.git
cd atman
cargo install --path crates/atman
```

Or build a container:

```bash
docker build -t atman:1.0.0 .
docker run --rm atman:1.0.0 --help
```

## Commands

```text
atman ingest-matrix      wide protein matrix + metadata → canonical TSVs
atman validate           check canonical TSV schema, keys, and sample support
atman qc                 apply QC masking rules
atman report qc          summarize QC, missingness, and condition support
atman matrix             canonical long TSV → per-panel wide NPX CSVs
atman fold-change        compute per-panel log2 fold-change tables
atman de                 paired, Welch, OLS, mixed, limma, msqrob, or ensemble DE
atman bootstrap protein  subject-level bootstrap intervals for protein effects
atman bootstrap module   subject-level bootstrap intervals for module effects
atman null               permutation/sign-flip null calibration for DE effects
atman asymmetry          compare matched contrast pairs
atman robustness         summarize rerun/LOO rank and sign stability
atman module-trajectory  score user-defined modules from per-subject deltas
atman module-de          aggregate proteins into modules and test at module level
atman score modules      score modules per sample from canonical measurements
atman enrich ora         over-representation analysis from DE hits
atman enrich gprofiler   live g:Profiler REST wrapper with cached responses
atman meta               combine DE results across cohorts
atman network influence  feature-covariance hub scoring
atman decompose ica      multi-seed FastICA with seed-stability reporting
atman decompose null     permutation-null calibration of archetype stability
atman decompose variance mixed-model variance partition per archetype
atman decompose unmix    VCA + FCLS compartmental unmixing
atman align programs     cross-cohort program alignment + sensitivity sweep
atman align bootstrap    subject-level bootstrap of alignment
atman run                execute a plan YAML/JSON and emit a hash manifest
```

**Data flow**

```
Input                  Atman stage                 Output
-------------------------------------------------------------------
raw NPX / matrix  ──▶  ingest / ingest-matrix  ──▶ canonical TSV
canonical TSV     ──▶  qc / validate / report  ──▶ QC'd TSV + report
QC'd TSV          ──▶  de / decompose / bench  ──▶ results.tsv
any stage         ──▶  (every command)         ──▶ <output>.run.json
                                                   (SHA-256 manifest)
```

## Layout

```text
crates/atman-core/       core data model and algorithms
crates/atman/            CLI, command orchestration, and file IO
adapters/                canonical TSV adapter helpers and templates
adapters/examples/       tiny synthetic matrix examples (DIA-NN, MaxQuant, SomaScan)
CITATION.cff             citation metadata for release archives
docs/tutorial.md         package tutorial using the bundled Dube fixture
docs/recipes.md          DE cookbook: every test path with examples
docs/reference.md        data model, sidecar schema, validation details, non-goals
docs/analytical-roadmap.md  implemented analytical capability roadmap
docs/release-checklist.md   standalone release checklist
example_data/            Dube et al. 2023 Olink Explore fixture data
```

## Further Reading

- **[Tutorial](docs/tutorial.md)** — end-to-end Dube fixture walkthrough
- **[DE Recipes](docs/recipes.md)** — OLS, mixed, limma, DEqMS, msqrob, TREAT, post-hoc, ensemble, bootstrap, null, modules, ORA, meta-analysis, decomposition, and `atman run`
- **[Reference](docs/reference.md)** — data model, sidecar schema, adapter responsibilities, validation details, non-goals

## Reference Dataset

Gagnon D, Barry H, Barhdadi A, Oussaid E, Mongrain I, Lemieux Perreault LP,
Dubé MP. *A dataset of proteomic changes during human heat stress and heat
acclimation.* Scientific Data (2023).
[10.1038/s41597-023-02809-5](https://doi.org/10.1038/s41597-023-02809-5).

## License

Atman is licensed under either [MIT](./LICENSE-MIT) or
[Apache-2.0](./LICENSE-APACHE), at your option.
