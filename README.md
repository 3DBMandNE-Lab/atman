# Atman

[![CI](https://github.com/kevinjoseph/atman/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/kevinjoseph/atman/actions/workflows/ci.yml)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](./LICENSE-MIT)

Atman is a standalone Rust command-line tool for reproducible proteomics
workflows. It has a built-in Olink Explore NPX ingest path, a simple canonical
TSV interchange format for other proteomics sources, and analysis commands for
differential abundance, module testing, asymmetry, and robustness summaries.

Atman does not require a service, database, notebook runtime, or workflow
system. The repository includes a Dube et al. 2023 Olink Explore fixture
dataset so the package can test its raw-to-result reproduction path end to end.

## Install

From source with stable Rust 1.75 or newer:

```bash
git clone https://github.com/kevinjoseph/atman.git
cd atman
cargo install --path crates/atman
```

Or build a container:

```bash
docker build -t atman:1.0.0 .
docker run --rm atman:1.0.0 --help
```

## Input Support

Atman has two input modes:

- **Built-in ingest:** `atman ingest --platform olink-explore-ngs` reads raw
  long-format Olink Explore NPX CSV files and writes Atman's canonical TSVs.
- **Canonical TSV adapters:** any script or converter can write
  `samples.tsv`, `proteins.tsv`, `measurements.tsv`, and optionally
  `qc_measurements.tsv`. Once those files exist, Atman's downstream commands
  (`de`, `module-de`, `robustness`, and related summaries) run the same way
  regardless of the original assay source.

The canonical TSV route is how Atman can be used with other proteomics
matrices, including log2 LFQ/intensity exports from MaxQuant-like workflows,
Spectronaut protein-group quantity tables, DIA-NN protein-group matrices,
SomaScan-style log2 abundance matrices, and already-normalized CSV/TSV
protein matrices. Those adapters are intentionally thin: normalize source
metadata, map samples and proteins, log-transform linear intensities when
needed, and emit Atman's TSV schema.

## Commands

```text
atman ingest             raw Olink Explore long CSV -> canonical TSVs
atman validate           check canonical TSV schema, keys, and sample support
atman qc                 apply QC masking rules
atman report qc          summarize QC, missingness, and condition support
atman matrix             canonical long TSV -> per-panel wide NPX CSVs
atman fold-change        compute per-panel log2 fold-change tables
atman de                 paired, moderated, Welch, or OLS differential abundance
atman asymmetry          compare matched contrast pairs
atman robustness         summarize rerun/LOO rank and sign stability
atman module-trajectory  score user-defined modules from per-subject deltas
atman module-de          aggregate proteins into modules and test at module level
```

## Quick Start

```bash
atman ingest \
    --platform olink-explore-ngs --parser dube \
    --output-dir out \
    example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
    example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv

atman qc --input-dir out --output-dir out --rule dube

atman validate \
    --input-dir out \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" \
    --min-pairs 5

atman report qc \
    --input-dir out \
    --output-dir out/report

atman matrix \
    --input-dir out --output-dir out \
    --format dube-wide --split-by panel

atman fold-change \
    --input-dir out --output-dir out \
    --groups "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1"

atman de \
    --input-dir out --output-dir out \
    --test paired-t --paired-by participant \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" \
    --min-pairs 5
```

Common outputs:

- `measurements.tsv`, `qc_measurements.tsv`, `samples.tsv`, `proteins.tsv`
- optional `validate_report.tsv`
- `report/qc_summary.tsv`, `report/sample_qc.tsv`, `report/protein_qc.tsv`,
  `report/condition_counts.tsv`
- `<panel>_npx.csv`
- `<panel>_log2_fc.csv`
- `de_results.tsv`, `de_report.tsv`

## Reproducibility Check

Run the package test suite:

```bash
cargo test --workspace --release
```

The integration tests execute `ingest -> qc -> matrix -> fold-change` against
the bundled Dube fixture data and compare outputs with the published reference
tables:

- filtered NPX files: byte-exact cell match
- log2 fold-change files: numerical match to floating-point precision

The DE sanity test also verifies that canonical heat-shock proteins increase
in both acute heat comparisons.

## Data Model

Atman writes a small canonical TSV dataset between commands. This schema is
also the adapter target for non-Olink sources:

- `samples.tsv`: sample IDs, subject IDs, conditions, controls, ingest order
- `proteins.tsv`: assay IDs, UniProt IDs, gene symbols, panel metadata
- `measurements.tsv`: one row per sample-assay abundance measurement
- `qc_measurements.tsv`: same schema after QC masking

Downstream commands operate on these files, so each stage can be inspected,
rerun, or replaced independently.

## Layout

```text
crates/atman-core/       core data model and algorithms
crates/atman/            CLI, command orchestration, and file IO
adapters/                canonical TSV adapter helpers and templates
docs/tutorial.md         package tutorial using the bundled Dube fixture
docs/analytical-roadmap.md
                         planned analytical capability build-out
example_data/dube_heat_2023/
                         Dube et al. 2023 Olink Explore fixture data
```

## Reference Dataset

Gagnon D, Barry H, Barhdadi A, Oussaid E, Mongrain I, Lemieux Perreault LP,
Dubé MP. *A dataset of proteomic changes during human heat stress and heat
acclimation.* Scientific Data (2023).
[10.1038/s41597-023-02809-5](https://doi.org/10.1038/s41597-023-02809-5).

## License

Atman is licensed under either [MIT](./LICENSE-MIT) or
[Apache-2.0](./LICENSE-APACHE), at your option.
