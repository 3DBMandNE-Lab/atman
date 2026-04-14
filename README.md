# karnaProteome

A local-first Rust engine for proteomics data, extending the karna bioinformatics platform into the protein-abundance modality.

**v0.1 scope:** Ingest, QC-filter, and reproduce the published Dube et al. *Scientific Data* 2023 Olink Explore NGS outputs — byte-for-byte on filtered NPX files and within 1e-4 on log2 fold-change files. No statistical testing, no normalization beyond what Olink already applied, no integrity gate (deferred to v0.2).

## Reproduction recipe

```bash
cargo test --workspace
```

The workspace integration test `dube_reproduction_end_to_end` runs the full pipeline against the raw NPX CSVs in `example_data/dube_heat_2023/` and diffs every output against the published Dube reference files. Observed result on the reference data:

- all 8 filtered NPX files: **byte-exact** cell match
- all 8 log2 fold-change files: **max delta 1.05e-15** (IEEE 754 last-bit drift)
- runtime: ~2.5 s

## Manual pipeline

```bash
karnaproteome ingest \
    --platform olink-explore-ngs --parser dube \
    --output-dir out/ \
    example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
    example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv

karnaproteome qc --input-dir out/ --output-dir out/ --rule dube
karnaproteome matrix --input-dir out/ --output-dir out/ --format dube-wide --split-by panel
karnaproteome fold-change --input-dir out/ --output-dir out/ \
    --groups "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1"
```

## Layout

- `crates/proteome-core/` — platform-agnostic library (types, ingest trait, QC, pivot, fold change)
- `crates/karnaproteome/` — engine binary with four subcommands
- `example_data/dube_heat_2023/` — Dube reference dataset (CC BY 4.0)
- `docs/superpowers/specs/` — design spec
- `docs/superpowers/plans/` — implementation plan

## Reference

Gagnon D, Barry H, Barhdadi A, Oussaid E, Mongrain I, Lemieux Perreault LP, Dubé MP.
*A dataset of proteomic changes during human heat stress and heat acclimation.* Scientific Data (2023).
[10.1038/s41597-023-02809-5](https://doi.org/10.1038/s41597-023-02809-5).
Data: CC BY 4.0.
