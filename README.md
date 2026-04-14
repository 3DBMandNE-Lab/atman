# karnaProteome

A local-first Rust engine for proteomics data, extending the karna bioinformatics platform into the protein-abundance modality.

## Scope

**v0.1 — reproduction base** (done): ingest, QC-filter, and reproduce the published Dube et al. *Scientific Data* 2023 Olink Explore NGS outputs byte-for-byte on filtered NPX files and within 1e-4 on log2 fold-change files. No normalization beyond what Olink already applied.

**v0.2 analysis — pass 1** (done): paired Student's t-test at subject level with BH-FDR per comparison family (`karnaproteome de`). Heat-shock instrument sanity check. Findings doc on the Dube heat-acclimation dataset.

**v0.2 analysis — pass 2** (done): pathway enrichment via Fisher's exact ORA against MSigDB Hallmark / KEGG / GO-BP, with universe correction to the Olink Explore measurable gene universe. Retracted the pass-1 headline after discovering the full MSigDB background inflated Olink-overlap-heavy pathways (EMT is 53% Olink-covered).

**v0.2 analysis — pass 3** (done, null): per-protein and per-pathway OLS regression against Dube's physiological covariates (Tcore, Tskin, LSR, LDF, CVC, HR, BP, SSNA). Genomewide phenotype regression doesn't work at n=10. External calibration against Dube's own per-subject delta NPX files matched at 1e-6 precision.

**Deferred (v0.2+):** karnadata integration, bridge normalization, additional ingest adapters (Olink Target qPCR, SomaScan, DIA-NN, Spectronaut).

## Reproduction

```bash
cargo test --workspace --release
```

50 tests pass. The workspace integration test `dube_reproduction_end_to_end` runs `ingest → qc → matrix → fold-change` against the raw NPX CSVs in `example_data/dube_heat_2023/` and diffs every output against the published Dube reference files (~2.5s):

- all 8 filtered NPX files: byte-exact cell match
- all 8 log2 fold-change files: max delta 1.05e-15 (IEEE 754 last-bit drift)

The workspace integration test `dube_de_heat_shock_sanity` runs `ingest → qc → de` against the same raw data and asserts canonical heat-shock proteins (HSPB1, HSPA1A, DNAJB1, HSPG2) go up in both acute heat comparisons. 4/4 up in PT1−PR1 and PT2−PR2.

## Manual pipeline

```bash
# Stage 1–4: reproduction base
karnaproteome ingest \
    --platform olink-explore-ngs --parser dube \
    --output-dir out/ \
    example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
    example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv
karnaproteome qc --input-dir out/ --output-dir out/ --rule dube
karnaproteome matrix --input-dir out/ --output-dir out/ --format dube-wide --split-by panel
karnaproteome fold-change --input-dir out/ --output-dir out/ \
    --groups "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1"

# Stage 5: analysis
karnaproteome de --input-dir out/ --output-dir out/ \
    --test paired-t --paired-by participant \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" --min-pairs 5

# Stage 6: pathway enrichment (two passes)
scripts/run_enrichment.sh out/de_results.tsv docs/findings/pathway_results/
scripts/enrich_restricted_universe.py out/de_results.tsv \
    docs/findings/pathway_results_restricted/

# Stage 7: phenotype regression (null, for the record)
scripts/phenotype_regression.py docs/findings/phenotype_results/
scripts/phenotype_pathway_regression.py out/de_results.tsv \
    docs/findings/phenotype_results/
```

The `scripts/run_enrichment.sh` Bash driver and the two Python scripts shell out to the sibling `genesets` Rust CLI for gene-set data (set `KARNA_GENESETS=/path/to/genesets` if it's not at the default).

## Layout

```
crates/proteome-core/       platform-agnostic library
  src/
    types.rs                Platform, AssayId, MeasurementRecord, ...
    errors.rs               IngestError
    sample_id.rs            SampleIdParser trait + DubeSampleIdParser
    ingest/olink_explore.rs Olink Explore NGS long-CSV adapter
    qc.rs                   Dube QC rule (mask on QC_Warning or Assay_Warning != PASS)
    matrix.rs               long → Dube-wide pivot (gene-symbol keyed, string values)
    fold_change.rs          Neumaier-compensated log2 FC
    de.rs                   paired Student's t + BH-FDR

crates/karnaproteome/       engine binary
  src/
    main.rs                 clap dispatch
    io.rs                   atomic write, long-TSV read/write, per-format CSV writers
    commands/
      ingest.rs             karnaproteome ingest
      qc.rs                 karnaproteome qc
      matrix.rs             karnaproteome matrix
      fold_change.rs        karnaproteome fold-change
      de.rs                 karnaproteome de
  tests/
    dube_reproduction.rs    acceptance: strict diff 8 filtered NPX + numeric diff 8 FC
    dube_de_sanity.rs       acceptance: canonical HSPs up in heat comparisons

scripts/                    Python + Bash analysis drivers (not part of the Rust crates)
  run_enrichment.sh         Fisher's exact ORA via sibling genesets CLI
  enrich_restricted_universe.py   universe-corrected ORA
  phenotype_regression.py         per-protein OLS vs Dube phenotypes
  phenotype_pathway_regression.py per-pathway OLS vs Dube phenotypes

example_data/dube_heat_2023/   all 5 Figshare items for project 163291 (CC BY 4.0)
  20212016_Dube_NPX_2021-11-30.csv               raw NPX, Explore 1536
  20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv   raw NPX, Expansion
  filtered_npx/npx/<panel>_npx.csv × 8           Dube's published filtered NPX
  fold_changes/fold_changes/<panel>_log2_fc.csv × 8    Dube's published log2 FC
  delta_npx/dNPX_{hd1,hd7,accpre,accpost}.csv    Dube's per-subject delta NPX
  Physiological_data.xlsx                        per-subject body temp / sweat / HR / BP / SSNA

docs/superpowers/
  specs/2026-04-14-karnaproteome-olink-reproduction-design.md
  plans/2026-04-14-karnaproteome-olink-reproduction.md

docs/findings/
  2026-04-14-dube-heat-acclimation-findings.md   per-protein + pathway + phenotype write-up
  pathway_results/                                ORA pass 1 JSONs (inflated universe)
  pathway_results_restricted/                     ORA pass 2 TSVs (Olink universe)
  phenotype_results/                              per-protein + per-pathway OLS results
```

## Reference

Gagnon D, Barry H, Barhdadi A, Oussaid E, Mongrain I, Lemieux Perreault LP, Dubé MP.
*A dataset of proteomic changes during human heat stress and heat acclimation.* Scientific Data (2023).
[10.1038/s41597-023-02809-5](https://doi.org/10.1038/s41597-023-02809-5).
Data CC BY 4.0 via Figshare project [163291](https://figshare.com/projects/163291).
