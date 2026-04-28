# Atman Adapters

Atman accepts proteomics data through three routes: platform-native
`ingest` (currently Olink Explore NPX), `ingest-matrix` for common wide
matrices (DIA-NN, Spectronaut, MaxQuant/LFQ, SomaScan, and other log2
abundance exports), and user-written adapters that emit the canonical TSV
files directly. Any proteomics source can feed Atman by producing:

- `samples.tsv`
- `proteins.tsv`
- `measurements.tsv`

The adapters in this directory are release-safe helpers for common matrix
layouts. They do not encode project-specific paths or manuscript assumptions.

## Canonical Schema

`samples.tsv`:

```text
sample_id    subject_id    condition    is_control    sample_type    ingest_order
```

Additional columns are allowed. `atman de --test ols --covariates ...` reads
extra sample columns directly from `samples.tsv`.

`proteins.tsv`:

```text
platform    assay_id    uniprot    gene_symbol    panel    panel_lot
```

`measurements.tsv`:

```text
platform    sample_id    assay_id    gene_symbol    panel    npx_source_str    abundance    abundance_raw    abundance_unit    qc_sample    qc_assay    detection_limit    below_lod    dropped_by_qc    plate_id    panel_lot    ingest_order
```

Well-known platform values (these have platform-specific code paths
elsewhere in atman, e.g. Olink NPX QC, MaxQuant LFQ matrix shape):

- `olink_explore_ngs`
- `olink_target_qpcr`
- `somascan`
- `maxquant_lfq`
- `maxquant_tmt`
- `diann_report`
- `spectronaut_report`

**Custom platforms.** Any non-empty platform string is accepted by the
canonical reader. Adapters for new platforms (e.g. `cptac_tmt_proteome`,
`bruker_timstof_diapasef`, `thermo_proteome_discoverer`) can declare
their own platform identifier without recompiling atman; the canonical
analytical paths (DE, decomposition, alignment, enrichment, signature
scoring) operate on the canonical TSV regardless of platform. Only the
well-known names get platform-specific dispatch.

Supported abundance units include:

- `log2_npx`
- `log2_intensity`
- `log2_lfq`
- `log2_rfu`
- `log2_pg_quantity`
- `log2_diann_pg_quantity`
- `ibaq`
- `raw`

Unknown units are accepted as `raw` by the Rust reader. For differential
abundance, use a scale where additive differences are meaningful. For most
linear intensity exports, log2-transform before writing the TSVs.

## Built-In Matrix Ingest

For common wide matrices, prefer the Rust CLI:

```bash
atman ingest-matrix \
  --matrix pg_matrix.tsv \
  --samples sample_metadata.tsv \
  --output-dir out_ms \
  --platform diann_report \
  --abundance-unit log2_diann_pg_quantity \
  --orientation proteins-rows \
  --assay-id-col Protein.Group \
  --gene-col Genes \
  --uniprot-col Protein.Ids \
  --sample-id-col sample_id \
  --condition-col diagnosis \
  --subject-id-col participant_id \
  --log2-transform
```

Then validate and run Atman:

```bash
atman validate --input-dir out_ms --groups Case-Control --min-pairs 2
atman de --input-dir out_ms --output-dir out_ms \
  --test welch-t --groups Case-Control --min-pairs 2
```

See `examples/` for tiny synthetic fixtures covering Spectronaut, DIA-NN,
MaxQuant/LFQ, and SomaScan-style matrices.

## Python Wide-Matrix Adapter

Use `generic/wide_matrix_to_atman.py` when you have:

- a protein-by-sample or sample-by-protein abundance matrix
- a sample metadata table
- protein metadata embedded in the matrix, or a separate protein metadata table

Examples:

```bash
python adapters/generic/wide_matrix_to_atman.py \
  --matrix pg_matrix.tsv \
  --samples sample_metadata.tsv \
  --out out_ms \
  --platform diann_report \
  --abundance-unit log2_diann_pg_quantity \
  --orientation proteins-rows \
  --assay-id-col Protein.Group \
  --gene-col Genes \
  --uniprot-col Protein.Ids \
  --sample-id-col sample_id \
  --condition-col diagnosis \
  --subject-id-col participant_id \
  --log2-transform
```

Then run Atman on the output directory:

```bash
atman validate --input-dir out_ms --groups Case-Control --min-pairs 5
atman de --input-dir out_ms --output-dir out_ms \
  --test welch-t --groups Case-Control --min-pairs 5
```

For sample-by-protein matrices, provide `--orientation samples-rows` and a
separate `--proteins` metadata table containing the assay ID and optional
gene/UniProt columns.

## CPTAC TMT Proteome Adapter

Use `cptac/cptac_tmt_proteome_to_atman.py` for CPTAC's protein-level TMT
proteome TSVs (the wide `<TUMOR>_proteome.tsv` files distributed by the
CPTAC Data Coordination Center). The adapter handles:

- Skipping the leading `Mean` / `Median` / `StdDev` summary rows.
- Discriminating `Log Ratio` (kept) from `Unshared Log Ratio` (dropped).
- Filtering non-biological columns (`TumorOnlyIR`, `NormalOnlyIR`,
  `QC*`, `Pool*`, `Reference`, `RefMix*`).
- Promoting NCBIGeneID into `assay_id` so cross-cohort joins are stable
  on the gene primary key, with `gene_symbol` carried as metadata.
- Recovering condition labels from sample-ID prefixes when the cohort
  encodes them inline (CPTAC HCC uses `T<N>` for tumor, `P<N>` for the
  paired non-tumor sample from the same subject).

Examples:

```bash
# CPTAC HCC: T<N> / P<N> condition recovery from sample-ID prefixes.
python3 adapters/cptac/cptac_tmt_proteome_to_atman.py \
  --proteome HCC_proteome.tsv \
  --tumor-tag HCC \
  --output-dir out/HCC \
  --sample-id-regex '^(?P<condition_key>[TP])(?P<subject>\d+)$' \
  --condition-key-map 'T=tumor,P=paired_non_tumor'

atman validate --input-dir out/HCC --groups tumor-paired_non_tumor --min-pairs 5
atman de --input-dir out/HCC --output-dir out/HCC \
  --test welch-t --groups tumor-paired_non_tumor --min-pairs 5

# CPTAC GBM / BRCA / COAD / HNSCC: opaque IDs, single-condition cohort.
python3 adapters/cptac/cptac_tmt_proteome_to_atman.py \
  --proteome GBM_proteome.tsv \
  --tumor-tag GBM \
  --output-dir out/GBM
```

The platform column is set to `cptac_tmt_proteome` and the abundance
unit to `log2_tmt_ratio`. Both are accepted by `atman validate`.

## Olink Explore NGS NPX Adapter

Use `generic/olink_explore_to_atman.py` for Olink Explore NGS NPX CSV
exports (the long-format file produced by Olink's analysis software):

```bash
python3 adapters/generic/olink_explore_to_atman.py \
  --output-dir out \
  path/to/file_one_NPX.csv path/to/file_two_NPX.csv
```

Sample IDs are parsed with the Dube et al. 2023 convention by default
(`SSNA-<subject>-<PR1|PR2|PT1|PT2>` for biological samples,
`CONTROL_SAMPLE_*` for technical controls). For datasets with a different
sample-ID layout, pass `--sample-id-regex` (must define named groups
`subject` and `condition`) and `--control-prefix`.

## Adapter Design Rules

- Preserve one biological replicate per `sample_id`.
- Use `subject_id` for paired designs. For unpaired designs, it can equal
  `sample_id`.
- Put the comparison label used by `atman de --groups A-B` in `condition`.
- Set `is_control` to `0` for biological samples and `1` for technical controls.
- Use `dropped_by_qc=1` for cells that should be excluded from analysis.
- Keep missing abundance cells blank; Atman will drop them per protein/test.
