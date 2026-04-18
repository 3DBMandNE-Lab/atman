#!/usr/bin/env bash
set -euo pipefail

out="target/adapter_examples/diann"
rm -rf "$out"
atman_bin="${ATMAN_BIN:-atman}"

"$atman_bin" ingest-matrix \
  --matrix adapters/examples/diann/protein_groups.tsv \
  --samples adapters/examples/diann/sample_metadata.tsv \
  --output-dir "$out" \
  --orientation proteins-rows \
  --platform diann_report \
  --abundance-unit log2_diann_pg_quantity \
  --sample-id-col sample_id \
  --subject-id-col subject_id \
  --condition-col condition \
  --assay-id-col Protein.Group \
  --gene-col Genes \
  --uniprot-col Protein.Ids \
  --panel diann \
  --log2-transform

"$atman_bin" validate --input-dir "$out" --groups Case-Control --min-pairs 2
