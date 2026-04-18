#!/usr/bin/env bash
set -euo pipefail

out="target/adapter_examples/spectronaut"
rm -rf "$out"
atman_bin="${ATMAN_BIN:-atman}"

"$atman_bin" ingest-matrix \
  --matrix adapters/examples/spectronaut/protein_groups.tsv \
  --samples adapters/examples/spectronaut/sample_metadata.tsv \
  --output-dir "$out" \
  --orientation proteins-rows \
  --platform spectronaut_report \
  --abundance-unit log2_intensity \
  --sample-id-col sample_id \
  --subject-id-col subject_id \
  --condition-col condition \
  --assay-id-col Protein.Group \
  --gene-col Genes \
  --uniprot-col Protein.Ids \
  --panel-col Panel

"$atman_bin" validate --input-dir "$out" --groups Case-Control --min-pairs 2
