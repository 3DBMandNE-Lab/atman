#!/usr/bin/env bash
set -euo pipefail

out="target/adapter_examples/maxquant"
rm -rf "$out"
atman_bin="${ATMAN_BIN:-atman}"

"$atman_bin" ingest-matrix \
  --matrix adapters/examples/maxquant/protein_groups.tsv \
  --samples adapters/examples/maxquant/sample_metadata.tsv \
  --output-dir "$out" \
  --orientation proteins-rows \
  --platform maxquant_lfq \
  --abundance-unit log2_lfq \
  --sample-id-col sample_id \
  --subject-id-col subject_id \
  --condition-col condition \
  --assay-id-col "Protein IDs" \
  --gene-col "Gene names" \
  --uniprot-col "Protein IDs" \
  --panel maxquant_lfq

"$atman_bin" validate --input-dir "$out" --groups Case-Control --min-pairs 2
