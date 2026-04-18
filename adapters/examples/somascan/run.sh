#!/usr/bin/env bash
set -euo pipefail

out="target/adapter_examples/somascan"
rm -rf "$out"
atman_bin="${ATMAN_BIN:-atman}"

"$atman_bin" ingest-matrix \
  --matrix adapters/examples/somascan/abundance.tsv \
  --samples adapters/examples/somascan/sample_metadata.tsv \
  --proteins adapters/examples/somascan/feature_metadata.tsv \
  --output-dir "$out" \
  --orientation samples-rows \
  --platform somascan \
  --abundance-unit log2_rfu \
  --sample-id-col sample_id \
  --subject-id-col subject_id \
  --condition-col condition \
  --assay-id-col aptamer_id \
  --gene-col gene_symbol \
  --uniprot-col uniprot \
  --panel somascan

"$atman_bin" validate --input-dir "$out" --groups Case-Control --min-pairs 2
