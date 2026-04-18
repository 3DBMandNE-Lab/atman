# Adapter Examples

These are tiny synthetic fixtures for common non-Olink matrix layouts. They are
not biological reference datasets; they exist to show source-to-Atman mappings.

Each example can be converted with `atman ingest-matrix` and then checked with
`atman validate`.

```bash
bash adapters/examples/spectronaut/run.sh
bash adapters/examples/diann/run.sh
bash adapters/examples/maxquant/run.sh
bash adapters/examples/somascan/run.sh
```

The scripts write into `target/adapter_examples/`, which is ignored by git.
