# Atman Release Checklist

Use this checklist before cutting a standalone GitHub release.

## Required Checks

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo package -p atman-core
```

`cargo package -p atman` requires `atman-core` to already exist in the registry
because the CLI crate depends on `atman-core = 1.0.0`. For crates.io releases,
publish and verify in this order:

```bash
cargo publish -p atman-core
cargo package -p atman
cargo publish -p atman
```

For a GitHub-only release, the package tree is still validated by the workspace
checks above and by the Docker build below.

## GitHub Release

1. Confirm `CHANGELOG.md` has the final release date and version.
2. Confirm `README.md`, `docs/tutorial.md`, and `adapters/README.md` describe
   Atman as a standalone tool.
3. Confirm no hidden system files are present:

   ```bash
   find . -name .DS_Store -print
   ```

4. Build the container image:

   ```bash
   docker build -t atman:1.0.0 .
   ```

5. Run the bundled Dube reproduction path. Atman has no native NPX
   ingest; the Olink Explore adapter is a Python script in `adapters/`.
   Run it on the host (or in a Python container), then run Atman on the
   resulting canonical TSVs:

   ```bash
   python3 adapters/generic/olink_explore_to_atman.py \
       --output-dir out \
       example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
       example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv

   docker run --rm -v "$(pwd)/out:/out" atman:1.0.0 \
     qc --input-dir /out --output-dir /out --rule mask-warn-fail
   ```

6. Tag the release commit:

   ```bash
   git tag -a v1.0.0 -m "Atman v1.0.0"
   ```
