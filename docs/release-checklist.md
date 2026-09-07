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
     qc --input-dir /out --output-dir /out
   ```

6. Tag the release commit:

   ```bash
   git tag -a v1.0.0 -m "Atman v1.0.0"
   ```

## Publishing a polished tree, and the provenance join

If the public release is a cleaned or squashed history, its commit sha
will not be the sha recorded in any analysis that used this tree. Every
run sidecar carries `atman_git_sha`, and a deposited analysis pointing at
a commit that does not exist publicly is a broken provenance link — in a
tree whose whole claim is that the links hold.

Close it one of two ways.

**Preferred: re-run the analysis at the published commit.** Tell the
analysis owner the sha once the release exists. This is the clean option
and it needs no argument in the methods.

**Fallback: prove source identity.** If a re-run is not possible, show
that the published tree builds the same numerics as the recorded commit.

```bash
# sources only
find crates -name '*.rs' -type f | sort | xargs shasum -a 256 | shasum -a 256

# sources + manifests + lockfile + toolchain  (prefer this one)
{ find crates -name '*.rs' -type f; ls Cargo.toml Cargo.lock rust-toolchain.toml; \
  find crates -name 'Cargo.toml'; } | sort -u | xargs shasum -a 256 | shasum -a 256
```

Use the second. The first covers only `crates/**/*.rs`, which does
include both `build.rs` files, but omits three things that change the
numbers: `Cargo.lock` pins dependency versions, the `Cargo.toml` files
carry features and profile settings that change codegen, and
`rust-toolchain.toml` pins the compiler. A digest narrower than the claim
it is cited for is the failure this project spent a release removing.

Recorded for commit `a92f380` (version 1.2.0):

| Scope | Digest | Use |
|---|---|---|
| sources only | `952cb1a7243060ec83c72bf78bb4460b5b651845994f36aa5fa2986c08c8f2b0` | **INSUFFICIENT — do not cite** |
| sources + manifests + lockfile + toolchain | `adcbad75e7850197967cdd4d20bf41a87b3f8d891c8a3a56580bd7832cd54269` | **Use this one** |

The narrow digest is kept only so a reader who computes it can tell which
of the two they have. It must not be cited as evidence of identical
numerics: two trees matching on it can still differ in `Cargo.lock`, the
manifests or the toolchain, and would then produce different numbers
while the digest called them identical. Labelling it rather than deleting
it is deliberate — an unlabelled digest in a file is something someone
picks up and uses.

Both were computed on a clean tree. `target/` is excluded by
construction, and the digest is path-sensitive, so the published tree
must keep the same layout for the comparison to mean anything.

## What a polished release must preserve

Trimming development history is fine. These are not history:

1. **The version string.** Sidecars record `atman_version`; if the
   release disagrees, that is the first thing a reader checks and the
   first thing that looks wrong.
2. **`bench/adapters/` and the bench fixtures.** The scikit-learn
   comparison is a published figure. A reader checking it needs the
   adapters, which are about twenty lines each.
3. **The platform-divergence documentation.** `README.md` "Determinism
   and Your Operating System" and the `docs/reference.md` section it
   points to. Reproducibility statements cite it, and it is the honest
   form of the determinism claim rather than the overclaim it replaced.
