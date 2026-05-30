# Atman 1.1.0 — pre-release checklist

Status of the gating items between the cut `v1.1.0` tag and a confident public
release. `[x]` = done in this repo; `[ ]` = needs your action, a decision, or
another machine. See `docs/release-checklist.md` for the generic publish flow.

## 1. Done in this pass (verify, don't redo)

- [x] Version bumped to 1.1.0 — `Cargo.toml`, `crates/atman/Cargo.toml` pin,
      `Cargo.lock`, `CITATION.cff` (+ date 2026-05-30), README docker tags.
- [x] `CHANGELOG.md` `[1.1.0]` section dated; fresh `[Unreleased]` opened.
- [x] Annotated tag `v1.1.0` created (local).
- [x] Kanban orchestration scaffolding removed from the product tree.
- [x] Docs reconciled with code (README command list, ensemble grading,
      recover-plex `--infer-pairs`, detectability flags).
- [x] Full suite green locally (490 passed / 0 failed) and the byte-identity
      determinism gates pass in release profile — on **darwin only**.

## 2. Local verification — DONE (2026-05-30)

- [x] `cargo fmt --all -- --check` — clean (the tree was not fmt-clean; ran
      `cargo fmt --all`, layout-only, no behavior change).
- [x] `cargo clippy --workspace --all-targets -- -D warnings` — **clean (exit 0)**.
      ~40 pre-existing style lints (28 in `atman-core` numeric kernels, ~12 in
      the command layer + test targets) were fixed: real iterator/slice/closure
      rewrites where behavior-identical, scoped `#[allow(...)]` with comments on
      load-bearing lockstep-index loops. Verified byte-identical via the
      determinism + golden suites; reviewed for numeric preservation.
- [x] `cargo test --workspace --release` — 490 passed, 0 failed.
- [x] `cargo package -p atman-core` — clean. `cargo package -p atman` fails
      only on the documented publish-order constraint (needs `atman-core` on
      crates.io first); moot for a local/non-crates.io release.
- [x] `docker build -t atman:1.1.0 .` + `--version`/`--help` — success. Builds
      on `rust:1.94-slim-bookworm`, so this also validated a **Linux release
      build** (see §3).

## 3. Cross-platform / CI (needs another machine or a remote)

This is the biggest open risk: everything so far was validated on one darwin
host, and the repo is local-only so CI has never actually run.

- [ ] Build + test on **Linux** (the CI target / Dockerfile base, rust 1.94).
- [ ] Build + test on **Windows** if it's a supported target (the build.rs
      git-dir resolution and BLAS/thread-pinning are the likely sore spots).
- [ ] If/when a remote is added: push and confirm `.github/workflows/ci.yml`
      goes green with the 1.94 toolchain + thread-pinning env.

## 4. Decisions only you can make

- [ ] **Domain sign-off on the four new commands' docs.** Their CHANGELOG and
      `docs/recipes.md` descriptions (`detectability`, `absence-topology`,
      `recover-plex`, `robust-paired`) were written from reading the code —
      confirm they match intent before publishing.
- [ ] **Release-note the behavior changes** so users aren't surprised:
      - RNG unification changes absolute outputs of `null`, `bootstrap`,
        `ratio`, `decompose null`, `decompose unmix --n-boot`, GSEA permutation,
        `align bootstrap`, Dunnett-Hsu MC (same-seed reruns still
        byte-identical) → regenerate any cached outputs from those commands.
      - Strict QC parsing now errors on unrecognized `qc_sample`/`qc_assay`.
      - `recover-plex --infer-pairs` now hard-fails on a coherent plex.
      - `de --test ensemble` grading is stricter (method agreement, not the
        combined p) — fewer VALIDATED calls.
- [ ] **clippy denial policy** (see §2).

## 5. Not done — out of scope for the fix campaign

- [ ] Security review of input parsing / file handling.
- [ ] Fuzzing of the TSV readers.
- [ ] SBOM generation + dependency/license-compliance audit (note: `statrs`
      pulls `nalgebra` + `rand 0.8` transitively, both unused — documented and
      consciously retained; see `docs/analytical-roadmap.md` "Dependency notes").
- [ ] Confirm `README.pdf` / `mdpdf.log` (gitignored, untracked) are not
      shipped in any release artifact.

## 6. Publish mechanics (when 2–4 are satisfied)

- [ ] `git push` + `git push --tags` (once a remote exists).
- [ ] Create the GitHub release for `v1.1.0` (paste the CHANGELOG section).
- [ ] If publishing to crates.io: `cargo publish -p atman-core` first, then
      `cargo publish -p atman` (the CLI depends on the published core).
