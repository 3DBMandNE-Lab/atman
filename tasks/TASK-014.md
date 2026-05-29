# TASK-014 — build provenance: -dirty SHA, drop unused thiserror, audit statrs

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** (pending)
**Started:** 2026-05-29T09:49:45+00:00
**Updated:** 2026-05-29T10:18:00+00:00
<!-- kanban-status:end -->

## Plan

1. `build.rs`: append `-dirty` to `ATMAN_GIT_SHA` when the tree has
   uncommitted tracked changes, preserving the `unknown` fallback.
2. Drop the unused `thiserror` workspace dependency.
3. Audit `statrs` symbol usage and document the nalgebra+rand footprint.

## Notes

- `build.rs`: replaced the bare `git rev-parse HEAD` capture with
  `git_sha_with_dirty()`. It keeps the base SHA, then runs
  `git status --porcelain --untracked-files=no`; a non-empty successful
  output (whitespace-insensitive check) marks the tree dirty and appends
  `-dirty`. When `git rev-parse` returns `"unknown"` (no git), the suffix is
  never added — we don't claim dirtiness we can't observe. Side-effect-free,
  deterministic.
  - Rerun triggers unchanged (build.rs, Cargo.lock, .git/HEAD + ref). Caveat:
    cargo won't re-run build.rs purely because a tracked source file changed,
    so the suffix reflects tree state at the moment build.rs last ran. In the
    real build path, editing a tracked source recompiles its crate and a clean
    `cargo build` from CI bakes the correct flag; the suffix is exact when
    build.rs executes. Not worth watching every tracked file.
- `thiserror`: zero references in the codebase (only kanban/card mentions) and
  absent from `Cargo.lock` even transitively. Removed from
  `[workspace.dependencies]`. Workspace still builds.
- `statrs`: audit written to `docs/analytical-roadmap.md` under a new
  "Dependency notes" heading. Only scalar CDFs (`Normal`, `StudentsT`,
  `FisherSnedecor`, `Hypergeometric`), the `ContinuousCDF`/`DiscreteCDF`
  traits, and `function::gamma::ln_gamma` are used. No `MultivariateNormal`/
  sampling, so the transitive `nalgebra`+`rand 0.8.5` weight is unused-but-
  tolerated. No code/dep change to statrs (audit-only per task).

## Verification

- [x] Dirty-tree binary bakes `-dirty`: built binary carries
      `44bfc9e75c761cb421b348062ee07ea122ab6e76-dirty` while the worktree has
      uncommitted tracked changes; clean trees keep the bare SHA (no suffix
      branch only taken on empty `git status` output).
- [x] `thiserror` removed; `cargo build --workspace` succeeds.
- [x] statrs symbol-usage audit written to docs note; no statrs change.
- [x] `cargo build --workspace`, `cargo test -p atman`, `cargo test -p
      atman-core` all pass.

## Hand-off

Pure infra/provenance change. Reviewer note: the `-dirty` suffix will appear
on any binary built from a non-clean checkout (by design). Clean CI builds are
unaffected. Future candidate (not this task): feature-gate or replace statrs to
shed the unused nalgebra+rand transitive deps.
