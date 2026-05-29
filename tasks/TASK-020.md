# TASK-020 — split god-module decompose.rs (3,349 LOC)

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** worktree-agent-aa4e0cbbbb8fd332d
**Started:** 2026-05-29T21:07:03+00:00
**Updated:** 2026-05-29T21:31:11Z
<!-- kanban-status:end -->

## Plan

Pure refactor of `commands/decompose.rs` (3,349 LOC) into a `commands/decompose/`
submodule tree, mirroring the `commands/de/` split (private `mod` + selective `use`).
One file per subcommand; shared types stay with their natural owner.

## Notes

Resulting layout (`crates/atman/src/commands/decompose/`):

| file               | LOC | contents |
|--------------------|-----|----------|
| `mod.rs`           |   65 | `Args`, `Command` enum, `run()` dispatch, shared `program_name()` |
| `ica.rs`           | 1022 | `IcaArgs`, `TransformArg`/`StabilityMetric`/`MissingnessModel`, `run_ica`, `resolve_k`, `parse_k_selection`, `derive_cohort`, `AbundanceMatrix`/`SampleInfo`/`AssayInfo`, `load_matrix`, `load_matrix_mnar`, `apply_compositional`, `canonicalize`, writers |
| `nmf.rs`           |  722 | `NmfArgs`, `nmf_run` |
| `null.rs`          |  243 | `NullModeArg`, `NullArgs`, `run_null` |
| `variance.rs`      |  544 | `VarianceArgs`, `run_variance`, formula parse, design build, TSV writers |
| `unmix.rs`         |  669 | `UnmixArgs`, `run_unmix`, marker/k-selection/endmember/abundance writers, `load_subject_protein_matrix` |
| `counterfactual.rs`|  162 | `CounterfactualArgs`, `run_counterfactual` |

Code moved verbatim; only module wiring changed. Cross-module sharing:
`null.rs` reuses `ica::{load_matrix, IcaArgs, MissingnessModel, StabilityMetric,
TransformArg}` (IcaArgs fields made `pub(super)` so the null adapter can construct
it directly); `program_name` lives in `mod.rs` as `pub(super)` and is used by
`ica`/`nmf`/`null`. `AbundanceMatrix` and its `samples`/`assays`/`data` fields made
`pub(super)` for `null.rs`. No logic, constant, ordering, RNG, or formatting change.

## Verification

- `cargo build --workspace`: clean, zero warnings.
- All decompose_* tests pass: decompose_ica, decompose_nmf, decompose_unmix,
  decompose_null, decompose_variance, decompose_counterfactual,
  decompose_ica_transforms, decompose_variance_type3, decompose_null_dube,
  align_programs_nmf.
- `determinism_new_methods`: 3/3 pass (byte-identity for nmf, missingness-ica,
  --adjust-for).
- `cargo test --workspace`: 41 test binaries, all `ok`, 0 failed (exit 0).

## Hand-off

Pure structural refactor; output byte-identical. `commands::decompose::run` and
`commands::decompose::Args` signatures preserved (only external references, from
`main.rs`). No follow-ups required.
