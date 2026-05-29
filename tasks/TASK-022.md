# TASK-022 — byte-identity determinism tests across all RNG-bearing commands

<!-- kanban-status:start -->
**Status:** done
**Owner:** orchestrator
**Branch:** worktree-agent-ac545f58ac2987a8a
**Started:** 2026-05-29T21:07:03+00:00
**Updated:** 2026-05-29T21:42:28+00:00
<!-- kanban-status:end -->

## Plan

Extend byte-identity determinism coverage (currently only nmf / missingness-ica /
--adjust-for in `determinism_new_methods.rs`) to every RNG/iterative command, now
that RNG is unified (`atman_core::rng`, TASK-017). New sibling test file
`crates/atman/tests/determinism_rng_commands.rs` mirrors the existing harness:
pin BLAS/OMP threads to 1, run each command twice with a fixed `--seed` into two
separate output dirs, byte-compare the data TSVs, and compare the sidecar SHA
fields (inputs_sha256 + output_files) — which excludes the non-deterministic
sidecar fields (run_uuid, started_at, finished_at). Fixtures reused inline from
each command's own integration test.

## Notes

TEST-ONLY task — no command source changed.

Commands covered (each: two runs, fixed seed, byte-identical outputs):
- `null` (welch-t permutation null, n=100, seed=99) — null_summary.tsv + empirical_p.tsv
- `bootstrap protein` (welch-t, n=200, seed=123)
- `ratio` (wilcoxon, n-bootstrap=500, seed=7)
- `align bootstrap` (k=2, n-boot=10, seed=20260418)
- `de --test ensemble` (welch-t,limma) — de_results.tsv + de_ensemble.tsv
- `decompose ica` (standard FastICA, n-seeds=5, seed=20260418) — loadings/activations/stability
- `robust-paired` (leave-one-pair-out; deterministic by construction — no RNG/seed,
  covered anyway to guard against accidental nondeterminism e.g. map iteration order)

All seven are byte-identical across runs. No command found non-deterministic.

## Verification

- [x] byte-identity tests added for null, bootstrap, ratio, align-bootstrap, de-ensemble, decompose-ica, robust-paired
- [x] each test runs the command twice with a fixed seed and asserts identical outputs, excluding run_uuid/started_at/finished_at
- [x] the tests PASS (`cargo test -p atman --test determinism_rng_commands`: 7 passed)
- [x] `cargo build --workspace` succeeds
- [x] `cargo test --workspace` passes (no failures across the suite)

## Hand-off

New file: `crates/atman/tests/determinism_rng_commands.rs`. No source changes.
Follow-up risk: BLAS-thread pinning is required for ICA/ensemble determinism;
tests set BLAS_NUM_THREADS/OMP_NUM_THREADS/OPENBLAS_NUM_THREADS=1 (matching the
existing harness). robust-paired has no `--seed` (LOSO diagnostic, no RNG); if a
seed flag is ever added it should be threaded through and this test extended.

### Codex adversarial-review fix (orchestrator)

Codex flagged two test-hardening MEDIUMs:
- robust-paired only directly byte-compared its primary output; the derived `<stem>_summary.<ext>` companion was covered only via the sidecar SHA helper. Added a direct `assert_byte_identical` on the summary file (via a `summary_sibling` helper).
- `assert_sidecar_shas_match` could pass vacuously if `output_files` were empty. Added a non-empty guard so an absent/empty output set fails loudly instead of comparing two empty maps.
All 7 determinism_rng_commands tests still pass.
