# TASK-004 — null: use numerically stable two-sided p-value

<!-- kanban-status:start -->
**Status:** done
**Owner:** orchestrator
**Branch:** worktree-agent-a7f43436fe4ec7450
**Started:** 2026-05-29T06:43:48+00:00
**Updated:** 2026-05-29T07:37:11+00:00
<!-- kanban-status:end -->

## Plan

1. Read `t_from_diffs` (null path), `observed_test` (observed path), and the core
   `de` two-sided-p helper to find the canonical implementation.
2. Replace the unstable `2.0 * (1.0 - dist.cdf(t.abs()))` in `t_from_diffs` with
   the canonical core helper `atman_core::de::two_sided_t_p_value`, matching the
   observed path which routes through `paired_t`/`welch_t` (both use the stable
   survival function).
3. Drop the now-unused `StudentsT`/`ContinuousCDF` import from `null.rs`.
4. Dedup the duplicated `"measurements.tsv"` entry in the
   `hash_canonical_inputs` provenance list.
5. Add a unit test exercising a large-|t| case that asserts a non-zero p where
   the old formula returned exactly 0.0.
6. Build the workspace and run the null + core test suites.

## Notes

- The canonical helper already exists: `atman_core::de::two_sided_t_p_value`
  (crates/atman-core/src/de.rs:28), which computes `2 * dist.sf(|t|)` clamped to
  [0,1] and handles `df = Inf` via the normal distribution. It is the same helper
  the observed path reaches through `paired_t`/`welch_t`, so reusing it puts the
  null and observed p-values on identical numerical footing. Core's existing
  regression test `two_sided_t_p_value_uses_stable_tail_for_extreme_t` covers it.
- `t_from_diffs` was the only consumer of `StudentsT`/`ContinuousCDF` in
  `null.rs`; removing the formula made that import unused, so it was deleted (no
  unused-import warning on build).
- Determinism impact: null p-values in the extreme tail change from 0.0 to a
  small positive value (the intended fix). This flows through
  `null_p -> bh_fdr -> null_q -> mean_null_hits/empirical_fdr`, so null-summary
  outputs in the tail will differ from pre-fix runs. This is the correct,
  observed-path-consistent behavior.
- Provenance impact: deduping `measurements.tsv` changes the input list passed to
  `hash_canonical_inputs`. The hash itself is over file *contents* keyed by name;
  removing a duplicate name removes one redundant hashed entry, so the recorded
  `inputs_sha256` in the run sidecar may change versus pre-fix runs. Noted as a
  benign, one-time provenance shift.

## Verification

Acceptance criteria:
- [x] null path uses the stable sf-based two-sided p-value via the canonical core
      helper `atman_core::de::two_sided_t_p_value` (same helper as observed path).
- [x] duplicate `measurements.tsv` removed from `null.rs` provenance hash list.
- [x] test `t_from_diffs_large_t_returns_nonzero_p` exercises a large-|t| case,
      asserts the naive `2*(1-cdf(|t|))` underflows to 0.0 for the fixture, and
      asserts the new path returns a strictly positive p equal to the core helper.
- [x] cargo build + tests pass.

Commands:
- `cargo build --workspace` -> Finished, clean (no warnings).
- `cargo test -p atman --lib commands::null` -> 1 passed (new unit test).
- `cargo test -p atman --test null` -> 3 passed.
- `cargo test -p atman-core` -> all passed (incl. extreme-t regression test).
- `cargo test -p atman` (full) -> all suites passed, 0 failed.

## Hand-off

Single-file code change in `crates/atman/src/commands/null.rs`:
- `t_from_diffs` now calls `two_sided_t_p_value(t, df)` instead of the unstable
  cdf-subtraction; returns `Skipped` if the helper yields `None`.
- Removed unused `statrs` import.
- Deduped provenance input list.
- Added `#[cfg(test)] mod tests` with the large-|t| regression test.

Risks/follow-ups:
- Pre-fix null-summary fixtures (if any are checked in as golden files) with
  extreme-tail p-values will need regeneration; none surfaced as failing tests.
- Reviewers may want a note in the changelog about the null-path numerical fix
  and the provenance-hash input-list change.
