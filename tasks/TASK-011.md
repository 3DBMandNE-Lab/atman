# TASK-011 — align project: use matrix sample-order, not re-derived order

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** worktree-agent-a22741557ce5801e4
**Started:** 2026-05-29T08:09:39+00:00
**Updated:** 2026-05-29T09:22:13+00:00
<!-- kanban-status:end -->

## Plan

Make the matrix's own sample order the single source of truth for labeling
projected activations, instead of re-deriving sample IDs from samples.tsv in
file order.

- `load_cohort_matrix` now returns `(CohortMatrix, Vec<String>)` where the
  second element is `sample_order` — the same first-seen, non-QC, finite-row
  order the matrix columns were built from.
- `run_project` labels activations with that `sample_order` directly.
- `run_bootstrap` (the other caller) destructures and discards the order; it
  aligns programs across cohorts and never labels individual samples.

## Notes

- `CohortMatrix` lives in `atman-core` and is shared with the bootstrap path,
  so I did NOT add a field to it (would ripple into bootstrap fixtures /
  construction sites). Returning a tuple keeps the change local to the CLI.
- The previous `len() != len()` guard is retained against `transformed.len()`.
- Added a complementary consistency guard: every sample the matrix was built
  from must be declared in samples.tsv (catches the inverse error — a measured
  sample missing from the sample sheet — without letting samples.tsv dictate
  activation labels).
- Output is unchanged for already-agreeing inputs: when samples.tsv file order
  equals the measurement first-seen order (the common case, and what all
  pre-existing tests/fixtures use), the labels are byte-identical. Verified by
  the 4 pre-existing align_project tests still passing unchanged.

## Verification

- [x] activations keyed by the same sample order the matrix columns were built
      from (single source of truth) — `subject_ids = sample_order`.
- [x] test where samples.tsv order differs from measurement first-seen order
      attaches activations to the correct samples:
      `align_project_keys_activations_by_measurement_order_not_samples_tsv`
      (samples.tsv = SUB01..04 ascending; measurement first-seen =
      SUB02, SUB04, SUB01, SUB03; counts coincide so the length guard cannot
      catch it). Confirmed it FAILS against the old samples.tsv-order logic
      ("SUB01: A0001 = 0, expected 1") and PASSES with the fix.
- [x] `cargo build --workspace` — clean.
- [x] `cargo test -p atman --test align_project` — 5 passed.
- [x] `cargo test -p atman --test align_programs` — 2 passed (bootstrap path
      unaffected).

## Hand-off

Code: `crates/atman/src/commands/align.rs` (load_cohort_matrix signature,
run_project labeling + guard, run_bootstrap destructure). Test:
`crates/atman/tests/align_project.rs` (new regression test). No public API or
core-crate changes; bootstrap output and existing projection outputs unchanged.
