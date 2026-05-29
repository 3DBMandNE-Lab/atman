# TASK-018 — validate: dedup the provenance input-hash list

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** worktree-agent-af93479cd8e11bb2f
**Started:** 2026-05-29T13:08:21+00:00
**Updated:** 2026-05-29T13:09:39Z
<!-- kanban-status:end -->

## Plan

Remove the duplicate `measurements.tsv` from the `hash_canonical_inputs`
call in `validate.rs` (~line 122), leaving the intended distinct set
`["measurements.tsv", "samples.tsv", "proteins.tsv"]`.

## Notes

- `hash_canonical_inputs(input_dir, filenames)` (run_sidecar.rs:57) maps
  each filename to `input_dir.join(name)` then calls `hash_labeled_inputs`,
  which inserts into a `BTreeMap` keyed by label (basename). The duplicate
  `measurements.tsv` was silently overwritten, so the hashed set was already
  effectively `{measurements, samples, proteins}` — the fix only removes the
  redundant entry; it does not drop any genuinely-distinct input.
- Only one `hash_canonical_inputs` call exists in validate.rs; no other
  duplicate provenance hash entries in this file.
- Determinism: the recorded provenance input list no longer contains the
  duplicate. Since the duplicate was overwritten in the BTreeMap, the
  resulting `inputs_sha256` map is unchanged. Benign, no sidecar hash change.

## Verification

- [x] validate.rs provenance input-hash list has no duplicate entries
- [x] hashed input set is the intended {measurements, samples, proteins}
- [x] `cargo build --workspace` — clean
- [x] `cargo test -p atman --test validate` — 5 passed, 0 failed

## Hand-off

Single-line change in `crates/atman/src/commands/validate.rs`. No behavior
change to recorded hashes (duplicate was already deduped by BTreeMap key);
input list is now accurate. Ready for review.
