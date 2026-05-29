# TASK-019 — de/mod: fix stale doc comments + dedup provenance hash list

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** (pending)
**Started:** 2026-05-29T13:08:21+00:00
**Updated:** 2026-05-29T13:14:13+00:00
<!-- kanban-status:end -->

## Plan

Doc/comment + provenance-list fixes in `crates/atman/src/commands/de/mod.rs`. No DE behavior or numeric-output changes.

1. Update `--paired-by` doc: it claimed the key is "fixed to subject_id", but the override to any samples.tsv column is implemented (`override_subject_id`, dispatched ~line 582).
2. Update `--post-hoc` doc + the unknown-method bail message: both claimed tukey/dunnett were "pending"/"rejected until DEBT-5/DEBT-6"; both are fully implemented (posthoc.rs) and dispatched (`run_posthoc_tukey`/`run_posthoc_dunnett`, ~lines 568-572).
3. Dedup the duplicated `"measurements.tsv"` in the single `hash_canonical_inputs` call (~line 1268).

## Notes

- Only one `hash_canonical_inputs` call exists in de/mod.rs (other DE methods own theirs in submodules; out of scope per task). It listed `measurements.tsv` twice; now `["measurements.tsv", "samples.tsv", "proteins.tsv"]`.
- Provenance output is byte-identical before/after the dedup: `hash_canonical_inputs` returns an `InputHashes = BTreeMap<String, String>` keyed by basename, so the duplicate `measurements.tsv` already collapsed to one map entry. The dedup is a source-level cleanup with zero effect on the sidecar bytes.
- Also fixed the stale `other => bail!("...sidak (tukey / dunnett pending)")` message to `"supported: sidak, tukey, dunnett"`, since it was part of the same stale-doc surface.
- No DE numeric output touched.

## Verification

- `cargo build --workspace`: clean (Finished dev profile).
- `cargo test -p atman --test de_paired_by --test de_posthoc_tukey --test de_posthoc_dunnett --test de_limma`: 4 + 3 + 3 + 7 = all pass.
- Full `cargo test -p atman`: all suites pass, 0 failed.
- Acceptance checklist:
  - [x] `--paired-by` / `--post-hoc` help/doc reflect current behavior (no false "fixed to subject_id" / "rejected until DEBT" claims)
  - [x] duplicate `measurements.tsv` removed from de/mod.rs provenance hash list
  - [x] DE numeric output unchanged (doc/comment + benign provenance-list dedup only)
  - [x] cargo test passes

## Hand-off

Determinism: only the source-level provenance input list changed; sidecar bytes are unchanged because the BTreeMap already deduped the repeated basename. No follow-ups required. Other DE submodules (posthoc.rs, limma.rs, msqrob.rs) own their own `hash_*_inputs` calls and were not in scope.
