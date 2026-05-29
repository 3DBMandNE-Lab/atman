# TASK-010 — robustness: deterministic top-k tie-break + atomic writes

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** worktree-agent-a90acf904996fab92
**Started:** 2026-05-29T08:09:39+00:00
**Updated:** 2026-05-29T09:42:13+00:00
<!-- kanban-status:end -->

## Plan

- Add a deterministic secondary key to `top_k_set`'s sort and make the primary
  comparator a total order (NaN handled), so tied `|mean_diff|` at the k-boundary
  resolves to a fixed lexicographic feature-id order.
- Route the three robustness output writes through `crate::io::atomic_write`.
- Add an integration test with tied effect sizes asserting byte-reproducibility
  across two runs.

## Notes

Nondeterminism documented and fixed in `crates/atman/src/commands/robustness.rs`:

1. **`top_k_set` (primary target).** Pairs `(feature_id, |mean_diff|)` are
   collected from a `HashMap<String, DeLite>` and were sorted by magnitude only
   with `partial_cmp` and no secondary key. Tied magnitudes at the k-boundary
   were resolved by HashMap iteration order → non-reproducible top-k sets →
   non-reproducible overlap/Jaccard in `rank_stability.tsv`. Fix: sort by
   `b.total_cmp(&a)` (total order incl. NaN) `.then_with(|| a.0.cmp(b.0))`
   (ascending feature id).

2. **`loo_sign_stability.tsv` row order.** The per-contrast loop iterated
   `base_rows` (a HashMap) directly, so emitted row order varied run-to-run. The
   tie-break test surfaced this. Fixed by iterating sorted keys.

3. **`stability_ranked.tsv` source order.** Same HashMap-iteration issue feeding
   the `Vec<StabRow>`; the rank sorts are stable `sort_by`, so a deterministic
   feature-id initial order makes tied scores resolve deterministically. Fixed
   by iterating sorted keys.

Writes now go through `crate::io::atomic_write` (tempfile + atomic rename) for
crash-safety, consistent with the rest of the toolset.

## Verification

- [x] top-k sort has a deterministic secondary key; total-order primary
      comparator (`total_cmp` handles NaN)
- [x] tied-boundary selection reproducible
- [x] all three robustness outputs written via `crate::io::atomic_write`
- [x] new test `crates/atman/tests/robustness_topk_ties.rs` — tied effect sizes,
      byte-identical `rank_stability.tsv` / `loo_sign_stability.tsv` /
      `stability_ranked.tsv` across two runs; passes
- [x] `cargo build --workspace` clean
- [x] `cargo test -p atman` → 216 passed, 0 failed

## Hand-off

Changes confined to `crates/atman/src/commands/robustness.rs` (top_k_set +
two sorted-key iterations + atomic_write) and new test
`crates/atman/tests/robustness_topk_ties.rs`. No CLI/flag/output-schema changes;
only ordering/atomicity. Branch name is the generic worktree name; rename on
merge if a conventional branch name is desired.

### Codex adversarial-review fix (orchestrator)

Codex flagged (MEDIUM): `b.1.total_cmp(&a.1)` sorts NaN magnitudes ahead of all finite values in descending order, so a non-finite `mean_diff` (parse_opt_f64 accepts "NaN"/"inf") could hijack the top-k set. Fixed: `top_k_set` now filters non-finite `|mean_diff|` out of the candidate set before ranking (a non-finite effect is not a meaningful top hit; strict ranking). Added unit test `top_k_set_excludes_non_finite_and_breaks_ties_by_id` (NaN/inf excluded; tie resolves to ascending id). Integration tie-break byte-identity test still green.
