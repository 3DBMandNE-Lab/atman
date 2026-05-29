# TASK-012 — robust-paired: zero-sign handling, infer-pairs hard-fail, dedup cluster helpers

<!-- kanban-status:start -->
**Status:** done
**Owner:** orchestrator
**Branch:** worktree-agent-a0b4782ffd7d62c41
**Started:** 2026-05-29T08:09:39+00:00
**Updated:** 2026-05-29T09:40:10+00:00
<!-- kanban-status:end -->

## Plan

1. robust_paired: handle `mean_diff == 0` explicitly for both `target_sign`
   and per-LOSO-replicate sign (Rust `f64::signum(0.0) == +1.0` is wrong here).
2. robust_paired summary: route the two float fields through the empty-cell
   formatter so an undefined median emits an empty cell, not literal `NaN`.
3. recover_plex: make `--infer-pairs` inference failure on a coherent plex a
   hard error; preserve the plex-incoherent SKIP path.
4. Extract the truly-identical cluster helpers into one shared module and use it
   from both `absence_topology.rs` and `recover_plex.rs`.

## Notes

- Behavioral changes (documented in code):
  - **sign-zero:** an exactly-zero `mean_diff` has no direction. `target_sign`
    of a zero full effect is `0.0` (not `+1.0`); a zero LOSO replicate is NOT
    counted as sign-stable against any target. Honest choice: zero is excluded
    from the "matching" count, so a no-direction effect yields sign_stability 0.
  - **NaN→empty:** `fraction_robust_of_significant` and
    `median_sign_stability_of_significant` now go through `opt_f`, so the
    undefined median (no significant proteins) is an empty cell.
  - **infer-pairs hard-fail:** on a `plex_coherent` plex, an inference failure
    now propagates via `?` (was warn + exit 0). The plex-incoherent skip
    (diagnostic != plex_coherent) is unchanged and still exits 0.
- Shared module: `crates/atman/src/commands/jaccard_cluster.rs` holds
  `find`, `union`, `jaccard_distance`, `median_f`, `format_f`, and a
  `quality_path_for(out, default_stem)` (parameterized fallback stem — the only
  difference between the two prior copies). `classify` was NOT shared: the plex
  and protein variants have different signatures and emit different diagnostic
  strings, so each command keeps its own (the task's "confirm truly identical
  before extracting" gate). The shared-helper extraction is byte-identical;
  `median_f` differed only in a local variable name between copies.

## Verification

- `cargo build --workspace` — clean, no warnings.
- `cargo test -p atman --test robust_paired --test absence_topology --test recover_plex`
  — 7/7 pass (robust_paired 2, absence_topology 2, recover_plex 3).
- `cargo clippy -p atman` — no new warnings reference `jaccard_cluster.rs`;
  remaining warnings are pre-existing across the crate.
- No determinism/e2e fixtures reference these three commands, so the refactor
  has no golden-output regression surface; output values are asserted by the
  per-command tests.

## Hand-off

- Acceptance checklist all met. Terminal status: review.
- Follow-up (optional, out of scope): the per-crate clippy warnings predate
  this task and could be swept separately.
