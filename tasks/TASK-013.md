# TASK-013 — atomic_write sweep for remaining fs::write callers

<!-- kanban-status:start -->
**Status:** done
**Owner:** orchestrator
**Branch:** worktree-agent-ac478b220efa31296
**Started:** 2026-05-29T09:49:45+00:00
**Updated:** 2026-05-29T10:19:38+00:00
<!-- kanban-status:end -->

## Plan

Route the remaining direct `std::fs::write` command-output calls in three
command modules through `crate::io::atomic_write` (tempfile + atomic rename),
matching the rest of the toolset (robustness.rs converted in TASK-010). Pure
mechanism change; output bytes unchanged.

## Notes

Converted one output write per file (each had exactly one `std::fs::write`):
- `asymmetry.rs` (~line 177): asymmetry metrics TSV → `crate::io::atomic_write`
- `module_de.rs` (~line 259): `module_de_results.tsv` → `crate::io::atomic_write`
- `module_trajectory.rs` (~line 123): module trajectory TSV → `crate::io::atomic_write`

Used the fully-qualified `crate::io::atomic_write` (no import changes needed).
`&args.output` / `&output_path` (`&PathBuf`) coerces to the `&Path` parameter.
The `out` byte buffer passed is identical to before, so output is byte-for-byte
unchanged. No other `std::fs::write` calls remain in these three files.

## Verification

- `grep "fs::write"` on the three files: no matches remaining.
- `cargo build --workspace`: clean.
- `cargo test -p atman`: all suites pass, 0 failures.

## Hand-off

Mechanism-only change; safe to merge. No behavioral or output changes.
