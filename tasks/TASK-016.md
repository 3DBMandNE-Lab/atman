# TASK-016 — delete the stray github/ nested clone

<!-- kanban-status:start -->
**Status:** done
**Owner:** orchestrator
**Branch:** (none — orchestrator cleanup)
**Started:** 2026-05-29T13:07:44+00:00
**Updated:** 2026-05-29T13:07:44+00:00
<!-- kanban-status:end -->

## Plan

(filled in by subagent)

## Notes

Orchestrator cleanup, no code branch: github/ was an untracked full duplicate checkout (own .git, stale CHANGELOG, crates/, build.rs) — a `git clone github` typo for .github. Untracked (git ls-files github/ = 0), no tracked source referenced it. Removed via rm -rf in the main working tree. Nothing to merge.

## Verification

(filled in by subagent)

## Hand-off

(filled in by subagent)
