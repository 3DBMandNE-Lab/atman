# Atman feature requests

No open requests. All filed items from the 2026-04-19 CSF cross-disease
batch shipped the same day.

## How to implement a feature request

Per-feature loop — do not batch. Ship one, then pick up the next.

1. **Read the whole request.** Command sketch, "what it does", "why we
   need it", scope hints. The CLI contract matters more than the
   internal algorithm.

2. **Flag constraint conflicts before coding.** If the request crosses
   a global constraint (e.g. commandment 8: "Rust. SQLite. Local.
   Deterministic. No cloud except PubMed API."), surface it and get
   explicit sign-off before adding the dependency. Don't rationalize.

3. **Match existing conventions.** Read one or two sibling commands
   under `crates/atman/src/commands/` before writing a new one. Keep
   argument naming, TSV column naming, and error-handling style
   consistent. Pure-Rust, pure-deterministic, TSV in / TSV out.

4. **Put generic algorithms in `atman-core`.** CLI wiring, file I/O, and
   path handling live in `crates/atman/src/commands/<name>.rs`; the
   algorithmic core (anything that could be unit-tested without a
   filesystem) goes in `crates/atman-core/src/<name>.rs` and is
   re-exported from `lib.rs`.

5. **Determinism is a correctness requirement.** Seed every randomized
   routine from a user-visible `--seed` (default `20260418`). Use the
   in-tree Xoshiro256++ where a PRNG is needed so byte-for-byte
   reproducibility holds across runs and across machines.

6. **Write tests at two levels.**
   - Unit tests inside the core module (`#[cfg(test)] mod tests`).
   - One integration test per command under `crates/atman/tests/
     <command>.rs`, invoking the built binary via
     `env!("CARGO_BIN_EXE_atman")` against a tempdir fixture.

7. **Run the full workspace before committing.**
   `cargo test --workspace --release` must pass with zero failures.
   No skipped tests, no `#[ignore]` without explicit justification.

8. **Smoke-test on real data when possible.** For commands that consume
   canonical Atman TSVs, run the new command against the Dube fixture
   (`example_data/dube_heat_2023`) after `ingest` + `qc`. This catches
   schema assumptions the synthetic test fixture didn't.

9. **Update user-facing docs in the same commit.**
   - `README.md` commands list
   - `CHANGELOG.md` under the `## [Unreleased]` heading
   - Any relevant section in `docs/` (tutorial, analytical-roadmap)

10. **Clear the shipped request from this file** — do not mark it as
    shipped and leave it in place. Delete it. Keep this file reflecting
    only what is still open.

11. **Commit cleanly.** One feature per commit. Tools don't get
    authorship — no `Co-Authored-By` trailers (commandment 9). Imperative
    subject line under 70 chars, body explains *why* in 2–4 sentences,
    end with `Closes <priority> feature <n> …` if applicable.

12. **Validate before moving on.** Before starting the next request,
    confirm the last one is (a) tested, (b) committed, (c) removed from
    this file. Only then advance.

## Notes for atman maintainers

- The CSF_CrossDisease manuscript is the first external consumer of atman
  v1.0.0 beyond the Dube NPX reproduction path. Any feature work here
  benefits other proteomics-analysis consumers of the canonical TSV
  contract.
- Contact: see `/Users/kevinjoseph/Cursor/CSF_CrossDisease/manuscript/`
  for the downstream analytic context each feature targets.
- Seed: all CSF manuscript analyses use `20260418` for direct
  compatibility with existing Python output; keeping that as the default
  seed in bootstrap/ICA/null commands simplifies cross-validation.
