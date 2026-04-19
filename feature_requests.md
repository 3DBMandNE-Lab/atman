# Atman feature requests

No open requests. All filed items from the 2026-04-19 CSF cross-disease
batch shipped the same day.

## How to implement a feature request

Work one feature end-to-end before starting the next.

1. **Read the entire request.** Command sketch, "what it does", "why we
   need it", and scope hints. Treat the CLI contract as load-bearing;
   the internal algorithm is a choice.

2. **Surface constraint conflicts first.** If the request touches a
   global constraint (e.g. commandment 8: "Rust. SQLite. Local.
   Deterministic. No cloud except PubMed API."), raise it and get
   explicit sign-off before adding the dependency.

3. **Match existing conventions.** Skim one or two sibling commands in
   `crates/atman/src/commands/` first, then mirror their argument
   naming, TSV column naming, and error-handling style. Keep
   everything pure Rust, deterministic, and TSV-in / TSV-out.

4. **Split algorithm and CLI.** Put generic algorithms in
   `crates/atman-core/src/<name>.rs` and re-export from `lib.rs`.
   Keep CLI wiring, file I/O, and path handling in
   `crates/atman/src/commands/<name>.rs`.

5. **Reuse shared helpers.** Before writing any utility, grep
   `atman-core::stats` (`mean`, `ranks`, `pearson`, `spearman`,
   `cosine`, `jaccard_top_n`, `top_abs_indices`) and `atman::io`
   (`need_col`, `find_col`, `optional_cell`, `format_float`,
   `escape_tsv`, `sha256_hex`, `atomic_write`). Extend the shared
   module when something is genuinely missing; add the local copy
   only when the semantics are truly command-specific.

6. **Seed every randomized routine.** Expose `--seed` (default
   `20260418`) and draw from the in-tree Xoshiro256++ so runs are
   byte-for-byte reproducible across machines.

7. **Write tests at two levels.**
   - Unit tests inside the core module (`#[cfg(test)] mod tests`).
   - One integration test per command under `crates/atman/tests/
     <command>.rs`, invoking the built binary via
     `env!("CARGO_BIN_EXE_atman")` against a tempdir fixture.

8. **Run the full workspace before committing.**
   `cargo test --workspace --release` should finish green with every
   test executed.

9. **Smoke-test on real data.** For commands that consume canonical
   Atman TSVs, run the new command against the Dube fixture
   (`example_data/dube_heat_2023`) after `ingest` + `qc`. This catches
   schema assumptions the synthetic test fixture misses.

10. **Update user-facing docs in the same commit.**
    - `README.md` commands list
    - `CHANGELOG.md` under the `## [Unreleased]` heading
    - Any relevant section in `docs/` (tutorial, analytical-roadmap)

11. **Delete shipped requests from this file.** Keep this file
    reflecting only what is still open.

12. **Commit cleanly.** One feature per commit. Commit as yourself
    alone (commandment 9). Imperative subject line under 70 chars,
    body explains *why* in 2–4 sentences, end with `Closes <priority>
    feature <n> …` when applicable.

13. **Validate before advancing.** Confirm the last feature is (a)
    tested, (b) committed, (c) removed from this file. Then pick up
    the next.

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
