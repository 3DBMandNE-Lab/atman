# atman fuzzing

Coverage-guided fuzzing of atman's untrusted-input parsers, via
[`cargo-fuzz`](https://github.com/rust-fuzz/cargo-fuzz) (libFuzzer). The targets
exercise the canonical TSV readers on arbitrary bytes; the invariant under test
is **no panic / OOB / hang — only `Ok` or `Err`**.

This is the deeper complement to the stable-toolchain
`crates/atman/tests/adversarial_input.rs` regression tests.

## Requirements

- A **nightly** Rust toolchain (libFuzzer needs nightly): `rustup toolchain add nightly`
- `cargo install cargo-fuzz`

## Run

```bash
# from the repo root
cargo +nightly fuzz run read_measurements
cargo +nightly fuzz run read_samples
cargo +nightly fuzz run read_proteins

# time-boxed (e.g. CI): 60 seconds
cargo +nightly fuzz run read_measurements -- -max_total_time=60
```

## Targets

| Target | Parser under test |
|---|---|
| `read_measurements` | `atman::io::read_measurements_long` |
| `read_samples` | `atman::io::read_samples` |
| `read_proteins` | `atman::io::read_proteins` |

A crash writes a reproducer to `fuzz/artifacts/<target>/`; re-run it with
`cargo +nightly fuzz run <target> fuzz/artifacts/<target>/<id>`.

This crate is its own Cargo workspace and is excluded from the stable
`cargo build --workspace` / `cargo test --workspace`.
