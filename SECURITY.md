# Security Policy

## Reporting a vulnerability

Please report suspected vulnerabilities privately via the repository's GitHub
Security Advisories ("Report a vulnerability"), or to the maintainer listed in
`CITATION.cff`. Please do not open a public issue for a security report.

## Threat model

`atman` is a **local, single-user command-line analysis tool**, not a network
service. It reads scientific data files (TSV/CSV/JSON), writes result files, and
makes at most one outbound network call (see below). Calibrate severity
accordingly: a panic on malformed input is a robustness bug, not a remote
exploit.

- **`unsafe` is confined to three FFI sites, on macOS only.** The workspace
  contained zero `unsafe` blocks until the NMF and ICA kernels were routed
  through Accelerate's BLAS. There are now exactly three, all in
  `crates/atman-core/src/blas.rs`, all behind `#[cfg(target_os = "macos")]`,
  and all reachable only through `gemm` and `syrk`:

  | Site | What it does | Why it is sound |
  |---|---|---|
  | `set_single_threaded_once` | `dlsym` + `transmute` to call `BLASSetThreading` | Null-checked before use; the transmuted signature matches Apple's `vecLib/thread_api.h`; `Once` bounds it to one call per process |
  | `gemm_accelerate` | `cblas_dgemm` | Private, reachable only via `gemm`, which asserts every slice length against its dimensions |
  | `syrk` | `cblas_dsyrk` | Same, via `syrk`'s own assertions |

  Each call site carries a `SAFETY` comment naming the invariant it depends
  on. `assert_fits_c_int` guards the one way those length checks could be
  defeated: CBLAS takes 32-bit signed dimensions, so a `usize` above
  `i32::MAX` would truncate and the library would read outside the slices
  whose lengths were validated in `usize`. atman's matrices are far below
  that bound, so the assertion never fires in practice; it exists because
  the alternative to firing is memory unsafety rather than a wrong number.

  No other platform compiles any `unsafe`: every other target uses the pure
  Rust `gemm_fallback`, which has none. Accelerate is a system framework, so
  this adds no dependency to the SBOM.

- **No manual memory management** outside those three sites.
- **Untrusted input (data files).** TSV/CSV/JSON readers are header-indexed and
  reject ragged rows (the `csv` reader runs in strict, non-flexible mode), so a
  malformed file produces a descriptive error rather than a crash or
  out-of-bounds access. Matrix dimensions are derived from data actually read,
  never from an integer claimed in a header, so a small malicious file cannot
  induce a large allocation. A non-finite (`NaN`) value in an abundance column
  is **intentional** — it is the canonical MNAR / below-LOD missingness marker
  that missingness-aware methods model; finite-assuming consumers go through
  `MeasurementRecord::effective_abundance()`, which excludes it.
- **Output paths.** `--output` / `--output-dir` are operator-supplied. Where a
  filename component is derived from *data* (the `panel` column drives the
  per-panel CSV names in `atman matrix` / `atman fold-change`), that component is
  sanitized and rejected if it contains a path separator or `..` (1.1.0+), so a
  crafted data file cannot write outside the intended directory. Writes are
  atomic (tempfile + rename), with no predictable-name temp race.

## Operator-controlled functionality (trust boundaries)

Two commands execute or load operator-supplied programs. Treat their inputs with
the same trust you would a shell script:

- **`atman run`** executes each stage of a pre-registered plan file via `sh -c`,
  as the invoking user, with full shell semantics. **Plan files must be
  operator-authored and trusted.** Do not run a plan from an untrusted source.
- **`atman bench`** runs adapter scripts from `--adapters-dir` (arguments are
  passed as separate argv entries, not through a shell). The adapters directory
  must be operator-controlled.

## Network surface

- The only outbound network call is **`atman enrich gprofiler`** (the g:Profiler
  REST API). It uses HTTPS with TLS certificate validation (rustls +
  webpki-roots), bounds the response size, and JSON-escapes all request fields.
- **`--offline`** disables the network entirely (cache hit required, else it
  fails). The on-disk response cache (`--cache-dir`) is keyed by a SHA-256 of the
  canonical request and stores raw upstream JSON; **treat `--cache-dir` as
  trusted** and do not share it across trust boundaries.

## Dependencies & supply chain

- **Licenses:** all transitive dependencies are permissively licensed
  (MIT / Apache-2.0 / BSD / ISC / Zlib / Unicode-3.0 and equivalents), compatible
  with atman's `MIT OR Apache-2.0`. No copyleft-only dependencies.
  - **`docs/sbom-1.2.0.cdx.json`** (CycloneDX 1.5, `cargo-cyclonedx` 0.5.9)
  describes the released `aarch64-apple-darwin` binary: its normal and build
  dependencies, **129** crates. The lockfile holds **187** entries. The
  **57** it lists and the SBOM does not fall into three classes, each
  checked with `cargo tree`: **21** other-platform crates (`windows-*`,
  `wasi*`, `linux-raw-sys` and their pins), **14** dev-only crates
  (`proptest` and its tree), and **22** optional dependencies of
  transitive crates that no feature in this workspace turns on
  (`wasm-bindgen`, `wit-bindgen` and their trees). Cargo locks optional
  dependencies whether or not they are activated. None of the 57 is
  compiled into the released binary. The two workspace crates carry the
  public git source as their identity, not the build machine's path.
- `docs/sbom-1.1.0.cdx.json` is the previous release's record. It listed
  every lockfile entry, including the three classes above, so the two
  files are not comparable by count.
- **Advisories (`cargo audit`):** `cargo audit` runs in CI (the `audit` job) and
  scans `Cargo.lock` against the RustSec database on every build.
  - **Fixed in 1.2.0:** `RUSTSEC-2026-0190` (unsoundness in
    `anyhow::Error::downcast_mut()`) — resolved by updating `anyhow` from
    1.0.102 to 1.0.104, which is past the advisory's `>= 1.0.103` fix. atman
    never calls `downcast_mut` (it calls `downcast_ref` once, in
    `axes/contrast.rs`), so the advisory was not reachable — but a patch
    existed, and taking it is better than carrying a reachability argument.
  - **Fixed in 1.1.0:** `RUSTSEC-2026-0104` (reachable panic in
    `rustls-webpki` CRL parsing, on the `enrich gprofiler` TLS path) — resolved
    by updating `rustls-webpki` to 0.103.13.
  - **Scoped ignores** (in `.cargo/audit.toml`) — two advisories that are
    transitive via `statrs → nalgebra` and **not reachable from atman's code
    paths**:
    - `RUSTSEC-2026-0097` — `rand 0.8.5` "unsound with a custom logger using
      `rand::rng()`". atman never calls `rand`'s sampling APIs (it ships its own
      deterministic `atman_core::rng`); `rand` is pulled in only by `statrs`'s
      distribution machinery, which atman uses for scalar CDFs/quantiles only.
    - `RUSTSEC-2024-0436` — `paste` unmaintained (transitive via `nalgebra`).
    - Both are **currently inert**: `rand 0.8.5` and `paste 1.0.15` are still in
      the lockfile, but neither advisory fires as of 1.2.0 and `cargo audit`
      reports zero warnings even under `--deny warnings`. The entries are kept
      rather than deleted, because a suppression that stops matching is not
      evidence the risk is gone. Re-check at each release.
  The ignores are scoped to those specific IDs, so any **new** advisory still
  fails CI. They stem from the unused `statrs → nalgebra/rand` weight. `statrs`
  is used only for closed-form CDFs, quantiles and special functions, never
  for matrix algebra or sampling, so the `nalgebra` and `rand` it pulls in
  serve no code path here (symbol audit, 2026-05-29). Dropping or replacing
  `statrs` would let us remove both ignores.
- **Per release:** regenerate the CycloneDX SBOM (`docs/sbom-<version>.cdx.json`).
