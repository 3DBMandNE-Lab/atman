# `atman bench decompose` — cross-tool benchmark harness

Score `atman decompose ica` against named reference tools on a shared
synthetic planted-archetype fixture. Produces a canonical
`bench_results.tsv` with one row per `(tool, planted_archetype)`.

## Fixture layout

Each fixture directory under `bench/` is self-contained:

- `planted_loadings.tsv` — columns `archetype_id`, `protein`,
  `loading`. Ground-truth loading vectors, one row per
  (archetype, protein) cell.
- `abundance.tsv` — first column `sample_id`, remaining columns one
  protein each. Synthetic abundance matrix generated from the
  planted loadings plus calibrated noise.
- `cohorts.tsv` *(optional)* — columns `sample_id`, `cohort`. Used
  by future cross-cohort benchmark extensions; single-cohort fixtures
  omit it.

The bundled `planted_archetypes_v1/` fixture plants 2 archetypes on
20 proteins across 20 subjects. Archetype 1 loads on proteins
`G001..G010`; archetype 2 loads on `G011..G020`; proteins
`G021..G030` are always zero. Subject activations are cubed i.i.d.
Gaussians (deterministic via LCG seeds 77, 88, 99), which breaks
Gaussianity so FastICA can recover the sources.

## Running atman natively

```bash
atman bench decompose \
  --fixture bench/planted_archetypes_v1 \
  --tools atman \
  --seed 20260420 \
  --top-n 10 \
  --output bench/atman_only_results.tsv
```

Atman runs FastICA (`fast_ica` + `canonicalize_ica` from
`atman-core`) on the fixture's abundance matrix and scores the
recovered loading matrix against the planted one via
`atman_core::bench_decompose::score_against_planted`:

- `archetype_correlation = |pearson(planted_row, recovered_row)|`
- `recovery_jaccard` = Jaccard over each side's top-N loading sets
  by `|loading|`
- Matching is reciprocal-best absolute cosine; archetypes with no
  reciprocal-best partner score zeros.
- `determinism_score = 1.0` when two repeated runs at the same seed
  yield byte-equal recovered loadings (the default).

## Adapter contract for external tools

Tools other than `atman` are invoked via shell adapters:

```bash
atman bench decompose \
  --fixture bench/planted_archetypes_v1 \
  --tools atman,fastica-icasso,mofa \
  --adapters-dir bench/adapters \
  --output bench/all_tools.tsv
```

Each adapter script lives at `bench/adapters/<tool>.sh` and must:

1. Accept **two positional arguments**: the fixture directory and an
   output directory that it must fill.
2. Read `ATMAN_BENCH_SEED` and `ATMAN_BENCH_K` from the environment
   if the underlying tool exposes them; otherwise ignore.
3. Write `recovered_loadings.tsv` into the output directory, with
   the **same three-column contract as `planted_loadings.tsv`**
   (`archetype_id`, `protein`, `loading`). No other files are
   required.
4. Exit 0 on success, non-zero on failure.

Missing adapter scripts or non-zero exits produce a single
`tool_not_available = 1` row in the results TSV rather than
aborting the harness. This keeps `atman bench decompose --tools
atman,fastica-icasso` usable on machines where `fastica-icasso`
isn't installed — the atman row still lands.

Adapters themselves live outside atman's analytical path
(commandment 8): they are thin wrappers around whichever R or
Python implementation the reference tool provides, stitched to the
canonical TSV-in / TSV-out contract. We don't ship any here —
writing one is ~20 lines per tool and stays on the benchmark
operator's machine.

## Two ways to measure the wrong thing

Both of these produced a published-looking number that was wrong, on this
project, in 2026-09.

### A microbenchmark measures the kernel, not the kernel as used

`cblas_dgemm` at NMF's shapes (n=110, p=10000, k=10) ran at 0.179 ms and
0.304 ms in a tight 200-repeat loop. The same calls inside the solve ran
at 0.743 ms and 0.864 ms — three to four times slower.

Back-to-back calls keep Accelerate's worker threads hot. Calls separated
by other work pay to wake them again. The isolated benchmark could not
see that, and it was the whole remaining gap to scikit-learn.

Time the phase inside the real loop before believing a kernel number.
Instrumenting the loop found this after two rounds of estimating did not.

### Check that the thing you measured actually ran

Two measurements in this project were invalid because the command under
test did not do what the timing implied, and in both cases the number
looked right.

One compared a release binary against itself. Another session had edited
the source in the same working tree, `cargo build` reported
`Finished in 0.15s` — no recompilation — and a golden-output comparison
across that build returned BYTE-IDENTICAL on every file. It was one
binary compared to itself.

The other timed a command that was failing. `--k-selection fixed=30` is
not a supported rule, and validation used to happen after the input was
read, so the run cost 3.5 s on a 123 MB cohort and wrote nothing. Four
of those at different `--max-iter` produced a flat slope, which was read
as evidence that FastICA iterations are free. (They are free — confirmed
later on runs that wrote their output — but the measurement behind the
claim was of a command that did nothing.)

Both preconditions are one line:

- did the build actually recompile? Look for `Compiling`, not `Finished`.
- did the run actually write its output? Check the file, or the exit
  code, or a line count.

**In both cases the wrong answer was the one we were hoping for.** One
confirmed a fix had worked; the other confirmed a hypothesis formed
minutes earlier. A check that fails silently is dangerous in proportion
to how much you want its result, and a silent failure that lands on the
side you expected will not feel like one.

Validation that runs after the expensive load has this shape generally:
a failure then costs about what success costs, and becomes timeable. In
`decompose ica` that flag is now rejected in 0.01 s regardless of input
size, so a timing containing a failure is obviously wrong.

### An adapter measures the interpreter, not the tool

`bench/adapters/*.sh` wrap external tools as subprocesses, so an adapter
timing includes interpreter startup, TSV parsing and TSV writing. On a
110 x 10000 fixture that overhead was 2.31 s of a 2.436 s measurement.

Quoted from adapter timings, atman looked 12.9x faster than scikit-learn
on ICA and 7.3x slower on NMF. Timed kernel to kernel, it was 1.4x slower
on ICA and 16x slower on NMF. Both published figures were artifacts, in
opposite directions.

Compare kernels to kernels. If you quote an adapter number, say so.

### And check the iteration count before optimising

A wall-clock comparison between two solvers is meaningless without their
iteration counts. `decompose ica` records `ica_n_iterations`,
`ica_converged` and `ica_final_tol` in its run sidecar; `bench decompose`
does not carry them in the row. FastICA's iterations run on the
k-dimensional whitened data and are effectively free — 2000 of them cost
the same wall clock as 10 — so for ICA the kernel is the whitening, not
the loop.

## Building a non-negative fixture for NMF

Make it non-negative BY CONSTRUCTION (`W @ H` with `W, H >= 0`), not by
shifting (`X - X.min() + eps`). A shift adds a constant offset that NMF
then has to model, and on such a fixture atman scores correlation 0.3155
and recovers 4 of 10 archetypes. It looks broken and is not.
