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
