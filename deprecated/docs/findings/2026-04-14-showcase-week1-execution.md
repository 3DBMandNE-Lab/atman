# Atman Showcase Week 1 Execution (2026-04-14)

## Run Scope

Executed from raw NPX inputs:

1. `atman ingest`
2. `atman qc`
3. `atman matrix`
4. `atman fold-change`
5. `atman de`

Output workspace: `/tmp/karna_showcase`

## Pipeline Outcome

1. Ingest: `122,135` rows, `2,943` assays, `43` samples.
2. QC: `3,725` masked, `118,410` passed.
3. DE: `11,752` rows, `11,752` computed, `0` skipped (`min_pairs=5`).

## DE Discovery Counts (BH q < 0.05)

1. `PT1-PR1`: `36`
2. `PT2-PR2`: `37`
3. `PT2-PT1`: `7`
4. `PR2-PR1`: `2`

Interpretation: strongest signal remains in acute heat contrasts (`PT1-PR1`, `PT2-PR2`).

## Canonical Heat-Shock Sanity Check

Observed positive acute effects in canonical heat-stress proteins:

1. `HSPB1`: positive in `PT1-PR1` and `PT2-PR2`.
2. `HSPA1A`: positive in `PT1-PR1` and `PT2-PR2`.
3. `DNAJB1`: positive in `PT2-PR2`.
4. `HSPG2`: positive in `PT2-PR2`.

These are directionally coherent with the acute heat response.

## Top Acute Signals by q-value

`PT1-PR1` includes strong positive changes such as `CALCA`, `FGFBP3`, `ANGPTL4`, and `PRL`.

`PT2-PR2` includes strong positive changes such as `EDN1`, `CXCL17`, `TIMP4`, and `SMOC1` (with some negative shifts including `MOG`, `PTPRR`).

## Week 1 Conclusion

1. Reproduction and core DE stack are operational from raw files.
2. Acute heat proteomic response is robust and biologically coherent.
3. Baseline is ready for Week 2 robustness profiling without changing public CLI surface.

