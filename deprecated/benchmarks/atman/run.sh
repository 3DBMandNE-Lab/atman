#!/usr/bin/env bash
# Atman baseline harness.
# Relies on ./out/ being already populated by scripts/reproduce_paper.sh.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"
export PATH="$REPO_ROOT/target/release:$PATH"

OUT_BENCH="benchmarks/out"
mkdir -p "$OUT_BENCH"

[[ -f out/de_results.tsv ]] || {
    echo "out/de_results.tsv missing. Run scripts/reproduce_paper.sh first." >&2
    exit 1
}

START=$(python3 -c 'import time; print(time.time())')

# Re-run DE in isolation to measure atman's own timing.
atman de --input-dir out --output-dir benchmarks/out/atman \
    --test paired-t --paired-by participant \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" --min-pairs 5

END=$(python3 -c 'import time; print(time.time())')
RUNTIME=$(python3 -c "print(round($END - $START, 3))")

python3 - "$RUNTIME" <<'PY'
import sys, pandas as pd
runtime = float(sys.argv[1])
src = pd.read_csv("benchmarks/out/atman/de_results.tsv", sep="\t")
out = pd.DataFrame({
    "contrast":    src["comparison"],
    "gene":        src["gene_symbol"],
    "mean_diff":   src["mean_diff"],
    "pvalue":      src["p_value"],
    "qvalue":      src["bh_q"],
    "tool":        "atman",
    "runtime_s":   runtime,
    "peak_rss_mb": pd.NA,
})
out.to_csv("benchmarks/out/atman_de.tsv", sep="\t", index=False)
print(f"atman DE rows: {len(out):,}")
PY
