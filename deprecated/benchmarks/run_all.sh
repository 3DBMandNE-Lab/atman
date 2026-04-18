#!/usr/bin/env bash
#
# Drive all three comparators on the Dube dataset and aggregate into Table 1.
# Each per-tool script is self-contained; failure in one does not abort the
# others (benchmark output reports per-tool status).

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

OUT="benchmarks/out"
mkdir -p "$OUT"

echo "[bench 1/4] atman baseline"
bash benchmarks/atman/run.sh   2>&1 | tee "$OUT/atman.log"   || echo "atman harness failed"

echo "[bench 2/4] OlinkAnalyze"
Rscript benchmarks/olinkanalyze/run.R  2>&1 | tee "$OUT/olinkanalyze.log"  || echo "OlinkAnalyze harness failed"

echo "[bench 3/4] limma"
Rscript benchmarks/limma/run.R         2>&1 | tee "$OUT/limma.log"         || echo "limma harness failed"

echo "[bench 4/4] aggregate"
python3 benchmarks/aggregate.py        2>&1 | tee "$OUT/aggregate.log"     || echo "aggregator failed"

echo ""
echo "Table 1: $OUT/table1.tsv"
echo "UpSet:   $OUT/upset_de_overlap.pdf"
