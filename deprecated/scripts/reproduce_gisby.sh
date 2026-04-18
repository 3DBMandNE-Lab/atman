#!/usr/bin/env bash
# Reproduce the Gisby 2021 Olink Target 96 generalizability demonstration.
# Assumes `atman` on PATH (run `cargo install --path crates/atman` first)
# and Python 3.11+ with pandas.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
export PATH="$REPO_ROOT/target/release:$PATH"

DATA=example_data/gisby_covid_2021
OUT=out_gisby

# Fetch raw data if missing.
if [[ ! -f "$DATA/plasma_npx_level.csv" ]]; then
    mkdir -p "$DATA"
    for f in plasma_npx_level.csv plasma_sample_level.csv serum_npx_level.csv serum_sample_level.csv; do
        curl -sL -o "$DATA/$f" \
            "https://raw.githubusercontent.com/jackgisby/longitudinal_olink_proteomics/master/data/$f"
    done
fi

echo "[1/4] ingest (Python adapter)"
python3 scripts/ingest_gisby.py "$DATA" "$OUT"

echo "[2/4] fold-change"
atman fold-change --input-dir "$OUT" --output-dir "$OUT" --groups "late-early"

echo "[3/4] de (paired-t, late - early)"
atman de --input-dir "$OUT" --output-dir "$OUT" \
    --test paired-t --paired-by participant --groups "late-early" --min-pairs 5

echo "[4/4] LOO robustness (one DE run per excluded patient)"
rm -rf out_gisby_loo
mkdir -p out_gisby_loo
PATIENTS=$(python3 -c "import pandas as pd; s=pd.read_csv('$OUT/samples.tsv', sep='\t'); print(' '.join(sorted(s.subject_id.unique())))")
LOO_PATHS=""
for pid in $PATIENTS; do
    d="out_gisby_loo/loo_${pid}"
    mkdir -p "$d"
    python3 - "$pid" "$OUT" "$d" <<'PY'
import sys, pandas as pd
pid, src, dst = sys.argv[1], sys.argv[2], sys.argv[3]
samples = pd.read_csv(f'{src}/samples.tsv', sep='\t')
bad = set(samples.loc[samples.subject_id==pid, 'sample_id'])
for f in ('measurements','qc_measurements','samples','proteins'):
    d = pd.read_csv(f'{src}/{f}.tsv', sep='\t')
    if 'subject_id' in d.columns:
        d = d[d.subject_id != pid]
    elif 'sample_id' in d.columns:
        d = d[~d.sample_id.isin(bad)]
    d.to_csv(f'{dst}/{f}.tsv', sep='\t', index=False)
PY
    atman de --input-dir "$d" --output-dir "$d" \
        --test paired-t --paired-by participant \
        --groups "late-early" --min-pairs 5 >/dev/null
    LOO_PATHS="${LOO_PATHS}${d}/de_results.tsv,"
done
LOO_PATHS=${LOO_PATHS%,}

atman robustness --baseline "$OUT/de_results.tsv" --loo "$LOO_PATHS" \
    --top-k 20 --output-dir "$OUT/robustness"

echo ""
echo "Done. Summary:"
python3 - <<'PY'
import pandas as pd
de = pd.read_csv("out_gisby/de_results.tsv", sep="\t")
print(f"  DE hits q<0.05: {(de.bh_q<0.05).sum()} / {len(de)}")
print(f"  DE hits q<0.10: {(de.bh_q<0.10).sum()} / {len(de)}")
sign = pd.read_csv("out_gisby/robustness/loo_sign_stability.tsv", sep="\t")
rank = pd.read_csv("out_gisby/robustness/rank_stability.tsv", sep="\t")
print(f"  Mean sign-match rate (LOO): {sign['sign_match_rate'].mean():.4f}")
print(f"  Mean top-20 Jaccard (LOO):  {rank['jaccard_index'].mean():.4f}")
PY
