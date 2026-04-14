#!/usr/bin/env bash
# Pathway enrichment pipeline for de_results.tsv → per-(comparison,
# direction, collection) enrichment JSON files. Shells out to the
# external `genesets` Rust CLI from ../Rust_bioinfo_cli/genesets.
#
# Usage:
#   scripts/run_enrichment.sh <de_results.tsv> <output-dir>
#
# Or set KARNA_GENESETS=/path/to/genesets to override the default
# location. Expects MSigDB Hallmark / KEGG / GO-BP gene sets already
# bundled in the genesets binary.
#
# Foreground lists are capped at top-100 genes per (comparison, dir)
# ranked by |mean_diff| among q < Q_THRESHOLD (default 0.10). The cap
# works around a panic in genesets on very large input lists and also
# keeps the enrichment from being dominated by a handful of large sets.

set -euo pipefail

DE_RESULTS="${1:-}"
OUTDIR="${2:-}"
Q_THRESHOLD="${Q_THRESHOLD:-0.10}"
TOP_N="${TOP_N:-100}"
GENESETS="${KARNA_GENESETS:-/Users/kevinjoseph/Cursor/Rust_bioinfo_cli/genesets/target/release/genesets}"

if [[ -z "$DE_RESULTS" || -z "$OUTDIR" ]]; then
    echo "usage: $0 <de_results.tsv> <output-dir>" >&2
    exit 2
fi
if [[ ! -f "$DE_RESULTS" ]]; then
    echo "de_results.tsv not found: $DE_RESULTS" >&2
    exit 2
fi
if [[ ! -x "$GENESETS" ]]; then
    echo "genesets binary not found at $GENESETS" >&2
    echo "set KARNA_GENESETS=/path/to/genesets to override" >&2
    exit 2
fi

mkdir -p "$OUTDIR"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Discover comparisons present in the file.
COMPARISONS=$(awk -F'\t' 'NR>1 {print $5}' "$DE_RESULTS" | sort -u)

for cmp in $COMPARISONS; do
    for dir in up down; do
        if [[ "$dir" = "up" ]]; then
            sign_expr='$9+0>0'
            sort_key='$9'
        else
            sign_expr='$9+0<0'
            sort_key='(-$9)'
        fi

        LIST="$WORK/${cmp}_${dir}.txt"
        awk -F'\t' -v c="$cmp" -v q="$Q_THRESHOLD" \
            "NR>1 && \$5==c && \$13!=\"\" && \$13+0<q && $sign_expr {print $sort_key\"\\t\"\$3}" \
            "$DE_RESULTS" \
            | sort -k1,1 -gr \
            | head -n "$TOP_N" \
            | cut -f2 \
            | sort -u \
            > "$LIST"

        n=$(wc -l < "$LIST" | tr -d ' ')
        if [[ "$n" -eq 0 ]]; then
            continue
        fi

        for coll in hallmark kegg go_bp; do
            out="$OUTDIR/${cmp}_${dir}_${coll}.json"
            if "$GENESETS" enrich --input "$LIST" --collection "$coll" --json > "$out" 2>/dev/null; then
                printf "  %-10s %-4s %-8s n=%-3d → %s\n" "$cmp" "$dir" "$coll" "$n" "$(basename "$out")"
            else
                # genesets panics on some inputs; record the failure and move on.
                rm -f "$out"
                printf "  %-10s %-4s %-8s n=%-3d → FAILED (likely NaN in Fisher)\n" \
                    "$cmp" "$dir" "$coll" "$n"
            fi
        done
    done
done

echo
echo "done. JSON outputs in $OUTDIR"
