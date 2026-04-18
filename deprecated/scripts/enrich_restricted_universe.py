#!/usr/bin/env python3
"""Universe-corrected pathway enrichment for Atman de_results.tsv.

genesets' built-in `enrich` subcommand uses the full MSigDB gene universe as
the background, but our measurement universe is only the ~2,943 proteins on
Olink Explore 1536 + Explore Expansion. This over-counts enrichment for
pathways with high Olink coverage.

This script:

    1. Loads gene set definitions from genesets.json.gz (all collections).
    2. Reads de_results.tsv and collects the set of gene symbols actually
       tested — this is our universe.
    3. For each (comparison, direction):
         foreground = top-100 genes at q<0.10, ranked by |mean_diff|
       For each gene set:
         set_in_universe   = pathway gene list ∩ universe
         hits              = foreground ∩ set_in_universe
         Fisher's exact 2x2 on:
             [hits,                    |foreground| - hits                 ]
             [|set_in_universe| - hits, |universe| - |foreground| - |set_in_universe| + hits]
       BH-FDR within the (comparison, direction, collection) family.
    4. Emits one TSV per (comparison, direction, collection) to the output
       directory, plus a combined summary.

Usage:
    scripts/enrich_restricted_universe.py \
        <de_results.tsv> \
        <output_dir> \
        [--genesets /path/to/genesets.json.gz]

Defaults:
    --genesets: ~/Cursor/Rust_bioinfo_cli/genesets/data/genesets.json.gz
    --q-threshold: 0.10
    --top-n: 100
    --min-set-size: 3
"""

import argparse
import csv
import gzip
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path


def fisher_exact_right_tail(a: int, b: int, c: int, d: int) -> float:
    """Right-tailed Fisher's exact p-value for the 2x2 contingency table:

        |       | in_set        | not_in_set     |
        | fg    | a             | b              |
        | bg    | c             | d              |

    Uses the hypergeometric distribution:
        P(X >= a) = sum_{k=a}^{min(a+b, a+c)} hypergeom(k; N, K, n)
    where N = a+b+c+d, K = a+c (set size in universe), n = a+b (foreground size).
    """
    N = a + b + c + d
    K = a + c
    n = a + b
    if K == 0 or n == 0 or K > N or n > N:
        return 1.0
    k_max = min(n, K)
    log_nCk = _log_binom_table(N)
    log_total = log_nCk(N, n)
    p = 0.0
    for k in range(a, k_max + 1):
        log_p_k = log_nCk(K, k) + log_nCk(N - K, n - k) - log_total
        p += math.exp(log_p_k)
    return min(max(p, 0.0), 1.0)


def _log_binom_table(max_n: int):
    """Return a function log_binom(n, k) using a shared log-factorial cache."""
    # Grow on demand.
    log_fact: list[float] = [0.0, 0.0]

    def ensure(n: int) -> None:
        while len(log_fact) <= n:
            log_fact.append(log_fact[-1] + math.log(len(log_fact)))

    def log_binom(n: int, k: int) -> float:
        if k < 0 or k > n:
            return -math.inf
        ensure(n)
        return log_fact[n] - log_fact[k] - log_fact[n - k]

    ensure(max_n)
    return log_binom


def bh_fdr(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction. Returns q-values aligned with input."""
    m = len(p_values)
    if m == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    # Compute sorted q-values and apply reverse monotonicity sweep.
    q_sorted = [p_values[o] * m / (r + 1) for r, (o, _) in enumerate(indexed)]
    q_sorted = [min(1.0, v) for v in q_sorted]
    for i in range(len(q_sorted) - 2, -1, -1):
        if q_sorted[i + 1] < q_sorted[i]:
            q_sorted[i] = q_sorted[i + 1]
    out = [0.0] * m
    for rank, (orig, _) in enumerate(indexed):
        out[orig] = q_sorted[rank]
    return out


def load_genesets(path: Path) -> list[dict]:
    """Load gene set definitions from genesets.json.gz."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as f:
            data = json.load(f)
    else:
        with open(path) as f:
            data = json.load(f)
    return data["sets"]


def load_de(path: Path):
    """Return (universe, hits_by_comp_dir) where hits_by_comp_dir maps
    (comparison, direction) → list of top-N gene symbols."""
    universe = set()
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row["gene_symbol"]
            if gene:
                universe.add(gene)
            rows.append(row)
    return universe, rows


def top_hits(rows, comparison: str, direction: str, q_threshold: float, top_n: int) -> list[str]:
    """Extract top-N foreground gene symbols for a given (comparison, direction)."""
    filtered = []
    for r in rows:
        if r["comparison"] != comparison:
            continue
        if not r["bh_q"] or not r["mean_diff"]:
            continue
        q = float(r["bh_q"])
        if q >= q_threshold:
            continue
        md = float(r["mean_diff"])
        if direction == "up" and md <= 0:
            continue
        if direction == "down" and md >= 0:
            continue
        filtered.append((abs(md), r["gene_symbol"]))
    filtered.sort(reverse=True)
    seen = set()
    out = []
    for _, gene in filtered:
        if gene not in seen:
            seen.add(gene)
            out.append(gene)
            if len(out) >= top_n:
                break
    return out


def enrich_one(
    foreground: list[str],
    universe: set[str],
    gene_sets: list[dict],
    min_set_size: int,
) -> list[dict]:
    """Fisher's exact on each set, restricted to `universe`. Returns list of dicts."""
    fg = set(foreground) & universe
    universe_size = len(universe)
    fg_size = len(fg)
    results = []
    for gs in gene_sets:
        set_in_univ = set(gs["genes"]) & universe
        set_size = len(set_in_univ)
        if set_size < min_set_size:
            continue
        hits = fg & set_in_univ
        a = len(hits)
        b = fg_size - a
        c = set_size - a
        d = universe_size - fg_size - c
        if d < 0:
            continue
        p = fisher_exact_right_tail(a, b, c, d)
        results.append({
            "name": gs["name"],
            "collection": gs["collection"],
            "set_size_raw": len(gs["genes"]),
            "set_size_universe": set_size,
            "foreground_size": fg_size,
            "universe_size": universe_size,
            "hits": a,
            "pvalue": p,
            "hit_genes": sorted(hits),
        })
    # BH-FDR within this family.
    p_values = [r["pvalue"] for r in results]
    q_values = bh_fdr(p_values)
    for r, q in zip(results, q_values):
        r["fdr"] = q
    results.sort(key=lambda r: (r["pvalue"], -r["hits"]))
    return results


def write_tsv(path: Path, results: list[dict]) -> None:
    with open(path, "w") as f:
        f.write(
            "pathway\tcollection\thits\tset_size_universe\tset_size_raw\t"
            "foreground_size\tuniverse_size\tpvalue\tfdr\thit_genes\n"
        )
        for r in results:
            f.write(
                f"{r['name']}\t{r['collection']}\t{r['hits']}\t"
                f"{r['set_size_universe']}\t{r['set_size_raw']}\t"
                f"{r['foreground_size']}\t{r['universe_size']}\t"
                f"{r['pvalue']:.6e}\t{r['fdr']:.6e}\t"
                f"{','.join(r['hit_genes'])}\n"
            )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("de_results", help="atman de_results.tsv")
    ap.add_argument("output_dir", help="output directory")
    ap.add_argument(
        "--genesets",
        default="/Users/kevinjoseph/Cursor/Rust_bioinfo_cli/genesets/data/genesets.json.gz",
    )
    ap.add_argument("--q-threshold", type=float, default=0.10)
    ap.add_argument("--top-n", type=int, default=100)
    ap.add_argument("--min-set-size", type=int, default=3)
    args = ap.parse_args()

    de_path = Path(args.de_results)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_sets = load_genesets(Path(args.genesets))
    by_collection: dict[str, list[dict]] = defaultdict(list)
    for gs in gene_sets:
        by_collection[gs["collection"]].append(gs)
    print(
        f"loaded {sum(len(v) for v in by_collection.values())} gene sets "
        f"across {len(by_collection)} collections: "
        f"{', '.join(f'{k}={len(v)}' for k, v in sorted(by_collection.items()))}",
        file=sys.stderr,
    )

    universe, rows = load_de(de_path)
    print(f"Olink universe: {len(universe)} unique gene symbols", file=sys.stderr)

    comparisons = sorted({r["comparison"] for r in rows})
    for cmp in comparisons:
        for direction in ("up", "down"):
            fg = top_hits(rows, cmp, direction, args.q_threshold, args.top_n)
            if not fg:
                continue
            for collection in sorted(by_collection):
                results = enrich_one(
                    fg, universe, by_collection[collection], args.min_set_size
                )
                out_path = out_dir / f"{cmp}_{direction}_{collection}_restricted.tsv"
                write_tsv(out_path, results)
                n_sig = sum(1 for r in results if r["fdr"] < 0.10)
                print(
                    f"  {cmp:10} {direction:4} {collection:8} n_fg={len(fg):3} "
                    f"n_tested_sets={len(results):5} n_fdr_lt_0.10={n_sig:4}",
                    file=sys.stderr,
                )

    print("done.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
