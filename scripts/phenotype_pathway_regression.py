#!/usr/bin/env python3
"""Pathway-level OLS regression of physiological heat-response endpoints on
pathway scores aggregated from per-subject delta NPX.

Per-protein regression (scripts/phenotype_regression.py) is severely
underpowered at n=10 — to pass BH-FDR<0.05 across 2,940 protein tests,
each protein needs R²>0.85, which only happens by overfitting on outliers.
Pathway aggregation reduces dimensionality from 2,940 to ~50-700 per
collection, giving meaningful statistical resolution.

For each pathway P in {hallmark, kegg, go_bp}:
    P_score(subject) = mean of dNPX_<contrast>(gene, subject)
                       over all genes in P that are in the Olink universe
For each (endpoint, pathway):
    OLS:  endpoint ~ α + β · P_score
    n = 10 subjects, df = 8.
BH-FDR within (endpoint × collection).

Reads the same dNPX files and Physiological_data.xlsx as phenotype_regression.py.
Restricts pathways to the Olink universe (gene_symbol set from de_results.tsv).

Usage:
    scripts/phenotype_pathway_regression.py <de_results.tsv> <output_dir>
"""

import csv
import gzip
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import t as student_t

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = REPO_ROOT / "example_data" / "dube_heat_2023"
DELTA_DIR = DATA_ROOT / "delta_npx"
GENESETS_PATH = Path(
    "/Users/kevinjoseph/Cursor/Rust_bioinfo_cli/genesets/data/genesets.json.gz"
)

# Reuse loaders from the per-protein script.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phenotype_regression import (
    PARTICIPANTS,
    PHENO_VARS,
    bh_fdr,
    build_endpoints,
    load_dnpx,
    load_phenotypes,
    ols_slope,
)

ENDPOINT_TO_DNPX = {"d1_": "hd1", "d7_": "hd7", "acc_": "accpre"}
MIN_PATHWAY_SIZE = 5  # in Olink universe


def load_genesets() -> List[dict]:
    with gzip.open(GENESETS_PATH, "rt") as f:
        return json.load(f)["sets"]


def load_olink_universe(de_results_path: Path) -> set:
    universe = set()
    with open(de_results_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r.get("gene_symbol", "")
            if g:
                universe.add(g)
    return universe


def pathway_scores(
    pathway_genes: set,
    dnpx: Dict[Tuple[str, str, str], Dict[str, float]],
) -> Dict[str, float]:
    """Per-subject pathway score = mean of dNPX over genes in pathway that
    have a measurement for that subject. Returns NaN if fewer than 2 genes
    contribute for a subject."""
    # Build gene -> {subject: value} from dnpx.
    gene_vals: Dict[str, Dict[str, float]] = defaultdict(dict)
    for (_oid, gene, _panel), per_subject in dnpx.items():
        if gene in pathway_genes:
            for s, v in per_subject.items():
                gene_vals[gene][s] = v
    scores: Dict[str, float] = {}
    for s in PARTICIPANTS:
        contributions = [
            v.get(s) for v in gene_vals.values() if s in v and math.isfinite(v[s])
        ]
        if len(contributions) >= 2:
            scores[s] = float(np.mean(contributions))
        else:
            scores[s] = math.nan
    return scores


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: phenotype_pathway_regression.py <de_results.tsv> <output_dir>", file=sys.stderr)
        return 2
    de_results_path = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    universe = load_olink_universe(de_results_path)
    print(f"Olink universe: {len(universe)} unique genes", file=sys.stderr)

    sets = load_genesets()
    by_collection: Dict[str, List[dict]] = defaultdict(list)
    for gs in sets:
        in_univ = set(gs["genes"]) & universe
        if len(in_univ) >= MIN_PATHWAY_SIZE:
            by_collection[gs["collection"]].append(
                {"name": gs["name"], "genes": in_univ, "size_raw": len(gs["genes"])}
            )
    print(
        "pathways with >= "
        f"{MIN_PATHWAY_SIZE} Olink genes: "
        + ", ".join(f"{k}={len(v)}" for k, v in sorted(by_collection.items())),
        file=sys.stderr,
    )

    pheno = load_phenotypes()
    endpoints = build_endpoints(pheno)
    dnpx_files = {
        "hd1": load_dnpx(DELTA_DIR / "dNPX_hd1.csv"),
        "hd7": load_dnpx(DELTA_DIR / "dNPX_hd7.csv"),
        "accpre": load_dnpx(DELTA_DIR / "dNPX_accpre.csv"),
    }

    all_rows: List[dict] = []

    for ep_name, ep_vals in endpoints.items():
        dnpx_tag = next(
            (tag for prefix, tag in ENDPOINT_TO_DNPX.items() if ep_name.startswith(prefix)),
            None,
        )
        if dnpx_tag is None:
            continue
        dnpx = dnpx_files[dnpx_tag]
        y = np.array([ep_vals.get(p, math.nan) for p in PARTICIPANTS], dtype=float)
        if np.sum(np.isfinite(y)) < 5:
            continue

        for collection, pathways in sorted(by_collection.items()):
            family_rows: List[dict] = []
            family_p: List[float] = []
            for pw in pathways:
                scores = pathway_scores(pw["genes"], dnpx)
                x = np.array([scores.get(p, math.nan) for p in PARTICIPANTS], dtype=float)
                n, beta, se, t_stat, p_val, r2 = ols_slope(x, y)
                row = {
                    "endpoint": ep_name,
                    "dnpx_source": dnpx_tag,
                    "collection": collection,
                    "pathway": pw["name"],
                    "set_size_raw": pw["size_raw"],
                    "set_size_universe": len(pw["genes"]),
                    "n": n,
                    "beta": beta,
                    "se": se,
                    "t": t_stat,
                    "p_value": p_val,
                    "r_squared": r2,
                }
                family_rows.append(row)
                family_p.append(p_val)
            qs = bh_fdr(family_p)
            for r, q in zip(family_rows, qs):
                r["bh_q"] = q
            all_rows.extend(family_rows)

        n_q05 = sum(
            1
            for r in all_rows
            if r["endpoint"] == ep_name
            and not math.isnan(r.get("bh_q", math.nan))
            and r["bh_q"] < 0.05
        )
        n_q10 = sum(
            1
            for r in all_rows
            if r["endpoint"] == ep_name
            and not math.isnan(r.get("bh_q", math.nan))
            and r["bh_q"] < 0.10
        )
        print(
            f"  {ep_name:12} ({dnpx_tag:6}) q<0.05={n_q05:3} q<0.10={n_q10:3}",
            file=sys.stderr,
        )

    # Sort: endpoint, collection, q ascending.
    def sort_key(r: dict):
        q = r.get("bh_q", math.nan)
        q_sort = q if not math.isnan(q) else math.inf
        return (r["endpoint"], r["collection"], q_sort)

    all_rows.sort(key=sort_key)

    out_path = out_dir / "phenotype_pathway_regression.tsv"
    with open(out_path, "w") as f:
        f.write(
            "endpoint\tcollection\tpathway\tset_size_universe\tset_size_raw\t"
            "dnpx_source\tn\tbeta\tse\tt\tp_value\tr_squared\tbh_q\n"
        )
        for r in all_rows:
            f.write(
                f"{r['endpoint']}\t{r['collection']}\t{r['pathway']}\t"
                f"{r['set_size_universe']}\t{r['set_size_raw']}\t{r['dnpx_source']}\t"
                f"{r['n']}\t{_fmt(r['beta'])}\t{_fmt(r['se'])}\t{_fmt(r['t'])}\t"
                f"{_fmt(r['p_value'])}\t{_fmt(r['r_squared'])}\t{_fmt(r['bh_q'])}\n"
            )
    print(f"\nwrote {out_path} ({len(all_rows)} rows)", file=sys.stderr)
    return 0


def _fmt(v) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    return f"{v}"


if __name__ == "__main__":
    sys.exit(main())
