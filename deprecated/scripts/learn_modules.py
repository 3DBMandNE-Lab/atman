#!/usr/bin/env python3
"""Data-driven protein-module learning.

Reads a per-subject-gene-delta table (columns `subject_id`, `gene_symbol`,
one delta column per contrast) and emits a `modules.tsv` in the format
`atman module-de` consumes (columns `module`, `gene_symbol`).

Method
------
For each protein, concatenate its per-subject delta values across all
contrast columns into a single feature vector. Compute pairwise distance
$1 - |\mathrm{Pearson}(\mathbf{v}_i, \mathbf{v}_j)|$ between proteins.
Absolute correlation is used so that proteins moving in anti-correlated
directions cluster together when that reflects a single up/down regulon.
Ward linkage is applied to the condensed distance vector and the tree is
cut at `--k-modules` (default 20). Modules are renamed to `module_01`,
`module_02`, ... in descending size order.

Missing values: proteins with any NA in their feature vector are dropped
with a warning; proteins with zero variance across subjects are also
dropped (distance to every other protein is undefined).

Usage
-----
    scripts/learn_modules.py <deltas.tsv> <modules.tsv> [--k-modules 20]
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deltas_tsv", type=Path)
    ap.add_argument("modules_tsv", type=Path)
    ap.add_argument("--k-modules", type=int, default=20)
    ap.add_argument("--contrast-cols", default=None,
                    help="Comma-separated delta column names; default = all numeric cols except subject_id/gene_symbol/panel")
    ap.add_argument("--summary", type=Path, default=None,
                    help="Optional path for module_summary.tsv with size + mean |within-module r|")
    args = ap.parse_args()

    d = pd.read_csv(args.deltas_tsv, sep="\t")
    meta_cols = {"subject_id", "gene_symbol", "panel", "assay_id", "uniprot"}
    if args.contrast_cols:
        ccols = [c.strip() for c in args.contrast_cols.split(",") if c.strip()]
    else:
        ccols = [c for c in d.columns if c not in meta_cols and pd.api.types.is_numeric_dtype(d[c])]
    print(f"deltas: {d.shape}  using {len(ccols)} contrast columns: {ccols}", file=sys.stderr)

    # Pivot to protein × (subject, contrast) feature matrix.
    # Unique protein key = gene_symbol (panel ignored — same gene on multiple panels
    # gets collapsed; this is the standard module-analysis convention).
    long = d.melt(id_vars=[c for c in ("subject_id","gene_symbol") if c in d.columns],
                  value_vars=ccols, var_name="contrast", value_name="delta")
    long["feature"] = long["subject_id"].astype(str) + "__" + long["contrast"]
    wide = long.pivot_table(index="gene_symbol", columns="feature", values="delta", aggfunc="mean")
    n_before = len(wide)

    # Drop proteins with any NA or zero variance.
    wide = wide.dropna()
    var = wide.var(axis=1, ddof=0)
    wide = wide[var > 0]
    n_after = len(wide)
    print(f"proteins: {n_before} → {n_after} after dropping NA / zero-var", file=sys.stderr)

    # Pairwise signed-correlation distance $d_{ij} = (1 - r_{ij})/2 \in [0,1]$.
    # Signed (not absolute) correlation keeps anti-correlated proteins far apart,
    # which produces compact up-regulon and down-regulon modules rather than
    # lumping them together. Complete linkage is used in place of average
    # linkage because average linkage exhibits severe chaining on high-P
    # proteomics correlation matrices (two giant modules absorbing most genes).
    X = wide.to_numpy(dtype=float)
    Xc = X - X.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(Xc, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    Xn = Xc / norms
    corr = np.clip(Xn @ Xn.T, -1.0, 1.0)
    dist_sq = (1.0 - corr) / 2.0
    np.fill_diagonal(dist_sq, 0.0)
    cond = squareform(dist_sq, checks=False)
    Z = linkage(cond, method="complete")
    labels = fcluster(Z, t=args.k_modules, criterion="maxclust")
    wide["_cluster"] = labels

    # Rename clusters in descending-size order: module_01 is the largest.
    size = wide.groupby("_cluster").size().sort_values(ascending=False)
    rename = {old: f"module_{i+1:02d}" for i, old in enumerate(size.index)}
    wide["module"] = wide["_cluster"].map(rename)

    # Write modules.tsv (atman module-de format).
    out = wide[["module"]].reset_index()[["module", "gene_symbol"]].sort_values(["module", "gene_symbol"])
    out.to_csv(args.modules_tsv, sep="\t", index=False)
    print(f"wrote {args.modules_tsv} — {out.module.nunique()} modules, {len(out)} gene assignments", file=sys.stderr)

    # Optional summary with within-module mean |correlation|.
    if args.summary is not None:
        rows = []
        for m, sub in wide.groupby("module"):
            idx = [wide.index.get_loc(g) for g in sub.index]
            if len(idx) < 2:
                mean_r = float("nan")
            else:
                block = corr[np.ix_(idx, idx)]
                iu = np.triu_indices(len(idx), k=1)
                mean_r = float(np.abs(block[iu]).mean())
            rows.append({"module": m, "n_genes": len(idx), "mean_abs_within_r": mean_r})
        summ = pd.DataFrame(rows).sort_values("n_genes", ascending=False)
        summ.to_csv(args.summary, sep="\t", index=False)
        print(f"wrote {args.summary}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
