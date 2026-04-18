#!/usr/bin/env python3
"""Floor 1 — null calibration for the full module-DE pipeline.

For each of B sign-flipping permutations of the paired condition labels,
we re-learn modules on the permuted per-subject deltas and re-run
module-level paired-t + BH-FDR. Under the null (no systematic A−B effect),
a calibrated pipeline should reject ≤α·K modules on average at nominal
BH threshold α.

This tests the *full pipeline*: permutations that re-learn modules also
capture any "clustering-on-noise" inflation in the multiscale claim.

Inputs
------
    --measurements PATH    out/qc_measurements.tsv (or equivalent)
    --samples      PATH    out/samples.tsv
    --contrast     A-B     e.g. "PT2-PR2" or "late-early"
    --k-modules    INT     K passed to learn_modules (default 20)
    --n-perm       INT     number of sign-flip permutations (default 1000)
    --output       PATH    permutation_q_values.tsv (per perm × module)
    --summary      PATH    summary TSV with empirical rejection rate at α ∈ {0.01, 0.05, 0.10}

Usage
-----
    scripts/calibrate_module_de.py \
        --measurements out/qc_measurements.tsv \
        --samples out/samples.tsv \
        --contrast PT2-PR2 \
        --k-modules 20 \
        --n-perm 1000 \
        --output out/robustness/module_de_null_PT2-PR2.tsv \
        --summary out/robustness/module_de_null_PT2-PR2_summary.tsv
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import t as student_t


# ---- permutation + cluster core ------------------------------------------

def learn_modules_from_deltas(deltas: pd.DataFrame, k: int) -> pd.Series:
    """Per-gene module assignment via complete linkage on (1-r)/2 distance.

    `deltas` is long-format (subject_id, gene_symbol, delta). Returns a
    pandas Series indexed by gene_symbol with integer cluster labels.
    """
    wide = deltas.pivot_table(index="gene_symbol", columns="subject_id",
                              values="delta", aggfunc="mean").dropna()
    var = wide.var(axis=1, ddof=0)
    wide = wide[var > 0]
    if len(wide) < 2:
        return pd.Series([], dtype=int)
    X = wide.to_numpy(dtype=float)
    Xc = X - X.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(Xc, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    Xn = Xc / norms
    corr = np.clip(Xn @ Xn.T, -1.0, 1.0)
    dist = (1.0 - corr) / 2.0
    np.fill_diagonal(dist, 0.0)
    cond = squareform(dist, checks=False)
    Z = linkage(cond, method="complete")
    labels = fcluster(Z, t=min(k, len(wide)), criterion="maxclust")
    return pd.Series(labels, index=wide.index)


def paired_t_per_module(module_scores: pd.DataFrame) -> pd.DataFrame:
    """`module_scores` wide: index = subject_id, columns = (module, condition).

    Returns one-row-per-module DataFrame with mean_diff, t, p_value.
    """
    modules = sorted(set(col[0] for col in module_scores.columns))
    rows = []
    for mod in modules:
        a = module_scores.get((mod, "A"))
        b = module_scores.get((mod, "B"))
        if a is None or b is None:
            continue
        pair = pd.concat([a, b], axis=1).dropna()
        if len(pair) < 2:
            continue
        d = pair.iloc[:, 0] - pair.iloc[:, 1]
        n = len(d)
        mean_diff = d.mean()
        sd = d.std(ddof=1)
        if sd == 0 or not np.isfinite(sd):
            continue
        tstat = mean_diff / (sd / np.sqrt(n))
        p = 2.0 * (1.0 - student_t.cdf(abs(tstat), df=n - 1))
        rows.append({"module": mod, "n": n, "mean_diff": mean_diff, "t": tstat, "p_value": p})
    return pd.DataFrame(rows)


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR, returns q-values aligned with input order."""
    m = len(pvals)
    if m == 0:
        return pvals
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    q_sorted = sorted_p * m / np.arange(1, m + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    q = np.empty_like(pvals)
    q[order] = q_sorted
    return q


# ---- one permutation -----------------------------------------------------

def one_permutation(ab_long: pd.DataFrame, k: int, rng: np.random.Generator) -> pd.DataFrame:
    """Given long-format (subject, gene, condition_A, condition_B) values,
    sign-flip each subject's (A,B) pair with prob 0.5, re-learn modules,
    run module paired-t, return DataFrame of per-module q-values.
    """
    # Flip 50% of subjects
    subjects = ab_long["subject_id"].unique()
    flip_mask = rng.random(len(subjects)) < 0.5
    flip_map = dict(zip(subjects, flip_mask))
    a = ab_long["A"].copy()
    b = ab_long["B"].copy()
    flipped = ab_long["subject_id"].map(flip_map).to_numpy()
    ab_long = ab_long.copy()
    ab_long["A_p"] = np.where(flipped, b, a)
    ab_long["B_p"] = np.where(flipped, a, b)
    ab_long["delta_p"] = ab_long["A_p"] - ab_long["B_p"]
    # learn modules on permuted deltas
    labels = learn_modules_from_deltas(
        ab_long[["subject_id","gene_symbol"]].assign(delta=ab_long["delta_p"]), k)
    if len(labels) < 2:
        return pd.DataFrame()
    # build (subject, module, condition) scores under permuted labels.
    # score = mean of member gene abundances per subject per condition.
    lab_df = labels.rename("module").reset_index()  # gene_symbol, module
    joined = ab_long.merge(lab_df, on="gene_symbol")
    # Per (subject, module) mean of A_p and B_p across member genes
    agg = joined.groupby(["subject_id","module"])[["A_p","B_p"]].mean().reset_index()
    # wide scores: index=subject_id, cols=(module, 'A'/'B')
    wide = agg.pivot(index="subject_id", columns="module")
    # restructure columns to (module, 'A'/'B')
    ab_cols = {}
    for col in wide.columns:
        kind = "A" if col[0] == "A_p" else ("B" if col[0] == "B_p" else None)
        if kind is None:
            continue
        ab_cols[(col[1], kind)] = wide[col]
    scores = pd.DataFrame(ab_cols)
    scores.columns = pd.MultiIndex.from_tuples(scores.columns)
    res = paired_t_per_module(scores)
    if res.empty:
        return res
    res["bh_q"] = bh_fdr(res["p_value"].to_numpy())
    return res


# ---- main ----------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--measurements", required=True, type=Path)
    ap.add_argument("--samples",      required=True, type=Path)
    ap.add_argument("--contrast",     required=True, help="A-B, e.g. 'PT2-PR2'")
    ap.add_argument("--k-modules",    type=int, default=20)
    ap.add_argument("--n-perm",       type=int, default=1000)
    ap.add_argument("--seed",         type=int, default=42)
    ap.add_argument("--output",       required=True, type=Path)
    ap.add_argument("--summary",      required=True, type=Path)
    args = ap.parse_args()
    cond_a, cond_b = args.contrast.split("-")

    m = pd.read_csv(args.measurements, sep="\t")
    s = pd.read_csv(args.samples, sep="\t")
    m = m[m["dropped_by_qc"] == 0]
    m = m.merge(s[["sample_id","subject_id","condition"]], on="sample_id")
    # Wide per (subject, gene) × condition
    wide = m[m.condition.isin([cond_a, cond_b])] \
        .pivot_table(index=["subject_id","gene_symbol"], columns="condition",
                     values="abundance", aggfunc="mean").reset_index()
    wide = wide.dropna(subset=[cond_a, cond_b])
    ab_long = wide.rename(columns={cond_a: "A", cond_b: "B"})
    print(f"[calibrate] {len(ab_long):,} subject×gene pairs, "
          f"{ab_long.subject_id.nunique()} subjects, "
          f"{ab_long.gene_symbol.nunique()} genes", file=sys.stderr)

    rng = np.random.default_rng(args.seed)
    rows = []
    t0 = time.time()
    for b in range(args.n_perm):
        res = one_permutation(ab_long, args.k_modules, rng)
        if res.empty: continue
        res["perm"] = b
        rows.append(res[["perm","module","n","mean_diff","p_value","bh_q"]])
        if (b + 1) % 50 == 0 or b == args.n_perm - 1:
            elapsed = time.time() - t0
            print(f"[calibrate] perm {b+1}/{args.n_perm}  elapsed {elapsed:.1f}s  "
                  f"({elapsed/(b+1)*1000:.0f} ms/perm)", file=sys.stderr)

    all_rows = pd.concat(rows, ignore_index=True)
    all_rows.to_csv(args.output, sep="\t", index=False)

    # Empirical rejection rate at nominal α ∈ {0.01, 0.05, 0.10}
    per_perm = all_rows.groupby("perm").agg(
        k=("module", "count"),
        rej_0_01=("bh_q", lambda x: (x < 0.01).sum()),
        rej_0_05=("bh_q", lambda x: (x < 0.05).sum()),
        rej_0_10=("bh_q", lambda x: (x < 0.10).sum()),
    ).reset_index()
    summary = pd.DataFrame({
        "alpha": [0.01, 0.05, 0.10],
        "mean_rejections":     [per_perm.rej_0_01.mean(), per_perm.rej_0_05.mean(), per_perm.rej_0_10.mean()],
        "fraction_any":        [(per_perm.rej_0_01 > 0).mean(),
                                (per_perm.rej_0_05 > 0).mean(),
                                (per_perm.rej_0_10 > 0).mean()],
        "expected_if_calibrated": [0.01 * per_perm.k.mean(),
                                   0.05 * per_perm.k.mean(),
                                   0.10 * per_perm.k.mean()],
    })
    summary.to_csv(args.summary, sep="\t", index=False)
    print("\n==== Null calibration summary ====", file=sys.stderr)
    print(summary.to_string(index=False), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
