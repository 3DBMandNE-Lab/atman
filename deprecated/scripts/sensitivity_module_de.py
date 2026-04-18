#!/usr/bin/env python3
"""Floor 1b — sensitivity of module-DE to linkage and K.

Sweeps linkage ∈ {single, average, complete, ward} and K ∈ {5,10,15,20,25,30,40,50}
on a chosen dataset/contrast. For each (linkage, K):
    - compute observed-data hit count at q<0.05 and q<0.10
    - run B=100 sign-flip permutations and report empirical rejection rate

Plus a modularity/silhouette score vs K (at fixed linkage=complete) to
expose the elbow.

Output: single TSV with all cells.
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform, pdist
from scipy.stats import t as student_t
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)  # matmul BLAS noise


def learn(wide: pd.DataFrame, method: str, k: int) -> pd.Series:
    X = wide.to_numpy(dtype=float)
    Xc = X - X.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(Xc, axis=1, keepdims=True); norms[norms == 0] = 1.0
    Xn = Xc / norms
    corr = np.clip(Xn @ Xn.T, -1.0, 1.0)
    dist = (1.0 - corr) / 2.0
    np.fill_diagonal(dist, 0.0)
    cond = squareform(dist, checks=False)
    # Ward requires squared Euclidean; we approximate by passing the
    # precomputed distance via scipy's default linkage semantics.
    try:
        Z = linkage(cond, method=method)
    except Exception:
        return pd.Series([], dtype=int)
    labels = fcluster(Z, t=min(k, len(wide)), criterion="maxclust")
    return pd.Series(labels, index=wide.index)


def silhouette_on_dist(dist: np.ndarray, labels: np.ndarray) -> float:
    """Silhouette mean; dist is a square distance matrix; labels aligned."""
    n = dist.shape[0]
    if len(set(labels)) < 2 or len(set(labels)) >= n:
        return float("nan")
    sils = []
    for i in range(n):
        same = (labels == labels[i]) & (np.arange(n) != i)
        other = (labels != labels[i])
        if same.sum() == 0:
            continue
        a = dist[i, same].mean()
        # mean distance to nearest other cluster
        other_labels = set(labels[other])
        if not other_labels:
            continue
        b = min(dist[i, labels == lab].mean() for lab in other_labels)
        sils.append((b - a) / max(a, b) if max(a, b) > 0 else 0.0)
    return float(np.mean(sils)) if sils else float("nan")


def module_de(ab: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    """Given long df with (subject, gene, A, B) and gene→module labels,
    compute per-module paired-t + BH q."""
    lab_df = labels.rename("module").reset_index()
    joined = ab.merge(lab_df, on="gene_symbol")
    agg = joined.groupby(["subject_id","module"])[["A","B"]].mean().reset_index()
    wide = agg.pivot(index="subject_id", columns="module")
    rows = []
    for mod in labels.unique():
        try:
            a = wide[("A", mod)]
            b = wide[("B", mod)]
        except KeyError:
            continue
        pair = pd.concat([a, b], axis=1).dropna()
        if len(pair) < 2:
            continue
        d = pair.iloc[:, 0] - pair.iloc[:, 1]
        sd = d.std(ddof=1)
        if not np.isfinite(sd) or sd == 0: continue
        t = d.mean() / (sd / np.sqrt(len(d)))
        p = 2.0 * (1.0 - student_t.cdf(abs(t), df=len(d) - 1))
        rows.append({"module": mod, "mean_diff": d.mean(), "t": t, "p_value": p})
    res = pd.DataFrame(rows)
    if res.empty: return res
    # BH
    p = res.p_value.to_numpy()
    m = len(p)
    order = np.argsort(p)
    q = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty_like(p); out[order] = q
    res["bh_q"] = out
    return res


def one_perm(ab: pd.DataFrame, wide: pd.DataFrame, method: str, k: int,
             rng: np.random.Generator):
    subjects = ab["subject_id"].unique()
    flip = dict(zip(subjects, rng.random(len(subjects)) < 0.5))
    mask = ab["subject_id"].map(flip).to_numpy()
    ab_p = ab.copy()
    ab_p["A"] = np.where(mask, ab["B"], ab["A"])
    ab_p["B"] = np.where(mask, ab["A"], ab["B"])
    # rebuild wide deltas for module learning
    w = ab_p.assign(delta=ab_p["A"] - ab_p["B"])\
            .pivot_table(index="gene_symbol", columns="subject_id",
                         values="delta", aggfunc="mean").dropna()
    var = w.var(axis=1, ddof=0); w = w[var > 0]
    if len(w) < 2: return None
    labels = learn(w, method, k)
    if labels.empty: return None
    return module_de(ab_p, labels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--measurements", required=True, type=Path)
    ap.add_argument("--samples", required=True, type=Path)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--n-perm", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    cond_a, cond_b = args.contrast.split("-")

    m = pd.read_csv(args.measurements, sep="\t")
    s = pd.read_csv(args.samples, sep="\t")
    m = m[m.dropped_by_qc == 0].merge(s[["sample_id","subject_id","condition"]], on="sample_id")
    wide = m[m.condition.isin([cond_a, cond_b])]\
        .pivot_table(index=["subject_id","gene_symbol"], columns="condition",
                     values="abundance", aggfunc="mean").reset_index().dropna(subset=[cond_a, cond_b])
    ab = wide.rename(columns={cond_a: "A", cond_b: "B"})
    # precompute delta-wide for observed module learning
    obs_wide = ab.assign(delta=ab["A"] - ab["B"])\
        .pivot_table(index="gene_symbol", columns="subject_id",
                     values="delta", aggfunc="mean").dropna()
    obs_wide = obs_wide[obs_wide.var(axis=1, ddof=0) > 0]

    # Precomputed observed-data distance (for silhouette)
    X = obs_wide.to_numpy(dtype=float)
    Xc = X - X.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(Xc, axis=1, keepdims=True); norms[norms == 0] = 1.0
    corr = np.clip((Xc / norms) @ (Xc / norms).T, -1.0, 1.0)
    dist = (1.0 - corr) / 2.0
    np.fill_diagonal(dist, 0.0)

    rng = np.random.default_rng(args.seed)
    methods = ["single", "average", "complete", "ward"]
    Ks      = [5, 10, 15, 20, 25, 30, 40, 50]
    rows = []
    t0 = time.time()
    for method in methods:
        for k in Ks:
            labels = learn(obs_wide, method, k)
            if labels.empty: continue
            sil = silhouette_on_dist(dist, labels.to_numpy())
            obs = module_de(ab, labels)
            n05 = int((obs.bh_q < 0.05).sum()) if not obs.empty else 0
            n10 = int((obs.bh_q < 0.10).sum()) if not obs.empty else 0
            # permutation FDR
            null_05 = 0; null_10 = 0; null_tests = 0
            for _ in range(args.n_perm):
                p = one_perm(ab, obs_wide, method, k, rng)
                if p is None or p.empty: continue
                null_05 += int((p.bh_q < 0.05).sum())
                null_10 += int((p.bh_q < 0.10).sum())
                null_tests += len(p)
            mean_rej_05 = null_05 / args.n_perm
            mean_rej_10 = null_10 / args.n_perm
            rows.append({
                "linkage": method, "K": k, "K_realised": labels.nunique(),
                "silhouette": round(sil, 4) if np.isfinite(sil) else None,
                "obs_q005": n05, "obs_q010": n10,
                "null_mean_rej_q005": round(mean_rej_05, 4),
                "null_mean_rej_q010": round(mean_rej_10, 4),
                "expected_if_calibrated_q005": round(0.05 * labels.nunique(), 3),
                "expected_if_calibrated_q010": round(0.10 * labels.nunique(), 3),
            })
            elapsed = time.time() - t0
            print(f"  {method:9s} K={k:3d} obs[q05={n05:2d},q10={n10:2d}] null_mean[q05={mean_rej_05:.3f},q10={mean_rej_10:.3f}] sil={sil:.3f}  {elapsed:.0f}s",
                  file=sys.stderr)
    df = pd.DataFrame(rows)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"\nwrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
