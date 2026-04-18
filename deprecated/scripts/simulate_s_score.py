#!/usr/bin/env python3
"""Floor 2 — simulation study for the stability-weighted S-score.

Generates synthetic paired proteomics matrices with known ground truth
and benchmarks three top-K ranking rules (by q, by |d|, by S) under:

    n ∈ {8, 12, 16, 20, 24}                  subject count
    frac_TP ∈ {0.05, 0.15, 0.30}             fraction of true positives
    effect ∈ {strong, mixed, weak}           |d| distribution for TPs
    var_hetero ∈ {homo, hetero}              subject-level variance heterogeneity
    reps = 50 per cell

For each rep we run:
    1. paired-t per protein → q
    2. leave-one-subject-out reruns → SSR per protein
    3. S = |d| · SSR · (1 − q)
    4. for K = true-TP count, score top-K precision and recall under each rule

Output: one TSV row per (cell, rep, rule).
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student_t


def paired_t_vec(D: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised paired-t across rows of D (P × n of per-subject deltas).

    Returns (t_stat, p_value) as length-P vectors.
    """
    n = D.shape[1]
    mean = D.mean(axis=1)
    # ddof=1 sample sd
    sd = D.std(axis=1, ddof=1)
    sd = np.where(sd == 0, np.nan, sd)
    t = mean / (sd / np.sqrt(n))
    p = 2.0 * (1.0 - student_t.cdf(np.abs(t), df=n - 1))
    return t, p


def bh_fdr(p: np.ndarray) -> np.ndarray:
    valid = ~np.isnan(p)
    m = valid.sum()
    q = np.full_like(p, np.nan)
    if m == 0:
        return q
    idx = np.where(valid)[0]
    order = np.argsort(p[idx])
    sorted_p = p[idx][order]
    ranks = np.arange(1, m + 1)
    q_sorted = sorted_p * m / ranks
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    q_out = np.empty(m)
    q_out[order] = q_sorted
    q[idx] = q_out
    return q


def simulate_one(n: int, P: int, frac_tp: float, effect: str, hetero: bool,
                 fragile_frac: float,
                 rng: np.random.Generator) -> pd.DataFrame:
    """One synthetic dataset → rep-level metrics for q / |d| / S rankings.

    `fragile_frac` controls what fraction of TPs are outlier-driven
    (single-subject spike). Fragile TPs have |d| that looks large but SSR
    drops substantially — the regime where S should beat |d|.
    """
    n_tp = int(round(frac_tp * P))
    is_tp = np.zeros(P, dtype=bool)
    tp_idx = rng.choice(P, size=n_tp, replace=False)
    is_tp[tp_idx] = True
    # Split TPs: robust (whole-cohort signal) vs fragile (single-subject outlier)
    n_fragile = int(round(fragile_frac * n_tp))
    fragile_idx = tp_idx[:n_fragile] if n_fragile > 0 else np.array([], dtype=int)
    robust_idx  = tp_idx[n_fragile:]

    # True robust effects
    if effect == "strong":
        d_true_robust = rng.normal(1.5, 0.3, size=len(robust_idx)) * rng.choice([-1, 1], size=len(robust_idx))
    elif effect == "mixed":
        d_true_robust = rng.normal(0.5, 0.5, size=len(robust_idx)) * rng.choice([-1, 1], size=len(robust_idx))
    elif effect == "weak":
        d_true_robust = rng.normal(0.3, 0.1, size=len(robust_idx)) * rng.choice([-1, 1], size=len(robust_idx))
    else:
        raise ValueError(effect)

    d_true = np.zeros(P)
    d_true[robust_idx] = d_true_robust
    # Fragile TPs: no systematic signal, but inject a large outlier into one random subject
    # (the outlier-to-be-injected is added to the noise matrix later)

    # Subject-level variance
    if hetero:
        sigma_s = rng.uniform(0.4, 1.2, size=n)
    else:
        sigma_s = np.full(n, 0.7)

    # Generate per-subject paired deltas: D[i, s] = d_true[i] + noise
    noise = rng.normal(0.0, 1.0, size=(P, n)) * sigma_s[None, :]
    D = d_true[:, None] + noise

    # Inject outlier spikes for fragile TPs (1–2 subjects carry the signal)
    for i in fragile_idx:
        k = rng.integers(0, n)                 # chosen subject
        sign = rng.choice([-1, 1])
        spike = sign * rng.uniform(4.0, 7.0)   # large enough to drive |d_bar| up
        D[i, k] += spike

    # Baseline paired-t
    t, p = paired_t_vec(D)
    q = bh_fdr(p)
    d_bar = D.mean(axis=1)

    # LOO SSR (vectorised across leave-out subjects)
    n_loo = n
    sign_match = np.zeros(P, dtype=int)
    n_valid = np.zeros(P, dtype=int)
    base_sign = np.sign(d_bar)
    for k in range(n):
        D_loo = np.delete(D, k, axis=1)
        d_loo = D_loo.mean(axis=1)
        loo_sign = np.sign(d_loo)
        valid = (base_sign != 0) & (loo_sign != 0)
        sign_match += ((base_sign == loo_sign) & valid).astype(int)
        n_valid += valid.astype(int)
    ssr = np.where(n_valid > 0, sign_match / np.maximum(n_valid, 1), 0.0)

    # S-score
    q_for_s = np.nan_to_num(q, nan=1.0)
    S = np.abs(d_bar) * ssr * (1 - q_for_s)

    # Rank by each criterion and compute top-K precision/recall at K=n_tp
    K = max(n_tp, 1)
    # Ranks: smallest-q highest rank, largest-|d| highest rank, largest-S highest rank
    def topk(scores_desc):
        idx = np.argsort(-scores_desc)[:K]
        tp_hits = is_tp[idx].sum()
        return tp_hits / K, tp_hits / max(n_tp, 1)  # precision, recall

    prec_q, rec_q = topk(-np.nan_to_num(q, nan=2.0))  # smaller q = higher rank → use -q
    prec_d, rec_d = topk(np.abs(d_bar))
    prec_s, rec_s = topk(S)

    rows = []
    for rule, prec, rec in [("q", prec_q, rec_q), ("abs_d", prec_d, rec_d), ("S", prec_s, rec_s)]:
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
        rows.append({"rule": rule, "precision_at_K": prec, "recall_at_K": rec, "f1_at_K": f1})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--P", type=int, default=500, help="number of proteins")
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    ns        = [8, 12, 16, 20, 24]
    fracs     = [0.05, 0.15, 0.30]
    effects   = ["strong", "mixed", "weak"]
    heteros   = [False, True]
    fragiles  = [0.0, 0.25, 0.5]  # fraction of TPs that are outlier-driven

    rows = []
    t0 = time.time()
    total_cells = len(ns) * len(fracs) * len(effects) * len(heteros) * len(fragiles)
    cell_i = 0
    for n in ns:
        for f in fracs:
            for e in effects:
                for h in heteros:
                    for fr in fragiles:
                        cell_i += 1
                        cell_df = []
                        for r in range(args.reps):
                            res = simulate_one(n, args.P, f, e, h, fr, rng)
                            res["rep"] = r
                            cell_df.append(res)
                        cdf = pd.concat(cell_df, ignore_index=True)
                        cdf["n"] = n; cdf["frac_tp"] = f; cdf["effect"] = e
                        cdf["hetero"] = h; cdf["fragile"] = fr
                        rows.append(cdf)
                        elapsed = time.time() - t0
                        if cell_i % 30 == 0 or cell_i == total_cells:
                            print(f"[sim] cell {cell_i}/{total_cells}  elapsed {elapsed:.1f}s", file=sys.stderr)
    df = pd.concat(rows, ignore_index=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {args.output}  ({len(df)} rows across {total_cells} cells × {args.reps} reps)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
