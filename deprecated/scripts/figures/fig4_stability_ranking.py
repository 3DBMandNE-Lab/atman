#!/usr/bin/env python3
"""Fig 5 — stability-aware ranking (legacy S-score characterisation).

Two standalone panels:
  fig5a_ranking_loo_jaccard  top-K LOO Jaccard under three rankings (q, |d|, S)
  fig5b_q_vs_s_rank          rank-by-q vs rank-by-S on Dube PT1-PR1
"""
from __future__ import annotations
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _style import (
    apply_style, save, set_title,
    PRIMARY, ACCENT, MUTED, HIGHLIGHT, SEQUENTIAL,
    PANEL_WIDE, PANEL_SQUARE, REPO,
)

apply_style()


def loo_jaccards(base_path, loo_glob, ssr_path, contrast, K=None):
    base = pd.read_csv(base_path, sep="\t")
    base["key"] = base["panel"] + "::" + base["assay_id"]
    ssr = pd.read_csv(ssr_path, sep="\t")
    ssr["key"] = ssr["panel"] + "::" + ssr["assay_id"]
    ssr_map = {(r.comparison, r.key): r.sign_match_rate
               for _, r in ssr.iterrows()
               if pd.notna(r.sign_match_rate)}
    b = base[base.comparison == contrast].copy()
    b["ssr"] = [ssr_map.get((contrast, k), 0.0) for k in b.key]
    b["abs_d"] = b.mean_diff.abs()
    b["S"] = b["abs_d"] * b["ssr"] * (1 - b.bh_q)
    if K is None:
        K = int((b.bh_q < 0.05).sum())
    full = {
        "q":   set(b.nsmallest(K, "bh_q").key),
        "|d|": set(b.nlargest(K, "abs_d").key),
        "S":   set(b.nlargest(K, "S").key),
    }
    over = {c: [] for c in full}
    for lp in sorted(glob.glob(loo_glob)):
        l = pd.read_csv(lp, sep="\t")
        l["key"] = l["panel"] + "::" + l["assay_id"]
        l = l[l.comparison == contrast].copy()
        l["ssr"] = [ssr_map.get((contrast, k), 0.0) for k in l.key]
        l["abs_d"] = l.mean_diff.abs()
        l["S"] = l["abs_d"] * l["ssr"] * (1 - l.bh_q)
        tops = {
            "q":   set(l.nsmallest(K, "bh_q").key),
            "|d|": set(l.nlargest(K, "abs_d").key),
            "S":   set(l.nlargest(K, "S").key),
        }
        for c in full:
            inter = len(full[c] & tops[c])
            union = len(full[c] | tops[c])
            over[c].append(inter / union if union else float("nan"))
    return K, {c: (float(np.mean(v)), float(np.std(v))) for c, v in over.items()}


def fig5a():
    configs = [
        ("Dube\nPT1-PR1",     REPO/"out/de_results.tsv",
         REPO/"out_dube_loo/loo_*/de_results.tsv",
         REPO/"out/robustness/loo_sign_stability.tsv", "PT1-PR1"),
        ("Dube\nPT2-PR2",     REPO/"out/de_results.tsv",
         REPO/"out_dube_loo/loo_*/de_results.tsv",
         REPO/"out/robustness/loo_sign_stability.tsv", "PT2-PR2"),
        ("Dube\nPT2-PT1",     REPO/"out/de_results.tsv",
         REPO/"out_dube_loo/loo_*/de_results.tsv",
         REPO/"out/robustness/loo_sign_stability.tsv", "PT2-PT1"),
        ("Gisby\nlate-early", REPO/"out_gisby/de_results.tsv",
         REPO/"out_gisby_loo/loo_*/de_results.tsv",
         REPO/"out_gisby/robustness/loo_sign_stability.tsv", "late-early"),
    ]
    labels, data = [], []
    for label, base, loo, ssr, contrast in configs:
        K, r = loo_jaccards(str(base), str(loo), str(ssr), contrast)
        labels.append(f"{label}\n(K={K})")
        data.append(r)

    fig, ax = plt.subplots(figsize=PANEL_WIDE)
    x = np.arange(len(labels))
    w = 0.26
    crit_colors = {"q": MUTED, "|d|": ACCENT, "S": PRIMARY}
    for i, crit in enumerate(("q", "|d|", "S")):
        means = [data[j][crit][0] for j in range(len(labels))]
        sds   = [data[j][crit][1] for j in range(len(labels))]
        ax.bar(x + (i - 1) * w, means, w, yerr=sds, capsize=2.5,
               label=f"by {crit}", color=crit_colors[crit],
               edgecolor="white", linewidth=0.6,
               error_kw=dict(elinewidth=0.7, ecolor="#374151"))
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("top-K Jaccard under LOO\n(mean ± sd)")
    ax.set_ylim(0, 1.18)  # headroom for legend above bars
    ax.axhline(1.0, ls=":", lw=0.5, color="#9CA3AF", zorder=0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.0),
              ncol=3, columnspacing=1.4)
    set_title(ax, "Ranking reproducibility under LOO resampling")
    save(fig, "fig4a_ranking_loo_jaccard")
    plt.close(fig)


def fig5b():
    st  = pd.read_csv(REPO / "out/robustness/stability_ranked.tsv", sep="\t")
    sub = st[st.comparison == "PT1-PR1"].dropna(
        subset=["baseline_bh_q", "sign_match_rate", "stability_score"]).copy()
    sub["abs_d"] = sub["baseline_mean_diff"].abs()
    keep = sub[(sub.rank_by_q <= 100) | (sub.rank_by_stability <= 100)].copy()
    lim = 200
    keep["rq"] = keep.rank_by_q.clip(upper=lim)
    keep["rs"] = keep.rank_by_stability.clip(upper=lim)

    fig, ax = plt.subplots(figsize=PANEL_SQUARE)
    ax.plot([0, lim], [0, lim], ls="--", lw=0.7, color="#4B5563", zorder=0)
    sc = ax.scatter(keep.rq, keep.rs, c=keep.abs_d, cmap=SEQUENTIAL,
                    s=22, edgecolor="white", linewidth=0.4,
                    vmin=0, vmax=keep.abs_d.quantile(0.95))
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_aspect("equal", adjustable="box")
    # Do not invert: rank 1 = top is at the origin (bottom-left); the identity
    # diagonal runs from origin to (lim, lim). Points below the diagonal are
    # q-rank > S-rank → promoted by S; points above → demoted by S.
    ax.set_xlabel(r"rank by $q$   (1 at origin)")
    ax.set_ylabel(r"rank by $S$   (1 at origin)")
    cbar = plt.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label(r"$|\bar d|$   ($\log_2$ NPX)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_linewidth(0)
    set_title(ax, r"Dube PT1-PR1: $q$-rank vs $S$-rank")
    save(fig, "fig4b_q_vs_s_rank")
    plt.close(fig)


def main():
    fig5a()
    fig5b()


if __name__ == "__main__":
    main()
