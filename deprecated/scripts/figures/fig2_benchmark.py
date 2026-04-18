#!/usr/bin/env python3
"""Fig 2 — benchmark against OlinkAnalyze and limma on the Dube contrasts.

Three standalone panels (one PDF+PNG per panel):
  fig2a_benchmark_yield       per-tool DE hit counts at q<0.05
  fig2b_benchmark_containment  Atman contains ≥94% of comparator hits at q<0.10
  fig2c_benchmark_concordance within-Atman paired-t vs moderated effect sizes
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _style import (
    apply_style, save, set_title,
    PRIMARY, ACCENT, MUTED, HIGHLIGHT,
    TOOL_COLORS, TOOL_LABELS, CONTRAST_COLORS,
    PANEL_SIZE, PANEL_SQUARE, PANEL_WIDE, REPO,
)

apply_style()
CONTRASTS = ["PT1-PR1", "PT2-PR2", "PT2-PT1", "PR2-PR1"]


def fig2a():
    t1 = pd.read_csv(REPO / "benchmarks/out/table1_counts.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=PANEL_SIZE)

    x = np.arange(len(CONTRASTS))
    w = 0.27
    for i, tool in enumerate(("atman", "OlinkAnalyze_ttest", "limma")):
        vals = [int(t1[(t1.tool == tool) & (t1.contrast == c)].n_q_lt_005.iloc[0])
                for c in CONTRASTS]
        bars = ax.bar(x + (i - 1) * w, vals, w,
                      color=TOOL_COLORS[tool],
                      edgecolor="white", linewidth=0.6,
                      label=TOOL_LABELS[tool])
        for b, v in zip(bars, vals):
            if v > 0:
                ax.annotate(str(v),
                            xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                            xytext=(0, 3), textcoords="offset points",
                            ha="center", fontsize=8, color="#374151")

    ax.set_xticks(x)
    ax.set_xticklabels(CONTRASTS)
    ax.set_ylabel("DE hits at q < 0.05")
    ax.set_ylim(0, max(ax.get_ylim()[1] * 1.12, 1))
    ax.legend(loc="upper right")
    set_title(ax, "Per-tool DE yield across four Dube contrasts")

    save(fig, "fig2a_benchmark_yield")
    plt.close(fig)


def fig2b():
    ov = pd.read_csv(REPO / "benchmarks/out/table1_overlap.tsv", sep="\t")
    ov = ov[(ov.tool_a == "atman") &
            (ov.contrast.isin(["PT1-PR1", "PT2-PR2"]))].copy()
    rows = [(r.contrast, r.tool_b, int(r.n_b), int(r.n_shared))
            for _, r in ov.iterrows()]

    fig, ax = plt.subplots(figsize=PANEL_SIZE)
    y = np.arange(len(rows))
    comparator = [r[2] for r in rows]
    shared     = [r[3] for r in rows]

    ax.barh(y, comparator, color=MUTED, alpha=0.35,
            edgecolor="none", height=0.65)
    ax.barh(y, shared, color=PRIMARY, edgecolor="white",
            linewidth=0.6, height=0.65)

    for i, (_, _, ntot, nsh) in enumerate(rows):
        pct = 100 * nsh / ntot if ntot else 0
        ax.annotate(f"  {nsh}/{ntot}  ({pct:.0f}%)",
                    xy=(ntot, i), xytext=(4, 0),
                    textcoords="offset points", va="center",
                    fontsize=9, color="#374151")

    labels = [f"{TOOL_LABELS[t]}\n{c}" for (c, t, _, _) in rows]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("comparator hits at q < 0.10")
    ax.set_xlim(0, max(comparator) * 1.70)

    # Legend above the plot to avoid collision with bar-end annotations.
    p_comparator = plt.Rectangle((0, 0), 1, 1, facecolor=MUTED, alpha=0.35)
    p_shared     = plt.Rectangle((0, 0), 1, 1, facecolor=PRIMARY)
    ax.legend([p_comparator, p_shared],
              ["comparator total", "shared with Atman"],
              loc="upper center", bbox_to_anchor=(0.5, -0.18),
              ncol=2, frameon=False)

    set_title(ax, "Atman contains ≥94% of comparator hits (q < 0.10)")
    save(fig, "fig2b_benchmark_containment")
    plt.close(fig)


def fig2c():
    p = pd.read_csv(REPO / "out/de_results.tsv", sep="\t")
    m = pd.read_csv(REPO / "out_mod/de_results.tsv", sep="\t")
    p = p.rename(columns={"mean_diff": "mean_diff_p"})
    m = m.rename(columns={"mean_diff": "mean_diff_m"})
    merged = p.merge(
        m[["panel", "assay_id", "comparison", "mean_diff_m"]],
        on=["panel", "assay_id", "comparison"], how="inner"
    ).dropna(subset=["mean_diff_p", "mean_diff_m"])

    fig, ax = plt.subplots(figsize=PANEL_SQUARE)
    lim = max(merged.mean_diff_p.abs().max(),
              merged.mean_diff_m.abs().max()) * 1.05

    # Identity line, drawn first so points sit on top.
    ax.plot([-lim, lim], [-lim, lim], ls="--", lw=0.7,
            color="#4B5563", zorder=0)

    for c in CONTRASTS:
        sub = merged[merged.comparison == c]
        ax.scatter(sub.mean_diff_p, sub.mean_diff_m,
                   s=6, alpha=0.45,
                   color=CONTRAST_COLORS[c],
                   edgecolor="none",
                   label=c)

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(r"paired-$t$ effect  ($\bar d$)")
    ax.set_ylabel(r"moderated effect  ($\bar d$)")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="upper left", markerscale=2, labelspacing=0.3)
    set_title(ax, r"Within-Atman mode concordance  ($\rho = 1.00$)")

    save(fig, "fig2c_benchmark_concordance")
    plt.close(fig)


def main():
    fig2a()
    fig2b()
    fig2c()


if __name__ == "__main__":
    main()
