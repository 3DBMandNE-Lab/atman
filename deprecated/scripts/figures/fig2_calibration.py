#!/usr/bin/env python3
"""Fig 2 — calibration and null behaviour.

Four panels, 2x2 at double-column width:
  a. QQ plot of module-DE p-values under sign-flip null (three datasets).
  b. Empirical module-DE FDR vs nominal BH-FDR cutoff α, same datasets.
  c. Discovery power (recall at K = true-TP count) vs n, three rankings,
     from the 270-cell synthetic simulation (baseline regime:
     frac_tp = 0.15, mixed effect, homoscedastic, no fragile TPs).
  d. Empirical protein-DE FPR at nominal q<0.05 by method, from the
     Floor 3 permutation benchmark (200 sign-flip perms x 4 Dube
     contrasts). Every method sits 1-4 orders of magnitude below nominal;
     the paper's "paired-t not inflated" claim made visual.

All numbers come from existing TSVs; no re-computation.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _style import (
    apply_style, save, panel_letter, clean_axes,
    PRIMARY, ACCENT, SECOND, LEGACY, MUTED, STRUCTURE, REPO, GRID_2x2,
)

apply_style()


# ------------------------------------------------------------------
# panel a — QQ of null p-values
# ------------------------------------------------------------------
def panel_a(ax):
    files = [
        ("Dube PT1-PR1",  REPO / "out/robustness/module_de_null_PT1-PR1.tsv",      PRIMARY),
        ("Dube PT2-PR2",  REPO / "out/robustness/module_de_null_PT2-PR2.tsv",      ACCENT),
        ("Gisby late-early", REPO / "out_gisby/robustness/module_de_null_late-early.tsv", SECOND),
    ]
    ax.plot([0, 1], [0, 1], color=MUTED, lw=0.5, ls="--", zorder=0)
    for label, path, col in files:
        df = pd.read_csv(path, sep="\t")
        p = np.sort(df["p_value"].dropna().to_numpy())
        n = len(p)
        # Theoretical uniform quantiles (i - 0.5) / n
        theo = (np.arange(1, n + 1) - 0.5) / n
        # Thin to ~300 points for file size / clarity while preserving shape
        step = max(1, n // 300)
        ax.plot(theo[::step], p[::step],
                color=col, lw=0.8, label=label)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    clean_axes(ax,
               xlabel="theoretical quantile",
               ylabel="empirical $p$-value")
    ax.legend(loc="lower right", handlelength=1.6)
    panel_letter(ax, "a", x=-0.22)


# ------------------------------------------------------------------
# panel b — empirical FDR vs nominal α
# ------------------------------------------------------------------
def panel_b(ax):
    files = [
        ("Dube PT1-PR1",  REPO / "out/robustness/module_de_null_PT1-PR1.tsv",      20, PRIMARY),
        ("Dube PT2-PR2",  REPO / "out/robustness/module_de_null_PT2-PR2.tsv",      20, ACCENT),
        ("Gisby late-early", REPO / "out_gisby/robustness/module_de_null_late-early.tsv", 15, SECOND),
    ]
    alphas = np.linspace(0.005, 0.20, 40)
    ax.plot(alphas, alphas, color=MUTED, lw=0.5, ls="--", zorder=0,
            label=r"nominal ($y = x$)")
    for label, path, K, col in files:
        df = pd.read_csv(path, sep="\t")
        B = df["perm"].nunique()
        # Empirical per-test FPR = mean rejections across perms / K
        emp = [
            (df[df["bh_q"] < a].groupby("perm").size().reindex(range(B), fill_value=0).mean()) / K
            for a in alphas
        ]
        ax.plot(alphas, emp, color=col, lw=0.9, label=label)
    ax.set_xlim(0, 0.20); ax.set_ylim(0, 0.20)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([0, 0.05, 0.10, 0.15, 0.20])
    ax.set_yticks([0, 0.05, 0.10, 0.15, 0.20])
    clean_axes(ax,
               xlabel=r"nominal BH-FDR  $\alpha$",
               ylabel="empirical per-module rejection rate")
    # Skip the legend here — it's redundant with panel a
    panel_letter(ax, "b", x=-0.22)


# ------------------------------------------------------------------
# panel c — discovery power vs n, by ranking rule
# ------------------------------------------------------------------
def panel_c(ax):
    sim = pd.read_csv(REPO / "out/simulations/s_score_regime_map.tsv", sep="\t")
    sub = sim[(sim.frac_tp == 0.15)
              & (sim.effect == "mixed")
              & (sim.hetero == False)
              & (sim.fragile == 0.0)].copy()
    g = sub.groupby(["n", "rule"])["recall_at_K"].agg(["mean", "sem"]).reset_index()

    rule_style = [
        ("q",     r"$q$",              LEGACY,  "o"),
        ("abs_d", r"$|\bar d|$",        ACCENT,  "s"),
        ("S",     r"$S$",              PRIMARY, "D"),
    ]
    for rule, label, col, mk in rule_style:
        r = g[g.rule == rule].sort_values("n")
        ax.errorbar(r.n, r["mean"], yerr=r["sem"],
                    fmt="-" + mk, color=col, mec=col, mfc=col,
                    markersize=3.2, lw=0.9, elinewidth=0.5, capsize=1.8,
                    label=label)
    ax.set_xticks([8, 12, 16, 20, 24])
    ax.set_ylim(0, 1.02)
    ax.set_box_aspect(1.0)
    clean_axes(ax,
               xlabel=r"subjects  $n$",
               ylabel="recall at $K$")
    ax.legend(loc="lower right", title="ranking by",
              title_fontsize=6.5, handlelength=1.6)
    panel_letter(ax, "c", x=-0.22)


# ------------------------------------------------------------------
# panel d — empirical protein-DE FPR at q<0.05 by method
# ------------------------------------------------------------------
def panel_d(ax):
    df = pd.read_csv(REPO / "benchmarks/out/null_fpr_summary.tsv", sep="\t")
    df = df[df["q_cutoff"] == 0.05].copy()
    # average/strip across the 4 Dube contrasts
    methods = [
        ("atman", "paired-t",  "Atman\npaired-t",  PRIMARY),
        ("atman", "moderated", "Atman\nmoderated", PRIMARY),
        ("atman", "welch-t",   "Atman\nWelch-t",   PRIMARY),
        ("limma", "eBayes",    "limma\neBayes",    SECOND),
    ]
    nominal = 0.05

    rng = np.random.default_rng(0)
    xs_all, ys_all, cs_all = [], [], []
    mean_xs, mean_ys = [], []
    for i, (tool, mode, _lab, col) in enumerate(methods):
        sub = df[(df.tool == tool) & (df["mode"] == mode)]
        vals = sub["empirical_fpr"].to_numpy()
        # replace exact zeros with a small floor so they show on log axis
        vals = np.where(vals > 0, vals, 5e-7)
        jitter = rng.uniform(-0.08, 0.08, size=len(vals))
        xs_all.extend(i + jitter)
        ys_all.extend(vals)
        cs_all.extend([col] * len(vals))
        mean_xs.append(i)
        mean_ys.append(np.mean(vals))

    # per-contrast dots
    ax.scatter(xs_all, ys_all, s=12, c=cs_all,
               edgecolors="white", linewidths=0.4, zorder=3)
    # mean marker: short horizontal bar
    for mx, my in zip(mean_xs, mean_ys):
        ax.plot([mx - 0.22, mx + 0.22], [my, my],
                color=STRUCTURE, lw=1.1, solid_capstyle="butt", zorder=4)

    # nominal reference
    ax.axhline(nominal, color=MUTED, lw=0.6, ls="--", zorder=1)
    ax.text(3.48, nominal * 1.15, r"nominal $q = 0.05$",
            ha="right", va="bottom", fontsize=6.3, color="#4B5563")

    ax.set_yscale("log")
    ax.set_ylim(3e-7, 0.2)
    ax.set_xlim(-0.5, 3.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels([m[2] for m in methods])
    ax.set_box_aspect(1.0)
    clean_axes(ax,
               xlabel=None,
               ylabel=r"empirical FPR at $q<0.05$")
    panel_letter(ax, "d", x=-0.22)


# ------------------------------------------------------------------
def main():
    fig, axes = plt.subplots(2, 2, figsize=GRID_2x2,
                             constrained_layout=True)
    panel_a(axes[0, 0])
    panel_b(axes[0, 1])
    panel_c(axes[1, 0])
    panel_d(axes[1, 1])
    save(fig, "fig2_calibration")
    plt.close(fig)


if __name__ == "__main__":
    main()
