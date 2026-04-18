"""Nature Methods / Bioinformatics figure style.

Conventions copied from the canonical figure vocabulary:
  - Helvetica throughout (fallback Arial, DejaVu Sans).
  - Typography hierarchy: 6pt ticks, 7pt axis labels, 8pt annotations,
    9pt bold lowercase panel letters.
  - 0.5pt spines and ticks (outward, 2.5pt long). No grid.
  - Muted two-colour contrasts. Colour only carries information.
  - No legend frames. No background washes. No rounded corners.

Figure sizing follows OUP Bioinformatics column widths:
  SINGLE_COL = 85 mm  (3.35 in)
  DOUBLE_COL = 180 mm (7.09 in)
"""
from __future__ import annotations
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt


# ---- palette (muted; used only where information requires it) -----------
PRIMARY    = "#1F4E79"   # Oxford blue — Atman, bootstrap, method
ACCENT     = "#C06B2A"   # burnt orange — OlinkAnalyze, comparator 1
SECOND     = "#6E9570"   # sage          — limma, comparator 2
LEGACY     = "#7A7A7A"   # mid grey      — legacy S / SSR, demoted
STRUCTURE  = "#1A1A1A"   # near-black    — spines, axis labels, text
MUTED      = "#9CA3AF"   # cool grey     — reference lines, no-signal
NEG        = "#A63A50"   # deep red      — used only when emphasising a negative
DIVERGING  = "RdBu_r"
SEQUENTIAL = "viridis"

# Per-tool colours for benchmark plots.
TOOL_COLORS = {
    "atman":              PRIMARY,
    "OlinkAnalyze_ttest": ACCENT,
    "limma":              SECOND,
}
TOOL_LABELS = {
    "atman":              "Atman",
    "OlinkAnalyze_ttest": "OlinkAnalyze",
    "limma":              "limma",
}

CONTRAST_COLORS = {
    "PT1-PR1":    PRIMARY,
    "PT2-PR2":    ACCENT,
    "PT2-PT1":    LEGACY,
    "PR2-PR1":    SECOND,
    "late-early": PRIMARY,
}

# ---- canvas sizes (inches) ----------------------------------------------
SINGLE_COL = 3.35        # 85 mm
DOUBLE_COL = 7.09        # 180 mm
PANEL_1    = (SINGLE_COL, 2.6)
PANEL_SQ   = (2.8, 2.8)
ROW_3      = (DOUBLE_COL, 2.4)     # three-panel row at double-column width
ROW_4      = (DOUBLE_COL, 2.1)     # four-panel row at double-column width
GRID_2x2   = (DOUBLE_COL, 4.6)
GRID_2x3   = (DOUBLE_COL, 4.4)


def apply_style():
    mpl.rcParams.update({
        # typography
        "font.family":         ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size":            7,
        "axes.titlesize":       8,
        "axes.labelsize":       7,
        "xtick.labelsize":      6,
        "ytick.labelsize":      6,
        "legend.fontsize":      6.5,
        "axes.titlepad":        4,
        "axes.labelpad":        2.5,
        # spines
        "axes.spines.top":      False,
        "axes.spines.right":    False,
        "axes.linewidth":       0.5,
        "axes.edgecolor":       STRUCTURE,
        "axes.labelcolor":      STRUCTURE,
        "text.color":           STRUCTURE,
        # ticks
        "xtick.direction":      "out",
        "ytick.direction":      "out",
        "xtick.major.size":     2.5,
        "ytick.major.size":     2.5,
        "xtick.major.width":    0.5,
        "ytick.major.width":    0.5,
        "xtick.minor.size":     1.3,
        "ytick.minor.size":     1.3,
        "xtick.minor.width":    0.4,
        "ytick.minor.width":    0.4,
        "xtick.color":          STRUCTURE,
        "ytick.color":          STRUCTURE,
        "xtick.major.pad":      2.5,
        "ytick.major.pad":      2.5,
        # grid — always off unless a panel explicitly enables it
        "axes.grid":            False,
        # legend
        "legend.frameon":       False,
        "legend.handlelength":  1.6,
        "legend.handletextpad": 0.5,
        "legend.borderpad":     0.2,
        "legend.columnspacing": 1.0,
        # lines / markers
        "lines.linewidth":      0.8,
        "lines.markersize":     3.0,
        "patch.linewidth":      0.5,
        # figure
        "figure.facecolor":     "white",
        "figure.dpi":           150,
        "savefig.dpi":          400,
        "savefig.bbox":         "tight",
        "savefig.pad_inches":   0.02,
        "pdf.fonttype":         42,
        "ps.fonttype":          42,
    })


REPO = Path(__file__).resolve().parents[2]
OUT  = REPO / "docs" / "findings" / "figures"


def save(fig, stem: str):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png")
    print(f"wrote {OUT / stem}.{{pdf,png}}")


def panel_letter(ax, letter: str, *, x=-0.18, y=1.08):
    """Bold lowercase panel letter in the top-left of an axes, outside the
    plotting area. Positions are in axes-fraction coordinates so they
    align consistently across panels regardless of axis limits."""
    ax.text(x, y, letter, transform=ax.transAxes,
            ha="left", va="bottom",
            fontsize=9, fontweight="bold",
            color=STRUCTURE)


def clean_axes(ax, *, xlabel=None, ylabel=None):
    """Apply final per-panel polish."""
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("bottom", "left"):
        ax.spines[sp].set_linewidth(0.5)
        ax.spines[sp].set_color(STRUCTURE)
