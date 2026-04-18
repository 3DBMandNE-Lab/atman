#!/usr/bin/env python3
"""Aggregate per-tool DE outputs into Table 1 and an UpSet figure.

Reads `benchmarks/out/{atman,olinkanalyze,limma}_de.tsv`, each with
the shared schema (contrast, gene, mean_diff, pvalue, qvalue, tool,
runtime_s, peak_rss_mb). Emits:
    benchmarks/out/table1.tsv
    benchmarks/out/upset_de_overlap.pdf
"""

from __future__ import annotations

import itertools
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent / "out"
TOOLS = ("atman", "OlinkAnalyze_ttest", "limma")
TOOL_FILES = {
    "atman":      OUT / "atman_de.tsv",
    "OlinkAnalyze_ttest": OUT / "olinkanalyze_de.tsv",
    "limma":              OUT / "limma_de.tsv",
}
CONTRASTS = ("PT1-PR1", "PT2-PR2", "PT2-PT1", "PR2-PR1")


def load_tool(tool: str) -> pd.DataFrame:
    path = TOOL_FILES[tool]
    if not path.exists():
        print(f"[warn] {tool}: {path} missing; skipping")
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    df["tool"] = tool
    return df


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b) if a | b else float("nan")


def main() -> None:
    frames = {t: load_tool(t) for t in TOOLS}
    frames = {t: f for t, f in frames.items() if not f.empty}

    # ---- Table 1: DE counts + runtime + overlap ---------------------------
    rows = []
    hitsets: dict[tuple[str, str], set[str]] = {}
    for tool, df in frames.items():
        runtime = float(df["runtime_s"].iloc[0])
        for c in CONTRASTS:
            sub = df[df["contrast"] == c]
            hits05 = set(sub.loc[sub["qvalue"] < 0.05, "gene"])
            hits10 = set(sub.loc[sub["qvalue"] < 0.10, "gene"])
            hitsets[(tool, c)] = hits10
            rows.append({
                "tool":           tool,
                "contrast":       c,
                "n_tested":       len(sub),
                "n_q_lt_005":     len(hits05),
                "n_q_lt_010":     len(hits10),
                "runtime_s":      runtime,
            })
    t1 = pd.DataFrame(rows)
    t1.to_csv(OUT / "table1_counts.tsv", sep="\t", index=False)
    print("wrote table1_counts.tsv")

    # Pairwise Jaccard per contrast, q<0.10.
    jrows = []
    for c in CONTRASTS:
        for a, b in itertools.combinations(frames.keys(), 2):
            jrows.append({
                "contrast": c,
                "tool_a":   a,
                "tool_b":   b,
                "jaccard":  jaccard(hitsets[(a, c)], hitsets[(b, c)]),
                "n_a":      len(hitsets[(a, c)]),
                "n_b":      len(hitsets[(b, c)]),
                "n_shared": len(hitsets[(a, c)] & hitsets[(b, c)]),
            })
    jt = pd.DataFrame(jrows)
    jt.to_csv(OUT / "table1_overlap.tsv", sep="\t", index=False)
    print("wrote table1_overlap.tsv")

    # Combined Table 1 — wide-form for manuscript insertion.
    wide = t1.pivot(index="tool",
                    columns="contrast",
                    values=["n_q_lt_005", "n_q_lt_010"])
    wide.columns = [f"{c}_{m}" for m, c in wide.columns]
    runtimes = t1.drop_duplicates("tool").set_index("tool")["runtime_s"]
    wide["runtime_s"] = runtimes
    wide.to_csv(OUT / "table1.tsv", sep="\t")
    print("wrote table1.tsv")

    # ---- UpSet figure over PT1-PR1 DE hits --------------------------------
    try:
        from upsetplot import from_contents, UpSet
        import matplotlib.pyplot as plt
    except ImportError:
        print("[warn] upsetplot / matplotlib not installed — skipping figure")
        return

    contents = {t: hitsets[(t, "PT1-PR1")] for t in frames.keys()}
    if all(len(s) == 0 for s in contents.values()):
        print("[warn] all hit sets empty for PT1-PR1 — skipping figure")
        return

    data = from_contents(contents)
    fig = plt.figure(figsize=(7, 4))
    UpSet(data, subset_size="count", show_counts=True).plot(fig=fig)
    fig.suptitle("DE hit overlap — PT1-PR1, q<0.10", fontsize=10)
    fig.savefig(OUT / "upset_de_overlap.pdf", bbox_inches="tight", dpi=300)
    print("wrote upset_de_overlap.pdf")


if __name__ == "__main__":
    main()
