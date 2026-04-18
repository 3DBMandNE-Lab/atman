#!/usr/bin/env python3
"""Floor 3 — null-FPR benchmark across Atman modes and limma.

For each of B sign-flip permutations, per contrast:
    - Sign-flip each subject's paired (A, B) labels with p=0.5
    - Run atman de with --test paired-t, moderated, welch-t
    - Run limma eBayes (via the companion R script)
    - Record rejection counts at nominal q ∈ {0.01, 0.05, 0.10}

OlinkAnalyze paired olink_ttest is mathematically equivalent to Atman
paired-t (paired Student, exact df, BH); we document that in the
supplementary and do not separately permute OlinkAnalyze here.

Output: benchmarks/out/null_fpr_results.tsv (one row per perm × contrast ×
tool × q_cutoff) and benchmarks/out/null_fpr_summary.tsv (aggregates).
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
ATMAN_BIN = REPO_ROOT / "target" / "release" / "atman"
Q_CUTOFFS = (0.01, 0.05, 0.10)


def sign_flip_samples(samples: pd.DataFrame, cond_a: str, cond_b: str,
                      rng: np.random.Generator) -> pd.DataFrame:
    """Swap (A, B) labels per subject independently with p=0.5.

    Operates on the `participant` column (atman's paired-by key). Caller
    must rename `subject_id` → `participant` upstream if the source file
    uses the former.
    """
    out = samples.copy()
    subset = out[out["condition"].isin([cond_a, cond_b])].copy()
    subjects = subset["participant"].unique()
    flip = {s: bool(rng.random() < 0.5) for s in subjects}
    mapped = subset.apply(
        lambda r: (cond_b if r["condition"] == cond_a else cond_a)
        if flip[r["participant"]] else r["condition"],
        axis=1,
    )
    out.loc[subset.index, "condition"] = mapped
    return out


def run_atman_de(tempdir: Path, contrast: str, mode: str) -> pd.DataFrame:
    out_dir = tempdir / f"de_{mode}"
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(ATMAN_BIN), "de",
        "--input-dir", str(tempdir),
        "--output-dir", str(out_dir),
        "--test", mode,
        "--groups", contrast,
        "--min-pairs", "5",
    ]
    if mode in ("paired-t", "moderated"):
        cmd += ["--paired-by", "participant"]
    subprocess.run(cmd, check=True, capture_output=True)
    return pd.read_csv(out_dir / "de_results.tsv", sep="\t")


def count_rejections(df: pd.DataFrame) -> dict[float, int]:
    q = pd.to_numeric(df["bh_q"], errors="coerce")
    return {c: int((q < c).sum()) for c in Q_CUTOFFS}


def run_atman_sweep(input_dir: Path, contrasts: list[str], modes: list[str],
                    n_perm: int, rng: np.random.Generator) -> pd.DataFrame:
    samples = pd.read_csv(input_dir / "samples.tsv", sep="\t")
    # atman's paired-by key is hardcoded to "participant" in v0.1.
    if "participant" not in samples.columns and "subject_id" in samples.columns:
        samples = samples.rename(columns={"subject_id": "participant"})
    rows = []
    with tempfile.TemporaryDirectory() as td:
        tempdir = Path(td)
        # Stage invariant inputs once.
        shutil.copy(input_dir / "qc_measurements.tsv", tempdir / "qc_measurements.tsv")
        shutil.copy(input_dir / "proteins.tsv", tempdir / "proteins.tsv")
        t0 = time.time()
        total = n_perm * len(contrasts) * len(modes)
        i = 0
        for contrast in contrasts:
            cond_a, cond_b = contrast.split("-")
            for p in range(n_perm):
                perm_samples = sign_flip_samples(samples, cond_a, cond_b, rng)
                perm_samples.to_csv(tempdir / "samples.tsv", sep="\t", index=False)
                for mode in modes:
                    de = run_atman_de(tempdir, contrast, mode)
                    counts = count_rejections(de)
                    for qc, nrej in counts.items():
                        rows.append({
                            "tool": "atman", "mode": mode, "contrast": contrast,
                            "perm": p, "q_cutoff": qc, "n_rejected": nrej,
                            "n_tested": int(de["bh_q"].notna().sum()),
                        })
                    i += 1
                if (p + 1) % 25 == 0:
                    el = time.time() - t0
                    print(f"[atman] {contrast} perm {p+1}/{n_perm}  elapsed {el:.1f}s",
                          file=sys.stderr)
    return pd.DataFrame(rows)


def run_limma_sweep(input_dir: Path, contrasts: list[str], n_perm: int,
                    seed: int) -> pd.DataFrame:
    """Shell out once to R; the R script does the perm loop internally
    using the same seed scheme."""
    r_script = REPO_ROOT / "benchmarks" / "null_fpr" / "run_null_fpr_limma.R"
    out_file = Path(tempfile.mkstemp(suffix=".tsv")[1])
    cmd = [
        "Rscript", str(r_script),
        "--input-dir", str(input_dir),
        "--contrasts", ",".join(contrasts),
        "--n-perm", str(n_perm),
        "--seed", str(seed),
        "--output", str(out_file),
    ]
    print(f"[limma] launching R with B={n_perm}", file=sys.stderr)
    t0 = time.time()
    subprocess.run(cmd, check=True)
    print(f"[limma] done in {time.time()-t0:.1f}s", file=sys.stderr)
    df = pd.read_csv(out_file, sep="\t")
    out_file.unlink()
    return df


def summarise(df: pd.DataFrame, n_tested: int) -> pd.DataFrame:
    """Per (tool, mode, contrast, q_cutoff) → mean FP count, empirical FPR."""
    g = df.groupby(["tool", "mode", "contrast", "q_cutoff"], as_index=False).agg(
        mean_n_rejected=("n_rejected", "mean"),
        median_n_rejected=("n_rejected", "median"),
        sd_n_rejected=("n_rejected", "std"),
        p95_n_rejected=("n_rejected", lambda x: float(np.percentile(x, 95))),
        n_perm=("perm", "count"),
    )
    # Empirical per-test FPR = mean_n_rejected / n_tested.
    g["empirical_fpr"] = g["mean_n_rejected"] / n_tested
    g["expected_fpr_if_uniform"] = g["q_cutoff"]
    return g


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, default=REPO_ROOT / "out")
    ap.add_argument("--contrasts", type=str,
                    default="PT1-PR1,PT2-PR2,PT2-PT1,PR2-PR1")
    ap.add_argument("--modes", type=str, default="paired-t,moderated,welch-t")
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--skip-limma", action="store_true",
                    help="Skip limma R subprocess (Atman-only sweep).")
    ap.add_argument("--output", type=Path,
                    default=REPO_ROOT / "benchmarks" / "out" / "null_fpr_results.tsv")
    ap.add_argument("--summary", type=Path,
                    default=REPO_ROOT / "benchmarks" / "out" / "null_fpr_summary.tsv")
    args = ap.parse_args()

    contrasts = args.contrasts.split(",")
    modes = args.modes.split(",")
    rng = np.random.default_rng(args.seed)

    atman_df = run_atman_sweep(args.input_dir, contrasts, modes, args.n_perm, rng)

    if not args.skip_limma:
        limma_df = run_limma_sweep(args.input_dir, contrasts, args.n_perm, args.seed)
        combined = pd.concat([atman_df, limma_df], ignore_index=True)
    else:
        combined = atman_df

    combined.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {args.output}  ({len(combined)} rows)", file=sys.stderr)

    # n_tested inferred from atman output (protein count ~ 2938).
    n_tested = int(atman_df["n_tested"].median())
    summary = summarise(combined, n_tested)
    summary.to_csv(args.summary, sep="\t", index=False)
    print(f"wrote {args.summary}  ({len(summary)} rows)", file=sys.stderr)

    # Headline print to stderr.
    print("\n=== Null-FPR headline (q<0.05) ===", file=sys.stderr)
    head = summary[summary["q_cutoff"] == 0.05][
        ["tool", "mode", "contrast", "mean_n_rejected", "empirical_fpr"]
    ].sort_values(["contrast", "tool", "mode"])
    print(head.to_string(index=False), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
