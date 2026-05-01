"""
Generate bench/planted_archetypes_nmf_v1/ fixture for NMF benchmarking.

Design
------
- 3 planted non-negative archetypes
- 80 features x 120 samples
- Sample x archetype mixing weights: Gamma(2.0, 1.0)  [shape=2, scale=1]
- Archetype x feature loadings: Gamma(1.0, 1.0), localised to ~1/3 of
  features per archetype (sparse biological-like signatures)
- Additive non-negative noise: Gamma(0.5, 0.2)
- All values non-negative by construction

Outputs (all in this directory):
  abundance.tsv          -- sample_id x gene_symbol observed matrix
  planted_loadings.tsv   -- archetype_id / protein / loading triples
  planted_activations.tsv -- sample_id / program / activation triples

Usage: python3 generate.py

Python >= 3.8 required. Stdlib + numpy only.
numpy version pinned implicitly by environment; recorded in header comments.
"""

import csv
import os
import sys

# --- version metadata recorded in comments ---
import numpy as np

RNG_SEED = 20260429
N_ARCHETYPES = 3
N_FEATURES = 80
N_SAMPLES = 120
# Each archetype localised to a contiguous block of ~1/3 of features.
BLOCK_SIZE = N_FEATURES // N_ARCHETYPES  # 26 per archetype; remainder shared

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    rng = np.random.default_rng(RNG_SEED)

    # ------------------------------------------------------------------ #
    # 1.  Ground-truth loadings: H_true  [K x P]
    # ------------------------------------------------------------------ #
    # Each archetype has a dominant block of features (Gamma-distributed)
    # plus near-zero background outside its block.
    H_true = np.zeros((N_ARCHETYPES, N_FEATURES))
    for k in range(N_ARCHETYPES):
        start = k * BLOCK_SIZE
        end = start + BLOCK_SIZE
        # Strong loadings in the own block.
        H_true[k, start:end] = rng.gamma(shape=1.0, scale=1.0, size=BLOCK_SIZE)
        # Weak loadings outside the block (sparse leakage, biologically plausible).
        other_size = N_FEATURES - BLOCK_SIZE
        leakage = rng.gamma(shape=0.3, scale=0.1, size=other_size)
        mask = np.ones(N_FEATURES, dtype=bool)
        mask[start:end] = False
        H_true[k, mask] = leakage

    # ------------------------------------------------------------------ #
    # 2.  Sample activations: W_true  [N x K]
    # ------------------------------------------------------------------ #
    W_true = rng.gamma(shape=2.0, scale=1.0, size=(N_SAMPLES, N_ARCHETYPES))

    # ------------------------------------------------------------------ #
    # 3.  Observed matrix: X = W H + noise  [N x P], all non-negative
    # ------------------------------------------------------------------ #
    X_signal = W_true @ H_true  # [N x P]
    noise = rng.gamma(shape=0.5, scale=0.2, size=(N_SAMPLES, N_FEATURES))
    X = X_signal + noise
    # Verify non-negativity.
    assert X.min() >= 0.0, "BUG: negative values in generated matrix"

    # ------------------------------------------------------------------ #
    # 4.  Labels
    # ------------------------------------------------------------------ #
    sample_ids = [f"S{i:03d}" for i in range(N_SAMPLES)]
    gene_ids = [f"G{j:03d}" for j in range(N_FEATURES)]
    archetype_ids = [f"P_A{k+1:02d}" for k in range(N_ARCHETYPES)]

    # ------------------------------------------------------------------ #
    # 5.  Write abundance.tsv
    # ------------------------------------------------------------------ #
    abundance_path = os.path.join(SCRIPT_DIR, "abundance.tsv")
    with open(abundance_path, "w", newline="") as f:
        header = ["sample_id"] + gene_ids
        f.write("\t".join(header) + "\n")
        for i, sid in enumerate(sample_ids):
            row = [sid] + [f"{X[i, j]:.6f}" for j in range(N_FEATURES)]
            f.write("\t".join(row) + "\n")
    print(f"wrote {abundance_path}  ({N_SAMPLES} samples x {N_FEATURES} features)")

    # ------------------------------------------------------------------ #
    # 6.  Write planted_loadings.tsv
    # ------------------------------------------------------------------ #
    loadings_path = os.path.join(SCRIPT_DIR, "planted_loadings.tsv")
    with open(loadings_path, "w", newline="") as f:
        f.write("archetype_id\tprotein\tloading\n")
        for k, aid in enumerate(archetype_ids):
            for j, gid in enumerate(gene_ids):
                f.write(f"{aid}\t{gid}\t{H_true[k, j]:.6f}\n")
    print(f"wrote {loadings_path}  ({N_ARCHETYPES} archetypes x {N_FEATURES} features)")

    # ------------------------------------------------------------------ #
    # 7.  Write planted_activations.tsv
    # ------------------------------------------------------------------ #
    activations_path = os.path.join(SCRIPT_DIR, "planted_activations.tsv")
    with open(activations_path, "w", newline="") as f:
        f.write("sample_id\tprogram\tactivation\n")
        for i, sid in enumerate(sample_ids):
            for k, aid in enumerate(archetype_ids):
                f.write(f"{sid}\t{aid}\t{W_true[i, k]:.6f}\n")
    print(
        f"wrote {activations_path}  ({N_SAMPLES} samples x {N_ARCHETYPES} programs)"
    )

    print(
        f"\nFixture summary:\n"
        f"  RNG seed:      {RNG_SEED}\n"
        f"  Archetypes:    {N_ARCHETYPES}\n"
        f"  Features:      {N_FEATURES}\n"
        f"  Samples:       {N_SAMPLES}\n"
        f"  X min:         {X.min():.4f}\n"
        f"  X max:         {X.max():.4f}\n"
        f"  numpy version: {np.__version__}\n"
        f"  python:        {sys.version}\n"
    )


if __name__ == "__main__":
    main()
