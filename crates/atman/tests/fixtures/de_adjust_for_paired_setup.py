#!/usr/bin/env python3
# crates/atman/tests/fixtures/de_adjust_for_paired_setup.py
#
# Generates the paired-design fixture for K5 (mixed-model) and K7 (paired-t)
# --adjust-for parity tests.
#
# Design:
#   12 subjects × 2 conditions (a, b) = 24 samples
#   Subject_k has samples  PA{k:02d} (condition a) and PB{k:02d} (condition b)
#   30 proteins, 2 numeric NMF covariates per sample
#
#   Abundance model:
#     abundance[g, i] = BASELINE
#                     + beta_g * condition_indicator_i
#                     + u_g_subject_i               (subject random intercept, SD 0.8)
#                     + gamma_g * nmf1_i
#                     + delta_g * nmf2_i
#                     + noise[g, i]                 (SD 0.3)
#
#   - Proteins 1-10:  beta_g = 2.0  (true DE)
#   - Proteins 11-20: gamma_g = 3.0 (nmf1-driven)
#   - Proteins 21-30: delta_g = 2.5 (nmf2-driven)
#
# Seeds are pinned so output is deterministic across Python versions.
# Usage:
#   python3 de_adjust_for_paired_setup.py <paired_input_dir> <paired_covariates_tsv>

import sys
import os
import math
import random

SEED_NOISE    = 20260501
SEED_NMF      = 20260502
SEED_SUBJECT  = 20260503

N_SUBJECTS  = 12
N_PROTEINS  = 30
BASELINE    = 8.0
SIGMA_NOISE = 0.30
SIGMA_SUBJ  = 0.80     # per-subject random intercept SD

def gauss(rng, mu, sigma):
    """Box-Muller normal deviate (no dependency on random.gauss internals)."""
    while True:
        u = rng.random()
        v = rng.random()
        if u > 0.0:
            break
    z = math.sqrt(-2.0 * math.log(u)) * math.cos(2.0 * math.pi * v)
    return mu + sigma * z

def write_tsv(path, header, rows, comment=None):
    with open(path, "w") as f:
        if comment:
            f.write("# " + comment + "\n")
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")

def main():
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} <paired_input_dir> <paired_covariates_tsv>",
              file=sys.stderr)
        sys.exit(1)

    paired_dir    = sys.argv[1]
    cov_tsv       = sys.argv[2]
    os.makedirs(paired_dir, exist_ok=True)

    rng_nmf     = random.Random(SEED_NMF)
    rng_subj    = random.Random(SEED_SUBJECT)
    rng_noise   = random.Random(SEED_NOISE)

    # ── subject random intercepts (per protein) ──────────────────────────────
    # u[g][k] ~ N(0, SIGMA_SUBJ)  — k = 0..N_SUBJECTS-1
    u = [[gauss(rng_subj, 0.0, SIGMA_SUBJ) for _ in range(N_SUBJECTS)]
         for _ in range(N_PROTEINS)]

    # ── NMF covariates ────────────────────────────────────────────────────────
    # 24 samples: PA01..PA12 then PB01..PB12
    # nmf values are independent of condition
    n_samples = N_SUBJECTS * 2
    nmf1 = [rng_nmf.random() for _ in range(n_samples)]
    nmf2 = [rng_nmf.random() for _ in range(n_samples)]

    # ── protein effect coefficients ───────────────────────────────────────────
    betas  = [0.0] * N_PROTEINS
    gammas = [0.0] * N_PROTEINS
    deltas = [0.0] * N_PROTEINS
    for g in range(N_PROTEINS):
        if g < 10:
            betas[g]  = 2.0
        elif g < 20:
            gammas[g] = 3.0
        else:
            deltas[g] = 2.5

    # ── build sample list ─────────────────────────────────────────────────────
    # First N_SUBJECTS samples: condition a (PA01..PA12)
    # Next N_SUBJECTS samples: condition b (PB01..PB12)
    # sample index i = 0..N_SUBJECTS-1 → condition a, subject k=i
    # sample index i = N_SUBJECTS..2*N_SUBJECTS-1 → condition b, subject k=i-N_SUBJECTS
    sample_ids   = ([f"PA{k+1:02d}" for k in range(N_SUBJECTS)] +
                    [f"PB{k+1:02d}" for k in range(N_SUBJECTS)])
    subject_ids  = ([f"SUB{k+1:02d}" for k in range(N_SUBJECTS)] +
                    [f"SUB{k+1:02d}" for k in range(N_SUBJECTS)])
    conditions   = ["a"] * N_SUBJECTS + ["b"] * N_SUBJECTS

    # ── abundance matrix ──────────────────────────────────────────────────────
    abundance = []   # [protein][sample]
    for g in range(N_PROTEINS):
        row = []
        for i in range(n_samples):
            cond_ind = 1.0 if conditions[i] == "b" else 0.0
            k = i if i < N_SUBJECTS else i - N_SUBJECTS  # subject index
            mu = (BASELINE
                  + betas[g] * cond_ind
                  + u[g][k]
                  + gammas[g] * nmf1[i]
                  + deltas[g] * nmf2[i])
            row.append(gauss(rng_noise, mu, SIGMA_NOISE))
        abundance.append(row)

    # ── write samples.tsv ─────────────────────────────────────────────────────
    samples_header = ["sample_id", "subject_id", "condition", "is_control",
                      "sample_type", "ingest_order"]
    samples_rows = []
    for i, (sid, subj, cond) in enumerate(zip(sample_ids, subject_ids, conditions)):
        samples_rows.append([sid, subj, cond, 0, "plasma", i + 1])
    write_tsv(os.path.join(paired_dir, "samples.tsv"), samples_header, samples_rows)

    # ── write proteins.tsv ────────────────────────────────────────────────────
    proteins_header = ["platform", "assay_id", "uniprot", "gene_symbol", "panel", "panel_lot"]
    proteins_rows = []
    for g in range(N_PROTEINS):
        proteins_rows.append([
            "olink_explore_ngs",
            f"A{g+1:03d}",
            f"Q{g+1:05d}",
            f"PROT{g+1:03d}",
            "P1",
            "",
        ])
    write_tsv(os.path.join(paired_dir, "proteins.tsv"), proteins_header, proteins_rows)

    # ── write measurements.tsv ────────────────────────────────────────────────
    meas_header = [
        "platform", "sample_id", "assay_id", "gene_symbol", "panel",
        "npx_source_str", "abundance", "abundance_raw", "abundance_unit",
        "qc_sample", "qc_assay", "detection_limit", "below_lod",
        "dropped_by_qc", "plate_id", "panel_lot", "ingest_order",
    ]
    meas_rows = []
    order = 0
    for g in range(N_PROTEINS):
        for i, sid in enumerate(sample_ids):
            order += 1
            v = abundance[g][i]
            meas_rows.append([
                "olink_explore_ngs",
                sid,
                f"A{g+1:03d}",
                f"PROT{g+1:03d}",
                "P1",
                f"{v:.6f}",
                f"{v:.6f}",
                f"{v:.6f}",
                "log2_npx",
                "PASS",
                "PASS",
                "",
                0,
                0,
                "",
                "",
                order,
            ])
    write_tsv(os.path.join(paired_dir, "measurements.tsv"), meas_header, meas_rows)

    # ── write covariates.tsv ──────────────────────────────────────────────────
    cov_header = ["sample_id", "nmf_program_01", "nmf_program_02"]
    cov_rows = []
    for i, sid in enumerate(sample_ids):
        cov_rows.append([sid, f"{nmf1[i]:.8f}", f"{nmf2[i]:.8f}"])
    import sys as _sys
    comment = (f"generated by de_adjust_for_paired_setup.py "
               f"Python {_sys.version.split()[0]}")
    write_tsv(cov_tsv, cov_header, cov_rows, comment=comment)

    print(f"wrote {paired_dir}/samples.tsv ({n_samples} samples, {N_SUBJECTS} subjects × 2 conditions)")
    print(f"wrote {paired_dir}/proteins.tsv ({N_PROTEINS} proteins)")
    print(f"wrote {paired_dir}/measurements.tsv ({order} rows)")
    print(f"wrote {cov_tsv} ({n_samples} rows)")

if __name__ == "__main__":
    main()
