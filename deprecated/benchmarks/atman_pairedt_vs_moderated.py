#!/usr/bin/env python3
"""Within-Atman comparison: paired-t vs moderated-t on the Dube contrasts.

Reviewer request: show that switching inferential mode within Atman produces
consistent direction and rank signal, so the benchmark's paired-t-centred
comparison against OlinkAnalyze and limma is not an artefact of mode choice.

Reads:
    out/de_results.tsv       atman --test paired-t
    out_mod/de_results.tsv   atman --test moderated

Writes:
    benchmarks/out/atman_pairedt_vs_moderated.tsv

Per-contrast metrics:
    n_tested, n_common
    sign_concordance             fraction of common proteins with matching
                                  sign of mean_diff across both modes
    spearman_logfc                rank correlation of mean_diff
    paired_n_q_lt_005             paired-t hits at q<0.05
    moderated_n_q_lt_005          moderated hits at q<0.05
    hits_overlap_count            shared hits at q<0.05
    jaccard_q_lt_005              Jaccard of q<0.05 hit sets
    containment_paired_in_mod     paired hits recovered by moderated
    containment_mod_in_paired     moderated hits recovered by paired
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parent.parent
paired = pd.read_csv(REPO / "out" / "de_results.tsv", sep="\t")
moder = pd.read_csv(REPO / "out_mod" / "de_results.tsv", sep="\t")

CONTRASTS = ("PT1-PR1", "PT2-PR2", "PT2-PT1", "PR2-PR1")

rows = []
for c in CONTRASTS:
    p = paired[paired["comparison"] == c].set_index("assay_id")
    m = moder [moder ["comparison"] == c].set_index("assay_id")
    common_ids = p.index.intersection(m.index)
    cp = p.reindex(common_ids)
    cm = m.reindex(common_ids)
    # direction concordance
    sign_match = ((cp["mean_diff"] > 0) == (cm["mean_diff"] > 0)).sum()
    sign_conc  = sign_match / len(common_ids)
    # rank correlation of effect sizes
    rho, _ = spearmanr(cp["mean_diff"], cm["mean_diff"])
    # hits
    hits_p = set(cp.index[cp["bh_q"] < 0.05])
    hits_m = set(cm.index[cm["bh_q"] < 0.05])
    shared = hits_p & hits_m
    union  = hits_p | hits_m
    jac    = len(shared) / len(union) if union else float("nan")
    cont_p_in_m = len(shared) / len(hits_p) if hits_p else float("nan")
    cont_m_in_p = len(shared) / len(hits_m) if hits_m else float("nan")
    rows.append({
        "contrast":                     c,
        "n_common":                     len(common_ids),
        "sign_concordance":             round(sign_conc, 4),
        "spearman_logfc":               round(rho, 4),
        "paired_n_q_lt_005":            len(hits_p),
        "moderated_n_q_lt_005":         len(hits_m),
        "hits_overlap_count":           len(shared),
        "jaccard_q_lt_005":             round(jac, 4) if not pd.isna(jac) else float("nan"),
        "containment_paired_in_mod":    round(cont_p_in_m, 4) if not pd.isna(cont_p_in_m) else float("nan"),
        "containment_mod_in_paired":    round(cont_m_in_p, 4) if not pd.isna(cont_m_in_p) else float("nan"),
    })

out = pd.DataFrame(rows)
out.to_csv(REPO / "benchmarks" / "out" / "atman_pairedt_vs_moderated.tsv", sep="\t", index=False)
print(out.to_string(index=False))
