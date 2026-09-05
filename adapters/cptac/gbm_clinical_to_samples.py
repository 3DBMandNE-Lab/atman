#!/usr/bin/env python3
"""CPTAC GBM clinical ingest: Wang 2021 Table S1 + PDC biospecimen → atman samples.

Reads
  --samples          runs/GBM/samples.tsv          (canonical, aliquot-keyed)
  --proteins         runs/GBM/proteins.tsv
  --measurements     runs/GBM/measurements.tsv
  --table-s1         wang2021_tableS1_mmc2.xlsx    (sheets clinical_data, additional_annotations)
  --pdc-json         pdc_biospecimen_PDC000204.json
Writes
  extended samples.tsv (in place, + .run.json sidecar)
  --contrast-dir     runs/GBM_subtype  (condition = mesenchymal|other)
  stratum dirs       <contrast-dir stem>_mesenchymal / _other → runs/GBM_mesenchymal, runs/GBM_other
  --report           runs/GBM/clinical_ingest_report.tsv (+ sidecar)

Deterministic: sorted iteration everywhere, no timestamps in TSV bodies.
"""
import argparse, hashlib, json, os, subprocess, sys
from datetime import datetime, timezone

import openpyxl


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with open(path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]
    return rows[0], rows[1:]


def write_tsv(path, header, rows):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def sidecar(out_path, command, inputs, params):
    doc = {
        "atman_adapter": "gbm_clinical_to_samples.py",
        "command": command,
        "git_commit": subprocess.run(
            ["git", "-C", os.path.dirname(os.path.abspath(__file__)), "rev-parse", "HEAD"],
            capture_output=True, text=True).stdout.strip(),
        "wall_clock_utc": datetime.now(timezone.utc).isoformat(),
        "inputs_sha256": {os.path.basename(p): "sha256:" + sha256(p) for p in inputs},
        "outputs_sha256": {os.path.basename(out_path): "sha256:" + sha256(out_path)},
        "args": params,
    }
    with open(out_path + ".run.json", "w") as f:
        json.dump(doc, f, indent=2, sort_keys=True)
        f.write("\n")


def load_s1(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    def sheet_rows(name):
        rows = [[("" if v is None else str(v)) for v in r]
                for r in wb[name].iter_rows(values_only=True)]
        hdr, data = rows[0], rows[1:]
        return [dict(zip(hdr, r)) for r in data if any(c != "" for c in r)]
    clin = {r["case_id"]: r for r in sheet_rows("clinical_data")}
    ann = {r["case"]: r for r in sheet_rows("additional_annotations")}
    return clin, ann


APPEND_COLS = ["case_id", "wang_subtype", "multiomic_subtype", "idh1_alter_status",
               "mgmt_methyl_stp27_prediction", "is_gcimp", "xcell_immune_score",
               "xcell_stroma_score", "xcell_microenvironment_score",
               "s1_sample_type", "age", "gender"]


def main():
    ap = argparse.ArgumentParser()
    for flag in ["samples", "proteins", "measurements", "table_s1", "pdc_json",
                 "contrast_dir", "report"]:
        ap.add_argument("--" + flag.replace("_", "-"), required=True)
    a = ap.parse_args()

    pdc = json.load(open(a.pdc_json))["data"]["biospecimenPerStudy"]
    aliquot_to_case = {}
    for r in sorted(pdc, key=lambda r: r["aliquot_submitter_id"]):
        aliquot_to_case.setdefault(r["aliquot_submitter_id"], r["case_submitter_id"])

    clin, ann = load_s1(a.table_s1)
    hdr, samples = read_tsv(a.samples)
    if hdr[-len(APPEND_COLS):] == APPEND_COLS:
        sys.exit("samples.tsv already extended; refusing to double-append")

    def annotate(row):
        sid = row[hdr.index("sample_id")]
        case = aliquot_to_case.get(sid, "NA")
        c, n = clin.get(case, {}), ann.get(case, {})
        sub = n.get("rna_wang_cancer_cell_2017", "NA") or "NA"
        return row + [case, sub, n.get("multiomic", "NA") or "NA",
                      n.get("IDH1_alter_status", "NA") or "NA",
                      n.get("mgmt_methyl_stp27_prediction", "NA") or "NA",
                      n.get("is_gcimp", "NA") or "NA",
                      n.get("xcell_immune_score", "NA") or "NA",
                      n.get("xcell_stroma_score", "NA") or "NA",
                      n.get("xcell_microenvironment_score", "NA") or "NA",
                      n.get("sample_type", "NA") or "NA",
                      c.get("age", "NA") or "NA", c.get("gender", "NA") or "NA"]

    ext_hdr = hdr + APPEND_COLS
    ext = [annotate(r) for r in samples]
    write_tsv(a.samples, ext_hdr, ext)
    ins = [a.table_s1, a.pdc_json, a.proteins, a.measurements]
    sidecar(a.samples, "adapter gbm-clinical extend-samples", ins,
            {"append_cols": APPEND_COLS})

    # ── contrast + stratum dirs ────────────────────────────────────────────
    i_sid, i_cond = ext_hdr.index("sample_id"), ext_hdr.index("condition")
    i_sub = ext_hdr.index("wang_subtype")
    cond_of = {"Mesenchymal": "mesenchymal", "Classical": "other", "Proneural": "other"}
    kept, excluded = [], []
    for r in ext:
        m = cond_of.get(r[i_sub])
        if m is None:
            excluded.append((r[i_sid], r[i_sub]))
        else:
            r2 = list(r); r2[i_cond] = m; kept.append(r2)

    mhdr, meas = read_tsv(a.measurements)
    mi_sid = mhdr.index("sample_id")

    def emit(dirpath, rows):
        os.makedirs(dirpath, exist_ok=True)
        write_tsv(os.path.join(dirpath, "samples.tsv"), ext_hdr, rows)
        keep = {r[i_sid] for r in rows}
        write_tsv(os.path.join(dirpath, "measurements.tsv"), mhdr,
                  [m for m in meas if m[mi_sid] in keep])
        phdr, prot = read_tsv(a.proteins)
        write_tsv(os.path.join(dirpath, "proteins.tsv"), phdr, prot)
        sidecar(os.path.join(dirpath, "samples.tsv"),
                "adapter gbm-clinical contrast-dir", ins,
                {"n_samples": len(rows), "dir": os.path.basename(dirpath)})

    emit(a.contrast_dir, kept)
    base = os.path.dirname(os.path.normpath(a.contrast_dir))
    emit(os.path.join(base, "GBM_mesenchymal"),
         [r for r in kept if r[i_cond] == "mesenchymal"])
    emit(os.path.join(base, "GBM_other"), [r for r in kept if r[i_cond] == "other"])

    # ── report ─────────────────────────────────────────────────────────────
    from collections import Counter
    subs = Counter(r[i_sub] for r in ext)
    rep = [["metric", "value"],
           ["n_samples_total", str(len(ext))],
           ["n_mapped_to_case", str(sum(1 for r in ext if r[ext_hdr.index('case_id')] != 'NA'))],
           ["n_contrast_mesenchymal", str(sum(1 for r in kept if r[i_cond] == 'mesenchymal'))],
           ["n_contrast_other", str(sum(1 for r in kept if r[i_cond] == 'other'))]]
    rep += [["n_subtype_" + k.replace(" ", "_"), str(v)] for k, v in sorted(subs.items())]
    rep += [["excluded_" + s, sub] for s, sub in sorted(excluded)]
    write_tsv(a.report, rep[0], rep[1:])
    sidecar(a.report, "adapter gbm-clinical report", ins, {})
    print("clinical ingest done:", dict(rep[1:5]))


if __name__ == "__main__":
    main()
