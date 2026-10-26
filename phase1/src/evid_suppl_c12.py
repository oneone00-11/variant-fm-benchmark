"""E20 C12 -- the in-frame attribution of the true positives, beside the false positives'.

E6 (evid_inframe) asked how many of a tool's false positives -- variants it scores at or
above Walker's PP3 cut point that the assay calls normal -- are predicted to change the
transcript in frame, and E2.7 asked the same on DDX3X and TP53. This gives the true
positives the same attribution -- the same cut point, the assay calling the variant
damaging -- so that the false positives' share stands beside a reference.

  * True positives: the event records of E20.1 (scripts/stage_E20/run_c12_events.py, E6's
    own scorers in E6's environments), on the 3-10 and 11-50 nt bands (DDX3X's variants
    reach 25 nt, TP53's 12 nt); the anchors those runs carried are left out. Every
    variant is classified by evid_inframe.classify (SpliceAI) or classify_pangolin
    (Pangolin), unchanged, on evid_inframe's transcript maps. DDX3X and TP53 are read on
    SpliceAI's events only; each label definition is counted on its own, as in E2.7.
  * False positives: E6's and E2.7's per-variant attribution files as they stand
    (inframe_attribution.csv, inframe_attribution_ddx3x.csv, inframe_attribution_tp53.csv),
    not recomputed, counted the way their summaries count them; every count must equal
    the summary's, cell by cell, or the stage stops. E2.7's DDX3X Pangolin rows, read on
    SpliceAI's false positives, have no true-positive counterpart and are not carried.

Per gene set x label definition x tool x band x group (true_positive, false_positive):
in-frame, out-of-frame and undetermined counts, n, the in-frame share of n with E15's
Wilson 95% interval, and the share among the resolved, in-frame / (in-frame +
out-of-frame); shares rounded to four decimals, as E6's summaries round them.

Output, in reports/evidence/suppl_analyses/, with a header (lines beginning '#'):
  c12_inframe_tp_fp.csv

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_c12
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import evid_common as K
from . import evid_inframe as F
from . import evid_practical_metrics as P

OUT_DIR = Path("reports/evidence/suppl_analyses")
OUT = OUT_DIR / "c12_inframe_tp_fp.csv"
REPORT_DIR = Path("reports/evidence")
EV = Path("data/evidence")
BANDS = ["s3_10", "s11_50"]
FRAMES = ["in_frame", "out_of_frame", "undetermined"]
SEVEN = "seven genes"
LABEL_SEVEN = "y_assay"
# gene set: (tool, subset, events, classifier, label columns or None)
TRUE_POSITIVES = [
    (SEVEN, "spliceai", EV / "inframe_subset_tp.parquet", EV / "inframe_events_tp.parquet",
     F.classify, None),
    (SEVEN, "pangolin", EV / "inframe_subset_pangolin_tp.parquet",
     EV / "inframe_events_pangolin_tp.parquet", F.classify_pangolin, None),
    ("DDX3X", "spliceai", EV / "external/ddx3x_inframe_subset_tp.parquet",
     EV / "external/ddx3x_inframe_events_tp.parquet", F.classify, F.EXTERNAL["ddx3x"]["label_cols"]),
    ("TP53", "spliceai", EV / "external/tp53_inframe_subset_12nt_tp.parquet",
     EV / "external/tp53_inframe_events_12nt_tp.parquet", F.classify, F.EXTERNAL["tp53"]["label_cols"]),
]
# gene set: (per-variant file, summary, label columns or None)
FALSE_POSITIVES = [
    (SEVEN, REPORT_DIR / "inframe_attribution.csv", REPORT_DIR / "inframe_attribution_summary.csv", None),
    ("DDX3X", REPORT_DIR / "inframe_attribution_ddx3x.csv",
     REPORT_DIR / "inframe_attribution_ddx3x_summary.csv", F.EXTERNAL["ddx3x"]["label_cols"]),
    ("TP53", REPORT_DIR / "inframe_attribution_tp53.csv",
     REPORT_DIR / "inframe_attribution_tp53_summary.csv", F.EXTERNAL["tp53"]["label_cols"]),
]
EXTERNAL_TOOLS = ["spliceai"]


def _count(frames: pd.DataFrame) -> dict:
    c = frames["frame"].value_counts()
    k = {f: int(c.get(f, 0)) for f in FRAMES}
    k["n"] = int(len(frames))
    return k


def _row(gene_set, label, tool, band, group, k) -> dict:
    n, inf, out = k["n"], k["in_frame"], k["out_of_frame"]
    lo, hi = P._wilson(inf, n)
    return {"gene_set": gene_set, "label_definition": label, "tool": tool, "stratum": band,
            "group": group, **{f: k[f] for f in FRAMES}, "n": n,
            "in_frame_share": round(inf / n, 4) if n else np.nan,
            "resolved_in_frame_share": round(inf / (inf + out), 4) if inf + out else np.nan,
            "in_frame_share_lo": round(lo, 4) if n else np.nan,
            "in_frame_share_hi": round(hi, 4) if n else np.nan}


def true_positives() -> list[dict]:
    parts = []
    for gene_set, tool, subset_path, events_path, classifier, labels in TRUE_POSITIVES:
        sub = pd.read_parquet(subset_path)
        ev = pd.read_parquet(events_path)
        ev = ev[ev["variant_id"].isin(sub["variant_id"])]
        if len(ev) != len(sub):
            raise SystemExit(f"[C12] {events_path}: {len(sub) - len(ev)} subset variants have "
                             "no event record: stopped")
        keep = ["variant_id", "stratum"] + (labels or [])
        ev = ev.merge(sub[keep], on="variant_id", how="left", validate="1:1")
        parts.append((gene_set, tool, ev, classifier, labels))
    genes = sorted({g for *_x, ev, _c, _l in parts for g in ev["gene"]} & set(F.GENE_TRANSCRIPT))
    maps = F._transcript_maps(genes)
    rows = []
    for gene_set, tool, ev, classifier, labels in parts:
        attributed = []
        for r in ev.to_dict("records"):
            if r["gene"] not in maps:
                attributed.append(r | {"frame": "undetermined",
                                       "reason": "no transcript model for this gene"})
            else:
                attributed.append(r | classifier(r, maps[r["gene"]]))
        out = pd.DataFrame(attributed)
        for label in (labels or [LABEL_SEVEN]):
            o = out if labels is None else out[out[label] == 1]
            for band in BANDS:
                rows.append(_row(gene_set, label, tool, band, "true_positive",
                                 _count(o[o["stratum"] == band])))
    return rows


def false_positives() -> list[dict]:
    rows = []
    for gene_set, per_variant, summary_path, labels in FALSE_POSITIVES:
        att = K.read_back(per_variant)
        summ = K.read_back(summary_path)
        tools = sorted(att["tool"].unique()) if labels is None else EXTERNAL_TOOLS
        for label in (labels or [LABEL_SEVEN]):
            o = att if labels is None else att[att[label] == 0]
            s = summ if labels is None else summ[summ["label_definition"] == label]
            for tool in tools:
                for band in BANDS:
                    k = _count(o[(o["tool"] == tool) & (o["stratum"] == band)])
                    cell = s[(s["tool"] == tool) & (s["stratum"] == band)]
                    want = ({f: int(cell[f].iloc[0]) for f in FRAMES + ["n"]} if len(cell)
                            else {f: 0 for f in FRAMES + ["n"]})
                    if k != want:
                        raise SystemExit(f"[C12] {gene_set} {label} {tool} {band}: the false "
                                         f"positives counted from {per_variant} ({k}) differ "
                                         f"from {summary_path} ({want}): stopped")
                    rows.append(_row(gene_set, label, tool, band, "false_positive", k))
    return rows


def header() -> list[str]:
    return [
        "# E20 C12: the in-frame attribution of the true positives, beside the false positives' (src/evid_suppl_c12.py).",
        "# In frame or not is inferred from SpliceAI's or Pangolin's own output, not measured on transcripts: the dominant "
        "event's position read against the transcript's exon boundaries (an exon skipped, a cryptic site's shift from the "
        "nearest natural site); nonsense-mediated decay, partial use of a new site and a variant that changes more than "
        "one site are not considered (evid_inframe, E6).",
        "# True positives: scored at or above Walker's PP3 cut point, 0.2 (SpliceAI: spliceai_walker; Pangolin: the same "
        "cut on pangolin, so that the two tools' sets are defined alike), and called damaging by the assay (y_assay == 1; "
        "DDX3X and TP53 under each of their label definitions); the event records of E20.1, E6's scorers in E6's "
        "environments, the anchors left out.",
        "# False positives: E6's and E2.7's attribution files as they stand, counted as their summaries count them (equal "
        "in every cell). E6 drew its false positives on 2026-09-17, before E16d switched PALB2 and BARD1 to the archived "
        "classification: under the current labels 6 PALB2 variants at 3-10 nt (SpliceAI) and 7 (Pangolin) would be "
        "added; they are not in E6's files. The true positives are the same under either version of the labels.",
        "# TP53's label (y_tp53_control_anchored, as E2.7 derives it) is its primary label, urn:mavedb:00001277-a-1's "
        "investigator-provided ranges, variant by variant (checked in E18: 288 of 288).",
        "# Bands: 3-10 and 11-50 nt (DDX3X 11-25 nt, TP53 11-12 nt: their variants reach 25 and 12 nt); +/-1-2 nt is "
        "outside PP3 and BP4 and is not counted. DDX3X and TP53 are read on SpliceAI's events only; E2.7's DDX3X "
        "Pangolin rows have no true-positive counterpart and are not carried.",
        "# in_frame_share = in_frame / n, with its Wilson 95% interval (E15's); resolved_in_frame_share = in_frame / "
        "(in_frame + out_of_frame); rounded to four decimals, as E6's summaries are.",
    ]


def main() -> None:
    rows = true_positives() + false_positives()
    out = pd.DataFrame(rows)
    order = {g: i for i, g in enumerate([SEVEN, "DDX3X", "TP53"])}
    # label definitions in E2.7's order, the primary label first
    labels = [LABEL_SEVEN] + F.EXTERNAL["ddx3x"]["label_cols"] + F.EXTERNAL["tp53"]["label_cols"]
    out["_g"] = out["gene_set"].map(order)
    out["_l"] = out["label_definition"].map({lab: i for i, lab in enumerate(labels)})
    out["_b"] = out["stratum"].map({b: i for i, b in enumerate(BANDS)})
    out["_grp"] = out["group"].map({"true_positive": 0, "false_positive": 1})
    out = out.sort_values(["_g", "_l", "tool", "_b", "_grp"],
                          kind="mergesort").drop(columns=["_g", "_l", "_b", "_grp"])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        f.write("\n".join(header()) + "\n")
        out.to_csv(f, index=False)
    print(f"[C12] wrote {OUT} ({len(out)} rows)")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
