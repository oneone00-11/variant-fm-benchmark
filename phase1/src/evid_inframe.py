"""E6 -- what the tool's false positives are predicting.

A variant the tool scores above the PP3 cut point and the assay calls normal is not
necessarily a wrong splicing prediction. A predicted exon skip whose exon length is
a multiple of three, or a cryptic site shifted by a multiple of three, removes or
adds amino acids without shifting the reading frame, and a cell-survival assay can
score that as normal function while the splicing prediction was right. Separating
those from genuine false positives says something different about the tool than a
specificity figure does.

This is an INFERENCE from SpliceAI's own output, not a measurement of transcripts:

  * the dominant event is the largest of the four delta scores, and its DP field
    gives the position of the lost or gained site relative to the variant;
  * for an acceptor or donor LOSS, the position is matched to a natural exon
    boundary (within 2 bp) and the predicted consequence is taken to be skipping of
    that exon -- in-frame when the exon's length is a multiple of three;
  * for a GAIN, the position is a cryptic site and the consequence is taken to be
    the shift between it and the nearest natural site of the same kind -- in-frame
    when that shift is a multiple of three.

Both readings are the standard first-pass rule a PVS1 decision tree applies; neither
accounts for nonsense-mediated decay, for partial usage of the new site, or for a
variant that changes more than one site. Cells that cannot be resolved are reported
as `undetermined` rather than assigned.

Pangolin is NOT covered. The atlas persisted only its maximum score, its per-event
components and positions were never written, and its scorer is materially slower
than SpliceAI's; a re-score of the needed subset is a separate job. This is recorded
rather than worked around.

Run (PYTHONPATH=phase1):  python -m src.evid_inframe --subset   # write the variant list
                          python -m src.evid_inframe --attribute
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import config as C
from . import evid_common as K

ATLAS_REPO = Path(os.environ.get(
    "evidence-strength_ATLAS_REPO", "/Users/cliffzhang/work/functional-standard-atlas"))
evidence-strength_DIR = Path("data/evidence")
REPORT_DIR = Path("reports/evidence")
SUBSET = evidence-strength_DIR / "inframe_subset.parquet"
EVENTS = evidence-strength_DIR / "inframe_events.parquet"
OUT = REPORT_DIR / "inframe_attribution.csv"
CONFIG_PATH = Path("config/walker2023.yaml")

BOUNDARY_TOL = 2          # bp; a predicted loss must sit on a real boundary
EVENT_NAMES = {"ds_ag": "acceptor_gain", "ds_al": "acceptor_loss",
               "ds_dg": "donor_gain", "ds_dl": "donor_loss"}

# The assay transcript each gene's c. notation is written against; the atlas maps
# every gene to its MANE Select transcript, so these are the MANE RefSeq accessions.
GENE_TRANSCRIPT = {
    "BAP1": "NM_004656.4", "BARD1": "NM_000465.4", "BRCA1": "NM_007294.4",
    "BRCA2": "NM_000059.4", "PALB2": "NM_024675.4", "RAD51C": "NM_058216.3",
    "VHL": "NM_000551.4",
}


def write_subset() -> None:
    """The variants E6 is about: called by the tool, called normal by the assay."""
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    pp3 = cfg["thresholds"]["pp3"]["value"]
    df = K.load_set()
    col = "spliceai_walker" if "spliceai_walker" in df.columns else "spliceai"
    sub = df[(df[col] >= pp3) & (df["y_assay"] == 0)].copy()
    sub["score_column"] = col
    sub[["variant_id", "gene", "chrom", "pos", "ref", "alt", "stratum",
         "hgvs_c", col, "score_column"]].to_parquet(SUBSET, index=False)
    print(f"[E6] {len(sub):,} variants scored >= {pp3} on `{col}` and labelled normal")
    print(sub.groupby("stratum", observed=True).size().to_string())
    print(f"     wrote {SUBSET}")
    print("\n[E6] next, in the SpliceAI environment:")
    print("  <atlas>/models/spliceai/.venv/bin/python phase1/src/evid_score_spliceai_events.py \\")
    print(f"      --variants {SUBSET} --distance 4999 --out {EVENTS} --run")


# ---------------------------------------------------------------------------
# frame inference
# ---------------------------------------------------------------------------
def _transcript_maps(genes: list[str]) -> dict:
    sys.path.insert(0, str(ATLAS_REPO / "src"))
    from atlas import mapping as M
    mane = M.load_mane_records(genes)
    maps = {}
    for g in genes:
        tmap, _ = M.build_transcript_map(GENE_TRANSCRIPT[g], mane[g])
        maps[g] = tmap
    return maps


def _boundaries(tmap):
    """(acceptors, donors, exon_lengths) in genomic coordinates.

    The first exon has no acceptor and the last no donor; both are excluded, since
    neither can be lost to produce a skip.
    """
    acc, don, length = [], [], {}
    n = len(tmap.exons)
    for i, e in enumerate(tmap.exons):
        ln = e.t_end - e.t_start
        if tmap.strand_sign == 1:
            a, d = e.g_start, e.g_end
        else:
            a, d = e.g_end, e.g_start
        length[a] = length[d] = ln
        if i > 0:
            acc.append((a, ln))
        if i < n - 1:
            don.append((d, ln))
    return acc, don


def classify(row, tmap) -> dict:
    ds = {k: row[f"spliceai_{k}"] for k in EVENT_NAMES}
    if not np.isfinite(list(ds.values())).any():
        return {"event": "none", "frame": "undetermined",
                "reason": "no SpliceAI record at this variant"}
    dom = max(ds, key=lambda k: (ds[k] if np.isfinite(ds[k]) else -np.inf))
    dp = row[f"spliceai_{dom.replace('ds_', 'dp_')}"]
    if not np.isfinite(dp):
        return {"event": EVENT_NAMES[dom], "dominant_score": ds[dom],
                "frame": "undetermined", "reason": "no position for the dominant event"}
    g_event = int(row["pos"]) + int(dp)
    acc, don = _boundaries(tmap)
    kind = EVENT_NAMES[dom]
    base = {"event": kind, "dominant_score": float(ds[dom]),
            "event_position": g_event, "dp": int(dp)}

    if kind.endswith("_loss"):
        sites = acc if kind.startswith("acceptor") else don
        hit = [(g, ln) for g, ln in sites if abs(g - g_event) <= BOUNDARY_TOL]
        if not hit:
            return base | {"frame": "undetermined",
                           "reason": "predicted loss does not sit on an internal "
                                     "exon boundary of this transcript"}
        g, ln = hit[0]
        return base | {"exon_length": ln, "shift": None,
                       "frame": "in_frame" if ln % 3 == 0 else "out_of_frame",
                       "reason": f"skipping an exon of {ln} bp"}

    sites = acc if kind.startswith("acceptor") else don
    if not sites:
        return base | {"frame": "undetermined", "reason": "no natural site of this kind"}
    g, ln = min(sites, key=lambda t: abs(t[0] - g_event))
    shift = abs(g_event - g)
    if shift == 0:
        return base | {"frame": "undetermined",
                       "reason": "gain predicted at the natural site itself"}
    return base | {"exon_length": ln, "shift": int(shift),
                   "frame": "in_frame" if shift % 3 == 0 else "out_of_frame",
                   "reason": f"cryptic site {shift} bp from the nearest natural "
                             f"{'acceptor' if kind.startswith('acceptor') else 'donor'}"}


def attribute() -> None:
    if not EVENTS.exists():
        raise SystemExit(f"[E6] {EVENTS} not found -- run the event scorer first "
                         "(see --subset output)")
    ev = pd.read_parquet(EVENTS)
    genes = sorted(set(ev["gene"]) & set(GENE_TRANSCRIPT))
    maps = _transcript_maps(genes)
    sub = pd.read_parquet(SUBSET)[["variant_id", "stratum"]]
    ev = ev.merge(sub, on="variant_id", how="left")

    rows = []
    for r in ev.to_dict("records"):
        if r["gene"] not in maps:
            rows.append(r | {"frame": "undetermined",
                             "reason": "no transcript model for this gene"})
            continue
        rows.append(r | classify(r, maps[r["gene"]]))
    out = pd.DataFrame(rows)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    keep = ["variant_id", "gene", "stratum", "event", "dominant_score",
            "event_position", "dp", "exon_length", "shift", "frame", "reason"]
    out[[c for c in keep if c in out.columns]].to_csv(OUT, index=False)

    summary = (out.groupby(["stratum", "frame"], observed=True).size()
                  .unstack(fill_value=0))
    for c in ("in_frame", "out_of_frame", "undetermined"):
        if c not in summary:
            summary[c] = 0
    summary["n"] = summary.sum(axis=1)
    summary["in_frame_share"] = (summary["in_frame"] / summary["n"]).round(4)
    summary["resolved_in_frame_share"] = (
        summary["in_frame"] / (summary["in_frame"] + summary["out_of_frame"])).round(4)
    summary.to_csv(REPORT_DIR / "inframe_attribution_summary.csv")

    print(f"[E6] wrote {OUT} ({len(out)} variants) and the per-stratum summary\n")
    print(summary.to_string())
    print("\n--- dominant event type ---")
    print(out.groupby(["stratum", "event"], observed=True).size()
             .unstack(fill_value=0).to_string())
    und = out[out.frame == "undetermined"]
    if len(und):
        print(f"\n--- why {len(und)} are undetermined ---")
        print(und["reason"].value_counts().to_string())
    print("\n[E6] Pangolin is not covered; see the module docstring for why.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--subset", action="store_true")
    ap.add_argument("--attribute", action="store_true")
    args = ap.parse_args()
    if args.subset:
        write_subset()
    elif args.attribute:
        attribute()
    else:
        ap.error("pass --subset or --attribute")


if __name__ == "__main__":
    main()
