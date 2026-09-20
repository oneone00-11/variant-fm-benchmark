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

Pangolin is covered too, on the seven frozen genes. The atlas persisted only its
maximum score, so its gain/loss values and their positions are re-scored for the
same kind of subset (`src/evid_score_pangolin_events.py`). Its output is coarser
than SpliceAI's: one gain and one loss per transcript record rather than four typed
events, so a predicted loss cannot be read as donor-side or acceptor-side from the
score alone and is matched against both boundary sets. On the external gene Pangolin
is not covered at all -- its gffutils database covers the frozen genes' chromosomes
only.

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
    "EVID_ATLAS_REPO", "/Users/cliffzhang/work/functional-standard-atlas"))
EVID_DIR = Path("data/evidence")
REPORT_DIR = Path("reports/evidence")
SUBSET = EVID_DIR / "inframe_subset.parquet"
EVENTS = EVID_DIR / "inframe_events.parquet"
PANG_SUBSET = EVID_DIR / "inframe_subset_pangolin.parquet"
PANG_EVENTS = EVID_DIR / "inframe_events_pangolin.parquet"
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
    # external genes, on their own MANE Select transcripts
    "TP53": "NM_000546.6", "DDX3X": "ENST00000644876.2",
}

# E2.7: the same attribution on the external genes. Each needs its own false-
# positive subset, because the cut point is applied to that gene's own scores.
EXTERNAL = {
    "ddx3x": {
        "scored": EVID_DIR / "external/ddx3x_scored.parquet",
        "labels": EVID_DIR / "external/ddx3x_splice.parquet",
        "label_cols": ["y_control_anchored", "y_fdr"],
        "spliceai_col": "spliceai_walker",
        "fasta": EVID_DIR / "refs/chrX.fa",
        "events": EVID_DIR / "external/ddx3x_inframe_events.parquet",
        "subset": EVID_DIR / "external/ddx3x_inframe_subset.parquet",
        "pangolin_events": EVID_DIR / "external/ddx3x_pangolin.parquet",
    },
    "tp53": {
        # TP53's columns come from the published study's frozen matrix, so its
        # SpliceAI value is the distance-50 one unless the Walker-basis column has
        # been merged; whichever is present is used and recorded in the subset.
        "scored": Path("data/external/tp53_splice_scored_v2.parquet"),
        "labels": Path("data/external/tp53_splice_scored_v2.parquet"),
        "label_cols": ["y_tp53_control_anchored"],
        "spliceai_col": "spliceai_walker",
        "fasta": None,          # chr17 is in the atlas subset fasta
        "events": EVID_DIR / "external/tp53_inframe_events.parquet",
        "subset": EVID_DIR / "external/tp53_inframe_subset.parquet",
        "pangolin_events": EVID_DIR / "external/tp53_pangolin.parquet",
    },
}


def write_subset() -> None:
    """The variants E6 is about: called by the tool, called normal by the assay.

    The cut point is Walker's, which is calibrated for SpliceAI only. It is applied
    to Pangolin as well so the two tools' false-positive sets are defined the same
    way; that is a comparison device, not a claim that 0.2 is Pangolin's threshold.
    """
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    pp3 = cfg["thresholds"]["pp3"]["value"]
    df = K.load_set()
    for tool, out_path in (("spliceai", SUBSET), ("pangolin", PANG_SUBSET)):
        col = "spliceai_walker" if (tool == "spliceai"
                                    and "spliceai_walker" in df.columns) else tool
        sub = df[(df[col] >= pp3) & (df["y_assay"] == 0)].copy()
        sub["score_column"] = col
        sub[["variant_id", "gene", "chrom", "pos", "ref", "alt", "stratum",
             "hgvs_c", col, "score_column"]].to_parquet(out_path, index=False)
        print(f"[E6] {tool}: {len(sub):,} variants scored >= {pp3} on `{col}` "
              "and labelled normal")
        print(sub.groupby("stratum", observed=True).size().to_string())
        print(f"     wrote {out_path}")
    print("\n[E6] next, in the model environments:")
    print("  <atlas>/models/spliceai/.venv/bin/python phase1/src/evid_score_spliceai_events.py \\")
    print(f"      --variants {SUBSET} --distance 4999 --out {EVENTS} --run")
    print("  <atlas>/models/pangolin/.venv/bin/python phase1/src/evid_score_pangolin_events.py \\")
    print(f"      --variants {PANG_SUBSET} --out {PANG_EVENTS} --run")


# ---------------------------------------------------------------------------
# frame inference
# ---------------------------------------------------------------------------
def _transcript_maps(genes: list[str]) -> dict:
    sys.path.insert(0, str(ATLAS_REPO / "src"))
    from atlas import mapping as M
    # atlas.mapping resolves its reference cache RELATIVE to the working directory
    # ("data/raw/reference"). These stages run from phase1/, so every cached
    # Mutalyzer model and Ensembl lookup would be missed and re-fetched -- slowly,
    # and with a second copy of the cache appearing under phase1/. Point it at the
    # atlas's own cache.
    M.REF_CACHE = ATLAS_REPO / "data/raw/reference"
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


def classify_pangolin(row, tmap) -> dict:
    gain, loss = row["pangolin_gain"], row["pangolin_loss"]
    if not (np.isfinite(gain) or np.isfinite(loss)):
        return {"event": "none", "frame": "undetermined",
                "reason": "no Pangolin record at this variant"}
    g = gain if np.isfinite(gain) else -np.inf
    l = -loss if np.isfinite(loss) else -np.inf
    is_gain = g >= l
    dp = row["pangolin_gain_pos"] if is_gain else row["pangolin_loss_pos"]
    if not np.isfinite(dp):
        return {"event": "gain" if is_gain else "loss", "frame": "undetermined",
                "reason": "no position for the dominant event"}
    g_event = int(row["pos"]) + int(dp)
    acc, don = _boundaries(tmap)
    base = {"event": "gain" if is_gain else "loss",
            "dominant_score": float(g if is_gain else l),
            "event_position": g_event, "dp": int(dp)}
    sites = acc + don          # Pangolin does not say which kind it is
    if not sites:
        return base | {"frame": "undetermined", "reason": "no internal boundaries"}
    if not is_gain:
        hit = [(gg, ln) for gg, ln in sites if abs(gg - g_event) <= BOUNDARY_TOL]
        if not hit:
            return base | {"frame": "undetermined",
                           "reason": "predicted loss does not sit on an internal "
                                     "exon boundary of this transcript"}
        gg, ln = hit[0]
        return base | {"exon_length": ln,
                       "frame": "in_frame" if ln % 3 == 0 else "out_of_frame",
                       "reason": f"skipping an exon of {ln} bp"}
    # Pangolin reports one gain without saying whether the gained site is a donor
    # or an acceptor. The reading-frame shift is only meaningful against a natural
    # site of the SAME kind, so measuring against whichever boundary happens to be
    # nearest turns a third of these into a coin flip. Both readings are taken; the
    # call stands only where they agree.
    reads = {}
    for kind, group in (("acceptor", acc), ("donor", don)):
        if not group:
            continue
        gg, ln = min(group, key=lambda tt: abs(tt[0] - g_event))
        sh = abs(g_event - gg)
        reads[kind] = (sh, ln)
    if not reads:
        return base | {"frame": "undetermined", "reason": "no natural site of either kind"}
    if any(sh == 0 for sh, _ in reads.values()):
        return base | {"frame": "undetermined",
                       "reason": "gain predicted at a natural site itself"}
    frames = {k: ("in_frame" if sh % 3 == 0 else "out_of_frame")
              for k, (sh, _) in reads.items()}
    detail = ", ".join(f"{k} {reads[k][0]} bp -> {frames[k]}" for k in sorted(reads))
    if len(set(frames.values())) > 1:
        return base | {"frame": "undetermined",
                       "reason": f"Pangolin does not type the gained site and the "
                                 f"two readings disagree ({detail})"}
    kind = sorted(reads)[0]
    sh, ln = reads[kind]
    return base | {"exon_length": ln, "shift": int(sh),
                   "frame": frames[kind],
                   "reason": f"cryptic site, both readings agree ({detail})"}


def _attribute_one(events_path: Path, subset_path: Path, tool: str,
                   classifier) -> pd.DataFrame | None:
    if not events_path.exists():
        print(f"[E6] {tool}: {events_path} not found -- skipped")
        return None
    ev = pd.read_parquet(events_path)
    genes = sorted(set(ev["gene"]) & set(GENE_TRANSCRIPT))
    maps = _transcript_maps(genes)
    sub = pd.read_parquet(subset_path)[["variant_id", "stratum"]]
    ev = ev.merge(sub, on="variant_id", how="left")
    rows = []
    for r in ev.to_dict("records"):
        if r["gene"] not in maps:
            rows.append(r | {"frame": "undetermined",
                             "reason": "no transcript model for this gene"})
            continue
        rows.append(r | classifier(r, maps[r["gene"]]))
    out = pd.DataFrame(rows)
    out["tool"] = tool
    return out


def attribute() -> None:
    parts = [p for p in (
        _attribute_one(EVENTS, SUBSET, "spliceai", classify),
        _attribute_one(PANG_EVENTS, PANG_SUBSET, "pangolin", classify_pangolin),
    ) if p is not None]
    if not parts:
        raise SystemExit("[E6] no event table found -- run the event scorers first "
                         "(see --subset output)")
    out = pd.concat(parts, ignore_index=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    keep = ["tool", "variant_id", "gene", "stratum", "event", "dominant_score",
            "event_position", "dp", "exon_length", "shift", "frame", "reason"]
    out[[c for c in keep if c in out.columns]].to_csv(OUT, index=False)

    summary = (out.groupby(["tool", "stratum", "frame"], observed=True).size()
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
    print(out.groupby(["tool", "stratum", "event"], observed=True).size()
             .unstack(fill_value=0).to_string())
    und = out[out.frame == "undetermined"]
    if len(und):
        print(f"\n--- why {len(und)} are undetermined ---")
        print(und["reason"].value_counts().to_string())
    tools = sorted(out["tool"].unique())
    print(f"\n[E6] tools covered: {', '.join(tools)}")


def external_subset(gene: str) -> None:
    """The external gene's false positives: its own score over the cut point, its
    own assay calling the variant normal. A variant is taken when EITHER label
    definition calls it normal and the tool calls it, and the definition is carried
    so the attribution can be read either way."""
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    pp3 = cfg["thresholds"]["pp3"]["value"]
    spec = EXTERNAL[gene]
    lab = pd.read_parquet(spec["labels"])
    sc = pd.read_parquet(spec["scored"])
    df = lab if spec["labels"] == spec["scored"] else lab.merge(
        sc, on="variant_id", how="left")
    if gene == "tp53":
        df = df[df["is_splice"]].copy()
        df["intron_offset_abs"] = df["intron_offset"].abs()
        df["stratum"] = df["intron_offset_abs"].map(
            lambda o: "pm12" if o <= 2 else ("s3_10" if o <= 10 else "s11_50"))
        # the published study's control-anchored definition: RFS is anchored at -1
        # for synonymous and +1 for nonsense, so damaging is RFS > 0
        df["y_tp53_control_anchored"] = (
            df["func_pathogenicity"].to_numpy(dtype=float) > 0).astype(float)
        w = EVID_DIR / "external/tp53_spliceai_walker.parquet"
        if w.exists():
            ww = pd.read_parquet(w)
            df = df.merge(ww[["variant_id", "spliceai_walker"]],
                          on="variant_id", how="left")
        else:
            df["spliceai_walker"] = df["spliceai"]
    col = spec["spliceai_col"]
    normal = np.zeros(len(df), dtype=bool)
    for c in spec["label_cols"]:
        normal |= (df[c] == 0).to_numpy()
    sub = df[(df[col] >= pp3) & normal].copy()
    sub["score_column"] = col
    keep = ["variant_id", "gene", "chrom", "pos", "ref", "alt", "stratum", "hgvs_c",
            col, "score_column"] + spec["label_cols"]
    sub[[k for k in keep if k in sub.columns]].to_parquet(spec["subset"], index=False)
    print(f"[E2.7] {gene.upper()}: {len(sub):,} variants over {pp3} on `{col}` and "
          f"called normal by at least one label definition")
    print(sub.groupby("stratum", observed=True).size().to_string())
    print(f"     wrote {spec['subset']}")
    print("\n[E2.7] next, in the SpliceAI environment:")
    fasta = f" \\\n      --fasta {spec['fasta']}" if spec["fasta"] else ""
    print("  <atlas>/models/spliceai/.venv/bin/python phase1/src/evid_score_spliceai_events.py \\")
    print(f"      --variants {spec['subset']} --distance 4999{fasta} \\")
    print(f"      --out {spec['events']} --run")


def attribute_external(gene: str) -> None:
    spec = EXTERNAL[gene]
    parts = []
    sai = _attribute_one(spec["events"], spec["subset"], "spliceai", classify)
    if sai is not None:
        parts.append(sai)
    pg = spec["pangolin_events"]
    if pg.exists():
        sub_ids = set(pd.read_parquet(spec["subset"])["variant_id"])
        ev = pd.read_parquet(pg)
        ev = ev[ev["variant_id"].isin(sub_ids)]
        if len(ev):
            tmp = EVID_DIR / f"external/_{gene}_pangolin_subset.parquet"
            ev.to_parquet(tmp, index=False)
            p = _attribute_one(tmp, spec["subset"], "pangolin", classify_pangolin)
            tmp.unlink(missing_ok=True)
            if p is not None:
                parts.append(p)
    if not parts:
        raise SystemExit(f"[E2.7] {gene}: no event table -- see --external-subset")
    out = pd.concat(parts, ignore_index=True)
    out["gene_set"] = gene.upper()
    path = REPORT_DIR / f"inframe_attribution_{gene}.csv"
    keep = ["tool", "gene_set", "variant_id", "gene", "stratum", "event",
            "dominant_score", "event_position", "dp", "exon_length", "shift",
            "frame", "reason"]
    out[[c for c in keep if c in out.columns]].to_csv(path, index=False)
    summary = (out.groupby(["tool", "stratum", "frame"], observed=True).size()
                  .unstack(fill_value=0))
    for c in ("in_frame", "out_of_frame", "undetermined"):
        if c not in summary:
            summary[c] = 0
    summary["n"] = summary.sum(axis=1)
    summary["in_frame_share"] = (summary["in_frame"] / summary["n"]).round(4)
    summary["resolved_in_frame_share"] = (
        summary["in_frame"] / (summary["in_frame"] + summary["out_of_frame"])).round(4)
    summary.to_csv(REPORT_DIR / f"inframe_attribution_{gene}_summary.csv")
    print(f"[E2.7] wrote {path} ({len(out)} variants)\n")
    print(summary.to_string())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--subset", action="store_true")
    ap.add_argument("--attribute", action="store_true")
    ap.add_argument("--external-subset", choices=sorted(EXTERNAL))
    ap.add_argument("--external-attribute", choices=sorted(EXTERNAL))
    args = ap.parse_args()
    if args.subset:
        write_subset()
    elif args.attribute:
        attribute()
    elif args.external_subset:
        external_subset(args.external_subset)
    elif args.external_attribute:
        attribute_external(args.external_attribute)
    else:
        ap.error("pass --subset, --attribute, --external-subset or "
                 "--external-attribute")


if __name__ == "__main__":
    main()
