"""E20.1 (C12): the true-positive subsets, and the score inputs that carry anchor 2's
variants, written beside E6's and E2.7's false-positive files.

The definition is E6's and E2.7's (evid_inframe.write_subset and external_subset) with
"called normal by the assay" replaced by "called damaging", on the bands PP3 and BP4
apply to (3-10 and 11-50 nt; DDX3X's variants reach 25 nt and TP53's 12 nt):

  * the seven training genes: the analysis set's spliceai_walker, and separately its
    pangolin, at or above Walker's PP3 cut point (config/walker2023.yaml, 0.2), with
    y_assay == 1 (as in E6, the cut point is applied to Pangolin only so that the two
    tools' subsets are defined the same way);
  * DDX3X and TP53, SpliceAI only: each gene's scored table, the same column and cut
    point, damaging under any of its label definitions (the union; DDX3X four, TP53 its
    control-anchored one, derived exactly as external_subset derives it).

The same code with "called normal" (and every band) must give E2.7's stored DDX3X and
TP53 subsets back exactly, and must contain E6's stored subsets of the seven genes; the
variants it finds beyond E6's are recorded (E6 drew its subsets on 2026-09-17, before
E16d switched PALB2 and BARD1 to the archived classification).

Anchor 2 (E20.1): 10 variants from each of E6's event files (SpliceAI, Pangolin), sorted
by variant_id and taken at equal steps, positions i * n // 10 for i = 0 ... 9. They are
scored in the same run as the seven genes' true positives: the score inputs are the
subset's rows and the anchors', each with its role.

Writes, under phase1/data/evidence:
  inframe_subset_tp.parquet, inframe_subset_pangolin_tp.parquet    the seven genes
  inframe_score_input_tp.parquet, inframe_score_input_pangolin_tp.parquet
  external/ddx3x_inframe_subset_tp.parquet, external/tp53_inframe_subset_12nt_tp.parquet

Run in the reference image (修改单 §0 rule 9), offline:
    docker run ... IMAGE bash scripts/stage_E20/in_image.sh scripts/stage_E20/write_c12_subsets.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

PHASE1 = REPO / "phase1"
sys.path.insert(0, str(PHASE1))
os.chdir(PHASE1)
from src import evid_common as K  # noqa: E402
from src import evid_inframe as F  # noqa: E402

BANDS = ["s3_10", "s11_50"]
EV = Path("data/evidence")
OUT = {
    "spliceai": EV / "inframe_subset_tp.parquet",
    "pangolin": EV / "inframe_subset_pangolin_tp.parquet",
    "spliceai_input": EV / "inframe_score_input_tp.parquet",
    "pangolin_input": EV / "inframe_score_input_pangolin_tp.parquet",
    "ddx3x": EV / "external/ddx3x_inframe_subset_tp.parquet",
    "tp53": EV / "external/tp53_inframe_subset_12nt_tp.parquet",
}
E6_EVENTS = {"spliceai": F.EVENTS, "pangolin": F.PANG_EVENTS}
E6_SUBSETS = {"spliceai": F.SUBSET, "pangolin": F.PANG_SUBSET}
SCORER_COLS = ["variant_id", "gene", "chrom", "pos", "ref", "alt"]
N_ANCHORS = 10


def pp3_cut() -> float:
    return yaml.safe_load(F.CONFIG_PATH.read_text())["thresholds"]["pp3"]["value"]


def seven_genes(df: pd.DataFrame, col: str, y: int, bands: list[str] | None) -> pd.DataFrame:
    """write_subset's selection, with the assay's call y and, if given, the bands."""
    keep = (df[col] >= pp3_cut()) & (df["y_assay"] == y)
    if bands is not None:
        keep &= df["stratum"].isin(bands)
    sub = df[keep].copy()
    sub["score_column"] = col
    return sub[["variant_id", "gene", "chrom", "pos", "ref", "alt", "stratum",
                "hgvs_c", col, "score_column"]]


def external(gene: str, y: int, bands: list[str] | None) -> pd.DataFrame:
    """external_subset's selection, with "called y by at least one label definition"."""
    spec = F.EXTERNAL[gene]
    lab = pd.read_parquet(spec["labels"])
    sc = pd.read_parquet(spec["scored"])
    df = lab if spec["labels"] == spec["scored"] else lab.merge(sc, on="variant_id", how="left")
    if gene == "tp53":
        df = df[df["is_splice"]].copy()
        df["intron_offset_abs"] = df["intron_offset"].abs()
        df["stratum"] = df["intron_offset_abs"].map(
            lambda o: "pm12" if o <= 2 else ("s3_10" if o <= 10 else "s11_50"))
        df["y_tp53_control_anchored"] = (
            df["func_pathogenicity"].to_numpy(dtype=float) > 0).astype(float)
        if "spliceai_walker" not in df.columns:
            raise SystemExit("TP53's scored table has no spliceai_walker column: stopped")
    col = spec["spliceai_col"]
    called = np.zeros(len(df), dtype=bool)
    for c in spec["label_cols"]:
        called |= (df[c] == y).to_numpy()
    keep = (df[col] >= pp3_cut()).to_numpy() & called
    if bands is not None:
        keep &= df["stratum"].isin(bands).to_numpy()
    sub = df[keep].copy()
    sub["score_column"] = col
    cols = ["variant_id", "gene", "chrom", "pos", "ref", "alt", "stratum", "hgvs_c",
            col, "score_column"] + spec["label_cols"]
    return sub[[c for c in cols if c in sub.columns]]


def anchors(tool: str) -> pd.DataFrame:
    ev = pd.read_parquet(E6_EVENTS[tool]).sort_values("variant_id", kind="mergesort")
    ev = ev.reset_index(drop=True)
    n = len(ev)
    picks = [i * n // N_ANCHORS for i in range(N_ANCHORS)]
    out = ev.loc[picks, SCORER_COLS].copy()
    out["role"] = "anchor"
    return out.reset_index(drop=True)


def same_frame(a: pd.DataFrame, b: pd.DataFrame) -> bool:
    a = a.reset_index(drop=True)
    b = b.reset_index(drop=True)
    if list(a.columns) != list(b.columns) or len(a) != len(b):
        return False
    for c in a.columns:
        x, y = a[c].to_numpy(dtype=object), b[c].to_numpy(dtype=object)
        both_na = pd.isna(x) & pd.isna(y)
        if not (both_na | (x == y)).all():
            return False
    return True


def main() -> None:
    with RunLog(__file__) as log:
        log.input(F.CONFIG_PATH)
        log.input(K.SET_PATH)
        df = K.load_set()
        counts = {}
        for tool, col in (("spliceai", "spliceai_walker"), ("pangolin", "pangolin")):
            # the definition, checked against E6's stored false positives first
            fp_now = seven_genes(df, col, 0, None)
            log.input(E6_SUBSETS[tool])
            stored = pd.read_parquet(E6_SUBSETS[tool])
            extra = fp_now[~fp_now["variant_id"].isin(stored["variant_id"])]
            missing = stored[~stored["variant_id"].isin(fp_now["variant_id"])]
            if len(missing) or set(extra["gene"]) - {"PALB2"}:
                raise SystemExit(f"{tool}: the false-positive definition does not contain E6's "
                                 f"stored subset, or differs from it beyond PALB2: stopped")
            log.note(f"{tool}_false_positive_definition_check", {
                "stored": len(stored), "now": len(fp_now),
                "stored_contained": True,
                "beyond_stored": extra[["variant_id", "gene", "stratum", "hgvs_c"]]
                .to_dict("records")})
            tp = seven_genes(df, col, 1, BANDS)
            tp.to_parquet(OUT[tool], index=False)
            log.output(OUT[tool])
            log.input(E6_EVENTS[tool])
            anc = anchors(tool)
            if anc["variant_id"].isin(tp["variant_id"]).any():
                raise SystemExit(f"{tool}: an anchor is also a true positive: stopped")
            inp = pd.concat([tp[SCORER_COLS].assign(role="true_positive"), anc],
                            ignore_index=True)
            inp.to_parquet(OUT[f"{tool}_input"], index=False)
            log.output(OUT[f"{tool}_input"])
            counts[f"seven_genes_{tool}"] = tp.groupby("stratum").size().to_dict()
            log.note(f"{tool}_anchors", {"from": str(E6_EVENTS[tool]), "rule":
                     "sorted by variant_id, positions i * n // 10, i = 0 ... 9",
                     "n_events": int(len(pd.read_parquet(E6_EVENTS[tool]))),
                     "variant_ids": anc["variant_id"].tolist()})
        for gene in ("ddx3x", "tp53"):
            spec = F.EXTERNAL[gene]
            for p in dict.fromkeys([spec["labels"], spec["scored"], spec["subset"]]):
                log.input(p)
            fp_now = external(gene, 0, None)
            stored = pd.read_parquet(spec["subset"])
            if not same_frame(fp_now, stored):
                raise SystemExit(f"{gene}: the false-positive definition does not give E2.7's "
                                 "stored subset back: stopped")
            log.note(f"{gene}_false_positive_definition_check",
                     {"stored": len(stored), "now": len(fp_now), "identical": True})
            tp = external(gene, 1, BANDS)
            tp.to_parquet(OUT[gene], index=False)
            log.output(OUT[gene])
            counts[gene] = tp.groupby("stratum").size().to_dict()
            counts[f"{gene}_by_label"] = {c: int((tp[c] == 1).sum()) for c in spec["label_cols"]}
        log.note("true_positive_counts", counts)
        print(counts)


if __name__ == "__main__":
    main()
