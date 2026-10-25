"""E19.10 -- E3's in-sample thresholds refitted on twelve genes.

The genes: the seven training genes and CTCF, SFPQ, RAD51D, XRCC2 and SBDS. DDX3X and
TP53 are left out: their cell lines and label mechanisms differ (修改单 §2).

The data: the seven genes' rows of the analysis set, with its labels; the five genes'
rows of E18's newgenes_scored.parquet (all at 3-50 nt), read as E7 reads them
(evid_external._load_external), each under its primary label (the investigator-provided
ranges of its score set for four, SBDS's class column for the fifth). The rows are
taken in that order, the analysis set's own and then the five genes' in
evid_external.NEW_GENES order and each file's order, because the cluster bootstrap
draws genes in their order of first appearance.

The columns: twelve, the panel without NT. NT has no column of the analysis set's
variable on the new genes (NT on the original GPU did not reproduce its anchors; the
author's decision of 2026-10-23), so it gets no twelve-gene fit. The fusion is not
refitted either: NT is among its inputs, and it is a trained model, outside this stage.
Their rows carry the reason.

The fit: E3's own in-sample function (evid_interval_lr.thresholds_for, through
evid_suppl_c3.fit_cells), its gene-clustered bootstrap, its number of draws and its seed
key (tool, stratum, "full"), the frozen fit's, on the bands 3-10, 11-50 and 3-50 nt. The
same code on the seven genes alone must give evidence_thresholds.csv back exactly, or
the stage stops. No leave-one-gene-out fit (E19.10). No frozen file changes.

Output, in reports/evidence/suppl_analyses/:
  twelve_gene_refit.csv   per column x band x side x tier: the frozen seven-gene fit and
                          the twelve-gene fit (genes, damaging and normal counts,
                          threshold, whether the tier is reached), the absolute and
                          relative change of the threshold (twelve minus seven; relative
                          to the seven-gene threshold's magnitude), and tier_change
                          ('gained', 'lost', 'unchanged')

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_twelve_gene
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import evid_common as K
from . import evid_external as X
from . import evid_interval_lr as E3
from . import evid_suppl_c3 as C3

OUT_DIR = Path("reports/evidence/suppl_analyses")
OUT = OUT_DIR / "twelve_gene_refit.csv"
NOT_REFITTED = {
    "nt": ("not refitted: NT has no column of the analysis set's variable on the new genes "
           "(on the original GPU it did not reproduce its anchors; the author's decision of "
           "2026-10-23)"),
    K.FUSION: ("not refitted: the fusion is a trained model, outside this stage, and NT is "
               "among its inputs")}


def twelve_gene_frame(df: pd.DataFrame, tools: list[str]) -> pd.DataFrame:
    cols = ["gene", "stratum", "y_assay"] + tools
    parts = [df[cols]]
    for g in X.NEW_GENES:
        e, labels = X._load_external(g)
        lab = X.NEW_GENE_LABEL[g]
        if labels != [lab] or not set(e["stratum"]) <= {"s3_10", "s11_50"}:
            raise SystemExit(f"[12] {g}: not the primary label at 3-50 nt that E7 reads")
        missing = [t for t in tools if t not in e.columns]
        if missing:
            raise SystemExit(f"[12] {g}: no column {missing}")
        parts.append(e.assign(y_assay=e[lab].astype("Float64"))[cols])
    out = pd.concat(parts, ignore_index=True)
    out["y_assay"] = out["y_assay"].astype("Float64")
    return out


def main() -> None:
    cfg = yaml.safe_load(E3.CONFIG_PATH.read_text())
    path_cuts, ben_cuts = K.acmg_bands(cfg)
    cuts = {"path": path_cuts, "ben": ben_cuts}
    df = K.load_set()
    tools = [t for t in K.panel_of(df) if t not in NOT_REFITTED]
    if len(tools) != 12:
        raise SystemExit(f"[12] {len(tools)} columns, not twelve: {tools}")
    frozen = C3.frozen_cells()
    seven = C3.fit_cells(df, df["y_assay"].astype("Float64"), tools, cuts)
    C3.check_reproduces_frozen(seven, frozen[frozen["tool"].isin(tools)])
    print(f"[12] the same fit on the seven genes gives {C3.FROZEN} back exactly")

    d12 = twelve_gene_frame(df, tools)
    fit = C3.fit_cells(d12, d12["y_assay"], tools, cuts)
    keys = ["tool", "stratum", "side", "tier"]
    n_genes = {(t, st): int(d12.loc[K.stratum_frame(d12, st).index]
                            .dropna(subset=["y_assay", t])["gene"].nunique())
               for t in tools for st in C3.BANDS}
    n7 = {(t, st): int(df.loc[K.stratum_frame(df, st).index]
                       .dropna(subset=["y_assay", t])["gene"].nunique())
          for t in tools for st in C3.BANDS}
    m = frozen[frozen["tool"].isin(tools)].rename(columns={
        "n_pos": "n_pos_7", "n_neg": "n_neg_7", "threshold": "threshold_7",
        "reached": "reached_7", "status": "status_7"}).merge(
        fit.drop(columns=["cut_lr"]).rename(columns={
            "n_pos": "n_pos_12", "n_neg": "n_neg_12", "threshold": "threshold_12",
            "reached": "reached_12", "status": "status_12"}),
        on=keys, validate="one_to_one")
    m["n_genes_7"] = [n7[(t, s)] for t, s in zip(m.tool, m.stratum)]
    m["n_genes_12"] = [n_genes[(t, s)] for t, s in zip(m.tool, m.stratum)]
    m["abs_change"] = m["threshold_12"] - m["threshold_7"]
    m["rel_change"] = m["abs_change"] / m["threshold_7"].abs()
    m["tier_change"] = np.select(
        [m.reached_12 & ~m.reached_7, ~m.reached_12 & m.reached_7], ["gained", "lost"],
        "unchanged")
    rest = frozen[frozen["tool"].isin(list(NOT_REFITTED))].rename(columns={
        "n_pos": "n_pos_7", "n_neg": "n_neg_7", "threshold": "threshold_7",
        "reached": "reached_7", "status": "status_7"})
    rest = rest.assign(status_12="not refitted", reason=rest["tool"].map(NOT_REFITTED))
    out = pd.concat([m, rest], ignore_index=True)[
        ["tool", "stratum", "side", "tier", "cut_lr",
         "n_genes_7", "n_pos_7", "n_neg_7", "threshold_7", "reached_7", "status_7",
         "n_genes_12", "n_pos_12", "n_neg_12", "threshold_12", "reached_12", "status_12",
         "abs_change", "rel_change", "tier_change", "reason"]]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    ch = out[out.tier_change.isin(["gained", "lost"])]
    print(f"[12] {len(ch)} of {len(m)} cells change tier"
          + "".join(f"\n      {r.tool} {r.stratum} {r.side} {r.tier} {r.tier_change}"
                    for r in ch.itertuples()))


if __name__ == "__main__":
    main()
