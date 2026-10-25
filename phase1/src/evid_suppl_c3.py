"""E19 C3 -- label sensitivity: the in-sample thresholds refitted under other labels.

Every refit is E3's in-sample fit (evid_interval_lr.thresholds_for: the gene-clustered
bootstrap, N_BOOT_GENE draws, seeded from the cell key (tool, stratum, "full") as E3
seeds it), on the bands 3-10, 11-50 and 3-50 nt, for every column E3 fits: the panel
and the fusion. Only the labels change. The scores, the tier cuts and the seed keys do
not; nor does the fusion's out-of-fold column, which is trained on the within-gene
rank of func_pathogenicity, not on these labels (evid_common.logo_fusion), and is
computed here exactly as E3 computes it. Each refit is set against the frozen
in-sample thresholds of evidence_thresholds.csv, cell by cell. Before any refit, the
same fit on the analysis set's own labels must give every frozen threshold back
exactly; otherwise the stage stops.

  1. intermediate classes (E19.5 item 1). The variants at 3-50 nt that the label
     files class as intermediate, indeterminate or enriched, and so leave without a
     label (y_assay_source 'indeterminate', 337), are (a) all taken as damaging and
     (b) all taken as normal.
  2. superseded labels (item 2). PALB2's and BARD1's labels are taken from the files
     in data/assay_labels/superseded/, the ones the label source audit replaced on
     2026-10-21, joined exactly as evid_build_set.attach_labels joins the current
     files (same key, class column and class map); the other five genes keep theirs.
  4. without BARD1, without PALB2, without both (item 4). Nothing is refitted. E11's
     threshold_dilution.csv has fitted every subset of the genes with the same
     function, so the in-sample threshold, and whether the tier is reached, are read
     there for the training subset that leaves out BARD1, PALB2 or both, beside the
     fit on every gene. E11 fits four columns (spliceai_walker, pangolin, alphagenome,
     avi), on the PP3 side, at Moderate and Strong; those are the cells this gives.
     In the 11-50 nt band the genes are six (BRCA2 has no variant there).

Item 3, a second BRCA2 dataset, is a stage of its own (E19 C3, BRCA2).

Columns of the two refit files, one row per analysis x column x band x side x tier:
  n_added, n_removed, n_switched   the band's variants (whatever their scores) that
                                   gain a label, lose one, or change damaging <->
                                   normal, against the analysis set's labels
  frozen_*, the rest               the frozen fit and this refit: damaging and normal
                                   counts among the column's scored variants, the
                                   threshold (NaN: the tier is not reached), whether
                                   it is reached
  tier_change                      'unchanged'; 'gained' (reached now, not frozen);
                                   'lost' (reached frozen, not now)
  threshold_changed                the threshold differs from the frozen one (both
                                   NaN counts as the same)

Outputs, in reports/evidence/suppl_analyses/:
  c3_intermediate.csv          analyses 'intermediate as damaging', 'intermediate as
                               normal'
  c3_superseded_labels.csv     analysis 'superseded PALB2 and BARD1 labels'
  c3_without_bard1_palb2.csv   per column x band x tier x subset: the subset's genes,
                               its threshold and whether the tier is reached, beside
                               every gene's

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_c3
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import config as C
from . import evid_common as K
from . import evid_interval_lr as E3
from .evid_build_set import _norm_hgvs

OUT_DIR = Path("reports/evidence/suppl_analyses")
OUT_INTERMEDIATE = OUT_DIR / "c3_intermediate.csv"
OUT_SUPERSEDED = OUT_DIR / "c3_superseded_labels.csv"
OUT_WITHOUT = OUT_DIR / "c3_without_bard1_palb2.csv"
FROZEN = E3.REPORT_DIR / "evidence_thresholds.csv"
DILUTION = E3.REPORT_DIR / "threshold_dilution.csv"
DILUTION_SUMMARY = E3.REPORT_DIR / "threshold_dilution_summary.csv"
SUPERSEDED_DIR = C.ASSAY_LABEL_DIR / "superseded"
SUPERSEDED_GENES = ["PALB2", "BARD1"]
BANDS = ["s3_10", "s11_50", "s3_50"]
SIDES = [("pp3", "path"), ("bp4", "ben")]
INDETERMINATE = "indeterminate"            # evid_build_set.attach_labels' source name
SUBSETS = [("without BARD1", {"BARD1"}), ("without PALB2", {"PALB2"}),
           ("without BARD1 and PALB2", {"BARD1", "PALB2"})]
ALL_GENES = "every gene"


def _same(a: float, b: float) -> bool:
    return (np.isnan(a) and np.isnan(b)) or a == b


def label_changes(df: pd.DataFrame, y: pd.Series) -> dict[str, dict[str, int]]:
    """Per band: variants gaining a label, losing one, switching side."""
    out = {}
    for stratum in BANDS:
        idx = K.stratum_frame(df, stratum).index
        old, new = df.loc[idx, "y_assay"], y.loc[idx]
        out[stratum] = {
            "n_added": int((old.isna() & new.notna()).sum()),
            "n_removed": int((old.notna() & new.isna()).sum()),
            "n_switched": int((old.notna() & new.notna() & (old != new)).sum())}
    return out


def fit_cells(df: pd.DataFrame, y: pd.Series, tools: list[str], cuts: dict) -> pd.DataFrame:
    """E3's in-sample fit for every column x band, under the labels y (aligned to df).
    The rows are taken in the analysis set's order, as E3 takes them: the cluster
    bootstrap draws genes in their order of first appearance."""
    rows = []
    for tool in tools:
        for stratum in BANDS:
            sub = K.stratum_frame(df, stratum)
            keep = y.loc[sub.index].notna() & sub[tool].notna()
            yy = y.loc[sub.index][keep].to_numpy(dtype=float)
            s = sub.loc[keep, tool].to_numpy(dtype=float)
            genes = sub.loc[keep, "gene"].astype(str).to_numpy()
            base = {"tool": tool, "stratum": stratum,
                    "n_pos": int((yy == 1).sum()), "n_neg": int((yy == 0).sum())}
            full, reason = None, ""
            if base["n_pos"] < K.MIN_POS or base["n_neg"] < K.MIN_NEG:
                reason = f"fewer than {K.MIN_POS} labelled variants on one side"
            else:
                full = E3.thresholds_for(yy, s, genes, cuts["path"], cuts["ben"],
                                         n_boot=E3.N_BOOT_GENE, cluster=True,
                                         cell=(tool, stratum, "full"))
                if full is None:
                    reason = "score takes too few distinct values"
            for side, key in SIDES:
                for tier in E3.TIERS:
                    thr = full[key][tier] if full is not None else np.nan
                    rows.append(base | {"side": side, "tier": tier,
                                        "cut_lr": cuts[key][tier],
                                        "status": "ok" if full is not None else "not evaluable",
                                        "reason": reason, "threshold": thr,
                                        "reached": bool(np.isfinite(thr))})
    return pd.DataFrame(rows)


def frozen_cells() -> pd.DataFrame:
    """evidence_thresholds.csv's in-sample fit, in fit_cells' long form."""
    e = K.read_back(FROZEN)
    e = e[e["stratum"].isin(BANDS)]
    parts = []
    for side, _ in SIDES:
        parts.append(e.assign(side=side).rename(columns={
            f"{side}_threshold_insample": "threshold", f"{side}_threshold_reachable": "reached",
            "path_cut_lr" if side == "pp3" else "ben_cut_lr": "cut_lr"})[
            ["tool", "stratum", "side", "tier", "cut_lr", "n_pos", "n_neg", "status",
             "threshold", "reached"]])
    return pd.concat(parts, ignore_index=True)


def check_reproduces_frozen(fit: pd.DataFrame, frozen: pd.DataFrame) -> None:
    keys = ["tool", "stratum", "side", "tier"]
    m = fit.merge(frozen, on=keys, suffixes=("", "_frozen"), how="outer", indicator=True)
    bad = m[(m["_merge"] != "both")
            | ~m.apply(lambda r: _same(r["threshold"], r["threshold_frozen"]), axis=1)
            | (m["reached"] != m["reached_frozen"]) | (m["n_pos"] != m["n_pos_frozen"])
            | (m["n_neg"] != m["n_neg_frozen"]) | (m["status"] != m["status_frozen"])]
    if len(bad):
        raise SystemExit(f"[C3] the in-sample fit on the analysis set's own labels does not "
                         f"give {FROZEN} back in {len(bad)} cell(s), e.g. "
                         f"{bad[keys].head(3).to_dict('records')}: stopped")


def compare(analysis: str, fit: pd.DataFrame, frozen: pd.DataFrame,
            changes: dict) -> pd.DataFrame:
    keys = ["tool", "stratum", "side", "tier"]
    m = frozen.rename(columns={c: f"frozen_{c}" for c in
                               ("n_pos", "n_neg", "threshold", "reached", "status")}).merge(
        fit.drop(columns=["cut_lr"]), on=keys, how="inner", validate="one_to_one")
    if len(m) != len(frozen):
        raise SystemExit(f"[C3] {analysis}: the refit does not cover the frozen cells")
    m["tier_change"] = np.select(
        [m["reached"] & ~m["frozen_reached"], ~m["reached"] & m["frozen_reached"]],
        ["gained", "lost"], "unchanged")
    m["threshold_changed"] = [not _same(a, b) for a, b in
                              zip(m["threshold"], m["frozen_threshold"])]
    for col in ("n_added", "n_removed", "n_switched"):
        m[col] = m["stratum"].map({s: c[col] for s, c in changes.items()})
    m.insert(0, "analysis", analysis)
    return m[["analysis", "tool", "stratum", "side", "tier", "cut_lr",
              "n_added", "n_removed", "n_switched",
              "frozen_n_pos", "frozen_n_neg", "n_pos", "n_neg",
              "frozen_threshold", "threshold", "frozen_reached", "reached",
              "tier_change", "threshold_changed", "frozen_status", "status", "reason"]]


def intermediate_labels(df: pd.DataFrame, as_value: float) -> pd.Series:
    """The analysis set's labels, with every unlabelled intermediate-class variant at
    3-50 nt given as_value."""
    y = df["y_assay"].astype("Float64").copy()
    y[(df["y_assay_source"] == INDETERMINATE) & (df["stratum"] != "pm12")] = as_value
    return y


def labels_from_file(df: pd.DataFrame, gene: str, path: Path) -> pd.Series:
    """One gene's labels from a label file, joined as evid_build_set.attach_labels
    joins its gene's current file: the normalised c. notation against the file's key
    column, the class column through the gene's class map (an unmapped class, or a
    variant the file lacks, gives no label)."""
    spec = C.ASSAY_LABELS[gene]
    lab = pd.read_csv(path, sep=None, engine="python")
    lut = dict(zip(lab[spec["key"]].map(_norm_hgvs), lab[spec["class_col"]].map(spec["map"])))
    mask = df["gene"] == gene
    return df.loc[mask, "hgvs_c"].map(_norm_hgvs).map(lut).astype("Float64")


def superseded_labels(df: pd.DataFrame) -> pd.Series:
    y = df["y_assay"].astype("Float64").copy()
    for gene in SUPERSEDED_GENES:
        old = labels_from_file(df, gene, SUPERSEDED_DIR / C.ASSAY_LABELS[gene]["file"])
        y.loc[old.index] = old
    return y


def without_bard1_palb2() -> pd.DataFrame:
    d = K.read_back(DILUTION)
    summ = K.read_back(DILUTION_SUMMARY)
    rows = []
    for (tool, stratum, tier), g in d.groupby(["tool", "stratum", "tier"], sort=False):
        genes = set(g["heldout_gene"])
        for t in g["training_genes"]:
            genes |= set(t.split(";"))
        ref = summ[(summ["tool"] == tool) & (summ["stratum"] == stratum)
                   & (summ["tier"] == tier)]["threshold_all_genes"].unique()
        if len(ref) != 1:
            raise SystemExit(f"[C3] {tool} {stratum} {tier}: no single every-gene threshold")
        ref = float(ref[0])
        every = {"tool": tool, "stratum": stratum, "tier": tier,
                 "path_cut_lr": float(g["path_cut_lr"].iloc[0]), "subset": ALL_GENES,
                 "n_genes": len(genes), "training_genes": ";".join(sorted(genes)),
                 "threshold": ref, "tier_reached": bool(np.isfinite(ref)),
                 "threshold_every_gene": ref, "tier_reached_every_gene": bool(np.isfinite(ref)),
                 "tier_change": "unchanged", "threshold_changed": False}
        rows.append(every)
        for name, drop in SUBSETS:
            want = genes - drop
            hit = g[g["training_genes"].map(lambda t: set(t.split(";")) == want)]
            thr = hit["threshold"].to_numpy(dtype=float)
            reached = set(hit["tier_reached"])
            if len(hit) == 0 or not all(_same(x, thr[0]) for x in thr) or len(reached) != 1:
                raise SystemExit(f"[C3] {tool} {stratum} {tier} {name}: threshold_dilution.csv "
                                 "gives no single fit for this subset: stopped")
            thr, reached = float(thr[0]), reached.pop()
            if bool(reached) != bool(np.isfinite(thr)):
                raise SystemExit(f"[C3] {tool} {stratum} {tier} {name}: tier_reached "
                                 "disagrees with the threshold")
            rows.append(every | {
                "subset": name, "n_genes": len(want), "training_genes": ";".join(sorted(want)),
                "threshold": thr, "tier_reached": bool(reached),
                "tier_change": ("gained" if reached and not np.isfinite(ref) else
                                "lost" if not reached and np.isfinite(ref) else "unchanged"),
                "threshold_changed": not _same(thr, ref)})
    return pd.DataFrame(rows)


def main() -> None:
    cfg = yaml.safe_load(E3.CONFIG_PATH.read_text())
    path_cuts, ben_cuts = K.acmg_bands(cfg)
    cuts = {"path": path_cuts, "ben": ben_cuts}
    df = K.load_set()
    df[K.FUSION] = K.logo_fusion(df, [t for t in K.FUSION_FEATURES if t in df.columns])
    tools = K.panel_of(df) + [K.FUSION]
    frozen = frozen_cells()
    if tools != list(dict.fromkeys(frozen["tool"])):
        raise SystemExit(f"[C3] the columns E3 fits now ({tools}) are not those of {FROZEN}")

    check_reproduces_frozen(fit_cells(df, df["y_assay"].astype("Float64"), tools, cuts), frozen)
    print(f"[C3] the in-sample fit on the set's own labels gives {FROZEN} back exactly")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    parts = []
    for name, value in (("intermediate as damaging", 1.0), ("intermediate as normal", 0.0)):
        y = intermediate_labels(df, value)
        parts.append(compare(name, fit_cells(df, y, tools, cuts), frozen,
                             label_changes(df, y)))
    inter = pd.concat(parts, ignore_index=True)
    inter.to_csv(OUT_INTERMEDIATE, index=False)

    y = superseded_labels(df)
    sup = compare("superseded PALB2 and BARD1 labels", fit_cells(df, y, tools, cuts), frozen,
                  label_changes(df, y))
    sup.to_csv(OUT_SUPERSEDED, index=False)

    wo = without_bard1_palb2()
    wo.to_csv(OUT_WITHOUT, index=False)

    for name, t in (("intermediate", inter), ("superseded", sup)):
        ch = t[t["tier_change"] != "unchanged"]
        print(f"[C3] {name}: {len(ch)} of {len(t)} cells change tier"
              + "".join(f"\n      {r.analysis}: {r.tool} {r.stratum} {r.side} {r.tier} "
                        f"{r.tier_change}" for r in ch.itertuples()))
    ch = wo[wo["tier_change"] != "unchanged"]
    print(f"[C3] without BARD1/PALB2: {len(ch)} of {len(wo[wo.subset != ALL_GENES])} subset "
          "cells change tier" + "".join(f"\n      {r.subset}: {r.tool} {r.stratum} {r.tier} "
                                        f"{r.tier_change}" for r in ch.itertuples()))


if __name__ == "__main__":
    main()
