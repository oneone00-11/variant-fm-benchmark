"""E17 -- the two AlphaGenome definitions side by side.

E17 replaced the study's `alphagenome` column. It was the companion atlas's
merged-quantile composite (client 0.7.0, a 16-kb window). It is now AlphaGenome's own
definition: the merged splicing score of AlphaGenome's splicing notebook, on raw scores
over a 1 Mb window (client 0.10.0; data/evidence/alphagenome_official/). The old column
stays in the analysis set as `alphagenome_atlas_quantile`, outside the panel, and this
table is the one place it appears, beside the new one.

One row per stratum (s3_10, s11_50, s3_50) x tier (Supporting, Moderate, Strong) x side
(PP3, BP4). For each definition it gives only quantities the study reports for every
predictor elsewhere, computed by the functions that report them there:

  * the in-sample threshold: E3's (evid_interval_lr.run_tool_stratum), as Tables S6
    and S7 print it;
  * the band likelihood ratio at that threshold: E15's (evid_practical_metrics), as the
    practical metrics give it;
  * held out, on the PP3 side: Table 3's "k/n (m)", the held-out genes whose ratio
    clears the tier's cut out of the folds, with the folds not evaluable in brackets,
    shown where the tier is reached in-sample (E3's per-fold thresholds, counted as
    E2.2 counts them). The BP4 side has no count, since E3 writes per-fold thresholds
    for the PP3 side only; E3's median held-out ratio is given instead, as Table S7
    gives it;
  * the stratum's pooled AUROC: E4's (evid_territory_metrics.per_gene_and_pool, every
    ClinVar arm).

The old definition is run under the name it had, `alphagenome`, so that E3's bootstrap
is seeded as it was and its numbers are the ones the study printed before E17;
scripts/stage_E17/check_alphagenome_only.py checks them against that commit's outputs.
The new definition's numbers are checked here against what E3, E2.2, E4 and E15 wrote
for `alphagenome` in the same run, and a disagreement stops the stage.

Run (PYTHONPATH=phase1):  python -m src.evid_alphagenome_definitions
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import evid_common as K
from . import evid_interval_lr as L
from . import evid_practical_metrics as P
from . import evid_territory_metrics as T
from .evid_tables import _folds

REPORT_DIR = Path("reports/evidence")
OUT = REPORT_DIR / "alphagenome_definition_comparison.csv"
CONFIG_PATH = Path("config/walker2023.yaml")
STRATA = ["s3_10", "s11_50", "s3_50"]
TIERS = ["supporting", "moderate", "strong"]
TOOL = "alphagenome"
OLD = "alphagenome_atlas_quantile"
LR_COLUMN = {"pp3": "band_lr_plus", "bp4": "bp4_band_lr"}


def one_definition(df: pd.DataFrame, cfg: dict, ce) -> pd.DataFrame:
    """Every quantity of the table for the column df[TOOL], by stratum, tier and side."""
    path_cuts, ben_cuts = K.acmg_bands(cfg)
    rows = []
    for st in STRATA:
        res, _, logo = L.run_tool_stratum(df, TOOL, st, path_cuts, ben_cuts, cfg,
                                          L.N_BOOT_VARIANT)
        e = {r["tier"]: r for r in res}
        frame = K.stratum_frame(df, st)
        auc = T.per_gene_and_pool(K.arm_frame(frame, "all").dropna(subset=["y_assay", TOOL]),
                                  TOOL, ce)
        for tier in TIERS:
            r = e[tier]
            if r["status"] != "ok":
                rows += [{"stratum": st, "tier": tier, "side": side, "status": "not evaluable"}
                         for side in ("pp3", "bp4")]
                continue
            # E2.2's counts over E3's per-fold thresholds, formatted as Table 3 formats them
            f = (logo[(logo.tier == tier) & (logo.side == "pp3")] if logo is not None
                 and len(logo) else pd.DataFrame(columns=["heldout_lr"]))
            lrs = f["heldout_lr"].astype(float)
            t = pd.DataFrame([{
                "tool": TOOL, "stratum": st, "tier": tier, "status": "ok",
                "in_sample_tier_reached": bool(r["pp3_threshold_reachable"]),
                "n_folds": len(f), "folds_with_heldout_lr": int(lrs.notna().sum()),
                "folds_heldout_lr_above_cut": int((lrs >= r["path_cut_lr"]).sum())}])
            for side in ("pp3", "bp4"):
                tau = float(r[f"{side}_threshold_insample"])
                reached = bool(r[f"{side}_threshold_reachable"])
                lr = (P._values(P._cell(frame, TOOL, "y_assay", tau), side,
                                (path_cuts, ben_cuts))[LR_COLUMN[side]]
                      if reached else np.nan)
                rows.append({
                    "stratum": st, "tier": tier, "side": side, "status": "ok",
                    "reached_in_sample": reached,
                    "in_sample_threshold": tau,
                    "band_lr_at_threshold": lr,
                    "heldout_k_of_n": _folds(t, TOOL, st, tier) if side == "pp3" else "",
                    "heldout_lr_median": (float(r["bp4_logo_heldout_lr_median"])
                                          if side == "bp4" and reached else np.nan),
                    "pooled_auroc": auc.get("auroc", np.nan),
                })
    return pd.DataFrame(rows)


def check_against_study(new: pd.DataFrame, df: pd.DataFrame) -> None:
    """The new definition's numbers must be the ones E3, E2.2, E4 and E15 wrote for
    `alphagenome` in this run, float for float."""
    ev = K.read_back(REPORT_DIR / "evidence_thresholds.csv")
    tl = K.read_back(REPORT_DIR / "tier_logo_table.csv")
    tm = K.read_back(REPORT_DIR / "territory_metrics.csv")
    pm = K.read_back(P.OUT_TABLE)
    bad = []

    def same(a, b) -> bool:
        return (pd.isna(a) and pd.isna(b)) or a == b

    for r in new.itertuples():
        if r.status != "ok":
            continue
        e = ev[(ev.tool == TOOL) & (ev.stratum == r.stratum) & (ev.tier == r.tier)].iloc[0]
        checks = [(r.in_sample_threshold, e[f"{r.side}_threshold_insample"]),
                  (r.reached_in_sample, bool(e[f"{r.side}_threshold_reachable"]))]
        if r.side == "bp4" and r.reached_in_sample:
            checks.append((r.heldout_lr_median, e["bp4_logo_heldout_lr_median"]))
        if r.side == "pp3":
            t = tl[(tl.tool == TOOL) & (tl.stratum == r.stratum) & (tl.tier == r.tier)]
            checks.append((r.heldout_k_of_n, _folds(t, TOOL, r.stratum, r.tier)))
        a = tm[(tm.tool == TOOL) & (tm.stratum == r.stratum) & (tm.clinvar_arm == "all")]
        checks.append((r.pooled_auroc, a["auroc"].iloc[0]))
        if r.reached_in_sample and (r.side == "pp3" or r.tier in P.BP4_TIERS):
            kind = r.tier if r.side == "pp3" else f"bp4_{r.tier}"
            p = pm[(pm.tool == TOOL) & (pm.stratum == r.stratum) & (pm.side == r.side)
                   & (pm.threshold_kind == kind) & (pm.threshold_source == "in_sample")]
            checks.append((r.band_lr_at_threshold, p[LR_COLUMN[r.side]].iloc[0]))
        for a_, b_ in checks:
            if not same(a_, b_):
                bad.append(f"{r.stratum} {r.tier} {r.side}: {a_!r} here, {b_!r} in the study")
    if bad:
        raise SystemExit("[E17] the new definition's numbers differ from the study's own "
                         "files:\n  " + "\n  ".join(bad))


def build() -> pd.DataFrame:
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    ce = T._atlas_clinical_evidence()
    df = K.load_set()
    if OLD not in df.columns:
        raise SystemExit(f"[E17] the analysis set has no {OLD} column; run src.evid_build_set")
    new = one_definition(df, cfg, ce)
    check_against_study(new, df)
    # the old column under the name it had, so that E3's bootstrap is seeded as it was
    old = one_definition(df.assign(**{TOOL: df[OLD]}), cfg, ce)
    path_cuts, ben_cuts = K.acmg_bands(cfg)
    key = ["stratum", "tier", "side"]
    out = new[key].copy()
    out["tier_lr_cut"] = [(path_cuts if s == "pp3" else ben_cuts)[t]
                          for t, s in zip(out.tier, out.side)]
    for name, d in ((TOOL, new), (OLD, old)):
        if list(d[key].itertuples(index=False)) != list(new[key].itertuples(index=False)):
            raise SystemExit("[E17] the two definitions' rows are not in the same order")
        for c in d.columns:
            if c not in key:
                out[f"{name}_{c}"] = d[c].to_numpy()
    out["side"] = out["side"].str.upper()
    out["note"] = np.where(out.side == "BP4", P.BP4_NO_LOGO, "")
    return out


def main() -> None:
    out = build()
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"[E17] wrote {OUT} ({len(out)} rows)\n")
    show = ["stratum", "tier", "side"] + [f"{n}_{c}" for n in (TOOL, OLD) for c in
                                          ("in_sample_threshold", "band_lr_at_threshold",
                                           "heldout_k_of_n")]
    print(out[show].to_string(index=False))


if __name__ == "__main__":
    main()
