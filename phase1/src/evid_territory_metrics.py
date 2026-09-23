"""E4 -- classification metrics by territory and by ClinVar arm.

One table, two questions.

  1. Is 11-50 bp the worst-served territory on classification metrics too, or only
     on rank correlation? The companion atlas located the genuine prediction
     failure there on rho against a reliability ceiling. AUROC and PR-AUC are a
     different question about the same variants, and a territory map that holds on
     one metric and not the other is worth knowing before either is built on.
  2. Does evidence fall monotonically across the three ClinVar arms -- classified,
     recorded but unclassified, not recorded at all? The published study could only
     contrast the first two, because every variant in its set carried a record.

AUROC is pooled across genes exactly as the companion atlas pools it (logit scale,
Hanley-McNeil variances, DerSimonian-Laird), by importing that code rather than
reimplementing it, so the two papers' AUROC columns are the same estimator.
Spearman rho is carried as a secondary column, not as the endpoint.

Run (PYTHONPATH=phase1):  python -m src.evid_territory_metrics
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import evid_common as K

from .evid_common import ATLAS_REPO as _atlas_repo_default  # noqa: E402
ATLAS_REPO = _atlas_repo_default  # EVID_ATLAS_REPO overrides; see evid_common
REPORT_DIR = Path("reports/evidence")
OUT = REPORT_DIR / "territory_metrics.csv"

MIN_POS = MIN_NEG = 10
TIER_RANK = {"below supporting": 0, "Supporting": 1, "Moderate": 2,
             "Strong": 3, "Very strong": 4, "not evaluable": -1}


def _atlas_clinical_evidence():
    """Import the atlas's clinical-evidence module by path, so AUROC pooling here
    is literally the same code as the companion paper's Figure 6."""
    src = str(ATLAS_REPO / "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    spec = importlib.util.spec_from_file_location(
        "atlas_clinical_evidence", ATLAS_REPO / "src/atlas/clinical_evidence.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def per_gene_and_pool(sub: pd.DataFrame, tool: str, ce) -> dict:
    per_gene = []
    for gene, g in sub.groupby("gene", observed=True):
        y = g["y_assay"].to_numpy(dtype=float)
        s = g[tool].to_numpy(dtype=float)
        m = np.isfinite(y) & np.isfinite(s)
        y, s = y[m], s[m]
        n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
        if n_pos < MIN_POS or n_neg < MIN_NEG:
            continue
        auc = ce.roc_auc(y, s)
        per_gene.append({"gene": gene, "n_pos": n_pos, "n_neg": n_neg,
                         "auroc": auc, "auprc": ce.average_precision(y, s),
                         "rho": spearmanr(s, y).statistic,
                         "var": ce.hanley_mcneil_var(auc, n_pos, n_neg)})
    if len(per_gene) < 2:
        return {"k_genes": len(per_gene), "status": "not evaluable",
                "reason": "fewer than two genes clear the minimum count"}
    pg = pd.DataFrame(per_gene)
    z = pg["auroc"].map(ce._logit)
    vz = pg["var"] / (pg["auroc"] * (1 - pg["auroc"])) ** 2
    p = ce.dl_pool(z.tolist(), vz.tolist())
    inv = lambda x: 1 / (1 + math.exp(-x))  # noqa: E731
    # Genes below the minimum count are dropped from the estimate. The row's `n`
    # and `positive_rate` describe the whole cell, so they are NOT the sample the
    # estimate was computed on; the retained figures are reported beside them rather
    # than leaving a cell size next to an estimate that did not use it.
    return {
        "k_genes": len(pg),
        "n_pos": int(pg["n_pos"].sum()), "n_neg": int(pg["n_neg"].sum()),
        "n_used": int(pg["n_pos"].sum() + pg["n_neg"].sum()),
        "positive_rate_used": float(pg["n_pos"].sum()
                                    / (pg["n_pos"].sum() + pg["n_neg"].sum())),
        "auroc": inv(p["pooled"]), "auroc_lo": inv(p["ci_lo"]),
        "auroc_hi": inv(p["ci_hi"]), "auroc_i2_pct": p["i2_pct"],
        "auprc_median": float(pg["auprc"].median()),
        "rho_median": float(pg["rho"].median()),
        "status": "ok", "reason": "",
    }


def attach_best_tier(out: pd.DataFrame) -> pd.DataFrame:
    """Highest evidence tier E3 finds reachable, joined on tool x stratum.

    E3 determines thresholds on the whole stratum, not per ClinVar arm, so the
    column repeats down the arms and is labelled as a stratum-level property.
    """
    path = REPORT_DIR / "evidence_thresholds.csv"
    if not path.exists():
        out["e3_best_tier_pp3"] = "E3 not run"
        return out
    ev = pd.read_csv(path)
    ev = ev[(ev.get("status") == "ok") & ev["pp3_threshold_reachable"].fillna(False)]
    order = {"supporting": 1, "moderate": 2, "strong": 3, "very_strong": 4}
    ev["rank"] = ev["tier"].map(order)
    best = (ev.sort_values("rank").groupby(["tool", "stratum"], observed=True)
              .tail(1)[["tool", "stratum", "tier", "pp3_threshold_insample",
                        "pp3_logo_heldout_lr_median"]]
              .rename(columns={"tier": "e3_best_tier_pp3",
                               "pp3_threshold_insample": "e3_threshold_at_best_tier",
                               "pp3_logo_heldout_lr_median": "e3_logo_heldout_lr_median"}))
    return out.merge(best, on=["tool", "stratum"], how="left")


def main() -> None:
    ce = _atlas_clinical_evidence()
    df = K.load_set()
    tools = K.panel_of(df)
    if K.FUSION not in df.columns:
        df[K.FUSION] = K.logo_fusion(df, [t for t in K.FUSION_FEATURES if t in df.columns])
    tools = tools + [K.FUSION]

    rows = []
    for tool in tools:
        for stratum in K.STRATA:
            ss = K.stratum_frame(df, stratum)
            for arm in K.ARMS:
                sub = K.arm_frame(ss, arm).dropna(subset=["y_assay", tool])
                base = {"tool": tool, "stratum": stratum, "clinvar_arm": arm,
                        "n_cell": int(len(sub)),
                        "positive_rate_cell": (float((sub["y_assay"] == 1).mean())
                                               if len(sub) else np.nan)}
                rows.append(base | per_gene_and_pool(sub, tool, ce))

    out = pd.DataFrame(rows)
    out = attach_best_tier(out)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)

    ok = out[out.status == "ok"]
    print(f"[E4] wrote {OUT} ({len(out)} rows; {len(ok)} evaluable)\n")
    print("--- AUROC by stratum, all arms pooled ---")
    print(ok[ok.clinvar_arm == "all"].pivot_table(
        index="tool", columns="stratum", values="auroc").round(3).to_string())
    print("\n--- AUROC by ClinVar arm, all_1_50 ---")
    print(ok[ok.stratum == "all_1_50"].pivot_table(
        index="tool", columns="clinvar_arm", values="auroc").round(3).to_string())
    print("\n--- PR-AUC (median over genes) by stratum, all arms ---")
    print(ok[ok.clinvar_arm == "all"].pivot_table(
        index="tool", columns="stratum", values="auprc_median").round(3).to_string())
    ne = out[out.status != "ok"]
    print(f"\n[E4] not evaluable: {len(ne)} cells")
    if len(ne):
        print(ne.groupby(["stratum", "clinvar_arm"], observed=True).size()
                .reset_index(name="n").to_string(index=False))


if __name__ == "__main__":
    main()
