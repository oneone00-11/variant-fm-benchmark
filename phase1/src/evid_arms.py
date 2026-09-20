"""E2.3 / E2.4 -- the ClinVar arms, with the gene confound taken out two ways.

The pooled arm comparison is confounded with gene composition and the size of the
confound is not small: the unrecorded arm holds no BRCA1 at all, because BRCA1's
saturation editing data were deposited into ClinVar, and BRCA1 supplies more than
half of the recorded-but-unclassified arm. A difference between arms is therefore
also a difference between genes.

Two ways out, both reported:

  1. drop BRCA1 and redo the pooled comparison. Cheap, and it answers whether the
     ordering survives without the gene that causes the imbalance.
  2. compare arms WITHIN each gene. Only this removes the confound rather than
     reducing it, because every gene then contributes to every arm it has.

E2.4 also makes the 588 ClinVar records that carry no clinical assertion visible.
They are ClinVar entries with an empty CLNSIG -- an allele record nobody has
assessed -- and the collapse rule files them under "Other", which the three-arm
rule then reads as recorded-but-unclassified. 576 of the 588 are BRCA1, 85% of that
gene's middle arm. The three arms the brief specifies are unchanged; the strict
split is reported beside them.

Run (PYTHONPATH=phase1):  python -m src.evid_arms
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import evid_common as K

ATLAS_REPO = Path(os.environ.get(
    "EVID_ATLAS_REPO", "/Users/cliffzhang/work/functional-standard-atlas"))
REPORT_DIR = Path("reports/evidence")
CONFIG_PATH = Path("config/walker2023.yaml")
MIN_POS = MIN_NEG = 10
ARMS = ["classified", "recorded_unclassified", "unrecorded"]
STRATA = ["s3_10", "s11_50", "s3_50", "all_1_50"]


def _ce():
    src = str(ATLAS_REPO / "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    spec = importlib.util.spec_from_file_location(
        "atlas_clinical_evidence", ATLAS_REPO / "src/atlas/clinical_evidence.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pooled_auroc(sub: pd.DataFrame, tool: str, ce) -> dict:
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
        per_gene.append({"gene": gene, "n_pos": n_pos, "n_neg": n_neg, "auroc": auc,
                         "auprc": ce.average_precision(y, s),
                         "var": ce.hanley_mcneil_var(auc, n_pos, n_neg)})
    if len(per_gene) < 2:
        return {"k_genes": len(per_gene), "status": "not evaluable"}
    pg = pd.DataFrame(per_gene)
    z = pg["auroc"].map(ce._logit)
    vz = pg["var"] / (pg["auroc"] * (1 - pg["auroc"])) ** 2
    p = ce.dl_pool(z.tolist(), vz.tolist())
    inv = lambda x: 1 / (1 + math.exp(-x))          # noqa: E731
    return {"k_genes": len(pg), "n_pos": int(pg.n_pos.sum()), "n_neg": int(pg.n_neg.sum()),
            "auroc": inv(p["pooled"]), "auroc_lo": inv(p["ci_lo"]),
            "auroc_hi": inv(p["ci_hi"]), "auprc_median": float(pg.auprc.median()),
            "status": "ok"}


def walker_lr(sub: pd.DataFrame, tool: str, pp3: float, bp4: float) -> dict:
    y = sub["y_assay"].to_numpy(dtype=float)
    s = sub[tool].to_numpy(dtype=float)
    m = np.isfinite(y) & np.isfinite(s)
    y, s = y[m], s[m]
    n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
    if n_pos < MIN_POS or n_neg < MIN_NEG:
        return {"n_pos": n_pos, "n_neg": n_neg, "status": "not evaluable"}
    return {"n_pos": n_pos, "n_neg": n_neg, "status": "ok",
            "sens_pp3": float((s[y == 1] >= pp3).mean()),
            "spec_pp3": float((s[y == 0] < pp3).mean()),
            "lr_pp3": K.band_lr(y, s, pp3, "upper"),
            "lr_bp4": K.band_lr(y, s, bp4, "lower"),
            "lr_bp4_for_tier": K.band_lr_benign_bound(y, s, bp4)}


def main() -> None:
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    pp3 = cfg["thresholds"]["pp3"]["value"]
    bp4 = cfg["thresholds"]["bp4"]["value"]
    path_bands, ben_bands = K.acmg_bands(cfg)
    ce = _ce()
    df = K.load_set()
    tools = K.panel_of(df)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    # ---- E2.4: where the no-assertion records are ---------------------------
    strict = (df.groupby(["clinvar_arm", "clinvar_arm_strict"], observed=True)
                .size().reset_index(name="n"))
    strict.to_csv(REPORT_DIR / "clinvar_arm_no_assertion.csv", index=False)
    noass = df[df["clinvar_arm_strict"] == "recorded_no_assertion"]
    print(f"[E2.4] {len(noass)} ClinVar records carry no clinical assertion "
          f"(empty CLNSIG). They sit in the `recorded_unclassified` arm, which is "
          f"{len(df[df.clinvar_arm == 'recorded_unclassified'])} variants, so they "
          f"are {len(noass) / max(len(df[df.clinvar_arm == 'recorded_unclassified']), 1):.0%} of it.")
    print("     by gene:", noass.gene.value_counts().to_dict())
    print("     by stratum:", noass.stratum.value_counts().to_dict())

    # ---- E2.3a: pooled, with and without BRCA1 ------------------------------
    rows = []
    for excl, label in ((False, "all_genes"), (True, "no_BRCA1")):
        base = df[df.gene != "BRCA1"] if excl else df
        for stratum in STRATA:
            ss = K.stratum_frame(base, stratum)
            for arm in ARMS:
                sa = ss[ss.clinvar_arm == arm]
                for tool in tools:
                    sub = sa.dropna(subset=["y_assay", tool])
                    r = {"gene_set": label, "stratum": stratum, "clinvar_arm": arm,
                         "tool": tool, "n": int(len(sub))}
                    rows.append(r | pooled_auroc(sub, tool, ce)
                                | {f"walker_{k}": v
                                   for k, v in walker_lr(sub, tool, pp3, bp4).items()
                                   if k not in ("n_pos", "n_neg", "status")})
    arms_pooled = pd.DataFrame(rows)
    arms_pooled["walker_tier_pp3"] = arms_pooled.get("walker_lr_pp3", pd.Series(dtype=float)).map(
        lambda v: K.tier_of(v, path_bands, "pathogenic"))
    arms_pooled.to_csv(REPORT_DIR / "territory_metrics_arms_noBRCA1.csv", index=False)

    # ---- E2.3b: within gene -------------------------------------------------
    rows = []
    for gene, gdf in df.groupby("gene", observed=True):
        for stratum in STRATA:
            ss = K.stratum_frame(gdf, stratum)
            for arm in ARMS:
                sa = ss[ss.clinvar_arm == arm]
                for tool in tools:
                    sub = sa.dropna(subset=["y_assay", tool])
                    y = sub["y_assay"].to_numpy(dtype=float)
                    s = sub[tool].to_numpy(dtype=float)
                    n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
                    r = {"gene": gene, "stratum": stratum, "clinvar_arm": arm,
                         "tool": tool, "n": int(len(sub)),
                         "n_pos": n_pos, "n_neg": n_neg}
                    if n_pos < MIN_POS or n_neg < MIN_NEG:
                        rows.append(r | {"status": "not evaluable"})
                        continue
                    rows.append(r | {"status": "ok", "auroc": ce.roc_auc(y, s),
                                     "auprc": ce.average_precision(y, s)}
                                | {f"walker_{k}": v
                                   for k, v in walker_lr(sub, tool, pp3, bp4).items()
                                   if k not in ("n_pos", "n_neg", "status")})
    within = pd.DataFrame(rows)
    within.to_csv(REPORT_DIR / "arms_within_gene.csv", index=False)

    # ---- console ------------------------------------------------------------
    ok = arms_pooled[arms_pooled.status == "ok"]
    print("\n[E2.3] pooled AUROC by arm, in-scope pool s3_50")
    for label in ("all_genes", "no_BRCA1"):
        t = ok[(ok.gene_set == label) & (ok.stratum == "s3_50")]
        if len(t):
            print(f"  --- {label} ---")
            print(t.pivot_table(index="tool", columns="clinvar_arm",
                                values="auroc").round(3).to_string())
    w = within[within.status == "ok"]
    print(f"\n[E2.3] within-gene: {len(w)} evaluable cells; "
          f"genes contributing all three arms in s3_50: "
          f"{sorted(w[w.stratum == 's3_50'].groupby('gene')['clinvar_arm'].nunique().pipe(lambda s: s[s == 3]).index)}")
    print(f"\n[E2.3] wrote territory_metrics_arms_noBRCA1.csv ({len(arms_pooled)} rows), "
          f"arms_within_gene.csv ({len(within)} rows), "
          f"clinvar_arm_no_assertion.csv")


if __name__ == "__main__":
    main()
