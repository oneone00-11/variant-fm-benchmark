"""
Phase 1 -- directionality gate (control-anchored, label-independent).

Every SGE assay is anchored so that nonsense/LoF variants are damaging and
synonymous variants are ~neutral. Since func_pathogenicity = -func_score is
supposed to be "larger = more damaging", we can verify orientation PER GENE
without any Method B labels:

    median func_pathogenicity( nonsense )  MUST be  >  median( synonymous )

A gene that fails this needs to go into FLIP_GENES. A gene whose controls barely
separate (AUROC near 0.5) is a data-quality flag worth investigating before it
enters the cross-gene meta-analysis.

Run:  python -m src.phase1_directionality_check
Reads the frozen matrix; writes reports/phase1/directionality_check.csv
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from . import config as C

# region_class values that act as controls. Extend if your annotation differs.
LOF_CONTROLS     = {"nonsense", "stop_gained", "frameshift", "stop_gained_variant"}
NEUTRAL_CONTROLS = {"synonymous", "synonymous_variant", "silent"}

AUROC_STRONG = 0.80   # controls separate cleanly
AUROC_WEAK   = 0.65   # below this: investigate the gene


def _verdict(delta: float, auroc: float) -> str:
    if not np.isfinite(auroc):
        return "NO_CONTROLS"
    if auroc < 0.50 or delta < 0:
        return "FLIP_NEEDED"
    if auroc < AUROC_WEAK:
        return "WEAK_SEPARATION"
    return "OK"


def check(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene, g in df.groupby("gene", observed=True):
        r = g["region"].astype(str).str.lower()
        lof = g.loc[r.isin(LOF_CONTROLS), "func_pathogenicity"].dropna()
        neu = g.loc[r.isin(NEUTRAL_CONTROLS), "func_pathogenicity"].dropna()
        med_lof = float(lof.median()) if len(lof) else np.nan
        med_neu = float(neu.median()) if len(neu) else np.nan
        delta = med_lof - med_neu

        auroc = np.nan
        if len(lof) >= 5 and len(neu) >= 5:
            y = np.r_[np.ones(len(lof)), np.zeros(len(neu))]
            s = np.r_[lof.values, neu.values]
            auroc = roc_auc_score(y, s)   # >0.5 => higher score for LoF (correct)

        rows.append({
            "gene": gene, "n_lof": len(lof), "n_syn": len(neu),
            "median_path_lof": round(med_lof, 4), "median_path_syn": round(med_neu, 4),
            "delta_lof_minus_syn": round(delta, 4) if np.isfinite(delta) else np.nan,
            "control_auroc": round(auroc, 4) if np.isfinite(auroc) else np.nan,
            "verdict": _verdict(delta, auroc),
        })
    return pd.DataFrame(rows).sort_values("gene").reset_index(drop=True)


def main() -> None:
    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    if not p.exists():
        sys.exit(f"[dircheck] frozen matrix not found at {p} -- run phase1 build first")
    df = pd.read_parquet(p)

    out = check(df)
    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(C.REPORT_DIR / "directionality_check.csv", index=False)

    print(out.to_string(index=False))
    flips = out.loc[out["verdict"] == "FLIP_NEEDED", "gene"].tolist()
    weak  = out.loc[out["verdict"].isin(["WEAK_SEPARATION", "NO_CONTROLS"]), "gene"].tolist()
    print("\n[dircheck] GATE SUMMARY")
    print(f"  FLIP_NEEDED -> add to FLIP_GENES and re-freeze: {flips or 'none'}")
    print(f"  investigate before Phase 2: {weak or 'none'}")
    if not flips and not weak:
        print("  ALL GREEN -- orientation consistent; safe to freeze and proceed.")


if __name__ == "__main__":
    main()
