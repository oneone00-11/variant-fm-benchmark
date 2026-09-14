"""(h) Per-gene head-to-head differences behind the "six of seven genes" sentences.

The manuscript says which held-out genes favour the fusion over the best single
tool on each axis (3.1: ranking, Δρ per gene; 3.3: calibration, ΔBrier per gene in
the primary condition). Those counts and the name of the dissenting gene were read
off Figure 1 and the sign-flip inputs; nothing tracked them. This stage writes them.

    delta_rho    = per-gene Spearman ρ, fusion (M1, LOGO out-of-fold) − best single tool
    delta_brier  = per-gene Brier, best single − fusion, after the identical out-of-gene
                   isotonic calibration on the functional standard, BRCA1 included
                   (positive favours the fusion on both axes)

Same features, fusion, calibration and best-single criterion as phases 2/3/5/7/8.

Usage (PYTHONPATH=phase1):  python -m src.phase8_per_gene_deltas
Writes reports/phase1/phase8_per_gene_deltas.csv
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from . import config as C
from .phase2_model import (usable_features, rank_within_gene, rank_target_within_gene,
                           logo_oof, per_gene_rho, dl_pool)
from .phase3_calibration import logo_calibrate, brier


def run():
    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    if not p.exists():
        sys.exit(f"[phase8h] frozen matrix not found at {p}")
    df = pd.read_parquet(p)
    df = df[df["is_splice"] & df["func_pathogenicity"].notna()].reset_index(drop=True)
    feats = usable_features(df)
    genes = df["gene"].astype(str)
    gv = genes.to_numpy()
    Xrn = rank_within_gene(df, feats)
    ytr = (rank_target_within_gene(df)
           if getattr(C, "TRAIN_TARGET", "rank_within_gene") == "rank_within_gene"
           else df["func_pathogenicity"])
    y = df["func_pathogenicity"].astype(float).to_numpy()

    oof = logo_oof(Xrn, ytr, genes, feats, "enet")
    single_pg = {f: per_gene_rho(Xrn[f].to_numpy(), y, gv) for f in feats}
    pool = {f: dl_pool(single_pg[f])["rho"] for f in feats}
    best = max(pool, key=lambda k: pool[k] if np.isfinite(pool[k]) else -9)
    fus_pg = per_gene_rho(oof, y, gv).set_index("gene")
    bs_pg = single_pg[best].set_index("gene")

    # min-max as phase3 does before calibrating
    def mm(s):
        lo, hi = np.nanmin(s), np.nanmax(s)
        return (s - lo) / (hi - lo) if hi > lo else s
    lab = df["y_assay"].astype(float).to_numpy()
    pf = logo_calibrate(mm(oof), lab, gv, method="isotonic")
    ps = logo_calibrate(mm(Xrn[best].to_numpy()), lab, gv, method="isotonic")

    rows = []
    for g in sorted(fus_pg.index):
        m = (gv == g) & ~np.isnan(pf) & ~np.isnan(ps) & ~np.isnan(lab)
        bf, bsg = brier(pf[m], lab[m]), brier(ps[m], lab[m])
        rows.append({"gene": g, "best_single": best,
                     "n_rho": int(fus_pg.loc[g, "n"]),
                     "rho_fusion": float(fus_pg.loc[g, "rho"]),
                     "rho_best_single": float(bs_pg.loc[g, "rho"]),
                     "delta_rho": float(fus_pg.loc[g, "rho"] - bs_pg.loc[g, "rho"]),
                     "n_brier": int(m.sum()),
                     "brier_fusion": float(bf), "brier_best_single": float(bsg),
                     "delta_brier": float(bsg - bf)})
    out = pd.DataFrame(rows)
    out["rho_favours_fusion"] = out["delta_rho"] > 0
    out["brier_favours_fusion"] = out["delta_brier"] > 0
    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(C.REPORT_DIR / "phase8_per_gene_deltas.csv", index=False)
    print(f"[phase8h] {len(out)} genes; best single {best}; "
          f"ranking favours fusion in {int(out.rho_favours_fusion.sum())} "
          f"(dissenting: {', '.join(out.gene[~out.rho_favours_fusion]) or 'none'}); "
          f"calibration favours fusion in {int(out.brier_favours_fusion.sum())} "
          f"(dissenting: {', '.join(out.gene[~out.brier_favours_fusion]) or 'none'})")
    print(out.round(4).to_string(index=False))


if __name__ == "__main__":
    run()
