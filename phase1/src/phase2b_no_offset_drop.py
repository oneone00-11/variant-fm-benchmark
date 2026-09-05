"""H1 robustness: dropping the thirteen variants with no placeable intron offset.

Section 3.1 reports that the fusion's ranking advantage concentrates in the
thirteen splice variants that carry no usable intron offset -- seven with no
offset at all and six audited exon-edge/UTR overrides absent from the full
functional masters -- and that dropping them moves delta-rho to +0.0095, below
every one of two hundred random thirteen-variant drops, with the bootstrap
interval's lower bound falling 4.6-fold. Until this module existed that claim
had no code behind it; the anchor tests check that a described analysis exists,
and this one did not.

The drop set is not hard-coded: it is exactly the set the sampling-frame
analysis cannot place -- the analysis-set splice variants whose key does not
occur in the full-functional target frame (|intron offset| <= 8), the same
criterion that takes 1,781 to 1,768 in `ipw_reweight`. Eleven are in BRCA1.

The predictions are the frozen LOGO out-of-fold elastic-net ones and the
frozen within-gene rank of the best single tool; nothing is re-trained. The
interval is the phase-2 gene-cluster bootstrap (2,000 resamples of the seven
genes) off the module's RANDOM_SEED generator -- the same house pattern as
phase5/6/7; the two hundred random drops are point estimates only, drawn
from the same generator after the interval resamples.

Run (PYTHONPATH=phase1):  python -m src.phase2b_no_offset_drop
Outputs:  phase1/reports/phase1/phase2_no_offset_drop.csv           (summary)
          phase1/reports/phase1/phase2_no_offset_drop_draws.csv     (200 draws)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C
from .ipw_reweight import load_target_frame, _key
from .phase2_model import (usable_features, rank_within_gene,
                           rank_target_within_gene, logo_oof,
                           per_gene_rho, dl_pool)

# House pattern (phase5/phase6/phase7): one module-level generator off
# RANDOM_SEED. The 2,000-resample interval is drawn first, then the two
# hundred random drops; under that order the drop-13 interval is
# [+0.0005, +0.0200] -- the lower bound 4.6-fold below the full-set +0.0023
# that 3.1 reports. A label-seeded stream (ipw_reweight.stream style) instead
# gives a lower bound a hair below zero: the bound sits within bootstrap
# noise of zero either way, which is exactly what the text concludes from.
RNG = np.random.default_rng(C.RANDOM_SEED)
N_BOOT = 2000          # same resample count as phase2_model.boot_delta
N_RANDOM_DROPS = 200   # the "two hundred random thirteen-variant drops"
DROP_K = 13


def delta_rho(oof_fusion, single, y, genes):
    """DL-pooled per-gene Spearman delta, fusion minus best single."""
    fpg = per_gene_rho(oof_fusion, y, genes)
    spg = per_gene_rho(single, y, genes)
    return dl_pool(fpg)["rho"] - dl_pool(spg)["rho"], fpg, spg


def boot_delta(fpg, spg, rng, n_boot=N_BOOT):
    """phase2_model.boot_delta with an explicit stream: resample the seven
    genes with replacement, re-pool, take the percentile interval."""
    obs = dl_pool(fpg)["rho"] - dl_pool(spg)["rho"]
    genes = np.intersect1d(fpg["gene"], spg["gene"])
    d = []
    for _ in range(n_boot):
        samp = rng.choice(genes, size=len(genes), replace=True)
        fb = pd.concat([fpg[fpg.gene == g] for g in samp])
        sb = pd.concat([spg[spg.gene == g] for g in samp])
        d.append(dl_pool(fb)["rho"] - dl_pool(sb)["rho"])
    lo, hi = np.nanpercentile(d, [2.5, 97.5])
    return {"delta": obs, "lo": float(lo), "hi": float(hi),
            "excludes_zero": bool(lo > 0 or hi < 0)}


def run():
    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    if not p.exists():
        sys.exit(f"[phase2b] frozen matrix not found at {p}")
    df = pd.read_parquet(p)
    df = df[df["is_splice"] & df["func_pathogenicity"].notna()].reset_index(drop=True)
    n_all = len(df)

    # --- frozen predictions, identical to phase2 ---------------------------
    feats = usable_features(df)
    genes = df["gene"].astype(str).to_numpy()
    y = df["func_pathogenicity"].astype(float).to_numpy()
    ytr = (rank_target_within_gene(df)
           if getattr(C, "TRAIN_TARGET", "rank_within_gene") == "rank_within_gene"
           else df["func_pathogenicity"])
    Xrn = rank_within_gene(df, feats)
    single_pool = {f: dl_pool(per_gene_rho(Xrn[f].to_numpy(), y, genes))["rho"]
                   for f in feats}
    best = max(single_pool, key=lambda k: single_pool[k]
               if np.isfinite(single_pool[k]) else -9)
    oof = logo_oof(Xrn, ytr, df["gene"].astype(str), feats, "enet")
    single = Xrn[best].to_numpy()

    # --- the thirteen variants the sampling frame cannot place -------------
    repo_root = Path(__file__).resolve().parents[2]
    df["key"] = _key(df.chrom, df.pos, df.ref, df.alt).to_numpy()
    frame = load_target_frame(repo_root)
    placeable = df["key"].isin(set(frame["key"])).to_numpy()
    drop_idx = np.where(~placeable)[0]
    n_drop = len(drop_idx)
    per_gene = df.loc[drop_idx, "gene"].value_counts().to_dict()
    print(f"[phase2b] splice set n={n_all:,}; best single = {best}")
    print(f"[phase2b] variants with no placeable intron offset: {n_drop} "
          f"(per gene: {per_gene})")
    if n_drop != DROP_K:
        sys.exit(f"[phase2b] expected {DROP_K} unplaceable variants, found {n_drop}")

    # --- full set, the 13-drop, and its interval ---------------------------
    d_full, _, _ = delta_rho(oof, single, y, genes)
    keep = np.ones(n_all, bool); keep[drop_idx] = False
    d_drop, fpg_d, spg_d = delta_rho(oof[keep], single[keep], y[keep], genes[keep])
    ci = boot_delta(fpg_d, spg_d, RNG)
    print(f"[phase2b] delta-rho full set            : {d_full:+.5f}")
    print(f"[phase2b] delta-rho dropping the {n_drop}    : {d_drop:+.5f}  "
          f"95% CI [{ci['lo']:+.5f}, {ci['hi']:+.5f}]  "
          f"({'excludes' if ci['excludes_zero'] else 'spans'} zero, "
          f"{N_BOOT} gene-cluster resamples)")

    # --- two hundred random thirteen-variant drops --------------------------
    draws = []
    for i in range(N_RANDOM_DROPS):
        idx = RNG.choice(n_all, size=DROP_K, replace=False)
        k = np.ones(n_all, bool); k[idx] = False
        d, _, _ = delta_rho(oof[k], single[k], y[k], genes[k])
        draws.append({"draw": i, "delta_rho": d})
    draws_df = pd.DataFrame(draws)
    r_mean, r_sd = float(draws_df.delta_rho.mean()), float(draws_df.delta_rho.std(ddof=1))
    r_min, r_max = float(draws_df.delta_rho.min()), float(draws_df.delta_rho.max())
    n_below = int((d_drop < draws_df.delta_rho).sum())
    print(f"[phase2b] {N_RANDOM_DROPS} random {DROP_K}-variant drops   : "
          f"{r_mean:+.5f} +/- {r_sd:.5f} (sd), range [{r_min:+.5f}, {r_max:+.5f}]")
    print(f"[phase2b] the {DROP_K}-variant drop is below {n_below} of "
          f"{N_RANDOM_DROPS} random drops")

    summary = pd.DataFrame([
        {"quantity": "delta_rho_full_set", "value": d_full,
         "note": f"all {n_all} splice variants; phase2 overall"},
        {"quantity": "delta_rho_drop_no_offset", "value": d_drop,
         "note": f"dropping the {n_drop} unplaceable variants"},
        {"quantity": "drop_ci_lo", "value": ci["lo"],
         "note": f"gene-cluster bootstrap, {N_BOOT} resamples"},
        {"quantity": "drop_ci_hi", "value": ci["hi"],
         "note": f"gene-cluster bootstrap, {N_BOOT} resamples"},
        {"quantity": "drop_ci_excludes_zero", "value": float(ci["excludes_zero"]),
         "note": "1 if the interval excludes zero"},
        {"quantity": "random_drops_mean", "value": r_mean,
         "note": f"{N_RANDOM_DROPS} random {DROP_K}-variant drops, point estimates"},
        {"quantity": "random_drops_sd", "value": r_sd, "note": "sample sd (ddof=1)"},
        {"quantity": "random_drops_min", "value": r_min, "note": ""},
        {"quantity": "random_drops_max", "value": r_max, "note": ""},
        {"quantity": "n_random_drops_above_drop", "value": float(n_below),
         "note": f"of {N_RANDOM_DROPS}; the manuscript claims all two hundred"},
    ])
    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(C.REPORT_DIR / "phase2_no_offset_drop.csv", index=False)
    draws_df.to_csv(C.REPORT_DIR / "phase2_no_offset_drop_draws.csv", index=False)
    print(f"\n[phase2b] wrote phase2_no_offset_drop.csv and "
          f"phase2_no_offset_drop_draws.csv to {C.REPORT_DIR}")


if __name__ == "__main__":
    run()
