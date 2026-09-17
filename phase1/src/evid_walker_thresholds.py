"""E2 -- the ClinGen SVI (Walker 2023) fixed cut points, applied to seven genes.

The cut points are not ours and are not fitted here: SpliceAI max raw delta >= 0.2
for PP3 and <= 0.1 for BP4, with the interval between them declared uninformative.
phase1/config/walker2023.yaml holds them with their source. What this stage
measures is what those fixed cut points deliver on a functional standard, by
intronic distance and by ClinVar arm.

Two things to read carefully in the output.

`lr_pp3` and `lr_bp4` are BAND likelihood ratios -- P(band | damaging) divided by
P(band | normal) -- because that is the quantity Walker's own table reports and the
quantity a laboratory applies when it says "this variant scored <= 0.1". For the
>= 0.2 band this coincides with the classical LR+ at a 0.2 threshold. For the
<= 0.1 band it does NOT coincide with the classical LR- at 0.2, which would fold
the uninformative middle in with the benign call; that classical value is carried
alongside as `lr_minus_at_pp3` so both readings are visible.

Prior sensitivity is reported as a posterior probability, not as a second set of
bands. Tavtigian et al. publish a solution for the odds constant at a prior of 0.10
only (OPVSt = 350; their Figure 1 shows the admissible region for other priors but
tabulates no values), so re-banding at 0.05 or 0.20 would mean fitting a constant
ourselves. The likelihood ratio does not depend on the prior; the posterior does,
and it is computed directly.

Run (PYTHONPATH=phase1):  python -m src.evid_walker_thresholds
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import config as C

SET_PATH = Path("data/evidence/analysis_set_v1.parquet")
CONFIG_PATH = Path("config/walker2023.yaml")
REPORT_DIR = Path("reports/evidence")
OUT = REPORT_DIR / "walker_thresholds.csv"

N_BOOT = 2000
MIN_POS = MIN_NEG = 10
# A band also needs occupants before it gets a ratio. Without this a band holding a
# single variant yields a finite likelihood ratio and an evidence tier: pm12's BP4
# band holds one variant, zero of them damaging, which is arithmetically a ratio of
# zero and reads as Very strong benign evidence from one observation.
MIN_BAND = 10
STRATUM_ORDER = ["pm12", "s3_10", "s11_50", "s3_50", "all_1_50"]
# The recommendation excludes the canonical +/-1,2 dinucleotides, so a pool that
# contains them is not a place to read a Walker-threshold result off. On this set
# pm12 supplies half the positives of the full pool.
IN_SCOPE = {"s3_10", "s11_50", "s3_50"}
ARM_ORDER = ["all", "classified", "recorded_unclassified", "unrecorded"]
PRIORS = [0.05, 0.10, 0.20]

def _rng_for(*key) -> np.random.Generator:
    """A generator seeded from the cell's identity, not from call order.

    A single module-level generator makes every row depend on how many rows ran
    before it, so one cell cannot be reproduced without re-running the whole table
    in the same order. Deriving the seed from the cell key makes each row
    independently reproducible and leaves the table unchanged if a row is added.
    """
    h = hashlib.sha256(("|".join(str(k) for k in key)).encode()).digest()
    return np.random.default_rng([C.RANDOM_SEED, int.from_bytes(h[:8], "big")])


def load_config() -> dict:
    return yaml.safe_load(CONFIG_PATH.read_text())


def tier_pathogenic(lr: float, bands: dict) -> str:
    if lr is None or not np.isfinite(lr):
        return "not evaluable"
    for name, thresh in [("Very strong", bands["very_strong"]), ("Strong", bands["strong"]),
                         ("Moderate", bands["moderate"]), ("Supporting", bands["supporting"])]:
        if lr >= thresh:
            return name
    return "below supporting"


def tier_benign(lr: float, bands: dict) -> str:
    if lr is None or not np.isfinite(lr):
        return "not evaluable"
    for name, thresh in [("Very strong", bands["very_strong"]), ("Strong", bands["strong"]),
                         ("Moderate", bands["moderate"]), ("Supporting", bands["supporting"])]:
        if lr <= thresh:
            return name
    return "below supporting"


def posterior(lr: float, prior: float) -> float:
    if lr is None or not np.isfinite(lr):
        return np.nan
    return float(lr * prior / ((lr - 1.0) * prior + 1.0))


def _band_stats(y: np.ndarray, s: np.ndarray, pp3: float, bp4: float) -> dict:
    """Band proportions and band likelihood ratios. No guarding, no clipping --
    the caller decides what an empty cell means."""
    pos, neg = s[y == 1], s[y == 0]
    n_pos, n_neg = len(pos), len(neg)
    if n_pos < MIN_POS or n_neg < MIN_NEG:
        return {}
    hi_p, hi_n = float((pos >= pp3).mean()), float((neg >= pp3).mean())
    lo_p, lo_n = float((pos <= bp4).mean()), float((neg <= bp4).mean())
    mid_p = 1.0 - hi_p - lo_p
    mid_n = 1.0 - hi_n - lo_n
    # A zero denominator is bounded away from zero rather than reported as an
    # infinite ratio, the convention phase5 already uses. A zero numerator is left
    # at zero, and an EMPTY band is not given a ratio at all: flooring the numerator
    # would report a likelihood ratio -- and an evidence tier -- computed from the
    # class sizes alone, with no variant from the band in it.
    n_hi = int((s >= pp3).sum())
    n_lo = int((s <= bp4).sum())
    lr_pp3 = (hi_p / max(hi_n, 1.0 / (n_neg + 1))) if n_hi >= MIN_BAND else np.nan
    lr_bp4 = (lo_p / max(lo_n, 1.0 / (n_neg + 1))) if n_lo >= MIN_BAND else np.nan
    spec = 1.0 - hi_n
    lr_minus = (1.0 - hi_p) / spec if spec > 0 else np.nan
    return {
        "n_pos": n_pos, "n_neg": n_neg,
        "n_in_pp3_band": n_hi, "n_in_bp4_band": n_lo,
        "frac_pp3": float((s >= pp3).mean()), "frac_grey": float(((s > bp4) & (s < pp3)).mean()),
        "frac_bp4": float((s <= bp4).mean()),
        "sens_pp3": hi_p, "spec_pp3": spec,
        "lr_pp3": lr_pp3, "lr_bp4": lr_bp4, "lr_minus_at_pp3": lr_minus,
        "prop_grey_pos": mid_p, "prop_grey_neg": mid_n,
    }


def _boot_ci(sub: pd.DataFrame, score: str, pp3: float, bp4: float,
             keys: tuple[str, ...], cell: tuple = ()) -> dict:
    """Gene-clustered bootstrap, 2,000 resamples of the genes, as phase5 does.

    Where a stratum sits inside one gene there is no cluster structure to resample
    and the interval is reported as not evaluable rather than as a variant-level
    interval that would understate the uncertainty.
    """
    rng = _rng_for(score, *cell)
    genes = sub["gene"].to_numpy()
    ug = pd.unique(genes)
    out = {f"{k}_lo": np.nan for k in keys} | {f"{k}_hi": np.nan for k in keys}
    out["n_genes"] = len(ug)
    if len(ug) < 2:
        out["ci_basis"] = "single gene -- no cluster structure"
        return out
    idx_by_gene = {g: np.where(genes == g)[0] for g in ug}
    y_all = sub["y_assay"].to_numpy(dtype=float)
    s_all = sub[score].to_numpy(dtype=float)
    draws = {k: [] for k in keys}
    for _ in range(N_BOOT):
        idx = np.concatenate([idx_by_gene[g] for g in rng.choice(ug, len(ug), replace=True)])
        st = _band_stats(y_all[idx], s_all[idx], pp3, bp4)
        for k in keys:
            draws[k].append(st.get(k, np.nan))
    for k in keys:
        v = np.asarray(draws[k], dtype=float)
        v = v[np.isfinite(v)]
        out[f"{k}_n_boot_used"] = int(v.size)
        if v.size >= N_BOOT * 0.5:
            lo, hi = np.nanpercentile(v, [2.5, 97.5])
            out[f"{k}_lo"], out[f"{k}_hi"] = float(lo), float(hi)
    # A draw is dropped when the resampled gene set fell below the minimum count on
    # one side, so the interval is a percentile of the draws that CLEARED the gate,
    # not of all N_BOOT. Where many are dropped the interval is conditioned on having
    # drawn the genes that carry the scarce side, which narrows it; the used count is
    # reported per quantity so that is visible rather than implied.
    used = min(int(out.get(f"{k}_n_boot_used", 0)) for k in keys)
    out["ci_basis"] = (f"gene-clustered bootstrap over {len(ug)} genes; "
                       f"{used} of {N_BOOT} draws evaluable")
    return out


def evaluate(df: pd.DataFrame, score: str, cfg: dict) -> pd.DataFrame:
    pp3 = cfg["thresholds"]["pp3"]["value"]
    bp4 = cfg["thresholds"]["bp4"]["value"]
    path_bands = cfg["evidence_bands"]["pathogenic_lr_thresholds"]
    ben_bands = cfg["evidence_bands"]["benign_lr_thresholds"]
    keys = ("lr_pp3", "lr_bp4", "sens_pp3", "spec_pp3", "lr_minus_at_pp3")

    rows = []
    for stratum in STRATUM_ORDER:
        ss = (df if stratum == "all_1_50"
              else df[df["stratum"] != "pm12"] if stratum == "s3_50"
              else df[df["stratum"] == stratum])
        for arm in ARM_ORDER:
            sa = ss if arm == "all" else ss[ss["clinvar_arm"] == arm]
            scopes = [("pooled", None)] + [("gene", g) for g in sorted(sa["gene"].unique())]
            for scope, gene in scopes:
                sub = sa if gene is None else sa[sa["gene"] == gene]
                sub = sub.dropna(subset=["y_assay", score])
                base = {
                    "score_column": score, "stratum": stratum, "clinvar_arm": arm,
                    "scope": scope, "gene": gene or "(pooled)",
                    "n": int(len(sub)),
                    "walker_in_scope": stratum in IN_SCOPE,
                }
                st = _band_stats(sub["y_assay"].to_numpy(dtype=float),
                                 sub[score].to_numpy(dtype=float), pp3, bp4)
                if not st:
                    rows.append(base | {
                        "n_pos": int((sub["y_assay"] == 1).sum()),
                        "n_neg": int((sub["y_assay"] == 0).sum()),
                        "status": "not evaluable",
                        "reason": f"fewer than {MIN_POS} labelled variants on one side",
                    })
                    continue
                ci = (_boot_ci(sub, score, pp3, bp4, keys, (stratum, arm))
                      if scope == "pooled"
                      else {"ci_basis": "per-gene point estimate; no interval"})
                row = base | st | ci | {"status": "ok", "reason": ""}
                row["tier_pp3"] = tier_pathogenic(st["lr_pp3"], path_bands)
                row["tier_bp4"] = tier_benign(st["lr_bp4"], ben_bands)
                for p in PRIORS:
                    row[f"posterior_pp3_prior{p:g}"] = posterior(st["lr_pp3"], p)
                    row[f"posterior_bp4_prior{p:g}"] = posterior(st["lr_bp4"], p)
                rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    cfg = load_config()
    df = pd.read_parquet(SET_PATH)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    cols = ["spliceai_walker", "spliceai"]
    cols = [c for c in cols if c in df.columns]
    if "spliceai_walker" not in cols:
        print("[E2] WARNING: spliceai_walker not in the set; "
              "running on the distance-50 column only")
    out = pd.concat([evaluate(df, c, cfg) for c in cols], ignore_index=True)
    out.to_csv(OUT, index=False)

    print(f"[E2] cut points PP3 >= {cfg['thresholds']['pp3']['value']} "
          f"({cfg['thresholds']['pp3']['evidence_strength']}), "
          f"BP4 <= {cfg['thresholds']['bp4']['value']} "
          f"({cfg['thresholds']['bp4']['evidence_strength']}); "
          f"grey zone {cfg['thresholds']['grey_zone']['interval']}")
    print(f"[E2] wrote {OUT} ({len(out)} rows)\n")
    show = out[(out.scope == "pooled") & (out.status == "ok")]
    for c in cols:
        print(f"--- {c} (pooled) ---")
        print(show[show.score_column == c][
            ["stratum", "clinvar_arm", "n", "n_pos", "n_neg", "frac_pp3", "frac_grey",
             "frac_bp4", "sens_pp3", "spec_pp3", "lr_pp3", "lr_pp3_lo", "lr_pp3_hi",
             "tier_pp3", "lr_bp4", "tier_bp4"]].round(4).to_string(index=False))
        print()
    ne = out[out.status == "not evaluable"]
    print(f"[E2] not evaluable: {len(ne)} cells")
    if len(ne):
        print(ne.groupby(["score_column", "stratum", "clinvar_arm", "scope"], observed=True)
                .size().reset_index(name="n").to_string(index=False))


if __name__ == "__main__":
    main()
