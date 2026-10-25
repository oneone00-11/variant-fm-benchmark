"""E19 C4 -- prior sensitivity: E3's in-sample thresholds at the priors 0.05, 0.10 and
0.185.

The priors are fixed values (E19.6): 0.10 is the frozen fit's; 0.05 and 0.185 are round
values of the damaging fractions of the distal and proximal bands, taken by the
literature check (its question 6) and fixed in advance, not recomputed from the present
fractions.

The tier edges move with the prior so that each keeps the posterior its frozen edge
gives at 0.10 (fixed posterior; E19.6):

    LR_p = LR_0.10 x [0.10 / (1 - 0.10)] / [p / (1 - p)]

All four tiers on both sides, Very strong included, are converted and refitted (the
author's decision of 2026-10-25), in one-to-one correspondence with the frozen fit:
every column E3 fits (the panel and the fusion, its out-of-fold column computed as E3
computes it), every stratum of evidence_thresholds.csv, with E3's own function
(evid_interval_lr.thresholds_for), bootstrap and seed key (tool, stratum, "full"). The
edges change the thresholds only, never the local likelihood ratios or their bounds.
The run at 0.10 must give evidence_thresholds.csv's in-sample thresholds back exactly,
or the stage stops.

Outputs, in reports/evidence/suppl_analyses/, each with the formula and the three
priors' edges in its header (lines beginning '#'):
  c4_band_edges.csv          per prior x side x tier: the prior's odds, the edge at 0.10
                             and at the prior, and the posterior both keep
  c4_prior_sensitivity.csv   per prior x column x stratum x side x tier: the edge, the
                             counts, the threshold and whether the tier is reached,
                             beside the same cell at 0.10; tier_change ('gained',
                             'lost', 'unchanged') and threshold_changed against 0.10; and
                             for the column, stratum and side, whether Strong is reached
                             at 0.10 and at the prior (strong_still_reached: blank where
                             Strong is not reached at 0.10)

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_c4
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import evid_common as K
from . import evid_interval_lr as E3

OUT_DIR = Path("reports/evidence/suppl_analyses")
OUT_EDGES = OUT_DIR / "c4_band_edges.csv"
OUT_SENS = OUT_DIR / "c4_prior_sensitivity.csv"
FROZEN = E3.REPORT_DIR / "evidence_thresholds.csv"
BASE_PRIOR = 0.10
PRIORS = [0.05, BASE_PRIOR, 0.185]
SIDES = [("pp3", "path"), ("bp4", "ben")]
FORMULA = "LR_p = LR_0.10 x [0.10 / (1 - 0.10)] / [p / (1 - p)]"


def _odds(p: float) -> float:
    return p / (1 - p)


def edges_at(prior: float, cuts: dict) -> dict:
    """The edges at a prior that keep the posteriors the 0.10 edges give."""
    if prior == BASE_PRIOR:
        return {k: dict(v) for k, v in cuts.items()}
    f = _odds(BASE_PRIOR) / _odds(prior)
    return {k: {t: v * f for t, v in side.items()} for k, side in cuts.items()}


def header(cuts: dict) -> str:
    lines = ["E19 C4: prior sensitivity of E3's in-sample thresholds (src/evid_suppl_c4.py).",
             "Fixed posterior: at prior p each tier's edge keeps the posterior its edge at 0.10 gives:",
             f"  {FORMULA}",
             "Priors 0.05, 0.10 (the frozen fit) and 0.185, fixed in advance (E19.6).",
             "Edges (supporting, moderate, strong, very_strong):"]
    for p in PRIORS:
        e = edges_at(p, cuts)
        for side, key in SIDES:
            lines.append(f"  p = {p}: {side.upper()} " + ", ".join(repr(e[key][t]) for t in E3.TIERS))
    return "".join(f"# {x}\n" for x in lines)


def band_edges(cuts: dict) -> pd.DataFrame:
    rows = []
    for p in PRIORS:
        e = edges_at(p, cuts)
        for side, key in SIDES:
            for t in E3.TIERS:
                o = _odds(BASE_PRIOR) * cuts[key][t]
                rows.append({"prior": p, "prior_odds": _odds(p), "side": side, "tier": t,
                             "lr_at_0_10": cuts[key][t], "lr": e[key][t],
                             "posterior": o / (1 + o)})
    return pd.DataFrame(rows)


def fit_all(df: pd.DataFrame, tools: list[str], strata: list[str], cuts: dict) -> pd.DataFrame:
    """E3's in-sample fit for every column x stratum at the given edges, rows taken as
    E3 takes them (run_tool_stratum)."""
    rows = []
    for tool in tools:
        for stratum in strata:
            sub = K.stratum_frame(df, stratum).dropna(subset=["y_assay", tool])
            y = sub["y_assay"].to_numpy(dtype=float)
            s = sub[tool].to_numpy(dtype=float)
            genes = sub["gene"].astype(str).to_numpy()
            base = {"tool": tool, "stratum": stratum,
                    "n_pos": int((y == 1).sum()), "n_neg": int((y == 0).sum())}
            full, reason = None, ""
            if base["n_pos"] < K.MIN_POS or base["n_neg"] < K.MIN_NEG:
                reason = f"fewer than {K.MIN_POS} labelled variants on one side"
            else:
                full = E3.thresholds_for(y, s, genes, cuts["path"], cuts["ben"],
                                         n_boot=E3.N_BOOT_GENE, cluster=True,
                                         cell=(tool, stratum, "full"))
                if full is None:
                    reason = "score takes too few distinct values"
            for side, key in SIDES:
                for t in E3.TIERS:
                    thr = full[key][t] if full is not None else np.nan
                    rows.append(base | {"side": side, "tier": t, "cut_lr": cuts[key][t],
                                        "status": "ok" if full is not None else "not evaluable",
                                        "reason": reason, "threshold": thr,
                                        "reached": bool(np.isfinite(thr))})
    return pd.DataFrame(rows)


def _same(a: float, b: float) -> bool:
    return (np.isnan(a) and np.isnan(b)) or a == b


def check_reproduces_frozen(fit: pd.DataFrame) -> None:
    e = K.read_back(FROZEN)
    bad = []
    for r in fit.itertuples():
        f = e[(e.tool == r.tool) & (e.stratum == r.stratum) & (e.tier == r.tier)]
        if len(f) != 1:
            bad.append((r.tool, r.stratum, r.tier, "no frozen row"))
            continue
        f = f.iloc[0]
        if not (_same(r.threshold, f[f"{r.side}_threshold_insample"])
                and r.reached == bool(f[f"{r.side}_threshold_reachable"])
                and (r.n_pos, r.n_neg) == (f.n_pos, f.n_neg) and r.status == f.status):
            bad.append((r.tool, r.stratum, r.side, r.tier))
    if bad or len(fit) != 2 * len(e):
        raise SystemExit(f"[C4] the fit at prior 0.10 does not give {FROZEN} back "
                         f"({len(bad)} cell(s), e.g. {bad[:3]}): stopped")


def main() -> None:
    cfg = yaml.safe_load(E3.CONFIG_PATH.read_text())
    path_cuts, ben_cuts = K.acmg_bands(cfg)
    cuts = {"path": path_cuts, "ben": ben_cuts}
    df = K.load_set()
    df[K.FUSION] = K.logo_fusion(df, [t for t in K.FUSION_FEATURES if t in df.columns])
    tools = K.panel_of(df) + [K.FUSION]
    strata = list(dict.fromkeys(K.read_back(FROZEN)["stratum"]))

    fits = {p: fit_all(df, tools, strata, edges_at(p, cuts)) for p in PRIORS}
    check_reproduces_frozen(fits[BASE_PRIOR])
    print(f"[C4] the fit at prior {BASE_PRIOR} gives {FROZEN} back exactly")

    keys = ["tool", "stratum", "side", "tier"]
    base = fits[BASE_PRIOR][keys + ["cut_lr", "threshold", "reached"]].rename(columns={
        "cut_lr": "cut_lr_at_0_10", "threshold": "threshold_at_0_10",
        "reached": "reached_at_0_10"})
    strong = fits[BASE_PRIOR][fits[BASE_PRIOR].tier == "strong"][
        ["tool", "stratum", "side", "reached"]].rename(columns={"reached": "strong_reached_at_0_10"})
    parts = []
    for p in PRIORS:
        m = fits[p].merge(base, on=keys, validate="one_to_one")
        m["tier_change"] = np.select(
            [m.reached & ~m.reached_at_0_10, ~m.reached & m.reached_at_0_10],
            ["gained", "lost"], "unchanged")
        m["threshold_changed"] = [not _same(a, b) for a, b in zip(m.threshold, m.threshold_at_0_10)]
        st = fits[p][fits[p].tier == "strong"][["tool", "stratum", "side", "reached"]].rename(
            columns={"reached": "strong_reached"})
        m = m.merge(strong, on=["tool", "stratum", "side"]).merge(st, on=["tool", "stratum", "side"])
        m["strong_still_reached"] = m.strong_reached.where(m.strong_reached_at_0_10, other=pd.NA)
        m.insert(0, "prior", p)
        parts.append(m)
    out = pd.concat(parts, ignore_index=True)[
        ["prior", "tool", "stratum", "side", "tier", "cut_lr", "cut_lr_at_0_10", "n_pos", "n_neg",
         "status", "reason", "threshold", "reached", "threshold_at_0_10", "reached_at_0_10",
         "tier_change", "threshold_changed", "strong_reached_at_0_10", "strong_reached",
         "strong_still_reached"]]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    head = header(cuts)
    with open(OUT_EDGES, "w", newline="") as f:
        f.write(head)
        band_edges(cuts).to_csv(f, index=False)
    with open(OUT_SENS, "w", newline="") as f:
        f.write(head)
        out.to_csv(f, index=False)

    for p in PRIORS:
        o = out[out.prior == p]
        lost = o[(o.tier == "strong") & o.strong_reached_at_0_10 & ~o.strong_reached]
        print(f"[C4] prior {p}: {int((o.tier_change != 'unchanged').sum())} of {len(o)} cells "
              f"change tier; Strong no longer reached in {len(lost)} column-stratum-side cell(s)"
              + "".join(f"\n      {r.tool} {r.stratum} {r.side}" for r in lost.itertuples()))


if __name__ == "__main__":
    main()
