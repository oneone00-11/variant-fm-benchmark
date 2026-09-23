"""E3 -- score-to-evidence intervals, per tool and per territory (Pejaver 2022).

The published calibration study asked whether fusing predictors improves a proper
score. This stage asks the question a laboratory actually has: at what score does
this tool start to be worth Supporting, Moderate or Strong evidence in this
territory, and does that threshold hold on a gene it was not chosen on.

Method, following Pejaver et al. 2022 rather than inventing one:

  * local likelihood ratio. For each observed score s, a window [s-eps, s+eps] is
    taken around it, with eps the smallest value that puts at least
    MIN_WINDOW labelled variants inside. lr(s) is the fraction of damaging
    variants falling in the window divided by the fraction of normal variants
    falling in it. This is the density ratio, which is what Pejaver's local
    posterior reduces to once the prior is divided out, so it needs no prior and
    is directly comparable across strata with different positive rates.
  * thresholds. Pejaver's Equations 8 and 9: the threshold for a tier is the
    smallest score at or above which EVERY local likelihood ratio's one-sided 95%
    bound clears the tier's cut (and symmetrically, the largest score at or below
    which every bound stays under the reciprocal cut). The "for every score above"
    requirement is what makes a threshold a threshold rather than a lucky bin.
  * the bound is taken in the stringent direction, as Pejaver does: the 5th
    percentile of the bootstrap distribution on the pathogenic side, the 95th on
    the benign side. Pejaver writes this as lr - Delta with Delta from 10,000
    bootstrap iterations; the percentile is the same bound without assuming the
    interval is symmetric.

Two deliberate departures, both recorded here because they change numbers:

  1. Pejaver resamples variants. The primary interval here resamples GENES, because
     these seven assays are the cluster structure of this dataset and a
     variant-level interval on clustered data is anti-conservative. The
     variant-level bootstrap is run as well, at Pejaver's 10,000 iterations, and
     reported beside it.
  2. The tier cuts are the Tavtigian point system at a prior of 0.10 (2.08 / 4.33 /
     18.7 / 350), which is what both companion manuscripts and Walker et al. use.
     Pejaver's own cuts (2.406 / 5.79 / 33.53 / 1124) are derived at their prior of
     0.0441 and are not used.

The fusion's LOGO columns carry a leakage the single tools' do not, and it is
flagged in the output rather than left for a reader to notice. `logo_fusion`
produces ONE out-of-fold column: gene g's score comes from a model trained on the
other six. When E3 then holds out gene g and fits a threshold on the six training
genes' scores, each of those six scores came from a model whose training set
INCLUDED g -- so g's labels have already influenced the scores the threshold is
fitted on. The threshold is therefore not free of the held-out gene, so its
held-out likelihood ratio is not an independent estimate and is read as optimistic.
Where one could be computed it nearly always cleared the tier's boundary; that is
not evidence that a fusion threshold transfers, and an earlier version of this
docstring, which said the fusion's thresholds did not carry across genes, described
the evaluable folds wrongly. A clean estimate would need a nested design -- refit
the fusion inside each outer fold -- which is not run here.
Rows for the fusion carry `logo_threshold_leakage`.

A fitted fusion threshold and a single tool's threshold are not the same kind of
object, and the LOGO columns show it. A single tool emits the same score scale in
every gene, so a threshold chosen on six genes means the same thing on the seventh.
The elastic-net fusion is refitted in every fold, so its out-of-fold score is on a
per-fold scale and a threshold carried across folds does not land in the same place.
`evidence_thresholds_logo_folds.csv` holds the per-fold thresholds and held-out
likelihood ratios so that instability is visible rather than averaged away.

The brief describes the BP4 side in terms of LR-. It is computed here as the local
likelihood ratio on the low-score side, which is Pejaver's own quantity and
Walker's; the classical LR- at a threshold is a different number and is reported by
E2, not here.

Run (PYTHONPATH=phase1):  python -m src.evid_interval_lr [--quick]
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import config as C
from . import evid_common as K

CONFIG_PATH = Path("config/walker2023.yaml")
REPORT_DIR = Path("reports/evidence")

MIN_WINDOW = 100          # Pejaver: at least 100 pathogenic and benign combined
N_BOOT_GENE = 2000        # project convention for gene-clustered intervals
N_BOOT_VARIANT = 10000    # Pejaver's iteration count
TIERS = ["supporting", "moderate", "strong", "very_strong"]

def _rng_for(*key) -> np.random.Generator:
    """Seeded from the cell's identity rather than from call order, so a single
    tool x stratum x fold is reproducible without re-running the whole table."""
    h = hashlib.sha256(("|".join(str(k) for k in key)).encode()).digest()
    return np.random.default_rng([C.RANDOM_SEED, int.from_bytes(h[:8], "big")])


# ---------------------------------------------------------------------------
# local likelihood ratio
# ---------------------------------------------------------------------------
# The window is [s-eps, s+eps] closed. Rebuilding those endpoints by subtracting a
# stored half-width does not round-trip exactly in binary floating point, and since
# the interval is closed that silently drops the whole tie group sitting on the
# boundary -- measured at 102 grid points across the delivered curves, the worst
# holding 88 observations instead of 100. The endpoints are therefore widened by a
# relative tolerance rather than trusted to reconstruct.
_BOUND_TOL = 1e-12


def _window_bounds(grid: np.ndarray, eps: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Closed window endpoints, widened by a relative tolerance (see _BOUND_TOL)."""
    lo, hi = grid - eps, grid + eps
    return (lo - _BOUND_TOL * np.maximum(np.abs(lo), 1.0),
            hi + _BOUND_TOL * np.maximum(np.abs(hi), 1.0))


def _eps_for_grid(scores_sorted: np.ndarray, grid: np.ndarray,
                  min_window: int) -> np.ndarray:
    """Smallest eps putting >= min_window observations inside [s-eps, s+eps].

    At the ends of the range the symmetric window runs off the data, so the window
    is clipped and holds >= min_window observations on the side that exists. Pejaver
    scales the requirement down at the ends instead; the count is kept fixed here
    and the clipping is what gives way, which is the more conservative of the two.
    """
    n = len(scores_sorted)
    eps = np.empty(len(grid))
    for j, s in enumerate(grid):
        i = int(np.searchsorted(scores_sorted, s))
        lo, hi = i - 1, i
        d = 0.0
        cnt = 0
        while cnt < min_window and (lo >= 0 or hi < n):
            dl = s - scores_sorted[lo] if lo >= 0 else np.inf
            dr = scores_sorted[hi] - s if hi < n else np.inf
            if dl <= dr:
                d = dl
                lo -= 1
            else:
                d = dr
                hi += 1
            cnt += 1
        eps[j] = d if np.isfinite(d) else (scores_sorted[-1] - scores_sorted[0])
    return eps


def _local_lr(pos_sorted: np.ndarray, neg_sorted: np.ndarray,
              lo: np.ndarray, hi: np.ndarray, *, counts: bool = False):
    """Vectorised density ratio over fixed window endpoints.

    A window holding no damaging variant has a local likelihood ratio of ZERO, and
    that is what is returned. An earlier version floored the numerator at
    1/(npos+1) -- the same defect that was removed from band_lr: it invents a
    positive ratio out of a window with no observation of the class in the
    numerator. On the benign side that floor was load-bearing and wrong. It sat
    ABOVE the Strong and Very strong benign cuts, so those tiers could never be
    reached, and "no BP4 evidence in this territory" was an artefact of the
    flooring rather than a property of the data. What a zero-count window does need
    is an upper bound; that is supplied by _rule_of_three_upper, not by moving the
    estimate.
    """
    npos, nneg = len(pos_sorted), len(neg_sorted)
    if npos == 0 or nneg == 0:
        nan = np.full(len(lo), np.nan)
        return (nan, np.zeros(len(lo), int), np.zeros(len(lo), int)) if counts else nan
    cp = (np.searchsorted(pos_sorted, hi, side="right")
          - np.searchsorted(pos_sorted, lo, side="left"))
    cn = (np.searchsorted(neg_sorted, hi, side="right")
          - np.searchsorted(neg_sorted, lo, side="left"))
    fn = np.maximum(cn / nneg, 1.0 / (nneg + 1))
    lr = (cp / npos) / fn
    return (lr, cp, cn) if counts else lr


def _rule_of_three_upper(cp: np.ndarray, cn: np.ndarray,
                         npos: int, nneg: int) -> np.ndarray:
    """One-sided 95% upper bound on the local ratio where a window holds no
    damaging variant.

    The bootstrap cannot bound such a window: resampling the same variants can never
    put a damaging variant into a score range that contains none, so every resample
    returns the same zero and the interval collapses onto it -- which would declare
    Very strong benign evidence from no observation at all. The rule of three
    (Hanley and Lippman-Hand 1983) is the standard one-sided 95% upper bound on a
    rate observed as zero in n trials, 3/n; dividing by the window's normal fraction
    gives the matching bound on the ratio. Applied only where cp == 0; everywhere
    else the bootstrap bound stands.
    """
    fn = np.maximum(cn / nneg, 1.0 / (nneg + 1))
    return np.where(cp == 0, (3.0 / npos) / fn, -np.inf)


def local_lr_curve(y: np.ndarray, s: np.ndarray, min_window: int = MIN_WINDOW):
    """(grid, eps, lo, hi, lr, cp, cn) on the observed score values."""
    m = np.isfinite(s) & np.isfinite(y)
    y, s = y[m], s[m]
    grid = np.unique(s)
    all_sorted = np.sort(s)
    eps = _eps_for_grid(all_sorted, grid, min_window)
    lo, hi = _window_bounds(grid, eps)
    lr, cp, cn = _local_lr(np.sort(s[y == 1]), np.sort(s[y == 0]), lo, hi, counts=True)
    return grid, eps, lo, hi, lr, cp, cn


def bootstrap_bounds(y: np.ndarray, s: np.ndarray, genes: np.ndarray,
                     wlo: np.ndarray, whi: np.ndarray, *, cluster: bool,
                     n_boot: int, cell: tuple = ()) -> tuple[np.ndarray, np.ndarray]:
    """One-sided 5th/95th percentile bands of the local likelihood ratio.

    The grid and the window endpoints are held at their full-data values, so every
    resample is evaluated at the same places; only the sample changes.
    """
    m = np.isfinite(s) & np.isfinite(y)
    y, s, genes = y[m], s[m], genes[m]
    rng = _rng_for("cluster" if cluster else "variant", n_boot, *cell)
    draws = np.empty((n_boot, len(wlo)))
    if cluster:
        ug = pd.unique(genes)
        if len(ug) < 2:
            return (np.full(len(wlo), np.nan),) * 2
        idx_by_gene = {g: np.where(genes == g)[0] for g in ug}
        for b in range(n_boot):
            idx = np.concatenate([idx_by_gene[g]
                                  for g in rng.choice(ug, len(ug), replace=True)])
            yb, sb = y[idx], s[idx]
            draws[b] = _local_lr(np.sort(sb[yb == 1]), np.sort(sb[yb == 0]), wlo, whi)
    else:
        n = len(y)
        for b in range(n_boot):
            idx = rng.integers(0, n, n)
            yb, sb = y[idx], s[idx]
            draws[b] = _local_lr(np.sort(sb[yb == 1]), np.sort(sb[yb == 0]), wlo, whi)
    with np.errstate(invalid="ignore"):
        lo = np.nanpercentile(draws, 5, axis=0)
        hi = np.nanpercentile(draws, 95, axis=0)
    return lo, hi


# ---------------------------------------------------------------------------
# thresholds (Pejaver Eq. 8 / Eq. 9)
# ---------------------------------------------------------------------------
def pathogenic_threshold(grid: np.ndarray, lower: np.ndarray, cut: float) -> float:
    """Smallest tau such that every score >= tau has lower bound >= cut."""
    ok = np.isfinite(lower) & (lower >= cut)
    tau = np.nan
    for j in range(len(grid) - 1, -1, -1):
        if not ok[j]:
            break
        tau = grid[j]
    return tau


def benign_threshold(grid: np.ndarray, upper: np.ndarray, cut: float) -> float:
    """Largest tau such that every score <= tau has upper bound <= cut."""
    ok = np.isfinite(upper) & (upper <= cut)
    tau = np.nan
    for j in range(len(grid)):
        if not ok[j]:
            break
        tau = grid[j]
    return tau


def thresholds_for(y, s, genes, path_cuts, ben_cuts, *, n_boot, cluster=True,
                   min_window=MIN_WINDOW, cell=()):
    grid, eps, wlo, whi, lr, cp, cn = local_lr_curve(y, s, min_window)
    if len(grid) < 2 or not np.isfinite(lr).any():
        return None
    boot_lo, boot_hi = bootstrap_bounds(y, s, genes, wlo, whi, cluster=cluster,
                                        n_boot=n_boot, cell=cell)
    ok = np.isfinite(s) & np.isfinite(y)
    npos, nneg = int((y[ok] == 1).sum()), int((y[ok] == 0).sum())
    # where a window holds no damaging variant the bootstrap carries no information;
    # the rule-of-three bound takes over there, on the benign side only
    rot = _rule_of_three_upper(cp, cn, npos, nneg)
    hi_eff = np.maximum(boot_hi, rot)
    out = {"grid": grid, "eps": eps, "lr": lr, "lo": boot_lo, "hi": hi_eff,
           "hi_boot": boot_hi, "cp": cp, "cn": cn, "path": {}, "ben": {}}
    for tier in TIERS:
        out["path"][tier] = pathogenic_threshold(grid, boot_lo, path_cuts[tier])
        out["ben"][tier] = benign_threshold(grid, hi_eff, ben_cuts[tier])
    return out


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def run_tool_stratum(df, tool, stratum, path_cuts, ben_cuts, walker, n_boot_variant):
    sub = K.stratum_frame(df, stratum).dropna(subset=["y_assay", tool])
    y = sub["y_assay"].to_numpy(dtype=float)
    s = sub[tool].to_numpy(dtype=float)
    genes = sub["gene"].astype(str).to_numpy()
    n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
    base = {"tool": tool, "stratum": stratum, "n": int(len(sub)),
            "n_pos": n_pos, "n_neg": n_neg, "n_genes": int(len(pd.unique(genes)))}
    if n_pos < K.MIN_POS or n_neg < K.MIN_NEG:
        return [base | {"tier": t, "status": "not evaluable",
                        "reason": f"fewer than {K.MIN_POS} labelled variants on one side"}
                for t in TIERS], None, None

    full = thresholds_for(y, s, genes, path_cuts, ben_cuts,
                          n_boot=N_BOOT_GENE, cluster=True,
                          cell=(tool, stratum, "full"))
    if full is None:
        return [base | {"tier": t, "status": "not evaluable",
                        "reason": "score takes too few distinct values"} for t in TIERS], None, None
    var = thresholds_for(y, s, genes, path_cuts, ben_cuts,
                         n_boot=n_boot_variant, cluster=False,
                         cell=(tool, stratum, "full"))

    # LOGO: the threshold is chosen on six genes and reported on the seventh
    logo = {t: [] for t in TIERS}
    logo_ben = {t: [] for t in TIERS}
    for g in pd.unique(genes):
        tr, te = genes != g, genes == g
        if (y[tr] == 1).sum() < K.MIN_POS or (y[tr] == 0).sum() < K.MIN_NEG:
            continue
        fit = thresholds_for(y[tr], s[tr], genes[tr], path_cuts, ben_cuts,
                             n_boot=N_BOOT_GENE, cluster=True,
                             cell=(tool, stratum, "logo", g))
        if fit is None:
            continue
        for t in TIERS:
            logo[t].append({"gene": g, "tau": fit["path"][t],
                            "heldout_lr": K.band_lr(y[te], s[te], fit["path"][t], "upper")})
            logo_ben[t].append({"gene": g, "tau": fit["ben"][t],
                                "heldout_lr": K.band_lr(y[te], s[te], fit["ben"][t], "lower")})

    logo_detail = []
    for t in TIERS:
        for d in logo[t]:
            logo_detail.append({"tool": tool, "stratum": stratum, "tier": t,
                                "side": "pp3", "heldout_gene": d["gene"],
                                "threshold_from_6_genes": d["tau"],
                                "heldout_lr": d["heldout_lr"]})

    rows = []
    for t in TIERS:
        lo_vals = [d["heldout_lr"] for d in logo[t] if np.isfinite(d["heldout_lr"])]
        lo_taus = [d["tau"] for d in logo[t] if np.isfinite(d["tau"])]
        bn_vals = [d["heldout_lr"] for d in logo_ben[t] if np.isfinite(d["heldout_lr"])]
        rows.append(base | {
            "tier": t,
            "status": "ok",
            "reason": "",
            # see the module docstring: the fusion is refitted in every fold, so the
            # six training genes' out-of-fold scores were produced by models that saw
            # the held-out gene. Its LOGO numbers are optimistic by construction.
            "logo_threshold_leakage": (
                "fusion out-of-fold scores in the training folds were produced by "
                "models trained on the held-out gene; held-out LR is optimistic"
                if tool == K.FUSION else ""),
            "path_cut_lr": path_cuts[t],
            "pp3_threshold_insample": full["path"][t],
            "pp3_threshold_reachable": bool(np.isfinite(full["path"][t])),
            "pp3_threshold_variant_boot": var["path"][t] if var else np.nan,
            "pp3_logo_folds_reaching": len(lo_taus),
            "pp3_logo_folds_with_heldout_lr": len(lo_vals),
            "pp3_logo_threshold_spread": (float(np.max(lo_taus) - np.min(lo_taus))
                                          if lo_taus else np.nan),
            "pp3_logo_threshold_median": float(np.median(lo_taus)) if lo_taus else np.nan,
            "pp3_logo_heldout_lr_median": float(np.median(lo_vals)) if lo_vals else np.nan,
            "pp3_logo_heldout_lr_min": float(np.min(lo_vals)) if lo_vals else np.nan,
            "pp3_logo_heldout_lr_max": float(np.max(lo_vals)) if lo_vals else np.nan,
            "ben_cut_lr": ben_cuts[t],
            "bp4_threshold_insample": full["ben"][t],
            "bp4_threshold_reachable": bool(np.isfinite(full["ben"][t])),
            "bp4_threshold_variant_boot": var["ben"][t] if var else np.nan,
            "bp4_logo_heldout_lr_median": float(np.median(bn_vals)) if bn_vals else np.nan,
        })

    curve = pd.DataFrame({"tool": tool, "stratum": stratum, "score": full["grid"],
                          "eps": full["eps"], "n_damaging_in_window": full["cp"],
                          "n_normal_in_window": full["cn"], "local_lr": full["lr"],
                          "lr_lo_gene_boot": full["lo"],
                          "lr_hi_gene_boot": full["hi"],
                          "lr_hi_bootstrap_only": full["hi_boot"]})
    if var is not None:
        curve["lr_lo_variant_boot"] = var["lo"]
        curve["lr_hi_variant_boot"] = var["hi"]

    # where Walker's fixed cut points land on this curve
    if tool in ("spliceai", "spliceai_walker"):
        for name, v in (("pp3", walker["thresholds"]["pp3"]["value"]),
                        ("bp4", walker["thresholds"]["bp4"]["value"])):
            j = int(np.argmin(np.abs(full["grid"] - v)))
            for r in rows:
                r[f"walker_{name}_local_lr"] = float(full["lr"][j])
                # Both tails, and the STRINGENT one named. For PP3 that is the 5th
                # percentile; for BP4 it is the 95th. Reporting the 5th on both
                # sides -- which an earlier version did -- makes Walker's 0.1 cut
                # point look like Supporting benign evidence when the bound that
                # actually applies is nowhere near it.
                r[f"walker_{name}_local_lr_lo"] = float(full["lo"][j])
                r[f"walker_{name}_local_lr_hi"] = float(full["hi"][j])
                r[f"walker_{name}_stringent_bound"] = float(
                    full["lo"][j] if name == "pp3" else full["hi"][j])
                r[f"walker_{name}_nearest_grid_score"] = float(full["grid"][j])
    return rows, curve, pd.DataFrame(logo_detail)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick", action="store_true",
                    help="cut the variant-level bootstrap to 1,000 for a smoke run")
    ap.add_argument("--tools", default=None, help="comma-separated subset")
    args = ap.parse_args()

    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    path_cuts, ben_cuts = K.acmg_bands(cfg)
    n_boot_variant = 1000 if args.quick else N_BOOT_VARIANT

    df = K.load_set()
    tools = K.panel_of(df)
    df[K.FUSION] = K.logo_fusion(df, [t for t in K.FUSION_FEATURES if t in df.columns])
    tools = tools + [K.FUSION]
    if args.tools:
        tools = [t for t in args.tools.split(",") if t in df.columns or t == K.FUSION]

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    all_rows, curves, details = [], [], []
    for tool in tools:
        for stratum in K.STRATA:
            rows, curve, detail = run_tool_stratum(df, tool, stratum, path_cuts,
                                                   ben_cuts, cfg, n_boot_variant)
            all_rows += rows
            if detail is not None and len(detail):
                details.append(detail)
            if curve is not None:
                curves.append(curve)
                curve.to_csv(REPORT_DIR / f"interval_lr_{tool}_{stratum}.csv", index=False)
            print(f"[E3] {tool:16s} {stratum:9s} "
                  f"n={rows[0]['n']:>5} {rows[0].get('status','')}", flush=True)

    out = pd.DataFrame(all_rows)
    out.to_csv(REPORT_DIR / "evidence_thresholds.csv", index=False)
    if details:
        pd.concat(details, ignore_index=True).to_csv(
            REPORT_DIR / "evidence_thresholds_logo_folds.csv", index=False)
    print(f"\n[E3] wrote evidence_thresholds.csv ({len(out)} rows) and "
          f"{len(curves)} curve files")

    ok = out[out.status == "ok"]
    print("\n[E3] highest tier reached in-sample (PP3 side):")
    best = (ok[ok.pp3_threshold_reachable]
            .assign(rank=lambda d: d.tier.map({t: i for i, t in enumerate(TIERS)}))
            .sort_values("rank").groupby(["tool", "stratum"], observed=True).tail(1))
    print(best.pivot_table(index="tool", columns="stratum", values="tier",
                           aggfunc="first").to_string())


if __name__ == "__main__":
    main()
