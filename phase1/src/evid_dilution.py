"""E11 -- does a fitted threshold depend on how many genes it is fitted on, and which?

E3 fits each tool's pathogenic threshold on six genes and reports it on the seventh.
That is one training-set size and one composition per held-out gene. The natural
follow-up is whether the threshold, and whether it carries to a gene it was not
chosen on, would change if fewer genes, or different ones, had been available to fit
it on. This stage answers that by refitting on every subset of the genes rather than
arguing from the seven folds, and by testing every refitted threshold on every gene
the subset left out.

Design, and why each piece is the way it is:

  * tools. spliceai_walker, pangolin, alphagenome and avi: the splice-aware columns,
    with SpliceAI on Walker's distance-4999 basis rather than the distance-50 panel
    column. The fusion is left out on purpose. Its out-of-fold score is refitted in
    every fold (see E3's docstring), so a subset design would need the fusion itself
    refitted inside every subset, which is a different experiment. The four columns
    here are fixed per variant, so a threshold fitted on some genes and read on
    another carries no leakage of the held-out gene's labels.
  * strata. s3_10, s11_50 and s3_50, the in-scope pools; pm12 goes through the PVS1
    decision tree, not PP3. Within each stratum the genes are those with at least
    one scored, labelled variant, read off the data rather than fixed here, so the
    number of genes G can differ between strata.
  * training subsets. Every subset of size k = 2, ..., G-1. k = 1 is left out
    because the gene-clustered bootstrap needs two clusters to produce an interval
    at all (E3's bootstrap_bounds returns none with one), and k = G leaves no gene to
    hold out. The k = G fit is still made, once per tool and stratum, and carried in
    the summary as `threshold_all_genes`: it is the reference every subset's
    threshold is read against.
  * the fit. E3's own thresholds_for, with E3's MIN_WINDOW, E3's gene-cluster
    bootstrap count and the tier cuts from config/walker2023.yaml, on the training
    genes' variants only. A subset with fewer than MIN_POS damaging or MIN_NEG
    normal variants is not fitted, as in E3's folds. Only the pathogenic side is
    recorded, and only Moderate and Strong; the fit returns every tier at no extra
    cost, so adding one is a change to TIERS and nothing else.
  * the evaluation. Each fitted threshold is applied to every gene outside its
    training subset, one gene at a time, and the likelihood ratio of the band at
    and above it is measured there with evid_common.band_lr: the same point
    estimate E3 reports per fold, under the same rule that a band with fewer than
    MIN_BAND occupants, or a held-out gene with fewer than MIN_POS damaging or
    MIN_NEG normal variants, is not evaluated. The band's occupant, damaging and
    normal counts are recorded beside the ratio so the reason a pair is not
    evaluable can be read off the row.

Seeds. Every fit's bootstrap is seeded from its cell's identity through E3's
_rng_for, so a single fit can be reproduced on its own and the tables do not depend
on the order the fits run in or on how many worker processes run them (--jobs only
changes the wall-clock time). The subsets of size G-1 are exactly E3's
leave-one-gene-out training sets, and they are given E3's own cell key, (tool,
stratum, "logo", held-out gene), rather than a key built from the training genes.
That is deliberate: it makes the k = G-1 rows a recomputation of E3's per-fold table
rather than a second, differently seeded estimate of it, so they can be required to
agree with evidence_thresholds_logo_folds.csv exactly, and a disagreement means one
of the two stages has drifted from the other. The k = G reference reuses E3's
(tool, stratum, "full") key for the same reason and reproduces its in-sample
threshold. Every other subset is keyed (tool, stratum, "subset", then its training
genes in sorted order). The training rows are also taken in E3's order, the analysis
set's own order filtered to the subset, because the cluster bootstrap draws genes in
their order of first appearance and a reordered frame would draw different
resamples from the same seed.

One property of the method shapes the small-k rows and should be understood before
reading them. With k training genes a bootstrap resample is one gene drawn k times
with probability k^(1-k): one half at k = 2, one ninth at k = 3. At k = 2 each
single-gene resample holds a quarter of the draws, far more than the 5% tail the
pathogenic bound is read from, so the bound is at most the weaker gene's local ratio
and a score clears a tier's cut only where each training gene's own ratio, in the
same window, clears it too. That is the correct behaviour for an interval that
treats the gene as the unit of replication -- two genes are two observations -- but
it means a low reach rate at small k measures the interval's width as well as any
shift in the curve. Part of the spread of thresholds across subsets is also the
bootstrap's own Monte Carlo noise, since each subset has its own seed. Where the
tier is reached, the threshold's position across subsets and the held-out ratio it
delivers are the quantities that speak to composition.

The held-out genes are balanced across k. At size k each gene lies outside
C(G-1, k) subsets, the same number for every gene, so every gene is the held-out
gene equally often at every k and a trend across k is not a change in which genes
are being tested. Evaluability then removes pairs in two ways: a gene with too few
damaging or normal variants in a stratum is never evaluable, at any k, and a band
with too few occupants depends on where the threshold landed.

What the fractions are taken over. `frac_subsets_reached` is over all subsets of a
size, so a subset that could not be fitted counts as not reaching the tier;
`n_subsets_fitted` says how often that happened. `frac_pairs_cleared` and the
held-out ratio's median and quartiles are over the (subset, held-out gene) pairs
with an evaluable ratio only, and `n_pairs_evaluable` is their number.

Outputs, in reports/evidence/:

  threshold_dilution.csv            one row per tool x stratum x training subset x
                                    tier x held-out gene
  threshold_dilution_summary.csv    one row per tool x stratum x tier x k
  supplement/tableS8_dilution.csv   the summary, rounded for print

Run (PYTHONPATH=phase1):  python -m src.evid_dilution [--jobs N]
"""
from __future__ import annotations

import argparse
import itertools
import os
import time
from concurrent.futures import ProcessPoolExecutor
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import evid_common as K
from . import evid_interval_lr as E

REPORT_DIR = K.REPORT_DIR
OUT = REPORT_DIR / "threshold_dilution.csv"
OUT_SUMMARY = REPORT_DIR / "threshold_dilution_summary.csv"
OUT_SUPP = REPORT_DIR / "supplement" / "tableS8_dilution.csv"

TOOLS = ["spliceai_walker", "pangolin", "alphagenome", "avi"]
STRATA = ["s3_10", "s11_50", "s3_50"]
TIERS = ["moderate", "strong"]
MIN_K = 2                      # a gene-clustered interval needs two clusters
N_BOOT = E.N_BOOT_GENE         # E3's count for the gene-clustered interval
MIN_WINDOW = E.MIN_WINDOW      # E3's local-window size

# columns of the printed supplement table, and how each is rounded
_SUPP_COLS = ["tool", "stratum", "tier", "path_cut_lr", "n_genes", "k", "n_subsets",
              "n_subsets_reached", "frac_subsets_reached", "threshold_median",
              "threshold_q25", "threshold_q75", "threshold_all_genes",
              "n_pairs_evaluable", "n_pairs_cleared", "frac_pairs_cleared",
              "heldout_lr_median"]
# Thresholds are rounded to decimal places in score units, not to significant
# figures: alphagenome's thresholds sit close together near the top of its scale,
# and four significant figures merge its quartiles into one printed value.
_DP_THRESHOLD = 4
_SIG_LR = 3
_DP_FRACTION = 2


# ---------------------------------------------------------------------------
# the cells
# ---------------------------------------------------------------------------
def cell_arrays(df: pd.DataFrame, tool: str, stratum: str):
    """(y, s, genes) for one tool and stratum, in E3's row order."""
    sub = K.stratum_frame(df, stratum).dropna(subset=["y_assay", tool])
    return (sub["y_assay"].to_numpy(dtype=float), sub[tool].to_numpy(dtype=float),
            sub["gene"].astype(str).to_numpy())


def genes_of(genes: np.ndarray) -> list[str]:
    return sorted(pd.unique(genes).tolist())


def cell_key(tool: str, stratum: str, train: tuple, all_genes: list[str]) -> tuple:
    """The bootstrap's seed key: E3's own where the fit is one E3 also makes."""
    left_out = [g for g in all_genes if g not in train]
    if not left_out:
        return (tool, stratum, "full")
    if len(left_out) == 1:
        return (tool, stratum, "logo", left_out[0])
    return (tool, stratum, "subset", *sorted(train))


def fit_subset(y, s, genes, train, path_cuts, ben_cuts, cell) -> tuple[str, dict]:
    """Pathogenic thresholds fitted on the training genes' variants only.

    Returns (status, {tier: threshold}); a threshold is NaN where the tier is not
    reached, and the dict is empty where no fit was possible.
    """
    tr = np.isin(genes, list(train))
    yt, st, gt = y[tr], s[tr], genes[tr]
    if (yt == 1).sum() < K.MIN_POS or (yt == 0).sum() < K.MIN_NEG:
        return (f"not fitted: fewer than {K.MIN_POS} labelled variants on one side",
                {})
    fit = E.thresholds_for(yt, st, gt, path_cuts, ben_cuts, n_boot=N_BOOT,
                           cluster=True, min_window=MIN_WINDOW, cell=cell)
    if fit is None:
        return "not fitted: score takes too few distinct values", {}
    return "fitted", {t: float(fit["path"][t]) for t in TIERS}


def heldout_eval(y, s, genes, gene: str, tau: float) -> dict:
    """The band at and above tau, measured in one held-out gene."""
    te = genes == gene
    yh, sh = y[te], s[te]
    out = {"heldout_n_damaging": int((yh == 1).sum()),
           "heldout_n_normal": int((yh == 0).sum()),
           "heldout_n_in_band": None, "heldout_n_damaging_in_band": None,
           "heldout_n_normal_in_band": None, "heldout_lr": np.nan,
           "evaluable": False, "cleared_cut": None, "not_evaluable_reason": ""}
    if not np.isfinite(tau):
        out["not_evaluable_reason"] = "tier not reached on the training genes"
        return out
    band = sh >= tau
    out["heldout_n_in_band"] = int(band.sum())
    out["heldout_n_damaging_in_band"] = int((band & (yh == 1)).sum())
    out["heldout_n_normal_in_band"] = int((band & (yh == 0)).sum())
    lr = K.band_lr(yh, sh, tau, "upper")
    out["heldout_lr"] = lr
    if np.isfinite(lr):
        out["evaluable"] = True
    elif out["heldout_n_damaging"] < K.MIN_POS or out["heldout_n_normal"] < K.MIN_NEG:
        out["not_evaluable_reason"] = (f"held-out gene has fewer than {K.MIN_POS} "
                                       f"damaging or {K.MIN_NEG} normal variants")
    else:
        out["not_evaluable_reason"] = f"fewer than {K.MIN_BAND} variants in the band"
    return out


# ---------------------------------------------------------------------------
# the work list, run serially or across processes with identical results
# ---------------------------------------------------------------------------
_STATE: dict = {}


def _state():
    """The analysis set and the tier cuts, loaded once per process."""
    if not _STATE:
        cfg = yaml.safe_load(E.CONFIG_PATH.read_text())
        _STATE["cuts"] = K.acmg_bands(cfg)
        _STATE["df"] = K.load_set()
    return _STATE


def _run_task(task: tuple) -> tuple[str, dict]:
    tool, stratum, train, cell = task
    st = _state()
    key = ("arrays", tool, stratum)
    if key not in st:
        st[key] = cell_arrays(st["df"], tool, stratum)
    y, s, genes = st[key]
    path_cuts, ben_cuts = st["cuts"]
    return fit_subset(y, s, genes, train, path_cuts, ben_cuts, cell)


def plan(df: pd.DataFrame) -> list[tuple]:
    """Every fit to make: (tool, stratum, G, k, training genes, cell key)."""
    tasks = []
    for tool in TOOLS:
        for stratum in STRATA:
            _, _, genes = cell_arrays(df, tool, stratum)
            ug = genes_of(genes)
            for k in list(range(MIN_K, len(ug))) + [len(ug)]:
                for train in itertools.combinations(ug, k):
                    tasks.append((tool, stratum, len(ug), k, train,
                                  cell_key(tool, stratum, train, ug)))
    return tasks


def run_fits(tasks: list[tuple], jobs: int) -> list[tuple[str, dict]]:
    work = [(t[0], t[1], t[4], t[5]) for t in tasks]
    block_end = {}
    for i, t in enumerate(tasks):
        block_end[(t[0], t[1])] = i
    t0 = time.perf_counter()

    def progress(i):
        t = tasks[i]
        if block_end[(t[0], t[1])] == i:
            n = sum(1 for u in tasks if (u[0], u[1]) == (t[0], t[1]))
            print(f"[E11] {t[0]:16s} {t[1]:7s} G={t[2]} fits={n:>4}  "
                  f"({i + 1}/{len(tasks)} done, {time.perf_counter() - t0:.0f}s)",
                  flush=True)

    results = []
    if jobs <= 1:
        for i, w in enumerate(work):
            results.append(_run_task(w))
            progress(i)
    else:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for i, r in enumerate(ex.map(_run_task, work, chunksize=1)):
                results.append(r)
                progress(i)
    return results


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------
def build_tables(df: pd.DataFrame, tasks: list[tuple], results: list[tuple],
                 path_cuts: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(fits, detail, summary). `fits` is one row per subset x tier, including the
    k = G reference; `detail` and `summary` hold k = 2..G-1 only."""
    arrays = {}
    fit_rows, detail_rows = [], []
    for (tool, stratum, G, k, train, cell), (status, taus) in zip(tasks, results):
        if (tool, stratum) not in arrays:
            arrays[(tool, stratum)] = cell_arrays(df, tool, stratum)
        y, s, genes = arrays[(tool, stratum)]
        ug = genes_of(genes)
        tr = np.isin(genes, list(train))
        base = {"tool": tool, "stratum": stratum, "n_genes": G, "k": k,
                "training_genes": ";".join(train),
                "training_n_damaging": int((y[tr] == 1).sum()),
                "training_n_normal": int((y[tr] == 0).sum()),
                "seed_key": "|".join(cell[2:3]), "fit_status": status}
        for tier in TIERS:
            tau = taus.get(tier, np.nan)
            row = base | {"tier": tier, "path_cut_lr": path_cuts[tier],
                          "threshold": tau, "tier_reached": bool(np.isfinite(tau))}
            fit_rows.append(row)
            if k == G:
                continue
            for g in ug:
                if g in train:
                    continue
                ev = heldout_eval(y, s, genes, g, tau)
                if ev["evaluable"]:
                    ev["cleared_cut"] = bool(ev["heldout_lr"] >= path_cuts[tier])
                detail_rows.append(row | {"heldout_gene": g} | ev)

    fits = pd.DataFrame(fit_rows)
    detail = pd.DataFrame(detail_rows)
    for c in ("heldout_n_in_band", "heldout_n_damaging_in_band",
              "heldout_n_normal_in_band"):
        detail[c] = detail[c].astype("Int64")
    detail["cleared_cut"] = detail["cleared_cut"].astype("boolean")
    detail = detail[["tool", "stratum", "tier", "path_cut_lr", "n_genes", "k",
                     "training_genes", "training_n_damaging", "training_n_normal",
                     "seed_key", "fit_status", "threshold", "tier_reached",
                     "heldout_gene", "heldout_n_damaging", "heldout_n_normal",
                     "heldout_n_in_band", "heldout_n_damaging_in_band",
                     "heldout_n_normal_in_band", "heldout_lr", "evaluable",
                     "cleared_cut", "not_evaluable_reason"]]
    return fits, detail, summarise(fits, detail)


def _q(x: pd.Series, q: float) -> float:
    return float(np.quantile(x, q)) if len(x) else np.nan


def summarise(fits: pd.DataFrame, detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tool in TOOLS:
        for stratum in STRATA:
            for tier in TIERS:
                f = fits[(fits.tool == tool) & (fits.stratum == stratum)
                         & (fits.tier == tier)]
                if not len(f):
                    continue
                G = int(f["n_genes"].iloc[0])
                ref = f[f.k == G]["threshold"]
                tau_all = float(ref.iloc[0]) if len(ref) else np.nan
                for k in range(MIN_K, G):
                    fk = f[f.k == k]
                    taus = fk.loc[fk.tier_reached, "threshold"]
                    d = detail[(detail.tool == tool) & (detail.stratum == stratum)
                               & (detail.tier == tier) & (detail.k == k)]
                    ev = d[d.evaluable]
                    lrs = ev["heldout_lr"]
                    n_clear = int(ev["cleared_cut"].sum())
                    rows.append({
                        "tool": tool, "stratum": stratum, "tier": tier,
                        "path_cut_lr": float(fk["path_cut_lr"].iloc[0]),
                        "n_genes": G, "k": k,
                        "n_subsets": int(len(fk)),
                        "n_subsets_fitted": int((fk.fit_status == "fitted").sum()),
                        "n_subsets_reached": int(len(taus)),
                        "frac_subsets_reached": len(taus) / len(fk),
                        "threshold_median": _q(taus, 0.5),
                        "threshold_q25": _q(taus, 0.25),
                        "threshold_q75": _q(taus, 0.75),
                        "threshold_min": float(taus.min()) if len(taus) else np.nan,
                        "threshold_max": float(taus.max()) if len(taus) else np.nan,
                        "threshold_all_genes": tau_all,
                        "n_pairs": int(len(d)),
                        "n_pairs_evaluable": int(len(ev)),
                        "n_pairs_cleared": n_clear,
                        "frac_pairs_cleared": (n_clear / len(ev)) if len(ev) else np.nan,
                        "heldout_lr_median": _q(lrs, 0.5),
                        "heldout_lr_q25": _q(lrs, 0.25),
                        "heldout_lr_q75": _q(lrs, 0.75),
                    })
    return pd.DataFrame(rows)


def _sig(x, n: int):
    return float(f"{x:.{n}g}") if pd.notna(x) else np.nan


def for_print(summary: pd.DataFrame) -> pd.DataFrame:
    """The summary rounded for the printed supplement. The unrounded table stays the
    source of record; this is a view of it, recomputed rather than hand-copied."""
    t = summary[_SUPP_COLS].copy()
    for c in ("threshold_median", "threshold_q25", "threshold_q75",
              "threshold_all_genes"):
        t[c] = t[c].round(_DP_THRESHOLD)
    t["heldout_lr_median"] = t["heldout_lr_median"].map(lambda v: _sig(v, _SIG_LR))
    for c in ("frac_subsets_reached", "frac_pairs_cleared"):
        t[c] = t[c].round(_DP_FRACTION)
    return t


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def _print_readout(summary: pd.DataFrame, fits: pd.DataFrame) -> None:
    for tier in TIERS:
        s = summary[summary.tier == tier]
        print(f"\n--- {tier}: fraction of subsets reaching the tier "
              f"| fraction of evaluable held-out pairs clearing its cut, by k ---")
        for (tool, stratum), d in s.groupby(["tool", "stratum"], sort=False):
            cells = "  ".join(
                f"k={r.k}: {r.frac_subsets_reached:.2f}|"
                + (f"{r.frac_pairs_cleared:.2f}" if pd.notna(r.frac_pairs_cleared)
                   else " na ") + f" (n={r.n_pairs_evaluable})"
                for r in d.itertuples())
            print(f"  {tool:16s} {stratum:7s} {cells}")

    # composition: at each k, the reach rate of subsets that include a gene minus
    # that of subsets that do not, averaged over k with equal weight
    print("\n--- composition: reach rate with the gene in training minus without, "
          "averaged over k ---")
    sub = fits[fits.k < fits.n_genes]
    for (tool, stratum, tier), f in sub.groupby(["tool", "stratum", "tier"],
                                                sort=False):
        genes = sorted(set(";".join(f.training_genes).split(";")))
        parts = []
        for g in genes:
            has = f.training_genes.str.split(";").map(lambda gs, g=g: g in gs)
            diffs = [f[has & (f.k == k)].tier_reached.mean()
                     - f[~has & (f.k == k)].tier_reached.mean()
                     for k in sorted(f.k.unique())
                     if (has & (f.k == k)).any() and (~has & (f.k == k)).any()]
            parts.append(f"{g}={np.mean(diffs):+.2f}")
        print(f"  {tool:16s} {stratum:7s} {tier:9s} " + "  ".join(parts))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1),
                    help="worker processes; the outputs do not depend on it")
    args = ap.parse_args()

    st = _state()
    df, (path_cuts, _) = st["df"], st["cuts"]
    tasks = plan(df)
    print(f"[E11] {len(tasks)} fits across {len(TOOLS)} tools x {len(STRATA)} "
          f"strata, {N_BOOT} gene-clustered resamples each, jobs={args.jobs}",
          flush=True)
    results = run_fits(tasks, args.jobs)
    fits, detail, summary = build_tables(df, tasks, results, path_cuts)

    # the subset counts are fixed by G; a shortfall means a fit went missing
    for r in summary.itertuples():
        assert r.n_subsets == comb(r.n_genes, r.k), (r.tool, r.stratum, r.k)

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    detail.to_csv(OUT, index=False)
    summary.to_csv(OUT_SUMMARY, index=False)
    # The printed Supplementary Table S8 is written by evid_supp_tables from the
    # summary, so that it carries the same labels and the same threshold precision
    # as every other supplementary table.
    print(f"\n[E11] wrote {OUT} ({len(detail)} rows) and {OUT_SUMMARY} "
          f"({len(summary)} rows)")
    _print_readout(summary, fits)


if __name__ == "__main__":
    main()
