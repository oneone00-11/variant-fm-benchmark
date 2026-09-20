"""E2.2 -- does a tier hold when the threshold is chosen on other genes?

Two quantities that are easy to conflate, so they are named apart and put in
separate columns:

  in_sample_tier_reached  Pejaver's rule on the whole stratum: the tier is reached
                          when there is a score at or above which EVERY local
                          likelihood ratio's one-sided 95% lower bound clears the
                          tier's cut. This is the claim the manuscript leads with.
  heldout_lr              a point estimate, nothing more: the threshold fitted on
                          six genes, applied to the seventh, and the likelihood
                          ratio of the whole band above it measured there. No
                          bound, no bootstrap.

A tier can be reached in-sample and still produce a held-out ratio below its own
cut, and that is not a contradiction -- one is a bounded statement about the
stratum, the other an unbounded one about one gene.

The fusion is excluded. Its threshold does not carry across genes (held-out ratios
spanning 3.5 to 106.7 in the published table) and its LOGO fit is contaminated: the
training folds' out-of-fold scores came from models that saw the held-out gene. The
per-fold table stays in the supplement as the evidence for that, and is not
repeated here as if it were a transferable threshold.

Run (PYTHONPATH=phase1):  python -m src.evid_tier_logo
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import evid_common as K

REPORT_DIR = Path("reports/evidence")
OUT = REPORT_DIR / "tier_logo_table.csv"
TIERS = ["supporting", "moderate", "strong"]
STRATA = ["s3_10", "s11_50", "s3_50"]


def build() -> pd.DataFrame:
    ev = pd.read_csv(REPORT_DIR / "evidence_thresholds.csv")
    folds = pd.read_csv(REPORT_DIR / "evidence_thresholds_logo_folds.csv")
    ev = ev[ev["tool"] != K.FUSION]
    folds = folds[folds["tool"] != K.FUSION]

    rows = []
    for tool in sorted(ev["tool"].unique()):
        for stratum in STRATA:
            for tier in TIERS:
                e = ev[(ev.tool == tool) & (ev.stratum == stratum)
                       & (ev.tier == tier)]
                if not len(e) or e["status"].iloc[0] != "ok":
                    rows.append({"tool": tool, "stratum": stratum, "tier": tier,
                                 "status": "not evaluable"})
                    continue
                e = e.iloc[0]
                f = folds[(folds.tool == tool) & (folds.stratum == stratum)
                          & (folds.tier == tier) & (folds.side == "pp3")]
                per_gene = f.sort_values("heldout_gene")
                taus = per_gene["threshold_from_6_genes"]
                lrs = per_gene["heldout_lr"]
                rows.append({
                    "tool": tool, "stratum": stratum, "tier": tier, "status": "ok",
                    # Pejaver's rule on the whole stratum -- the primary claim
                    "in_sample_tier_reached": bool(e["pp3_threshold_reachable"]),
                    "in_sample_threshold": e["pp3_threshold_insample"],
                    "tier_lr_cut": e["path_cut_lr"],
                    # leave-one-gene-out, per fold, as point estimates
                    "n_folds": int(len(per_gene)),
                    "folds_reached": int(taus.notna().sum()),
                    "folds_with_heldout_lr": int(lrs.notna().sum()),
                    "heldout_lr_median": float(lrs.median()) if lrs.notna().any() else np.nan,
                    "heldout_lr_min": float(lrs.min()) if lrs.notna().any() else np.nan,
                    "heldout_lr_max": float(lrs.max()) if lrs.notna().any() else np.nan,
                    "folds_heldout_lr_above_cut": int(
                        (lrs >= e["path_cut_lr"]).sum()) if lrs.notna().any() else 0,
                    "threshold_from_6_genes_per_fold": "; ".join(
                        f"{g}={t:.4g}" if np.isfinite(t) else f"{g}=none"
                        for g, t in zip(per_gene["heldout_gene"], taus)),
                    "heldout_lr_per_fold": "; ".join(
                        f"{g}={v:.4g}" if np.isfinite(v) else f"{g}=na"
                        for g, v in zip(per_gene["heldout_gene"], lrs)),
                })
    return pd.DataFrame(rows)


def main() -> None:
    out = build()
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    ok = out[out.status == "ok"]
    print(f"[E2.2] wrote {OUT} ({len(out)} rows; fusion excluded by design)\n")
    print("--- Strong: reached in-sample, and how many folds its threshold reached ---")
    s = ok[ok.tier == "strong"]
    print(s[["tool", "stratum", "in_sample_tier_reached", "in_sample_threshold",
             "folds_reached", "n_folds", "heldout_lr_median",
             "folds_heldout_lr_above_cut"]].round(4).to_string(index=False))
    print("\n--- in-sample reached by tier and stratum ---")
    print(ok.pivot_table(index="tool", columns=["stratum", "tier"],
                         values="in_sample_tier_reached",
                         aggfunc="first").to_string())


if __name__ == "__main__":
    main()
