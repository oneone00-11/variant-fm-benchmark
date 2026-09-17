"""E8 -- the two diagnostic tables that used to be written by the test suite.

`walker_cutpoint_local_vs_band.csv` and `local_lr_monotonicity.txt` are cited in
the report, so they are pipeline outputs and belong to a module that an entry point
runs. They were originally emitted as a side effect of two tests, which meant
nothing regenerated them and running pytest rewrote tracked files.

The tests now read these files instead of writing them.

Run (PYTHONPATH=phase1):  python -m src.evid_diagnostics
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPORT_DIR = Path("reports/evidence")


def cutpoint_local_vs_band() -> pd.DataFrame | None:
    """At Walker's cut point: the ratio over the band above it, against the density
    ratio in a window around it. They are different quantities and they diverge."""
    w, e = REPORT_DIR / "walker_thresholds.csv", REPORT_DIR / "evidence_thresholds.csv"
    if not (w.exists() and e.exists()):
        return None
    e2, e3 = pd.read_csv(w), pd.read_csv(e)
    if "walker_pp3_local_lr" not in e3.columns:
        return None
    e3 = e3[(e3.status == "ok") & e3["walker_pp3_local_lr"].notna()]
    rows = []
    for tool in e3.tool.unique():
        for stratum in e3[e3.tool == tool].stratum.unique():
            a = e2[(e2.score_column == tool) & (e2.stratum == stratum)
                   & (e2.scope == "pooled") & (e2.clinvar_arm == "all")
                   & (e2.status == "ok")]
            b = e3[(e3.tool == tool) & (e3.stratum == stratum)]
            if not len(a) or not len(b):
                continue
            rows.append({
                "tool": tool, "stratum": stratum,
                "walker_in_scope": bool(a.walker_in_scope.iloc[0]),
                "frac_called_pp3": round(float(a.frac_pp3.iloc[0]), 4),
                "band_lr_above_cut": round(float(a.lr_pp3.iloc[0]), 3),
                "local_lr_at_cut": round(float(b.walker_pp3_local_lr.iloc[0]), 3),
                "local_lr_lower_bound": round(float(b.walker_pp3_local_lr_lo.iloc[0]), 3),
            })
    return pd.DataFrame(rows) if rows else None


def local_lr_monotonicity() -> str | None:
    """The threshold rule reads the local likelihood ratio as rising with the score.
    Where a curve dips, the dip is written down rather than assumed away."""
    curves = sorted(REPORT_DIR.glob("interval_lr_*.csv"))
    if not curves:
        return None
    lines = ["tool_stratum n_grid n_drops drop_frac end_minus_start"]
    for path in curves:
        c = pd.read_csv(path).dropna(subset=["local_lr"])
        if len(c) < 3:
            continue
        d = np.diff(c["local_lr"].to_numpy())
        drops = int((d < 0).sum())
        lines.append(" ".join(str(x) for x in (
            path.stem, len(c), drops, round(drops / max(len(d), 1), 3),
            float(c["local_lr"].iloc[-1] - c["local_lr"].iloc[0]))))
    return "\n".join(lines) + "\n" if len(lines) > 1 else None


def main() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    cv = cutpoint_local_vs_band()
    if cv is not None:
        cv.to_csv(REPORT_DIR / "walker_cutpoint_local_vs_band.csv", index=False)
        print(f"[E8] wrote walker_cutpoint_local_vs_band.csv ({len(cv)} rows)")
        print(cv.to_string(index=False))
    else:
        print("[E8] walker_cutpoint_local_vs_band.csv: inputs missing, not written")
    mono = local_lr_monotonicity()
    if mono is not None:
        (REPORT_DIR / "local_lr_monotonicity.txt").write_text(mono)
        print(f"\n[E8] wrote local_lr_monotonicity.txt ({len(mono.splitlines()) - 1} curves)")
    else:
        print("[E8] local_lr_monotonicity.txt: no curve files, not written")


if __name__ == "__main__":
    main()
