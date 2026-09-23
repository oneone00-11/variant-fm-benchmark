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
import yaml

from . import evid_common as K

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
    _main_concordance()
    dc = depth_counts()
    dc.to_csv(REPORT_DIR / "depth_counts.csv", index=False)
    print(f"\n[E8] wrote depth_counts.csv ({len(dc)} rows)")
    print(dc[dc.gene == "all"].to_string(index=False))


# ---------------------------------------------------------------------------
# E8c -- how deep into the intron the labelled data reach
# ---------------------------------------------------------------------------
# The deposits stop at different depths, and the distal band's damaging variants
# sit almost entirely in its first twenty nucleotides. A statement about "11-50 bp"
# is only as deep as the data under it, so the depth profile is written here and
# every depth count the text quotes is read from it.
DEPTH_BANDS = [("3-10", 3, 10), ("11-20", 11, 20), ("21-30", 21, 30),
               ("31-40", 31, 40), ("41-50", 41, 50), ("11-30", 11, 30),
               ("31-50", 31, 50), ("11-50", 11, 50), ("3-50", 3, 50)]


def depth_counts() -> pd.DataFrame:
    cfg = yaml.safe_load(Path("config/walker2023.yaml").read_text())
    pp3 = cfg["thresholds"]["pp3"]["value"]
    df = K.load_set()
    rows = []
    for gene in ["all"] + sorted(df["gene"].unique()):
        g = df if gene == "all" else df[df["gene"] == gene]
        for name, lo, hi in DEPTH_BANDS:
            b = g[(g["intron_offset_abs"] >= lo) & (g["intron_offset_abs"] <= hi)]
            lab = b[b["y_assay"].notna()]
            called = lab["spliceai_walker"] >= pp3
            rows.append({
                "gene": gene, "depth_nt": name, "n": int(len(b)),
                "labelled": int(len(lab)),
                "damaging": int((lab["y_assay"] == 1).sum()),
                "normal": int((lab["y_assay"] == 0).sum()),
                "at_or_above_pp3_cut": int(called.sum()),
                "damaging_at_or_above_pp3_cut": int((called & (lab["y_assay"] == 1)).sum()),
                "max_offset_in_gene": int(g["intron_offset_abs"].max()),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# E8b -- concordance between columns that measure related things
# ---------------------------------------------------------------------------
# The manuscript states how far the two SpliceAI bases agree (the atlas's -D 50
# column against the Walker-basis -D 4999 re-score) and how far the Atlas's
# precomputed AVI columns agree with the model-computed alphagenome column. Those
# numbers were computed ad hoc for the execution reports; a number the text cites
# has to come from a pipeline output, so they are written here.
_PAIRS = [
    ("spliceai", "spliceai_walker"),
    ("avi", "alphagenome"),
    ("avi_splice_sites", "alphagenome"),
    ("avi_splice_site_usage", "alphagenome"),
    ("avi_splice_junctions", "alphagenome"),
    ("avi", "alphagenome_v061"),
]


def column_concordance() -> pd.DataFrame | None:
    from scipy.stats import spearmanr
    import yaml
    set_path = Path("data/evidence/analysis_set_v1.parquet")
    cfg_path = Path("config/walker2023.yaml")
    if not (set_path.exists() and cfg_path.exists()):
        return None
    df = pd.read_parquet(set_path)
    pp3 = float(yaml.safe_load(cfg_path.read_text())["thresholds"]["pp3"]["value"])
    rows = []
    for a, b in _PAIRS:
        if a not in df.columns or b not in df.columns:
            continue
        x = df[a].to_numpy(dtype=float)
        y = df[b].to_numpy(dtype=float)
        m = np.isfinite(x) & np.isfinite(y)
        row = {"column_a": a, "column_b": b, "n": int(m.sum()),
               "spearman_rho": float(spearmanr(x[m], y[m]).statistic)}
        if a == "spliceai" and b == "spliceai_walker":
            # the same tool on two bases: how many variants change side of the
            # PP3 cut point, in each direction
            row["cut_point"] = pp3
            row["n_a_below_b_at_or_above"] = int(((x[m] < pp3) & (y[m] >= pp3)).sum())
            row["n_a_at_or_above_b_below"] = int(((x[m] >= pp3) & (y[m] < pp3)).sum())
        rows.append(row)
    return pd.DataFrame(rows) if rows else None


def _main_concordance() -> None:
    cc = column_concordance()
    if cc is not None:
        cc.to_csv(REPORT_DIR / "column_concordance.csv", index=False)
        print(f"\n[E8] wrote column_concordance.csv ({len(cc)} pairs)")
        print(cc.to_string(index=False))
    else:
        print("[E8] column_concordance.csv: inputs missing, not written")


if __name__ == "__main__":
    main()
