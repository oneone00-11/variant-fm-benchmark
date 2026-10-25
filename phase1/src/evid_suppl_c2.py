"""E19 C2 -- how ClinVar covers the analysis set, position by position.

The denominator is every SNV of the analysis set, 8,853, the 1-2 nt canonical
dinucleotides included: the assays measured all of them, so the share ClinVar records
at each position is measured on an unselected set. Positions are the signed intronic
offset, -50 ... -1 on the acceptor side and +1 ... +50 on the donor side, read from
the c. notation by the companion atlas's own offset pattern (atlas/evaluate.py,
_OFFSET), with the sign kept; its magnitude is the set's intron_offset_abs.

ClinVar's fields are the analysis set's own, from the VCF of 2026-06-15
(evid_build_set.annotate_clinvar): `clinvar_arm` (unrecorded, recorded but not
classified P/LP or B/LB, classified) and `clnsig_class`, the project's collapse of
the germline classification CLNSIG into P/LP, B/LB, VUS, Conflicting and Other (Other:
a record without a germline classification, or not_provided). Nothing is fitted and no
threshold is read.

Outputs, in reports/evidence/suppl_analyses/:
  c2_clinvar_coverage_by_position.csv   per gene and for the seven genes together, at
                                        each signed offset: the variants, those with
                                        any ClinVar record, those classified P/LP or
                                        B/LB, and P/LP, B/LB, VUS, Conflicting and
                                        Other, each as a count and a share of the
                                        variants at that offset
  c2_clinvar_coverage_by_stratum.csv    the same genes, by band (1-2, 3-10, 11-50 nt):
                                        unrecorded, recorded but unclassified, and
                                        classified, as counts and shares
  c2_clinvar_coverage_draft.png         a draft to look at, in evid_fig_data's manner;
                                        the figure's layout is left to the revision

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_c2
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

from . import evid_common as K            # noqa: E402
from . import evid_figures as F           # noqa: E402  (the typeface)

OUT_DIR = Path("reports/evidence/suppl_analyses")
OUT_POSITION = OUT_DIR / "c2_clinvar_coverage_by_position.csv"
OUT_STRATUM = OUT_DIR / "c2_clinvar_coverage_by_stratum.csv"
OUT_DRAFT = OUT_DIR / "c2_clinvar_coverage_draft.png"

# the atlas's offset pattern (atlas/evaluate.py, _OFFSET), the sign captured as well
SIGNED_OFFSET = re.compile(r"^c\.[-*]?\d+([+-])(\d+)")
ALL = "Seven genes"
CLASSES = [("P/LP", "plp"), ("B/LB", "blb"), ("VUS", "vus"),
           ("Conflicting", "conflicting"), ("Other", "other")]
# the bands of the set, in the order the text reads them; pm12 is the 1-2 nt band
BANDS = [("pm12", "1–2"), ("s3_10", "3–10"), ("s11_50", "11–50")]
ARMS = ["unrecorded", "recorded_unclassified", "classified"]


def signed_offset(hgvs_c: str) -> int:
    m = SIGNED_OFFSET.search(hgvs_c)
    if m is None:
        raise SystemExit(f"[C2] no intronic offset in {hgvs_c!r}")
    return int(m.group(2)) * (1 if m.group(1) == "+" else -1)


def _share(k: int, n: int) -> float:
    return k / n if n else np.nan


def _cell(d: pd.DataFrame) -> dict:
    n = int(len(d))
    recorded = int((d["clinvar_arm"] != "unrecorded").sum())
    classified = int((d["clinvar_arm"] == "classified").sum())
    row = {"n_variants": n,
           "n_recorded": recorded, "share_recorded": _share(recorded, n),
           "n_classified": classified, "share_classified": _share(classified, n)}
    for cls, short in CLASSES:
        k = int((d["clnsig_class"] == cls).sum())
        row |= {f"n_{short}": k, f"share_{short}": _share(k, n)}
    return row


def by_position(df: pd.DataFrame) -> pd.DataFrame:
    positions = [*range(-50, 0), *range(1, 51)]
    rows = []
    groups = [(g, df[df["gene"] == g]) for g in sorted(df["gene"].unique())] + [(ALL, df)]
    for gene, d in groups:
        for pos in positions:
            rows.append({"gene": gene, "signed_offset": pos,
                         "side": "donor" if pos > 0 else "acceptor"}
                        | _cell(d[d["signed_offset"] == pos]))
    return pd.DataFrame(rows)


def by_stratum(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    groups = [(g, df[df["gene"] == g]) for g in sorted(df["gene"].unique())] + [(ALL, df)]
    for gene, d in groups:
        for stratum, band in BANDS:
            s = d[d["stratum"] == stratum]
            n = int(len(s))
            row = {"gene": gene, "stratum": stratum, "band_nt": band, "n_variants": n}
            for arm in ARMS:
                k = int((s["clinvar_arm"] == arm).sum())
                row |= {f"n_{arm}": k, f"share_{arm}": _share(k, n)}
            rows.append(row)
    return pd.DataFrame(rows)


def draft(pos: pd.DataFrame) -> None:
    F._use_font()
    d = pos[pos["gene"] == ALL].set_index("signed_offset")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax1.plot(d.index, d["share_recorded"], "o-", ms=3, label="any ClinVar record")
    ax1.plot(d.index, d["share_classified"], "s-", ms=3, label="classified P/LP or B/LB")
    ax1.set_ylabel("share of variants"); ax1.set_ylim(0, 1); ax1.legend(fontsize=7)
    ax1.set_title("C2 draft: ClinVar coverage by signed intronic offset, seven genes")
    bottom = np.zeros(len(d))
    for cls, short in CLASSES:
        ax2.bar(d.index, d[f"share_{short}"].fillna(0), bottom=bottom, width=0.8, label=cls)
        bottom += d[f"share_{short}"].fillna(0).to_numpy()
    ax2.set_ylabel("share of variants"); ax2.set_ylim(0, 1); ax2.legend(fontsize=7, ncol=5)
    for ax in (ax1, ax2):
        for x in (-10.5, -2.5, 2.5, 10.5):
            ax.axvline(x, ls=":", lw=0.8, color="grey")
    ax2.set_xlabel("signed intronic offset (acceptor side < 0 < donor side)")
    fig.tight_layout(); fig.savefig(OUT_DRAFT, dpi=110); plt.close(fig)


def main() -> None:
    df = K.load_set()
    df["signed_offset"] = df["hgvs_c"].map(signed_offset)
    if not (df["signed_offset"].abs() == df["intron_offset_abs"]).all():
        raise SystemExit("[C2] the signed offset's magnitude is not the set's intron_offset_abs")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pos = by_position(df)
    pos.to_csv(OUT_POSITION, index=False)
    st = by_stratum(df)
    st.to_csv(OUT_STRATUM, index=False)
    draft(pos)
    total = pos[pos["gene"] == ALL]["n_variants"].sum()
    print(f"[C2] wrote {OUT_POSITION} ({len(pos)} rows; {total:,} variants), {OUT_STRATUM} "
          f"({len(st)} rows) and {OUT_DRAFT}")
    s = st[st["gene"] == ALL]
    print(s[["band_nt", "n_variants"] + [f"share_{a}" for a in ARMS]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
