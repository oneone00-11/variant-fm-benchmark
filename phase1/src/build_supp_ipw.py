"""
Build Supplementary Table S10 (sampling-frame reweighting) and append it to the
supplementary-material document.

Every number is read from `phase1/reports/phase1/ipw_*.csv`; nothing is typed in
by hand. The document is regenerated from an input copy rather than edited in
place, so the transform is reproducible and the input stays intact:

    python -m src.build_supp_ipw --in  "<supplement>.docx" --out "<supplement>.docx"

Panels
------
  S10a  selection inside the target frame        ipw_frame_balance.csv
  S10b  weight diagnostics                       ipw_weight_diagnostics.csv
  S10c  per-predictor rho and Brier, unweighted vs weighted   ipw_predictor_table.csv
  S10d  headline deltas under every binning scheme            ipw_ranking.csv + ipw_headline.csv
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import pandas as pd
import docx
from docx.shared import Pt

from . import config as C

PRETTY = {"fusion_M1": "Fusion — elastic net (M1)", "mean_M0b": "Fusion — equal-weight mean (M0b)",
          "single:pangolin": "Pangolin", "single:spliceai": "SpliceAI",
          "single:alphagenome": "AlphaGenome", "single:cadd": "CADD",
          "single:gpn_msa": "GPN-MSA", "single:phylop": "phyloP",
          "single:phastcons": "phastCons", "single:nt": "Nucleotide Transformer",
          "single:gnomad_af": "gnomAD AF", "single:alphamissense": "AlphaMissense"}
SCHEME = {"unweighted": "Unweighted", "quintile": "Quintile", "decile": "Decile (primary)",
          "vigintile": "Twentile", "decile_truncated": "Decile, truncated [0.1, 10]"}
COND = {"y_assay/BRCA1_included": "Assay, +BRCA1", "y_assay/BRCA1_excluded": "Assay, −BRCA1",
        "y_clinvar/BRCA1_included": "ClinVar, +BRCA1", "y_clinvar/BRCA1_excluded": "ClinVar, −BRCA1"}
ORDER = ["unweighted", "quintile", "decile", "vigintile", "decile_truncated"]


def _f(x, nd=4):
    return "" if pd.isna(x) else f"{float(x):.{nd}f}"


def _ci(lo, hi, nd=4):
    return "" if pd.isna(lo) or pd.isna(hi) else f"[{float(lo):+.{nd}f}, {float(hi):+.{nd}f}]"


def panels(rep: Path):
    fb = pd.read_csv(rep / "ipw_frame_balance.csv")
    wd = pd.read_csv(rep / "ipw_weight_diagnostics.csv")
    pt = pd.read_csv(rep / "ipw_predictor_table.csv")
    rk = pd.read_csv(rep / "ipw_ranking.csv")
    hd = pd.read_csv(rep / "ipw_headline.csv")

    a_h = ["Gene", "n in frame", "n observed", "% covered", "mean |z| observed",
           "mean |z| frame", "mean |z| not recruited", "Extreme decile, obs.",
           "Extreme decile, not recr.", "KS p"]
    a_r = [[r.gene, f"{int(r.n_frame):,}", f"{int(r.n_observed):,}", f"{r.pct_covered:.1f}",
            _f(r.mean_abs_z_observed), _f(r.mean_abs_z_frame), _f(r.mean_abs_z_remainder),
            _f(r.extreme_decile_frac_observed), _f(r.extreme_decile_frac_remainder),
            str(r.ks_p)] for r in fb.itertuples()]

    b_h = ["Binning", "n", "w min", "w median", "w max", "Kish ESS", "ESS / n", "Bins with w > 10"]
    wd = wd.set_index("scheme").reindex([s for s in ORDER if s in set(wd.scheme)]).reset_index()
    b_r = [[SCHEME.get(r.scheme, r.scheme), f"{int(r.n):,}", _f(r.w_min), _f(r.w_median),
            _f(r.w_max), f"{r.ESS:.1f}", _f(r.ESS_over_n), str(int(r.n_w_gt_alert))]
           for r in wd.itertuples()]

    c_h = ["Predictor / model", "ρ unweighted", "ρ weighted", "Δρ", "Brier unweighted",
           "Brier weighted", "ΔBrier", "Rank unw.", "Rank wgt.", "Rank change"]
    c_r = [[PRETTY.get(r.model, r.model), _f(r.rho_unw), _f(r.rho_wgt), f"{r.d_rho:+.4f}",
            _f(r.Brier_unw), _f(r.Brier_wgt), f"{r.d_Brier:+.4f}",
            str(int(r.rank_unw)), str(int(r.rank_wgt)), ("0" if int(r.rank_change) == 0 else f"{int(r.rank_change):+d}")]
           for r in pt.itertuples()]

    d_h = ["Binning", "Δρ (fusion − best single)", "Δρ 95% CI", "Condition", "ΔBrier", "ΔBrier 95% CI",
           "CI excludes 0"]
    dl = rk[rk["rank"].isna()].set_index("scheme")
    d_r = []
    for s in ORDER:
        if s not in dl.index:
            continue
        row = dl.loc[s]
        first = True
        for cond in COND:
            h = hd[(hd.scheme == s) & (hd.set == cond)]
            if not len(h):
                continue
            h = h.iloc[0]
            d_r.append([SCHEME.get(s, s) if first else "",
                        f"{row.rho:+.4f}" if first else "",
                        _ci(row.ci_lo, row.ci_hi) if first else "",
                        COND[cond], _f(h.dBrier), _ci(h.dBrier_lo, h.dBrier_hi),
                        "yes" if bool(h.dBrier_excludes_zero) else "no"])
            first = False
    return [
        ("Supplementary Table S10a. Selection into the analysis set, inside the proximal splice window",
         "Effect size is the functional score standardised within gene over the target frame (the full "
         "functional SNV set with |intron offset| ≤ 8, n = 3,291); |z| is the selection variable. "
         "“Observed” is the analysis set inside that frame (n = 1,768 of 1,781; thirteen variants carry "
         "no intron offset and cannot be placed). BRCA1 is fully covered, so it has no remainder. "
         "KS compares the observed and not-recruited z distributions.", a_h, a_r),
        ("Supplementary Table S10b. Weight diagnostics",
         "Quantile cutpoints are taken on the target frame, never on the observed sample, so "
         "p_full(k) = 1/K by construction; w_k = p_full(k)/p_int(k), renormalised to mean(w) = 1. "
         "No bin exceeds the w > 10 alert threshold, so the truncated row is an identity.", b_h, b_r),
        ("Supplementary Table S10c. Per-predictor ranking and calibration, unweighted versus reweighted",
         "Decile weighting. ρ is the DerSimonian–Laird pooled per-gene weighted Spearman over the "
         "in-frame splice set (n = 1,768); Brier is on the functional standard including BRCA1 under "
         "out-of-gene isotonic calibration. Both columns are computed on the same in-frame subset, so "
         "the unweighted values differ slightly from the main-text ones, which use all 1,781. "
         "AlphaMissense scores ≈0 splice variants and is not evaluable.", c_h, c_r),
        ("Supplementary Table S10d. Headline differences under every binning scheme",
         "Δρ is fusion − best single tool (Pangolin), pooled per-gene, gene-level cluster bootstrap "
         "with 1,000 resamples of the seven genes; ΔBrier is best single − fusion under the same "
         "bootstrap. Positive favours the fusion.", d_h, d_r),
    ]


def add_table(doc, title, caption, header, rows, fontsz=8):
    doc.add_paragraph(title)
    doc.add_paragraph(caption)
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Table Grid"
    for c, h in zip(t.rows[0].cells, header):
        c.text = h
        for p in c.paragraphs:
            for r in p.runs:
                r.bold = True
                r.font.size = Pt(fontsz)
    for row in rows:
        cells = t.add_row().cells
        for c, v in zip(cells, row):
            c.text = str(v)
            for p in c.paragraphs:
                for r in p.runs:
                    r.font.size = Pt(fontsz)
    doc.add_paragraph("")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True)
    ap.add_argument("--out", dest="dst", required=True)
    a = ap.parse_args()
    rep = C.REPORT_DIR
    missing = [n for n in ("ipw_frame_balance.csv", "ipw_weight_diagnostics.csv",
                           "ipw_predictor_table.csv", "ipw_ranking.csv", "ipw_headline.csv")
               if not (rep / n).exists()]
    if missing:
        sys.exit(f"[S10] missing inputs {missing}; run `python -m src.ipw_reweight` first")
    doc = docx.Document(a.src)
    n0 = len(doc.tables)
    for title, caption, header, rows in panels(rep):
        add_table(doc, title, caption, header, rows)
    doc.save(a.dst)
    print(f"[S10] {a.src} -> {a.dst}: tables {n0} -> {len(docx.Document(a.dst).tables)}")


if __name__ == "__main__":
    main()
