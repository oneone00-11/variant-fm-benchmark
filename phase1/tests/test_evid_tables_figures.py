"""The printed tables and the figures are derived, not typed.

These tests hold the main tables to the upstream files they are built from, and
check that the arm comparison never reports the published SpliceAI cut point for
a column it is not defined on. The figure check is existence and determinism of
metadata only; the figures' content is covered by the tables they are drawn from.
"""
from __future__ import annotations

from pathlib import Path

import sys

import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
EV = PHASE1 / "reports/evidence"
TABLES = EV / "tables"
FIGS = EV / "figures"

need = pytest.mark.skipif(not (TABLES / "table3_evidence_tiers.csv").exists(),
                          reason="tables not built")


@need
def test_table3_fold_counts_match_the_tier_table():
    """Every 'k/n' cell of Table 3 is folds clearing over folds, as tier_logo_table says."""
    t3 = pd.read_csv(TABLES / "table3_evidence_tiers.csv", dtype=str).set_index("Predictor")
    t = pd.read_csv(EV / "tier_logo_table.csv")
    from src.evid_figures import NAME
    band = {"s3_10": "3–10 bp", "s11_50": "11–50 bp", "s3_50": "3–50 bp"}
    checked = 0
    for r in t[(t.status == "ok") & t.in_sample_tier_reached.astype(bool)].itertuples():
        if r.tool not in NAME:
            continue
        cell = t3.loc[NAME[r.tool], f"{band[r.stratum]}: {r.tier.capitalize()}, held out"] \
            if r.tier in ("moderate", "strong") else None
        if cell is None:
            continue
        k, n = cell.split(" ")[0].split("/")
        assert int(k) == r.folds_heldout_lr_above_cut and int(n) == r.n_folds, (r.tool, r.stratum, r.tier)
        checked += 1
    assert checked > 20


@need
def test_table1_totals_add_up_to_the_count_table():
    a = pd.read_csv(TABLES / "table1a_composition_by_gene.csv", dtype=str)
    sc = pd.read_csv(EV / "set_counts.csv")
    ins = sc[sc.stratum != "pm12"]
    total = a[a.Gene == "Total"].iloc[0]
    assert int(total["Variants, 3–50 bp"].replace(",", "")) == int(ins.n.sum())
    assert int(total["3–50 bp: labelled"].replace(",", "")) == int(ins.n_labelled.sum())


@need
def test_the_published_cut_point_is_only_applied_to_spliceai():
    a = pd.read_csv(EV / "arms_at_tool_threshold.csv")
    pub = a[a.threshold_basis == "published cut point"]
    assert set(pub.tool) <= {"spliceai", "spliceai_walker"}


@need
def test_table4a_matches_the_arm_file():
    t = pd.read_csv(TABLES / "table4a_spliceai_published_cut_by_arm.csv", dtype=str)
    a = pd.read_csv(EV / "arms_at_tool_threshold.csv")
    r = a[(a.tool == "spliceai_walker") & (a.threshold_basis == "published cut point")
          & (a.gene_set == "all_genes") & (a.stratum == "s3_50")
          & (a.clinvar_arm == "classified")].iloc[0]
    cell = t[(t.Genes == "Seven genes") & (t.Band == "3–50 bp")].iloc[0]["Classified"]
    assert float(cell) == round(r.band_lr, 1)


@pytest.mark.skipif(not FIGS.exists(), reason="figures not built")
def test_every_figure_exists_in_both_formats_without_a_timestamp():
    for stem in ("figure1", "figure2", "figure3", "figure4", "figureS1", "figureS2", "figureS3"):
        pdf, png = FIGS / f"{stem}.pdf", FIGS / f"{stem}.png"
        assert pdf.exists() and png.exists(), stem
        head = pdf.read_bytes()[:4000] + pdf.read_bytes()[-4000:]
        assert b"CreationDate" not in head, f"{stem}.pdf carries a creation date"
