"""The list of every analysed variant (src/evid_variant_list.py, E8.5; the author's
decision of 2026-11-07, D5): its six checks, its counts against the tables that count
the same variants, and its form."""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
from src import evid_variant_list as V  # noqa: E402

OUT = PHASE1 / "reports" / "evidence" / "layout"
LIST = OUT / "variants" / "analysed_variants.csv"
COUNTS = OUT / "variants" / "analysed_variants_counts.csv"
need_out = pytest.mark.skipif(not LIST.exists(), reason="the layout stage has not run")


@pytest.fixture(scope="module")
def built(request):
    mp = pytest.MonkeyPatch()
    mp.chdir(PHASE1)
    try:
        df, raw = V.build()
    finally:
        mp.undo()
    return df, raw


@pytest.fixture(scope="module")
def printed():
    return pd.read_csv(LIST, dtype=str, keep_default_na=False)


def _n(s: str) -> int:
    return int(str(s).replace(",", ""))


def test_the_six_checks_pass(built, monkeypatch):
    monkeypatch.chdir(PHASE1)
    df, raw = built
    assert {k: v[:3] for k, v in V.checks(df, raw).items() if v} == {}


def test_the_checks_catch_what_they_are_for(built, monkeypatch):
    """A wrong c., a distance off by one, a swapped base, a row twice, a stray transcript
    and a missing row each fail their check."""
    monkeypatch.chdir(PHASE1)
    df, raw = built
    bad = df.copy()
    raw2 = [dict(r) for r in raw]
    bad.loc[0, "HGVS c."] = "c.100_101del"
    bad.loc[1, "Signed intronic offset"] = str(int(bad.loc[1, "Signed intronic offset"]) + 1)
    bad.loc[2, "Alternate allele"] = bad.loc[2, "Reference allele"]
    raw2[3]["tx"] = "NM_000000.1"
    bad = pd.concat([bad, bad.iloc[[5]]], ignore_index=True)
    raw2.append(raw2[5])
    res = V.checks(bad, raw2)
    assert res["form"] and res["distance"] and res["bases"] and res["transcripts"] and res["unique"] and res["counts"]


@need_out
def test_the_list_is_what_build_gives_and_its_columns_are_all_there(built, printed):
    assert list(printed.columns) == V.COLUMNS
    assert printed.equals(built[0].astype(str))
    assert len(printed) == V.COUNTS["total"] == 15889
    assert set(printed.Set) == {"Analysis set", "External gene"}          # no second-set row
    assert (printed["HGVS g."].str.match(r"^NC_\d{6}\.\d+:g\.\d+[ACGT]>[ACGT]$")).all()
    assert (printed["Variant in the source"].str.match(r"^urn:mavedb:[\w-]+#\d+$")).all()
    assert (printed["Source (MaveDB score set)"] == printed["Variant in the source"].str.split("#").str[0]).all()
    assert not printed["GRCh38 position"].str.contains(",").any()          # a coordinate, as written


@need_out
def test_the_counts_file_counts_the_list(printed):
    c = pd.read_csv(COUNTS, dtype=str, keep_default_na=False)
    keys = ["Set", "Gene", "In the analysis range", "Distance band", "Functional class in the main analysis"]
    want = printed.groupby(keys, sort=False).size()
    got = c.set_index(keys)["Variants"].map(_n)
    assert got.to_dict() == want.to_dict() and got.sum() == len(printed)


@need_out
def test_the_training_genes_are_table_1a_and_the_analysis_set_by_status(printed):
    """In the range, by gene, band and class: Table 1(a)'s labelled and damaging counts,
    and the analysis set by ClinVar status (its rows summed over the statuses)."""
    t1 = pd.read_csv(OUT / "main" / "analysis_set_composition.by_gene.csv", dtype=str, keep_default_na=False)
    t1.columns = [c.replace("\xa0", " ") for c in t1.columns]          # the printed table's spaces
    a = printed[(printed.Set == "Analysis set") & (printed["In the analysis range"] == "Yes")]
    for _, r in t1.iterrows():
        gene = re.sub(r"\*", "", r.iloc[0])
        if gene not in V.TRAINING:
            continue
        g = a[a.Gene == gene]
        assert len(g) == _n(r["Variants at 3–50 bp"]), gene
        for band in ("3–10 bp", "11–50 bp"):
            x = g[g["Distance band"] == band]
            lab = x["Functional class in the main analysis"].isin(["damaging", "normal"]).sum()
            dmg = (x["Functional class in the main analysis"] == "damaging").sum()
            assert lab == _n(r[f"{band}: labelled"]), (gene, band)
            assert r[f"{band}: damaging (%)"].startswith(f"{dmg:,} ") or (dmg == 0 and r[f"{band}: damaging (%)"] in ("0", "–", "")), (gene, band, r[f"{band}: damaging (%)"])
    s3 = pd.read_csv(OUT / "supplementary" / "analysis_set_by_status.csv", dtype=str, keep_default_na=False)
    st = s3["Stratum (intronic offset, bp)"].str.replace("\xa0", " ")
    for gene in V.TRAINING:
        rows = s3[s3.Gene.str.replace("*", "", regex=False) == gene]
        g = printed[printed.Gene == gene]
        for band in ("3–10", "11–50"):
            x = g[g["Distance band"] == f"{band} bp"]
            r = rows[st[rows.index] == band]
            assert r["Variants"].map(_n).sum() == len(x), (gene, band)
            assert r["Damaging"].map(_n).sum() == (x["Functional class in the main analysis"] == "damaging").sum(), (gene, band)
            assert r["Normal"].map(_n).sum() == (x["Functional class in the main analysis"] == "normal").sum(), (gene, band)
        canon = rows[st[rows.index].str.startswith("Canonical")]
        assert canon["Variants"].map(_n).sum() == (g["In the analysis range"] == "No").sum(), gene


@need_out
def test_the_external_genes_are_their_tables(printed):
    """In the range, by band: the external genes' table (its primary label rows) gives
    the same damaging and normal counts; the new genes' table its labels by band."""
    from src import evid_external as X
    from src import evid_supp_tables as S
    e = pd.read_csv(OUT / "supplementary" / "external_genes.csv", dtype=str, keep_default_na=False)
    x = printed[(printed.Set == "External gene") & (printed["In the analysis range"] == "Yes")]
    checked = 0
    for gene in V.EXTERNAL:
        prim = S.label_name(X.PRIMARY_LABEL[gene], gene)
        rows = e[(e.Gene == gene) & (e["Label definition"] == prim)]
        for st in ("s3_10", "s11_50"):
            band = f"{S._ext_band(gene, st)} bp"
            r = rows[rows["Stratum (intronic offset, bp)"].str.replace(" bp", "") == band.replace(" bp", "")]
            if not len(r):
                continue
            g = x[(x.Gene == gene) & (x["Distance band"] == band)]
            cls = g["Functional class in the main analysis"]
            assert _n(r["Damaging"].iloc[0]) == (cls == "damaging").sum(), (gene, band)
            assert _n(r["Normal"].iloc[0]) == (cls == "normal").sum(), (gene, band)
            checked += 1
    assert checked == 14


@need_out
def test_the_new_genes_canonical_dinucleotides_are_their_counts_by_offset(printed):
    """At 1-2 bp, by gene and class: the new genes' classes by intronic offset."""
    c = pd.read_csv(OUT / "supplementary" / "new_external_genes.counts_by_offset.csv", dtype=str, keep_default_na=False)
    assert len(c)
    x = printed[(printed.Set == "External gene") & (printed["In the analysis range"] == "No")]
    for gene in ("CTCF", "SFPQ", "RAD51D", "XRCC2", "SBDS"):
        assert len(x[x.Gene == gene]) == V.COUNTS["1-2 bp"][gene]
