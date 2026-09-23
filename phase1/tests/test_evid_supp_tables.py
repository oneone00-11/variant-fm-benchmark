"""Guardrails on the supplementary tables.

The supplement is the place a number is most often misprinted without anyone
noticing, because by then nobody reruns anything: a table is regenerated from one
set of outputs while the figure beside it came from another, or a value stored at
four decimals is rounded again to three. These tests check the printed tables
against the files they claim to come from, not against themselves.

  * a printed cell reproduces the upstream value it was rounded from;
  * counts add up, and the grand total is the analysis set;
  * S10 never carries a published-cut-point ratio for a column the cut point was
    not defined on;
  * the rounding rule is the documented one (half-up on the decimal digits), and a
    rerun is byte-identical, both to itself and to the tables on disk.
"""
from __future__ import annotations

import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import config as C                    # noqa: E402
from src import evid_supp_tables as S          # noqa: E402

REPORTS = PHASE1 / "reports/evidence"
SUPP = REPORTS / "supplement"

needs_atlas = pytest.mark.skipif(
    not (S.ATLAS_COMPOSITION.exists() and S.ATLAS_RESOURCES.exists()),
    reason="companion atlas checkout not available (EVID_ATLAS_REPO)")


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(SUPP / name, dtype=str, keep_default_na=False)


# ---------------------------------------------------------------------------
# rounding
# ---------------------------------------------------------------------------
def test_rounding_is_half_up_on_the_decimal_digits():
    # format() and round() work on the binary value, which for these sits just
    # under the tie, and print 2.67 and 0.004; a reader rounding the stored digits
    # by hand writes 2.68 and 0.005
    assert format(2.675, ".3g") == "2.67" and S.fmt_sig(2.675, 3) == "2.68"
    assert f"{round(0.0045, 3):.3f}" == "0.004" and S.fmt_dp(0.0045, 3) == "0.005"
    assert S.fmt_dp(0.4005, 3) == "0.401"
    # significant figures keep their trailing zeros, carry over a power of ten, and
    # never fall into scientific notation
    assert S.fmt_sig(1.0, 3) == "1.00"
    assert S.fmt_sig(9.996, 3) == "10.0"
    assert S.fmt_sig(1124.3, 3) == "1120"
    assert S.fmt_sig(0.000123456, 3) == "0.000123"
    assert S.fmt_sig(0.0, 3) == "0" and S.fmt_sig(float("nan"), 3) == ""


# ---------------------------------------------------------------------------
# printed values against their sources
# ---------------------------------------------------------------------------
def test_s5_in_scope_pp3_ratio_is_the_walker_table_value():
    w = pd.read_csv(REPORTS / "walker_thresholds.csv")
    src = w[(w.score_column == "spliceai_walker") & (w.scope == "pooled")
            & (w.stratum == "s3_50") & (w.clinvar_arm == "all")].iloc[0]
    s5 = _read("tableS5_published_cut_points.csv")
    row = s5[(s5["Stratum (intronic offset, bp)"] == "3–50")
             & (s5["ClinVar record status"] == "All variants")
             & (s5["Gene"] == "All genes (pooled)")]
    assert len(row) == 1
    cell = row.filter(like="PP3 band").iloc[0, 0]
    lead, ci = cell.split(" ", 1)
    lo, hi = (float(x) for x in ci.strip("()").split("–"))
    # within half a unit of the third significant figure of the stored value
    for printed, stored in ((float(lead), src.lr_pp3), (lo, src.lr_pp3_lo),
                            (hi, src.lr_pp3_hi)):
        unit = 10 ** (np.floor(np.log10(abs(stored))) - 2)
        assert abs(printed - stored) <= unit / 2 + 1e-12
    assert row["PP3 tier"].iloc[0] == src.tier_pp3


def test_s2_total_is_the_analysis_set_and_rows_add_up():
    s2 = _read("tableS2_analysis_set.csv")
    n_set = len(pd.read_parquet(PHASE1 / "data/evidence/analysis_set_v1.parquet",
                                columns=["variant_id"]))
    assert n_set == 8853
    last = s2.iloc[-1]
    assert (last["Gene"], last["ClinVar record status"]) == ("All genes", "All variants")
    assert int(last["Variants"]) == n_set
    num = s2[["Variants", "Labelled", "Damaging", "Normal"]].astype(int)
    assert (num["Labelled"] == num["Damaging"] + num["Normal"]).all()
    assert (num["Labelled"] <= num["Variants"]).all()
    # the three stratum totals add up to the overall total
    tot = s2[(s2.Gene == "All genes") & (s2["ClinVar record status"] == "All variants")]
    per = tot[tot["Stratum (intronic offset, bp)"].isin(
        [S.STRATUM_LABEL[k] for k in ("pm12", "s3_10", "s11_50")])]
    assert per["Variants"].astype(int).sum() == n_set
    # and the gene rows of each stratum add up to that stratum's total
    for st, t in zip(per["Stratum (intronic offset, bp)"], per["Variants"].astype(int)):
        genes = s2[(s2["Stratum (intronic offset, bp)"] == st) & (s2.Gene != "All genes")]
        assert genes["Variants"].astype(int).sum() == t


def test_s10_published_cut_point_only_for_the_spliceai_columns():
    s10 = _read("tableS10_clinvar_status.csv")
    spliceai = {S.tool_labels()[t] for t in S.WALKER_TOOLS}
    pub = [c for c in s10.columns if c.startswith("Band LR at published cut point")]
    assert len(pub) == 1
    other = s10[~s10["Tool"].isin(spliceai)]
    assert (other[pub[0]] == "").all()
    # and the SpliceAI rows carry it wherever a band is defined
    band = s10[s10["Tool"].isin(spliceai) & (s10["Damaging"] != "")]
    assert (band[pub[0]] != "").all()
    # nothing from territory_metrics' walker_* columns: every published-cut value
    # is one of arms_at_tool_threshold's
    a = pd.read_csv(REPORTS / "arms_at_tool_threshold.csv")
    allowed = {S.fmt_lr(v) for v in a.loc[a.threshold_basis == "published cut point",
                                          "band_lr"]} | {"not evaluable"}
    assert set(band[pub[0]]) <= allowed
    assert not any("walker" in c.lower() for c in s10.columns)


def test_s6_per_fold_ratios_rebuild_the_fold_file():
    s6 = _read("tableS6_thresholds_pathogenic.csv")
    folds = pd.read_csv(REPORTS / "evidence_thresholds_logo_folds.csv")
    labels = {v: k for k, v in S.tool_labels().items()}
    strata = {v: k for k, v in S.STRATUM_LABEL.items()}
    checked = 0
    for _, r in s6[s6["Held-out LR per fold"] != ""].iterrows():
        tool = labels[r["Tool"]]
        st, tier = strata[r["Stratum (intronic offset, bp)"]], r["Tier"].lower()
        f = folds[(folds.tool == tool) & (folds.stratum == st) & (folds.tier == tier)
                  & (folds.side == "pp3")].set_index("heldout_gene")
        parts = dict(p.split("=") for p in r["Held-out LR per fold"].split("; "))
        assert set(parts) == set(f.index)
        assert int(r["Held-out genes"]) == len(f)
        cut = float(r["Tier LR cut"])
        assert int(r["Held-out genes clearing the cut"]) == int(
            (f.heldout_lr >= cut).sum())
        for g, v in parts.items():
            if v in ("no threshold", "n.e."):
                assert pd.isna(f.loc[g, "heldout_lr"])
            else:
                assert v == S.fmt_lr(f.loc[g, "heldout_lr"])
                checked += 1
    assert checked > 100


def test_s13_shares_come_from_the_counts():
    s13 = _read("tableS13_inframe_attribution.csv")
    n = s13[["In-frame", "Out-of-frame", "Undetermined", "Total"]].astype(int)
    assert (n["In-frame"] + n["Out-of-frame"] + n["Undetermined"] == n["Total"]).all()
    for (_, r), (_, c) in zip(s13.iterrows(), n.iterrows()):
        assert r["In-frame share"] == S.fmt_dp(c["In-frame"] / c["Total"])
        res = c["In-frame"] + c["Out-of-frame"]
        assert r["In-frame share among resolved"] == (
            S.fmt_dp(c["In-frame"] / res) if res else "")


def test_s1_counts_and_publication_names():
    s1 = _read("tableS1_functional_standards.csv").set_index("Gene")
    sc = pd.read_csv(REPORTS / "set_counts.csv").groupby("gene")[
        ["n", "n_labelled", "n_damaging", "n_normal"]].sum()
    for g in S.SEVEN:
        r = s1.loc[g]
        assert (int(r["Intron-side SNVs"]), int(r["Labelled"]), int(r["Damaging"]),
                int(r["Normal"])) == tuple(int(x) for x in sc.loc[g])
        # where the label file's name carries author and year, the printed short
        # name must agree with it
        m = re.fullmatch(r"[a-z0-9]+_([a-z]+)(\d{4})\.tsv", C.ASSAY_LABELS[g]["file"])
        if m:
            # a preprint carries a "(preprint)" suffix after author and year
            assert r["Assay publication"].split(" (")[0] == \
                f"{m.group(1).capitalize()} {m.group(2)}"
    assert set(s1.index) == set(S.SEVEN) | {"DDX3X", "TP53"}
    ext = pd.read_csv(REPORTS / "external_tp53.csv")
    t = ext[(ext.label_definition == "y_control_anchored") & (ext.stratum == "all_1_50")]
    assert int(s1.loc["TP53", "Damaging"]) == int(t.n_pos.iloc[0])
    assert int(s1.loc["TP53", "Normal"]) == int(t.n_neg.iloc[0])


def test_s3_orientation_and_coverage():
    s3 = _read("tableS3_predictor_columns.csv")
    assert list(s3["Score column"]) == S.TOOLS and len(s3) == 13
    fo = pd.read_csv(REPORTS / "feature_orientation.csv").set_index("feature")
    for _, r in s3.iterrows():
        printed = r["Spearman rho with functional pathogenicity"]
        assert printed == S.fmt_dp(fo.loc[r["Score column"], "spearman_vs_pathogenicity"])
        assert float(printed) > 0


def test_s12_external_rows_follow_the_basis_check():
    s12 = _read("tableS12_external_genes.csv")
    thr = [c for c in s12.columns if c.endswith("(median of leave-one-gene-out fits)")]
    assert len(thr) == 3
    no = s12[s12["Same variable as the fitted column"] == "No"]
    assert len(no) and (no["Why not comparable"] != "").all()
    assert (no[thr] == "").all().all()
    # the Walker-basis ratio at the published cut point is the external table's
    e = pd.read_csv(REPORTS / "external_ddx3x.csv")
    e = e[(e.tool == "spliceai_walker") & (e.stratum == "s3_50")
          & (e.label_definition == "y_control_anchored")].iloc[0]
    row = s12[(s12.Gene == "DDX3X") & (s12["Label definition"] == "Control-anchored")
              & (s12["Stratum (intronic offset, bp)"] == S._ext_band("DDX3X", "s3_50"))
              & (s12.Tool == S.tool_labels()["spliceai_walker"])]
    col = [c for c in s12.columns if c.startswith("LR at published cut point")][0]
    # the cell is the point estimate followed by its interval
    assert row[col].iloc[0].split(" (")[0] == S.fmt_lr(e.walker_lr_pp3)
    assert row[col].iloc[0].endswith(")")


# ---------------------------------------------------------------------------
# the manifest, determinism and the house rules
# ---------------------------------------------------------------------------
def test_every_manifest_table_exists_and_matches_its_record():
    man = _read(S.MANIFEST)
    # every table this stage writes, plus S9, which its own stage writes
    assert set(t[1] for t in S.TABLES) <= set(man["File"])
    assert "tableS9_fusion.csv" in set(man["File"])
    for _, r in man.iterrows():
        p = SUPP / r["File"]
        assert p.exists() and p.stat().st_size > 0
        df = pd.read_csv(p, dtype=str, keep_default_na=False)
        assert len(df) == int(r["Rows"]) > 0
        assert hashlib.sha256(p.read_bytes()).hexdigest() == r["File sha256"]
        # sources are named so that another machine can resolve them, and each
        # listed file has one hash entry, paired by position; only the data inputs
        # are hashed, since a code edit is not a change of input
        srcs = [x for x in r["Source files"].split("; ") if x]
        assert not any(x.startswith("/") for x in srcs)
        shas = [x for x in r["Source sha256"].split("; ") if x]
        if r["File"] == "tableS9_fusion.csv":
            continue
        assert len(srcs) == len(shas), r["File"]
        for src, sha in zip(srcs, shas):
            if src.endswith(".py"):
                assert sha == "n/a (code)", (r["File"], src)
            else:
                assert len(sha) == 64 and int(sha, 16) >= 0, (r["File"], src)


@needs_atlas
def test_rerun_is_byte_identical_and_matches_the_tables_on_disk(tmp_path):
    """Two separate processes, so nothing cached in this one is shared."""
    a, b = tmp_path / "a", tmp_path / "b"
    env = dict(os.environ, PYTHONPATH=str(PHASE1))
    for out in (a, b):
        subprocess.run([sys.executable, "-W", "ignore", "-m", "src.evid_supp_tables",
                        "--out-dir", str(out)], cwd=PHASE1, env=env, check=True,
                       capture_output=True)
    names = sorted(p.name for p in a.iterdir())
    assert names == sorted(p.name for p in b.iterdir())
    for n in names:
        assert (a / n).read_bytes() == (b / n).read_bytes(), n
        # a table on disk older than its inputs fails here
        assert (a / n).read_bytes() == (SUPP / n).read_bytes(), f"{n} is stale"



def test_s8_counts_are_the_dilution_summary():
    """S8 is the dilution summary in print form: every count it prints is the
    summary's."""
    summ = pd.read_csv(REPORTS / "threshold_dilution_summary.csv")
    s8 = pd.read_csv(SUPP / "tableS8_dilution.csv", dtype=str)
    assert len(s8) == len(summ)
    reached = s8["Subsets reaching the tier"].str.split(" ").str[0].astype(int)
    assert reached.sum() == int(summ.n_subsets_reached.sum())


def test_printed_thresholds_select_the_same_variants():
    """A printed threshold must select exactly the variants the fitted value
    selects. At four significant figures AlphaGenome's 2.19988 prints as 2.200,
    above every score in the set, which would empty the band."""
    ev = pd.read_csv(REPORTS / "evidence_thresholds.csv")
    ev = ev[(ev.status == "ok") & ev.pp3_threshold_reachable.astype(bool)
            & ev.tool.isin(S.TOOLS)]
    df = pd.read_parquet(PHASE1 / "data/evidence/analysis_set_v1.parquet")
    checked = 0
    for r in ev.itertuples():
        s = df[r.tool].dropna().to_numpy(dtype=float)
        printed = float(S.fmt_thr(r.pp3_threshold_insample, r.tool))
        assert (s >= printed).sum() == (s >= r.pp3_threshold_insample).sum(), (r.tool, r.stratum, r.tier)
        checked += 1
    assert checked > 50
