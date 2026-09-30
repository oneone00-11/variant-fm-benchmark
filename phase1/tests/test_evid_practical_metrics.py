"""E15: the classical metrics at the evidence thresholds are read off the study's own
thresholds and agree with every stage that already reports the same quantity.

The checks the stage was specified with:
  * at Walker's 0.2 cut point the published-basis SpliceAI column's sensitivity and
    specificity are E2's (walker_thresholds.csv, pooled rows);
  * every row's counts add up: called plus not called is the class total;
  * the class totals are those of the gene-pooled AUROC table the instruction calls
    fig3 data, once the genes that pool leaves out (fewer than ten variants on a side)
    are added back;
  * the pooled AUROC is the Mann-Whitney statistic computed by hand;
  * the curves are monotone;
  * a held-out row is the sum of its folds.
And the ones that tie the stage to the data and to E3 and E7: coverage and the curve
points are recomputed from the analysis set, the rebuilt fusion is E3's, every
in-sample threshold is E3's, every held-out fold reproduces E3's held-out ratio
exactly, and DDX3X and TP53 reproduce E7's counts and ratios.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import rankdata

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
EV = PHASE1 / "reports/evidence"
TABLE = EV / "practical_metrics_by_threshold.csv"
BY_GENE = EV / "practical_metrics_logo_by_gene.csv"
EXTERNAL = EV / "practical_metrics_external.csv"
ROC = EV / "practical_roc_curves.csv"
PR = EV / "practical_pr_curves.csv"
AUC = EV / "practical_curve_auc.csv"
SUMMARY = REPO / "docs/practical-metrics-summary.md"
FIGURE = EV / "figures/figure_practical"
SET = PHASE1 / "data/evidence/analysis_set_v1.parquet"
MAIN = ["spliceai_walker", "pangolin", "alphagenome", "avi"]

need = pytest.mark.skipif(not TABLE.exists(), reason="E15 not run")
need_set = pytest.mark.skipif(not (TABLE.exists() and SET.exists()),
                              reason="E15 not run, or no analysis set")


def _csv(path: Path) -> pd.DataFrame:
    """Read back exactly: pandas' default float parser can return a value one binary
    digit away from the one written."""
    return pd.read_csv(path, float_precision="round_trip")


def _same(a, b) -> bool:
    return (pd.isna(a) and pd.isna(b)) or a == b


@pytest.fixture(scope="module")
def table():
    return _csv(TABLE)


@pytest.fixture(scope="module")
def auc():
    return _csv(AUC)


@pytest.fixture(scope="module")
def curves():
    return _csv(ROC), _csv(PR)


@pytest.fixture(scope="module")
def analysis_set():
    return pd.read_parquet(SET)


def test_the_stage_runs_after_the_figures_and_before_the_delta():
    text = (REPO / "scripts/reproduce_evidence.py").read_text()
    mods = re.findall(r'\(\s*"src\.([a-z_0-9]+)"', re.search(r"STAGES\s*=\s*\[(.*?)\n\]",
                                                            text, re.S).group(1))
    i = mods.index("evid_practical_metrics")
    assert mods[i - 1] == "evid_figures" and mods[i + 1] == "evid_delta"
    assert "practical-metrics-summary.md" in text        # compared by --verify


def test_the_stage_writes_only_its_own_files():
    """No output of E1-E14 is written here: every target is named practical, and the
    module writes through exactly these three calls."""
    from src import evid_practical_metrics as P
    targets = [P.OUT_TABLE, P.OUT_BY_GENE, P.OUT_EXTERNAL, P.OUT_ROC, P.OUT_PR,
               P.OUT_AUC, P.OUT_MD, Path(P.FIGURE)]
    assert all("practical" in t.name for t in targets)
    src = (PHASE1 / "src/evid_practical_metrics.py").read_text()
    assert len(re.findall(r"\.to_csv\(|write_text\(|_save\(", src)) == 3


@need
def test_walker_cut_point_matches_e2(table):
    """E4.1: sensitivity and specificity at 0.2 are E2's. E2 writes specificity as
    1 - (false-positive rate), this stage as (normal variants below) / (normal
    variants); the two agree to the last binary digit or so, and the counts exactly."""
    w = _csv(EV / "walker_thresholds.csv")
    w = w[(w.score_column == "spliceai_walker") & (w.scope == "pooled")
          & (w.clinvar_arm == "all") & (w.status == "ok")].set_index("stratum")
    rows = table[(table.tool == "spliceai_walker") & (table.threshold_source == "published")]
    checked = 0
    for _, r in rows.iterrows():
        e = w.loc[r.stratum]
        assert (r.n_pos, r.n_neg) == (e.n_pos, e.n_neg)
        if r.threshold_kind == "walker_pp3":
            assert r.sensitivity == e.sens_pp3
            assert r.specificity == pytest.approx(e.spec_pp3, rel=1e-12, abs=1e-15)
            assert r.n_pos_above + r.n_neg_above == e.n_in_pp3_band
            assert r.coverage_labelled == e.frac_pp3
            assert _same(r.band_lr_plus, e.lr_pp3)
            assert r.tier_by_band_lr == e.tier_pp3
            post = e["posterior_pp3_prior0.1"]
        else:
            assert r.bp4_n_pos_at_or_below + r.bp4_n_neg_at_or_below == e.n_in_bp4_band
            assert _same(r.bp4_band_lr, e.lr_bp4)
            assert r.tier_by_band_lr == e.tier_bp4
            post = e["posterior_bp4_prior0.1"]
        got = r["posterior_at_prior_0.10"]
        assert (pd.isna(got) and pd.isna(post)) or got == pytest.approx(post, rel=1e-12)
        checked += 1
    assert checked == 8                       # two cut points x four strata


@need
def test_every_rows_counts_add_up(table):
    """E4.1: on each side, the variants in the band and those outside it make up the
    class. _cell counts both; the recomputation below checks them against the data."""
    ok = table[table.status == "ok"]
    pp3, bp4 = ok[ok.side == "pp3"], ok[ok.side == "bp4"]
    assert len(pp3) > 150 and len(bp4) > 40
    assert (pp3.n_pos_above + pp3.n_pos_below == pp3.n_pos).all()
    assert (pp3.n_neg_above + pp3.n_neg_below == pp3.n_neg).all()
    assert (bp4.bp4_n_pos_at_or_below + bp4.bp4_n_pos_above == bp4.n_pos).all()
    assert (bp4.bp4_n_neg_at_or_below + bp4.bp4_n_neg_above == bp4.n_neg).all()
    assert (ok.n_pos + ok.n_neg == ok.n).all()
    for path in (BY_GENE, EXTERNAL):
        x = _csv(path)
        x = x[x.status == "ok"]
        assert (x.n_pos_above + x.n_pos_below == x.n_pos).all()
        assert (x.n_neg_above + x.n_neg_below == x.n_neg).all()
    # a held-out row without counts carries no counts, rather than the stratum's
    logo = table[(table.threshold_source == "logo") & (table.status != "ok")]
    assert len(logo) and logo[["n_scored", "n", "n_pos", "n_neg"]].isna().all().all()


@need_set
def test_coverage_and_counts_recomputed_from_the_set(table, analysis_set):
    """Coverage counts every scored variant in the band, labelled or not; the rest
    counts labelled variants. Recomputed here from the analysis set for every
    in-sample and published row of the thirteen stored columns."""
    from src import evid_common as K
    rows = table[table.threshold_source.isin(["in_sample", "published"])
                 & (table.status == "ok") & (table.tool != K.FUSION)]
    for _, r in rows.iterrows():
        f = K.stratum_frame(analysis_set, r.stratum)
        s = f[r.tool].to_numpy(dtype=float)
        y = f["y_assay"].to_numpy(dtype=float)
        band = (s >= r.threshold) if r.side == "pp3" else (s <= r.threshold)
        lab = np.isfinite(y) & np.isfinite(s)
        assert r.n_scored == np.isfinite(s).sum() and r.n == lab.sum()
        assert r.n_scored_in_band == band.sum()
        assert r.coverage == band.sum() / np.isfinite(s).sum()
        assert r.coverage_labelled == (band & lab).sum() / lab.sum()
        k_pos, k_neg = int((band & (y == 1)).sum()), int((band & (y == 0)).sum())
        if r.side == "pp3":
            assert (r.n_pos_above, r.n_neg_above) == (k_pos, k_neg)
        else:
            assert (r.bp4_n_pos_at_or_below, r.bp4_n_neg_at_or_below) == (k_pos, k_neg)
    assert len(rows) > 120


@need_set
def test_class_totals_are_the_gene_pooled_tables_plus_the_genes_it_leaves_out(
        table, auc, analysis_set):
    """E4.1, as it can hold. The gene-pooled AUROC (fig3 data) drops a gene with fewer
    than ten variants on a side, and its n_pos and n_neg count the genes it keeps (at
    11-50 bp BRCA1, three damaging variants). Here nothing is dropped. So the totals
    are equal wherever that pool keeps every gene, and differ by exactly the dropped
    genes' counts, recomputed from the analysis set, where it does not."""
    from src import evid_common as K
    pooled = _csv(EV / "fig_data/fig3_evidence_by_territory.csv")
    pooled = pooled[pooled.arm == "all"].set_index(["tool", "stratum"])
    ins = table[(table.threshold_source == "in_sample") & (table.threshold_kind == "supporting")]
    direct = 0
    for r in auc.itertuples():
        t = ins[(ins.tool == r.tool) & (ins.stratum == r.stratum)].iloc[0]
        assert (t.n_pos, t.n_neg) == (r.n_pos, r.n_neg)
        f = pooled.loc[(r.tool, r.stratum)]
        cols = ["y_assay"] if r.tool == K.FUSION else ["y_assay", r.tool]
        sub = K.stratum_frame(analysis_set, r.stratum).dropna(subset=cols)
        per = sub.groupby("gene", observed=True).y_assay.agg(
            pos=lambda v: int((v == 1).sum()), neg=lambda v: int((v == 0).sum()))
        out = per[(per.pos < K.MIN_POS) | (per.neg < K.MIN_NEG)]
        assert "; ".join(sorted(map(str, out.index))) == (
            "" if pd.isna(r.genes_not_in_gene_pool) else r.genes_not_in_gene_pool)
        assert r.n_pos == f.n_pos + int(out.pos.sum())
        assert r.n_neg == f.n_neg + int(out.neg.sum())
        assert (r.n_pos_gene_pooled, r.n_neg_gene_pooled) == (f.n_pos, f.n_neg)
        if not len(out):
            assert (r.n_pos, r.n_neg) == (f.n_pos, f.n_neg)
            direct += 1
    assert direct >= 28                       # 3-10 and 3-50 bp, all fourteen columns


@need_set
def test_pooled_auroc_is_the_mann_whitney_statistic(auc, analysis_set):
    """The fusion is rebuilt, not stored; its curves are tied to E3's below."""
    from src import evid_common as K
    checked = 0
    for r in auc[auc.tool != K.FUSION].itertuples():
        sub = K.stratum_frame(analysis_set, r.stratum).dropna(subset=["y_assay", r.tool])
        y = sub.y_assay.to_numpy(dtype=float)
        s = sub[r.tool].to_numpy(dtype=float)
        rk = rankdata(s)
        n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
        u = rk[y == 1].sum() - n_pos * (n_pos + 1) / 2
        assert r.auroc_pooled_curve == pytest.approx(u / (n_pos * n_neg), rel=1e-12)
        # average precision, the step-wise sum over distinct scores
        order = np.argsort(-s, kind="mergesort")
        ss, yy = s[order], y[order]
        last = np.append(np.flatnonzero(np.diff(ss)), len(ss) - 1)
        tp, fp = np.cumsum(yy == 1)[last], np.cumsum(yy == 0)[last]
        prec, rec = tp / (tp + fp), tp / n_pos
        ap = float(np.sum(np.diff(np.concatenate([[0.0], rec])) * prec))
        assert r.auprc_pooled_curve == pytest.approx(ap, rel=1e-12)
        checked += 1
    assert checked == 39


@need_set
def test_the_same_gene_and_cross_gene_split_by_brute_force(auc, analysis_set):
    """Every damaging-normal pair counted directly, a win as 1 and a tie as 1/2, for
    the four main columns at 11-50 bp, where the split matters most."""
    from src import evid_common as K
    for tool in MAIN:
        r = auc[(auc.tool == tool) & (auc.stratum == "s11_50")].iloc[0]
        sub = K.stratum_frame(analysis_set, "s11_50").dropna(subset=["y_assay", tool])
        s = sub[tool].to_numpy(dtype=float)
        y = sub.y_assay.to_numpy(dtype=float)
        g = sub.gene.astype(str).to_numpy()
        pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
        diff = s[pos][:, None] - s[neg][None, :]
        wins = (diff > 0) + 0.5 * (diff == 0)
        same = g[pos][:, None] == g[neg][None, :]
        assert r.auroc_same_gene_pairs == pytest.approx(wins[same].mean(), rel=1e-12)
        assert r.auroc_cross_gene_pairs == pytest.approx(wins[~same].mean(), rel=1e-12)
        assert r.share_cross_gene_pairs == pytest.approx((~same).mean(), rel=1e-12)
        assert r.auroc_pooled_curve == pytest.approx(wins.mean(), rel=1e-12)


@need
def test_curves_are_monotone_and_complete(curves, auc):
    roc, pr = curves
    for (tool, st), c in roc.groupby(["tool", "stratum"], sort=False):
        thr = c.threshold.to_numpy()
        assert np.isinf(thr[0]) and (np.diff(thr[1:]) < 0).all()
        assert (np.diff(c.fpr.to_numpy()) >= 0).all() and (np.diff(c.tpr.to_numpy()) >= 0).all()
        assert (c.fpr.iloc[0], c.tpr.iloc[0]) == (0.0, 0.0)
        assert (c.fpr.iloc[-1], c.tpr.iloc[-1]) == (1.0, 1.0)
        p = pr[(pr.tool == tool) & (pr.stratum == st)]
        assert len(p) == len(c) and np.array_equal(p.threshold.to_numpy(), thr)
        assert (np.diff(p.recall.to_numpy()) >= 0).all()
        assert (p.recall.iloc[0], p.precision.iloc[0]) == (0.0, 1.0)
        a = auc[(auc.tool == tool) & (auc.stratum == st)].iloc[0]
        assert p.recall.iloc[-1] == 1.0
        assert p.precision.iloc[-1] == pytest.approx(a.prevalence, rel=1e-12)
    assert roc.groupby(["tool", "stratum"]).ngroups == len(auc) == 42


@need_set
def test_curve_points_belong_to_their_thresholds(curves, analysis_set):
    """Each row's rates are those of 'score >= its threshold', recomputed from the data
    at every 40th row of every stored curve and at its last row."""
    from src import evid_common as K
    roc, pr = curves
    checked = 0
    for (tool, st), c in roc[roc.tool != K.FUSION].groupby(["tool", "stratum"], sort=False):
        sub = K.stratum_frame(analysis_set, st).dropna(subset=["y_assay", tool])
        s = sub[tool].to_numpy(dtype=float)
        y = sub.y_assay.to_numpy(dtype=float)
        n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
        assert len(c) == len(np.unique(s)) + 1
        p = pr[(pr.tool == tool) & (pr.stratum == st)]
        for i in sorted(set(range(1, len(c), 40)) | {len(c) - 1}):
            t = c.threshold.iloc[i]
            tp, fp = int(((s >= t) & (y == 1)).sum()), int(((s >= t) & (y == 0)).sum())
            assert (c.tpr.iloc[i], c.fpr.iloc[i]) == (tp / n_pos, fp / n_neg)
            assert (p.recall.iloc[i], p.precision.iloc[i]) == (tp / n_pos, tp / (tp + fp))
            checked += 1
    assert checked > 3000


@need
def test_the_rebuilt_fusion_is_e3s(curves):
    """The fusion is rebuilt here, not stored: its distinct labelled scores in each band
    are exactly the grid E3 wrote for the same column."""
    roc, _ = curves
    for st in ("s3_10", "s11_50", "s3_50"):
        grid = _csv(EV / f"interval_lr_fusion_enet_{st}.csv")["score"].to_numpy()
        mine = roc[(roc.tool == "fusion_enet") & (roc.stratum == st)].threshold.to_numpy()[1:]
        assert np.array_equal(np.sort(mine), np.sort(grid))


@need
def test_a_held_out_row_is_the_sum_of_its_folds(table):
    """E4.2: counts added over the folds whose threshold reached the tier, never a
    mean of proportions; the combined Atlas score without BRCA1 and RAD51C."""
    g = _csv(BY_GENE)
    logo = table[(table.threshold_source == "logo") & (table.side == "pp3")]
    checked = 0
    for r in logo.itertuples():
        f = g[(g.tool == r.tool) & (g.stratum == r.stratum) & (g.tier == r.threshold_kind)]
        used = f[f.counted_in_logo_row.astype(bool)]
        if r.status == "not evaluable":
            continue
        assert r.folds_total == len(f) and r.folds_used == len(used)
        assert r.folds_with_band_lr == used.band_lr_plus.notna().sum()
        if r.status != "ok":
            assert len(used) == 0
            continue
        for col in ("n_scored", "n", "n_pos", "n_neg", "n_pos_above", "n_neg_above",
                    "n_pos_below", "n_neg_below", "n_scored_in_band"):
            assert getattr(r, col) == used[col].sum(), (r.tool, r.stratum, col)
        assert r.sensitivity == r.n_pos_above / r.n_pos
        assert r.coverage == r.n_scored_in_band / r.n_scored
        assert r.threshold_fold_min == used.threshold_from_6_genes.min()
        assert r.threshold_fold_max == used.threshold_from_6_genes.max()
        checked += 1
    assert checked > 90
    avi = g[(g.tool == "avi") & g.heldout_gene.isin(["BRCA1", "RAD51C"])]
    assert len(avi) and not avi.counted_in_logo_row.astype(bool).any()


@need
def test_every_held_out_fold_reproduces_e3s_ratio():
    """Each fold's threshold, applied to its held-out gene here, gives the held-out
    ratio E3 wrote for that fold, exactly: same gene, same band, same rule."""
    g = _csv(BY_GENE)
    g = g[g.status == "ok"]
    assert len(g) > 500
    for r in g.itertuples():
        assert _same(r.band_lr_plus, r.e3_heldout_lr), (r.tool, r.stratum, r.tier,
                                                        r.heldout_gene)
    folds = _csv(EV / "evidence_thresholds_logo_folds.csv")
    folds = folds[(folds.side == "pp3") & folds.stratum.isin(g.stratum.unique())]
    merged = g.merge(folds, on=["tool", "stratum", "tier", "heldout_gene"])
    assert len(merged) == len(g)
    assert (merged.threshold_from_6_genes_x == merged.threshold_from_6_genes_y).all()
    assert all(_same(a, b) for a, b in zip(merged.e3_heldout_lr, merged.heldout_lr))


@need
def test_in_sample_thresholds_are_e3s(table):
    ev = _csv(EV / "evidence_thresholds.csv").set_index(["tool", "stratum", "tier"])
    ins = table[table.threshold_source == "in_sample"]
    for r in ins.itertuples():
        tier = r.threshold_kind.removeprefix("bp4_")
        e = ev.loc[(r.tool, r.stratum, tier)]
        tau = e[f"{r.side}_threshold_insample"]
        if e.status != "ok":
            assert r.status == "not evaluable"
        elif pd.isna(tau):
            assert r.status == "not reached" and pd.isna(r.threshold)
        else:
            assert r.status == "ok" and r.threshold == tau
    # the BP4 side has no per-fold thresholds in E3, and its held-out rows say so
    bp4 = table[(table.threshold_source == "logo") & (table.side == "bp4")]
    assert len(bp4) and (bp4.status == "no held-out threshold").all()


@need
def test_external_genes_reproduce_e7():
    """E4.3: at the thresholds E7 carried over, the counts and the band ratio are
    E7's; at the published cut point, E7's sensitivity and specificity. Where E7
    carried nothing, the row takes E7's words and applies no threshold."""
    x = _csv(EXTERNAL)
    checked = 0
    for gene, label in (("DDX3X", "y_deposit"), ("TP53", "y_control_anchored")):
        e = _csv(EV / f"external_{gene.lower()}.csv")
        e = e[e.label_definition == label].set_index(["tool", "stratum"])
        for r in x[x.gene == gene].itertuples():
            f = e.loc[(r.tool, r.stratum)]
            if r.status != "ok":
                assert pd.isna(r.threshold) and pd.isna(r.sensitivity)
                if f.status != "ok":
                    assert (r.status, r.note) == (f.status, f.reason)
                elif r.threshold_kind != "walker_pp3" and bool(f.column_basis_comparable):
                    k = r.threshold_kind
                    assert (r.status, r.note) == (f[f"e3_{k}_tier_here"], f[f"e3_{k}_note"])
                continue
            if r.threshold_kind == "walker_pp3":
                assert r.sensitivity == f.walker_sens_pp3
                assert r.specificity == f.walker_spec_pp3
                assert r.n_pos_above == f.walker_pp3_n_pos_band
            else:
                k = r.threshold_kind
                assert r.threshold == f[f"e3_{k}_threshold"]
                assert r.n_pos_above == f[f"e3_{k}_n_pos_band"]
                assert r.n_neg_above == f[f"e3_{k}_n_neg_band"]
                assert _same(r.band_lr_plus, f[f"e3_{k}_lr_here"])
            checked += 1
    assert checked > 100


@need
def test_the_summary_is_written_from_the_tables(table, auc, curves):
    """The one-page summary is the stage's own output from the committed tables, not a
    hand edit; its cells carry the tables' values and its sentences at most two
    numbers each."""
    from src import evid_practical_metrics as P
    ext = _csv(EXTERNAL)
    roc, _ = curves
    xmax, _ = P._inset_range(roc, P._tier_marks(table))
    md = SUMMARY.read_text(encoding="utf-8")
    assert md == P.summary(table, auc, ext, roc, xmax)
    row = next(line for line in md.splitlines() if line.startswith("| Pangolin | 3–10 bp |"))
    cells = [c.strip() for c in row.strip("|").split("|")]

    def pick(source):
        return table[(table.tool == "pangolin") & (table.stratum == "s3_10")
                     & (table.threshold_kind == "strong")
                     & (table.threshold_source == source)].iloc[0]
    assert cells[8] == f"{pick('logo').sensitivity:.2f} ({pick('in_sample').sensitivity:.2f})"
    a = auc[(auc.tool == "pangolin") & (auc.stratum == "s3_10")].iloc[0]
    assert cells[-2] == f"{a.auroc_pooled_curve:.2f}"
    sentences = re.findall(r"^\d\. (.+)$", md, re.M)
    assert len(sentences) == 5
    for sentence in sentences:
        bare = re.sub(r"\d+–\d+ bp", "", sentence)
        numbers = re.findall(r"(?<![\w–-])\d+(?:\.\d+)?%?(?![\w–])", bare)
        assert len(numbers) <= 2, sentence


def test_summary_sentences_follow_the_data():
    """The sentence helpers on made-up tables: ties, no common leader, no gene."""
    from src import evid_practical_metrics as P
    assert P._direction(pd.Series([], dtype=float)) is None
    assert P._direction(pd.Series([0.0, 0.0])) == "is unchanged"
    assert P._direction(pd.Series([0.1, -0.1])) == "moves both ways"
    col = "sensitivity_at_specificity_0.95"
    tied = [{"tool": t, "stratum": st, col: v}
            for st in P.CURVE_STRATA
            for t, v in zip(P.MAIN_TOOLS, (0.6, 0.6, 0.5, 0.4))]
    assert "SpliceAI and Pangolin are tied" in P._sentence_matched(pd.DataFrame(tied))
    rotating = [{"tool": t, "stratum": st, col: 0.7 if t == lead else 0.5}
                for st, lead in zip(P.CURVE_STRATA, P.MAIN_TOOLS)
                for t in P.MAIN_TOOLS]
    assert P._sentence_matched(pd.DataFrame(rotating)).startswith(
        "At a specificity of 0.95 no tool is the most sensitive in every band")
    empty = pd.DataFrame(columns=["status", "external_test", "main_table", "stratum",
                                  "threshold_kind", "sensitivity_minus_seven_gene",
                                  "specificity_minus_seven_gene", "gene", "tool"])
    assert P._sentence_external(empty).startswith("No threshold carried to DDX3X")


def test_wilson_and_band_ratio_from_counts():
    """The interval against its closed form at the edges and a textbook value, and the
    count-based band ratio against evid_common.band_lr on the same data."""
    from src import evid_common as K
    from src import evid_practical_metrics as P
    z2 = P.Z95 ** 2
    lo, hi = P._wilson(0, 10)
    assert lo == 0.0 and hi == pytest.approx(z2 / (10 + z2), rel=1e-12)
    lo, hi = P._wilson(10, 10)
    assert hi == 1.0 and lo == pytest.approx(10 / (10 + z2), rel=1e-12)
    lo, hi = P._wilson(30, 100)
    assert lo == pytest.approx(0.2189, abs=5e-5) and hi == pytest.approx(0.3958, abs=5e-5)
    rng = np.random.default_rng(7)
    for _ in range(200):
        n = int(rng.integers(15, 300))
        y = (rng.random(n) < rng.uniform(0.05, 0.6)).astype(float)
        s = np.round(rng.random(n) + 0.4 * y, 2)
        thr = float(rng.choice(s))
        pos, neg = s[y == 1], s[y == 0]
        for side in ("upper", "lower"):
            band = (s >= thr) if side == "upper" else (s <= thr)
            k_pos, k_neg = int((band & (y == 1)).sum()), int((band & (y == 0)).sum())
            assert _same(P._band_lr(k_pos, len(pos), k_neg, len(neg)),
                         K.band_lr(y, s, thr, side))
        k_pos, k_neg = int(((s <= thr) & (y == 1)).sum()), int(((s <= thr) & (y == 0)).sum())
        assert _same(P._benign_bound(k_pos, len(pos), k_neg, len(neg)),
                     K.band_lr_benign_bound(y, s, thr))


@need
def test_the_figure_carries_no_timestamp():
    pdf = (FIGURE.with_suffix(".pdf")).read_bytes()
    assert b"/CreationDate" not in pdf and b"/ModDate" not in pdf
    png = (FIGURE.with_suffix(".png")).read_bytes()
    assert b"Software" not in png[:4096]
