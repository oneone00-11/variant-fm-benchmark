"""Main-text numbers that come from a supplementary table must be that table's
value rounded once.

The supplement prints each such quantity to more decimals than the main text, at a
precision chosen so that rounding the printed value half-up reproduces the main-text
cell (phase1/src/build_supp_tables.col_digits). Until the v2 review round nothing
checked that boundary, and the two documents could disagree: Table 1 said I^2 = 60
where the supplement printed 59.5, both from one unrounded 59.49.

Both files are needed; the tests skip when either is absent. The manuscript is found
as tests/test_reported_numbers.py finds it; the supplement through $VARIANT_FM_SUPPLEMENT
or a *supplementary*.docx beside it.
"""
from __future__ import annotations

import os
import re
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


def _newest(paths):
    paths = [p for p in paths if p.exists()]
    return max(paths, key=lambda q: q.stat().st_mtime) if paths else None


def _manuscript():
    env = os.environ.get("VARIANT_FM_MANUSCRIPT")
    if env:
        p = Path(env).expanduser()
        return p if p.exists() else None
    desk = Path.home() / "Desktop"
    return _newest([*desk.glob("calibration_draft*.docx"), *(desk / "calibration").glob("calibration_draft*.docx")])


def _supplement():
    env = os.environ.get("VARIANT_FM_SUPPLEMENT")
    if env:
        p = Path(env).expanduser()
        return p if p.exists() else None
    desk = Path.home() / "Desktop"
    return _newest([*desk.glob("*supplementary*.docx"), *(desk / "calibration").glob("*supplementary*.docx")])


def _docs():
    m, s = _manuscript(), _supplement()
    if m is None or s is None:
        pytest.skip("manuscript or supplement not available in this checkout")
    from docx import Document
    return Document(str(m)), Document(str(s))


def once(text, nd):
    """Round a printed number half-up to nd decimals; returns the printed form."""
    return str(Decimal(str(text).replace("−", "-").strip()).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))


def _table(doc, pred):
    ts = [t for t in doc.tables if pred(t)]
    assert len(ts) == 1
    return ts[0]


def _rows(t):
    return [[c.text.strip() for c in r.cells] for r in t.rows]


MAIN_TO_SUPP = {"Fusion — elastic net (M1)": "M1_enet", "Pangolin": "single:pangolin",
                "Fusion — gradient-boosted trees (M2)": "M2_gbt", "SpliceAI": "single:spliceai",
                "Fusion — equal-weight mean (M0b)": "M0b_mean", "AlphaGenome": "single:alphagenome",
                "CADD": "single:cadd", "GPN-MSA": "single:gpn_msa", "phyloP": "single:phylop",
                "phastCons": "single:phastcons", "Nucleotide Transformer": "single:nt",
                "gnomAD AF": "single:gnomad_af", "AlphaMissense": "single:alphamissense"}
CAL_KEYS = dict(MAIN_TO_SUPP, **{"Fusion — elastic net (M1)": "fusion_M1", "Fusion — equal-weight mean (M0b)": "mean_M0b"})


def test_table_1_is_the_supplement_leaderboard_rounded_once():
    main, supp = _docs()
    t1 = _table(main, lambda t: t.rows[0].cells[0].text.startswith("Predictor") and "Pooled" in t.rows[0].cells[2].text)
    s5b = _table(supp, lambda t: t.rows[0].cells[0].text.strip() == "Model / tool"
                 and t.rows[0].cells[3].text.strip().startswith("I²"))   # S16c shares the first two headers
    supp_rows = {r[0]: r for r in _rows(s5b)[1:]}
    bad = []
    for name, _, rho, ci, i2 in _rows(t1)[1:]:
        srow = supp_rows[MAIN_TO_SUPP[name]]
        if rho == "—":
            if srow[1] not in ("", "n/a"):
                bad.append(f"{name}: main prints — but the supplement has {srow[1]!r}")
            continue
        lo, hi = re.findall(r"[-+]?\d*\.?\d+", srow[2])
        want = (once(srow[1], 3), f"{once(lo, 3)} – {once(hi, 3)}", once(srow[3], 0))
        if (rho, ci, i2) != want:
            bad.append(f"{name}: main {(rho, ci, i2)} but the supplement rounds to {want}")
    assert not bad, "\n  ".join(["Table 1 cells that are not the supplement rounded once:", *bad])


def _s3b_primary(supp):
    s3b = _table(supp, lambda t: t.rows[0].cells[0].text.strip() == "Condition" and t.rows[0].cells[3].text.strip() == "ECE")
    return {r[2]: r for r in _rows(s3b)[1:] if r[0] == "Assay, +BRCA1" and r[1] == "isotonic"}


def test_table_3_is_the_supplement_calibration_table_rounded_once():
    main, supp = _docs()
    t3 = _table(main, lambda t: t.rows[0].cells[2].text.strip() == "ECE")
    supp_rows = _s3b_primary(supp)
    bad = []
    for name, _, ece, brier, frac, acc in _rows(t3)[1:]:
        srow = supp_rows[CAL_KEYS[name]]
        if ece == "—":
            if srow[3] not in ("", "n/a"):
                bad.append(f"{name}: main prints — but the supplement has {srow[3]!r}")
            continue
        want = (once(srow[3], 3), once(srow[4], 3), once(Decimal(srow[5]) * 100, 1), once(Decimal(srow[6]) * 100, 1))
        if (ece, brier, frac, acc) != want:
            bad.append(f"{name}: main {(ece, brier, frac, acc)} but the supplement rounds to {want}")
    assert not bad, "\n  ".join(["Table 3 cells that are not the supplement rounded once:", *bad])


def test_figure_2_caption_quotes_the_supplement_rounded_once():
    main, supp = _docs()
    cap = [p.text for p in main.paragraphs if p.text.startswith("Figure 2.")]
    assert len(cap) == 1
    m = re.search(r"attains ECE ([\d.]+) and Brier ([\d.]+) \(Pangolin ([\d.]+) and ([\d.]+);", cap[0])
    assert m, "the Figure 2 caption no longer quotes ECE and Brier in the expected form"
    rows = _s3b_primary(supp)
    f, p = rows["fusion_M1"], rows["single:pangolin"]
    assert (m.group(1), m.group(2), m.group(3), m.group(4)) == (once(f[3], 3), once(f[4], 4), once(p[3], 3), once(p[4], 4))


def test_section_3_4_quotes_the_supplement_tp53_table_rounded_once():
    main, supp = _docs()
    para = [p.text for p in main.paragraphs if p.text.startswith("The fusion's calibrated Brier score on held-out TP53")]
    assert len(para) == 1
    t = para[0]
    nums = {
        "brier": re.search(r": ([\d.]+) versus ([\d.]+) \(ΔBrier", t),
        "ece": re.search(r"\(ECE ([\d.]+) versus ([\d.]+)\)", t),
        "acc": re.search(r"high-confidence subset \(([\d.]+) versus ([\d.]+)\)", t),
        "share": re.search(r"\(([\d.]+)% of variants against ([\d.]+)%\)", t),
    }
    assert all(nums.values()), {k: bool(v) for k, v in nums.items()}
    s8a = _table(supp, lambda tb: tb.rows[0].cells[0].text.strip() == "Label")
    rows = {r[1]: r for r in _rows(s8a)[1:] if r[0].startswith("Control-anchored")}
    fus = [v for k, v in rows.items() if k.startswith("fusion")][0]
    best = [v for k, v in rows.items() if k.startswith("best_single")][0]
    assert nums["brier"].groups() == (once(fus[3], 3), once(best[3], 3))
    assert nums["ece"].groups() == (once(fus[2], 3), once(best[2], 3))
    assert nums["acc"].groups() == (once(fus[5], 3), once(best[5], 3))
    assert nums["share"].groups() == (once(Decimal(fus[4]) * 100, 1), once(Decimal(best[4]) * 100, 1))
