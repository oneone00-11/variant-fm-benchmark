"""The label source audit (stage E16c): the seven label files against the sources they
cite, read-only, and what its tables may claim.

The audit compares each label file in phase1/data/assay_labels/ with what the score set
Supplementary Table S1 cites gives under rules R1-R6 (src/evid_label_audit.py). These
tests hold it to reading only, to a row for every rule, to tables that add up, to
quotes that were checked against the papers, and to compact extracts that match the
files they were made from (those not committed are checked where the cache holds them).
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
AUDIT = PHASE1 / "data/evidence/label_audit"
REPORTS = PHASE1 / "reports/evidence"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
from src import config as C                 # noqa: E402
from src import evid_common as K            # noqa: E402
from src import evid_label_audit as A       # noqa: E402
from src import evid_label_sources as LS    # noqa: E402

GENES = list(C.ASSAY_LABELS)
RULES = ["R1", "R2", "R3", "R4", "R5", "R6"]
OUTPUTS = ["label_audit_summary.csv", "label_audit_diff.csv", "label_audit_impact.csv",
           "label_audit_unclassified.csv"]


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


@pytest.fixture(scope="module")
def summary() -> pd.DataFrame:
    return pd.read_csv(REPORTS / "label_audit_summary.csv")


@pytest.fixture(scope="module")
def diff() -> pd.DataFrame:
    return pd.read_csv(REPORTS / "label_audit_diff.csv")


def _manifests() -> list[Path]:
    return sorted(AUDIT.glob("*_raw/manifest.json"))


# --------------------------------------------------------------------------
# read-only, and reproducible
# --------------------------------------------------------------------------
def test_the_audit_reads_only_and_rewrites_its_own_files_exactly(tmp_path, monkeypatch):
    """The stage run again, its outputs sent to a scratch folder: the label files,
    config.py and the analysis set keep their sha256, and every output comes out byte
    for byte as committed."""
    if K.resolve_atlas()[1]:
        pytest.skip("the companion atlas is not here (its mapping module aligns the -a-1 versions)")
    keep = [C.ASSAY_LABEL_DIR / s["file"] for s in C.ASSAY_LABELS.values()]
    keep += [Path("src/config.py"), K.SET_PATH]
    monkeypatch.chdir(PHASE1)
    before = {p: _sha(p) for p in keep}
    for name, attr in zip(OUTPUTS, ("OUT_SUMMARY", "OUT_DIFF", "OUT_IMPACT", "OUT_UNCLASSIFIED")):
        monkeypatch.setattr(A, attr, tmp_path / name)
    monkeypatch.setattr(A, "DOC_PATH", tmp_path / "label-source-audit.md")
    A.main()
    assert {p: _sha(p) for p in keep} == before
    for name in OUTPUTS:
        assert (tmp_path / name).read_bytes() == (REPORTS / name).read_bytes(), name
    assert ((tmp_path / "label-source-audit.md").read_bytes()
            == (REPO / "docs/label-source-audit.md").read_bytes())


# --------------------------------------------------------------------------
# the tables
# --------------------------------------------------------------------------
def test_every_gene_has_a_row_for_every_rule(summary):
    for gene in GENES:
        have = set(summary.loc[summary["gene"] == gene, "rule"])
        assert set(RULES) <= have, (gene, sorted(set(RULES) - have))


def test_identical_and_different_add_up(summary):
    assert (summary["n_identical"] + summary["n_differ"] == summary["n_compared"]).all()
    live = summary[summary["n_only_ours"].notna()]
    assert (live["n_compared"] + live["n_only_ours"] == live["n_our_file"]).all()
    assert (summary.loc[summary["exact_match"] == "yes", "n_differ"] == 0).all()
    assert (summary.loc[summary["exact_match"] == "no", "n_differ"] > 0).all()


def test_the_diff_table_lists_every_difference_once(summary, diff):
    assert len(diff) == int(summary["n_differ"].sum())
    keys = ["gene", "rule", "version", "rule_detail"]
    per = diff.groupby(keys).size()
    for r in summary[summary["n_differ"] > 0].itertuples():
        assert per[(r.gene, r.rule, r.version, r.rule_detail)] == r.n_differ, r.rule_detail
    assert not diff.duplicated(keys + ["hgvs_nt"]).any()
    assert (diff["our_canonical"] != diff["source_canonical"]).all()


def test_impact_rows_only_for_genes_no_rule_reproduces(summary):
    impact = pd.read_csv(REPORTS / "label_audit_impact.csv")
    reproduced = set(summary.loc[summary["exact_match"] == "yes", "gene"])
    assert set(impact["gene"]) == set(GENES) - reproduced
    for gene, sub in impact.groupby("gene"):            # each listed gene, every stratum
        assert set(sub["stratum"]) == set(A.STRATA), gene


def test_the_impact_table_starts_from_the_counts_files_and_adds_up(monkeypatch):
    """Its 'now' columns are depth_counts.csv / set_counts.csv (which the stage checks
    against the analysis set); both sides add up within each stratum."""
    monkeypatch.chdir(PHASE1)
    impact = pd.read_csv(REPORTS / "label_audit_impact.csv")
    for r in impact.itertuples():
        now = A.current_counts(r.gene)[r.stratum]
        assert (r.n_variants, r.n_labelled_now, r.n_damaging_now, r.n_normal_now) == (
            now["n"], now["labelled"], now["damaging"], now["normal"]), (r.gene, r.stratum)
    for side in ("now", "rule"):
        lab, dmg, nor, exc = (impact[f"n_{k}_{side}"] for k in ("labelled", "damaging",
                                                                 "normal", "excluded"))
        assert (lab == dmg + nor).all() and (lab + exc == impact["n_variants"]).all(), side


# --------------------------------------------------------------------------
# what the stage found, as the author's decisions of 2026-10-21 use it
# --------------------------------------------------------------------------
def test_brca2_is_reproduced_by_its_supplementary_table(summary):
    r6 = summary[(summary["gene"] == "BRCA2") & (summary["rule"] == "R6")]
    assert len(r6) == 1 and r6.iloc[0]["exact_match"] == "yes"
    assert r6.iloc[0]["n_compared"] == r6.iloc[0]["n_our_file"] == 6959
    others = summary[(summary["gene"] != "BRCA2") & (summary["rule"] == "R6")]
    assert (others["exact_match"] == "not applicable").all()


def test_the_vhl_variants_without_a_class_are_counted_by_stratum(summary):
    unc = pd.read_csv(REPORTS / "label_audit_unclassified.csv")
    vhl = unc[unc["gene"] == "VHL"]
    r1 = summary[(summary["gene"] == "VHL") & (summary["rule"] == "R1")].iloc[0]
    assert vhl["n"].sum() == r1["n_only_source"] == 68
    assert set(vhl["stratum"]) >= set(A.UNCLASSIFIED_STRATA)
    assert vhl["n_in_analysis_set"].sum() == 0
    assert set(unc["gene"]) == {"VHL"}


def test_the_record_page_carries_the_decisions_and_the_bap1_median():
    text = (REPO / "docs/label-source-audit.md").read_text()
    for gene, decision in A.AUTHOR_DECISIONS.items():
        section = text.split(f"\n## {gene}\n", 1)[1].split("\n## ", 1)[0]
        assert decision in section, gene
    assert "-0.1260642" in A.AUTHOR_DECISIONS["BAP1"] and "c.1132del" in A.AUTHOR_DECISIONS["BAP1"]


# --------------------------------------------------------------------------
# sources: the quotes, the extracts, the cache
# --------------------------------------------------------------------------
def _numbers(text: str) -> set[float]:
    """Every number a quote writes, its own markup dropped and U+2212 read as a minus."""
    plain = re.sub(r"<[^>]+>", "", text)
    return {float(x.replace("−", "-")) for x in re.findall(r"[−-]?\d+(?:\.\d+)?", plain)}


@pytest.mark.parametrize("path", sorted(AUDIT.glob("*_raw/paper_rule.json")),
                         ids=lambda p: p.parent.name)
def test_each_paper_rule_was_checked_against_the_paper_and_uses_only_quoted_numbers(path):
    pr = json.loads(path.read_text())
    check = json.loads((path.parent / "paper_quote_check.json").read_text())
    assert check["paper_rule_sha256"] == _sha(path)
    assert check["quotes"] and all(q["found"] for q in check["quotes"])
    assert len(check["quotes"]) == len(pr["quotes"])
    quoted = set().union(*(_numbers(q["text"]) for q in pr["quotes"]))
    for rule in pr.get("rules", []):
        for c in rule["classes"]:
            for cond in c["conditions"]:
                assert cond["value"] == 0 or float(cond["value"]) in quoted, (rule["name"], cond)


def test_the_brca2_grouping_is_the_one_the_paper_writes():
    pr = json.loads((AUDIT / "brca2_raw/paper_rule.json").read_text())
    rule = pr["supplementary_rule"]
    quote = next(q["text"] for q in pr["quotes"] if q["where"] == rule["quote_where"])
    assert "combined benign" in quote and "combined pathogenic" in quote and "VUS" in quote
    cats = {r["category"] for r in LS.read_tsv(AUDIT / "brca2_raw/supplementary_table3_categories.tsv")}
    assert cats == set(rule["groups"])
    for cat, canon in rule["groups"].items():
        want = {"P": "abnormal", "B": "normal"}.get(cat.split()[0], "intermediate")
        assert canon == want, cat


@pytest.mark.parametrize("path", _manifests(), ids=lambda p: p.parent.name)
def test_the_compact_extracts_match_the_files_they_were_made_from(path):
    """Row for row. A source that is not committed is checked where the cache holds it
    (scripts/reproduce_evidence.py --fetch-inputs fills the cache and runs the same
    check, which may not skip)."""
    man = json.loads(path.read_text())
    differ, absent = LS.extract_differences(man, REPO)
    assert not differ, differ
    if absent:
        print(f"skipped: {', '.join(absent)}: its source is not in the cache")
        pytest.skip(f"skipped: source of {', '.join(absent)} not in the cache "
                    "(scripts/reproduce_evidence.py --fetch-inputs)")


@pytest.mark.parametrize("path", _manifests(), ids=lambda p: p.parent.name)
def test_every_listed_file_is_the_file_on_disk_and_the_cached_ones_stay_out_of_git(path):
    man = json.loads(path.read_text())
    cache = (PHASE1 / "data/evidence/label_audit_cache").relative_to(REPO).as_posix()
    for e in man["files"]:
        f = REPO / e["path"]
        if e.get("cached"):
            assert e["path"].startswith(cache + "/"), e["path"]
            if not f.exists():
                continue
        assert _sha(f) == e["sha256"], e["path"]
    if (REPO / ".git").exists():
        cached = [e["path"] for e in man["files"] if e.get("cached")]
        if cached:
            out = subprocess.run(["git", "-C", str(REPO), "ls-files", "--", *cached],
                                 capture_output=True, text=True, check=True).stdout
            assert out == "", out


def test_the_label_audit_sources_are_stored_byte_for_byte():
    if not (REPO / ".git").exists():
        pytest.skip("a release archive has no .git")
    paths = [e["path"] for m in _manifests() for e in json.loads(m.read_text())["files"]
             if not e.get("cached") and e["path"].startswith("phase1/data/evidence/label_audit/")]
    out = subprocess.run(["git", "check-attr", "text", "--", *paths], cwd=REPO,
                         capture_output=True, text=True, check=True).stdout.splitlines()
    assert len(out) == len(paths) and all(l.endswith(": text: unset") for l in out), out[:2]


def test_the_label_files_entered_the_repository_on_the_date_the_audit_uses():
    if not (REPO / ".git").exists():
        pytest.skip("a release archive has no .git")
    dates = set()
    for spec in C.ASSAY_LABELS.values():
        out = subprocess.run(["git", "-C", str(REPO), "log", "--diff-filter=A", "--format=%cs",
                              "--", f"phase1/data/assay_labels/{spec['file']}"],
                             capture_output=True, text=True, check=True).stdout.split()
        if not out:
            pytest.skip("the commit that added the label files is not in this clone (shallow)")
        dates.add(out[-1])
    assert dates == {A.LABEL_FILES_COMMITTED}


# --------------------------------------------------------------------------
# the entry point
# --------------------------------------------------------------------------
def _entry_point():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_entry_point_runs_the_audit_after_e16b_and_before_the_delta():
    E = _entry_point()
    rows = [(m, a) for m, _, a in E.STAGES]
    i = rows.index(("src.evid_label_audit", []))
    assert rows[i - 1] == ("src.evid_newgenes", ["--gene", "SBDS"])
    assert rows[i + 1][0] == "src.evid_delta"
    assert REPO / "docs" / "label-source-audit.md" in E.OUTPUT_ROOTS
    assert PHASE1 / "data" / "evidence" / "label_audit_cache" in E.NOT_OUTPUTS


def test_every_committed_source_of_the_audit_is_listed_as_an_input():
    E = _entry_point()
    inputs = E._run_inputs()
    need = {e["path"] for m in _manifests() for e in json.loads(m.read_text())["files"]
            if not e.get("cached")}
    need |= {str(p.relative_to(REPO)) for p in AUDIT.glob("*_raw/*") if p.is_file()}
    assert need <= inputs, sorted(need - inputs)[:5]
