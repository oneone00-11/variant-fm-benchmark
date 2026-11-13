"""The check of the variant list's HGVS descriptions with VariantValidator
(src/evid_hgvs_check.py, E9.4, stage 42): its four checks on a small list made here and
on altered copies of it, each alteration caught by the check that is for it; then the
committed files: the inputs hold together, the two files are what the stage writes from
them, and the stage is registered where the layout's are. Nothing here asks a service."""
from __future__ import annotations

import copy
import hashlib
import json
import re
import sys
from pathlib import Path

import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
from src import evid_hgvs_check as H  # noqa: E402

OUT = PHASE1 / H.OUT
need_out = pytest.mark.skipif(not (OUT / H.FILES[0]).exists(), reason="stage 42 has not run")
need_excerpt = pytest.mark.skipif(not (PHASE1 / H.EXCERPT).exists(), reason="the excerpt is not here")


# ----------------------------------------------------------------- a small list, made here
def _list() -> list[dict]:
    rows = []
    for gene, nm, enst, nc, c, pos, ref, alt in (
            ("BAP1", "NM_004656.4", "ENST00000460680.6", "NC_000003.12", "c.37+3G>A", 52409000, "C", "T"),
            ("BAP1", "NM_004656.4", "ENST00000460680.6", "NC_000003.12", "c.37+3G>C", 52409000, "C", "G"),
            ("BARD1", "NM_000465.4", "ENST00000260947.9", "NC_000002.12", "c.2002-3C>G", 214729011, "G", "C"),
            ("SBDS", "NM_016038.4", "ENST00000246868.7", "NC_000007.14", "c.129-3C>T", 66994000, "G", "A"),
            ("SBDS", "NM_016038.4", "ENST00000246868.7", "NC_000007.14", "c.129-4T>C", 66994001, "A", "G")):
        rows.append({"Set": "External gene" if gene == "SBDS" else "Analysis set", "Gene": gene,
                     "RefSeq transcript (MANE Select)": nm, "HGVS c.": c, "Ensembl transcript (MANE Select)": enst,
                     "HGVS g.": f"{nc}:g.{pos}{ref}>{alt}"})
    return rows


def _entry(description: str, interface: str, g: str, warnings=()) -> dict:
    return {"submitted_description": description, "interface": interface, "transcript_description": description,
            "genomic_description_grch38": g, "vcf_chr": "1", "vcf_pos": "1", "vcf_ref": "A", "vcf_alt": "C",
            "mane_select": "true", "flag": "gene_variant", "warnings": json.dumps(list(warnings))}


def _excerpt(rows: list[dict]) -> list[dict]:
    out = [_entry(f'{r["RefSeq transcript (MANE Select)"]}:{r["HGVS c."]}', H.REFSEQ, r["HGVS g."],
                  ["The current status of LRG_297 is pending"] if r["Gene"] == "BARD1" else [])
           for r in rows]
    out += [_entry(f'{r["Ensembl transcript (MANE Select)"]}:{r["HGVS c."]}', H.ENSEMBL, r["HGVS g."])
            for r in rows if r["Gene"] in H.SECOND]
    return out


def _manifest(rows, excerpt) -> dict:
    return {"list": {"rows": len(rows), "sha256": "list-sha"}, "excerpt": {"rows": len(excerpt), "sha256": "excerpt-sha"}}


def _run(rows, excerpt, manifest=None, list_sha="list-sha", excerpt_sha="excerpt-sha", spoil=None):
    check, problems = H.build(rows, excerpt)
    summed = H.summary(check)
    if spoil:
        spoil(check)
    failed = H.checks(rows, excerpt, manifest or _manifest(rows, excerpt), list_sha, excerpt_sha, check, summed, problems)
    return check, summed, {k: v for k, v in failed.items() if v}


def test_a_list_and_its_excerpt_hold_together_and_every_conclusion_is_yes():
    rows = _list()
    check, summed, failed = _run(rows, _excerpt(rows))
    assert failed == {}
    assert [r["HGVS c."] for r in check] == [r["HGVS c."] for r in rows]
    assert all(r[c] == "Yes" for r in check for c in H.CONCLUSIONS)
    assert [list(r) for r in check] == [H.COLUMNS] * len(rows)
    # the second query's columns are filled for the genes submitted a second time, and for no other
    assert [bool(r[H.AGAIN[0]]) for r in check] == [r["Gene"] in H.SECOND for r in rows]
    assert all(r[H.AGAIN[4]] == r[H.AGAIN[5]] == ("Yes" if r["Gene"] in H.SECOND else "") for r in check)
    total = summed[-1]
    assert (total["Set"], total["Gene"], total["Variants"]) == (H.ALL, H.ALL, 5)
    assert [(s["Gene"], s["Variants"], s[H.SUMMARY[5]], s[H.SUMMARY[6]]) for s in summed[:-1]] == \
        [("BAP1", 2, 0, 0), ("BARD1", 1, 1, 0), ("SBDS", 2, 0, 2)]
    assert check[2]["Warnings"] == "The current status of LRG_297 is pending"


def test_a_returned_description_that_differs_is_printed_as_returned_with_no():
    """What the tool changed is recorded, not corrected: the row says No, the summary counts
    it, and the excerpt is then no longer the one the manifest records."""
    rows = _list()
    excerpt = _excerpt(rows)
    excerpt[0]["transcript_description"] = "NM_004656.4:c.37+3="
    excerpt[3]["genomic_description_grch38"] = "NC_000007.14:g.1A>C"
    check, summed, failed = _run(rows, excerpt, excerpt_sha="another-sha")
    assert list(failed) == ["the_list_and_the_excerpt_are_the_manifests"]
    assert check[0][H.FIRST[1]] == "NM_004656.4:c.37+3=" and check[0][H.CONCLUSIONS[0]] == "No"
    assert check[0][H.CONCLUSIONS[1]] == "Yes"
    assert check[3][H.CONCLUSIONS[1]] == "No" and check[3][H.CONCLUSIONS[0]] == "Yes"
    assert summed[-1][H.SUMMARY[3]] == 4 and summed[-1][H.SUMMARY[4]] == 4


def test_a_row_of_the_excerpt_taken_away_is_caught():
    rows = _list()
    for k in (1, 5):                                   # an entry on a RefSeq transcript, one on an Ensembl transcript
        excerpt = _excerpt(rows)
        gone = excerpt.pop(k)
        _, _, failed = _run(rows, excerpt, manifest=_manifest(rows, _excerpt(rows)))
        assert failed["every_row_has_its_one_entry"][0][1:3] == (gone["interface"], gone["submitted_description"])
        assert "the_list_and_the_excerpt_are_the_manifests" in failed          # the manifest's count of rows


def test_a_row_too_many_in_the_excerpt_is_caught():
    rows = _list()
    stray = _excerpt(rows) + [_entry("NM_000059.4:c.68-7T>A", H.REFSEQ, "NC_000013.11:g.32316415T>A")]
    _, _, failed = _run(rows, stray, manifest=_manifest(rows, stray))
    assert failed == {"the_excerpt_has_no_other_row": [(H.REFSEQ, "NM_000059.4:c.68-7T>A")]}
    # an Ensembl description of a gene that was not submitted a second time is one too many, too
    extra = _excerpt(rows) + [_entry("ENST00000460680.6:c.37+3G>A", H.ENSEMBL, rows[0]["HGVS g."])]
    _, _, failed = _run(rows, extra, manifest=_manifest(rows, extra))
    assert failed == {"the_excerpt_has_no_other_row": [(H.ENSEMBL, "ENST00000460680.6:c.37+3G>A")]}
    twice = _excerpt(rows) + [_excerpt(rows)[2]]
    _, _, failed = _run(rows, twice, manifest=_manifest(rows, twice))
    assert failed["every_row_has_its_one_entry"] == [(3, H.REFSEQ, "NM_000465.4:c.2002-3C>G", 2)]


def test_a_description_must_be_the_rows_character_for_character():
    rows = _list()
    excerpt = _excerpt(rows)
    excerpt[0]["submitted_description"] = "NM_004656.4:c.37+3G>a"
    _, _, failed = _run(rows, excerpt)
    assert failed["every_row_has_its_one_entry"] == [(1, H.REFSEQ, "NM_004656.4:c.37+3G>A")]
    assert failed["the_excerpt_has_no_other_row"] == [(H.REFSEQ, "NM_004656.4:c.37+3G>a")]


def test_another_list_than_the_one_the_replies_were_fetched_for_is_caught():
    rows = _list()
    _, _, failed = _run(rows, _excerpt(rows), list_sha="the-sha-of-another-list")
    assert failed == {"the_list_and_the_excerpt_are_the_manifests": [("list", "the-sha-of-another-list", "list-sha")]}


def test_a_summary_that_is_not_the_rows_counted_again_is_caught():
    """The summary is held to a second count, made from the descriptions themselves: a
    description altered after the summary was made, or a conclusion that does not follow
    from its descriptions, shows."""
    rows = _list()

    def alter(check):
        check[1][H.FIRST[1]] = "NM_004656.4:c.37+3G>T"

    _, _, failed = _run(rows, _excerpt(rows), spoil=alter)
    assert list(failed) == ["the_summary_is_the_rows_counted_again"]
    assert failed["the_summary_is_the_rows_counted_again"][0][0]["Gene"] == "BAP1"

    def drop(check):
        del check[4]

    _, _, failed = _run(rows, _excerpt(rows), spoil=drop)
    assert list(failed) == ["the_summary_is_the_rows_counted_again"]


def test_several_warnings_are_joined_and_a_warning_that_holds_the_separator_stops_the_run():
    rows = _list()
    excerpt = _excerpt(rows)
    excerpt[0]["warnings"] = json.dumps(["first", "second"])
    assert H.build(rows, excerpt)[0][0]["Warnings"] == "first | second"
    excerpt[0]["warnings"] = json.dumps(["one | two"])
    with pytest.raises(SystemExit, match="separates warnings"):
        H.build(rows, excerpt)


def test_the_files_text_is_comment_lines_then_the_header_then_the_rows():
    rows = _list()
    check, _ = H.build(rows, _excerpt(rows))
    text = H.as_text(H.COLUMNS, check, ["first line", "second, with a comma"])
    lines = text.split("\n")
    assert lines[:2] == ["# first line", "# second, with a comma"]
    assert lines[2].startswith("Set,Gene,RefSeq transcript (MANE Select),HGVS c.,HGVS g.,Description submitted,")
    assert len(lines) == 2 + 1 + len(rows) + 1 and lines[-1] == ""


def test_the_manifest_of_the_layout_gains_the_two_files_once(tmp_path):
    """The layout's manifest.csv keeps every line it had; the two files' rows go to its end,
    and a second run replaces them."""
    a, b, src = tmp_path / "hgvs_check.csv", tmp_path / "hgvs_check_summary.csv", tmp_path / "source.csv"
    for p, t in ((a, "one\n"), (b, "two\n"), (src, "source\n")):
        p.write_text(t)
    m = tmp_path / "manifest.csv"
    before = ('file,bytes,sha256,number,content,sources,sources_sha256\n'
              'phase1/x.csv,3,abc,Table 1,"a, b",phase1/s.csv; phase1/t.csv,d; e\n'
              'phase1/contents.csv,5,def,,,,\n')
    m.write_text(before)
    H.add_to_manifest([a, b], "Variant list check", [src], m)
    once = m.read_text()
    assert once.startswith(before) and len(once.splitlines()) == 5
    row = once.splitlines()[3].split(",")
    assert row[0] == H.REPO_PREFIX + a.as_posix() and row[1] == "4" and row[2] == hashlib.sha256(b"one\n").hexdigest()
    assert row[3:5] == ["Variant list check", H.CONTENT] and row[6] == hashlib.sha256(b"source\n").hexdigest()
    a.write_text("one more\n")
    H.add_to_manifest([a, b], "Variant list check", [src], m)
    again = m.read_text()
    assert again.startswith(before) and len(again.splitlines()) == 5
    assert again.splitlines()[3].split(",")[1] == "9"


def test_the_stage_asks_no_service():
    source = (PHASE1 / "src" / "evid_hgvs_check.py").read_text()
    assert not re.search(r"^\s*(?:import|from)\s+(?:urllib|requests|socket|http|ssl|ftplib|subprocess)\b", source, re.M)


# ----------------------------------------------------------------- the committed files
def _committed():
    mp = pytest.MonkeyPatch()
    mp.chdir(PHASE1)
    try:
        rows, excerpt = H.read_csv(H.LIST), H.read_csv(H.EXCERPT)
        manifest = json.loads(H.MANIFEST.read_text())
        check, problems = H.build(rows, excerpt)
        summed = H.summary(check)
        failed = H.checks(rows, excerpt, manifest, H.sha256(H.LIST), H.sha256(H.EXCERPT), check, summed, problems)
        return rows, excerpt, manifest, check, summed, failed
    finally:
        mp.undo()


@need_excerpt
def test_the_committed_inputs_hold_together():
    rows, excerpt, manifest, check, summed, failed = _committed()
    assert {k: v[:3] for k, v in failed.items() if v} == {}
    assert len(rows) == 15889 == summed[-1]["Variants"]
    second = sum(r["Gene"] in H.SECOND for r in rows)
    assert len(excerpt) == len(rows) + second == manifest["descriptions"]["submitted"]
    assert manifest["descriptions"]["by_interface"] == {H.REFSEQ: len(rows), H.ENSEMBL: second}
    assert list(excerpt[0]) == manifest["excerpt"]["columns"]
    assert {e["interface"] for e in excerpt} == {H.REFSEQ, H.ENSEMBL}


@need_excerpt
def test_the_manifest_says_where_the_replies_came_from():
    m = json.loads((PHASE1 / H.MANIFEST).read_text())
    assert m["tool"] == "VariantValidator" and m["interface"]["address"] == "https://rest.variantvalidator.org"
    for k in ("hello_at_start", "hello_at_end"):
        assert set(json.loads(m[k]["text"])["metadata"]) == set(H.VERSION_KEYS)
    run = m["run"]
    assert re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", run["first_request_utc"]) and run["first_request_utc"] <= run["last_request_utc"]
    assert run["requests_stored"] == m["descriptions"]["requests"]
    assert run["requests_made"] == run["requests_stored"] + run["requests_made_again"]
    assert run["list_sha256_the_run_read"] == [m["list"]["sha256"]]
    archive = m["raw_replies_archive"]
    assert re.fullmatch(r"[0-9a-f]{64}", archive["sha256"]) and archive["bytes"] > 0 and archive["files"] > m["descriptions"]["requests"]
    again = m["second_asking"]
    assert again["rows"] == 50 and again["entries"] >= 50
    # what was not submitted, and the tool's own list that shows why
    not_asked = m["not_submitted"]
    held = (PHASE1 / H.MANIFEST).parent / f"transcripts_{not_asked['gene']}_ensembl.json"
    assert not_asked["transcript"].split(".")[0] not in held.read_text()
    assert not_asked["gene"] not in H.SECOND
    for gene in H.SECOND:
        assert ((PHASE1 / H.MANIFEST).parent / f"transcripts_{gene}_ensembl.json").is_file()


@need_out
def test_the_two_files_are_what_the_stage_writes_from_the_committed_inputs():
    rows, _, manifest, check, summed, _ = _committed()
    comments = H.comment_lines(manifest)
    assert (OUT / H.FILES[0]).read_text() == H.as_text(H.COLUMNS, check, comments)
    printed = [{k: (H.grouped(v) if isinstance(v, int) else v) for k, v in r.items()} for r in summed]
    assert (OUT / H.FILES[1]).read_text() == H.as_text(H.SUMMARY, printed, comments)
    # the list's five columns, as they stand, in its order
    assert [[r[c] for c in H.FROM_LIST] for r in check] == [[r[c] for c in H.FROM_LIST] for r in rows]
    # the comment lines name the tool, the five versions and the days of the queries
    head = "\n".join(comments)
    hello = json.loads(manifest["hello_at_start"]["text"])["metadata"]
    assert "VariantValidator" in head and all(f"{k} {hello[k]}" in head for k in H.VERSION_KEYS)
    assert manifest["run"]["first_request_utc"][:10] in head and manifest["not_submitted"]["transcript"] in head


@need_out
def test_the_stage_is_registered_where_the_layouts_are():
    import yaml
    sys.path.insert(0, str(REPO / "scripts"))
    import reproduce_evidence as E
    mods = [m for m, _, _ in E.STAGES]
    assert mods[-2:] == ["src.evid_layout", "src.evid_hgvs_check"] and len(mods) == 42
    item = [it for it in yaml.safe_load((PHASE1 / H.CONFIG).read_text())["variants"] if it["content"] == H.CONTENT]
    assert len(item) == 1 and item[0]["make"] == "later_stage"
    assert item[0]["files"] == [str(H.OUT / f) for f in H.FILES]
    assert sorted(item[0]["sources"]) == sorted([str(H.EXCERPT), str(H.MANIFEST)])
    import pandas as pd
    c = pd.read_csv(OUT.parent / "contents.csv", dtype=str, keep_default_na=False).set_index("content")
    assert c.loc[H.CONTENT, "files"] == "; ".join(H.REPO_PREFIX + str(H.OUT / f) for f in H.FILES)
    assert c.loc[H.CONTENT, "number"] == item[0]["number"]
    m = pd.read_csv(OUT.parent / "manifest.csv", dtype=str, keep_default_na=False)
    mine = m[m.content == H.CONTENT]
    assert list(mine.file) == [H.REPO_PREFIX + str(H.OUT / f) for f in H.FILES]
    assert list(m.file[-2:]) == list(mine.file)                       # added at its end, after the layout's own
    # the excerpt and its record are inputs the stages only read
    inputs = (REPO / "scripts" / "reproduce_inputs.txt").read_text().split("\n")
    for name in ("excerpt.csv", "manifest.json", "transcripts_DDX3X_ensembl.json", "transcripts_SBDS_ensembl.json"):
        assert f"phase1/data/evidence/hgvs_check/{name}" in inputs
    assert "(42 stages" in (REPO / "README.md").read_text()
