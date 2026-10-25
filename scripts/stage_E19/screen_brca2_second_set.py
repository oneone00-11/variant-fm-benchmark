"""E19 C3, BRCA2 (E19.5 item 3, as the author specified it on 2026-10-25): which of the
BRCA2 score sets is a second dataset, by the four conditions, from the records
scripts/stage_E19/fetch_brca2_second_set.py saved. Offline; run in the reference image
(修改单 §0 rule 9).

The conditions (the author's items 2 and 3), each read from a saved record:
  (a) an SGE score set: the catalog's own rule (scripts/20_catalog_mavedb.py,
      assay_type): 'saturation genome editing', or SGE as a word, in the score set's
      title, its experiment's abstract or method, or its own method (case folded). For
      every score set results/tables/mavedb_catalog.tsv lists, its is_sge must agree, or
      this stops.
  (b) the experiment comes from a laboratory other than Huang 2025's, judged on the
      authors of the papers, not on MaveDB's submitters or contributors: the papers
      MaveDB records on the score set (else on its experiment) do not include Huang
      2025, and no such paper shares its last author with Huang 2025 (the last author
      taken as the head of the laboratory) or names Huang 2025's last author among its
      authors, nor does Huang 2025 name its last author among its own. Authors are
      Crossref's, as served, compared on family name and first initial, case and
      accents folded. A score set with no paper recorded cannot be judged so: (b) is not
      met ('no paper').
  (c) it shares variants with the analysis set's BRCA2 intronic variants, which are
      551, all at 3-10 nt (BRCA2 has none at 11-50 nt): rows of the score table written
      on BRCA2's MANE Select pair (HGNC's record: ENST00000380152.8, NM_000059.4) whose
      c. notation, normalised as evid_build_set.attach_labels normalises it, is one of
      those 551. Rows on any other reference are counted, and not matched.
  (d) it has a per-variant experimental classification, taken by the label source
      audit's order (the author's item 4): R1, an investigator-provided calibration
      (evid_label_audit's R1); R2, a class column of the score table (evid_label_audit's
      CLASS_COLUMNS); R4, the rule the paper states (paper_rule.json in this folder,
      quoted and checked); R6, a per-variant table published with the paper, the
      experimental classes only. R1 and R2 are read for every score set; R4 and R6 only
      for one that meets (a) to (c), and until paper_rule.json is written they are
      pending.
The author's expectation (item 3) is checked against, not assumed: the condition each
score set fails is written beside it.

Written into phase1/data/evidence/second_set/brca2_raw/:
  screening.csv    one row per BRCA2 score set: the catalog's place for it, each
                   condition's facts and verdict, and the conditions it fails
  coverage_3_10.csv  the author's coverage record (review of 2026-10-25, item 3; recounted
                   from Supplementary Table 5 by the author's decision of 2026-10-26): the
                   analysis set's 551 BRCA2 variants at 3-10 nt by side (donor +3 ... +10,
                   acceptor -3 ... -10) and offset, how many of them Table 5 classes,
                   joined on the genomic coordinate after its row-by-row check
                   (evid_second_set.table5_rows), and how many the chosen score set's
                   MaveDB table carries, joined on c. as condition (c) counts them, with
                   the positions each covers; the MaveDB total must equal (c)'s count
  decision.json    the conditions as applied, the score sets meeting all four, the
                   outcome; and the author's item 7, recorded and acted on in no way:
                   whether urn:mavedb:00001279-a-1 has an investigator-provided
                   calibration, and whether its scores are 00001225-a-1's variant by
                   variant
Exactly one score set meeting all four is the second dataset; none means no usable
second set; more than one stops this script (the author's item 8).

    docker run ... IMAGE bash scripts/stage_E19/in_image.sh scripts/stage_E19/screen_brca2_second_set.py
"""
from __future__ import annotations

import csv
import io
import json
import re
import sys
import unicodedata
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

sys.path.insert(0, str(REPO / "phase1"))
from src.evid_label_audit import CLASS_COLUMNS  # noqa: E402
from src.phase1_method_b_labels import _norm_hgvs  # noqa: E402

RAW = REPO / "phase1/data/evidence/second_set/brca2_raw"
MANIFEST = RAW / "manifest.json"
PAPER_RULE = RAW / "paper_rule.json"
CATALOG = REPO / "results/tables/mavedb_catalog.tsv"
SET = REPO / "phase1/data/evidence/analysis_set_v1.parquet"
HUANG = {"doi": "10.1038/s41586-024-08388-8", "urn": "urn:mavedb:00001225-a-1"}
AFTER_CATALOG = "urn:mavedb:00001279-a-1"     # published after the catalog (item 1)
SCREENING_COLUMNS = [
    "urn", "title", "catalog", "license", "published_date", "experiment",
    "a_sge", "a_evidence", "a_catalog_is_sge",
    "b_papers", "b_last_authors", "b_same_paper_as_huang", "b_last_author_shared",
    "b_authors_shared_with_huang", "b_other_laboratory",
    "c_rows", "c_rows_on_mane_pair", "c_intron_side_snvs", "c_shared_3_10",
    "c_shared_1_2", "c_shared",
    "d_r1", "d_r2", "d_r4", "d_r6", "d_rule_used", "d_classification",
    "meets_all", "fails", "reason"]


def _fold(s: str) -> str:
    s = unicodedata.normalize("NFKD", s or "")
    return "".join(ch for ch in s if not unicodedata.combining(ch)).lower().strip()


def _person(a: dict) -> tuple[str, str] | None:
    """(family, first initial), folded, for a Crossref author who is a person."""
    if not a.get("family"):
        return None
    return _fold(a["family"]), _fold(a.get("given") or "")[:1]


def _crossref(doi: str) -> dict:
    path = RAW / "papers" / f"crossref_{doi.replace('/', '_')}.json"
    return json.loads(path.read_text())["message"]


def _papers(meta: dict, exp: dict) -> list[str]:
    def dois(rec):
        return [p["doi"].lower() for p in rec.get("primaryPublicationIdentifiers") or []
                if p.get("doi")]
    return dois(meta) or dois(exp)


def condition_a(meta: dict) -> tuple[bool, str]:
    exp = meta.get("experiment") or {}
    parts = {"title": meta.get("title") or "", "experiment abstract": exp.get("abstractText") or "",
             "experiment method": exp.get("methodText") or "", "method": meta.get("methodText") or ""}
    hits = [k for k, v in parts.items()
            if "saturation genome editing" in v.lower() or re.search(r"\bsge\b", v.lower())]
    return bool(hits), "; ".join(hits) or "none"


def condition_b(papers: list[str], huang: dict) -> dict:
    h_people = [p for p in (_person(a) for a in huang.get("author") or []) if p]
    h_last = h_people[-1]
    out = {"papers": papers, "last_authors": [], "same_paper": HUANG["doi"] in papers,
           "last_shared": False, "shared": 0}
    for doi in papers:
        people = [p for p in (_person(a) for a in _crossref(doi).get("author") or []) if p]
        out["last_authors"].append(" ".join(reversed(people[-1])) if people else "none")
        if people and (people[-1] == h_last or people[-1] in h_people or h_last in people):
            out["last_shared"] = True
        out["shared"] = max(out["shared"], len(set(people) & set(h_people)))
    out["met"] = bool(papers) and not out["same_paper"] and not out["last_shared"]
    return out


def score_rows(urn: str) -> list[dict]:
    path = RAW / "candidates" / f"{urn.replace(':', '_')}_scores.csv"
    return list(csv.DictReader(io.StringIO(path.read_text(encoding="utf-8"), newline="")))


def condition_c(rows: list[dict], mane: list[str], ours: dict[str, str]) -> dict:
    """Intron-side SNVs are counted on whatever reference a row is written on (one
    written without a reference is on the score set's own target sequence); only rows
    on the MANE Select pair are matched to the analysis set's variants."""
    on_pair, intronic, shared = 0, 0, {"s3_10": 0, "pm12": 0}
    for r in rows:
        nt = r.get("hgvs_nt") or ""
        tx, _, c = nt.partition(":") if ":" in nt else ("", "", nt)
        if re.fullmatch(r"c\.[-*]?\d+[+-]\d+[ACGT]>[ACGT]", c):
            intronic += 1
        if tx not in mane or not c.startswith("c."):
            continue
        on_pair += 1
        band = ours.get(_norm_hgvs(c))
        if band in shared:
            shared[band] += 1
    return {"rows": len(rows), "on_pair": on_pair, "intronic": intronic,
            "shared_3_10": shared["s3_10"], "shared_1_2": shared["pm12"],
            "met": shared["s3_10"] > 0}


SIGNED = re.compile(r"^c\.([-*]?\d+)([+-])(\d+)[ACGT]>[ACGT]$")


def coverage(table: list[dict], mane: list[str], b, table5: set) -> list[dict]:
    """By side and offset, the set's BRCA2 variants at 3-10 nt, Table 5's cover of them
    (genomic keys) and the MaveDB score table's (c. notation)."""
    ours = {}
    for c, st, key in zip(b["hgvs_c"], b["stratum"],
                          zip(b["chrom"].astype(str), b["pos"].astype(int), b["ref"], b["alt"])):
        m = SIGNED.match(c)
        if st == "s3_10" and m:
            ours.setdefault((m.group(2), int(m.group(3))), []).append((_norm_hgvs(c), m.group(1), key))
    mavedb = set()
    for r in table:
        tx, _, c = (r.get("hgvs_nt") or "").partition(":")
        if tx in mane:
            mavedb.add(_norm_hgvs(c))
    rows = []
    for sign, side in (("+", "donor"), ("-", "acceptor")):
        for off in range(3, 11):
            mine = ours.get((sign, off), [])
            in5 = [x for x in mine if x[2] in table5]
            inm = [x for x in mine if x[0] in mavedb]
            rows.append({"side": side, "offset": f"{sign}{off}", "n_ours": len(mine),
                         "n_ours_in_table5": len(in5), "n_ours_in_mavedb_score_set": len(inm),
                         "positions_ours": len({x[1] for x in mine}),
                         "positions_ours_in_table5": len({x[1] for x in in5})})
    return rows


def condition_d_r1_r2(meta: dict, rows: list[dict]) -> tuple[str, str]:
    cals = [c for c in meta.get("scoreCalibrations") or [] if c.get("investigatorProvided") is True
            and c.get("functionalClassifications")]
    r1 = ("usable: " + "; ".join(c.get("title") or "untitled" for c in cals)) if cals else \
        "none: no investigator-provided calibration"
    cols = [c for c in CLASS_COLUMNS if rows and c in rows[0]]
    r2 = ("usable: " + ", ".join(cols)) if cols else "none: no class column in the score table"
    return r1, r2


def condition_d_r4_r6() -> tuple[str, str]:
    if not PAPER_RULE.is_file():
        return ("pending: the paper's rule is not yet read (paper_rule.json)",
                "pending: R4 is not yet settled")
    pr = json.loads(PAPER_RULE.read_text())
    r4 = {"stated": "usable: the paper's rule, quoted (paper_rule.json)",
          "not_stated": f"not stated: {pr.get('not_applicable_because', '')}"}.get(
        pr["status"], f"unknown status {pr['status']!r}")
    sup = pr.get("supplementary_rule")
    r6 = ("not examined: R4 is usable" if pr["status"] == "stated" else
          f"usable: {sup['name']}" if sup else "none: no supplementary table named")
    return r4, r6


def main() -> None:
    import pandas as pd
    with RunLog(__file__) as log:
        log.input(MANIFEST)
        man = json.loads(MANIFEST.read_text())
        mane = man["brca2_identifiers"]["mane_select"]
        log.input(SET)
        s = pd.read_parquet(SET, columns=["gene", "hgvs_c", "stratum", "chrom", "pos", "ref", "alt"])
        b = s[s["gene"] == "BRCA2"]
        ours = dict(zip(b["hgvs_c"].map(_norm_hgvs), b["stratum"]))
        n_ours = int((b["stratum"] == "s3_10").sum())
        if set(b["stratum"]) - {"pm12", "s3_10"} or n_ours != 551:
            raise SystemExit("the analysis set's BRCA2 rows are not 551 at 3-10 nt and none "
                             "at 11-50 nt, as the author's item 2 reads them: stopped")
        log.input(CATALOG)
        with open(CATALOG, newline="") as f:
            catalog = {r["urn"]: r for r in csv.DictReader(f, delimiter="\t")
                       if r["gene"] == "BRCA2"}
        huang = _crossref(HUANG["doi"])
        r4, r6 = condition_d_r4_r6()
        rows, notes = [], {}
        for urn in man["candidates"]:
            safe = urn.replace(":", "_")
            meta = json.loads((RAW / "candidates" / f"{safe}.json").read_text())
            exp_urn = meta["experiment"]["urn"]
            exp = json.loads((RAW / "candidates" / f"{exp_urn.replace(':', '_')}.json").read_text())
            table = score_rows(urn)
            a, a_ev = condition_a(meta)
            cat = catalog.get(urn)
            if cat is not None and (cat["is_sge"] == "True") != a:
                raise SystemExit(f"{urn}: condition (a) gives {a}, the catalog's is_sge "
                                 f"{cat['is_sge']}: stopped")
            bb = condition_b(_papers(meta, exp), huang)
            cc = condition_c(table, mane, ours)
            d1, d2 = condition_d_r1_r2(meta, table)
            reach = a and bb["met"] and cc["met"]
            d4, d6 = (r4, r6) if reach else ("not examined", "not examined")
            used = next((name for name, v in (("R1", d1), ("R2", d2), ("R4", d4), ("R6", d6))
                         if v.startswith("usable")), "")
            pending = reach and not used and (d4.startswith("pending") or d6.startswith("pending"))
            if reach:
                d_state = f"met ({used})" if used else "pending" if pending else "not met"
            else:
                d_state = (f"met ({used})" if used in ("R1", "R2") else
                           "R1 and R2 none; R4 and R6 not examined")
            fails = [k for k, ok in (("(a)", a), ("(b)", bb["met"]), ("(c)", cc["met"]),
                                     ("(d)", d_state != "not met")) if not ok]
            why = []
            if not a:
                why.append("(a) not an SGE score set by the catalog's rule")
            if not bb["met"]:
                why.append("(b) " + ("no paper recorded" if not bb["papers"] else
                                     "Huang 2025's own paper" if bb["same_paper"] else
                                     "shares its laboratory head with Huang 2025"))
            if not cc["met"]:
                why.append(f"(c) no row is one of the 551 BRCA2 variants at 3-10 nt "
                           f"({cc['on_pair']:,} of {cc['rows']:,} rows on the MANE Select pair; "
                           f"{cc['intronic']:,} intron-side SNVs in all)")
            if reach and not used:
                why.append("(d) " + ("pending: the paper's rule" if pending else
                                     "no per-variant experimental classification"))
            rows.append({
                "urn": urn, "title": meta.get("title"),
                "catalog": ("in the catalog" if cat is not None else
                            "after the catalog, listed separately" if urn == AFTER_CATALOG else
                            "not in the catalog"),
                "license": (meta.get("license") or {}).get("shortName"),
                "published_date": meta.get("publishedDate"), "experiment": exp_urn,
                "a_sge": a, "a_evidence": a_ev,
                "a_catalog_is_sge": cat["is_sge"] if cat is not None else "",
                "b_papers": ";".join(bb["papers"]) or "none",
                "b_last_authors": ";".join(bb["last_authors"]) or "none",
                "b_same_paper_as_huang": bb["same_paper"],
                "b_last_author_shared": bb["last_shared"],
                "b_authors_shared_with_huang": bb["shared"], "b_other_laboratory": bb["met"],
                "c_rows": cc["rows"], "c_rows_on_mane_pair": cc["on_pair"],
                "c_intron_side_snvs": cc["intronic"], "c_shared_3_10": cc["shared_3_10"],
                "c_shared_1_2": cc["shared_1_2"], "c_shared": cc["met"],
                "d_r1": d1, "d_r2": d2, "d_r4": d4, "d_r6": d6, "d_rule_used": used,
                "d_classification": d_state,
                "meets_all": reach and bool(used), "fails": ", ".join(fails) or "none",
                "reason": "; ".join(why) or ("meets all four" if reach and used else "")})
            if urn == AFTER_CATALOG:
                ref = {_norm_hgvs(r["hgvs_nt"].partition(":")[2]): r["score"]
                       for r in score_rows(HUANG["urn"])}
                both = [(ref[k], r["score"]) for r in table
                        if (k := _norm_hgvs(r["hgvs_nt"].partition(":")[2])) in ref]
                same = sum(1 for x, y in both if float(x) == float(y))
                notes["item_7"] = {
                    "urn": urn,
                    "investigator_provided_calibration": d1.startswith("usable"),
                    "calibrations": [(c.get("title"), c.get("investigatorProvided"))
                                     for c in meta.get("scoreCalibrations") or []],
                    "variants_in_both": len(both), "scores_equal": same,
                    "same_scores_variant_by_variant": len(both) > 0 and same == len(both),
                    "compared_with": HUANG["urn"],
                    "line": (f"{urn}: investigator-provided calibration "
                             f"{'yes' if d1.startswith('usable') else 'no'}; its scores equal "
                             f"{HUANG['urn']}'s on {same:,} of the {len(both):,} variants both "
                             "carry (same c. notation); nothing is changed on this account")}
        meets = [r["urn"] for r in rows if r["meets_all"]]
        pend = [r["urn"] for r in rows if "pending" in r["reason"]]
        chosen = meets or pend
        if len(chosen) == 1 and r4.startswith("not stated") and r6.startswith("usable"):
            from src.evid_second_set import table5_rows
            t5, check = table5_rows(RAW)
            log.note("table5_check", check)
            cov = coverage(score_rows(chosen[0]), mane, b, {r["key"] for r in t5})
            shared = next(r["c_shared_3_10"] for r in rows if r["urn"] == chosen[0])
            if sum(r["n_ours"] for r in cov) != n_ours or \
                    sum(r["n_ours_in_mavedb_score_set"] for r in cov) != shared:
                raise SystemExit("the coverage record does not add up to (c)'s counts: stopped")
            with open(RAW / "coverage_3_10.csv", "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(cov[0]))
                w.writeheader()
                w.writerows(cov)
            log.output(RAW / "coverage_3_10.csv")
        if len(meets) > 1:
            raise SystemExit(f"more than one BRCA2 score set meets the four conditions: {meets}: "
                             "stopped (the author's item 8)")
        outcome = (f"use {meets[0]}" if meets else
                   f"pending: {', '.join(pend)} meets (a) to (c); (d) waits for the paper's rule"
                   if pend else "no usable second set: skipped")
        with open(RAW / "screening.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=SCREENING_COLUMNS)
            w.writeheader()
            w.writerows(rows)
        log.output(RAW / "screening.csv")
        decision = {"what": "the second BRCA2 dataset, by the author's four conditions "
                            "(2026-10-25)",
                    "written_by": rel(Path(__file__)), "conditions": __doc__.split(
                        "Written into")[0].strip(),
                    "analysis_set_brca2": {"s3_10": n_ours, "pm12": int((b["stratum"] == "pm12").sum()),
                                           "s11_50": 0},
                    "catalog_brca2": sorted(catalog), "candidates": [r["urn"] for r in rows],
                    "meeting_all": meets, "outcome": outcome, **notes}
        (RAW / "decision.json").write_text(json.dumps(decision, indent=1, ensure_ascii=False) + "\n")
        log.output(RAW / "decision.json")
        log.note("outcome", outcome)
        for r in rows:
            print(f"{r['urn']}: fails {r['fails']}; {r['reason']}")
        print(outcome)
        if "item_7" in notes:
            print(notes["item_7"]["line"])


if __name__ == "__main__":
    main()
