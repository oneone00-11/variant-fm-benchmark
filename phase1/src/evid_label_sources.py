"""The label audit's sources as served and the compact extracts committed in their place
(E16c, decided by the author on 2026-10-21 in E16d).

Three of MaveDB's per-variant memberships of class-based calibrations (BAP1's one,
RAD51C's two; 9-18 MB each) and the BRCA2 paper's Supplementary Table 3 (an xlsx) are
not committed. They are kept, as served, in a gitignored cache that
scripts/reproduce_evidence.py --fetch-inputs fills and checks by sha256; what is
committed is a compact extract of each (variant -> class), from which the audit runs
offline. VHL's membership (1.4 MB) is committed as served, with its extract beside it.
These functions read both forms the same way, so that the fetch script, the audit,
--fetch-inputs and the tests compare like with like. Standard library only.
"""
from __future__ import annotations

import csv
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

MEMBERSHIP_COLUMNS = ["variant_urn", "hgvs_nt", "functional_classification_id", "label"]
SUPPLEMENTARY_COLUMNS = ["hgvs_nt", "category", "score"]
_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
       "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def membership_rows(blobs: list, calibration: dict) -> list[dict]:
    """A calibration's per-variant membership, in the order served: one row per variant,
    with the functional classification it is served under and that class's label."""
    label = {c["id"]: c["label"] for c in calibration["functionalClassifications"]}
    rows = []
    for blob in blobs:
        for group in blob if isinstance(blob, list) else [blob]:
            fc = group["functionalClassificationId"]
            for v in group["variants"]:
                rows.append({"variant_urn": v["urn"], "hgvs_nt": v.get("hgvsNt") or "",
                             "functional_classification_id": str(fc), "label": label[fc]})
    return rows


def _col(ref: str) -> str:
    return re.match(r"[A-Z]+", ref).group()


def xlsx_sheet(path: Path, sheet: str) -> dict[int, dict[str, str]]:
    """Every cell of one sheet as text, {row number: {column letter: value}}, read with
    the standard library (shared and inline strings resolved; numbers as stored)."""
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rid = next(s.get(f"{{{_NS['r']}}}id") for s in wb.find("m:sheets", _NS)
                   if s.get("name") == sheet)
        target = next(r.get("Target") for r in rels if r.get("Id") == rid)
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", _NS):
                shared.append("".join(t.text or "" for t in si.iter(f"{{{_NS['m']}}}t")))
        root = ET.fromstring(z.read("xl/" + target.lstrip("/").removeprefix("xl/")))
    rows = {}
    for row in root.iter(f"{{{_NS['m']}}}row"):
        vals = {}
        for c in row.findall("m:c", _NS):
            v = c.find("m:v", _NS)
            if c.get("t") == "s" and v is not None:
                vals[_col(c.get("r"))] = shared[int(v.text)]
            elif c.get("t") == "inlineStr":
                vals[_col(c.get("r"))] = "".join(t.text or "" for t in c.iter(f"{{{_NS['m']}}}t"))
            elif v is not None:
                vals[_col(c.get("r"))] = v.text
        rows[int(row.get("r"))] = vals
    return rows


def supplementary_rows(path: Path, spec: dict) -> tuple[list[dict], dict]:
    """The rows of a supplementary table that carry a class: the c. HGVS, the class and
    the score, in sheet order; and the sheet's title cell and header row as written."""
    sheet = xlsx_sheet(path, spec["sheet"])
    header = sheet[spec["header_row"]]
    col = {name: next(k for k, v in header.items() if v.strip() == name)
           for name in (spec["hgvs_column"], spec["class_column"], spec["score_column"])}
    rows = []
    for i in range(spec["header_row"] + 1, max(sheet) + 1):
        r = sheet.get(i, {})
        if not (r.get(col[spec["class_column"]]) or "").strip():
            continue
        rows.append({"hgvs_nt": r[col[spec["hgvs_column"]]].strip(),
                     "category": r[col[spec["class_column"]]].strip(),
                     "score": (r.get(col[spec["score_column"]]) or "").strip()})
    return rows, {"title_cell": sheet[1].get("A"), "header": header}


def extract_differences(manifest: dict, root: Path) -> tuple[list[str], list[str]]:
    """Every compact extract a label_audit manifest lists, against the files it was made
    from (cached or committed), row for row. Returns (extracts that differ, extracts
    whose source files are not here); the caller decides whether an absent source is
    allowed."""
    import json
    differ, absent = [], []
    metas = {Path(e["path"]).name: root / e["path"] for e in manifest["files"]}
    for e in manifest["files"]:
        if "compact_of" not in e:
            continue
        srcs = [root / p for p in e["compact_of"]]
        if not all(p.is_file() for p in srcs):
            absent.append(e["path"])
            continue
        if "supplementary" in e:
            rows, _ = supplementary_rows(srcs[0], e["supplementary"])
        else:
            meta = json.loads(metas[f"scoreset_meta_{e['score_set'].split(':')[-1]}.json"]
                              .read_text())
            cal = next(c for c in meta["scoreCalibrations"] if c["urn"] == e["calibration_urn"])
            rows = membership_rows([json.loads(p.read_text()) for p in srcs], cal)
        if read_tsv(root / e["path"]) != rows:
            differ.append(f"{e['path']} differs from {', '.join(e['compact_of'])}")
    return differ, absent
