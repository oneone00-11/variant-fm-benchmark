"""E16d: build the seven training genes' label files from their recorded sources.

No script in the repository produced the label files in phase1/data/assay_labels/
(they entered it on 2026-07-08, commit 10d94bb). This one builds each from the source
the label source audit (stage E16c) fetched and pinned by sha256, by one rule per gene
that the author fixed on 2026-10-21 and that must not be changed:

  BRCA1   urn:mavedb:00000097-0-2   R1, the investigator-provided ranges
  BRCA2   Supplementary Table 3     R6, its functional category, grouped as the
          of the paper              paper's Results group it
  BARD1   urn:mavedb:00001250-a-2   R1, the investigator-provided ranges
  PALB2   urn:mavedb:00001259-a-2   R1, the investigator-provided ranges
  RAD51C  urn:mavedb:00000673-0-1   R2, the score table's functional_classification
  VHL     urn:mavedb:00000675-a-1   R1, MaveDB's membership of the class-based
                                    calibration: the 2,200 variants in a class
  BAP1    urn:mavedb:00000662-0-1   R4, the paper's rule (FDR < 0.01 and the sign;
                                    depleted split at the printed median -0.1260642;
                                    the one variant on it is strongly depleted, as the
                                    file in use has it: the author's decision)

The classes are written in the wording of the file in use, which config.ASSAY_LABELS
maps (WORDING, the archive's word to the file's); the columns, delimiter, form of the
HGVS and row order follow the file in use too (FORMAT). An unassigned variant, or a
score on a point two ranges share, stops the script.

Writes phase1/data/assay_labels/_generated/<label file>, never a label file in use,
and phase1/data/assay_labels/PROVENANCE.json: per gene the source (MaveDB version, or
paper supplementary table), its files and sha256, the calibration used with its id
and modification date, the rule, the wording, this script's sha256, the output's
sha256, and whether the label file in use is this output ("in_use"). The BRCA2 table
is read from the label audit's cache (scripts/reproduce_evidence.py --fetch-inputs
fills it); where it is not there the gene is skipped, with the reason printed, unless
--require-all is given.

Run from the repository root:
    .venv/bin/python scripts/stage_E16d/build_labels.py
    .venv/bin/python scripts/stage_E16d/build_labels.py --require-all
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of  # noqa: E402

PHASE1 = REPO / "phase1"
sys.path.insert(0, str(PHASE1))
from src import config as C                    # noqa: E402
from src import evid_label_audit as A          # noqa: E402
from src import evid_label_sources as LS       # noqa: E402
from src import evid_newgenes as N             # noqa: E402

AUDIT = PHASE1 / "data" / "evidence" / "label_audit"
LABELS = PHASE1 / "data" / "assay_labels"
OUT = LABELS / "_generated"
PROVENANCE = LABELS / "PROVENANCE.json"
INVESTIGATOR = "Investigator-provided functional classes"

# The rule of each gene (the author's, 2026-10-21) and the archive's word for each
# class mapped to the word of the file in use; the latter are config.ASSAY_LABELS keys.
RULES = {
    "BRCA1": {"rule": "R1", "version": "urn:mavedb:00000097-0-2", "calibration": INVESTIGATOR},
    "BRCA2": {"rule": "R6"},
    "BARD1": {"rule": "R1", "version": "urn:mavedb:00001250-a-2", "calibration": INVESTIGATOR},
    "PALB2": {"rule": "R1", "version": "urn:mavedb:00001259-a-2", "calibration": INVESTIGATOR},
    "RAD51C": {"rule": "R2", "version": "urn:mavedb:00000673-0-1",
               "column": "functional_classification"},
    "VHL": {"rule": "R1", "version": "urn:mavedb:00000675-a-1", "calibration": INVESTIGATOR},
    "BAP1": {"rule": "R4", "version": "urn:mavedb:00000662-0-1",
             "paper_rule": "Results and Methods: FDR and sign, depleted split at the printed median"},
}
WORDING = {
    "BRCA1": {"Functional": "functional", "Intermediate": "intermediate",
              "Non-functional": "non-functional"},
    "BRCA2": {"P strong": "abnormal", "P moderate": "abnormal", "P supporting": "abnormal",
              "VUS": "intermediate",
              "B supporting": "normal", "B moderate": "normal", "B strong": "normal"},
    "BARD1": {"Functionally Abnormal": "functionally abnormal", "Indeterminate": "indeterminate",
              "Functionally Normal": "normal"},
    "PALB2": {"Functionally Abnormal": "abnormal", "Indeterminate": "indeterminate",
              "Functionally Normal": "normal"},
    "RAD51C": {"fast depleted": "Fast depleted", "slow depleted": "Slow depleted",
               "unchanged": "Unchanged", "enriched": "Enriched"},
    "VHL": {"Neutral": "neutral", "Intermediate": "intermediate", "LOF1": "LOF1", "LOF2": "LOF2"},
    # the variant on the printed median: strongly depleted, the author's decision of
    # 2026-10-21 (the paper does not say on which side it falls)
    "BAP1": {"strongly depleted": "strongly depleted", "weakly depleted": "weakly depleted",
             "depleted (on the printed median; side not stated)": "strongly depleted",
             "unchanged": "unchanged", "enriched": "enriched"},
}
# How each file in use writes its rows: the HGVS with or without the transcript, and
# the order (the source's, or by c. position, intron offset and alternative base)
FORMAT = {
    "BRCA1": {"hgvs": "without transcript", "order": "c. position"},
    "BRCA2": {"hgvs": "as the source writes it", "order": "source"},
    "BARD1": {"hgvs": "as the source writes it", "order": "source"},
    "PALB2": {"hgvs": "as the source writes it", "order": "source"},
    "RAD51C": {"hgvs": "as the source writes it", "order": "source"},
    "VHL": {"hgvs": "without transcript", "order": "source"},
    "BAP1": {"hgvs": "as the source writes it", "order": "source"},
}
_C = re.compile(r"^c\.([-*]?)(\d+)([+-]\d+)?([ACGT])>([ACGT])$")


class Skip(Exception):
    """A source is not here; the gene is skipped (unless --require-all)."""


def _c_position(hgvs: str) -> tuple:
    pre, pos, off, _ref, alt = _C.match(hgvs.split(":")[-1]).groups()
    p = int(pos)
    p = -p if pre == "-" else (10 ** 6 + p if pre == "*" else p)
    return (p, int(off or 0), "ACGT".index(alt))


def _manifest(gene: str) -> dict:
    return json.loads((AUDIT / f"{gene.lower()}_raw" / "manifest.json").read_text())


def _path(man: dict, name: str) -> Path:
    return REPO / next(e["path"] for e in man["files"] if Path(e["path"]).name == name)


def _table(man: dict, urn: str) -> tuple[list[dict], Path]:
    path = _path(man, f"urn_mavedb_{urn.split(':')[-1]}.csv")
    with open(path, newline="") as f:
        return list(csv.DictReader(f)), path


def _calibration(man: dict, urn: str, title: str) -> tuple[dict, Path]:
    path = _path(man, f"scoreset_meta_{urn.split(':')[-1]}.json")
    meta = json.loads(path.read_text())
    hit = [c for c in meta["scoreCalibrations"] if c.get("investigatorProvided") is True
           and c.get("title") == title]
    if len(hit) != 1:
        raise SystemExit(f"{urn}: {len(hit)} investigator calibrations titled {title!r}")
    return hit[0], path


def source_classes(gene: str) -> tuple[list[tuple[str, str]], dict]:
    """[(hgvs as the source writes it, the archive's class)] in source order, and the
    record of the source for PROVENANCE.json."""
    spec, man = RULES[gene], _manifest(gene)
    if spec["rule"] == "R6":
        sup = next(e for e in man["files"] if "supplementary" in e)
        raw = next(e for e in man["files"] if e.get("cached") and e["path"] in sup["compact_of"])
        xlsx = REPO / raw["path"]
        if not xlsx.is_file():
            raise Skip(f"{raw['path']} is not in the cache (scripts/reproduce_evidence.py "
                       "--fetch-inputs)")
        if sha256_of(xlsx) != raw["sha256"]:
            raise SystemExit(f"{raw['path']} differs from the sha256 its manifest records")
        rows, layout = LS.supplementary_rows(xlsx, sup["supplementary"])
        return ([(r["hgvs_nt"], r["category"]) for r in rows],
                {"type": "paper supplementary table", "pmcid": raw["pmcid"],
                 "file": raw["zip_member"], "sheet": raw["sheet"], "title_cell": layout["title_cell"],
                 "files": [{"path": raw["path"], "sha256": raw["sha256"],
                            "source_url": raw["source_url"], "committed": False}],
                 "rule": ("R6: the table's 'Functional category' of each variant, grouped as the "
                          "paper's Results group it (label_audit/brca2_raw/paper_rule.json)")})
    table, tpath = _table(man, spec["version"])
    files = [{"path": rel(tpath), "sha256": sha256_of(tpath), "committed": True}]
    record = {"type": "MaveDB score set version", "urn": spec["version"]}
    if spec["rule"] == "R2":
        out = [(r["hgvs_nt"], r[spec["column"]]) for r in table]
        record |= {"rule": f"R2: the score table's column {spec['column']}"}
    elif spec["rule"] == "R4":
        prpath = AUDIT / f"{gene.lower()}_raw" / "paper_rule.json"
        pr = json.loads(prpath.read_text())
        rule = next(r for r in pr["rules"] if r["name"] == spec["paper_rule"])
        got = A.rule_paper({"table": table, "genomic": False}, rule)
        out = [(r["hgvs_nt"], got[A._norm_hgvs(r["hgvs_nt"])][0]) for r in table]
        files.append({"path": rel(prpath), "sha256": sha256_of(prpath), "committed": True})
        record |= {"rule": f"R4: the paper's rule '{rule['name']}' ({pr['source']['citation_ids']})"}
    else:
        cal, mpath = _calibration(man, spec["version"], spec["calibration"])
        files.append({"path": rel(mpath), "sha256": sha256_of(mpath), "committed": True})
        record |= {"calibration": {"title": cal["title"], "urn": cal["urn"], "id": cal["id"],
                                   "created": cal.get("creationDate", "not_provided"),
                                   "modified": cal.get("modificationDate", "not_provided")}}
        if A._ranges(cal):
            classes = cal["functionalClassifications"]
            out = []
            for r in table:
                label, why = N.classify(r["score"], classes)
                if label == N.UNASSIGNED:
                    raise SystemExit(f"{gene}: {r['hgvs_nt']} unassigned under the ranges ({why}); "
                                     "stopped (E16d.6)")
                out.append((r["hgvs_nt"], label))
            record |= {"rule": "R1: the score placed in the investigator-provided ranges"}
        else:
            entry = next(e for e in man["files"] if e.get("calibration_urn") == cal["urn"]
                         and "compact_of" in e)
            label = {str(c["id"]): c["label"] for c in cal["functionalClassifications"]}
            member = {m["variant_urn"]: label[m["functional_classification_id"]]
                      for m in LS.read_tsv(REPO / entry["path"])}
            out = [(r["hgvs_nt"], member[r["accession"]]) for r in table if r["accession"] in member]
            files.append({"path": entry["path"], "sha256": entry["sha256"], "committed": True})
            record |= {"rule": ("R1: MaveDB's membership of the class-based calibration; the "
                                f"{len(table) - len(out)} score-table variants in no class are left "
                                "out (the author's decision of 2026-10-21)")}
    return out, record | {"files": files}


def build(gene: str) -> tuple[bytes, dict]:
    spec = C.ASSAY_LABELS[gene]
    rows, record = source_classes(gene)
    words, fmt = WORDING[gene], FORMAT[gene]
    missing = sorted({c for _, c in rows} - set(words))
    if missing:
        raise SystemExit(f"{gene}: classes outside the wording table: {missing}")
    lines = [(h.split(":")[-1] if fmt["hgvs"] == "without transcript" else h, words[c])
             for h, c in rows]
    if fmt["order"] == "c. position":
        lines.sort(key=lambda x: _c_position(x[0]))
    if not {c for _, c in lines} <= set(spec["map"]):
        raise SystemExit(f"{gene}: a class config.ASSAY_LABELS does not map")
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t" if spec["file"].endswith(".tsv") else ",",
                   lineterminator="\n")
    w.writerow([spec["key"], spec["class_col"]])
    w.writerows(lines)
    return buf.getvalue().encode(), record | {"wording": words, "format": fmt, "rows": len(lines)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--require-all", action="store_true",
                    help="fail, rather than skip a gene, when a source is not here")
    args = ap.parse_args()
    with RunLog(__file__) as log:
        OUT.mkdir(parents=True, exist_ok=True)
        old = json.loads(PROVENANCE.read_text())["genes"] if PROVENANCE.exists() else {}
        genes, skipped = {}, {}
        for gene, spec in C.ASSAY_LABELS.items():
            try:
                blob, record = build(gene)
            except Skip as e:
                if args.require_all:
                    raise SystemExit(f"{gene}: {e}")
                skipped[gene] = str(e)
                print(f"{gene}: skipped, {e}")
                if gene in old:
                    genes[gene] = old[gene]
                continue
            path = OUT / spec["file"]
            path.write_bytes(blob)
            log.output(path)
            in_use = LABELS / spec["file"]
            log.input(in_use)
            genes[gene] = {"label_file": rel(in_use),
                           "generated": {"path": rel(path), "sha256": sha256_of(path),
                                         "rows": record.pop("rows")},
                           "source": record, "in_use": in_use.read_bytes() == blob}
            print(f"{gene}: {genes[gene]['generated']['rows']:,} rows; the file in use is "
                  f"{'this file, byte for byte' if genes[gene]['in_use'] else 'not this file'}")
        doc = {"written_by": "scripts/stage_E16d/build_labels.py",
               "generator_sha256": sha256_of(Path(__file__)),
               "rules_fixed_by": "the author, 2026-10-21 (stage E16d)",
               "in_use": ("true when the label file config.ASSAY_LABELS names is, byte for "
                          "byte, the generated file"),
               "genes": genes}
        PROVENANCE.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")
        log.output(PROVENANCE)
        log.note("skipped", skipped)


if __name__ == "__main__":
    main()
