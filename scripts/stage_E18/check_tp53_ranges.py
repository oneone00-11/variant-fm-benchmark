"""E18.4: TP53's published classification ranges against E7's label, variant by variant.

From committed files and the tracked reference cache only, no network:

  (i)  the score set urn:mavedb:00001277-a-1 (tp53_raw/, fetched by
       scripts/stage_E18/fetch_tp53_ranges.py) covers every TP53 row E7 reads, all 288:
       each row's hgvs_nt names one row of the score table, and that row's c.
       description maps, by E16's own route (evid_newgenes: the atlas mapping module
       through the Mutalyzer transcript model of NM_000546.6, anchored to GRCh38 by the
       Ensembl lookup of the MANE Select partner, the stated reference base checked
       against the cached Ensembl slice), to the row's existing chrom, pos, ref and alt;
  (ii) the score set's investigator-provided ranges, read from its metadata and applied
       by E16's own rule (evid_newgenes.label_source and classify: exactly one
       investigator-provided calibration, ranges that apply mechanically), give every row
       the same class as y_control_anchored.

A row that is missing, maps elsewhere or is classed otherwise stops it (E18.4). The
cross-check E16 also makes against MaveDB's mapped-variant records is not made: their
file is 40.4 MB, above the 25 MB limit, and by the author's decision of 2026-10-23 the
288/288 identical coordinates of (i) stand in for it. Writes tp53_raw/range_check.json
(the result) and tp53_raw/range_check.tsv (row by row).

Run from the repository root:
    .venv/bin/python scripts/stage_E18/check_tp53_ranges.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of  # noqa: E402

PHASE1 = REPO / "phase1"
RAW = PHASE1 / "data/evidence/external/tp53_raw"
URN = "urn:mavedb:00001277-a-1"
TRANSCRIPT = "NM_000546.6"


def main() -> None:
    os.chdir(PHASE1)
    sys.path.insert(0, str(PHASE1))
    import pandas as pd
    from src import evid_external as X
    from src import evid_newgenes as NG
    name = X.TP53_PRIMARY_LABEL_NAME
    with RunLog(__file__) as log:
        man = json.loads((RAW / "manifest.json").read_text())
        for f in man["files"]:
            p = REPO / f["path"]
            if sha256_of(p) != f["sha256"]:
                raise SystemExit(f"{f['path']} is not the file its manifest records")
            log.input(p)
        meta = json.loads((RAW / "scoreset_meta.json").read_text())
        table = list(csv.DictReader((RAW / f"urn_mavedb_{URN.split(':')[-1]}.csv").open(newline="")))
        cals = meta.get("scoreCalibrations") or []
        blob = {"investigator_provided": [c for c in cals if c.get("investigatorProvided") is True]}
        src = NG.label_source("TP53", table, blob)
        if src["kind"] != "ranges":
            raise SystemExit("TP53: the score set's classes are not investigator-provided ranges")
        classes, fc = src["classes"], src["fc"]
        by_hgvs = {r["hgvs_nt"]: r for r in table}
        if len(by_hgvs) != len(table):
            raise SystemExit("TP53: an hgvs_nt names two rows of the score table")

        M = NG._atlas_mapping()
        mane = M.load_mane_records(["TP53"])["TP53"]
        tmap, prov = M.build_transcript_map(TRANSCRIPT, mane)
        lo = min(e.g_start for e in tmap.exons) - M.REGION_PAD
        hi = max(e.g_end for e in tmap.exons) + M.REGION_PAD
        s0, _, seq = M.fetch_ensembl_region(mane.chr_name, lo, hi)
        rseq = M.RegionSeq(s0, seq)

        df, _ = X._load_external("TP53")
        log.input(PHASE1 / X.TP53_TABLE)
        rows = []
        for r in df.itertuples(index=False):
            s = by_hgvs.get(r.hgvs_nt)
            m = NG.HGVS_NT.search(r.hgvs_nt or "")
            row = {"variant_id": r.variant_id, "hgvs_nt": r.hgvs_nt,
                   "score_set_row": s["accession"] if s else "", "mapped": "", "why": ""}
            if s is None or m is None:
                row["why"] = "no row of the score table" if s is None else "not an intron-side SNV"
            else:
                try:
                    pos, ref, alt = M.map_variant(tmap, rseq, m.group(2))
                    row["mapped"] = f"{mane.chr_name}-{pos}-{ref}-{alt}"
                except M.MappingFailure as e:
                    row["why"] = f"{e.category}: {e}"
            cls, why = NG.classify(s["score"], classes) if s else (NG.UNASSIGNED, "")
            label = NG.LABEL_OF.get(fc.get(cls))
            row |= {"score": s["score"] if s else "", "range_class": cls,
                    "range_label": "" if label is None else int(label),
                    "y_control_anchored": int(r.y_control_anchored),
                    "same_coordinates": row["mapped"] == r.variant_id,
                    "same_label": label is not None and label == r.y_control_anchored}
            rows.append(row)
        out = pd.DataFrame(rows)
        out.to_csv(RAW / "range_check.tsv", sep="\t", index=False)
        log.output(RAW / "range_check.tsv")
        cal = src["calibration"]
        result = {
            "urn": URN, "label_name": name,
            "rows_e7_reads": int(len(out)),
            "covered": int((out["score_set_row"] != "").sum()),
            "same_coordinates": int(out["same_coordinates"].sum()),
            "same_label": int(out["same_label"].sum()),
            "classes": out["range_class"].value_counts().to_dict(),
            "calibration": {"title": cal.get("title"), "investigatorProvided": True,
                            "ranges": [{k: c.get(k) for k in ("label", "functionalClassification",
                                                              "range", "inclusiveLowerBound",
                                                              "inclusiveUpperBound")}
                                       for c in classes]},
            "mapping": ("evid_newgenes' route (E16): the atlas mapping module, Mutalyzer "
                        f"model of {TRANSCRIPT}, the MANE Select partner's Ensembl lookup, the "
                        "reference base checked against the cached Ensembl slice; offline"),
            "not_done": ("the cross-check against MaveDB's mapped-variant records: their file is "
                         "40.4 MB, above the 25 MB limit; replaced by 288/288 coordinates "
                         "identical (the author's decision of 2026-10-23)"),
            "file": rel(RAW / "range_check.tsv"),
        }
        (RAW / "range_check.json").write_text(json.dumps(result, indent=1, ensure_ascii=False) + "\n")
        log.output(RAW / "range_check.json")
        log.note("result", {k: result[k] for k in ("rows_e7_reads", "covered", "same_coordinates",
                                                   "same_label")})
        print(f"TP53: {result['covered']}/{len(out)} covered, {result['same_coordinates']} at the "
              f"same coordinates, {result['same_label']} with the same label")
        if not (result["covered"] == result["same_coordinates"] == result["same_label"] == len(out) == 288):
            raise SystemExit("TP53's investigator-provided ranges do not reproduce E7's label "
                             "row for row: "
                             "stopped (E18.4)")


if __name__ == "__main__":
    main()
