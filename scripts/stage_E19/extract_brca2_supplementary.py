"""E19 C3, BRCA2: the compact extract of the second set's Supplementary Table 5, committed in
place of the xlsx (which stays in the gitignored cache, sha256-pinned), as the label
audit committed Huang 2025's Supplementary Table 3 (evid_label_sources).

The cached file is checked against the sha256 label_source_manifest.json records, then
every row of its sheet that carries a classification is written, by the label audit's own
reader (evid_label_sources.supplementary_rows): the c. HGVS, the class, the function score
and the g. HGVS, each as written (the g. HGVS because three rows write their c. wrongly and
the join is on the genomic coordinate, the author's decision of 2026-10-26). The
extract's entry is added to the manifest in the label audit's form, so
scripts/reproduce_evidence.py --fetch-inputs checks the extract against the file row for
row (evid_label_sources.extract_differences). No class is mapped here; the mapping is
paper_rule.json's.

Run in the reference image (the extract is committed; 修改单 §0 rule 9), offline:
    docker run ... IMAGE bash scripts/stage_E19/in_image.sh scripts/stage_E19/extract_brca2_supplementary.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of  # noqa: E402

sys.path.insert(0, str(REPO / "phase1"))
from src import evid_label_sources as LS  # noqa: E402

RAW = REPO / "phase1/data/evidence/second_set/brca2_raw"
MANIFEST = RAW / "label_source_manifest.json"
EXTRACT = RAW / "supplementary_table5_classes.tsv"
SPEC = {"sheet": "Sheet 1", "header_row": 1, "hgvs_column": "c.nom",
        "class_column": "Classification", "score_column": "Function.score",
        "genomic_column": "g.nom"}
COLUMNS = LS.SUPPLEMENTARY_COLUMNS + ["hgvs_g"]


def main() -> None:
    with RunLog(__file__) as log:
        log.input(MANIFEST)
        man = json.loads(MANIFEST.read_text())
        src = next(e for e in man["files"] if e.get("cached"))
        path = REPO / src["path"]
        if sha256_of(path) != src["sha256"]:
            raise SystemExit(f"{src['path']} is not the file the manifest records: stopped")
        log.input(path)
        rows, layout = LS.supplementary_rows(path, SPEC)
        LS.write_tsv(EXTRACT, COLUMNS, rows)
        log.output(EXTRACT)
        entry = {"path": rel(EXTRACT), "bytes": EXTRACT.stat().st_size,
                 "sha256": sha256_of(EXTRACT), "compact_of": [src["path"]], "n_rows": len(rows),
                 "supplementary": SPEC, "written_by": rel(Path(__file__)),
                 "note": ("compact extract: every row of the sheet with a classification; the "
                          "ACMG column (ACMG.class.AVENGERS) and the probability are not taken")}
        man["files"] = [e for e in man["files"] if e.get("path") != rel(EXTRACT)] + [entry]
        MANIFEST.write_text(json.dumps(man, indent=1, ensure_ascii=False) + "\n")
        log.output(MANIFEST)
        print(f"{rel(EXTRACT)}: {len(rows):,} rows")


if __name__ == "__main__":
    main()
