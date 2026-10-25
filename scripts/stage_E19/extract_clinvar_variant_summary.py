"""E19 C1 (E19.4): the analysis set's rows of ClinVar's variant_summary archive.

The archive (variant_summary_2026-06.txt.gz, 439 MB) is not committed (E19.0 rule 4);
scripts/stage_E19/fetch_clinvar_variant_summary.py keeps it in the gitignored cache and
records its sha256. This writes the small table committed in its place: every row of
the archive on the GRCh38 assembly whose VariationID is one of the analysis set's, with
every column as served (text, unchanged). The archive is checked against the sha256
its manifest records before a byte is read.

Run in the reference image (the extract is committed; 修改单 §0 rule 9), offline:
    docker run ... IMAGE bash scripts/stage_E19/in_image.sh scripts/stage_E19/extract_clinvar_variant_summary.py
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of  # noqa: E402

OUT = REPO / "phase1/data/evidence/e19/clinvar_variant_summary"
MANIFEST = OUT / "manifest.json"
EXTRACT = OUT / "variant_summary_2026-06_analysis_set.tsv"
SET = REPO / "phase1/data/evidence/analysis_set_v1.parquet"


def main() -> None:
    import pandas as pd
    with RunLog(__file__) as log:
        man = json.loads(MANIFEST.read_text())
        archive = REPO / man["cached_at"]
        log.input(MANIFEST)
        if sha256_of(archive) != man["sha256"]:
            raise SystemExit(f"{rel(archive)} is not the file its manifest records: stopped")
        log.input(SET)
        ids = {str(v) for v in pd.read_parquet(SET, columns=["variation_id"])["variation_id"].dropna()}
        n_rows = n_kept = 0
        with gzip.open(archive, "rt", encoding="utf-8", newline="") as f, \
                open(EXTRACT, "w", encoding="utf-8", newline="") as g:
            header = f.readline().rstrip("\n").split("\t")
            if header != man["columns"]:
                raise SystemExit("the archive's columns are not the ones its manifest records")
            vid, asm = header.index("VariationID"), header.index("Assembly")
            g.write("\t".join(header) + "\n")
            for line in f:
                n_rows += 1
                p = line.rstrip("\n").split("\t")
                if p[asm] == "GRCh38" and p[vid] in ids:
                    g.write(line if line.endswith("\n") else line + "\n")
                    n_kept += 1
        log.output(EXTRACT)
        found = set(pd.read_csv(EXTRACT, sep="\t", dtype=str, keep_default_na=False,
                                quoting=csv.QUOTE_NONE)["VariationID"])
        summary = {"archive_rows": n_rows, "rows_kept": n_kept,
                   "variation_ids_in_set": len(ids),
                   "variation_ids_found": len(found & ids),
                   "variation_ids_not_found": sorted(ids - found)}
        man["extract"] = {"file": rel(EXTRACT), "sha256": sha256_of(EXTRACT),
                          "bytes": EXTRACT.stat().st_size, "written_by": rel(Path(__file__)),
                          "selection": "Assembly == GRCh38 and VariationID in the analysis set",
                          **{k: v for k, v in summary.items() if k != "variation_ids_not_found"},
                          "variation_ids_not_found": summary["variation_ids_not_found"]}
        MANIFEST.write_text(json.dumps(man, indent=1) + "\n")
        log.output(MANIFEST)
        log.note("summary", {k: v for k, v in summary.items() if k != "variation_ids_not_found"})
        print(f"kept {n_kept:,} of {n_rows:,} rows; {len(found & ids):,} of {len(ids):,} "
              "VariationIDs found")


if __name__ == "__main__":
    main()
