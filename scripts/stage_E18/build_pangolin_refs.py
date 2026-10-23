"""E18.2: Pangolin's gene database for the external genes' chromosomes outside the atlas's.

The atlas scored Pangolin on chromosomes 2, 3, 13, 16 and 17 with a gffutils database
built from the Ensembl release 112 GTF (models/pangolin/NOTES.md). The external genes
also lie on chromosomes 1 (SFPQ), 7 (XRCC2, SBDS) and X (DDX3X). This builds one
database per chromosome by the atlas's recipe, unchanged: the GTF rows of that
chromosome, then gffutils.create_db(subset, db, force=True, keep_order=True,
merge_strategy="create_unique", disable_infer_genes=True,
disable_infer_transcripts=True). The GTF is the atlas's own copy
(<atlas>/data/refs/Homo_sapiens.GRCh38.112.gtf.gz), its sha256 recorded. As a check
of the recipe, the chrX subset must equal, byte for byte, the one E2.7 built for
DDX3X (phase1/data/evidence/refs/chrX.gtf).

Writes phase1/data/e18_cache/refs/chr<N>.gtf and chr<N>.gtf.db (gitignored).

Run from the repository root, in the atlas's pinned Pangolin environment:
    <atlas>/models/pangolin/.venv/bin/python scripts/stage_E18/build_pangolin_refs.py
"""
from __future__ import annotations

import gzip
import hashlib
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

ATLAS = Path(os.environ.get("EVID_ATLAS_REPO", str(REPO.parent / "functional-standard-atlas")))
GTF = ATLAS / "data/refs/Homo_sapiens.GRCh38.112.gtf.gz"
OUT = REPO / "phase1/data/e18_cache/refs"
E27_CHRX = REPO / "phase1/data/evidence/refs/chrX.gtf"
CHROMS = ["1", "7", "X"]


def main() -> None:
    import gffutils
    with RunLog(__file__) as log:
        log.input(GTF)
        OUT.mkdir(parents=True, exist_ok=True)
        lines = {c: [] for c in CHROMS}
        with gzip.open(GTF, "rt") as f:
            for line in f:
                c = line.split("\t", 1)[0]
                if c in lines:
                    lines[c].append(line)
        for c in CHROMS:
            gtf = OUT / f"chr{c}.gtf"
            gtf.write_text("".join(lines[c]))
            db = OUT / f"chr{c}.gtf.db"
            gffutils.create_db(str(gtf), str(db), force=True, keep_order=True,
                               merge_strategy="create_unique", disable_infer_genes=True,
                               disable_infer_transcripts=True)
            log.output(gtf)
            log.output(db)
            print(f"chr{c}: {len(lines[c]):,} GTF rows -> {rel(db)}", flush=True)
        same = hashlib.sha256((OUT / "chrX.gtf").read_bytes()).hexdigest() == \
            hashlib.sha256(E27_CHRX.read_bytes()).hexdigest()
        log.note("chrX_subset_equals_E2.7", same)
        if not same:
            raise SystemExit("the chrX subset differs from E2.7's: the recipe is not the same")
        print("chrX subset: byte for byte E2.7's")


if __name__ == "__main__":
    main()
