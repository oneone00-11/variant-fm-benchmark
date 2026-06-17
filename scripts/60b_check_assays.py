"""Verify assay nature of borderline non-SGE candidate sets + confirm VHL/BAP1 SGE."""
import json
from config import DATA_RAW

CACHE = DATA_RAW / "mavedb" / "details"
for urn, tag in [("urn:mavedb:00000675-a-1", "VHL (SGE?)"),
                 ("urn:mavedb:00000662-0-1", "BAP1 (SGE?)"),
                 ("urn:mavedb:00000050-a-1", "MSH2"),
                 ("urn:mavedb:00001218-a-1", "MLH1"),
                 ("urn:mavedb:00000068-0-1", "TP53"),
                 ("urn:mavedb:00000102-0-1", "PTEN")]:
    p = CACHE / f"{urn.replace(':', '_')}.json"
    if not p.exists():
        print(f"{tag} {urn}: (not cached)\n"); continue
    d = json.loads(p.read_text(encoding="utf-8"))
    exp = d.get("experiment", {})
    ab = (exp.get("abstractText") or "")[:340]
    print(f"=== {tag}  {urn} | title: {d.get('title','')} ===")
    print(f"  {ab}\n")
