"""E20.3 (C13): the VCEP specifications of the seven training genes, as ClinGen's Criteria
Specification Registry (CSpec) serves them. Records only; no predictor score is read.

  1. the registry's list of every specification (api/svis), as served;
  2. for each training gene, every specification whose rule set names the gene (read
     from the list, not assumed), and for each: its document (api/SequenceVariant-
     Interpretation/id/<GN>) and the full record its assertionMethod points to (version,
     approval date, status history, release notes, every criterion code with its
     strengths, instructions to use and comments);
  3. for each specification the registry deposited on Zenodo (the full record's DOIs):
     the Zenodo record's metadata, which carries the licence (the deposit's files are not
     downloaded: they are the registry's own JSON and a PDF rendering of it);
  4. the terms of use the registry's pages link to (ClinGen's).

A response that is not HTTP 200 stops the run (E20.0 rule 4: nothing is filled in by hand).
Every file is written as served under phase1/data/evidence/vcep/raw/, and raw/manifest.json
records each one's URL, final URL, status, content type, size, sha256 and time, with the
genes searched and the specifications found for each.

Run from the repository root, with network access:
    .venv/bin/python scripts/stage_E20/fetch_vcep_specs.py
"""
from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, utc_now  # noqa: E402

RAW = REPO / "phase1/data/evidence/vcep/raw"
MANIFEST = RAW / "manifest.json"
GENES = ["BAP1", "BARD1", "BRCA1", "BRCA2", "PALB2", "RAD51C", "VHL"]
CSPEC = "https://cspec.genome.network/cspec"
LIST_URL = f"{CSPEC}/api/svis"
TERMS_URL = "https://www.clinicalgenome.org/docs/terms-of-use/"
ZENODO = "https://zenodo.org/api/records/"
AGENT = "variant-fm-benchmark E20 record fetch (urllib)"


def get(url: str) -> tuple[bytes, dict]:
    req = urllib.request.Request(url, headers={"User-Agent": AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        body = r.read()
        info = {"url": url, "final_url": r.geturl(), "http_status": r.status,
                "content_type": r.headers.get("Content-Type"), "fetched_utc": utc_now()}
    if info["http_status"] != 200:
        raise SystemExit(f"{url}: HTTP {info['http_status']}: stopped")
    return body, info


def save(name: str, url: str, files: list, log: RunLog) -> bytes:
    body, info = get(url)
    path = RAW / name
    path.write_bytes(body)
    log.output(path)
    files.append(info | {"file": rel(path), "bytes": len(body),
                         "sha256": hashlib.sha256(body).hexdigest()})
    return body


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    files: list = []
    with RunLog(__file__) as log:
        listing = json.loads(save("cspec_svis.json", LIST_URL, files, log))["data"]
        found = {g: [] for g in GENES}
        for entry in listing:
            genes = {x.get("label") for rs in entry.get("ruleSets", []) for x in rs.get("genes", [])}
            for g in sorted(genes & set(GENES)):
                found[g].append({"id": entry["@id"].rsplit("/", 1)[1], "status": entry.get("status"),
                                 "version": entry.get("version"),
                                 "vcep": entry["affiliation"]["label"], "url": entry.get("url")})
        specs = {}
        for g in GENES:
            for s in found[g]:
                specs[s["id"]] = s
        records = {}
        for gn in sorted(specs):
            doc = json.loads(save(f"{gn}.json", f"{CSPEC}/api/SequenceVariantInterpretation/id/{gn}",
                                  files, log))
            am = doc.get("assertionMethod", {}).get("url")
            if not am:
                raise SystemExit(f"{gn}: the document names no full record: stopped")
            rec = json.loads(save(f"{gn}_record.json", am, files, log))["data"]["entContent"]
            records[gn] = {"version": rec.get("version"), "approvedOn": rec.get("approvedOn"),
                           "doi": {k: (rec.get("doi") or {}).get(k) for k in ("conceptDoi", "docDoi")}}
        zenodo = {}
        for gn, r in sorted(records.items()):
            doi = r["doi"]["docDoi"] or r["doi"]["conceptDoi"]
            if not doi:
                continue
            recid = doi.rsplit(".", 1)[1]
            meta = json.loads(save(f"zenodo_{gn}_{recid}.json", ZENODO + recid, files, log))
            m = meta.get("metadata", {})
            zenodo[gn] = {"doi_followed": doi, "record": meta.get("id"), "version": m.get("version"),
                          "publication_date": m.get("publication_date"),
                          "license": (m.get("license") or {}).get("id"),
                          "files": [f.get("key") for f in meta.get("files", [])]}
        save("clingen_terms_of_use.html", TERMS_URL, files, log)
        manifest = {"what": "the VCEP specifications of the seven training genes in ClinGen's "
                            "Criteria Specification Registry, as served (E20.3)",
                    "written_by": rel(Path(__file__)), "genes_searched": GENES,
                    "rule": "a specification whose rule set names the gene, read from the "
                            "registry's list", "found": found, "records": records,
                    "zenodo": zenodo, "files": files}
        MANIFEST.write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
        log.output(MANIFEST)
        for g in GENES:
            print(g, [(s["id"], s["status"], s["version"]) for s in found[g]] or "none")


if __name__ == "__main__":
    main()
