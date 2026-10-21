"""E16c: fetch, as served, the MaveDB score set that each training gene's labels cite
(supplementary/Table_S1.tsv), with every version of its chain, for the read-only label
source audit (phase1/src/evid_label_audit.py). Nothing here touches a label file, the
configuration, the analysis set or the manuscript.

For each of the seven genes this writes, under phase1/data/evidence/label_audit/<gene>_raw/:
  scoreset_meta_<id>.json   each version's metadata, as served
  urn_mavedb_<id>.csv       each version's score table, byte for byte as served
  score_ranges_<id>.json    each version's scoreCalibrations, the investigator-provided
                            ones selected, entries unchanged
  experiment_meta.json      the experiment of the cited version
  calibration_variants_<calibration id>[_<class id>].json
                            for an investigator-provided calibration whose classes are
                            classes rather than score ranges, the per-variant membership
                            MaveDB serves (/score-calibrations/{urn}/variants, or one
                            response per class when the whole is above the size limit);
                            above LARGE_BYTES it goes to the gitignored cache instead
                            (the author's decision of 2026-10-21), as served
  calibration_classes_<calibration id>.tsv
                            the same membership as a compact extract (variant -> class),
                            committed; the audit reads this
  supplementary_table3_categories.tsv (BRCA2)
                            the c. HGVS, functional category and score of every classed
                            row of the paper's Supplementary Table 3, committed; the xlsx
                            itself, fetched from Europe PMC's supplementary files for
                            the paper, goes to the cache and must have the sha256 E16c
                            recorded (E16d.1)
  manifest.json             every file's path, size, sha256, source URL, fetch time and
                            API version; the cited URN and where it is cited; the version
                            chain, with each score set's and each calibration's own
                            creation and modification dates (not_provided where a field is
                            absent); the licence; any request that failed

The cache is phase1/data/evidence/label_audit_cache/<gene>/ (gitignored); each cached
file is listed in the manifest with "cached": true, its URL and sha256, and
scripts/reproduce_evidence.py --fetch-inputs fetches it again on another machine.

BARD1 and PALB2: E16b fetched every version of both chains already
(phase1/data/evidence/external/label_versions/). Those files are reused where they
are, listed in the manifest with the fetch record E16b wrote, and not downloaded
again; only the experiment metadata is fetched here.

The network is used here only; the audit reads what this writes.

Run from the repository root:
    .venv/bin/python scripts/stage_E16c/fetch_sources.py              # all seven genes
    .venv/bin/python scripts/stage_E16c/fetch_sources.py --gene BRCA1
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import sys
import zipfile
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

sys.path.insert(0, str(REPO / "phase1"))
from src import evid_label_sources as LS  # noqa: E402

BASE = "https://api.mavedb.org/api/v1"
TABLE_S1 = REPO / "supplementary" / "Table_S1.tsv"
OUT = REPO / "phase1" / "data" / "evidence" / "label_audit"
E16B = REPO / "phase1" / "data" / "evidence" / "external" / "label_versions"
CACHE = REPO / "phase1" / "data" / "evidence" / "label_audit_cache"
# a file above this size is not committed as served (E16.2); the script stops instead
MAX_BYTES = 25 * 1024 * 1024
# a calibration membership above this size is kept in the cache, not committed (BAP1's
# one and RAD51C's two, 9-18 MB; the author's decision of 2026-10-21); VHL's 1.4 MB is
# committed as served
LARGE_BYTES = 5 * 1024 * 1024
NOT_PROVIDED = "not_provided"
# The paper supplementary table the author made BRCA2's label source on 2026-10-21
# (E16d.1, rule R6): Europe PMC serves a paper's supplementary files as one zip; the
# xlsx in it must be the file E16c read (sha256 below).
SUPPLEMENTARY = {
    "BRCA2": {"pmcid": "PMC11821525", "member": "41586_2024_8388_MOESM3_ESM.xlsx",
              "sha256": "020f3d7acc04357410c0217447bbfcbc5bd3562e7aead39ae5726d16e8da61c0",
              "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11821525/supplementaryFiles",
              "sheet": "Table S3", "header_row": 2,
              "hgvs_column": "Coding sequence change (c.)", "class_column": "Functional category",
              "score_column": "Model based functional score",
              "compact": "supplementary_table3_categories.tsv"},
}

S = requests.Session()
S.trust_env = False


def _get(url: str) -> requests.Response:
    return S.get(url, timeout=600)


def _id(urn: str) -> str:
    return urn.split(":")[-1]


def _save(path: Path, blob: bytes, log: RunLog) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(blob)
    log.output(path)
    if len(blob) > MAX_BYTES:
        raise SystemExit(f"{rel(path)} is {len(blob):,} bytes, above the {MAX_BYTES:,}-byte "
                         "limit for committing a file as served; stopped (E16.2)")


def _entry(path: Path, url: str, version: str, status: int) -> dict:
    return {"path": rel(path), "bytes": path.stat().st_size, "sha256": sha256_of(path),
            "source_url": url, "fetched_utc": utc_now(), "api_version": version,
            "http_status": status}


def _ranges_extract(meta: dict, source: Path) -> dict:
    """As E16 extracts them (scripts/stage_E16/fetch_mavedb.py): entries unchanged."""
    cals = meta.get("scoreCalibrations") or []
    return {"extracted_from": rel(source), "key": "scoreCalibrations",
            "selection": "investigatorProvided is true",
            "all_calibration_titles": [c.get("title") for c in cals],
            "legacy_scoreRanges_key_present": meta.get("scoreRanges") is not None,
            "investigator_provided": [c for c in cals if c.get("investigatorProvided") is True]}


def _dates_inside(obj, found: list) -> list:
    """Every creationDate / modificationDate anywhere inside a calibration (its
    classifications carry their own)."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in ("creationDate", "modificationDate") and isinstance(v, str):
                found.append(v)
            else:
                _dates_inside(v, found)
    elif isinstance(obj, list):
        for v in obj:
            _dates_inside(v, found)
    return found


def version_record(meta: dict, n_rows: int, cited: str) -> dict:
    lic = meta.get("license") or {}
    cals = []
    for c in meta.get("scoreCalibrations") or []:
        inside = _dates_inside(c.get("functionalClassifications") or [], [])
        cals.append({
            "title": c.get("title", NOT_PROVIDED),
            "investigator_provided": c.get("investigatorProvided", NOT_PROVIDED),
            "research_use_only": c.get("researchUseOnly", NOT_PROVIDED),
            "created": c.get("creationDate", NOT_PROVIDED),
            "modified": c.get("modificationDate", NOT_PROVIDED),
            "n_classes": len(c.get("functionalClassifications") or []),
            "latest_date_inside_its_classes": max(inside) if inside else NOT_PROVIDED,
        })
    return {
        "urn": meta["urn"], "role": "cited" if meta["urn"] == cited else "other version",
        "created": meta.get("creationDate", NOT_PROVIDED),
        "modified": meta.get("modificationDate", NOT_PROVIDED),
        "published": meta.get("publishedDate", NOT_PROVIDED),
        "supersedes": (meta.get("supersededScoreSet") or {}).get("urn"),
        "superseded_by": (meta.get("supersedingScoreSet") or {}).get("urn"),
        "num_variants_metadata": meta.get("numVariants", NOT_PROVIDED),
        "n_rows_score_table": n_rows,
        "license": {k: lic.get(k, NOT_PROVIDED) for k in ("shortName", "version", "link")},
        "private": meta.get("private", NOT_PROVIDED),
        "score_calibrations": cals,
    }


def chain(urn: str) -> tuple[list[requests.Response], list[dict]]:
    """Every version of a score set's chain, oldest first: back through
    supersededScoreSet and forward through supersedingScoreSet. A linked version the
    API does not serve is recorded (URL, status, body) and the walk stops there in
    that direction; the cited version itself must be served."""
    seen, back, missing = {}, [], []

    def get(u: str) -> requests.Response | None:
        r = _get(f"{BASE}/score-sets/{u}")
        if r.status_code != 200 and u != urn:
            missing.append({"urn": u, "url": r.url, "http_status": r.status_code,
                            "body": r.text[:500], "at_utc": utc_now()})
            return None
        r.raise_for_status()
        return r

    cur = urn
    while cur and cur not in seen:
        r = get(cur)
        if r is None:
            break
        seen[cur] = r
        back.append(cur)
        cur = (r.json().get("supersededScoreSet") or {}).get("urn")
    fwd, cur = [], (seen[urn].json().get("supersedingScoreSet") or {}).get("urn")
    while cur and cur not in seen:
        r = get(cur)
        if r is None:
            break
        seen[cur] = r
        fwd.append(cur)
        cur = (r.json().get("supersedingScoreSet") or {}).get("urn")
    return [seen[u] for u in list(reversed(back)) + fwd], missing


def _rows(blob: bytes) -> int:
    return sum(1 for _ in csv.DictReader(io.StringIO(blob.decode("utf-8"), newline="")))


def _class_based(cal: dict) -> bool:
    fcs = cal.get("functionalClassifications") or []
    return bool(fcs) and all(fc.get("range") is None and fc.get("class") is not None
                             for fc in fcs)


def fetch_class_members(meta: dict, gene: str, out: Path, version: str, log: RunLog) -> list[dict]:
    """The per-variant membership of each investigator-provided calibration whose
    classes are given as classes rather than score ranges (no score places a variant
    in them): /score-calibrations/{urn}/variants as served or, when that response is
    above the size limit, one response per class; in the cache when above LARGE_BYTES.
    Beside it, committed, a compact extract of the same membership."""
    files = []
    for cal in meta.get("scoreCalibrations") or []:
        if cal.get("investigatorProvided") is not True or not _class_based(cal):
            continue
        cid = cal["urn"].split("calibration-")[-1]
        r = _get(f"{BASE}/score-calibrations/{cal['urn']}/variants")
        r.raise_for_status()
        parts = [(r, f"calibration_variants_{cid}.json", "all")]
        if len(r.content) > MAX_BYTES:
            parts = []
            for fc in cal["functionalClassifications"]:
                rr = _get(f"{BASE}/score-calibrations/{cal['urn']}/functional-classifications/"
                          f"{fc['id']}/variants")
                rr.raise_for_status()
                parts.append((rr, f"calibration_variants_{cid}_{fc['id']}.json", fc["id"]))
        cached = sum(len(resp.content) for resp, _, _ in parts) > LARGE_BYTES
        raws = []
        for resp, name, fcid in parts:
            path = (CACHE / gene.lower() / name) if cached else (out / name)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(resp.content)
            log.output(path)
            if not cached and len(resp.content) > MAX_BYTES:
                raise SystemExit(f"{rel(path)} is above {MAX_BYTES:,} bytes; stopped (E16.2)")
            raws.append(_entry(path, resp.url, version, resp.status_code)
                        | {"score_set": meta["urn"], "calibration_urn": cal["urn"],
                           "calibration_title": cal.get("title"),
                           "functional_classification_id": fcid, "cached": cached})
        compact = out / f"calibration_classes_{cid}.tsv"
        rows = LS.membership_rows([json.loads(resp.content) for resp, _, _ in parts], cal)
        LS.write_tsv(compact, LS.MEMBERSHIP_COLUMNS, rows)
        log.output(compact)
        files += raws + [{"path": rel(compact), "bytes": compact.stat().st_size,
                          "sha256": sha256_of(compact), "score_set": meta["urn"],
                          "calibration_urn": cal["urn"], "calibration_title": cal.get("title"),
                          "compact_of": [e["path"] for e in raws], "n_rows": len(rows),
                          "note": "compact extract of the served membership: variant URN, "
                                  "hgvs_nt, functional classification id and its label"}]
    return files


def fetch_supplementary(gene: str, out: Path, log: RunLog) -> list[dict]:
    """A paper supplementary table the author made a label source (SUPPLEMENTARY): the
    zip Europe PMC serves for the paper, the table taken out of it into the cache and
    checked against the sha256 E16c recorded, and a compact extract committed."""
    spec = SUPPLEMENTARY.get(gene)
    if spec is None:
        return []
    r = _get(spec["url"])
    r.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(r.content)) as z:
        blob = z.read(spec["member"])
    sha = hashlib.sha256(blob).hexdigest()
    if sha != spec["sha256"]:
        raise SystemExit(f"{gene}: {spec['member']} from {spec['url']} has sha256 {sha}, not the "
                         f"{spec['sha256']} E16c recorded; stopped (E16d.6)")
    path = CACHE / gene.lower() / spec["member"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(blob)
    log.output(path)
    rows, layout = LS.supplementary_rows(path, spec)
    compact = out / spec["compact"]
    LS.write_tsv(compact, LS.SUPPLEMENTARY_COLUMNS, rows)
    log.output(compact)
    raw = {"path": rel(path), "bytes": len(blob), "sha256": sha, "source_url": spec["url"],
           "zip_member": spec["member"], "zip_bytes": len(r.content), "fetched_utc": utc_now(),
           "http_status": r.status_code, "cached": True, "pmcid": spec["pmcid"],
           "sheet": spec["sheet"], "title_cell": layout["title_cell"],
           "header_row": layout["header"]}
    return [raw, {"path": rel(compact), "bytes": compact.stat().st_size,
                  "sha256": sha256_of(compact), "compact_of": [raw["path"]], "n_rows": len(rows),
                  "supplementary": {k: spec[k] for k in ("pmcid", "member", "sheet", "header_row",
                                                       "hgvs_column", "class_column",
                                                       "score_column")},
                  "note": "compact extract: every row of the sheet with a functional category"}]


def fetch_gene(gene: str, cited: dict, label_file: str, version: str, log: RunLog) -> dict:
    out = OUT / f"{gene.lower()}_raw"
    files, versions, missing = [], [], []
    e16b = E16B / f"{gene.lower()}_raw" / "manifest.json"
    if e16b.is_file():
        log.input(e16b)
        man = json.loads(e16b.read_text())
        if cited["urn"] not in [v["urn"] for v in man["versions"]]:
            raise SystemExit(f"{gene}: E16b's chain lacks the cited {cited['urn']}")
        by_name = {Path(e["path"]).name: e for e in man["files"]}
        for v in man["versions"]:
            vid = _id(v["urn"])
            for name in (f"scoreset_meta_{vid}.json", f"urn_mavedb_{vid}.csv",
                         f"score_ranges_{vid}.json"):
                e = by_name[name]
                path = REPO / e["path"]
                if sha256_of(path) != e["sha256"]:
                    raise SystemExit(f"{e['path']} differs from E16b's manifest")
                log.input(path)
                files.append(e | {"reused_from": "E16b", "e16b_manifest": rel(e16b)})
            meta = json.loads((REPO / by_name[f"scoreset_meta_{vid}.json"]["path"]).read_text())
            n = _rows((REPO / by_name[f"urn_mavedb_{vid}.csv"]["path"]).read_bytes())
            versions.append(version_record(meta, n, cited["urn"]))
            files += fetch_class_members(meta, gene, out, version, log)
        cited_meta = json.loads((REPO / by_name[f"scoreset_meta_{_id(cited['urn'])}.json"]
                                 ["path"]).read_text())
    else:
        cited_meta = None
        served, missing = chain(cited["urn"])
        for resp in served:
            meta, vid = resp.json(), _id(resp.json()["urn"])
            mpath = out / f"scoreset_meta_{vid}.json"
            _save(mpath, resp.content, log)
            files.append(_entry(mpath, resp.url, version, resp.status_code))
            r = _get(f"{BASE}/score-sets/{meta['urn']}/scores")
            r.raise_for_status()
            spath = out / f"urn_mavedb_{vid}.csv"
            _save(spath, r.content, log)
            files.append(_entry(spath, r.url, version, r.status_code))
            rpath = out / f"score_ranges_{vid}.json"
            _save(rpath, (json.dumps(_ranges_extract(meta, mpath), indent=1,
                                     ensure_ascii=False) + "\n").encode(), log)
            files.append({"path": rel(rpath), "bytes": rpath.stat().st_size,
                          "sha256": sha256_of(rpath), "source_url": resp.url,
                          "derived_from": rel(mpath),
                          "note": "extract of scoreCalibrations, entries unchanged"})
            versions.append(version_record(meta, _rows(r.content), cited["urn"]))
            files += fetch_class_members(meta, gene, out, version, log)
            if meta["urn"] == cited["urn"]:
                cited_meta = meta
    exp_urn = cited_meta["experiment"]["urn"]
    r = _get(f"{BASE}/experiments/{exp_urn}")
    r.raise_for_status()
    epath = out / "experiment_meta.json"
    _save(epath, r.content, log)
    files.append(_entry(epath, r.url, version, r.status_code))
    files += fetch_supplementary(gene, out, log)
    manifest = {
        "gene": gene, "label_file": f"phase1/data/assay_labels/{label_file}",
        "cited": {"urn": cited["urn"], "where": rel(TABLE_S1),
                  "license_in_table_s1": cited["license"]},
        "experiment_urn": exp_urn,
        "api": {"base": BASE, "version": version},
        "versions": versions,
        "versions_not_served": missing,
        "cache_dir": rel(CACHE / gene.lower()),
        "files": files,
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
    log.output(out / "manifest.json")
    print(f"{gene}: cited {cited['urn']}; chain " + ", ".join(
        f"{v['urn'].split('-', 1)[-1]} ({v['n_rows_score_table']:,} rows, "
        f"{sum(1 for c in v['score_calibrations'] if c['investigator_provided'] is True)} "
        f"investigator calibration(s))" for v in versions))
    return {"gene": gene, "cited": cited["urn"], "versions": [v["urn"] for v in versions]}


def main() -> None:
    import argparse
    sys.path.insert(0, str(REPO / "phase1"))
    from src.config import ASSAY_LABELS
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", nargs="+", choices=list(ASSAY_LABELS),
                    help="fetch these genes only (default: all seven)")
    genes = ap.parse_args().gene or list(ASSAY_LABELS)
    with RunLog(__file__) as log:
        log.input(TABLE_S1)
        with open(TABLE_S1, newline="") as f:
            s1 = {r["gene"]: r for r in csv.DictReader(f, delimiter="\t")}
        r = _get(f"{BASE}/api/version")
        r.raise_for_status()
        version = r.json()["version"]
        log.note("api_version", version)
        done, failures = [], []
        for gene in genes:
            spec = ASSAY_LABELS[gene]
            try:
                done.append(fetch_gene(gene, {"urn": s1[gene]["mavedb_urn"], **s1[gene]},
                                       spec["file"], version, log))
            except requests.RequestException as e:
                failures.append({"gene": gene, "urn": s1[gene]["mavedb_urn"], "error": str(e),
                                 "at_utc": utc_now()})
                print(f"{gene}: request failed ({e}); recorded")
        log.note("genes", done)
        log.note("failures", failures)


if __name__ == "__main__":
    main()
