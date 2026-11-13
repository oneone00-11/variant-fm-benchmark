"""E9.4: what stage 42 reads of VariantValidator's replies, taken from the replies as stored.

scripts/stage_E9/fetch_hgvs_check.py stored each reply as received, in a folder outside
the repository (--raw). This script reads them and writes, under
phase1/data/evidence/hgvs_check/:

    excerpt.csv     one row for each description submitted, in the order submitted (every
                    listed row on its RefSeq transcript, then the rows of SBDS on their
                    Ensembl transcript). Each field is copied from the reply's entry for
                    that description (found by its "submitted_variant") and nothing is
                    changed, compared or judged here:
                      submitted_description      the description as submitted
                      interface                  variantvalidator, or variantvalidator_ensembl
                      transcript_description     the entry's hgvs_transcript_variant
                      genomic_description_grch38 its primary_assembly_loci.grch38.hgvs_genomic_description
                      vcf_chr, vcf_pos, vcf_ref, vcf_alt     that locus's vcf
                      mane_select                its annotations.mane_select, as JSON prints it
                      flag                       the flag of the reply the entry came in
                      warnings                   its validation_warnings, the list as JSON
                    A field the entry does not hold, or holds as null, is left empty.
    manifest.json   where the replies came from: the interface's address, /hello/ as it
                    answered when the run started and when it ended, the run's first and last
                    request (UTC), the requests made and those made again, the descriptions,
                    the list's sha256, the second asking of 50 rows, and the bytes and sha256
                    of the archive of the raw replies.
    transcripts_DDX3X_ensembl.json, transcripts_SBDS_ensembl.json
                    the tool's own lists of the Ensembl transcripts it holds for the two
                    genes whose deposits are written on one, as received (why DDX3X's
                    descriptions on ENST00000644876.2 were not submitted).

and, beside --raw, the archive itself (--archive, not committed): every file of the raw
folder in a tar, the names sorted, each with the same owner, mode and time, compressed
with gzip without a name or a time in its header.

Stops, writing nothing, if a request's reply is not stored, if a reply lacks the entry of
a description, if the replies' versions are not all the same, or if the list is not the
one the run read (its sha256 in the journal).

    python scripts/stage_E9/extract_hgvs_check.py --raw <the folder of the replies> \\
        --archive <where the archive goes, outside the repository>
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
import tarfile
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _runlog import REPO, RunLog  # noqa: E402
import fetch_hgvs_check as F  # noqa: E402

LIST = REPO / "phase1" / "reports" / "evidence" / "layout" / "variants" / "analysed_variants.csv"
OUT = REPO / "phase1" / "data" / "evidence" / "hgvs_check"
COLUMNS = ["submitted_description", "interface", "transcript_description", "genomic_description_grch38",
           "vcf_chr", "vcf_pos", "vcf_ref", "vcf_alt", "mane_select", "flag", "warnings"]
NOT_ASKED = {
    "gene": "DDX3X", "transcript": "ENST00000644876.2",
    "why": "The deposits of DDX3X and SBDS write their c. descriptions on the Ensembl transcript of the gene's "
           "MANE Select pair. VariantValidator's interface for Ensembl transcripts holds SBDS's "
           "(ENST00000246868.7) and does not hold DDX3X's in any version: on 2026-11-12 three descriptions on "
           "ENST00000644876.2 came back \"The transcript ENST00000644876.2 is not in our database. Please check "
           "the transcript ID\", and the tool's list of the Ensembl transcripts it holds for DDX3X "
           "(transcripts_DDX3X_ensembl.json, as received) does not name it. DDX3X's descriptions were therefore "
           "submitted on the RefSeq transcript of the same MANE Select pair (NM_001356.5) only (the author's "
           "decision of 2026-11-12).",
}


def text(value) -> str:
    return "" if value is None else str(value)


def excerpt_row(item: dict, entry: dict, flag) -> dict:
    locus = (entry.get("primary_assembly_loci") or {}).get("grch38") or {}
    vcf = locus.get("vcf") or {}
    select = (entry.get("annotations") or {}).get("mane_select")
    return {"submitted_description": item["description"], "interface": item["interface"],
            "transcript_description": text(entry.get("hgvs_transcript_variant")),
            "genomic_description_grch38": text(locus.get("hgvs_genomic_description")),
            "vcf_chr": text(vcf.get("chr")), "vcf_pos": text(vcf.get("pos")),
            "vcf_ref": text(vcf.get("ref")), "vcf_alt": text(vcf.get("alt")),
            "mane_select": "" if select is None else json.dumps(select), "flag": text(flag),
            "warnings": json.dumps(entry.get("validation_warnings") or [], ensure_ascii=False)}


def read_replies(raw: Path, rows: list[dict]) -> tuple[list[dict], dict]:
    """The excerpt's rows, and the versions every reply gives."""
    out, seen = [], None
    for req in F.requests_of(F.descriptions(rows)):
        path = F.reply_path(raw, req)
        if not path.is_file():
            raise SystemExit(f"no reply is stored for rows {req[0]['row']}-{req[-1]['row']} ({req[0]['interface']})")
        body = path.read_bytes()
        problem = F.whole(body, req)
        if problem:
            raise SystemExit(f"{path.relative_to(raw)}: {problem}")
        versions = F.versions(body)
        seen = seen or versions
        if versions != seen:
            raise SystemExit(f"{path.relative_to(raw)}: its versions {versions} are not the first reply's {seen}")
        found, flag = F.entries(body), json.loads(body).get("flag")
        out += [excerpt_row(it, found[it["description"]], flag) for it in req]
    return out, seen


def journal_facts(raw: Path) -> dict:
    """What the run's journal records: the list it read, its requests, when."""
    lines = [json.loads(ln) for ln in (raw / "journal.jsonl").read_text().splitlines() if ln.strip()]
    requests = [r for r in lines if r.get("what") == "request"]
    stored = [r for r in requests if r.get("outcome") == "stored"]
    failed = [r for r in requests if r.get("outcome") != "stored"]
    by_status: dict[str, int] = {}
    for r in failed:
        by_status[str(r["status"])] = by_status.get(str(r["status"]), 0) + 1
    starts = [r for r in lines if r.get("what") == "run starts"]
    rechecks = [r for r in lines if r.get("what") == "recheck"]
    return {"list_sha256_the_run_read": sorted({r["list_sha256"] for r in starts}),
            "runs_started": len(starts), "first_request_utc": requests[0]["utc"], "last_request_utc": requests[-1]["utc"],
            "requests_made": len(requests), "requests_stored": len(stored), "requests_made_again": len(failed),
            "made_again_by_status": dict(sorted(by_status.items())),
            "pauses_of_ten_minutes": sum(1 for r in lines if r.get("what") == "pause"),
            "second_asking_requests": len(rechecks),
            "second_asking_last_utc": rechecks[-1]["utc"] if rechecks else None}


def write_archive(raw: Path, dest: Path) -> dict:
    """Every file under `raw` in a tar.gz that depends on the files alone: names sorted,
    owner 0, mode 0644, time 0; gzip's header without a name or a time."""
    files = sorted(p for p in raw.rglob("*") if p.is_file() and p.name != ".DS_Store")
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w", format=tarfile.PAX_FORMAT) as tar:
        for p in files:
            info = tarfile.TarInfo(name=f"hgvs_raw/{p.relative_to(raw).as_posix()}")
            data = p.read_bytes()
            info.size, info.mtime, info.mode, info.uid, info.gid, info.uname, info.gname = len(data), 0, 0o644, 0, 0, "", ""
            tar.addfile(info, io.BytesIO(data))
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as f, gzip.GzipFile(filename="", mode="wb", fileobj=f, mtime=0, compresslevel=9) as z:
        z.write(buffer.getvalue())
    return {"name": dest.name, "bytes": dest.stat().st_size, "sha256": hashlib.sha256(dest.read_bytes()).hexdigest(),
            "files": len(files), "tar_bytes": len(buffer.getvalue()), "tar_sha256": hashlib.sha256(buffer.getvalue()).hexdigest()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=Path)
    ap.add_argument("--archive", required=True, type=Path)
    args = ap.parse_args()
    raw, archive = args.raw.expanduser().resolve(), args.archive.expanduser().resolve()
    with RunLog(__file__) as log:
        log.input(LIST)
        rows = F.read_list(LIST)
        items = F.descriptions(rows)
        if len({(it["interface"], it["description"]) for it in items}) != len(items):
            raise SystemExit("a description is submitted twice to one interface")
        excerpt, versions = read_replies(raw, rows)
        facts = journal_facts(raw)
        list_sha = hashlib.sha256(LIST.read_bytes()).hexdigest()
        if facts["list_sha256_the_run_read"] != [list_sha]:
            raise SystemExit(f"the run read a list of sha256 {facts['list_sha256_the_run_read']}, the list is {list_sha}")
        hellos = sorted(raw.glob("hello_*.json"))
        if len(hellos) < 2:
            raise SystemExit("/hello/ is stored fewer than twice: the run has not ended")
        comparison = raw / "recheck" / "comparison.json"
        if not comparison.is_file():
            raise SystemExit("the second asking of 50 rows has not been made (fetch_hgvs_check.py recheck)")
        again = json.loads(comparison.read_text())

        OUT.mkdir(parents=True, exist_ok=True)
        with open(OUT / "excerpt.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator="\n")
            w.writeheader()
            w.writerows(excerpt)
        for gene in F.DEPOSITS_ON_ENSEMBL:
            (OUT / f"transcripts_{gene}_ensembl.json").write_bytes((raw / f"transcripts_{gene}_ensembl.json").read_bytes())
        packed = write_archive(raw, archive)
        manifest = {
            "what": "VariantValidator's replies to the HGVS descriptions of every listed variant (E9.4): where "
                    "excerpt.csv comes from. Written by scripts/stage_E9/extract_hgvs_check.py from the replies "
                    "scripts/stage_E9/fetch_hgvs_check.py stored.",
            "tool": "VariantValidator",
            "interface": {"address": F.BASE,
                          "paths": [f"/VariantValidator/{F.REFSEQ}/{F.BUILD}/<descriptions>/<transcript>",
                                    f"/VariantValidator/{F.ENSEMBL}/{F.BUILD}/<descriptions>/<transcript>"],
                          "genome_build": F.BUILD, "descriptions_per_request_at_most": F.PER_REQUEST,
                          "seconds_between_requests_at_least": F.PAUSE},
            "hello_at_start": {"file": hellos[0].name, "text": hellos[0].read_text()},
            "hello_at_end": {"file": hellos[-1].name, "text": hellos[-1].read_text()},
            "hello_files": [p.name for p in hellos],
            "versions_in_every_reply": versions,
            "list": {"path": LIST.relative_to(REPO).as_posix(), "rows": len(rows), "sha256": list_sha},
            "descriptions": {"submitted": len(items),
                             "by_interface": {k: sum(1 for it in items if it["interface"] == k) for k in (F.REFSEQ, F.ENSEMBL)},
                             "requests": len(F.requests_of(items))},
            "run": facts,
            "second_asking": {"seed": again["seed"], "rows": len(again["rows"]), "entries": len(again["entries"]),
                              "entries_the_same_as_stored": sum(1 for e in again["entries"] if e["same_entry"]),
                              "entries_that_differ": [e for e in again["entries"] if not e["same_entry"]]},
            "not_submitted": NOT_ASKED,
            "excerpt": {"path": (OUT / "excerpt.csv").relative_to(REPO).as_posix(), "rows": len(excerpt),
                        "sha256": hashlib.sha256((OUT / "excerpt.csv").read_bytes()).hexdigest(), "columns": COLUMNS},
            "raw_replies_archive": packed | {"kept": "outside the repository; to go into the release's archive of inputs"},
        }
        (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
        for name in ("excerpt.csv", "manifest.json", "transcripts_DDX3X_ensembl.json", "transcripts_SBDS_ensembl.json"):
            log.output(OUT / name)
        log.note("archive", packed)
        log.note("rows", len(excerpt))
        print(f"{len(excerpt)} rows in {(OUT / 'excerpt.csv').relative_to(REPO)}; archive {packed['name']}: "
              f"{packed['bytes']:,} bytes, sha256 {packed['sha256']}")


if __name__ == "__main__":
    main()
