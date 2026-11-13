"""E9.4: every listed variant's HGVS description, submitted to VariantValidator's public REST
interface, each reply stored as received.

The list is phase1/reports/evidence/layout/variants/analysed_variants.csv (stage 41). Each
row is submitted as its RefSeq transcript and its c. description, "NM_...:c....", on GRCh38,
with that transcript as the only one to report. The deposits of DDX3X and SBDS describe their
variants on an Ensembl transcript. The rows of SBDS are therefore submitted a second time as
"ENST...:c....", to the interface for Ensembl transcripts. The rows of DDX3X are not: the
interface does not hold ENST00000644876 in any version (on 2026-11-12 three of its
descriptions came back "The transcript ENST00000644876.2 is not in our database", and the
tool's own list of the Ensembl transcripts it holds for the gene, which a run stores, names
19 others), so DDX3X is checked on the RefSeq transcript of the same MANE Select pair only
(the author's decision of 2026-11-12). Nothing is interpreted here:
scripts/stage_E9/extract_hgvs_check.py reads the stored replies afterwards.

    python scripts/stage_E9/fetch_hgvs_check.py run     --list <the list> --raw <a folder outside the repository>
    python scripts/stage_E9/fetch_hgvs_check.py recheck --list <the list> --raw <the same folder>

run      Requests go one at a time, never in parallel: at most 10 descriptions a request,
         all on one transcript, joined by "|", and at least one second between two requests
         (the interface's notes recommend no more than 3 requests a second). A reply counts
         only if it is JSON holding an entry for every description submitted (found by the
         entry's "submitted_variant", never by its key: where the reference base does not
         match the genome the tool rewrites the description and keys the entry by the
         rewritten one). On 429 or 5xx, on a network error, or on a reply that is not that,
         the request is tried again after the reply's Retry-After, or else after a wait that
         doubles from 5 seconds; five failures in a row pause the run for 10 minutes, and a
         third such round ends it (exit 3). Any other status ends it (exit 4). The run also
         ends if the software or data versions change on the way (exit 5): each reply's
         "metadata" must be that of the first reply stored, and /hello/ must answer as it
         did the first time it was asked (a reply writes one of the five versions with a
         suffix /hello/ does not print, "VV_SR_2025_02/master", so replies are compared
         with replies and /hello/ with /hello/). It ends, too, if 30 requests in the whole
         would take more than 24 hours (exit 6). It can be started again at any time: a
         request whose reply is already stored is not made again.
recheck  After the run: 50 rows drawn with a fixed seed are submitted again, one row a
         request, and each entry is compared with the stored one.

Stored under --raw, nothing else written anywhere:
    hello_<UTC>.json                  /hello/ as received, when a run starts and when it ends
    transcripts_<gene>_ensembl.json   the Ensembl transcripts the tool holds for DDX3X and for SBDS
                                      (tools/gene2transcripts_v2/<gene>/all/ensembl/GRCh38), as received
    replies/<interface>/<transcript>/<first row>_<last row>.json     each reply as received
    journal.jsonl                     one line per request made: when, what was submitted, the
                                      status, the bytes and sha256 of the reply, the seconds it
                                      took, and what was done with it
    recheck/                          the second replies of the 50 rows and the comparison

Requests carry the descriptions and nothing else: no key, no name, no e-mail address; the
User-Agent is this repository's address. Not a stage: scripts/reproduce_evidence.py does not
run it, and stage 42 reads only what extract_hgvs_check.py wrote from these replies.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE = "https://rest.variantvalidator.org"
USER_AGENT = "https://github.com/oneone00-11/variant-fm-benchmark"
BUILD = "GRCh38"
REFSEQ, ENSEMBL = "variantvalidator", "variantvalidator_ensembl"       # the two interfaces
ENSEMBL_GENES = ("SBDS",)               # its deposit writes c. on the Ensembl transcript, which the tool holds
DEPOSITS_ON_ENSEMBL = ("DDX3X", "SBDS")  # the tool's list of Ensembl transcripts is stored for both
PER_REQUEST = 10
PAUSE = 1.0                 # seconds between two requests, at least
FIRST_WAIT, LONGEST_WAIT = 5.0, 300.0
IN_A_ROW, ROUND_PAUSE, ROUNDS = 5, 600.0, 3
TIMEOUT = 600.0
LIMIT_HOURS, JUDGE_AFTER = 24.0, 30
RECHECK_ROWS, RECHECK_SEED = 50, 20261112
VERSION_KEYS = ("variantvalidator_hgvs_version", "variantvalidator_version", "vvdb_version",
                "vvseqrepo_db", "vvta_version")


def utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_list(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def descriptions(rows: list[dict]) -> list[dict]:
    """What is submitted, in the list's order: every row on its RefSeq transcript, then the
    rows of SBDS again on their Ensembl transcript."""
    out = []
    for interface, column, genes in ((REFSEQ, "RefSeq transcript (MANE Select)", None),
                                     (ENSEMBL, "Ensembl transcript (MANE Select)", ENSEMBL_GENES)):
        for n, r in enumerate(rows, start=1):
            if genes is None or r["Gene"] in genes:
                out.append({"row": n, "gene": r["Gene"], "interface": interface,
                            "transcript": r[column], "description": f'{r[column]}:{r["HGVS c."]}'})
    return out


def requests_of(items: list[dict], size: int = PER_REQUEST) -> list[list[dict]]:
    """The items cut into requests: one interface and one transcript each, in order."""
    out: list[list[dict]] = []
    for it in items:
        if out and len(out[-1]) < size and (out[-1][0]["interface"], out[-1][0]["transcript"]) == (it["interface"], it["transcript"]):
            out[-1].append(it)
        else:
            out.append([it])
    return out


def url_of(req: list[dict], base: str = BASE) -> str:
    """The descriptions joined by "|", each character outside the unreserved ones percent-encoded
    ("+" %2B, ">" %3E, "|" %7C, ":" %3A)."""
    joined = urllib.parse.quote("|".join(it["description"] for it in req), safe="")
    return (f"{base}/VariantValidator/{req[0]['interface']}/{BUILD}/{joined}/"
            f"{urllib.parse.quote(req[0]['transcript'], safe='')}?content-type=application%2Fjson")


def reply_path(raw: Path, req: list[dict]) -> Path:
    return raw / "replies" / req[0]["interface"] / req[0]["transcript"] / f"{req[0]['row']:05d}_{req[-1]['row']:05d}.json"


def get(url: str) -> tuple[int, dict, bytes, float]:
    """One GET: the status (0 for a network error), the headers, the body as received, the seconds."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    start = time.monotonic()
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as r:
            return r.status, dict(r.headers), r.read(), time.monotonic() - start
    except urllib.error.HTTPError as e:
        return e.code, dict(e.headers or {}), e.read() or b"", time.monotonic() - start
    except Exception as e:      # a refused or reset connection, a timeout, a name that does not resolve
        return 0, {}, f"{type(e).__name__}: {e}".encode(), time.monotonic() - start


def entries(body: bytes) -> dict[str, dict]:
    """A reply's entries by the description that was submitted (their "submitted_variant")."""
    reply = json.loads(body)
    if not isinstance(reply, dict):
        raise ValueError("the reply is not a JSON object")
    return {v["submitted_variant"]: v for v in reply.values() if isinstance(v, dict) and "submitted_variant" in v}


def whole(body: bytes, req: list[dict]) -> str:
    """"" if the reply holds one entry for each description submitted, or else what is wrong."""
    try:
        found = entries(body)
    except Exception as e:
        return f"not the JSON expected ({type(e).__name__}: {e})"[:200]
    missing = [it["description"] for it in req if it["description"] not in found]
    return f"no entry for {len(missing)} of {len(req)} descriptions, the first {missing[0]}" if missing else ""


def versions(body: bytes) -> dict:
    meta = json.loads(body).get("metadata") or {}
    return {k: meta.get(k) for k in VERSION_KEYS}


class Journal:
    def __init__(self, raw: Path):
        self.path = raw / "journal.jsonl"

    def write(self, **record) -> None:
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(json.dumps({"utc": utc(), **record}, ensure_ascii=False, sort_keys=True) + "\n")


def hello(raw: Path, journal: Journal, base: str, when: str) -> dict:
    status, _, body, seconds = get(f"{base}/hello/")
    if status != 200:
        journal.write(what="hello", when=when, status=status, outcome="failed")
        sys.exit(f"/hello/ answered {status}: {body[:200]!r}")
    name = f"hello_{utc().replace(':', '').replace('-', '')}.json"
    (raw / name).write_bytes(body)
    journal.write(what="hello", when=when, status=status, file=name, bytes=len(body),
                  sha256=hashlib.sha256(body).hexdigest(), seconds=round(seconds, 3), outcome="stored")
    return versions(body)


def stored(raw: Path, req: list[dict]) -> bool:
    p = reply_path(raw, req)
    return p.is_file() and whole(p.read_bytes(), req) == ""


def first_hello(raw: Path) -> dict | None:
    """The versions /hello/ gave the first time this folder's run asked, if it has."""
    files = sorted(raw.glob("hello_*.json"))
    return versions(files[0].read_bytes()) if files else None


def run(args) -> None:
    raw = Path(args.raw).expanduser().resolve()
    raw.mkdir(parents=True, exist_ok=True)
    items = descriptions(read_list(Path(args.list).expanduser()))
    reqs = requests_of(items)
    todo = [r for r in reqs if not stored(raw, r)]
    journal = Journal(raw)
    journal.write(what="run starts", list=str(args.list), list_sha256=hashlib.sha256(Path(args.list).expanduser().read_bytes()).hexdigest(),
                  descriptions=len(items), requests=len(reqs), already_stored=len(reqs) - len(todo))
    print(f"{utc()} {len(items)} descriptions in {len(reqs)} requests; {len(reqs) - len(todo)} already stored", flush=True)
    if not todo:
        return
    before = first_hello(raw)
    at_start = hello(raw, journal, args.base, "start")
    print(f"{utc()} /hello/: {at_start}", flush=True)
    if before is not None and before != at_start:
        journal.write(what="run stops", why="/hello/ does not answer as it did the first time", first=before, now=at_start)
        sys.exit(5)
    for gene in DEPOSITS_ON_ENSEMBL:
        path = raw / f"transcripts_{gene}_ensembl.json"
        if not path.is_file():
            status, _, body, seconds = get(f"{args.base}/VariantValidator/tools/gene2transcripts_v2/{gene}/all/ensembl/{BUILD}?content-type=application%2Fjson")
            journal.write(what="transcripts held", gene=gene, status=status, bytes=len(body), sha256=hashlib.sha256(body).hexdigest(), seconds=round(seconds, 3))
            if status != 200:
                sys.exit(f"the tool's list of Ensembl transcripts for {gene} answered {status}")
            path.write_bytes(body)
            time.sleep(PAUSE)
    done = [r for r in reqs if stored(raw, r)]
    of_replies = versions(reply_path(raw, done[0]).read_bytes()) if done else None
    in_a_row = rounds = made = 0
    began = time.monotonic()
    for k, req in enumerate(todo, start=1):
        url, wait = url_of(req, args.base), FIRST_WAIT
        while True:
            status, headers, body, seconds = get(url)
            made += 1
            record = dict(what="request", interface=req[0]["interface"], transcript=req[0]["transcript"],
                          rows=[it["row"] for it in req], submitted=[it["description"] for it in req], status=status,
                          bytes=len(body), sha256=hashlib.sha256(body).hexdigest(), seconds=round(seconds, 3))
            problem = whole(body, req) if status == 200 else f"status {status}"
            if not problem:
                seen = versions(body)
                of_replies = of_replies or seen
                if seen != of_replies:
                    journal.write(**record, outcome="the reply's versions differ from the first stored reply's", versions=seen, first=of_replies)
                    (raw / "stopped_reply.json").write_bytes(body)
                    sys.exit(5)
                path = reply_path(raw, req)
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(body)
                journal.write(**record, outcome="stored", file=str(path.relative_to(raw)))
                in_a_row = rounds = 0
                break
            if status not in (0, 200, 429) and not 500 <= status < 600:
                journal.write(**record, outcome="stopped: a status that is not retried", body=body[:500].decode("utf-8", "replace"))
                sys.exit(4)
            in_a_row += 1
            after = headers.get("Retry-After") or headers.get("retry-after")
            pause = float(after) if after and after.replace(".", "", 1).isdigit() else wait
            wait = min(wait * 2, LONGEST_WAIT)
            journal.write(**record, outcome=f"failed ({problem}); again after {pause:g} s", retry_after=after,
                          body=body[:300].decode("utf-8", "replace"))
            if in_a_row >= IN_A_ROW:
                rounds += 1
                if rounds >= ROUNDS:
                    journal.write(what="run stops", why=f"{IN_A_ROW} failures in a row, {ROUNDS} rounds")
                    sys.exit(3)
                journal.write(what="pause", seconds=ROUND_PAUSE, round=rounds)
                time.sleep(ROUND_PAUSE)
                in_a_row = 0
            else:
                time.sleep(pause)
        elapsed = time.monotonic() - began
        ahead = elapsed / k * (len(todo) - k)
        if k % 20 == 0 or k == len(todo):
            print(f"{utc()} {k}/{len(todo)} requests of this run stored; {made - k} retried; about {ahead / 3600:.1f} h to go", flush=True)
        if k == JUDGE_AFTER and (elapsed + ahead) / 3600 > LIMIT_HOURS:
            journal.write(what="run stops", why=f"{JUDGE_AFTER} requests took {elapsed:.0f} s: the whole would take {(elapsed + ahead) / 3600:.1f} h")
            sys.exit(6)
        time.sleep(PAUSE)
    at_end = hello(raw, journal, args.base, "end")
    journal.write(what="run ends", requests_made=made, retried=made - len(todo), versions_as_at_start=at_end == at_start)
    print(f"{utc()} done: {made} requests made for {len(todo)} stored; /hello/ as at the start: {at_end == at_start}", flush=True)
    if at_end != at_start:
        sys.exit(5)


def recheck(args) -> None:
    raw = Path(args.raw).expanduser().resolve()
    rows = read_list(Path(args.list).expanduser())
    items = descriptions(rows)
    by_row: dict[int, list[dict]] = {}
    for it in items:
        by_row.setdefault(it["row"], []).append(it)
    where = {it["description"]: reply_path(raw, req) for req in requests_of(items) for it in req}
    drawn = sorted(random.Random(RECHECK_SEED).sample(range(1, len(rows) + 1), RECHECK_ROWS))
    journal, out = Journal(raw), raw / "recheck"
    out.mkdir(exist_ok=True)
    result = []
    for n in drawn:
        for it in by_row[n]:
            first = entries(where[it["description"]].read_bytes())[it["description"]]
            path = out / f"{it['interface']}_{n:05d}.json"
            if not (path.is_file() and whole(path.read_bytes(), [it]) == ""):
                status, _, body, seconds = get(url_of([it], args.base))
                journal.write(what="recheck", row=n, submitted=it["description"], status=status, bytes=len(body),
                              sha256=hashlib.sha256(body).hexdigest(), seconds=round(seconds, 3))
                if status != 200 or whole(body, [it]):
                    sys.exit(f"row {n}: {status} {whole(body, [it]) if status == 200 else body[:200]!r}")
                path.write_bytes(body)
                time.sleep(PAUSE)
            second = entries(path.read_bytes())[it["description"]]
            result.append({"row": n, "interface": it["interface"], "submitted": it["description"], "same_entry": first == second,
                           "differs_in": sorted(k for k in set(first) | set(second) if first.get(k) != second.get(k))})
    (out / "comparison.json").write_text(json.dumps({"seed": RECHECK_SEED, "rows": drawn, "entries": result}, indent=1))
    same = sum(r["same_entry"] for r in result)
    print(f"{len(drawn)} rows, {len(result)} entries asked again: {same} the same as stored, {len(result) - same} not")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("what", choices=["run", "recheck"])
    ap.add_argument("--list", required=True)
    ap.add_argument("--raw", required=True)
    ap.add_argument("--base", default=BASE)
    args = ap.parse_args()
    {"run": run, "recheck": recheck}[args.what](args)


if __name__ == "__main__":
    main()
