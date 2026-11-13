"""Every listed variant's HGVS description against what VariantValidator returned for it
(E9.4, stage 42): the check the list's descriptions were given with an external tool,
printed beside the list under layout/variants/.

What it reads, and nothing else:
- the list the layout stage wrote, reports/evidence/layout/variants/analysed_variants.csv;
- data/evidence/hgvs_check/excerpt.csv: for each description submitted to VariantValidator's
  public REST interface, what the tool returned, copied from its reply and left as it was
  (scripts/stage_E9/extract_hgvs_check.py; the replies were fetched once, by
  scripts/stage_E9/fetch_hgvs_check.py, which is not a stage);
- data/evidence/hgvs_check/manifest.json: where those replies came from (the interface,
  /hello/ at the start and at the end, when, the list's sha256, the excerpt's);
- config/layout.yaml, for the number the layout gives these files.
No service is asked and nothing is mapped, normalised or corrected here.

What was submitted. Each row of the list as its RefSeq transcript and its c. description,
"NM_...:c....", on GRCh38, with that transcript as the only one to report. The deposits of
DDX3X and SBDS write c. on the Ensembl transcript of the MANE Select pair: SBDS's rows were
submitted a second time as "ENST...:c....", to the tool's interface for Ensembl
transcripts; DDX3X's were not, because the tool does not hold its Ensembl transcript
(the manifest's "not_submitted", with the tool's own list of the transcripts it holds).

What it writes, under reports/evidence/layout/variants/:
- hgvs_check.csv, one row for each listed variant, in the list's order. Five columns of
  the list, as they stand; the description submitted, the transcript description and the
  GRCh38 genomic description returned, and the tool's warnings; for the rows submitted a
  second time, the same four of that query and its two conclusions; and, last, two
  conclusions, Yes or No: the transcript description returned is the one submitted, and
  the genomic description returned is the list's HGVS g., both character for character;
- hgvs_check_summary.csv, by set and gene: the rows, those with Yes in each conclusion,
  those with a warning, the same of the second query, and a last row for all.
Both carry, as comment lines above the header (the layout's way), the tool, the five
versions /hello/ gave, the dates of the queries and what was not submitted and why.
The two files are added to the layout's manifest.csv, which lists every file under
layout/ with its sha256 (the layout stage, which removes and rewrites the folder, lists
them in contents.csv from config/layout.yaml and leaves the writing to this stage).

Before it writes, four things must hold (checks(); the stage stops otherwise, and the tests
run the same on altered copies):
- every row of the list has exactly one entry of the excerpt on its RefSeq transcript,
  the description submitted being the row's transcript and c., character for character,
  and every row of a gene submitted a second time exactly one on its Ensembl transcript;
- the excerpt has no other row;
- the list is the one the replies were fetched for, and the excerpt the one the manifest
  records (their sha256);
- the summary's numbers are those counted again from the rows.
A description the tool changed, a genomic description that is not the list's, or an error
is not a failed check: it is printed as returned, with No.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

import yaml

LIST = Path("reports/evidence/layout/variants/analysed_variants.csv")
EXCERPT = Path("data/evidence/hgvs_check/excerpt.csv")
MANIFEST = Path("data/evidence/hgvs_check/manifest.json")
CONFIG = Path("config/layout.yaml")
OUT = Path("reports/evidence/layout/variants")
LAYOUT_MANIFEST = Path("reports/evidence/layout/manifest.csv")
REPO_PREFIX = "phase1/"                         # as the layout's contents and manifest name paths
CONTENT = "hgvs_check"                          # the item's name in config/layout.yaml
FILES = ["hgvs_check.csv", "hgvs_check_summary.csv"]
REFSEQ, ENSEMBL = "variantvalidator", "variantvalidator_ensembl"
SECOND = ("SBDS",)                              # the genes submitted a second time, on the Ensembl transcript
VERSION_KEYS = ("variantvalidator_version", "variantvalidator_hgvs_version", "vvta_version",
                "vvdb_version", "vvseqrepo_db")
JOIN = " | "                                    # between two warnings of one description

FROM_LIST = ["Set", "Gene", "RefSeq transcript (MANE Select)", "HGVS c.", "HGVS g."]
FIRST = ["Description submitted", "Transcript description returned",
         "Genomic description returned (GRCh38)", "Warnings"]
AGAIN = ["Description submitted on the Ensembl transcript",
         "Transcript description returned for the Ensembl transcript",
         "Genomic description returned for the Ensembl transcript (GRCh38)",
         "Warnings for the Ensembl transcript",
         "Ensembl transcript description same as submitted",
         "Genomic description for the Ensembl transcript same as the list's"]
CONCLUSIONS = ["Transcript description same as submitted", "Genomic description same as the list's"]
COLUMNS = FROM_LIST + FIRST + AGAIN + CONCLUSIONS
SUMMARY = ["Set", "Gene", "Variants", "Transcript description same as submitted",
           "Genomic description same as the list's", "Variants with a warning",
           "Submitted on the Ensembl transcript", "Ensembl transcript description same as submitted",
           "Genomic description for the Ensembl transcript same as the list's",
           "With a warning for the Ensembl transcript"]
ALL = "All"


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def yes(same: bool) -> str:
    return "Yes" if same else "No"


def _warnings(cell: str) -> str:
    texts = json.loads(cell) if cell else []
    if any(JOIN in t for t in texts):
        raise SystemExit(f"[hgvs] a warning holds {JOIN!r}, which separates warnings: {texts}")
    return JOIN.join(texts)


def match(rows: list[dict], excerpt: list[dict]) -> tuple[list[tuple[dict, dict | None, dict | None]], dict[str, list]]:
    """Each row of the list with its entry on the RefSeq transcript and, for a gene
    submitted a second time, its entry on the Ensembl transcript; and what does not pair:
    rows without an entry or with several, entries no row asks for."""
    by_key: dict[tuple[str, str], list[dict]] = {}
    for e in excerpt:
        by_key.setdefault((e["interface"], e["submitted_description"]), []).append(e)
    problems: dict[str, list] = {"row_without_its_entry": [], "row_with_several_entries": [], "entry_no_row_asks_for": []}
    used, out = set(), []
    for n, r in enumerate(rows, start=1):
        pair = []
        for interface, column, wanted in ((REFSEQ, "RefSeq transcript (MANE Select)", True),
                                          (ENSEMBL, "Ensembl transcript (MANE Select)", r["Gene"] in SECOND)):
            key = (interface, f'{r[column]}:{r["HGVS c."]}')
            found = by_key.get(key, []) if wanted else []
            if wanted and not found:
                problems["row_without_its_entry"].append((n, *key))
            if len(found) > 1:
                problems["row_with_several_entries"].append((n, *key, len(found)))
            if wanted:
                used.add(key)
            pair.append(found[0] if found else None)
        out.append((r, pair[0], pair[1]))
    problems["entry_no_row_asks_for"] = [k for k in by_key if k not in used]
    return out, problems


def build(rows: list[dict], excerpt: list[dict]) -> tuple[list[dict], dict[str, list]]:
    """hgvs_check.csv's rows, in the list's order, and what does not pair (match())."""
    pairs, problems = match(rows, excerpt)
    out = []
    for r, first, again in pairs:
        row = {c: r[c] for c in FROM_LIST}
        f = first or {}
        row |= {FIRST[0]: f.get("submitted_description", ""), FIRST[1]: f.get("transcript_description", ""),
                FIRST[2]: f.get("genomic_description_grch38", ""), FIRST[3]: _warnings(f.get("warnings", ""))}
        if again is None:
            row |= {c: "" for c in AGAIN}
        else:
            row |= {AGAIN[0]: again["submitted_description"], AGAIN[1]: again["transcript_description"],
                    AGAIN[2]: again["genomic_description_grch38"], AGAIN[3]: _warnings(again["warnings"]),
                    AGAIN[4]: yes(again["transcript_description"] == again["submitted_description"]),
                    AGAIN[5]: yes(again["genomic_description_grch38"] == r["HGVS g."])}
        row |= {CONCLUSIONS[0]: yes(bool(first) and row[FIRST[1]] == row[FIRST[0]]),
                CONCLUSIONS[1]: yes(bool(first) and row[FIRST[2]] == r["HGVS g."])}
        out.append(row)
    return out, problems


def _count(rows: list[dict]) -> list[int]:
    again = [r for r in rows if r[AGAIN[0]]]
    return [len(rows), sum(r[CONCLUSIONS[0]] == "Yes" for r in rows), sum(r[CONCLUSIONS[1]] == "Yes" for r in rows),
            sum(bool(r[FIRST[3]]) for r in rows), len(again), sum(r[AGAIN[4]] == "Yes" for r in again),
            sum(r[AGAIN[5]] == "Yes" for r in again), sum(bool(r[AGAIN[3]]) for r in again)]


def summary(check: list[dict]) -> list[dict]:
    """By set and gene, in the list's order, then all: the rows, the Yes of each conclusion,
    the rows with a warning, and the same of the second query."""
    groups: dict[tuple[str, str], list[dict]] = {}
    for r in check:
        groups.setdefault((r["Set"], r["Gene"]), []).append(r)
    out = [dict(zip(SUMMARY, [s, g, *_count(rs)])) for (s, g), rs in groups.items()]
    return out + [dict(zip(SUMMARY, [ALL, ALL, *_count(check)]))]


def recount(check: list[dict]) -> list[dict]:
    """The summary's numbers counted again from the rows' descriptions themselves, not
    from their Yes and No: what checks() holds the summary to."""
    keys = list(dict.fromkeys((r["Set"], r["Gene"]) for r in check)) + [(ALL, ALL)]
    out = []
    for s, g in keys:
        rs = [r for r in check if (s, g) == (ALL, ALL) or (r["Set"], r["Gene"]) == (s, g)]
        again = [r for r in rs if r[AGAIN[0]]]
        out.append(dict(zip(SUMMARY, [
            s, g, len(rs), sum(bool(r[FIRST[0]]) and r[FIRST[1]] == r[FIRST[0]] for r in rs),
            sum(bool(r[FIRST[0]]) and r[FIRST[2]] == r["HGVS g."] for r in rs), sum(r[FIRST[3]] != "" for r in rs),
            len(again), sum(r[AGAIN[1]] == r[AGAIN[0]] for r in again),
            sum(r[AGAIN[2]] == r["HGVS g."] for r in again), sum(r[AGAIN[3]] != "" for r in again)])))
    return out


def checks(rows: list[dict], excerpt: list[dict], manifest: dict, list_sha: str, excerpt_sha: str,
           check: list[dict], summed: list[dict], problems: dict[str, list]) -> dict[str, list]:
    """The four checks: what fails each, empty when it passes."""
    out = {"every_row_has_its_one_entry": problems["row_without_its_entry"] + problems["row_with_several_entries"],
           "the_excerpt_has_no_other_row": list(problems["entry_no_row_asks_for"]),
           "the_list_and_the_excerpt_are_the_manifests": [], "the_summary_is_the_rows_counted_again": []}
    want = len(rows) + sum(r["Gene"] in SECOND for r in rows)
    if len(excerpt) != want and not any(out[k] for k in ("every_row_has_its_one_entry", "the_excerpt_has_no_other_row")):
        out["the_excerpt_has_no_other_row"].append(("rows", len(excerpt), want))
    for name, got, recorded in (("list", list_sha, manifest["list"]["sha256"]),
                                ("excerpt", excerpt_sha, manifest["excerpt"]["sha256"])):
        if got != recorded:
            out["the_list_and_the_excerpt_are_the_manifests"].append((name, got, recorded))
    if manifest["list"]["rows"] != len(rows) or manifest["excerpt"]["rows"] != len(excerpt):
        out["the_list_and_the_excerpt_are_the_manifests"].append(("rows", len(rows), len(excerpt)))
    if [r["HGVS c."] for r in check] != [r["HGVS c."] for r in rows] or len(check) != len(rows):
        out["the_summary_is_the_rows_counted_again"].append(("order",))
    again = recount(check)
    out["the_summary_is_the_rows_counted_again"] += [(a, b) for a, b in zip(summed, again) if a != b]
    if len(summed) != len(again):
        out["the_summary_is_the_rows_counted_again"].append(("rows", len(summed), len(again)))
    return out


def comment_lines(manifest: dict) -> list[str]:
    """What the two files say above their header: the tool, its versions, when it was
    asked, what the conclusions compare, and what was not submitted."""
    hello = [json.loads(manifest[k]["text"])["metadata"] for k in ("hello_at_start", "hello_at_end")]
    versions = ["; ".join(f"{k} {h.get(k)}" for k in VERSION_KEYS) for h in hello]
    said = (f"as /hello/ gave them when the queries started and when they ended: {versions[0]}" if versions[0] == versions[1]
            else f"as /hello/ gave them when the queries started: {versions[0]}; when they ended: {versions[1]}")
    run, again = manifest["run"], manifest["second_asking"]
    days = sorted({run["first_request_utc"][:10], run["last_request_utc"][:10]})
    not_asked = manifest["not_submitted"]
    return [
        "A check of every listed variant's HGVS description with VariantValidator (src/evid_hgvs_check.py, stage 42): "
        "one row for each row of analysed_variants.csv, in its order.",
        f"Tool: {manifest['tool']}, its public REST interface at {manifest['interface']['address']}. Versions, {said}.",
        f"Queried on {' and '.join(days)} (UTC), from {run['first_request_utc']} to {run['last_request_utc']}: "
        f"{manifest['descriptions']['submitted']:,} descriptions in {manifest['descriptions']['requests']:,} requests. "
        f"Afterwards {again['rows']} rows drawn with a fixed seed were asked again, {again['entries']} descriptions: "
        f"{again['entries_the_same_as_stored']} came back the same as stored.",
        "Each row was submitted as its RefSeq transcript (MANE Select) and its HGVS c., on GRCh38, with that transcript "
        f"as the only one to report. \"{CONCLUSIONS[0]}\" is Yes where the transcript description returned is the one "
        f"submitted, character for character; \"{CONCLUSIONS[1]}\" is Yes where the GRCh38 genomic description returned "
        "is the list's HGVS g., character for character. Nothing returned was changed.",
        "The deposits of DDX3X and SBDS write c. on the Ensembl transcript of the MANE Select pair. The rows of "
        f"{' and '.join(SECOND)} were submitted a second time on that transcript, to the tool's interface for Ensembl "
        f"transcripts (the columns that name the Ensembl transcript; empty for the other genes). The rows of "
        f"{not_asked['gene']} were not: the tool does not hold {not_asked['transcript']} "
        f"(data/evidence/hgvs_check/transcripts_{not_asked['gene']}_ensembl.json is its own list of the Ensembl "
        f"transcripts it holds for the gene), so {not_asked['gene']} is checked on the RefSeq transcript of the same pair.",
        f"Warnings are the tool's own, as returned; several for one description are separated by \"{JOIN.strip()}\".",
    ]


def as_text(header: list[str], rows: list[dict], comments: list[str]) -> str:
    out = io.StringIO()
    for line in comments:
        out.write(f"# {line}\n")
    w = csv.DictWriter(out, fieldnames=header, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return out.getvalue()


def grouped(n: int) -> str:
    return f"{n:,}"


def _item() -> dict:
    """This stage's item in config/layout.yaml: its number, and the files it must name."""
    items = [it for it in yaml.safe_load(CONFIG.read_text()).get("variants", []) if it["content"] == CONTENT]
    if len(items) != 1 or [Path(f).name for f in items[0].get("files", [])] != FILES:
        raise SystemExit(f"[hgvs] {CONFIG} does not list {CONTENT} once with the files {FILES}")
    return items[0]


def add_to_manifest(paths: list[Path], number: str, sources: list[Path], manifest_csv: Path = LAYOUT_MANIFEST) -> None:
    """The layout's manifest.csv with these files' rows at its end, every other line as
    it stood (a rerun replaces the rows it added)."""
    text = manifest_csv.read_text()
    header = next(csv.reader(io.StringIO(text)))
    if header != ["file", "bytes", "sha256", "number", "content", "sources", "sources_sha256"]:
        raise SystemExit(f"[hgvs] {manifest_csv}: its columns are {header}")
    names = {REPO_PREFIX + p.as_posix() for p in paths}
    kept = [ln for ln in text.splitlines(keepends=True) if next(csv.reader([ln]))[0] not in names]
    out = io.StringIO()
    w = csv.writer(out, lineterminator="\n")
    for p in paths:
        w.writerow([REPO_PREFIX + p.as_posix(), p.stat().st_size, sha256(p), number, CONTENT,
                    "; ".join(REPO_PREFIX + s.as_posix() for s in sources), "; ".join(sha256(s) for s in sources)])
    manifest_csv.write_text("".join(kept) + out.getvalue())


def prepare() -> tuple[list[dict], list[dict], list[str]]:
    """The check's rows, the summary and the comment lines, after the four checks; the
    stage stops if one fails."""
    rows, excerpt, manifest = read_csv(LIST), read_csv(EXCERPT), json.loads(MANIFEST.read_text())
    check, problems = build(rows, excerpt)
    summed = summary(check)
    failed = {k: v[:5] for k, v in checks(rows, excerpt, manifest, sha256(LIST), sha256(EXCERPT),
                                          check, summed, problems).items() if v}
    if failed:
        raise SystemExit(f"[hgvs] the check's inputs do not hold together: {failed}")
    return check, summed, comment_lines(manifest)


def write(folder: Path = OUT, prepared: tuple | None = None) -> list[Path]:
    """The two files, from what prepare() gives."""
    check, summed, comments = prepared or prepare()
    folder.mkdir(parents=True, exist_ok=True)
    a, b = folder / FILES[0], folder / FILES[1]
    a.write_text(as_text(COLUMNS, check, comments))
    printed = [{k: (grouped(v) if isinstance(v, int) else v) for k, v in r.items()} for r in summed]
    b.write_text(as_text(SUMMARY, printed, comments))
    return [a, b]


def main() -> None:
    item = _item()
    prepared = prepare()
    paths = write(OUT, prepared)
    add_to_manifest(paths, item["number"], [LIST, EXCERPT, MANIFEST])
    total = prepared[1][-1]
    print(f"[hgvs] {total['Variants']:,} listed variants: transcript description same as submitted in "
          f"{total[SUMMARY[3]]:,}, genomic description same as the list's in {total[SUMMARY[4]]:,}, "
          f"{total[SUMMARY[5]]:,} with a warning; wrote {', '.join(p.as_posix() for p in paths)}")


if __name__ == "__main__":
    main()
