"""E9.3: docs/manuscript_bindings_v7.tsv, the binding table of the text as submitted.

v6 (scripts/stage_E8/build_bindings_v6.py, which stays as it is) bound the working draft
of 2026-11-10. The three figure legends were shortened afterwards and the manuscript was
assembled into the file for submission, so the numbers are bound again, this time in the
submitted files themselves, each given as an argument so that a later text is bound by
running this again:

- --manuscript, the manuscript for submission (.docx): its title page, abstract and
  keywords, main text, statements (Data availability to Competing interests) and figure
  legends; its reference list is not read;
- --tables, the four main tables' files (.docx), in order: each one's title, part titles
  and note;
- --supplement, the supplementary information (.pdf): the supplementary notes, read from
  the PDF's own tagged structure as poppler prints it (pdfinfo -struct-text): the heading
  "Supplementary Notes", each note's heading and each paragraph, a list item without its
  label, up to "References cited in the Supplementary Notes". Where poppler is not
  installed (the reference image), --structure gives that print, made where it is. Either
  way the paragraphs are then held to the PDF itself, read with pypdf: each stands in the
  pages' text in order, white space aside, and nothing lies between two of them but a
  list item's label.

The texts are written to a folder outside the repository (--texts), each paragraph on one
line, in the form the number gate reads (scripts/check_evidence_manuscript_numbers.py):
- v7_main.md: "# Abstract", the abstract and the keywords; "# Title page" and its lines
  (the gate reads from the heading "Abstract" on, so the title page's lines are written
  after the abstract); the main text and the statements under their headings, in the
  manuscript's order; "# Figure legends"; "# Tables"; then "# References", where the gate
  stops;
- v7_supplement.md: the notes, ending at "# References cited in the Supplementary Notes";
- the same two with a source tag after each sentence that holds a measured number
  (v7_main_keyed.md, v7_supplement_keyed.md), for the gate's bound mode and self-test.

The table. One row per number the gate reads, in reading order. A number that v6 bound or
declared keeps v6's entry (its cell, or its reason), found by the sentence it stands in;
the numbers read here for the first time, on the title page and in the statements, have
their entries below (NEW_REASONS), and the title page's counts are counted from the
manuscript before they are written down as what they are.

Checks before writing; any failure stops it:
- every number of the two texts has one entry, and every entry of v6's and of the new
  ones is used once;
- every excerpt is a piece of its sentence as the gate reads it;
- the title page's counts are the manuscript's own: words of the main text and of the
  abstract, figure legends, table files, references.

    .venv/bin/python scripts/stage_E9/build_bindings_v7.py \\
        --manuscript ~/Desktop/submission/Manuscript.docx \\
        --tables ~/Desktop/submission/Table_1.docx ~/Desktop/submission/Table_2.docx \\
                 ~/Desktop/submission/Table_3.docx ~/Desktop/submission/Table_4.docx \\
        --supplement ~/Desktop/submission/Supplementary_Information.pdf \\
        --texts ~/Desktop/整理文件/E9_工作文件/v7_texts

In the reference image, which holds no poppler and sees the repository only: copies of the
files in results/scratch/ (ignored by git, removed afterwards), and
--structure <the output of pdfinfo -struct-text on the same PDF>.
"""
from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _runlog import REPO, RunLog  # noqa: E402

sys.path.insert(0, str(REPO / "scripts" / "stage_E8"))
import build_bindings_v6 as V6  # noqa: E402  its entries, its tagging of the keyed texts

sys.path.insert(0, str(REPO / "scripts"))
import check_evidence_manuscript_numbers as gate  # noqa: E402

OUT = REPO / "docs" / "manuscript_bindings_v7.tsv"
SECTIONS = ["Abstract", "Introduction", "Materials and methods", "Results", "Discussion",
            "Data availability", "Code availability", "References", "Acknowledgements",
            "Author contributions", "Funding", "Ethical approval", "Competing interests",
            "Figure legends"]
NOTES, NOTES_END = "Supplementary Notes", "References cited in the Supplementary Notes"
COUNTED = "counted from the manuscript by scripts/stage_E9/build_bindings_v7.py: "

# The numbers read for the first time in v7: (document, a piece of the sentence, the number
# as printed, its place among the numbers so printed in the sentence, what it is). The five
# counts of the title page are filled in by counted_reasons(), from the manuscript.
NEW_REASONS: list[tuple[str, str, str, int, str]] = [
    ("main", "registry entries GN092 and GN097, version 1.2", "1.2", 1,
     "the version of the two ENIGMA specifications, as their registry entries give it (references 39 and 40; "
     "phase1/data/evidence/suppl_analyses/c13_vcep_sources/manifest.json)"),
    ("main", "ghcr.io/oneone00-11/variant-fm-benchmark-reference", "11", 1,
     "digits of the repository owner's name in the reference image's address (container/ENVIRONMENT.json)"),
]


# ------------------------------------------------------------------ the manuscript
def manuscript_parts(path: Path) -> dict:
    """The submitted manuscript's paragraphs by part: the title page's lines, then each
    Heading 1 section's (style, text) pairs."""
    import docx
    doc = docx.Document(str(path))
    paras = [(p.style.name, p.text.strip()) for p in doc.paragraphs if p.text.strip()]
    heads = [t for s, t in paras if s == "Heading 1"]
    if heads != SECTIONS:
        raise SystemExit(f"{path.name}: its parts are {heads}, not {SECTIONS}")
    parts: dict = {"Title page": []}
    here = "Title page"
    for style, text in paras:
        if style == "Heading 1":
            here = text
            parts[here] = []
        else:
            parts[here].append((style, text))
    return parts


def table_paragraphs(path: Path) -> list[str]:
    import docx
    return [p.text.strip() for p in docx.Document(str(path)).paragraphs if p.text.strip()]


def main_text(manuscript: Path, tables: list[Path]) -> tuple[str, dict]:
    parts = manuscript_parts(manuscript)
    line = lambda style, text: ("## " + text) if style == "Heading 2" else text        # noqa: E731
    out = ["# Abstract", ""]
    for _, t in parts["Abstract"]:
        out += [t, ""]
    out += ["# Title page", ""]
    for _, t in parts["Title page"]:
        out += [t, ""]
    for name in SECTIONS[1:]:
        if name == "References":
            continue                                    # the reference list is not read
        out += ["# " + name, ""]
        for style, t in parts[name]:
            out += [line(style, t), ""]
    out += ["# Tables", ""]
    for path in tables:
        for t in table_paragraphs(path):
            out += [t, ""]
    out += ["# References", ""]
    body = [t for name in SECTIONS[1:5] for s, t in parts[name] if s != "Heading 2"]
    counts = {"main text words": sum(len(t.split()) for t in body),
              "abstract words": len(parts["Abstract"][0][1].split()),
              "figures": sum(1 for _, t in parts["Figure legends"] if re.match(r"^Figure \d\. ", t)),
              "tables": len(tables),
              "references": len(parts["References"]),
              "numbers in the reference list": sum(len(V7_NUMBER.findall(t)) for _, t in parts["References"])}
    return "\n".join(out) + "\n", counts


V7_NUMBER = re.compile(r"\d+(?:[.,]\d+)*")


def counted_reasons(text: str, counts: dict) -> list[tuple[str, str, str, int, str]]:
    """The title page's line of counts, each number held to what this script counts in the
    manuscript itself; a count that differs stops the run."""
    m = re.search(r"^Word count: Main text ([\d,]+) words; abstract (\d+) words; (\d+) figures; (\d+) tables; "
                  r"(\d+) references$", text, re.M)
    if not m:
        raise SystemExit("the title page's line of counts is not in the form expected")
    piece = "Word count: Main text"
    names = ["main text words", "abstract words", "figures", "tables", "references"]
    what = ["words of the main text (Introduction to Discussion, headings aside)", "words of the abstract",
            "figure legends in the manuscript", "table files given", "entries of the reference list"]
    out = []
    for printed, name, w in zip(m.groups(), names, what):
        if int(printed.replace(",", "")) != counts[name]:
            raise SystemExit(f"the title page prints {printed} for the {name}, the manuscript has {counts[name]}")
        out.append(("main", piece, printed, 1, COUNTED + w))
    return out


# ------------------------------------------------------------------ the supplementary PDF
def structure_of(pdf: Path, structure: Path | None = None) -> str:
    """The PDF's tagged structure as poppler prints it: the file given, or pdfinfo run here."""
    if structure is not None:
        text = structure.read_text()
    else:
        r = subprocess.run(["pdfinfo", "-struct-text", str(pdf)], capture_output=True, text=True)
        if r.returncode:
            raise SystemExit(f"{pdf.name}: pdfinfo -struct-text failed ({r.stderr.strip()[:200]})")
        text = r.stdout
    if not text.lstrip().startswith("Document"):
        raise SystemExit(f"{pdf.name}: no tagged structure to read")
    return text


def struct_blocks(structure: str) -> list[tuple[str, str]]:
    """The text-holding blocks of a PDF's tagged structure, in order, as poppler prints
    it (pdfinfo -struct-text): (tag, text) for each heading (H1 to H6), paragraph (P) and
    list label (Lbl). A block's text is its lines' strings joined as they stand: a line
    ends with the space that follows its last word, and with none after a hyphen."""
    out: list[list[str]] = []
    for raw in structure.split("\n"):
        s = raw.strip()
        m = re.match(r"^(P|H[1-6]?|Lbl)\b(?: \((?:block|inline)\))?:?$", s)
        if m:
            out.append([m.group(1), ""])
        elif s.startswith('"') and s.endswith('"') and len(s) >= 2 and out:
            out[-1][1] += s[1:-1]
    return [(tag, text) for tag, text in out]


def squeeze(text: str) -> str:
    return re.sub(r"\s+", "", text)


def confirm_in_pdf(pdf: Path, text: str) -> dict:
    """Hold the notes' text to the PDF itself, read with pypdf (each page's text without
    the page number printed above it): every heading and paragraph stands in the pages'
    text in order, white space aside, and between two of them lies nothing, or the label
    of a list item. Stops otherwise; returns what it counted."""
    from pypdf import PdfReader
    pages = []
    for i, page in enumerate(PdfReader(str(pdf)).pages, 1):
        t = page.extract_text()
        first, _, rest = t.partition("\n")
        pages.append(rest if first.strip() == str(i) else t)
    whole = squeeze("".join(pages))
    paras = [re.sub(r"^#+ ", "", ln) for ln in text.split("\n") if ln.strip()]
    at, between = 0, []
    for k, p in enumerate(paras):
        found = whole.find(squeeze(p), at)
        if found < 0:
            raise SystemExit(f"{pdf.name}: its pages do not hold, after the paragraph before it: {p[:100]}")
        if k:
            between.append(whole[at:found])
        at = found + len(squeeze(p))
    other = [g for g in between if g and not re.fullmatch(r"\d+\.", g)]
    if other:
        raise SystemExit(f"{pdf.name}: text between two paragraphs that the structure did not give: {other[:3]}")
    return {"paragraphs and headings found in the pages' text in order": len(paras),
            "list labels between them": sum(1 for g in between if g)}


def supplement_text(pdf: Path, structure: Path | None = None) -> tuple[str, dict]:
    blocks = struct_blocks(structure_of(pdf, structure))
    starts = [i for i, (tag, t) in enumerate(blocks) if tag.startswith("H") and t.strip() == NOTES]
    ends = [i for i, (tag, t) in enumerate(blocks) if t.strip() == NOTES_END]
    if len(starts) != 1 or not ends or ends[-1] < starts[0]:
        raise SystemExit(f"{pdf.name}: the notes' heading or the heading of their references is not found once")
    start, end = starts[0], ends[-1]
    out, labels, headings = [], 0, []
    for tag, t in blocks[start:end]:
        t = " ".join(t.split())
        if tag == "Lbl":
            if not re.fullmatch(r"\d+\.", t):
                raise SystemExit(f"a list label that is not a number: {t!r}")
            labels += 1
        elif tag.startswith("H"):
            level = "##" if t == NOTES else "###"
            headings.append(t)
            out += [f"{level} {t}", ""]
        elif t:
            out += [t, ""]
    out += ["# " + NOTES_END, ""]
    refs = []
    for tag, t in blocks[end + 1:]:
        if tag.startswith("H") or not re.match(r"^\d+\. ", t.strip()):
            break
        refs.append(t)
    text = "\n".join(out) + "\n"
    seen = confirm_in_pdf(pdf, text)
    if seen["list labels between them"] != labels:
        raise SystemExit(f"{pdf.name}: {labels} list labels in the structure, {seen['list labels between them']} in the pages")
    return text, {"note headings": headings[1:], "list labels left out": labels,
                  "numbers in the notes' reference list": sum(len(V7_NUMBER.findall(t)) for t in refs),
                  "entries of the notes' reference list": len(refs), **seen}


# ------------------------------------------------------------------ the rows
def rows_for(doc: str, text: str, bindings: list, reasons: list) -> tuple[list[dict], list]:
    """Each number of a text with its entry, as v6 found it (build_bindings_v6.rows_for,
    with the entries given); the entries matched, for the check that each is used once."""
    out, used = [], []
    ranks: dict[tuple[str, str], int] = {}
    for t in gate.read_tokens(text):
        k = (t.sentence, t.printed)
        ranks[k] = ranks.get(k, 0) + 1
        bound, why = V6._key(bindings, doc, t, ranks[k]), V6._key(reasons, doc, t, ranks[k])
        if len(bound) + len(why) > 1:
            raise SystemExit(f"{doc} line {t.line}: {t.printed} has {len(bound) + len(why)} entries")
        clean = t.excerpt_space
        excerpt = clean[max(0, t.context.start - 30):t.context.end + 30].strip()
        if excerpt not in clean:
            raise SystemExit(f"{doc} line {t.line}: the excerpt is not a piece of the sentence")
        row = {"document": doc, "line": t.line, "occurrence": t.occurrence, "printed": t.printed, "excerpt": excerpt}
        blank = {"file": "", "rows": "", "column": "", "derive": "", "scale": ""}
        if bound:
            used.append(bound[0])
            row |= {"kind": "measured", "reason": "", **bound[0][4]}
        elif why:
            used.append(why[0])
            row |= {"kind": "constant", "reason": why[0][4], **blank}
        elif t.reason:
            row |= {"kind": "constant", "reason": t.reason, **blank}
        else:
            raise SystemExit(f"{doc} line {t.line}: {t.printed} is measured and has no binding ({t.sentence[:120]})")
        out.append(row)
    return out, used


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manuscript", required=True, type=Path)
    ap.add_argument("--tables", required=True, type=Path, nargs=4)
    ap.add_argument("--supplement", required=True, type=Path)
    ap.add_argument("--structure", type=Path, help="the output of pdfinfo -struct-text on --supplement, "
                                                   "where poppler is not installed")
    ap.add_argument("--texts", required=True, type=Path)
    args = ap.parse_args()
    ms, pdf = args.manuscript.expanduser().resolve(), args.supplement.expanduser().resolve()
    structure = args.structure.expanduser().resolve() if args.structure else None
    tables = [p.expanduser().resolve() for p in args.tables]
    folder = args.texts.expanduser().resolve()
    with RunLog(__file__) as log:
        for name, p in [("manuscript", ms), ("supplement", pdf)] + [(f"table_{i}", p) for i, p in enumerate(tables, 1)]:
            log.note(f"{name}_sha256", hashlib.sha256(p.read_bytes()).hexdigest())
        folder.mkdir(parents=True, exist_ok=True)
        main_md, counts = main_text(ms, tables)
        supp_md, notes = supplement_text(pdf, structure)
        log.note("structure_read_from", "the file given, sha256 " + hashlib.sha256(structure.read_bytes()).hexdigest()
                 if structure else "pdfinfo -struct-text, run here")
        (folder / "v7_main.md").write_text(main_md)
        (folder / "v7_supplement.md").write_text(supp_md)
        log.note("counted_in_the_manuscript", counts)
        log.note("read_from_the_supplementary_pdf", notes)
        reasons = V6.REASONS + NEW_REASONS + counted_reasons(main_md, counts)
        rows, used = [], []
        for doc, text in (("main", main_md), ("supplement", supp_md)):
            r, u = rows_for(doc, text, V6.BINDINGS, reasons)
            rows += r
            used += u
            (folder / f"v7_{doc}_keyed.md").write_text(V6.keyed(text, r))
            log.note(f"{doc}_sha256", hashlib.sha256((folder / f"v7_{doc}.md").read_bytes()).hexdigest())
        unused = [e[:4] for e in V6.BINDINGS + reasons if e not in used]
        twice = [e[:4] for e in used if used.count(e) > 1]
        if unused or twice:
            raise SystemExit(f"entries matched no number: {unused[:5]}; matched twice: {twice[:5]}")
        gate.write_bindings(OUT, rows)
        n = {k: sum(1 for r in rows if r["kind"] == k) for k in ("measured", "constant")}
        log.note("rows", n | {"total": len(rows), "main": sum(1 for r in rows if r["document"] == "main"),
                              "supplement": sum(1 for r in rows if r["document"] == "supplement")})
        log.output(OUT)
        print(f"{OUT.relative_to(REPO)}: {len(rows)} numbers, {n['measured']} measured, {n['constant']} constants; "
              f"texts in {folder}")


if __name__ == "__main__":
    main()
