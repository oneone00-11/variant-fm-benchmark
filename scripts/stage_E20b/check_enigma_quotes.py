"""E20b.2 (commit 2): the sentences of the ENIGMA specification that say how SpliceAI is
scored, checked verbatim against the cached files (which stay out of git: their licence is
unclear), and written to raw/enigma_quote_check.json for the rules file and the tests.

The files are the registry's (fetch_enigma_specification.py; their sha256 must be the
manifest's): for BRCA1 (GN092) and for BRCA2 (GN097), Specifications_V1.2 (Figure 1A) and
Appendix_V1.2 (Appendix J). Each file's text is read paragraph by paragraph from its Word
XML (the w:t runs of each w:p, joined), and every quote must stand in one paragraph of the
section it is cited from, in both specifications' copies: Appendix J runs from its heading
to the next appendix's. Figure 1A is drawn as images; its caption is checked as text, its
images are found through the document's relationships (those between its caption and
Figure 1B's), and the splicing key inside the first is transcribed, marked as read from the
image and not machine-checked. For each part of the setting E20b.2 asks about (the source or online tool,
the version, the annotation or genome build, the distance, the masking, which score), the
Appendix J paragraphs that speak to it are recorded; so is every paragraph of either file
that names SpliceAI, so that no other statement of a setting goes unseen.

Run in the reference image (修改单 §0 rule 9), offline:
    docker run ... IMAGE bash scripts/stage_E20b/in_image.sh scripts/stage_E20b/check_enigma_quotes.py
"""
from __future__ import annotations

import hashlib
import html
import json
import re
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of  # noqa: E402

RAW = REPO / "phase1/data/evidence/vcep/raw"
MANIFEST = RAW / "manifest.json"
OUT = RAW / "enigma_quote_check.json"
APPENDIX_J = "Appendix J: Bioinformatic Analysis and Code Application - incorporating PM1, PP3, BP1, BP4, BP7, and PS1"
NEXT_APPENDIX = re.compile(r"^Appendix K\b")
QUOTES = [
    "For splicing predictions, our recommendations are based on results from unpublished analyses "
    "conducted by the SVI Splicing subgroup.",
    "For variants outside of the donor/acceptor ±1,2 positions, the optimal Splice AI cutpoints (based on "
    "the maximum ∆ score of donor loss, acceptor loss, donor gain or acceptor gain) to predict "
    "spliceogenicity were:",
    "≤0.1 for BP4 for no predicted splicing impact",
    "≥0.2 for PP3 for predicted splicing impact",
    "We recommend setting the Max Distance to 10 000 and display ‘unmasked’ and REF & ALT scores when "
    "curating variants.",
    "Using ‘unmasked’ scores may require additional interpretation to determine predicted splicing impact.",
    "Calibration of predictors against functional and splicing assay results showed that the SpliceAI tool "
    "(Jaganathan et al 2019, PMID:\u00a030661751) outperformed multiple other splicing predictors (including "
    "MaxEntScan).",
]
# what E20b.2 asks of the setting, searched for in Appendix J; and every paragraph of either
# file that names SpliceAI, so that no other statement of a setting goes unseen
ASKED = {"source or online tool": r"(?i)look\s*-?up|broad|https?://|www\.",
         "version": r"(?i)\bversion\b|\bv\s?\d|\b1\.3(\.\d)?\b",
         "annotation or genome build": r"(?i)gencode|ensembl|refseq|\bmane\b|grch3\d|\bhg\d\d\b",
         "distance": r"(?i)distance",
         "masking": r"(?i)mask",
         "which score": r"(?i)\bdelta\b|∆|Δ|\bmax(imum)?\b"}
SPLICEAI = re.compile(r"(?i)splice\s?ai")
FIGURE_1A = "Figure 1A: Application of bioinformatic codes based on variant type and location."
FIGURE_1A_KEY = ("Splice Predictions, using SpliceAI / No impact ≤0.1 / Not informative (do not use) >0.1 "
                 "and <0.2 / Impact ≥0.2")
NEXT_FIGURE = re.compile(r"^Figure 1B\b")


def images_between(z: zipfile.ZipFile, start: int, stop: int) -> list[str]:
    """The media files drawn in paragraphs start..stop-1, through the document's relationships."""
    xml = z.read("word/document.xml").decode("utf-8")
    rels = dict(re.findall(r'Id="([^"]+)"[^>]*Target="([^"]+)"',
                           z.read("word/_rels/document.xml.rels").decode("utf-8")))
    paras = re.findall(r"<w:p[ >].*?</w:p>", xml, flags=re.S)
    return ["word/" + rels[r] for i in range(start, stop) for r in re.findall(r'r:embed="([^"]+)"', paras[i])]


def paragraphs(path: Path) -> tuple[list[str], zipfile.ZipFile]:
    z = zipfile.ZipFile(path)
    xml = z.read("word/document.xml").decode("utf-8")
    out = []
    for p in re.findall(r"<w:p[ >].*?</w:p>", xml, flags=re.S):
        out.append(html.unescape("".join(re.findall(r"<w:t(?: [^>]*)?>([^<]*)</w:t>", p))))
    return out, z


def section(paras: list[str], heading: str) -> tuple[int, int]:
    starts = [i for i, p in enumerate(paras) if p.strip() == heading]
    if not starts:
        raise SystemExit(f"the heading {heading!r} is not in the file: stopped")
    s = starts[-1]                      # the section itself, after the table of contents
    e = next((i for i in range(s + 1, len(paras)) if NEXT_APPENDIX.match(paras[i].strip())), len(paras))
    return s, e


def main() -> None:
    man = json.loads(MANIFEST.read_text())
    specs = man["enigma_specification"]["specifications"]
    record = {"what": "the ENIGMA specification's sentences on how SpliceAI is scored, checked verbatim "
                      "(E20b.2)", "written_by": rel(Path(__file__)), "files": {}, "quotes": [],
              "figure_1a": {}}
    with RunLog(__file__) as log:
        log.input(MANIFEST)
        texts = {}
        for gn, s in specs.items():
            for a in s["attachments"]:
                if a.get("cached_at"):
                    p = REPO / a["cached_at"]
                    log.input(p)
                    if sha256_of(p) != a["sha256"]:
                        raise SystemExit(f"{a['cached_at']} is not the file the manifest records: stopped")
                    texts[(gn, a["name"])] = paragraphs(p)
                    record["files"][f"{gn} {a['name']}"] = {"cached_at": a["cached_at"], "sha256": a["sha256"],
                                                           "filename": a["filename"]}
        for q in QUOTES:
            where = {}
            for gn in specs:
                paras, _ = texts[(gn, "Appendix_V1.2")]
                s, e = section(paras, APPENDIX_J)
                hits = [i for i in range(s, e) if q in paras[i]]
                where[gn] = {"found": bool(hits), "paragraph": hits[0] if hits else None,
                             "section_paragraphs": [s, e]}
            record["quotes"].append({"quote": q, "section": f"Appendix J ({APPENDIX_J})",
                                     "file": "Appendix_V1.2", "found": where})
        record["appendix_j_searched"], record["spliceai_paragraphs"] = {}, {}
        for gn in specs:
            paras, _ = texts[(gn, "Appendix_V1.2")]
            s, e = section(paras, APPENDIX_J)
            record["appendix_j_searched"][gn] = {
                what: [{"paragraph": i, "names_spliceai": bool(SPLICEAI.search(paras[i])), "text": paras[i]}
                       for i in range(s, e) if re.search(rx, paras[i])]
                for what, rx in ASKED.items()}
            for name in ("Specifications_V1.2", "Appendix_V1.2"):
                paras, _ = texts[(gn, name)]
                record["spliceai_paragraphs"][f"{gn} {name}"] = [
                    {"paragraph": i, "in_appendix_j": name == "Appendix_V1.2" and s <= i < e, "text": p}
                    for i, p in enumerate(paras) if SPLICEAI.search(p)]
        for gn in specs:
            paras, z = texts[(gn, "Specifications_V1.2")]
            hits = [i for i, p in enumerate(paras) if p.strip() == FIGURE_1A]
            stop = next((i for i in range(hits[0] + 1, len(paras)) if NEXT_FIGURE.match(paras[i].strip())),
                        len(paras)) if hits else 0
            media = images_between(z, hits[0], stop) if hits else []
            record["figure_1a"][gn] = {
                "caption": FIGURE_1A, "caption_found": bool(hits), "caption_paragraph": hits[0] if hits else None,
                "images": [{"media": m, "sha256": hashlib.sha256(z.read(m)).hexdigest()} for m in media],
                "key_transcribed": FIGURE_1A_KEY, "key_in": media[0] if media else None,
                "note": "the figure is drawn as images, read by eye: the first (synonymous and intronic "
                        "variants) carries the key, transcribed here and not machine-checked, which names "
                        "SpliceAI and its cut points and no setting; the second (missense and in-frame "
                        "variants) carries no splicing key"}
        missing = [q["quote"][:60] for q in record["quotes"] if not all(v["found"] for v in q["found"].values())]
        missing += [gn for gn, f in record["figure_1a"].items() if not (f["caption_found"] and f["images"])]
        record["all_found"] = not missing
        OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n")
        log.output(OUT)
        log.note("missing", missing)
        if missing:
            raise SystemExit(f"not found verbatim: {missing}: stopped")
        print(f"{len(QUOTES)} quotes found in both copies of Appendix J; Figure 1A's caption found")


if __name__ == "__main__":
    main()
