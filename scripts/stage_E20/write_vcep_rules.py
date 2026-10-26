"""E20.3 (C13): the splicing-prediction rules of the VCEP specifications, as a rules file.

Reads the records scripts/stage_E20/fetch_vcep_specs.py saved (phase1/data/evidence/
vcep/raw/) and writes phase1/data/evidence/vcep/vcep_rules.json: one entry per gene and
criterion (PP3, BP4, BP7) of each current, approved specification, and one per training
gene without one. No predictor score is read; the rules file is committed before any
score is joined to it (E20.3).

What is read where. Every quote is copied verbatim from a field of the specification's full
record (the criterion's strength text, its instructions to use, its comments) and must be
found there exactly, or this stops. The strengths a criterion allows are read from the
record (a strength whose applicability is "Applicable"). The rest -- the tool, the direction
and the cut point, the position scope and the extra conditions -- is this script's reading
of the quoted text, written out in READINGS below with the sentence it rests on; a reading
the text does not settle uniquely stops the work for the author (E20.6), and none is
guessed.

Score definition. The registry's text names SpliceAI and its cut points. Where it names
how SpliceAI's score is taken, that is recorded with its quote (PALB2: the SVI's thresholds,
Walker et al. 2023); where it does not, "not stated".

Run in the reference image (修改单 §0 rule 9), offline:
    docker run ... IMAGE bash scripts/stage_E20/in_image.sh scripts/stage_E20/write_vcep_rules.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

VCEP = REPO / "phase1/data/evidence/vcep"
RAW = VCEP / "raw"
OUT = VCEP / "vcep_rules.json"
GENES = ["BAP1", "BARD1", "BRCA1", "BRCA2", "PALB2", "RAD51C", "VHL"]
APPROVED = {"Released"}
NBSP = " "
STATED_NOT = ("not stated: the criterion names SpliceAI and its cut point only, and refers to "
              "the specification's Figure 1A and Appendix J, which the registry's record does "
              "not carry")
ENIGMA_PP3 = (f"Apply PP3 for predicted splicing (SpliceAI{NBSP}≥0.2){NBSP}for silent, "
              f"missense/in-frame{NBSP}(irrespective of location in clinically important "
              "functional domain) and for intronic variants outside of donor and acceptor 1,2 sites.")
ENIGMA_BP4 = ("Intronic variants outside of the native donor and acceptor splice sites (i.e. not "
              f"+/- 1,2 positions) AND no predicted impact{NBSP}via splicing (SpliceAI{NBSP}≤0.1).")
ENIGMA_BP7 = ("Intronic variants located outside conserved donor or acceptor motif positions (at "
              "or beyond positions +7/-21) IF BP4 met.")
ENIGMA_BP7_USE = ("Following convention, this code is applied in addition to BP4 (no splicing "
                  f"prediction, Splice AI{NBSP}≤0.1) to capture the low prior probability of "
                  "pathogenicity of silent variants.")
PALB2_SVI = ("This VCEP recommends SpliceAI thresholds set forth by the SVI in applying PP3 and "
             "BP4 to non-canonical splice variants: Apply PP3 for SpliceAI scores ≥0.2 and apply "
             "BP4 for SpliceAI scores ≤0.1 (Walker et al., 2023).")
VHL_ALONE = "The Splice AI score alone can be applied if VarSeak is unable to accept the variant type."


def enigma(gene: str, spec: str) -> list[dict]:
    return [
        {"gene": gene, "spec": spec, "criterion": "PP3", "strength": "Supporting", "quote": ENIGMA_PP3,
         "tool": "SpliceAI", "direction": ">=", "cut": 0.2,
         "position_scope": "intronic variants outside the donor and acceptor 1,2 sites",
         "extra_conditions": [],
         "further_quotes": [], "score_definition_stated": STATED_NOT},
        {"gene": gene, "spec": spec, "criterion": "BP4", "strength": "Supporting", "quote": ENIGMA_BP4,
         "tool": "SpliceAI", "direction": "<=", "cut": 0.1,
         "position_scope": "intronic variants outside the native donor and acceptor splice sites "
                           "(not the +/-1,2 positions)",
         "extra_conditions": [], "further_quotes": [], "score_definition_stated": STATED_NOT},
        {"gene": gene, "spec": spec, "criterion": "BP7", "strength": "Supporting", "quote": ENIGMA_BP7,
         "tool": "SpliceAI, through BP4", "direction": "<=", "cut": 0.1,
         "position_scope": "intronic variants at or beyond +7 (donor side) and -21 (acceptor side)",
         "extra_conditions": ["BP4 met (SpliceAI <= 0.1)"],
         "further_quotes": [("instructionsToUse", ENIGMA_BP7_USE)],
         "score_definition_stated": STATED_NOT},
    ]


READINGS = enigma("BRCA1", "GN092") + enigma("BRCA2", "GN097") + [
    {"gene": "PALB2", "spec": "GN077", "criterion": "PP3", "strength": "Supporting",
     "quote": (f"Splicing: Predicted impact via splicing (SpliceAI{NBSP}**≥0.2**){NBSP}for silent, "
               "missense/in-frame and for intronic variants outside of donor and acceptor 1,2 sites."),
     "tool": "SpliceAI", "direction": ">=", "cut": 0.2,
     "position_scope": "intronic variants outside the donor and acceptor 1,2 sites",
     "extra_conditions": ["not applied in addition to PVS1 or PVS1_Variable(RNA) (no PVS1 or RNA "
                          "evidence enters this analysis)"],
     "further_quotes": [("instructionsToUse", PALB2_SVI),
                        ("instructionsToUse", "NOTE: PP3 for splice predictions may not be applied in "
                                              "addition to PVS1 or PVS1\\_Variable(RNA) codes.")],
     "score_definition_stated": "the SVI's thresholds, Walker et al. 2023 (instructions to use)"},
    {"gene": "PALB2", "spec": "GN077", "criterion": "BP4", "strength": "Supporting",
     "quote": f"Splicing: No predicted impact via splicing (SpliceAI{NBSP}**≤0.1**). Do not apply for missense variants.",
     "tool": "SpliceAI", "direction": "<=", "cut": 0.1,
     "position_scope": "non-canonical splice variants (instructions to use), missense variants excluded",
     "extra_conditions": ["with RNA data, the instructions to use allow BP4 with BP7_Variable(RNA) and "
                          "the comments do not (the v1.3 release notes say it is allowed); no RNA "
                          "evidence enters this analysis"],
     "further_quotes": [("instructionsToUse", PALB2_SVI),
                        ("instructionsToUse", "NOTE: BP4 for splice predictions may be applied in "
                                              "conjunction with BP7\\_Variable(RNA) (a lack of observed RNA defect)"),
                        ("additionalComments", "NOTE: BP4 for splice predictions may not be applied in "
                                               "conjunction with BP7_Variable(RNA) (a lack of observed RNA defect)")],
     "score_definition_stated": "the SVI's thresholds, Walker et al. 2023 (instructions to use)"},
    {"gene": "PALB2", "spec": "GN077", "criterion": "BP7", "strength": "Supporting",
     "quote": ("BP7: Use for synonymous and deep intronic variants defined as further than (but not "
               "including) +7 and further than (but not including) -21 at donor and acceptor sites, "
               "respectively."),
     "tool": "none stated in the criterion's text", "direction": None, "cut": None,
     "position_scope": "deep intronic variants beyond +7 (donor side) and -21 (acceptor side), "
                       "+7 and -21 not included: +8 and beyond, -22 and beyond",
     "extra_conditions": ["no splicing-prediction condition stated in the criterion's text; the "
                          "instructions to use say BP4 may also be applied to reach Likely Benign"],
     "further_quotes": [("instructionsToUse", "May also apply BP4 to achieve Likely Benign")],
     "score_definition_stated": "not stated"},
    {"gene": "VHL", "spec": "GN078", "criterion": "PP3", "strength": "Supporting",
     "quote": "For splice variants, concordance of Splice AI (>0.5) and VarSeak (class 4 or class 5).",
     "tool": "SpliceAI", "direction": ">", "cut": 0.5,
     "position_scope": "splice variants (no position stated)",
     "extra_conditions": ["VarSeak class 4 or class 5 (not in this study's panel: not evaluated; "
                          "the instructions to use allow SpliceAI alone where VarSeak cannot take "
                          "the variant)"],
     "further_quotes": [("instructionsToUse", VHL_ALONE)],
     "score_definition_stated": "not stated"},
    {"gene": "VHL", "spec": "GN078", "criterion": "BP4", "strength": "Supporting",
     "quote": ("BP4 can be applied to assess lack of splicing impact, with concordance of Splice AI "
               "(≤0.1) and VarSeak (Class 1 or Class 2)."),
     "tool": "SpliceAI", "direction": "<=", "cut": 0.1,
     "position_scope": "variants assessed for splicing impact (no position stated)",
     "extra_conditions": ["VarSeak class 1 or class 2 (not in this study's panel: not evaluated; "
                          "the instructions to use allow SpliceAI alone where VarSeak cannot take "
                          "the variant)"],
     "further_quotes": [("instructionsToUse", VHL_ALONE)],
     "score_definition_stated": "not stated"},
    {"gene": "VHL", "spec": "GN078", "criterion": "BP7", "strength": "Supporting",
     "quote": ("If BP4 is met for lack of splice effect, BP7 can be applied to silent or intronic "
               "variants where the PhyloP score is ≤0.2."),
     "tool": "SpliceAI and VarSeak through BP4; PhyloP", "direction": "<=", "cut": 0.1,
     "position_scope": "silent or intronic variants (no position stated)",
     "extra_conditions": ["BP4 met", "PhyloP <= 0.2"],
     "further_quotes": [], "score_definition_stated": "not stated"},
]


def criterion(record: dict, label: str) -> dict:
    hits = [c["entContent"] for c in record["ld"]["CriteriaCode"] if c["entContent"].get("label") == label]
    if len(hits) != 1:
        raise SystemExit(f"{record['entId']}: {len(hits)} criteria labelled {label}: stopped")
    return hits[0]


def main() -> None:
    with RunLog(__file__) as log:
        man_path = RAW / "manifest.json"
        log.input(man_path)
        man = json.loads(man_path.read_text())
        fetched = {Path(f["file"]).name: f["fetched_utc"] for f in man["files"]}
        records, docs = {}, {}
        for gn in sorted({s["id"] for g in GENES for s in man["found"][g]}):
            for name, store in ((f"{gn}_record.json", records), (f"{gn}.json", docs)):
                log.input(RAW / name)
                body = json.loads((RAW / name).read_text())
                store[gn] = body["data"] if name.endswith("_record.json") else body
        genes = {}
        for g in GENES:
            found = man["found"][g]
            current = [s for s in found if s["status"] in APPROVED]
            if len(current) > 1:
                raise SystemExit(f"{g}: more than one current, approved specification: stopped")
            genes[g] = {"specification": current[0]["id"] if current else None,
                        "registry_entries": found,
                        "note": ("" if current else
                                 "no current, approved specification in the registry"
                                 + ("" if not found else
                                    "; the registry lists " + ", ".join(
                                        f"{s['id']} v{s['version']} ({s['status']})" for s in found)))}
        rules = []
        for r in READINGS:
            if genes[r["gene"]]["specification"] != r["spec"]:
                raise SystemExit(f"{r['gene']}: the reading names {r['spec']}, the registry "
                                 f"{genes[r['gene']]['specification']}: stopped")
            rec = records[r["spec"]]
            c = criterion(rec, r["criterion"])
            strengths = [s for s in c["strengthDescriptor"] if s.get("applicability") == "Applicable"]
            text = {s["strength"]: s.get("text") or "" for s in strengths}
            if r["quote"] not in text.get(r["strength"], ""):
                raise SystemExit(f"{r['gene']} {r['criterion']}: the quote is not in the "
                                 f"{r['strength']} text: stopped")
            further = []
            for field, q in r["further_quotes"]:
                if q not in (c.get(field) or ""):
                    raise SystemExit(f"{r['gene']} {r['criterion']}: a quote is not in {field}: stopped")
                further.append({"quote": q, "section": f"{r['criterion']}, {field}"})
            cut_text = {">=": "≥", "<=": "≤", ">": ">"}.get(r["direction"] or "", None)
            if r["cut"] is not None and r["criterion"] != "BP7":
                if f"{cut_text}{r['cut']}" not in r["quote"].replace("**", "").replace(" ", ""):
                    raise SystemExit(f"{r['gene']} {r['criterion']}: the cut point read is not in "
                                     "the quote: stopped")
            ent = rec["entContent"]
            rules.append({
                "gene": r["gene"], "vcep": docs[r["spec"]]["affiliation"]["label"],
                "source": {"cspec_id": r["spec"], "record": rel(RAW / f"{r['spec']}_record.json"),
                           "url": docs[r["spec"]]["assertionMethod"]["url"],
                           "doi": (ent.get("doi") or {}).get("docDoi") or (ent.get("doi") or {}).get("conceptDoi")},
                "version": ent.get("version"), "approved_on": ent.get("approvedOn"),
                "status": docs[r["spec"]].get("currentStatus"),
                "fetched_utc": fetched[f"{r['spec']}_record.json"],
                "criterion": r["criterion"], "tool": r["tool"],
                "score_definition_stated": r["score_definition_stated"],
                "direction": r["direction"], "cut": r["cut"],
                "strength_allowed": [s["strength"] for s in strengths],
                "strength_quoted": r["strength"],
                "position_scope": r["position_scope"], "extra_conditions": r["extra_conditions"],
                "quote": r["quote"], "section": f"{r['criterion']}, {r['strength']} (strength text)",
                "further_quotes": further,
            })
        out = {"what": "the splicing-prediction rules (PP3, BP4, BP7) of the current, approved VCEP "
                       "specifications of the seven training genes, as ClinGen's Criteria "
                       "Specification Registry serves them (E20.3)",
               "written_by": rel(Path(__file__)), "source_manifest": rel(man_path),
               "licence": "ClinGen's terms of use (raw/clingen_terms_of_use.html): curated content "
                          "under CC0 1.0, attribution and access date requested; the registry's "
                          "Zenodo deposits: cc-pddc",
               "genes": genes, "rules": rules}
        OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
        log.output(OUT)
        log.note("rules", len(rules))
        print(f"{rel(OUT)}: {len(rules)} rules")


if __name__ == "__main__":
    main()
