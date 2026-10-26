"""E20 C13 -- the VCEP specifications' splicing rules, applied to this study's data.

For every PP3 and BP4 rule of data/evidence/vcep/vcep_rules.json (committed before any
score was joined to it), within the rule's gene, on the bands 3-10, 11-50 and 3-50 nt:

  * scope: the rule's position scope against each band; every rule here covers the
    intronic positions from 3 nt on (stated per rule below), so every band is in scope;
  * band: PP3, a score at or above the rule's cut point; BP4, at or below it; VHL's PP3 is
    written "> 0.5", and the strict reading is applied, its count checked against "at or
    above" (equal on these data, or the stage stops);
  * score column: a rule that names how SpliceAI is scored gets the column scored that way;
    one that does not gets the published basis (spliceai_walker), as every cut-point result
    of this study, and the row says which (PALB2's specification names the SVI's
    thresholds, Walker et al. 2023, whose basis that column follows). ENIGMA's states its
    setting in its Appendix J, not in the criterion's text (E20b.2): the largest of the
    four deltas, 'unmasked', Max Distance 10 000. Neither SpliceAI column of this study
    has that distance (spliceai_walker 4,999, spliceai 50), so its rows keep the published
    basis, nothing is computed on another column, and the header records the setting;
  * the ratio, its gates and the zero-event rule: E2's (evid_walker_thresholds._band_stats:
    at least 10 damaging and 10 normal in the gene and band, at least 10 variants in the
    scored band, the denominator bounded, a BP4 tier read off the rule-of-three bound where
    the band holds no damaging variant);
  * the interval: the external genes' (evid_external._band_detail: the gene's variants
    resampled within each class, 2,000 draws, class sizes fixed, seeded from the cell's
    identity);
  * the tiers (Tavtigian, prior 0.10): at the point estimate and at the conservative end of
    the interval (PP3 its lower bound, BP4 its upper bound; a BP4 band without a damaging
    variant at its rule-of-three bound), beside the strengths the specification allows.

Conditions outside the panel (VarSeak, RNA evidence, PVS1) are not evaluated; the row says
so. A gene with no labelled variant in a band gets that as its reason, not E2's gates.

BP7, as each specification defines it (E20b.3), counted as C6 counts
(evid_suppl_c6_c7.labelled: the same labels, the signed offset read the same way, the
same score column), per gene; each row gives the variants at the positions (or in the
band), those meeting the condition, the damaging and normal among them, and the damaging
share:
  * ENIGMA (BRCA1, BRCA2): donor +7 and beyond, acceptor -21 and beyond, with BP4
    (SpliceAI at or below 0.1); per side and both. This is C6's condition (a), and the
    rows must equal c6_bp4_bp7.csv's cell for cell, or the stage stops;
  * PALB2: donor +8 and beyond, acceptor -22 and beyond (+7 and -21 not included); (i) as
    the criterion's text, with no prediction condition, (ii) with BP4 (SpliceAI at or
    below 0.1), the combination the instructions to use allow ("May also apply BP4 to
    achieve Likely Benign"); per side and both;
  * VHL: no position stated, so the study's range, on the bands 3-10, 11-50 and 3-50 nt;
    BP4's SpliceAI half (at or below 0.1; VarSeak not evaluated) and PhyloP at or below
    0.2. PhyloP is the analysis set's phylop (UCSC phyloP100way, hg38; the specification
    names no version); a variant without it is not counted, and the row gives their
    number.
The +/-1-2 nt positions are outside this study's PP3 and BP4 bands and are not counted.

Labels are the analysis set's SGE labels. The VHL panel's concern about using the VHL SGE
as functional evidence (circularity) is not evaluated here.

Output, in reports/evidence/suppl_analyses/, with a header (lines beginning '#'):
  c13_vcep_band_lr.csv   part 'band_lr' (PP3, BP4), part 'bp7', part 'no_specification'

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_c13
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import evid_common as K
from . import evid_external as X
from . import evid_suppl_c6_c7 as C6C7
from . import evid_walker_thresholds as W

OUT_DIR = Path("reports/evidence/suppl_analyses")
OUT = OUT_DIR / "c13_vcep_band_lr.csv"
RULES = Path("data/evidence/vcep/vcep_rules.json")
C6 = OUT_DIR / "c6_bp4_bp7.csv"
CONFIG_PATH = Path("config/walker2023.yaml")
BANDS = {"s3_10": ["s3_10"], "s11_50": ["s11_50"], "s3_50": ["s3_10", "s11_50"]}
PUBLISHED = "spliceai_walker"
# the rules' position scopes, as vcep_rules.json writes them, and the bands each covers
SCOPES = {
    "intronic variants outside the donor and acceptor 1,2 sites": list(BANDS),
    "intronic variants outside the native donor and acceptor splice sites (not the +/-1,2 positions)": list(BANDS),
    "non-canonical splice variants (instructions to use), missense variants excluded": list(BANDS),
    "splice variants (no position stated)": list(BANDS),
    "variants assessed for splicing impact (no position stated)": list(BANDS),
}
C6_CONDITION = "(a) SpliceAI at or below the published cut, 0.1"
# ENIGMA's setting, as vcep_rules.json records it from the specification's Appendix J (E20b.2)
ENIGMA_STATED = "stated in the specification's Appendix J"
ENIGMA_SETTING = ["the largest of the four delta scores", "Max Distance 10 000", "'unmasked' scores"]
NO_LABELS = "no labelled variants at this distance in this gene's assay"
# BP7 as each specification defines it (E20b.3): its scope as vcep_rules.json writes it, and
# the positions (signed offsets, C6's reading) or bands it covers
BP7_ENIGMA = "intronic variants at or beyond +7 (donor side) and -21 (acceptor side)"
BP7_PALB2 = ("deep intronic variants beyond +7 (donor side) and -21 (acceptor side), +7 and -21 not "
             "included: +8 and beyond, -22 and beyond")
BP7_VHL = "silent or intronic variants (no position stated)"
PALB2_SIDES = [("donor, +8 and beyond", lambda o: o >= 8),
               ("acceptor, -22 and beyond", lambda o: o <= -22),
               ("both sides", lambda o: (o >= 8) | (o <= -22))]
VHL_BANDS = {"s3_10": (3, 10), "s11_50": (11, 50), "s3_50": (3, 50)}
PHYLOP, PHYLOP_MAX = "phylop", 0.2
PHYLOP_SOURCE = ("the analysis set's phylop: UCSC phyloP100way, hg38 (docs/column-provenance.md); "
                 "the specification names no version")


def score_column(rule: dict) -> tuple[str, str]:
    stated = rule["score_definition_stated"]
    if ENIGMA_STATED in stated:
        missing = [s for s in ENIGMA_SETTING if s not in stated]
        if missing:
            raise SystemExit(f"[C13] {rule['gene']} {rule['criterion']}: the setting recorded from "
                             f"Appendix J lacks {missing}: stopped")
        return PUBLISHED, ("stated in the specification's Appendix J, not in the criterion's text: "
                           "Max Distance 10 000, which neither SpliceAI column of this study has (see "
                           "the header); the published basis, as every cut-point result of this study")
    if stated.startswith("not stated"):
        return PUBLISHED, ("not stated in the specification: the published basis, as every "
                           "cut-point result of this study (E20.4)")
    if "Walker et al. 2023" in stated:
        return PUBLISHED, ("the specification names the SVI's thresholds, Walker et al. 2023, "
                           "whose basis the published-basis column follows")
    raise SystemExit(f"[C13] {rule['gene']} {rule['criterion']}: a score definition this stage "
                     f"does not know ({stated!r}): stopped")


def band_rows(df: pd.DataFrame, rules: list[dict], bands_pp3: dict, bands_bp4: dict) -> list[dict]:
    rows = []
    for rule in rules:
        if rule["criterion"] not in ("PP3", "BP4"):
            continue
        if rule["position_scope"] not in SCOPES:
            raise SystemExit(f"[C13] {rule['gene']} {rule['criterion']}: a position scope this "
                             "stage does not know: stopped")
        if rule["direction"] not in (">=", "<=", ">"):
            raise SystemExit(f"[C13] {rule['gene']} {rule['criterion']}: direction "
                             f"{rule['direction']!r}: stopped")
        col, basis = score_column(rule)
        cut, side = float(rule["cut"]), ("upper" if rule["criterion"] == "PP3" else "lower")
        bands, direction = (bands_pp3, "pathogenic") if side == "upper" else (bands_bp4, "benign")
        for band, members in BANDS.items():
            base = {"part": "band_lr", "gene": rule["gene"], "vcep": rule["vcep"],
                    "cspec_id": rule["source"]["cspec_id"], "version": rule["version"],
                    "criterion": rule["criterion"], "stratum": band,
                    "rule": f"SpliceAI {rule['direction']} {rule['cut']}",
                    "strength_allowed": "; ".join(rule["strength_allowed"]),
                    "score_column": col, "score_column_basis": basis,
                    "not_evaluated": "; ".join(rule["extra_conditions"])}
            if band not in SCOPES[rule["position_scope"]]:
                rows.append(base | {"status": "not in the rule's scope",
                                    "reason": f"the rule's scope: {rule['position_scope']}"})
                continue
            sub = df[(df["gene"] == rule["gene"]) & df["stratum"].isin(members)
                     & df["y_assay"].notna() & df[col].notna()]
            y = sub["y_assay"].to_numpy(dtype=float)
            s = sub[col].to_numpy(dtype=float)
            st = W._band_stats(y, s, cut if side == "upper" else np.nan,
                               cut if side == "lower" else np.nan)
            n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
            in_band = (s > cut) if rule["direction"] == ">" else ((s >= cut) if side == "upper" else (s <= cut))
            if rule["direction"] == ">" and int(in_band.sum()) != int((s >= cut).sum()):
                raise SystemExit(f"[C13] {rule['gene']} {rule['criterion']} {band}: '> {cut}' and "
                                 f"'>= {cut}' select different variants: stopped")
            row = base | {"n_damaging": n_pos, "n_normal": n_neg,
                          "n_in_band": int(in_band.sum()),
                          "n_damaging_in_band": int((in_band & (y == 1)).sum()),
                          "n_normal_in_band": int((in_band & (y == 0)).sum())}
            if not st:
                rows.append(row | {"status": "not evaluable",
                                   "reason": NO_LABELS if n_pos + n_neg == 0 else
                                   "fewer than 10 damaging or 10 normal variants"})
                continue
            lr = st["lr_pp3"] if side == "upper" else st["lr_bp4"]
            if not np.isfinite(lr):
                rows.append(row | {"status": "not evaluable",
                                   "reason": "fewer than 10 variants in the band"})
                continue
            d = X._band_detail(y, s, cut, side, bands, direction,
                               (rule["gene"], rule["source"]["cspec_id"], rule["criterion"], band, col, "c13"))
            if not (d["lr"] == lr and d["n_pos_band"] == row["n_damaging_in_band"]):
                raise SystemExit(f"[C13] {rule['gene']} {rule['criterion']} {band}: E2's and E7's "
                                 "ratios differ: stopped")
            if side == "upper":
                point, bound = lr, d["lr_lo"]
            else:
                point = st["lr_bp4_for_tier"]
                bound = point if row["n_damaging_in_band"] == 0 else d["lr_hi"]
            rows.append(row | {"status": "ok", "reason": "",
                               "lr": lr, "lr_lo": d["lr_lo"], "lr_hi": d["lr_hi"],
                               "lr_for_tier": point,
                               "tier_at_point": K.tier_of(point, bands, direction),
                               "bound_for_tier": bound,
                               "tier_at_bound": K.tier_of(bound, bands, direction)})
    return rows


def _counts(g: pd.DataFrame, meet: pd.Series) -> dict:
    """C6's counts (evid_suppl_c6_c7.c6) for the variants g at the positions and those of
    them meeting the condition."""
    m = g[meet]
    n_d, n_n = int((m["y_assay"] == 1).sum()), int((m["y_assay"] == 0).sum())
    return {"n_at_positions": len(g),
            "n_damaging_at_positions": int((g["y_assay"] == 1).sum()),
            "n_normal_at_positions": int((g["y_assay"] == 0).sum()),
            "n_meeting": len(m), "n_damaging": n_d, "n_normal": n_n,
            "share_damaging": n_d / len(m) if len(m) else np.nan}


def _bp4_cut(rule: dict) -> float:
    if (rule["direction"], rule["cut"]) != ("<=", 0.1):
        raise SystemExit(f"[C13] {rule['gene']} {rule['criterion']}: BP4 is not SpliceAI at or below "
                         "0.1: stopped")
    return float(rule["cut"])


def bp7_rows(df: pd.DataFrame, rules: list[dict]) -> list[dict]:
    d = C6C7.labelled(df)
    by = {(r["gene"], r["criterion"]): r for r in rules}
    score = C6C7.SCORE
    rows = []
    for rule in rules:
        if rule["criterion"] != "BP7":
            continue
        g = d[d["gene"] == rule["gene"]]
        base = {"part": "bp7", "gene": rule["gene"], "vcep": rule["vcep"],
                "cspec_id": rule["source"]["cspec_id"], "version": rule["version"], "criterion": "BP7",
                "strength_allowed": "; ".join(rule["strength_allowed"]),
                "bp7_scope_in_specification": rule["position_scope"],
                "bp7_conditions_in_specification": "; ".join(rule["extra_conditions"]),
                "score_column": score}
        found = []
        if rule["position_scope"] == BP7_ENIGMA:
            cut = _bp4_cut(rule)
            for side, at in C6C7.C6_SIDES:
                x = g[at(g["signed_offset"])]
                found.append(base | {"side": side, "bp7_condition": f"BP4: SpliceAI at or below {cut!r}",
                                     "checked_against": f"c6_bp4_bp7.csv, {C6_CONDITION}"}
                             | _counts(x, x[score] <= cut))
        elif rule["position_scope"] == BP7_PALB2:
            cut = _bp4_cut(by[(rule["gene"], "BP4")])
            for reading in ("i", "ii"):
                for side, at in PALB2_SIDES:
                    x = g[at(g["signed_offset"])]
                    if reading == "i":
                        cond, meet = "(i) the criterion's text: no prediction condition", \
                            pd.Series(True, index=x.index)
                    else:
                        cond, meet = (f"(ii) with BP4, SpliceAI at or below {cut!r} (instructions to "
                                      "use: 'May also apply BP4 to achieve Likely Benign')"), x[score] <= cut
                    found.append(base | {"side": side, "bp7_condition": cond} | _counts(x, meet))
        elif rule["position_scope"] == BP7_VHL:
            if f"PhyloP <= {PHYLOP_MAX}" not in rule["extra_conditions"]:
                raise SystemExit(f"[C13] {rule['gene']} BP7: no 'PhyloP <= {PHYLOP_MAX}' condition: stopped")
            cut = _bp4_cut(rule)
            for band, (lo, hi) in VHL_BANDS.items():
                x = g[g["signed_offset"].abs().between(lo, hi)]
                if set(x.index) != set(g.index[g["stratum"].isin(BANDS[band])]):
                    raise SystemExit(f"[C13] {rule['gene']} BP7 {band}: the offsets and the set's "
                                     "stratum disagree: stopped")
                has = x[PHYLOP].notna()
                x = x[has]
                found.append(base | {"stratum": band,
                                     "bp7_condition": (f"BP4's SpliceAI half, at or below {cut!r}, and "
                                                       f"PhyloP at or below {PHYLOP_MAX!r}"),
                                     "phylop_source": PHYLOP_SOURCE, "n_phylop_missing": int((~has).sum()),
                                     "not_evaluated": "VarSeak (BP4's other half)"}
                             | _counts(x, (x[score] <= cut) & (x[PHYLOP] <= PHYLOP_MAX)))
        else:
            raise SystemExit(f"[C13] {rule['gene']} BP7: a scope this stage does not know: stopped")
        for r in found:
            r |= {"status": "ok", "reason": NO_LABELS if r["n_at_positions"] == 0 else ""}
        rows += found
    check_against_c6(rows)
    return rows


def check_against_c6(rows: list[dict]) -> None:
    """ENIGMA's BP7 is C6's condition (a): its rows must be c6_bp4_bp7.csv's, cell for cell."""
    c6 = K.read_back(C6, comment="#")
    c6 = c6[c6["condition"] == C6_CONDITION].set_index(["gene", "side"])
    mine = [r for r in rows if r.get("checked_against")]
    if c6.empty or not mine:
        raise SystemExit(f"[C13] {C6}: nothing to check ENIGMA's BP7 against: stopped")
    for r in mine:
        c = c6.loc[(r["gene"], r["side"])]
        for k in ("n_at_positions", "n_damaging_at_positions", "n_normal_at_positions", "n_meeting",
                  "n_damaging", "n_normal", "share_damaging"):
            a, b = float(r[k]), float(c[k])
            if not (a == b or (np.isnan(a) and np.isnan(b))):
                raise SystemExit(f"[C13] {r['gene']} BP7 {r['side']} {k}: {a!r} here, {b!r} in c6 (a): "
                                 "stopped (E20b.6)")


def header(spec: dict, rows: list[dict]) -> list[str]:
    lines = ["# E20 C13: the VCEP specifications' splicing rules applied to this study's data "
             "(src/evid_suppl_c13.py), from data/evidence/vcep/vcep_rules.json; BP7 as each "
             "specification defines it (E20b)."]
    seen = set()
    for r in spec["rules"]:
        key = (r["gene"], r["source"]["cspec_id"])
        if key in seen:
            continue
        seen.add(key)
        lines.append(f"# {r['gene']}: {r['vcep']}, ClinGen CSpec {r['source']['cspec_id']} v{r['version']} "
                     f"({r['status']}, approved {str(r['approved_on'])[:10]}; fetched {r['fetched_utc'][:10]}).")
    none = [g for g, v in spec["genes"].items() if v["specification"] is None]
    no_phylop = sum(r["n_phylop_missing"] for r in rows if r.get("stratum") == "s3_50"
                    and "n_phylop_missing" in r)
    lines += [
        f"# No current, approved specification: {', '.join(none)} (see part 'no_specification').",
        "# Score column: spliceai_walker (published basis: SpliceAI 1.3.1, distance 4,999, unmasked, the "
        "largest of the four deltas), for every rule. PALB2's specification names the SVI's thresholds "
        "(Walker et al. 2023); VHL's states no setting. ENIGMA's (BRCA1, BRCA2) states its setting in its "
        "Appendix J, not in the criterion's text (E20b; Appendices_V1.2_2024-11-18.docx, recorded in "
        "data/evidence/vcep/raw/manifest.json): the largest of the four deltas, 'unmasked' scores, Max "
        "Distance 10 000; SpliceAI is cited by its paper (Jaganathan et al. 2019), with no version, online "
        "tool or annotation named. Neither SpliceAI column of this "
        "study has that distance (spliceai_walker 4,999; spliceai 50), so no row is computed on another "
        "column. The rule's cut point is applied to the column as it stands.",
        "# BP7, as each specification defines it, counted as c6_bp4_bp7.csv counts (the SGE labels, the "
        "signed offset read from the c. notation, spliceai_walker). ENIGMA: donor +7 and beyond, acceptor "
        "-21 and beyond, with BP4 (SpliceAI <= 0.1); the same definition as c6's condition (a), and equal "
        "to it cell for cell. PALB2: donor +8 and beyond, acceptor -22 and beyond (+7 and -21 not "
        "included); (i) the criterion's text, no prediction condition; (ii) with BP4 (SpliceAI <= 0.1), "
        "as the instructions to use allow ('May also apply BP4 to achieve Likely Benign'). VHL: no "
        "position stated, so this study's 3-10, 11-50 and 3-50 nt, with BP4's SpliceAI half (<= 0.1) and "
        f"PhyloP <= {PHYLOP_MAX!r}; PhyloP is {PHYLOP_SOURCE}; variants without it are not counted "
        f"({no_phylop} in 3-50 nt). The +/-1-2 nt positions are outside this study's PP3 and BP4 bands "
        "and are not counted.",
        "# Not evaluated: VarSeak (VHL's PP3, BP4 and BP7: only the SpliceAI half is applied; the "
        "specification allows SpliceAI alone where VarSeak cannot take the variant), RNA evidence and "
        "PVS1 (PALB2). A gene with no labelled variant in a band or at the positions: "
        f"'{NO_LABELS}'.",
        "# Ratios, gates and the zero-event rule: E2's (evid_walker_thresholds._band_stats); interval: "
        "evid_external._band_detail (within-gene resampling by class, 2,000 draws, seeded from the cell); "
        "tiers: Tavtigian at a prior of 0.10, at the point estimate (BP4: the rule-of-three bound where the "
        "band holds no damaging variant) and at the conservative end (PP3 the lower bound, BP4 the upper).",
        "# Labels: this study's analysis set, the SGE labels. The VHL panel's concern about using the VHL "
        "SGE as functional evidence (circularity) is not evaluated here. No VCEP rule is taken as the "
        "truth for this study's thresholds, and no threshold is changed.",
    ]
    return lines


def main() -> None:
    spec = json.loads(RULES.read_text())
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    bands_pp3, bands_bp4 = K.acmg_bands(cfg)
    df = K.load_set()
    bp7 = bp7_rows(df, spec["rules"])
    rows = band_rows(df, spec["rules"], bands_pp3, bands_bp4) + bp7
    rows += [{"part": "no_specification", "gene": g, "status": "no specification evaluated",
              "reason": v["note"]} for g, v in spec["genes"].items() if v["specification"] is None]
    out = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        f.write("\n".join(header(spec, bp7)) + "\n")
        out.to_csv(f, index=False)
    print(f"[C13] wrote {OUT} ({len(out)} rows)")
    print(out[[c for c in ("part", "gene", "criterion", "stratum", "side", "n_damaging", "n_normal",
                           "n_in_band", "lr", "lr_lo", "lr_hi", "tier_at_point", "tier_at_bound",
                           "status") if c in out.columns]].to_string(index=False))


if __name__ == "__main__":
    main()
