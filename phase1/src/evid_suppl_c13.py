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
    thresholds, Walker et al. 2023, whose basis that column follows);
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

Conditions outside the panel (VarSeak, PhyloP, RNA evidence, PVS1) are not evaluated; the row
says so. BP7 is read from c6_bp4_bp7.csv, not recomputed: per gene and side, the variants
at Walker 2023's BP7 positions (donor +7 and beyond, acceptor -21 and beyond) at or below
SpliceAI 0.1 and the damaging among them, beside the specification's own BP7 scope.

Labels are the analysis set's SGE labels. The VHL panel's concern about using the VHL SGE
as functional evidence (circularity) is not evaluated here.

Output, in reports/evidence/suppl_analyses/, with a header (lines beginning '#'):
  c13_vcep_band_lr.csv   part 'band_lr' (PP3, BP4), part 'bp7_from_c6', part 'no_specification'

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


def score_column(rule: dict) -> tuple[str, str]:
    stated = rule["score_definition_stated"]
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
                                   "reason": "fewer than 10 damaging or 10 normal variants"})
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


def bp7_rows(rules: list[dict]) -> list[dict]:
    c6 = K.read_back(C6, comment="#")
    c6 = c6[c6["condition"] == C6_CONDITION]
    if c6.empty:
        raise SystemExit(f"[C13] {C6}: no rows for condition {C6_CONDITION!r}: stopped")
    rows = []
    for rule in rules:
        if rule["criterion"] != "BP7":
            continue
        for r in c6[c6["gene"] == rule["gene"]].itertuples():
            rows.append({"part": "bp7_from_c6", "gene": rule["gene"], "vcep": rule["vcep"],
                         "cspec_id": rule["source"]["cspec_id"], "version": rule["version"],
                         "criterion": "BP7", "side": r.side,
                         "strength_allowed": "; ".join(rule["strength_allowed"]),
                         "bp7_scope_in_specification": rule["position_scope"],
                         "bp7_conditions_in_specification": "; ".join(rule["extra_conditions"]),
                         "c6_positions": "Walker 2023: donor +7 and beyond, acceptor -21 and beyond",
                         "c6_condition": C6_CONDITION,
                         "n_at_positions": int(r.n_at_positions), "n_meeting": int(r.n_meeting),
                         "n_damaging": int(r.n_damaging), "n_normal": int(r.n_normal),
                         "share_damaging": r.share_damaging,
                         "status": "read from c6_bp4_bp7.csv, not recomputed"})
    return rows


def header(spec: dict) -> list[str]:
    lines = ["# E20 C13: the VCEP specifications' splicing rules applied to this study's data "
             "(src/evid_suppl_c13.py), from data/evidence/vcep/vcep_rules.json."]
    seen = set()
    for r in spec["rules"]:
        key = (r["gene"], r["source"]["cspec_id"])
        if key in seen:
            continue
        seen.add(key)
        lines.append(f"# {r['gene']}: {r['vcep']}, ClinGen CSpec {r['source']['cspec_id']} v{r['version']} "
                     f"({r['status']}, approved {str(r['approved_on'])[:10]}; fetched {r['fetched_utc'][:10]}).")
    none = [g for g, v in spec["genes"].items() if v["specification"] is None]
    lines += [
        f"# No current, approved specification: {', '.join(none)} (see part 'no_specification').",
        "# Score column: spliceai_walker (published basis: SpliceAI 1.3.1, distance 4,999, unmasked, the "
        "largest of the four deltas), for every rule: none of the specifications states another setting; "
        "PALB2's names the SVI's thresholds (Walker et al. 2023). The rule's cut point is applied to that "
        "column as it stands.",
        "# Not evaluated: VarSeak (VHL's PP3 and BP4: only the SpliceAI half is applied; the specification "
        "allows SpliceAI alone where VarSeak cannot take the variant), PhyloP (VHL's BP7), RNA evidence and "
        "PVS1 (PALB2). BP7 is read from c6_bp4_bp7.csv as it stands, at Walker 2023's BP7 positions; each "
        "row gives the specification's own BP7 scope beside it.",
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
    rows = band_rows(df, spec["rules"], bands_pp3, bands_bp4) + bp7_rows(spec["rules"])
    rows += [{"part": "no_specification", "gene": g, "status": "no specification evaluated",
              "reason": v["note"]} for g, v in spec["genes"].items() if v["specification"] is None]
    out = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        f.write("\n".join(header(spec)) + "\n")
        out.to_csv(f, index=False)
    print(f"[C13] wrote {OUT} ({len(out)} rows)")
    print(out[[c for c in ("part", "gene", "criterion", "stratum", "side", "n_damaging", "n_normal",
                           "n_in_band", "lr", "lr_lo", "lr_hi", "tier_at_point", "tier_at_bound",
                           "status") if c in out.columns]].to_string(index=False))


if __name__ == "__main__":
    main()
