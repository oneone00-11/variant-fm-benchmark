"""E19 C6 -- BP4 at the BP7 positions; C7 -- SpliceAI's published cut points by position.

Positions are the signed intronic offset (evid_suppl_c2.signed_offset: the companion
atlas's offset pattern on the c. notation, the sign kept). Labels are the analysis
set's (y_assay: damaging 1, normal 0); SpliceAI is the published basis
(spliceai_walker); the cut points are Walker 2023's, from config/walker2023.yaml (PP3
at or above 0.2, BP4 at or below 0.1).

C6 (E19.8). The BP7 positions, after Walker 2023: the donor side at +7 and beyond, the
acceptor side at -21 and beyond (to the set's 50 nt). At these positions, the SGE
damaging and normal variants whose score meets each of two conditions, and the
damaging share among them:
  (a) the published cut: SpliceAI at or below 0.1;
  (b) Table 3(b): SpliceAI at or below the BP4 Supporting threshold of the variant's
      band, +7 ... +10 in the 3-10 nt band and every other BP7 position in the 11-50
      nt band. The threshold is the fitted one, evidence_thresholds.csv's
      bp4_threshold_insample at full precision, read as E18's external test reads
      Table 3(b)'s thresholds (the author's decision of 2026-10-26). Table 3(b) prints
      it to three figures, rounded in the conservative direction (PP3 up, BP4 down).
      The file's header names the variants the printed values would leave out, as the
      stage finds them, and says the damaging counts are the same; were a damaging
      variant among them, or one between the two values rather than on the fitted one,
      the stage stops.
Per gene and for the seven genes, for each side and both together; counts and shares
only.

C7 (E19.8). Walker 2023's position groups: on the donor side +3...+6, +7...+10 and
+11...+50; on the acceptor side -3...-10, -11...-20 and -21...-50. Per group, pooled over
the seven genes: the damaging and normal variants, and the band likelihood ratio at
each published cut point (at or above 0.2, PP3; at or below 0.1, BP4) with its
interval, by E2's own functions (evid_walker_thresholds._band_stats and the
gene-clustered _boot_ci, 2,000 draws, seeded from the cell's identity). A group with
fewer than 10 damaging or 10 normal variants gets its counts only, and a band with
fewer than 10 occupants no ratio (E2's rules).

Outputs, in reports/evidence/suppl_analyses/:
  c6_bp4_bp7.csv       per condition x side x gene; its header (lines beginning '#')
                       gives the thresholds of (b) and the printed values' difference
  c7_by_position.csv   per position group

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_c6_c7
"""
from __future__ import annotations

from decimal import Decimal
from pathlib import Path

import numpy as np
import pandas as pd

from . import evid_common as K
from . import evid_suppl_c2 as C2
from . import evid_tables as T
from . import evid_walker_thresholds as W

OUT_DIR = Path("reports/evidence/suppl_analyses")
OUT_C6 = OUT_DIR / "c6_bp4_bp7.csv"
OUT_C7 = OUT_DIR / "c7_by_position.csv"
FROZEN = Path("reports/evidence/evidence_thresholds.csv")
TABLE3B = T.OUT / "table3b_thresholds.csv"
SCORE = "spliceai_walker"
ALL = C2.ALL
# Walker 2023's BP7 positions
BP7_DONOR_FROM, BP7_ACCEPTOR_FROM = 7, -21
C6_SIDES = [("donor, +7 and beyond", lambda o: o >= BP7_DONOR_FROM),
            ("acceptor, -21 and beyond", lambda o: o <= BP7_ACCEPTOR_FROM),
            ("both sides", lambda o: (o >= BP7_DONOR_FROM) | (o <= BP7_ACCEPTOR_FROM))]
C7_GROUPS = [("donor", "+3...+6", 3, 6), ("donor", "+7...+10", 7, 10),
             ("donor", "+11...+50", 11, 50), ("acceptor", "-3...-10", -10, -3),
             ("acceptor", "-11...-20", -20, -11), ("acceptor", "-21...-50", -50, -21)]


def labelled(df: pd.DataFrame) -> pd.DataFrame:
    d = df.dropna(subset=["y_assay", SCORE]).copy()
    d["signed_offset"] = d["hgvs_c"].map(C2.signed_offset)
    if (d["signed_offset"].abs() != d["intron_offset_abs"]).any():
        raise SystemExit("[C6/C7] a signed offset disagrees with the set's intron_offset_abs")
    return d


def bp4_supporting() -> dict[str, float]:
    """The fitted BP4 Supporting thresholds of the 3-10 and 11-50 nt bands, and
    Table 3(b)'s printed values of them."""
    e = K.read_back(FROZEN)
    t3b = K.read_back(TABLE3B, dtype=str).set_index(["Predictor", "Band"])
    out = {}
    for st in ("s3_10", "s11_50"):
        r = e[(e.tool == SCORE) & (e.stratum == st) & (e.tier == "supporting")].iloc[0]
        if not bool(r.bp4_threshold_reachable):
            raise SystemExit(f"[C6] {SCORE} {st}: no BP4 Supporting threshold in Table 3(b)")
        printed = t3b.loc[(T.NAME[SCORE], T.STRATUM_LABEL[st]), "BP4: Supporting"]
        out[st] = float(r.bp4_threshold_insample)
        out[f"{st}_printed_text"] = printed.removeprefix("≤ ").strip()
        out[f"{st}_printed"] = float(Decimal(out[f"{st}_printed_text"]))
    return out


def _band_of(offset: pd.Series) -> pd.Series:
    return np.where(offset.abs() <= 10, "s3_10", "s11_50")


def printed_difference(d: pd.DataFrame, sup: dict) -> list[str]:
    """The BP7-position variants the printed Table 3(b) values leave out and the fitted
    ones count: every one must be SGE normal and sit exactly on its band's fitted
    threshold, or the stage stops."""
    at = C6_SIDES[-1][1](d["signed_offset"])
    g = d[at].assign(band=lambda x: _band_of(x["signed_offset"]))
    fitted = g["band"].map({st: sup[st] for st in ("s3_10", "s11_50")})
    printed = g["band"].map({st: sup[st + "_printed"] for st in ("s3_10", "s11_50")})
    diff = g[(g[SCORE] <= fitted) & ~(g[SCORE] <= printed)]
    if (diff["y_assay"] != 0).any() or (diff[SCORE] != fitted[diff.index]).any():
        raise SystemExit("[C6] the printed values leave out a damaging variant, or one that is "
                         "not on the fitted threshold: stopped (E19.12)")
    label = {"s3_10": "3-10 nt", "s11_50": "11-50 nt"}
    return [f"{r.gene} {r.hgvs_c} ({label[r.band]})" for r in
            diff.sort_values(["band", "gene"], key=lambda c: c.map({"s3_10": 0, "s11_50": 1})
                             if c.name == "band" else c).itertuples()]


def header(bp4_cut: float, sup: dict, left_out: list[str]) -> str:
    lines = [
        "E19 C6: BP4 at Walker 2023's BP7 positions (src/evid_suppl_c6_c7.py).",
        f"(a) SpliceAI (published basis) at or below the published cut, {bp4_cut!r}.",
        "(b) at or below the band's BP4 Supporting threshold of Table 3(b), the fitted value at full "
        f"precision (evidence_thresholds.csv): 3-10 nt {sup['s3_10']!r}, 11-50 nt {sup['s11_50']!r}.",
        f"With Table 3(b)'s printed values instead (3-10 nt {sup['s3_10_printed_text']}, 11-50 nt "
        f"{sup['s11_50_printed_text']}), the {len(left_out)} SGE normal variants that sit exactly on the "
        f"fitted thresholds would not be counted: {'; '.join(left_out)}. The damaging counts are "
        "the same."]
    return "".join(f"# {x}\n" for x in lines)


def c6(d: pd.DataFrame, bp4_cut: float, sup: dict) -> pd.DataFrame:
    rows = []
    genes = sorted(d["gene"].unique())
    for cond in ("a", "b"):
        for side, at in C6_SIDES:
            for gene in genes + [ALL]:
                g = d if gene == ALL else d[d["gene"] == gene]
                g = g[at(g["signed_offset"])]
                if cond == "a":
                    cut = pd.Series(bp4_cut, index=g.index)
                else:
                    cut = pd.Series(_band_of(g["signed_offset"]), index=g.index).map(
                        {st: sup[st] for st in ("s3_10", "s11_50")})
                m = g[g[SCORE] <= cut]
                n_d, n_n = int((m["y_assay"] == 1).sum()), int((m["y_assay"] == 0).sum())
                rows.append({
                    "condition": cond, "side": side, "gene": gene,
                    "n_at_positions": len(g),
                    "n_damaging_at_positions": int((g["y_assay"] == 1).sum()),
                    "n_normal_at_positions": int((g["y_assay"] == 0).sum()),
                    "n_meeting": len(m), "n_damaging": n_d, "n_normal": n_n,
                    "share_damaging": n_d / len(m) if len(m) else np.nan})
    out = pd.DataFrame(rows)
    out["condition"] = out["condition"].map({
        "a": f"(a) SpliceAI at or below the published cut, {bp4_cut!r}",
        "b": (f"(b) SpliceAI at or below Table 3(b)'s BP4 Supporting threshold: 3-10 nt "
              f"{sup['s3_10']!r}, 11-50 nt {sup['s11_50']!r}")})
    return out.reset_index(drop=True)


def c7(d: pd.DataFrame, pp3: float, bp4: float) -> pd.DataFrame:
    rows = []
    keys = ("lr_pp3", "lr_bp4")
    for side, label, lo, hi in C7_GROUPS:
        g = d[(d["signed_offset"] >= lo) & (d["signed_offset"] <= hi)]
        y = g["y_assay"].to_numpy(dtype=float)
        s = g[SCORE].to_numpy(dtype=float)
        base = {"side": side, "group": label, "offset_from": lo, "offset_to": hi,
                "n": len(g), "n_damaging": int((y == 1).sum()), "n_normal": int((y == 0).sum())}
        st = W._band_stats(y, s, pp3, bp4)
        if not st:
            rows.append(base | {"status": "counts only",
                                "reason": f"fewer than {W.MIN_POS} damaging or normal variants"})
            continue
        ci = W._boot_ci(g, SCORE, pp3, bp4, keys, ("c7", label))
        rows.append(base | {"status": "ok", "reason": "",
                            "n_in_pp3_band": st["n_in_pp3_band"], "lr_pp3": st["lr_pp3"],
                            "lr_pp3_lo": ci["lr_pp3_lo"], "lr_pp3_hi": ci["lr_pp3_hi"],
                            "n_in_bp4_band": st["n_in_bp4_band"], "lr_bp4": st["lr_bp4"],
                            "lr_bp4_lo": ci["lr_bp4_lo"], "lr_bp4_hi": ci["lr_bp4_hi"],
                            "n_genes": ci["n_genes"], "ci_basis": ci.get("ci_basis", "")})
    return pd.DataFrame(rows)


def main() -> None:
    cfg = W.load_config()
    pp3 = cfg["thresholds"]["pp3"]["value"]
    bp4 = cfg["thresholds"]["bp4"]["value"]
    d = labelled(K.load_set())
    sup = bp4_supporting()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t6 = c6(d, bp4, sup)
    with open(OUT_C6, "w", newline="") as f:
        f.write(header(bp4, sup, printed_difference(d, sup)))
        t6.to_csv(f, index=False)
    t7 = c7(d, pp3, bp4)
    t7.to_csv(OUT_C7, index=False)
    s = t6[(t6.gene == ALL) & (t6.side == "both sides")]
    for r in s.itertuples():
        print(f"[C6] {r.condition}: {r.n_meeting} variants, {r.n_damaging} damaging "
              f"({r.share_damaging:.3f})")
    print(t7[["side", "group", "n_damaging", "n_normal", "status", "lr_pp3", "lr_bp4"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
