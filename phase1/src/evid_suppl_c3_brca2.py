"""E19 C3, item 3 -- BRCA2 against a second dataset: concordance and two refits.

The second dataset is the one scripts/stage_E19/screen_brca2_second_set.py chose by the
author's four conditions (2026-10-25; data/evidence/second_set/brca2_raw/decision.json).
Its per-variant classes are taken by the rule and mapping committed in that folder before
this stage existed (paper_rule.json; the label source audit's order R1, R2, R4, R6): R4
does not apply, so R6, the paper's Supplementary Table 5, its experimental
`Classification` column read from the committed compact extract
(supplementary_table5_classes.tsv), never its ACMG column. Each class goes to abnormal,
normal or intermediate by paper_rule.json's groups; before that, the table's classes,
grouped as the paper's text counts them, must be the text's counts in every class, or
the stage stops. The table writes each row twice over: on BRCA2's MANE Select RefSeq
transcript (NM_000059.4(BRCA2):c. ...) and on GRCh38 (NC_000013.11:g. ...). The classes
are joined to the analysis set's BRCA2 variants on the genomic coordinate (chromosome,
position, reference, alternate), by the author's decisions of 2026-10-26, after a check
of every row (evid_second_set.table5_rows): each row's c., mapped to GRCh38 as E16 maps
coordinates (the atlas mapping module: the Mutalyzer model of NM_000059.4, anchored by
Ensembl's lookup of its MANE Select partner, the reference base checked against the
cached Ensembl slice; all from the tracked reference cache), is set against the row's g.
Twenty-eight rows write a
c. that does not give their g. while their g. is GRCh38's: in the exon 25 donor site
c.9502-c.9505 for c.9501+1 ... +4 (so c.9505A>C/G/T stands twice in the table, once for
the exonic variants at g.32396901; MaveDB's urn:mavedb:00001242-a-1 has none of those
six rows), in the exon 23 acceptor site c.8953-4 ... c.8953 for c.8954-5 ... -1, and one
malformed c.8633-2-1A>T; they are joined by their g. Four rows, c.8610A>C/T and
c.8613A>C/T, write a reference base GRCh38 does not have in either notation, and are
dropped (none is one of the analysis set's variants). Any other disagreement stops the
stage. BRCA2's own classes are its label file's (Huang 2025, Supplementary Table 3,
grouped as the paper groups them), with the intermediate ones that leave a variant
unlabelled.

The second set covers only part of BRCA2's intronic variants at 3-10 nt (coverage_3_10.csv
in the raw folder): both files' headers give the range, as the shared variants span it
(the author's note of 2026-10-26); outside it BRCA2's labels are its own in both refits.

  concordance   the 3 x 3 table of the two classes on the variants both cover at 3-50
                nt (BRCA2 has none at 11-50 nt, so these are its 3-10 nt variants):
                counts and shares of the shared variants.
  refit 1       'only variants the two sets agree on': every BRCA2 variant the two
                classify differently loses its label (one side intermediate counts as
                different); a variant the second set does not cover keeps its own.
  refit 2       'the second set's labels': every BRCA2 variant the second set classes
                takes its class (abnormal 1, normal 0, intermediate no label), the ones
                BRCA2's own file calls intermediate included; one it does not cover
                keeps its own.
Both refits are E3's in-sample fit, cell for cell as evid_suppl_c3 runs it, on the
bands 3-10, 11-50 and 3-50 nt, beside the frozen fit; the 11-50 nt band holds no BRCA2
variant, so its cells cannot change, which the file shows rather than hides.

Outputs, in reports/evidence/suppl_analyses/, each with the second set and its range in
its header (lines beginning '#'):
  c3_brca2_concordance.csv    the 3 x 3 table, with each refit's counts of variants
                              removed and changed
  c3_brca2_second_set.csv     the two refits, in c3_intermediate.csv's form

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_c3_brca2
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import yaml

from . import config as C
from . import evid_common as K
from . import evid_label_sources as LS
from . import evid_second_set as SS
from . import evid_interval_lr as E3
from . import evid_suppl_c3 as C3
from .evid_build_set import _norm_hgvs

RAW = SS.RAW
CLASSES = SS.CLASSES
RULE = RAW / "paper_rule.json"
DECISION = RAW / "decision.json"
SIGNED = re.compile(r"^c\.[-*]?\d+([+-])(\d+)[ACGT]>[ACGT]$")
G_HGVS = SS.G_HGVS
OUT_DIR = C3.OUT_DIR
OUT_CONCORDANCE = OUT_DIR / "c3_brca2_concordance.csv"
OUT_REFITS = OUT_DIR / "c3_brca2_second_set.csv"
GENE = "BRCA2"
ORDER = ["abnormal", "normal", "intermediate"]
VALUE = {"abnormal": 1.0, "normal": 0.0, "intermediate": None}
AGREE = "only variants the two sets agree on"
SECOND = "the second set's labels"


def own_classes(df: pd.DataFrame) -> pd.Series:
    """BRCA2's classes as its label file gives them, on the set's BRCA2 rows: the class
    column through the gene's map, named abnormal / normal / intermediate."""
    spec = C.ASSAY_LABELS[GENE]
    lab = K.read_back(C.ASSAY_LABEL_DIR / spec["file"], sep="\t", dtype=str,
                      keep_default_na=False)
    names = {k: ("intermediate" if v is None else "abnormal" if v == 1.0 else "normal")
             for k, v in spec["map"].items()}
    lut = dict(zip(lab[spec["key"]].map(_norm_hgvs), lab[spec["class_col"]].map(names)))
    mask = df["gene"] == GENE
    return df.loc[mask, "hgvs_c"].map(_norm_hgvs).map(lut)


def second_classes(df: pd.DataFrame) -> tuple[pd.Series, dict]:
    """The second set's class of each of the set's BRCA2 variants (abnormal, normal,
    intermediate; none where the table has no row), by paper_rule.json's R6 rule, joined
    on the genomic coordinate; and table5_rows' record."""
    rule = json.loads(RULE.read_text())
    sup = rule.get("supplementary_rule") or {}
    if rule.get("status") != "not_stated" or sup.get("rule") != "R6":
        raise SystemExit(f"[C3 BRCA2] {RULE}: not the R6 rule this stage reads: stopped")
    every = LS.read_tsv(CLASSES)
    unknown = sorted({r["category"] for r in every} - set(sup["groups"]))
    if unknown or set(sup["groups"].values()) - set(ORDER):
        raise SystemExit(f"[C3 BRCA2] classes the mapping does not name: {unknown}: stopped")
    # the paper counts every SNV of the table, before any row is dropped
    paper = pd.Series([sup["paper_groups"][r["category"]] for r in every]).value_counts().to_dict()
    want = rule["reported_counts"]
    if paper != want["counts"] or len(every) != want["total"]:
        raise SystemExit(f"[C3 BRCA2] the table's classes, grouped as the paper counts them, are "
                         f"{paper}, not the paper's {want['counts']}: stopped")
    rows, check = SS.table5_rows()
    lut = {r["key"]: sup["groups"][r["category"]] for r in rows}
    b = df[df["gene"] == GENE]
    keys = zip(b["chrom"].astype(str), b["pos"].astype(int), b["ref"], b["alt"])
    return pd.Series([lut.get(k) for k in keys], index=b.index, dtype=object), check


def covered_range(df: pd.DataFrame, shared: pd.Series) -> dict[str, tuple[int, int]]:
    """The signed offsets the shared variants span, per side."""
    out = {}
    for c in df.loc[shared[shared].index, "hgvs_c"]:
        m = SIGNED.match(c)
        side = "donor" if m.group(1) == "+" else "acceptor"
        o = int(m.group(2))
        lo, hi = out.get(side, (o, o))
        out[side] = (min(lo, o), max(hi, o))
    return out


def _blocks(items: list[dict]) -> list[list[dict]]:
    """Rows in runs of adjacent genomic positions."""
    items = sorted(items, key=lambda x: (int(G_HGVS.match(x["hgvs_g"]).group(1)), x["hgvs_g"]))
    out = []
    for x in items:
        p = int(G_HGVS.match(x["hgvs_g"]).group(1))
        if out and p - int(G_HGVS.match(out[-1][-1]["hgvs_g"]).group(1)) <= 1:
            out[-1].append(x)
        else:
            out.append([x])
    return out


def header(urn: str, rng: dict, n_shared: int, n_ours: int, check: dict,
           n_mavedb: int, ours: dict) -> str:
    d, a = rng.get("donor"), rng.get("acceptor")

    def g_of(x):
        m = G_HGVS.match(x["hgvs_g"])
        return (str(13), int(m.group(1)), m.group(2), m.group(3))

    lines = [
        "E19 C3, BRCA2 against a second dataset (src/evid_suppl_c3_brca2.py).",
        f"Second set: {urn} (Sahu et al. 2025), classes by rule R6: Supplementary Table 5's "
        "Classification, abnormal / normal / intermediate by side (paper_rule.json).",
        f"Joined on the genomic coordinate (the author's decisions of 2026-10-26). Table 5's "
        f"{check['rows']:,} rows were checked one by one: each row's c., mapped as E16 maps "
        f"coordinates, gives its g. on {check['c_gives_g']:,}; {len(check['c_miswritten'])} write "
        "a c. that does not, their g. being GRCh38's, and are joined by g.; "
        f"{len(check['dropped'])} write a reference base GRCh38 does not have and are dropped."]
    for blk in _blocks(check["c_miswritten"]):
        gs = [G_HGVS.match(x["hgvs_g"]).group(1) for x in blk]
        where = f"g.{gs[0]}" if gs[0] == gs[-1] else f"g.{gs[0]}-g.{gs[-1]}"
        lines.append(f"  c. miswritten, {where} ({len(blk)} row{'s' if len(blk) > 1 else ''}): "
                     + "; ".join(
            f"{x['hgvs_nt'].partition(':')[2]} at {x['hgvs_g'].partition(':')[2]}"
            + (f" = the set's {ours[g_of(x)]}" if g_of(x) in ours else " (not in the set)")
            for x in blk))
    lines.append("  dropped, the reference base not GRCh38's: " + "; ".join(
        f"{x['hgvs_nt'].partition(':')[2]} at {x['hgvs_g'].partition(':')[2]} (GRCh38 "
        f"{x['grch38_base']})" for x in check["dropped"]))
    lines += [
        "c.9505A>C/G/T thus stands twice in Table 5 (for g.32394937, the set's c.9501+4A>C/G/T, and "
        f"for the exonic g.32396901); MaveDB's {urn} lacks all six rows.",
        f"Shared with the analysis set's {n_ours} BRCA2 variants at 3-10 nt: {n_shared} by Table 5 "
        f"(the screening record's {n_mavedb}, by MaveDB's score table, is unchanged).",
        f"The second set covers them only on the donor side +{d[0]} ... +{d[1]} and the acceptor "
        f"side -{a[0]} ... -{a[1]} (coverage_3_10.csv); elsewhere BRCA2's labels are its own in "
        "both refits.",
        "BRCA2 has no variant at 11-50 nt, so that band cannot change."]
    return "".join(f"# {x}\n" for x in lines)


def main() -> None:
    dec = json.loads(DECISION.read_text())
    if not dec["outcome"].startswith("use "):
        raise SystemExit(f"[C3 BRCA2] {DECISION}: {dec['outcome']!r}: no second set to use")
    cfg = yaml.safe_load(E3.CONFIG_PATH.read_text())
    path_cuts, ben_cuts = K.acmg_bands(cfg)
    cuts = {"path": path_cuts, "ben": ben_cuts}
    df = K.load_set()
    df[K.FUSION] = K.logo_fusion(df, [t for t in K.FUSION_FEATURES if t in df.columns])
    tools = K.panel_of(df) + [K.FUSION]
    frozen = C3.frozen_cells()
    C3.check_reproduces_frozen(C3.fit_cells(df, df["y_assay"].astype("Float64"), tools, cuts),
                               frozen)

    own = own_classes(df)
    sec, check = second_classes(df)
    if own.isna().any():
        raise SystemExit("[C3 BRCA2] a BRCA2 variant of the set is not in its label file")
    in_scope = df.loc[own.index, "stratum"] != "pm12"
    shared = sec.notna() & in_scope
    tab = pd.crosstab(own[shared], sec[shared]).reindex(index=ORDER, columns=ORDER, fill_value=0)

    y_agree = df["y_assay"].astype("Float64").copy()
    differ = shared & (own != sec)
    y_agree.loc[differ[differ].index] = pd.NA
    y_second = df["y_assay"].astype("Float64").copy()
    take = sec.notna()
    y_second.loc[take[take].index] = sec[take].map(VALUE).astype("Float64")

    parts = []
    for name, y in ((AGREE, y_agree), (SECOND, y_second)):
        parts.append(C3.compare(name, C3.fit_cells(df, y, tools, cuts), frozen,
                                C3.label_changes(df, y)))
    refits = pd.concat(parts, ignore_index=True)
    n_shared = int(shared.sum())
    rows = [{"huang_class": i, "second_set_class": j, "n": int(tab.loc[i, j]),
             "share_of_shared": tab.loc[i, j] / n_shared if n_shared else float("nan")}
            for i in ORDER for j in ORDER]
    conc = pd.DataFrame(rows)
    ch = refits.drop_duplicates(["analysis", "stratum"]).set_index(["analysis", "stratum"])
    conc["second_set"] = dec["outcome"].removeprefix("use ")
    conc["n_shared_3_50"] = n_shared
    conc["n_agree"] = int((shared & (own == sec)).sum())
    for name, short in ((AGREE, "agree"), (SECOND, "second")):
        for col in ("n_added", "n_removed", "n_switched"):
            conc[f"{short}_refit_{col}_3_50"] = int(ch.loc[(name, "s3_50"), col])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    n_ours = int(((df["gene"] == GENE) & (df["stratum"] == "s3_10")).sum())
    scr = K.read_back(RAW / "screening.csv", dtype=str, keep_default_na=False)
    n_mavedb = int(scr.loc[scr["urn"] == conc["second_set"].iloc[0], "c_shared_3_10"].iloc[0])
    b = df[df["gene"] == GENE]
    ours = {(str(c), int(p), r, a): h for c, p, r, a, h in
            zip(b["chrom"], b["pos"], b["ref"], b["alt"], b["hgvs_c"])}
    head = header(conc["second_set"].iloc[0], covered_range(df, shared), n_shared, n_ours,
                  check, n_mavedb, ours)
    for path, table in ((OUT_CONCORDANCE, conc), (OUT_REFITS, refits)):
        with open(path, "w", newline="") as f:
            f.write(head)
            table.to_csv(f, index=False)
    print(f"[C3 BRCA2] {n_shared} shared variants at 3-50 nt; {conc.n_agree.iloc[0]} agree")
    print(tab.to_string())
    for name in (AGREE, SECOND):
        c = refits[(refits.analysis == name) & (refits.tier_change != "unchanged")]
        print(f"[C3 BRCA2] {name}: {len(c)} cells change tier"
              + "".join(f"\n      {r.tool} {r.stratum} {r.side} {r.tier} {r.tier_change}"
                        for r in c.itertuples()))


if __name__ == "__main__":
    main()
