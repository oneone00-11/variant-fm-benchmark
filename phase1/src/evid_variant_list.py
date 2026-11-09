"""Every analysed variant, one row each, with its HGVS description (E8.5; the author's
decision of 2026-11-07, D5), written by the layout stage under layout/variants/.

What it holds: each intron-side single-nucleotide variant at 1-50 bp from the exon of the
seven training genes (the analysis set) and of the seven external genes, these to the
depth their deposits reach (DDX3X 25 bp, TP53 12 bp). 3-50 bp is the analysis range;
the variants at 1-2 bp, the canonical dinucleotides, are outside it but counted or
analysed in the supplementary tables, and are listed, marked so. A variant beyond 50 bp
(counted by offset only, in the new external genes' table) is not listed, nor is a
control (synonymous or nonsense) that only sets a label. The second BRCA2 dataset
(urn:mavedb:00001242-a-1) would be a group of its own for any variant of it outside the
analysis set that entered a refit; its classes are joined to the analysis set's BRCA2
variants (evid_suppl_c3_brca2), so none does, and the group is empty.

Everything is read from outputs and inputs already written, as they stand: the analysis
set (its c., coordinates, alleles, score set, the deposit's variant URN and the
functional standard's class), DDX3X's splice table, TP53's table as the external test
reads it (evid_external, its primary label) and the deposit copy it came from (each
variant's URN, joined on its c. as the deposit writes it), and each new gene's variant
table (its primary label); transcripts, strands and the chromosomes' accessions from the
MANE v1.5 summary. No coordinate is mapped again and no service is asked.

Before it writes, the list must pass six checks (checks(), the stage stops otherwise),
and they are the tests' too: every c. is an intron-side SNV; its offset's size is the
variant's distance as its source records it, in that distance's band; its bases are the
genome's on the gene's strand; no variant twice in a group; every transcript the MANE
summary's, RefSeq and Ensembl from one MANE Select row; and the counts the author set
(COUNTS), gene by gene, in and out of the range.
"""
from __future__ import annotations

import gzip
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from . import evid_external as X
from . import evid_supp_tables as S

ANALYSIS_SET = Path("data/evidence/analysis_set_v1.parquet")
EXT = Path("data/evidence/external")
DDX3X_TABLE = EXT / "ddx3x_splice.parquet"
NEW_GENE_TABLE = EXT / "{gene}_variants.parquet"
TP53_DEPOSIT = Path("../data/external/tp53_mavedb_scores.tsv")      # from phase1/
MANE = Path("data/evidence/reference_cache/mane/MANE.GRCh38.v1.5.summary.txt.gz")
TRAINING = ["BAP1", "BARD1", "BRCA1", "BRCA2", "PALB2", "RAD51C", "VHL"]
EXTERNAL = ["DDX3X", "TP53"] + X.NEW_GENES
SECOND_BRCA2 = "urn:mavedb:00001242-a-1"
SETS = ["Analysis set", "External gene", "Second BRCA2 dataset"]
INTRONIC_SNV = re.compile(r"^c\.([-*]?\d+)([+-])(\d+)([ACGT])>([ACGT])$")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
CLASS = {1.0: "damaging", 0.0: "normal"}
CANONICAL = "1–2 bp (canonical dinucleotides)"

# the counts the author set (E8.5): in the analysis range, and at 1-2 bp
COUNTS = {
    "in range": {"BAP1": 2268, "BARD1": 1931, "BRCA1": 891, "BRCA2": 551, "PALB2": 1041,
                 "RAD51C": 1152, "VHL": 163, "DDX3X": 1665, "TP53": 240, "CTCF": 1370,
                 "RAD51D": 1296, "SFPQ": 1019, "SBDS": 447, "XRCC2": 390},
    "1-2 bp": {"training genes": 856, "DDX3X": 192, "TP53": 48, "CTCF": 90, "RAD51D": 108,
               "SFPQ": 102, "SBDS": 45, "XRCC2": 24},
    "total": 15889,
}

COLUMNS = ["Set", "Gene", "In the analysis range", "Distance band", "Signed intronic offset",
           "RefSeq transcript (MANE Select)", "HGVS c.", "Ensembl transcript (MANE Select)",
           "GRCh38 chromosome", "GRCh38 position", "Reference allele", "Alternate allele",
           "HGVS g.", "Functional class in the main analysis", "Source (MaveDB score set)",
           "Variant in the source"]


def mane_select(path: Path = MANE) -> dict[str, dict]:
    """The MANE summary's MANE Select rows, by gene symbol."""
    with gzip.open(path, "rt") as f:
        header = f.readline().lstrip("#").rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]
    return {r["symbol"]: r for r in rows if r["MANE_status"] == "MANE Select"}


def _c_part(hgvs: str) -> str:
    """The c. of an HGVS description, without a transcript before it."""
    return hgvs.split(":", 1)[1] if ":" in hgvs else hgvs


def _class(v) -> str:
    return CLASS.get(float(v), "not labelled") if pd.notna(v) else "not labelled"


def _sources() -> list[dict]:
    """Each source's rows in one shape: gene, c., the transcript it is written on, the
    chromosome, position and alleles, the distance and band as the source records them,
    the class, the score set and the variant's URN in it."""
    out = []
    a = pd.read_parquet(ANALYSIS_SET)
    for r in a.itertuples(index=False):
        out.append({"set": SETS[0], "gene": r.gene, "c": r.hgvs_c, "tx": r.transcript,
                    "chrom": str(r.chrom), "pos": int(r.pos), "ref": r.ref, "alt": r.alt,
                    "distance": int(r.intron_offset_abs), "stratum": r.stratum,
                    "class": _class(r.y_assay), "score_set": r.urn, "urn": r.variant_id})
    d = pd.read_parquet(DDX3X_TABLE)
    for r in d.itertuples(index=False):
        out.append({"set": SETS[1], "gene": "DDX3X", "c": r.hgvs_c, "tx": None,
                    "chrom": str(r.chrom), "pos": int(r.pos), "ref": r.ref, "alt": r.alt,
                    "distance": int(r.intron_offset_abs), "stratum": r.stratum,
                    "class": _class(getattr(r, X.PRIMARY_LABEL["DDX3X"])),
                    "score_set": r.scoreset, "urn": r.mavedb_accession})
    t, _ = X._load_external("TP53")
    dep = pd.read_csv(TP53_DEPOSIT, sep="\t", dtype=str, keep_default_na=False)
    urn = dict(zip(dep["HGVS(cDNA)"], dep["accession"]))
    if dep["HGVS(cDNA)"].duplicated().any():
        raise SystemExit(f"[variants] {TP53_DEPOSIT}: a c. stands twice; its URNs cannot be joined on it")
    for r in t.itertuples(index=False):
        u = urn.get(r.hgvs_nt, "")
        out.append({"set": SETS[1], "gene": "TP53", "c": _c_part(r.hgvs_nt),
                    "tx": r.hgvs_nt.split(":", 1)[0], "chrom": str(r.chrom), "pos": int(r.pos),
                    "ref": r.ref, "alt": r.alt, "distance": int(abs(r.intron_offset)),
                    "stratum": r.stratum, "class": _class(getattr(r, X.PRIMARY_LABEL["TP53"])),
                    "score_set": u.split("#", 1)[0], "urn": u})
    for g in X.NEW_GENES:
        n = pd.read_parquet(str(NEW_GENE_TABLE).format(gene=g.lower()))
        n = n[n["offset"].abs() <= 50]
        for r in n.itertuples(index=False):
            out.append({"set": SETS[1], "gene": g, "c": _c_part(r.hgvs_nt), "tx": r.transcript,
                        "chrom": str(r.chrom), "pos": int(r.pos), "ref": r.ref, "alt": r.alt,
                        "distance": int(abs(r.offset)), "stratum": r.stratum,
                        "class": _class(r.label), "score_set": r.urn, "urn": r.variant_urn})
    return out


@lru_cache(maxsize=None)
def _band(gene: str, distance: int) -> str:
    if distance <= 2:
        return CANONICAL
    st = "s3_10" if distance <= 10 else "s11_50"
    if gene in EXTERNAL:
        return f"{S._ext_band(gene, st)} bp"        # cut at the depth the deposit reaches
    return f"{S.STRATUM_LABEL[st]} bp"


def build() -> tuple[pd.DataFrame, list[dict]]:
    """The list, one row per variant, and each source row's fields for the checks."""
    mane = mane_select()
    raw = _sources()
    rows = []
    for s in raw:
        m = mane[s["gene"]]
        c = INTRONIC_SNV.match(s["c"])
        offset = int(c.group(3)) * (1 if c.group(2) == "+" else -1) if c else None
        g_acc = m["GRCh38_chr"]
        rows.append({
            "Set": s["set"], "Gene": s["gene"],
            "In the analysis range": "Yes" if 3 <= s["distance"] <= 50 else "No",
            "Distance band": _band(s["gene"], s["distance"]),
            "Signed intronic offset": "" if offset is None else f"{offset:+d}",
            "RefSeq transcript (MANE Select)": m["RefSeq_nuc"], "HGVS c.": s["c"],
            "Ensembl transcript (MANE Select)": m["Ensembl_nuc"],
            "GRCh38 chromosome": s["chrom"], "GRCh38 position": str(s["pos"]),
            "Reference allele": s["ref"], "Alternate allele": s["alt"],
            "HGVS g.": f"{g_acc}:g.{s['pos']}{s['ref']}>{s['alt']}" if g_acc.startswith("NC_") else "",
            "Functional class in the main analysis": "–" if s["set"] == SETS[2] else s["class"],
            "Source (MaveDB score set)": s["score_set"], "Variant in the source": s["urn"]})
    df = pd.DataFrame(rows, columns=COLUMNS)
    order = {g: k for k, g in enumerate(TRAINING + EXTERNAL)}
    df["_s"], df["_g"] = df.Set.map(SETS.index), df.Gene.map(order)
    df["_p"] = df["GRCh38 position"].astype(int)
    keep = df.sort_values(["_s", "_g", "_p", "Alternate allele"], kind="mergesort").index
    return df.loc[keep].drop(columns=["_s", "_g", "_p"]).reset_index(drop=True), [raw[i] for i in keep]


def counts(df: pd.DataFrame) -> pd.DataFrame:
    """Variants by set, gene, range, band and functional class, in the list's order."""
    keys = ["Set", "Gene", "In the analysis range", "Distance band", "Functional class in the main analysis"]
    c = df.groupby(keys, sort=False).size().reset_index(name="Variants")
    return c.sort_values(keys[:4], kind="mergesort", key=lambda s: s.map(
        {v: k for k, v in enumerate(dict.fromkeys(df[s.name]))})).reset_index(drop=True)


def checks(df: pd.DataFrame, raw: list[dict]) -> dict[str, list]:
    """The six checks (E8.5): what fails each, empty when it passes."""
    mane = mane_select()
    out = {k: [] for k in ("counts", "form", "distance", "bases", "unique", "transcripts")}
    for (_, r), s in zip(df.iterrows(), raw):
        m = INTRONIC_SNV.match(r["HGVS c."])
        if not m:
            out["form"].append(r["HGVS c."])
            continue
        if int(m.group(3)) != s["distance"] or abs(int(r["Signed intronic offset"])) != s["distance"]:
            out["distance"].append((r["Gene"], r["HGVS c."], s["distance"]))
        want_band = ("pm12" if s["distance"] <= 2 else "s3_10" if s["distance"] <= 10 else "s11_50")
        if s["stratum"] not in (want_band, "out_of_scope_deep") or r["Distance band"] != _band(r["Gene"], s["distance"]):
            out["distance"].append((r["Gene"], r["HGVS c."], s["stratum"], r["Distance band"]))
        strand = mane[r["Gene"]]["chr_strand"]
        ref, alt = m.group(4), m.group(5)
        if strand == "-":
            ref, alt = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
        if (ref, alt) != (r["Reference allele"], r["Alternate allele"]):
            out["bases"].append((r["Gene"], r["HGVS c."], r["Reference allele"], r["Alternate allele"], strand))
        mm = mane[r["Gene"]]
        if (r["RefSeq transcript (MANE Select)"], r["Ensembl transcript (MANE Select)"]) != (mm["RefSeq_nuc"], mm["Ensembl_nuc"]) \
                or (s["tx"] is not None and s["tx"] not in (mm["RefSeq_nuc"], mm["Ensembl_nuc"])):
            out["transcripts"].append((r["Gene"], s["tx"]))
    dup = df.duplicated(["Set", "Gene", "GRCh38 position", "Reference allele", "Alternate allele"], keep=False)
    out["unique"] = df.loc[dup, ["Gene", "HGVS c."]].values.tolist()
    inr = df[df["In the analysis range"] == "Yes"].groupby("Gene").size().to_dict()
    pm = df[df["In the analysis range"] == "No"].groupby("Gene").size().to_dict()
    got = {"in range": inr,
           "1-2 bp": {"training genes": sum(pm.get(g, 0) for g in TRAINING)} | {g: pm.get(g, 0) for g in EXTERNAL},
           "total": len(df[df.Set != SETS[2]])}
    for k, want in COUNTS.items():
        if got[k] != want:
            out["counts"].append((k, got[k], want))
    if (df["Distance band"] == CANONICAL).sum() != len(df[df["In the analysis range"] == "No"]):
        out["counts"].append(("out of range without the canonical band",))
    return out


def write(folder: Path) -> list[Path]:
    """The list and its counts, after the six checks; the stage stops if one fails."""
    df, raw = build()
    failed = {k: v[:5] for k, v in checks(df, raw).items() if v}
    if failed:
        raise SystemExit(f"[variants] the list fails its checks: {failed}")
    folder.mkdir(parents=True, exist_ok=True)
    a, b = folder / "analysed_variants.csv", folder / "analysed_variants_counts.csv"
    df.to_csv(a, index=False, lineterminator="\n")
    c = counts(df)
    c["Variants"] = c["Variants"].map(lambda n: f"{n:,}")
    c.to_csv(b, index=False, lineterminator="\n")
    return [a, b]


SOURCES = [str(ANALYSIS_SET), str(DDX3X_TABLE), str(X.TP53_TABLE), str(TP53_DEPOSIT)] + \
    [str(NEW_GENE_TABLE).format(gene=g.lower()) for g in X.NEW_GENES] + [str(MANE)]
