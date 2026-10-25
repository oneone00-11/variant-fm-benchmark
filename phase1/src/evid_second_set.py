"""E19 C3, BRCA2: the second BRCA2 set's label source, Supplementary Table 5 of Sahu et al.
2025, read and checked row by row (rule R6, data/evidence/second_set/brca2_raw/).

The table writes each row twice over: on BRCA2's MANE Select RefSeq transcript
(NM_000059.4(BRCA2):c. ...) and on GRCh38 (NC_000013.11:g. ...). Its rows are joined to
the analysis set by the genomic coordinate (the author's decisions of 2026-10-26), after
the check table5_rows makes of every row. Read by the screening of the second set
(scripts/stage_E19/screen_brca2_second_set.py) and by the stage evid_suppl_c3_brca2; no
predictor score is read here.
"""
from __future__ import annotations

import re
from pathlib import Path

from . import evid_common as K
from . import evid_label_sources as LS

RAW = Path("data/evidence/second_set/brca2_raw")
CLASSES = RAW / "supplementary_table5_classes.tsv"
GENE = "BRCA2"
MANE_REFSEQ = "NM_000059.4"
G_HGVS = re.compile(r"^NC_000013\.11:g\.(\d+)([ACGT])>([ACGT])$")
# the author's decisions of 2026-10-26: Table 5's rows are joined by g.; 28 of them write a
# c. that does not give their g. (the g. being GRCh38's), and the 4 whose g. reference base
# is not GRCh38's are dropped (none is one of the analysis set's variants)
N_MISWRITTEN = 28
DROPPED = ["c.8610A>C", "c.8610A>T", "c.8613A>C", "c.8613A>T"]
REF_CACHE = Path("data/evidence/reference_cache")


def table5_rows(raw: Path = RAW) -> tuple[list[dict], dict]:
    """Supplementary Table 5's rows (the committed extract) with their genomic keys, after
    the row-by-row check of the module's docstring, and the check's record. Every row is
    joined by its g. (the author's decisions of 2026-10-26): a row whose g. reference base
    is not GRCh38's (the cached Ensembl slice) is dropped, and those must be exactly the
    DROPPED rows; a row whose c., mapped as E16 maps coordinates, does not give its g.
    while its g. is GRCh38's is kept, its c. miswritten, and those must be exactly
    N_MISWRITTEN; anything else (an unreadable g., a g. twice, a g. outside the slice)
    stops the stage."""
    import sys
    sys.path.insert(0, str(K.ATLAS_REPO / "src"))
    from atlas import mapping as M
    M.REF_CACHE = REF_CACHE if REF_CACHE.exists() else Path("phase1") / REF_CACHE
    rows = LS.read_tsv(raw / CLASSES.name)
    mane = M.load_mane_records([GENE])[GENE]
    tmap, _ = M.build_transcript_map(MANE_REFSEQ, mane)
    lo = min(e.g_start for e in tmap.exons) - M.REGION_PAD
    hi = max(e.g_end for e in tmap.exons) + M.REGION_PAD
    s0, _, seq = M.fetch_ensembl_region(mane.chr_name, lo, hi)
    rseq = M.RegionSeq(s0, seq)
    kept, dropped, miswritten, keys = [], [], [], set()
    for r in rows:
        tx, _, c = r["hgvs_nt"].partition(":")
        g = G_HGVS.match(r.get("hgvs_g", ""))
        if tx.split("(")[0] != MANE_REFSEQ or g is None:
            raise SystemExit(f"[C3 BRCA2] {r['hgvs_nt']} / {r.get('hgvs_g')}: not on {MANE_REFSEQ} "
                             "and NC_000013.11: stopped")
        pos, ref, alt = int(g.group(1)), g.group(2), g.group(3)
        if not s0 <= pos < s0 + len(seq):
            raise SystemExit(f"[C3 BRCA2] {r['hgvs_g']}: outside the cached slice: stopped")
        r["key"] = (str(mane.chr_name), pos, ref, alt)
        if r["key"] in keys:
            raise SystemExit(f"[C3 BRCA2] {r['hgvs_g']} twice in the table: stopped")
        keys.add(r["key"])
        if seq[pos - s0] != ref:
            dropped.append({"hgvs_nt": r["hgvs_nt"], "hgvs_g": r["hgvs_g"],
                            "grch38_base": seq[pos - s0]})
            continue
        try:
            mp = M.map_variant(tmap, rseq, c)
            gives = f"g.{int(mp[0])}{mp[1]}>{mp[2]}"
            same = (int(mp[0]), mp[1], mp[2]) == (pos, ref, alt)
        except M.MappingFailure as e:
            gives, same = f"not mapped ({e.category})", False
        if not same:
            miswritten.append({"hgvs_nt": r["hgvs_nt"], "hgvs_g": r["hgvs_g"], "c_gives": gives})
        kept.append(r)
    if sorted(x["hgvs_nt"].partition(":")[2] for x in dropped) != sorted(DROPPED):
        raise SystemExit(f"[C3 BRCA2] rows whose g. reference base is not GRCh38's: {dropped}; "
                         f"the author's decision of 2026-10-26 drops only {DROPPED}: stopped")
    if len(miswritten) != N_MISWRITTEN:
        raise SystemExit(f"[C3 BRCA2] {len(miswritten)} rows whose c. does not give their g., not "
                         f"the {N_MISWRITTEN} the author's decision of 2026-10-26 allows: stopped")
    return kept, {"rows": len(rows), "kept": len(kept), "c_gives_g": len(kept) - len(miswritten),
                  "c_miswritten": miswritten, "dropped": dropped}
