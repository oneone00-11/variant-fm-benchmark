"""
phase4a_score_tp53.py -- scoring orchestrator for TP53 splice variants (external validation).

Goal: score TP53's 192 intron-side splice SNVs (Funk et al. 2024 SGE, NM_000546.6)
and emit a row set in the SAME schema as the benchmark's score_matrix_final.tsv, so
it can pass through the identical phase-1 freeze to sit alongside the seven genes.

DESIGN RULE (the whole point of external validation): TP53 must be scored by the
EXACT SAME code, parameters and versions used for the seven genes. This file writes
NO new predictor logic — it only (A) loads the MaveDB data and (B) dispatches each of
the ten predictors to the benchmark's own scorer. If a scorer here diverges from the
seven-gene scorer, the external comparison is invalid.

There is NO single `score_all_predictors()` in this benchmark: the seven genes were
scored by ten heterogeneous mechanisms (per-predictor scripts + external SpliceAI /
Pangolin CLI + a GPU NT run + the AlphaGenome API + a precomputed GPN-MSA lookup).
HOOK B therefore dispatches to each of those exact entry points (see BENCHMARK_SCORERS).

Run (HOOK A only, prints schema + scorer provenance; does NOT score):
    python scripts/phase4a_score_tp53.py
"""
from __future__ import annotations
import csv
import re
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
MAVEDB_TSV = ROOT / "data" / "external" / "tp53_mavedb_scores.tsv"
OUT_TSV = ROOT / "data" / "external" / "tp53_splice_scored.tsv"   # written by HOOK B (later)

GENE = "TP53"
CHROM = "17"                       # TP53, GRCh38 (NC_000017.11); g.-coords are chr17 absolute
SOURCE_URN = "urn:mavedb:00001213-a-1"
TRANSCRIPT = "NM_000546.6"
SPLICE_REGION_MAX = 8              # |offset| <= 8 (matches the frozen-matrix splice window)

# The ten predictor columns of score_matrix_final.tsv (frozen-matrix raw input schema).
PREDICTOR_COLS = ["spliceai_ds", "pangolin_score", "alphagenome_splice", "gpn_msa_score",
                  "nucleotide_transformer", "cadd_phred", "alphamissense",
                  "phylop100way", "phastcons100way", "gnomad_af_global"]


# ===========================================================================
# HOOK A -- data: load the MaveDB TP53 splice variants into the raw schema
# ===========================================================================
# Column mapping CONFIRMED against data/external/tp53_mavedb_scores.tsv header:
#   accession, hgvs_nt, hgvs_splice, hgvs_pro, score, mut_ID,
#   HGVS(genomic), HGVS(cDNA), HGVS(protein), score_RFS_enrich2, SE_RFS_enrich2
# --> frozen `functional_score` <- MaveDB `score`
# --> frozen `hgvs_nt`          <- MaveDB `HGVS(cDNA)`  (NM_000546.6:c...)
#     NB: the MaveDB column literally named `hgvs_nt` is LOCAL n.-notation
#     (e.g. n.12903T>G) and is NOT the transcript key -- do not use it.
# --> chrom/pos/ref/alt         <- parsed from MaveDB `HGVS(genomic)` g.-notation
CDNA_COL, SCORE_COL, GENOMIC_COL = "HGVS(cDNA)", "score", "HGVS(genomic)"

_INTRON_SNV = re.compile(r"c\.\d+([+-]\d+)[ACGT]>[ACGT]$")      # intron-side substitution + signed offset
_G_SNV = re.compile(r"g\.(\d+)([ACGT])>([ACGT])$")             # genomic substitution (plus strand)


def load_tp53_splice() -> pd.DataFrame:
    """Load the 192 intron-side splice SNVs (|offset| <= 8) with a functional score,
    in the score_matrix_final.tsv raw schema. Predictor columns are left NA for HOOK B."""
    if not MAVEDB_TSV.exists():
        sys.exit(f"[phase4a] MaveDB TSV not found: {MAVEDB_TSV} (run HOOK-A download first)")
    df = pd.read_csv(MAVEDB_TSV, sep="\t", dtype=str)

    rows = []
    for _, r in df.iterrows():
        cdna = str(r[CDNA_COL]).split(":")[-1]           # strip 'NM_000546.6:'
        m = _INTRON_SNV.search(cdna)
        if not m or abs(int(m.group(1))) > SPLICE_REGION_MAX:
            continue                                     # keep only intron-side SNV, |offset| <= 8
        g = _G_SNV.search(str(r[GENOMIC_COL]))
        if not g:
            continue                                     # need clean genomic SNV coords
        score = str(r[SCORE_COL]).strip()
        if score in ("", "NA", "nan", "None"):
            continue                                     # must carry a functional score
        offset = int(m.group(1))
        pos, ref, alt = g.group(1), g.group(2), g.group(3)
        rows.append({
            "variant_id": f"{CHROM}-{pos}-{ref}-{alt}",
            "gene": GENE, "chrom": CHROM, "pos": pos, "ref": ref, "alt": alt,
            "hgvs_nt": str(r[CDNA_COL]),                 # transcript c.-notation (join key)
            "functional_score": score,
            "intron_offset": float(offset),
            "region_class": "splice",
            "splice_class": "splice_core" if abs(offset) <= 2 else "splice_region",
            "source_urn": SOURCE_URN, "assay_type": "SGE",
            **{c: pd.NA for c in PREDICTOR_COLS},         # <- filled by HOOK B
        })
    out = pd.DataFrame(rows)
    return out


# ===========================================================================
# HOOK B -- scoring: dispatch every predictor to the benchmark's OWN scorer
# ===========================================================================
# Each entry names the EXACT script + entry point that scored the seven genes.
# `mode`: "function" = import the sibling script and call this function on the
#         TP53 (chrom,pos,ref,alt) rows; "external" = build the same VCF and run
#         the same external CLI, then parse with the benchmark's parser; "lookup"
#         = query the same precomputed/remote source with the benchmark's code.
# NOTHING here re-implements scoring; it points at the seven-gene code.
BENCHMARK_SCORERS = {
    "spliceai_ds": dict(
        script="80_make_vcf.py (+ external SpliceAI CLI) -> 82_assemble_v2.py",
        entry="parse_spliceai()", mode="external",
        detail="VCF built exactly as 80_make_vcf.py; SpliceAI run externally on it; "
               "score = max of the four SpliceAI delta scores (parse_spliceai)."),
    "pangolin_score": dict(
        script="80_make_vcf.py (+ external Pangolin CLI) -> 82_assemble_v2.py",
        entry="parse_pangolin()", mode="external",
        detail="same VCF; Pangolin run externally; score = max(splice gain, |splice loss|)."),
    "alphagenome_splice": dict(
        script="91_score_alphagenome.py", entry="score_one(v)", mode="function",
        detail="alphagenome dna_client; SEQUENCE_LENGTH_1MB interval; variant_scorers = the "
               "SPLICE_* recommended scorers; aggregate = max |raw_score| (needs ALPHAGENOME_API_KEY)."),
    "gpn_msa_score": dict(
        script="90_score_gpn.py", entry="tabix query of songlab/gpn-msa-hg38-scores", mode="lookup",
        detail="pysam/tabix lookup of the precomputed scores.tsv.bgz by (chrom,pos,ref,alt) "
               "(same HuggingFace release + local .tbi used for the seven genes)."),
    "nucleotide_transformer": dict(
        script="92_score_nt.py", entry="load_model(); masked_llr(...)", mode="function",
        detail="masked-token LLR = logP(REF 6-mer) - logP(ALT 6-mer); WINDOW_BP = 6000; "
               "non-overlapping 6-mer tokenizer (needs GPU + GRCh38 FASTA)."),
    "cadd_phred": dict(
        script="75_cadd.py", entry="fetch(c,p,r,a)", mode="function",
        detail="CADD REST API, build GRCh38-v1.7 (PHRED)."),
    "alphamissense": dict(
        script="74_alphamissense.py", entry="download() + position lookup", mode="lookup",
        detail="AlphaMissense_hg38.tsv.gz (Zenodo 8208688) lookup; ~0 splice coverage expected."),
    "phylop100way": dict(
        script="73_conservation.py", entry="fetch_track('phyloP100way', ...)", mode="function",
        detail="UCSC getData/track API, phyloP100way."),
    "phastcons100way": dict(
        script="73_conservation.py", entry="fetch_track('phastCons100way', ...)", mode="function",
        detail="UCSC getData/track API, phastCons100way."),
    "gnomad_af_global": dict(
        script="72_gnomad.py", entry="fetch_gene('TP53')", mode="function",
        detail="gnomAD GraphQL API; global AF per variant (same DATASET as the seven genes)."),
}


def _load_benchmark_module(filename: str):
    """Import a benchmark scorer script by filename (numbered names aren't importable
    normally). scripts/ is put on sys.path so their `from config import ...` resolves."""
    import importlib.util
    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    path = SCRIPTS / filename
    spec = importlib.util.spec_from_file_location(f"bench_{path.stem}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # NOTE: only call for scripts whose module level is side-effect free
    return mod


def provenance() -> str:
    lines = ["HOOK B -- predictor -> benchmark scorer (the exact seven-gene code):"]
    for col, s in BENCHMARK_SCORERS.items():
        lines.append(f"  {col:24s} <- scripts/{s['script']} :: {s['entry']}  [{s['mode']}]")
        lines.append(f"      {s['detail']}")
    return "\n".join(lines)


# --- GPN-MSA pinned revision (byte-stable since 2024-03-02; verified 2026-07-10) --
GPN_PIN = "cf1718a926db35d078bcb08b9dcebc803494bf37"
GPN_RESOLVE = ("https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/"
               f"resolve/{GPN_PIN}/scores.tsv.bgz")
GPN_TBI = str(ROOT / "refs" / "gpn_scores.tsv.bgz.tbi")
AM_GZ = ROOT / "data" / "raw" / "scores" / "AlphaMissense_hg38.tsv.gz"

# External predictors: version-pinned tools that ran in WSL/API for the 7 genes.
# NT is EXCLUDED by decision (2026-07-10: original code lost, window unrecoverable).
# The remaining three are scored in that env and dropped here as cache TSVs
# (schema: chrom<TAB>pos<TAB>ref<TAB>alt<TAB><col>); read by _read_external().
# TP53 uses an 8-FEATURE set: NT and Pangolin are EXCLUDED (their columns stay NA).
# AlphaGenome is read from its API cache TSV; SpliceAI is parsed directly from the
# SpliceAI output VCF via 82.parse_spliceai (see _score_spliceai).
EXTERNAL_CACHE = {
    "alphagenome_splice": ROOT / "data" / "external" / "tp53_alphagenome.tsv",
}
SPLICEAI_VCF_TP53 = ROOT / "data" / "external" / "tp53_spliceai.vcf"
FROZEN = ROOT / "phase1" / "data" / "frozen" / "frozen_matrix_v1.parquet"


def _keys(df):
    return [(str(r.chrom), int(r.pos), str(r.ref), str(r.alt)) for r in df.itertuples()]


# ---- feasible-on-this-Mac scorers: reuse the EXACT seven-gene entry points ----
def _score_gnomad(df):                                     # 72_gnomad.fetch_gene
    m = _load_benchmark_module("72_gnomad.py")
    af = {}
    for v in m.fetch_gene(GENE):
        c, p, ref, alt = v["variant_id"].split("-")
        seqs = [x for x in (v.get("genome"), v.get("exome")) if x and x.get("af") is not None]
        af[(c, int(p), ref, alt)] = max(x["af"] for x in seqs) if seqs else None
    return {"gnomad_af_global": [af.get(k) for k in _keys(df)]}


def _score_conservation(df):                               # 73_conservation.fetch_track
    m = _load_benchmark_module("73_conservation.py")
    pos = df["pos"].astype(int)
    lo, hi = int(pos.min()), int(pos.max())
    cols = {}
    for track, col in (("phyloP100way", "phylop100way"), ("phastCons100way", "phastcons100way")):
        p2v = {}
        for it in m.fetch_track(track, CHROM, lo - 1, hi):
            v, s, e = it.get("value"), it.get("start"), it.get("end")
            if v is None or s is None:
                continue
            for p in range(int(s) + 1, int(e) + 1):        # 0-based half-open -> 1-based
                p2v[p] = v
        cols[col] = [p2v.get(p) for p in pos]
    return cols


def _score_cadd(df):                                       # 75_cadd.fetch (cached, gentle)
    import json
    from concurrent.futures import ThreadPoolExecutor
    m = _load_benchmark_module("75_cadd.py")
    cache_path = ROOT / "data" / "external" / "tp53_cadd_cache.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    keys = _keys(df)
    kstr = lambda k: f"{k[0]}:{k[1]}:{k[2]}:{k[3]}"

    def one(k):
        res = m.fetch(k[0], k[1], k[2], k[3])
        return k, (res[0] if isinstance(res, tuple) else None)

    # the CADD public API throttles bursts -> gentle cached passes until saturated
    for p in range(6):
        todo = [k for k in keys if cache.get(kstr(k)) is None]
        if not todo:
            break
        print(f"  [cadd] pass {p}: {len(todo)} to (re)fetch (gentle)", flush=True)
        with ThreadPoolExecutor(max_workers=3) as ex:
            for k, phred in ex.map(one, todo):
                cache[kstr(k)] = phred
        cache_path.write_text(json.dumps(cache))
        time.sleep(2)
    return {"cadd_phred": [cache.get(kstr(k)) for k in keys]}


def _score_alphamissense(df):                              # 74 stream-filter (exact logic)
    import gzip
    want = {f"{c}:{p}:{r}:{a}" for (c, p, r, a) in _keys(df)}
    pos = df["pos"].astype(int)
    lo, hi = int(pos.min()), int(pos.max())
    chrom = f"chr{CHROM}"
    got = {}
    with gzip.open(AM_GZ, "rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("CHROM"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[0] != chrom:
                continue
            pp = int(p[1])
            if pp < lo or pp > hi:
                continue
            if f"{CHROM}:{pp}:{p[2]}:{p[3]}" in want:
                got[(CHROM, pp, p[2], p[3])] = float(p[8])
    return {"alphamissense": [got.get(k) for k in _keys(df)]}


def _gpn_url():
    import subprocess
    return subprocess.check_output(
        ["curl", "-s", "-o", "/dev/null", "-w", "%{url_effective}", "-I", "-L", GPN_RESOLVE]
    ).decode().strip()


def _score_gpn(df):                                        # 90 tabix (pinned revision)
    import pysam
    url = _gpn_url()
    pos = df["pos"].astype(int)
    lo, hi = int(pos.min()), int(pos.max())
    keyset = set(_keys(df))
    score = {}
    for attempt in range(5):
        try:
            tbx = pysam.TabixFile(url, index=GPN_TBI)
            for line in tbx.fetch(CHROM, lo - 1, hi + 1):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                k = (parts[0], int(parts[1]), parts[2], parts[3])
                if k in keyset:
                    score[k] = float(parts[4])
            tbx.close()
            break
        except Exception as e:
            print(f"  [gpn] attempt {attempt} failed: {str(e)[:60]}; re-resolving", flush=True)
            time.sleep(3)
            try:
                url = _gpn_url()
            except Exception:
                pass
    return {"gpn_msa_score": [score.get(k) for k in _keys(df)]}


def _score_spliceai(df):                                   # 82.parse_spliceai on the SpliceAI VCF
    if not SPLICEAI_VCF_TP53.exists():
        return {"spliceai_ds": [None] * len(df)}
    m = _load_benchmark_module("82_assemble_v2.py")
    m.SPLICEAI_VCF = str(SPLICEAI_VCF_TP53)                 # point the exact seven-gene parser at the TP53 VCF
    parsed = m.parse_spliceai()                            # -> df: chrom,pos,ref,alt,spliceai_ds (+components)
    sc = {(str(r.chrom), int(r.pos), str(r.ref), str(r.alt)): float(r.spliceai_ds)
          for r in parsed.itertuples()}
    return {"spliceai_ds": [sc.get(k) for k in _keys(df)]}


def _read_external(df):
    cols = {}
    for col, path in EXTERNAL_CACHE.items():
        if path.exists():
            m = {}
            with open(path) as f:
                for d in csv.DictReader(f, delimiter="\t"):
                    val = d.get(col)
                    m[(d["chrom"], int(d["pos"]), d["ref"], d["alt"])] = (
                        float(val) if val not in (None, "", "NA", "nan") else None)
            cols[col] = [m.get(k) for k in _keys(df)]
            print(f"[phase4a] external {col}: {sum(x is not None for x in cols[col])}/{len(df)} "
                  f"from {path.name}", flush=True)
        else:
            cols[col] = [None] * len(df)
            print(f"[phase4a] external {col}: PENDING -> run in WSL, drop {path.name}", flush=True)
    return cols


FEASIBLE = [_score_gnomad, _score_conservation, _score_cadd, _score_alphamissense,
            _score_gpn, _score_spliceai]


def score_with_benchmark(df: pd.DataFrame) -> pd.DataFrame:
    """Fill the predictor columns for the TP53 rows by dispatching to the benchmark's
    OWN scorers. Nothing here re-implements scoring. Six run via their exact seven-gene
    entry points (gnomAD/conservation/CADD/AlphaMissense/GPN + SpliceAI parsed from its
    VCF by 82.parse_spliceai); AlphaGenome is read from its API cache TSV. NT and
    Pangolin are EXCLUDED (8-feature TP53 set) -> their columns stay NA."""
    df = df.copy()
    for col in PREDICTOR_COLS:
        df[col] = pd.NA
    for fn in FEASIBLE:
        for col, vals in fn(df).items():
            df[col] = pd.Series(vals, index=df.index, dtype="object")
        print(f"[phase4a] {fn.__name__}: done", flush=True)
    for col, vals in _read_external(df).items():
        df[col] = pd.Series(vals, index=df.index, dtype="object")
    return df                                              # nucleotide_transformer + pangolin_score stay NA (excluded)


# ===========================================================================
# Assemble into the frozen-matrix schema by reusing the EXACT phase-1 transforms
# ===========================================================================
def to_frozen_schema(scored: pd.DataFrame) -> pd.DataFrame:
    p1 = ROOT / "phase1"
    if str(p1) not in sys.path:
        sys.path.insert(0, str(p1))
    from src import phase1_build_frozen_matrix as P1
    from src import config as C1

    d = scored.copy()
    for src in ("mc_terms", "clnsig_class"):               # src cols TP53 lacks -> NA
        if src not in d.columns:
            d[src] = pd.NA
    numeric = ["pos", "functional_score", "cadd_phred", "alphamissense", "gnomad_af_global",
               "phylop100way", "phastcons100way", "spliceai_ds", "pangolin_score",
               "alphagenome_splice", "gpn_msa_score", "nucleotide_transformer", "intron_offset"]
    for c in numeric:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")

    rename = {src: canon for canon, src in C1.COLUMNS.items() if src is not None}
    d = d.rename(columns=rename)[list(dict.fromkeys(rename.values()))].copy()   # == load_raw, no gene filter
    d["gene"] = pd.Categorical(d["gene"])                  # dtype parity (category); TP53 kept

    d = P1.orient_functional(d)                            # func_pathogenicity = +func_score for TP53 (in FLIP_GENES)
    d = P1.derive_intron_offset(d)                         # offset_source, splice_side
    d = P1.assign_region(d)                                # region/is_splice/splice_subclass/bin/analysis_region
    d = P1.classes_from_clinvar(d)                         # clinvar all-NA -> y_clinvar NA
    d["y_assay"] = pd.array([pd.NA] * len(d), dtype="Float64")   # no official TP53 assay labels ingested
    d["y_assay_source"] = pd.NA
    d = P1.add_missingness_flags(d)

    frozen_cols = list(pd.read_parquet(FROZEN).columns)
    for c in frozen_cols:
        if c not in d.columns:
            d[c] = pd.NA
    return d[frozen_cols].sort_values("variant_id").reset_index(drop=True)


def schema_parity(tp53: pd.DataFrame):
    ref = pd.read_parquet(FROZEN)
    rc, tc = list(ref.columns), list(tp53.columns)
    dtype_mismatch = {c: (str(ref[c].dtype), str(tp53[c].dtype))
                      for c in rc if c in tp53 and str(ref[c].dtype) != str(tp53[c].dtype)}
    return {"columns_identical": rc == tc, "n_cols_ref": len(rc), "n_cols_tp53": len(tc),
            "missing_in_tp53": [c for c in rc if c not in tc],
            "extra_in_tp53": [c for c in tc if c not in rc],
            "dtype_mismatches": dtype_mismatch}


FEATURES8 = ["spliceai", "alphagenome", "gpn_msa", "cadd",
             "alphamissense", "phylop", "phastcons", "gnomad_af"]


def feature_parity(tp53: pd.DataFrame):
    """Schema parity restricted to the 8 active TP53 feature columns (name + dtype)."""
    ref = pd.read_parquet(FROZEN)
    rows, all_ok = [], True
    for c in FEATURES8:
        in_both = (c in ref.columns) and (c in tp53.columns)
        dt_ref = str(ref[c].dtype) if c in ref.columns else "MISSING"
        dt_tp = str(tp53[c].dtype) if c in tp53.columns else "MISSING"
        ok = in_both and dt_ref == dt_tp
        all_ok = all_ok and ok
        rows.append((c, in_both, dt_ref, dt_tp, ok))
    return rows, all_ok


# expected coverage on the 192 intron-side splice SNVs (8-feature TP53: NT + Pangolin excluded)
EXPECT = {"spliceai": "full", "alphagenome": "full", "gpn_msa": "full",
          "cadd": "full", "phylop": "full", "phastcons": "full",
          "alphamissense": "~zero", "gnomad_af": "partial",
          "pangolin": "excluded", "nt": "excluded"}
EXCLUDED_SHORT = {"nt", "pangolin"}
MUST_BE_FULL = {"spliceai", "alphagenome"}    # user's hard stop-condition


def coverage_report(tp53: pd.DataFrame):
    n = len(tp53)
    rows, flags = [], []
    for short in ["spliceai", "alphagenome", "gpn_msa", "cadd", "alphamissense",
                  "phylop", "phastcons", "gnomad_af", "pangolin", "nt"]:
        nn = int(tp53[short].notna().sum())
        frac = nn / n
        exp = EXPECT[short]
        status = "excluded" if short in EXCLUDED_SHORT else "ok"
        if short in MUST_BE_FULL and frac < 1.0:
            flags.append(f"STOP {short}: {frac:.0%} -- must be 100% for the TP53 8-feature matrix")
        elif short not in EXCLUDED_SHORT and exp == "full" and frac < 0.80:
            flags.append(f"{short}: {frac:.0%} (expected ~full) -- possible scoring-convention mismatch")
        if exp == "~zero" and frac > 0.20:
            flags.append(f"{short}: {frac:.0%} (expected ~0 on intron-side splice) -- check")
        rows.append((short, nn, n, frac, exp, status))
    return rows, flags


def write_tp53_vcf(df, path):
    """Minimal VCF (no-chr contig '17', matching the 7-gene convention) for the
    WSL SpliceAI/Pangolin/AlphaGenome step."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=17>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in df.sort_values("pos").itertuples():
            f.write(f"{r.chrom}\t{r.pos}\t.\t{r.ref}\t{r.alt}\t.\t.\t.\n")
    return path


# ===========================================================================
def main() -> None:
    tp = load_tp53_splice()
    print(f"[phase4a] HOOK A: {len(tp)} TP53 intron-side splice SNVs "
          f"({tp['splice_class'].value_counts().to_dict()})", flush=True)
    write_tp53_vcf(tp, ROOT / "data" / "external" / "tp53_variants.vcf")

    scored = score_with_benchmark(tp)
    frozen = to_frozen_schema(scored)

    out = ROOT / "data" / "external" / "tp53_splice_scored.parquet"
    frozen.to_parquet(out, index=False)
    print(f"\n[phase4a] (a) wrote {out}  shape={frozen.shape}", flush=True)

    par = schema_parity(frozen)
    print("\n=== (b) SCHEMA PARITY vs frozen_matrix_v1 ===")
    print(f"  columns identical (name+order): {par['columns_identical']} "
          f"({par['n_cols_tp53']} vs {par['n_cols_ref']})")
    print(f"  missing_in_tp53: {par['missing_in_tp53']}")
    print(f"  extra_in_tp53:   {par['extra_in_tp53']}")
    print(f"  dtype_mismatches: {par['dtype_mismatches']}")
    frows, fok = feature_parity(frozen)
    print(f"  8-feature parity (active TP53 features): {'ALL MATCH' if fok else 'MISMATCH'}")
    for c, in_both, dr, dt, ok in frows:
        print(f"    {c:14} present={in_both}  frozen={dr:8} tp53={dt:8}  {'ok' if ok else 'MISMATCH'}")

    rows, flags = coverage_report(frozen)
    print("\n=== (c) PER-PREDICTOR COVERAGE on 192 TP53 splice variants ===")
    print(f"  {'predictor':13} {'scored':>8} {'cov':>6}  {'expected':9}  status")
    for short, nn, n, frac, exp, status in rows:
        print(f"  {short:13} {nn:4}/{n:<3}   {frac:5.0%}  {exp:9}  {status}")
    print("\n  FLAGS:", flags or "none (feasible tools match the expected pattern)")


if __name__ == "__main__":
    main()
