"""
phase4a_tp53_external.py -- produce the 3 WSL/API predictor caches for TP53.

These are the version-pinned predictors that did NOT run on the Mac session and
must run in the same env as the 7-gene run (macOS lacks the SpliceAI/Pangolin
CLIs and the AlphaGenome client, and the AG key is env-only). Each subcommand
reuses the EXACT seven-gene code and writes a cache TSV that
`phase4a_score_tp53.py` reads:  chrom<TAB>pos<TAB>ref<TAB>alt<TAB><col>

FLOW (in the WSL checkout, after copying data/external/tp53_variants.vcf there):

  # 1) SpliceAI (conda spliceai) -- exactly as refs/run_spliceai_full.sh, but on the TP53 VCF:
  conda activate spliceai
  spliceai -I data/external/tp53_variants.vcf -O data/external/tp53_spliceai_out.vcf \
           -R refs/grch38_subset.fa -A grch38 -D 50
  python scripts/phase4a_tp53_external.py parse-spliceai data/external/tp53_spliceai_out.vcf

  # 2) Pangolin (conda pangolin) -- exactly as refs/run_pangolin_full.sh, CSV mode:
  conda activate pangolin
  printf 'CHROM,POS,REF,ALT\n' > data/external/tp53_pangolin_in.csv
  grep -v '^#' data/external/tp53_variants.vcf | awk 'BEGIN{OFS=","}{print $1,$2,$4,$5}' \
       >> data/external/tp53_pangolin_in.csv
  pangolin -d 50 data/external/tp53_pangolin_in.csv refs/grch38_subset.fa \
           refs/grch38_subset.gtf.db data/external/tp53_pangolin_out
  python scripts/phase4a_tp53_external.py parse-pangolin data/external/tp53_pangolin_out.csv

  # 3) AlphaGenome (conda alphagenome, client 0.6.1, key in env) -- 91_score_alphagenome logic:
  conda activate alphagenome
  pip install 'alphagenome==0.6.1'          # PIN to the 7-gene client (drift-verified)
  export ALPHAGENOME_API_KEY=...            # env only; never printed/cached
  python scripts/phase4a_tp53_external.py alphagenome

Then copy the 3 tp53_{spliceai,pangolin,alphagenome}.tsv back beside the Mac's
data/external/ and re-run:  python scripts/phase4a_score_tp53.py
(the assembler picks them up automatically and fills the 9-feature matrix).
"""
from __future__ import annotations
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
VCF = ROOT / "data" / "external" / "tp53_variants.vcf"
OUT_SPLICEAI = ROOT / "data" / "external" / "tp53_spliceai.tsv"
OUT_PANGOLIN = ROOT / "data" / "external" / "tp53_pangolin.tsv"
OUT_ALPHAGENOME = ROOT / "data" / "external" / "tp53_alphagenome.tsv"


def _load_82():
    import importlib.util
    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location("bench_82", SCRIPTS / "82_assemble_v2.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _variants_from_vcf():
    rows = []
    with open(VCF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            rows.append((c[0], int(c[1]), c[3], c[4]))
    return rows


def parse_spliceai(out_vcf):
    m = _load_82()
    m.SPLICEAI_VCF = out_vcf                       # reuse the exact seven-gene parser
    df = m.parse_spliceai()[["chrom", "pos", "ref", "alt", "spliceai_ds"]]
    df.to_csv(OUT_SPLICEAI, sep="\t", index=False)
    print(f"[external] wrote {OUT_SPLICEAI}: {df.spliceai_ds.notna().sum()}/{len(df)} scored")


def parse_pangolin(out_csv):
    m = _load_82()
    m.PANGOLIN_CSV = out_csv                       # reuse the exact seven-gene parser
    df = m.parse_pangolin()[["chrom", "pos", "ref", "alt", "pangolin_score"]]
    df.to_csv(OUT_PANGOLIN, sep="\t", index=False)
    print(f"[external] wrote {OUT_PANGOLIN}: {df.pangolin_score.notna().sum()}/{len(df)} scored")


def alphagenome():
    """score_one COPIED VERBATIM from 91_score_alphagenome.py (client pinned 0.6.1)."""
    import os, time
    import alphagenome
    if getattr(alphagenome, "__version__", "?") != "0.6.1":
        sys.exit(f"alphagenome=={alphagenome.__version__}; pin the 7-gene client: pip install 'alphagenome==0.6.1'")
    if not os.environ.get("ALPHAGENOME_API_KEY"):
        sys.exit("ALPHAGENOME_API_KEY not set (env only).")
    from alphagenome.models import dna_client, variant_scorers
    from alphagenome.data import genome
    model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
    recs = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    SPLICE = [recs[n] for n in recs if "SPLICE" in str(n).upper()]

    def score_one(v):
        chrom, pos, ref, alt = v
        var = genome.Variant(chromosome="chr" + str(chrom), position=pos,
                             reference_bases=ref, alternate_bases=alt)
        iv = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        for attempt in range(5):
            try:
                out = model.score_variant(interval=iv, variant=var, variant_scorers=SPLICE)
                df = variant_scorers.tidy_scores(out)
                if df is None or len(df) == 0 or "raw_score" not in df.columns:
                    return ""
                return f"{df['raw_score'].abs().max():.6f}"
            except Exception:
                time.sleep(2 * (attempt + 1))
        return None

    rows = []
    for v in _variants_from_vcf():
        val = score_one(v)
        rows.append({"chrom": v[0], "pos": v[1], "ref": v[2], "alt": v[3],
                     "alphagenome_splice": (val if val not in (None, "") else "")})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_ALPHAGENOME, sep="\t", index=False)
    print(f"[external] wrote {OUT_ALPHAGENOME}: "
          f"{(df.alphagenome_splice.astype(str) != '').sum()}/{len(df)} scored")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "parse-spliceai":
        parse_spliceai(sys.argv[2])
    elif cmd == "parse-pangolin":
        parse_pangolin(sys.argv[2])
    elif cmd == "alphagenome":
        alphagenome()
    else:
        sys.exit("usage: phase4a_tp53_external.py {parse-spliceai <out.vcf> | "
                 "parse-pangolin <out.csv> | alphagenome}")
