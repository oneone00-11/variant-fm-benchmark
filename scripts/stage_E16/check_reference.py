"""E16: check every reference base of the four new genes' variant tables against the
Ensembl release-112 chromosome FASTA.

The FASTA (phase1/data/evidence/refs/, gitignored) is fetched by fetch_reference.py
or by scripts/reproduce_evidence.py --fetch-inputs. It is a different source from
the Ensembl REST sequence slice the mapping checked each base against, so this check
does not lean on the mapping's own inputs. Its sha256 must match the one recorded in
the gene's raw manifest.

Writes phase1/data/evidence/external/newgenes_reference_check.json: an input the
pipeline only reads, listed in scripts/reproduce_inputs.txt and held to a 100% match
by phase1/tests/test_evid_newgenes.py, which also repeats the check wherever the
FASTA is present. Exits non-zero, after writing the record, if any base differs.

Run from the repository root, after the stage has written the variant tables:
    .venv/bin/python scripts/stage_E16/check_reference.py               # E16: newgenes_reference_check.json
    .venv/bin/python scripts/stage_E16/check_reference.py --gene SBDS   # E16b: sbds_reference_check.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

EXT = REPO / "phase1" / "data" / "evidence" / "external"
GENES = ["CTCF", "SFPQ", "RAD51D", "XRCC2"]


def check_gene(gene: str, log: RunLog) -> dict:
    man_path = EXT / f"{gene.lower()}_raw" / "manifest.json"
    man = json.loads(man_path.read_text())
    ref = man["reference_fasta"]
    fa = REPO / ref["file_fa"]
    if sha256_of(fa) != ref["sha256_fa"]:
        raise SystemExit(f"{rel(fa)} differs from the sha256 its manifest records")
    table = EXT / f"{gene.lower()}_variants.parquet"
    v = pd.read_parquet(table)
    for path in (man_path, fa, table):
        log.input(path)
    got = []
    with pysam.FastaFile(str(fa)) as f:
        for chrom, pos in zip(v["chrom"].astype(str), v["pos"].astype(int)):
            got.append(f.fetch(chrom, pos - 1, pos).upper())
    v["fasta_base"] = got
    bad = v[v["fasta_base"] != v["ref"]]
    return {
        "gene": gene, "chrom": ref["chrom"],
        "fasta": {"url": ref["url"], "ensembl_release": ref["ensembl_release"],
                  "sha256_fa": ref["sha256_fa"], "sha256_gz": ref["sha256_gz"]},
        "variants_file": rel(table), "variants_sha256": sha256_of(table),
        "n": int(len(v)), "matches": int(len(v) - len(bad)), "mismatches": int(len(bad)),
        "mismatch_examples": bad[["hgvs_nt", "chrom", "pos", "ref", "fasta_base"]]
                             .head(10).to_dict(orient="records"),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", nargs="+")
    genes = ap.parse_args().gene
    out = EXT / (f"{genes[0].lower()}_reference_check.json" if genes else "newgenes_reference_check.json")
    with RunLog(__file__) as log:
        results = [check_gene(g, log) for g in (genes or GENES)]
        record = {"checked_utc": utc_now(), "rule": "every ref equals the FASTA base "
                  "at its 1-based position; 100% required (E16.3)", "genes": results}
        out.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n")
        log.output(out)
        log.note("summary", {r["gene"]: f"{r['matches']}/{r['n']}" for r in results})
        for r in results:
            print(f"{r['gene']:7s} chr{r['chrom']}: {r['matches']:,}/{r['n']:,} reference "
                  f"bases match the Ensembl release-{r['fasta']['ensembl_release']} FASTA")
        if any(r["mismatches"] for r in results):
            raise SystemExit("reference bases differ; stopped (E16.10)")


if __name__ == "__main__":
    main()
