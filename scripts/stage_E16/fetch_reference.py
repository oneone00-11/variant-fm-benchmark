"""E16: transcript models, sequence slices and reference chromosomes for the four
new genes.

The coordinates are made the way the DDX3X deposit's were (evid_external
prepare_ddx3x), because MaveDB holds no GRCh38 postMapped record for any intron-side
SNV of these four score sets (scripts/stage_E16/run_log/probe_raw_*.json; route
chosen by the author on 2026-10-20). For each gene this fetches, through the atlas
mapping module, which caches and registers every file itself in
phase1/data/evidence/reference_cache/:
  - the MANE Select record (the cached MANE v1.5 summary; no download);
  - the Mutalyzer transcript model of the deposit's RefSeq transcript;
  - the Ensembl lookup of its MANE Select Ensembl partner (the genomic exon anchor);
  - an Ensembl sequence slice over the gene with the module's REGION_PAD each side.

It also fetches the Ensembl release-112 single-chromosome FASTA of each gene's
chromosome, the source refs/download_refs.sh and refs/chrX.fa come from, into
phase1/data/evidence/refs/ (gitignored), unpacks it, and records its URL, release
and sha256 (packed and unpacked) in the gene's raw manifest, where
scripts/reproduce_evidence.py --fetch-inputs reads them. The chromosome is used only
to check every reference base independently of the Ensembl REST slice
(scripts/stage_E16/check_reference.py and the tests).

Run from the repository root, after fetch_mavedb.py:
    .venv/bin/python scripts/stage_E16/fetch_reference.py               # E16's four genes
    .venv/bin/python scripts/stage_E16/fetch_reference.py --gene SBDS   # E16b
A chromosome already fetched for another gene is reused, not downloaded again.
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of  # noqa: E402

sys.path.insert(0, str(REPO / "phase1"))
from src.evid_common import resolve_atlas  # noqa: E402

ATLAS, _missing = resolve_atlas()
sys.path.insert(0, str(ATLAS / "src"))
from atlas import mapping as M  # noqa: E402

EVID = REPO / "phase1" / "data" / "evidence"
EXT = EVID / "external"
REFS = EVID / "refs"
M.REF_CACHE = EVID / "reference_cache"     # the tracked cache, as evid_external uses
GENES = ["CTCF", "SFPQ", "RAD51D", "XRCC2"]
ENSEMBL_RELEASE = 112
FASTA_URL = ("https://ftp.ensembl.org/pub/release-{rel}/fasta/homo_sapiens/dna/"
             "Homo_sapiens.GRCh38.dna.chromosome.{chrom}.fa.gz")


def _entry_point():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def chromosome(chrom: str, log: RunLog) -> dict:
    url = FASTA_URL.format(rel=ENSEMBL_RELEASE, chrom=chrom)
    gz = REFS / url.rsplit("/", 1)[1]
    fa = REFS / f"chr{chrom}.fa"
    if not gz.exists():
        print(f"downloading {url}", flush=True)
        _entry_point()._download(url, gz)
    if not fa.exists():
        with gzip.open(gz, "rb") as src, open(fa, "wb") as dst:   # gzip checks its CRC
            shutil.copyfileobj(src, dst, 1 << 20)
    import pysam
    with pysam.FastaFile(str(fa)) as f:                           # writes the .fai
        names = list(f.references)
    if chrom not in names:
        raise SystemExit(f"{rel(fa)} holds {names}, not chromosome {chrom}")
    log.output(gz)
    log.output(fa)
    return {"chrom": chrom, "ensembl_release": ENSEMBL_RELEASE, "url": url,
            "file_gz": rel(gz), "sha256_gz": sha256_of(gz), "bytes_gz": gz.stat().st_size,
            "file_fa": rel(fa), "sha256_fa": sha256_of(fa), "bytes_fa": fa.stat().st_size,
            "tracked": False, "note": "phase1/data/evidence/refs/ is gitignored; "
                                      "reproduce_evidence.py --fetch-inputs fetches it"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", nargs="+", default=GENES)
    genes = ap.parse_args().gene
    with RunLog(__file__) as log:
        log.note("atlas", rel(ATLAS))
        mane = M.load_mane_records(genes)
        done = {}
        for gene in genes:
            raw = EXT / f"{gene.lower()}_raw"
            man_path = raw / "manifest.json"
            log.input(man_path)
            log.input(raw / "scoreset_meta.json")
            man = json.loads(man_path.read_text())
            meta = json.loads((raw / "scoreset_meta.json").read_text())
            target = meta["targetGenes"][0]["targetAccession"]["accession"]
            rec = mane[gene]
            tmap, prov = M.build_transcript_map(target, rec)
            lo = min(e.g_start for e in tmap.exons) - M.REGION_PAD
            hi = max(e.g_end for e in tmap.exons) + M.REGION_PAD
            M.fetch_ensembl_region(rec.chr_name, lo, hi)
            cache = {
                "mane_summary": "mane/MANE.GRCh38.v1.5.summary.txt.gz",
                "mutalyzer_model": f"mutalyzer/reference_model_{prov['mutalyzer_model_accession']}.json",
                "ensembl_lookup": f"ensembl/lookup_{prov['ensembl_anchor']}.json",
                "region": f"ensembl/seq_GRCh38_chr{rec.chr_name}_{lo}_{hi}.fa",
            }
            for key, path in cache.items():
                full = M.REF_CACHE / path
                if not full.is_file():
                    raise SystemExit(f"{gene}: {key} not cached at {rel(full)}")
                log.output(full)
            man["transcript_model"] = {
                "route": ("atlas mapping module, as for DDX3X: Mutalyzer transcript model "
                          "anchored to GRCh38 by the Ensembl lookup of the MANE Select "
                          "partner; not MaveDB postMapped"),
                "mane_select": {"symbol": rec.symbol, "refseq_nuc": rec.refseq_nuc,
                                "ensembl_nuc": rec.ensembl_nuc,
                                "grch38_chr_accession": rec.grch38_chr_accession,
                                "chr_name": rec.chr_name, "strand": rec.strand},
                "deposit_target_accession": target,
                "deposit_target_is_mane_select_refseq": target == rec.refseq_nuc,
                "provenance": prov,
                "region": {"chrom": rec.chr_name, "start": lo, "end": hi,
                           "pad": M.REGION_PAD},
                "cache_files": {k: {"path": rel(M.REF_CACHE / p),
                                    "sha256": sha256_of(M.REF_CACHE / p)}
                                for k, p in cache.items()},
            }
            man["reference_fasta"] = chromosome(rec.chr_name, log)
            man_path.write_text(json.dumps(man, indent=1, ensure_ascii=False) + "\n")
            log.output(man_path)
            done[gene] = {"target": target, "mane_refseq": rec.refseq_nuc,
                          "mane_ensembl": rec.ensembl_nuc, "chrom": rec.chr_name,
                          "strand": rec.strand, "n_exons": prov["n_exons"],
                          "exon_anchor_adjustments": prov["exon_anchor_adjustments"]}
            print(f"{gene}: {target} (MANE {rec.refseq_nuc}, {rec.ensembl_nuc}), chr"
                  f"{rec.chr_name} {rec.strand}, {prov['n_exons']} exons, region {lo}-{hi}")
        log.output(M.REF_CACHE / "MANIFEST.json")
        log.note("genes", done)


if __name__ == "__main__":
    main()
