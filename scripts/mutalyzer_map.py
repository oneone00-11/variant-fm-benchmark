"""Map transcript-c. HGVS to GRCh38 genomic coords via Mutalyzer, using the
chromosomal-reference + RefSeq-transcript-selector form (the one that resolves
intronic offsets). Thread-safe session, on-disk cache, retries.

ENST-based score sets are mapped through their MANE-Select RefSeq equivalent
(c. numbering identical under MANE; validated separately).
"""
import json
import re
import threading
import time

import requests

from config import DATA_RAW

# gene -> (chromosome NC accession GRCh38, RefSeq transcript for c. numbering)
GENE_REF = {
    "BRCA1":  ("NC_000017.11", "NM_007294.3"),
    "BARD1":  ("NC_000002.12", "NM_000465.4"),
    "BRCA2":  ("NC_000013.11", "NM_000059.4"),   # MANE equiv of ENST00000380152.8
    "PALB2":  ("NC_000016.10", "NM_024675.4"),
    "RAD51C": ("NC_000017.11", "NM_058216.3"),   # MANE equiv of ENST00000337432.9
}
NC_CHROM = {
    **{f"NC_0000{n:02d}": str(n) for n in range(1, 23)},
    "NC_000023": "X", "NC_000024": "Y", "NC_012920": "MT",
}

CACHE_PATH = DATA_RAW / "mavedb" / "mutalyzer_cache.json"
_GHGVS = re.compile(r"(NC_\d+)\.\d+:g\.(\d+)([ACGTN]+)>([ACGTN]+)")
_local = threading.local()
_cache_lock = threading.Lock()


def _session():
    if not hasattr(_local, "s"):
        s = requests.Session()
        s.trust_env = False
        _local.s = s
    return _local.s


def load_cache():
    if CACHE_PATH.exists():
        return json.loads(CACHE_PATH.read_text(encoding="utf-8"))
    return {}


def save_cache(cache):
    with _cache_lock:
        CACHE_PATH.write_text(json.dumps(cache), encoding="utf-8")


def ctoken_from_hgvs(hgvs_nt):
    """'ENST00000380152.8:c.7436-10T>A' -> 'c.7436-10T>A' (None if not c.)."""
    if ":c." not in hgvs_nt:
        return None
    return hgvs_nt.split(":", 1)[1]


def build_query(gene, ctoken):
    nc, nm = GENE_REF[gene]
    return f"{nc}({nm}):{ctoken}"


def map_one(gene, ctoken, retries=3):
    """Return (chrom, pos, ref, alt) GRCh38 SNV, or None. Raises only on give-up."""
    q = build_query(gene, ctoken)
    last = None
    for attempt in range(retries):
        try:
            r = _session().get(f"https://mutalyzer.nl/api/normalize/{q}", timeout=90)
            if r.status_code == 200:
                d = r.json()
                eq = (d.get("equivalent_descriptions") or {}).get("g", [])
                for g in eq:
                    m = _GHGVS.match(g["description"])
                    if m:
                        ncacc, pos, ref, alt = m.groups()
                        chrom = NC_CHROM.get(ncacc)
                        if chrom and len(ref) == 1 and len(alt) == 1:
                            return (chrom, int(pos), ref, alt)
                return None  # mapped but not a clean SNV / no g
            elif r.status_code in (429, 500, 502, 503):
                last = f"HTTP {r.status_code}"
                time.sleep(1.5 * (attempt + 1))
                continue
            else:
                return None  # 4xx parse error -> unmappable, don't retry
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
            time.sleep(1.0 * (attempt + 1))
    raise RuntimeError(f"give up {q}: {last}")
