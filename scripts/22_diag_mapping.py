"""Diagnostics: (1) what do RAD51C post-mapped expressions look like (why few SNVs);
(2) are unmapped variants concentrated in splice/intronic (headline-bias check)?"""
import json
import re
from collections import Counter

import pandas as pd
import requests

from config import DATA_RAW, DATA_PROCESSED

CACHE = DATA_RAW / "mavedb" / "mapped"
S = requests.Session(); S.trust_env = False

# --- (1) RAD51C expression shapes ---
data = json.loads((CACHE / "urn_mavedb_00000673-0-1.json").read_text(encoding="utf-8"))
cur = [m for m in data if m.get("current") and m.get("postMapped")]
synt = Counter()
op = Counter()
samples = []
GHGVS = re.compile(r"NC_\d+\.\d+:g\.(.+)$")
for m in cur:
    for e in m["postMapped"].get("expressions", []):
        synt[e.get("syntax")] += 1
        if e.get("syntax") == "hgvs.g":
            g = e["value"]
            mm = GHGVS.match(g)
            tail = mm.group(1) if mm else g
            if ">" in tail and re.match(r"^\d+[ACGTN]>[ACGTN]$", tail):
                op["SNV"] += 1
            elif "delins" in tail:
                op["delins"] += 1
            elif "del" in tail:
                op["del"] += 1
            elif "ins" in tail:
                op["ins"] += 1
            elif "dup" in tail:
                op["dup"] += 1
            elif ">" in tail:
                op["MNV/other_sub"] += 1
                if len(samples) < 8:
                    samples.append(g)
            else:
                op["other"] += 1
print("RAD51C postMapped expression syntaxes:", dict(synt))
print("RAD51C hgvs.g operation types:", dict(op))
print("sample non-simple substitutions:", samples)

# --- (2) unmapped-vs-region from the master ---
df = pd.read_csv(DATA_PROCESSED / "functional_scores_master.tsv", sep="\t")
print("\n=== mapped (SNV coords) vs region_class, per gene ===")
for gene, sub in df.groupby("gene"):
    tab = sub.assign(has_snv=sub["is_snv"]).pivot_table(
        index="region_class", columns="has_snv", values="variant_urn",
        aggfunc="count", fill_value=0)
    tab.columns = ["no_coord" if c is False else "snv_coord" for c in tab.columns]
    print(f"\n--- {gene} ---")
    print(tab.to_string())
