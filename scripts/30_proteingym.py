"""Check ProteinGym DMS coverage for the 6 genes (Milestone 2, subtask 5).

ProteinGym is missense / protein-level only (substitution DMS + clinical sets),
so it CANNOT serve the splice/non-coding line. We just record which genes have
a substitution-DMS assay.
"""
import io
import requests
import pandas as pd

from config import TARGET_GENES, RESULTS_TABLES, DATA_RAW

S = requests.Session(); S.trust_env = False
UNIPROT = {
    "BRCA1": "P38398", "BRCA2": "P51587", "BARD1": "Q99728",
    "RAD51C": "O43502", "PALB2": "Q86YC2", "ATM": "Q13315",
}
URLS = [
    "https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/reference_files/DMS_substitutions.csv",
    "https://marks.hms.harvard.edu/proteingym/reference_files/DMS_substitutions.csv",
]

text = None
for url in URLS:
    try:
        r = S.get(url, timeout=60)
        if r.status_code == 200 and len(r.content) > 1000:
            text = r.text
            print(f"fetched reference from {url} ({len(text)} bytes)")
            break
        print(f"[{r.status_code}] {url}")
    except Exception as e:
        print(f"[ERR] {url}: {type(e).__name__}: {str(e)[:80]}")

if text is None:
    print("\nCould not fetch ProteinGym reference directly (likely network).")
    print("Fallback: will rely on WebFetch / manual entry.")
    raise SystemExit(0)

(DATA_RAW / "proteingym").mkdir(parents=True, exist_ok=True)
(DATA_RAW / "proteingym" / "DMS_substitutions.csv").write_text(text, encoding="utf-8")

df = pd.read_csv(io.StringIO(text))
print("\nreference columns:", list(df.columns)[:12])
# Identify the gene/uniprot column
cols = {c.lower(): c for c in df.columns}
id_col = cols.get("dms_id") or cols.get("dms id") or df.columns[0]
uni_col = next((df.columns[i] for i, c in enumerate(df.columns)
                if "uniprot" in c.lower()), None)

hits = []
for gene, up in UNIPROT.items():
    mask = df[id_col].astype(str).str.contains(gene, case=False, na=False)
    if uni_col:
        mask = mask | df[uni_col].astype(str).str.contains(up, na=False)
    sub = df[mask]
    for _, row in sub.iterrows():
        hits.append({
            "gene": gene, "uniprot": up,
            "DMS_id": row[id_col],
            "n_mutants": row.get(cols.get("dms_number_single_mutants"))
                         or row.get("DMS_number_single_mutants", ""),
            "assay": "protein-level DMS (missense)",
        })

res = pd.DataFrame(hits)
res.to_csv(RESULTS_TABLES / "proteingym_coverage.tsv", sep="\t", index=False)
print("\n=== ProteinGym DMS coverage (PROTEIN-LEVEL / MISSENSE ONLY) ===")
if res.empty:
    print("  no matches")
else:
    print(res.to_string(index=False))
for gene in TARGET_GENES:
    n = (res.gene == gene).sum() if not res.empty else 0
    print(f"  {gene}: {n} ProteinGym DMS assay(s)")
