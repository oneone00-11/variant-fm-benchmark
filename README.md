# Variant-FM-Benchmark

Benchmarking language models for clinical interpretation of cancer gene variants,
using **functional gold standards** (SGE / MAVE) as the primary evaluation, with a
focus on **splice / non-coding variants** that protein-based methods cannot score.

Target genes (DNA-repair / hereditary breast-ovarian cancer):
**BRCA1, BRCA2, BARD1, RAD51C, PALB2, ATM**

See `..\变异解读FM基准测试_项目计划.md` (Desktop) for the full project plan.

## Layout
```
scripts/        pipeline scripts, numbered by phase
  config.py                 shared paths, gene list, ClinVar classification maps
  01_download_clinvar.py    fetch ClinVar GRCh38 VCF
  02_clinvar_inventory.py   Milestone 1: per-gene P/B/VUS counts + strata
data/
  raw/          downloaded sources (git-ignored)
  interim/      intermediate (git-ignored)
  processed/    cleaned tables (the variant master table lives here)
results/
  tables/       output TSVs
  figures/      output figures
docs/           notes, data provenance, decisions
notebooks/      exploratory analysis
```

## Environment
- Python 3.9, CPU-only (no GPU). Strategy: rely on **precomputed scores**
  (AlphaMissense, GPN-MSA, SpliceAI, CADD, REVEL) rather than local LM inference.
- Windows note: PyPI is reached via the Tsinghua mirror; the local system proxy
  is bypassed (`NO_PROXY=*` for pip, `trust_env=False` in download scripts).

### Setup
```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
$env:NO_PROXY="*"
pip install -r requirements.txt
```

### Run Milestone 1 (ClinVar inventory)
```powershell
python scripts\01_download_clinvar.py
python scripts\02_clinvar_inventory.py
```

## Status
- [x] Phase 0: repo skeleton + environment
- [ ] Milestone 1: ClinVar inventory (in progress)
- [ ] MaveDB / functional gold-standard availability check
- [ ] gnomAD frequencies
- [ ] Model score collection
- [ ] Evaluation engine
