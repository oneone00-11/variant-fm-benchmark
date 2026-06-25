# Milestone 9 — Reproducible repository for submission

Turned the working project into a clean, public-ready, one-command-reproducible
repository (target journal: *Bioinformatics*). **Pure organization + documentation:
no model was re-scored and no result number changed** (verified — see §5).

## 1. What was added / changed

| Item | Action |
|---|---|
| `README.md` | Full rewrite: summary, headline result, layout, install, **one-command reproduce**, data-acquisition table, model versions, gold-standard URN/license table, credentials note |
| `scripts/reproduce_main.py` | **New.** One command runs Phase-3 (100→103) from `score_matrix_final.tsv`, regenerates all tables + 3 figures, prints the headline splice ranking |
| `scripts/README.md` | **New.** Maps the 50 numbered scripts to milestones (layout kept flat — see §3) |
| `requirements.txt` | Pinned to the tested versions (pandas 2.3.3 … matplotlib 3.9.4); re-scoring deps marked optional |
| `environment.yml` | **New.** Conda alternative |
| `LICENSE` | **New.** MIT for code + explicit third-party data-license note |
| `.gitignore` | Expanded: `refs/` large binaries, root scratch copies, all `*cache*` files |
| `.gitattributes` | **New.** `* text=auto eol=lf` — removes the Windows→Mac CRLF churn |
| `docs/milestone9_repo.md` | This brief |
| scoring scripts | Hardcoded paths → repo-relative (see §4) |

## 2. Credential security check — PASS (no leak)

Highest-priority guardrail. Swept **working tree + entire git history**:
- AlphaGenome API key is read **only** from `os.environ["ALPHAGENOME_API_KEY"]`
  (`scripts/91_score_alphagenome.py`; the script aborts if unset). Never hard-coded,
  printed, cached, or logged.
- No API keys / tokens / HF tokens / bearer secrets anywhere in `*.py`, `*.sh`,
  `*.md`, `*.json`, `*.log` (all grep hits were false positives — HGVS "ctoken"
  helpers, the word "tokens" in comments).
- `git log -p --all` secret scan: **clean** — nothing sensitive was ever committed.
- The author's own `refs/leak_check.sh` corroborates this.
- Belt-and-suspenders: cache files (`*cache*.json/tsv`, `ag_cache.tsv`) are
  git-ignored even though verified key-free.

## 3. Repository-structure decisions

- **`scripts/` kept flat (not foldered by milestone).** 34 scripts import shared
  modules via `from config import …` / `from phase3_lib import …`; moving them into
  subdirectories would break those imports and risk the results. The numeric prefixes
  already encode phase order, and `scripts/README.md` documents the milestone grouping.
- **Git history preserved.** It is clean (no secrets; largest tracked blob 13 MB) and
  the milestone commits are valuable provenance. Milestone 9 is a new commit on `main`.
- **Personal / internal notes left untracked** (not committed): `HANDOFF.md`,
  `变异解读FM基准测试_项目计划.md`, `大致结构.txt`. They stay local; commit them later if desired.
- **Redundant root copies removed from tracking:** `score_matrix_final_with_nt.tsv`
  (byte-identical to `data/processed/score_matrix_final.tsv`) and `nt_cache2.tsv` are
  git-ignored. The canonical final matrix is `data/processed/score_matrix_final.tsv`.

## 4. Path portability (clone-and-run)

| Script | Before | After |
|---|---|---|
| `90_score_gpn.py` | `/mnt/d/variant-fm-benchmark/…` | `Path(__file__).resolve().parents[1] / …` |
| `91_score_alphagenome.py` | `/mnt/d/…` | repo-relative |
| `95_check_nt.py` | `C:\Users\张宁一\Desktop\…` (personal name) | repo-relative `data/processed/score_matrix_final.tsv` |
| `74_alphamissense.py` | docstring "to D:" | "to data/raw/scores" |

The Phase-3 analysis scripts (`100`–`103`, `config.py`) were already repo-relative.
`refs/*.sh` are WSL/conda **provenance** scripts (illustrative `/mnt/d` paths, not run on clone).

## 5. Large files NOT committed (with acquisition route)

Git-ignored; README "Data acquisition" has full details:

| Excluded | Size | Get it from |
|---|--:|---|
| `refs/grch38_subset.gtf.db` / `.fa` | 799 MB / 740 MB | rebuilt from Ensembl release-112 |
| `data/raw/scores/AlphaMissense_hg38.tsv.gz` | 642 MB | Zenodo |
| `refs/Homo_sapiens.GRCh38.*` FASTA/GTF | ~55–72 MB each | Ensembl release-112 |
| `data/raw/clinvar_GRCh38.vcf.gz` | 192 MB | NCBI ClinVar FTP (GRCh38, 2026-06-15) |
| `data/raw/mavedb/…` score sets | up to 70 MB | MaveDB API (URNs in README) |
| GPN-MSA scores | (81 GB remote) | HuggingFace `songlab/gpn-msa-hg38-scores`; only the 2.65 MB `.tbi` is local |

The **final evaluation matrix ships in the repo** (`data/processed/score_matrix_final.tsv`,
8 MB), so reproduction needs none of the above.

## 6. One-command reproduction — VERIFIED

```bash
pip install -r requirements.txt
python scripts/reproduce_main.py
```
Ran clean (exit 0) on this machine. Regenerated all 19 result tables + 3 figures.
`git diff --ignore-cr-at-eol` vs the committed outputs: **all tables byte-identical**;
`score_matrix_final.tsv` unchanged. Determinism via `np.random.default_rng(20260619)`.

## 7. Remaining manual steps for the author (out of scope here)

- Set the GitHub remote and push (the repo is committed locally on `main`).
- Confirm the LICENSE copyright name (currently "Zhang Ningyi" from git config) and
  license choice (MIT) — swap if a different one is preferred.
- Decide whether to also publish the data (e.g. Zenodo DOI) alongside the code.
