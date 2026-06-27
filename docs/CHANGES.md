# Milestone 10 — Manuscript content revision (CHANGES)

Content-level revision of the draft. **No result number was changed, no figure/table
datum altered, no citation/version fabricated.** Format/layout untouched (the v2 docx
is produced by the same generator with text edits only).

- **Input:** `Variant-FM-Benchmark_manuscript_draft.docx` (Desktop) — preserved, unchanged.
- **Output:** `Variant-FM-Benchmark_manuscript_v2.docx` (Desktop).
- **Generator:** `scratchpad/build_manuscript_v2.py` (copy of `build_manuscript.py` with the edits below).
- **Companion:** `docs/OPEN_ITEMS.md` (author to-verify / to-decide list).

Every version number / parameter below was taken from the cited milestone report,
scoring script, or `results/tables/` file — not from memory.

## Part 1 — Numerical & internal consistency
| # | Change | Where | Source / basis |
|---|---|---|---|
| 1 | "eleven predictors" → **"ten predictors"** (4 occurrences: Abstract Results, Intro §1, Methods 2.3, Discussion resource para) | Abstract, Intro, 2.3, 4 | Actual scored predictors = 10 (SpliceAI, Pangolin, AlphaGenome, GPN-MSA, NT, AlphaMissense, CADD, phyloP, phastCons, gnomAD AF); Table 2 has 10 rows. gnomAD global+popmax = two summaries of **one** baseline — stated explicitly in 2.3. |
| 2 | Gene order standardized to **BRCA1, BRCA2, BARD1, PALB2, RAD51C, VHL, BAP1** | Intro §1 (was BRCA1,BRCA2,PALB2,RAD51C,BARD1,…) | Matches Abstract, Methods 2.1, Table 1, `phase3_lib.GENES`. |

## Part 2 — Method completeness & reproducibility (versions/params)
| # | Change | Source |
|---|---|---|
| 3 | Added exact versions/params for every tool in 2.3: SpliceAI **v1.3.1**, distance **−D 50**, max of 4 deltas; Pangolin (**tkzeng/Pangolin**), **−d 50**, default mask, max(gain,\|loss\|); reference **Ensembl GRCh38 release-112**; AlphaGenome client **v0.6.1**, scorers **SPLICE_SITES + SPLICE_SITE_USAGE + SPLICE_JUNCTIONS**, **1-Mb (SEQUENCE_LENGTH_1MB)**, max\|raw_score\|; GPN-MSA source **songlab/gpn-msa-hg38-scores**; NT **InstaDeep v2-500m-multi-species**, masked ALT/REF log-likelihood ratio; AlphaMissense **hg38, Zenodo rec. 8208688**; CADD **GRCh38-v1.7**; phyloP/phastCons **UCSC 100-way**; gnomAD **v4**. | milestone6 (SpliceAI/Pangolin/ref), milestone7 (AlphaGenome/GPN-MSA), milestone5 (CADD/AlphaMissense/conservation/gnomAD). **NT checkpoint is author-stated** (score_nt2.py not in repo) → see OPEN_ITEMS. |
| 4 | 2.2: added "≈2% of functional variants could not be mapped after Mutalyzer retries and were excluded; codon-level MNV/delins out of scope." | milestone2 ("142 variants ~2% gave up"), milestone3 (9,220 non-SNV/unmapped held out). |
| 5 | 2.1: explicit honest disclosure that the splice/intron sub-class uses the **functional intron offset**, so **exon-side splice-region variants (last 1–3 exonic nt) are classified as missense/synonymous, not splice** → the splice subset is the intron-side splice region (conservative). Cross-referenced in Limitations. | milestone3 "Region rule" + documented caveat. |

## Part 3 — Added caveats / limitations
| # | Change | Source / basis |
|---|---|---|
| 6 | Limitations: added that "DNA language model" conclusions rest on **only two** models (GPN-MSA, NT); generalize cautiously (folded together with the existing Evo2/ESM-not-scored note). | Scope of this study (2.3); milestone7/8b (NT/Evo2/ESM status). |
| 7 | 3.2: added gnomAD-baseline caveat — ρ computed only on **526/1,781 (≈30%)** splice variants present in gnomAD v4; low ρ partly reflects a frequency-biased, range-restricted subset; baseline only. | `full_model_coverage.tsv` (gnomAD splice scored=526), `spearman_meta.tsv` (gnomad splice total_n=526). |
| 8 | Softened over-claim: "the largest assembled against a functional standard for these genes" → "**sizeable … and adequately powered for the pooled continuous analysis by a priori power calculations**". | milestone3 power analysis (splice pooled sufficient). |

## Part 4 — Auxiliary metrics & statistical transparency
| # | Change | Source |
|---|---|---|
| 10 | 3.4: **AUPRC now reported** (splice, excl. BRCA1: Pangolin 0.999, SpliceAI 0.998, AlphaGenome 0.998, GPN-MSA 0.983, NT 0.910) — resolves the 2.4 promise of AUROC/AUPRC. | `results/tables/auroc_clinvar.tsv` (auprc column; 101_auroc.py computes `average_precision_score`). |
| 11 | 2.4: stated **no multiple-comparison correction** was applied; pairwise CIs interpreted descriptively. | `scripts/102_clinical_pairwise.py` — raw percentile bootstrap CIs, no Bonferroni/FDR. |
| 12 | 2.4: noted pairwise test uses the **unweighted mean of per-gene ρ** whereas the headline uses **inverse-variance-weighted DerSimonian–Laird**, so point estimates differ marginally. | `102_clinical_pairwise.py` (`mean_gene_rho` = unweighted) vs `phase3_lib.dl_meta` (DL). |

## Part 5 — References (only verifiable additions; rest → OPEN_ITEMS)
| # | Change | Note |
|---|---|---|
| 13 | Added **phyloP (Pollard et al., 2010)** and **phastCons (Siepel et al., 2005)** references; CADD (Rentzsch) and gnomAD already cited. | Marked "[bibliographic details to verify]" — do not treat page numbers as confirmed. |
| 14 | Flagged the **gnomAD version mismatch** in the reference list: text uses v4 but the listed Karczewski 2020 is the v2 paper. | Resolution deferred to OPEN_ITEMS (cite gnomAD v4 paper or clarify). Not auto-edited to avoid fabricating a citation. |
| 15 | All "to verify / Preprint / details" reference entries consolidated into OPEN_ITEMS. | See OPEN_ITEMS §A. |

## Part 6 — Wording & clarity
| # | Change |
|---|---|
| 16 | 3.1 & Abstract: "all four DNA/splice models … and the Nucleotide Transformer …" → "**all five DNA-reading models (SpliceAI, Pangolin, AlphaGenome, GPN-MSA, NT)** scored every splice variant" — removes the four-plus-one ambiguity. |
| 17 | AlphaGenome terminology unified to "**DNA sequence model**" (Intro previously "multi-modal sequence model"; Table 2 already "DNA seq. model"). |
| 18 | Abstract/Intro cancer scope broadened from "breast, ovarian and renal" to "**breast, ovarian, renal and other**" (BAP1 → mesothelioma/uveal melanoma noted in Intro). |

## Part 7 — Placeholders
| # | Change |
|---|---|
| 19 | No placeholders invented. All ([author], [affiliation], [contact], GitHub [user], Funding, Associate Editor, Table S1) left as-is and listed in OPEN_ITEMS §B. |

## Part 8 — Design decision (manuscript side)
| # | Change |
|---|---|
| 9/20 | 2.1: added a neutral one-line rationale that the primary set is the ClinVar∩functional intersection, with the full-functional-set extension flagged to the Discussion. **No re-analysis performed or results invented** — the decision to actually re-run on the full ~34,582-SNV set is left to the author (OPEN_ITEMS §C). |

## Not changed (guardrails)
- All ρ, AUROC, AUPRC, Δρ, CI, I², N values — unchanged (verified identical to `results/tables/`).
- Figures and tables — data untouched.
- Layout, fonts, styles — identical generator settings.

---

# Milestone 11 — NT reproducibility + Supplementary Table S1

Resolves OPEN_ITEMS C-3 (NT scoring not reproducible) and B (Table S1 missing).
**No model was re-run, no score changed, no PMID/version/NT-window fabricated.**

- **Manuscript:** produced `Variant-FM-Benchmark_manuscript_v3.docx` (v2 and the
  original draft preserved). Only edit: Methods 2.3 NT sentence now states the full
  checkpoint (InstaDeepAI/nucleotide-transformer-v2-500m-multi-species), the precise
  scoring (mask the non-overlapping 6-mer token at the variant; reference-minus-
  alternate masked-token log-probability; higher = pathogenic), and points to the
  archived script. **No context-window number was inserted** (unrecoverable → OPEN_ITEMS).
  Generator: `scratchpad/build_manuscript_v3.py`.

| # | Change | Source / basis |
|---|---|---|
| C-3a | **Searched** the machine + repo for the NT script (`score_nt*.py`, `*nucleotide*`, notebooks) and cache. Result: **no script anywhere**; only the cache `nt_cache2.tsv` (root, git-ignored) survives. → rebuild case. | filesystem search |
| C-3b | **Rebuilt `scripts/92_score_nt.py`** as a clearly-labelled archival reconstruction of the NT scoring flow (read variants → GRCh38 reference window → mask the variant's 6-mer token → reference-minus-alternate masked log-prob → per-variant score), repo-relative paths, resumable, regenerates `nt_cache2.tsv`. Header states it is a flow archive (original ran on a single cloud GPU), **not a re-run**, with **no fabricated scores**. | docs/milestone8b; orientation verified from `score_matrix_final.tsv` (splice P/LP NT mean ≈+2.9 vs B/LB ≈+0.5); cache format read from `nt_cache2.tsv` (key=chrom:pos:ref:alt, nt_score) |
| C-3c | NT orientation confirmed **higher = pathogenic** (matches `phase3_lib.ORIENT["nucleotide_transformer"]=+1`); documented in the script. | `score_matrix_final.tsv`, `phase3_lib.py` |
| 5 | NT **context window**: not inserted into the manuscript. The original run's window is unrecorded; the rebuilt script uses a documented default `WINDOW_BP=6000` (~1000 tokens, << NT v2's 2048-token max) explicitly labelled "to be confirmed", and the manuscript points to the script rather than asserting a number. | header of `92_score_nt.py`; flagged in OPEN_ITEMS C-3 |
| 6 | **README.md** and **scripts/README.md**: NT moved into the reproducible-models list (needs GPU), with how-to-rerun and the cache-regeneration note; Status line updated. | this milestone |
| S1 | **Supplementary Table S1** generated: `supplementary/Table_S1.tsv` + `Table_S1.docx`, all 7 genes — gene, MaveDB URN, assay (SGE), reference transcript, MANE RefSeq used for mapping, assembly (GRCh38), reference (PMID/source), licence. | README gold-standard table + docs/milestone2 + milestone4 |
| S1-note | PMIDs in Table S1 carry `[verify]`; **BARD1 and PALB2 PMIDs are "not recorded"** in project materials (not invented). | docs/milestone2 (no PMID for BARD1/PALB2) |

## Not changed (Milestone 11 guardrails)
- No model re-scored; `nt_cache2.tsv`, `score_matrix_final.tsv`, all result tables/figures unchanged.
- No PMID, citation, version, or NT context-window value fabricated.
- Cache / large files not committed (`nt_cache2.tsv` stays git-ignored; the script regenerates it).
