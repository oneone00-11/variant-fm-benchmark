# Milestone 4 — Candidate panel-expansion genes (SGE gold-standard check)

Reuse of the Milestone-2 pipeline on candidate tumour-suppressor genes, same
caliber and gates. **No model scoring; nothing merged into analysis_ready yet** —
awaiting your include decision.

Candidates checked (in order): **VHL, BAP1, MSH2**; fallbacks TP53, PTEN, MLH1, NF1.

## 1. SGE availability (gate ①)

| Gene | SGE? | Primary set | N var | Transcript | License | PMID/source |
|---|---|---|--:|---|---|---|
| **BAP1** | ✅ SGE | urn:mavedb:00000662-0-1 | 18108 | ENST00000460680.6 | CC BY 4.0 | 38969833 (Waters 2024) — "CDS, splice-site intron, 3′UTR" |
| **VHL** | ✅ SGE | urn:mavedb:00000675-a-1 | 2268 | ENST00000256474.3 | CC BY 4.0 | 38969834 (Buckley 2023) |
| MSH2 | ❌ | urn:mavedb:00000050-a-1 | 17746 | — | CC0 | 33357406 — missense LOF screen (HAP1), **not SGE, no splice/intron** |
| TP53 | ❌ | (many) | — | — | CC0 | transactivation/DN assays, not SGE |
| PTEN | ❌ | urn:mavedb:00000102 | — | — | CC0 | VAMP-seq abundance, not SGE |
| MLH1 | ❌ | urn:mavedb:00001218 | — | — | CC0 | DHFR-PCA abundance (C-terminal, missense), not SGE |
| NF1 | ❌ | — | 0 | — | — | no on-target score set |

Only **VHL and BAP1** have SGE; both explicitly cover splice/intron. Full list:
`results/tables/candidate_mavedb_catalog.tsv`.

## 2–4. Mapping, region coverage, ClinVar overlap (gates ②–④)

Coordinate mapping identical to M2 (post-mapped VRS + Mutalyzer NC(NM) for
splice/intronic; ENST sets via MANE-Select RefSeq NM_000551.4 / NM_004656.4).

`results/tables/candidate_gene_coverage.tsv`:

| Gene | region | functional SNV | ClinVar ≥1★ overlap | of which VUS |
|---|---|--:|--:|--:|
| VHL | splice_core | 24 | 19 | 0 |
| VHL | splice_region | 69 | 38 | 8 |
| VHL | intronic | 481 | 137 | 87 |
| VHL | coding | 1585 | 880 | 317 |
| VHL | utr | 109 | 11 | 2 |
| BAP1 | splice_core | 192 | 77 | 2 |
| BAP1 | splice_region | 576 | 191 | 65 |
| BAP1 | intronic | 1692 | 361 | 14 |
| BAP1 | coding | 6501 | 2363 | 1168 |
| BAP1 | utr | 480 | 26 | 15 |

ClinVar SNVs for VHL/BAP1 were extracted from the full GRCh38 VCF (release
2026-06-15): VHL 225 P/LP, 539 B/LB, 770 VUS; BAP1 170 P/LP, 1296 B/LB, 1298 VUS (≥1★).

## Mandatory validation (gate ②, passed)

- **ENST→NM equivalence:** VHL **12/12** and BAP1 **12/12** coding variants re-mapped
  via the MANE RefSeq NM gave coordinates identical to MaveDB's own post-mapping.
- **Known pathogenic splice → ClinVar + score direction:** both genes' canonical
  splice P/LP variants land in ClinVar and sit in the **most-damaging tail** of
  the assay's score distribution (VHL c.340+1G>T −3.21, ~min of VHL's −3.7..+1.5
  scale; BAP1 c.37+1G>T −0.21, ~1st-percentile of BAP1's compressed −0.27..+0.03
  scale). **BAP1's scale is compressed — must normalise per assay** (rank-based
  metrics like Spearman/AUROC are scale-invariant, so this is not a problem).
- **Give-up / unmapped:** VHL **0**, BAP1 **0** single-substitution variants
  unmapped (VHL 2268/2268 SNV; BAP1 9441 SNV — the rest are codon-level MNV/delins,
  excluded from SNV analysis exactly as for RAD51C/BARD1).

## 5. Circularity (gate ⑤, measured not assumed)

Fraction of functional SNVs already in ClinVar:
- **VHL: 1103/2268 = 49% → circularity_flag = FALSE**
- **BAP1: 3035/9441 = 32% → circularity_flag = FALSE**

Neither approaches BRCA1's ~100% — their SGE was **not** wholesale-deposited to
ClinVar, so their ClinVar labels are **independent**. This is a real bonus: both
genes *strengthen* the auxiliary binary (ClinVar) evaluation rather than adding
circular labels.

## Splice-N comparison (key decision table)

`results/tables/candidate_splice_comparison.tsv` — splice (core+region) functional SNVs:

| Gene | splice_core | splice_region | splice_total | group |
|---|--:|--:|--:|---|
| BRCA1 | 137 | 399 | 536 | existing |
| BRCA2 | 138 | 414 | 552 | existing |
| BARD1 | 107 | 323 | 430 | existing |
| PALB2 | 144 | 384 | 528 | existing |
| RAD51C | 96 | 288 | 384 | existing |
| **BAP1** | **192** | **576** | **768** | **CANDIDATE** |
| **VHL** | **24** | **69** | **93** | **CANDIDATE** |

Pooled splice N (existing 5) = **2430**.
- **+ BAP1 → 3198 (+32%)** — bigger single-gene splice contribution than any current gene.
- + VHL → 2523 (+4%).
- + both → **3291 (+35%)**.

## Recommendations (per gene)

| Gene | Verdict | Rationale |
|---|---|---|
| **BAP1** | **INCLUDE** | SGE ✓ (splice+intron+3′UTR); validated 12/12, 0 give-up; CC BY 4.0; **+32% pooled splice N**; non-circular (independent ClinVar labels); broadens panel to a chromatin/deubiquitinase tumour suppressor. Caveat: compressed score scale → normalise per assay. |
| **VHL** | **INCLUDE (breadth)** | SGE ✓; validated 12/12, 0 give-up; CC BY 4.0; non-circular; strong splice signal. Splice N gain small (+4%), but adds a **non-DNA-repair (hypoxia/VHL-HIF) tumour suppressor** → improves generalisability & journal scope. Include for breadth, not for N. |
| **MSH2** | **EXCLUDE** | No SGE — only a missense LOF screen (HAP1); zero splice/intron coverage, the project's banner region. |
| TP53 / PTEN / MLH1 / NF1 | **EXCLUDE** | No SGE (transactivation / VAMP-seq / abundance assays / none); cannot serve the splice/non-coding line. |

**Bottom line:** adding **BAP1 + VHL** takes the panel from 5 → 7 genes, raises
pooled splice N by **+35% (2430 → 3291)**, adds two independent-label (non-circular)
genes, and broadens beyond DNA-repair — all at SGE gold-standard quality matching
the existing panel. Recommend including both; BAP1 is the high-value add, VHL is a
low-cost generalisability add.

## Deliverables
- `results/tables/candidate_mavedb_catalog.tsv`, `candidate_gene_coverage.tsv`,
  `candidate_splice_comparison.tsv`
- `data/processed/candidate_functional_master.tsv` (20,376 variants; 11,709 SNV-mapped),
  `candidate_clinvar_variants.tsv` (VHL/BAP1 ClinVar)
- **Not yet merged into analysis_ready** — awaiting your include decision.

## Caveats
- BAP1 functional scores are on a compressed scale (−0.27..+0.03); use rank-based
  or per-assay-normalised metrics (the eval design already does).
- BAP1, like RAD51C/BARD1, has many codon-level MNV/delins (8,667 of 18,108) held
  out of the SNV analysis.
- proximal/deep intronic and splice thresholds identical to analysis_ready
  (|off|≤2 core, 3–8 region, >8 intronic).
