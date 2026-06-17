# Milestone 2 — Functional gold-standard (SGE/MAVE) availability & region coverage

**Question answered:** which of the 6 genes have downloadable functional gold
standards, *which regions those assays actually cover* (especially splice /
non-coding), and how many of our ClinVar target variants can be evaluated
against functional data. This decides the final gene scope and whether the
banner angle stays "non-coding" or narrows to "splice region".

Data sources (cached under `data/raw/mavedb/`, on D: via junction):
- **MaveDB** API `https://api.mavedb.org/api/v1` (score sets identified by URN;
  scores CSV + post-mapped VRS variants). ClinVar release for intersection: 2026-06-15.
- **ProteinGym** DMS substitutions reference (GitHub `OATML-Markslab/ProteinGym`).

---

## 1. Availability — SGE/MAVE per gene

5 of 6 genes have **Saturation Genome Editing (SGE)** score sets; **ATM has none**
in MaveDB (all 45 text hits were off-target "ATM pathway" mentions). Primary
SGE function-score set chosen per gene:

| Gene | URN | N variants | Transcript | Assembly | PMID | License |
|---|---|--:|---|---|---|---|
| BRCA1 | urn:mavedb:00000097-0-2 (Findlay 2018) | 3893 | NM_007294.3 | GRCh38 | 30209399 | CC0 |
| BARD1 | urn:mavedb:00001250-a-2 | 10911 | NM_000465.4 | GRCh38 | — | CC0 |
| BRCA2 | urn:mavedb:00001225-a-1 (2025) | 6959 | ENST00000380152.8 | GRCh38 | 39779857 | CC0 |
| PALB2 | urn:mavedb:00001259-a-2 | 12851 | NM_024675.4 | GRCh38 | — | CC0 |
| RAD51C | urn:mavedb:00000673-0-1 | 9188 | ENST00000337432.9 | GRCh38 | 39299233 | CC BY 4.0 |
| **ATM** | **none** | **0** | — | — | — | — |

Additional non-SGE MAVE assays exist (BRCA1: Starita yeast E3-ligase, HDR assays;
BRCA2: HDR) — see `results/tables/mavedb_catalog.tsv` (59 on-target score sets).
All selected sets are openly licensed (CC0 / CC BY 4.0).

## 2. The critical coordinate-mapping finding

MaveDB auto-maps variants against the **mature transcript** sequence (NM_/ENST),
which contains no introns. Therefore **every splice / intronic variant fails
MaveDB's post-mapping** (`errorMessage: "Variant is intronic and cannot be
processed"`, no genomic coords, no ClinGen ID) — i.e. *the project's banner
region is exactly what MaveDB cannot place*. Coding + UTR variants map fine.

Fix applied: splice/intronic transcript-HGVS → GRCh38 via **Mutalyzer**
(`NC_xxx(NM_xxx):c....` form). ENST-based sets (BRCA2, RAD51C) mapped through
their MANE-Select RefSeq equivalent (NM_000059.4 / NM_058216.3).

**Validation (do not skip — this was the #1 risk):**
- ENST→NM MANE equivalence: 12/12 coding BRCA2 & RAD51C variants re-mapped via
  the NM gave coordinates **identical** to MaveDB's own post-mapping.
- Known BRCA1 intronic SGE variants mapped to GRCh38 and **landed in the ClinVar
  master at matching coordinates** with consistent labels (deep-intronic +17–20
  → B/LB, near-zero functional scores).
- 142 variants (~2%) gave up after Mutalyzer retries; RAD51C/BARD1/PALB2 also
  carry codon-level **delins/MNV** variants (not SNVs) that are intentionally
  not mapped. Net: **34,582 / 43,802** functional variants have GRCh38 SNV coords.

## 3. Region coverage & ClinVar overlap (deliverable: `goldstandard_coverage_matrix.tsv`)

Functional variants do cover **splice (±1/2) and intronic** positions — SGE tiles
the targeted introns. Overlap = functional SNV that is also a ClinVar SNV.

**Splice + intronic + non-coding subset (the gap), overlapping ClinVar ≥1★:**

| Gene | functional SNVs (gap) | overlap ClinVar ≥1★ | of which **VUS** |
|---|--:|--:|--:|
| BRCA1 | 1090 | 466 | 46 |
| PALB2 | 1519 | 428 | 80 |
| BARD1 | 2328 | 424 | 78 |
| RAD51C | 1353 | 371 | 64 |
| BRCA2 | 689 | 322 | 56 |

Splice-specifically (ClinVar SO `splice`, ≥1★, with functional score): BRCA1 113,
BRCA2 108, PALB2 99, BARD1 71, RAD51C 53 → **~444 splice variants with an
independent functional gold standard** across 5 genes.

**Totals:** 34,582 functional SNVs → **15,106 overlap ClinVar ≥1★**, of which
**6,089 are VUS** (the reclassification opportunity). Coding overlap by ClinVar
SO label is caliber-consistent with Milestone 1 (missense/synonymous/nonsense
split in the matrix's secondary view).

## 4. ProteinGym (deliverable: `proteingym_coverage.tsv`)

**Missense / protein-level ONLY — cannot serve the splice/non-coding line.**
Only 2 of 6 genes present: BRCA1 (`BRCA1_HUMAN_Findlay_2018`, 1837 single
mutants — the missense slice of the same SGE) and BRCA2
(`BRCA2_HUMAN_Erwood_2022`, 265). BARD1/RAD51C/PALB2/ATM: absent. Use ProteinGym
only as a convenience missense benchmark for protein LMs, not for the gap.

## 5. Recommendations (the decisions this milestone informs)

1. **Final gene scope = BRCA1, BRCA2, BARD1, PALB2, RAD51C** (all have SGE
   covering coding + splice + proximal intronic, openly licensed, GRCh38-mappable).
2. **ATM: drop from the functional-gold-standard arm.** Keep only for an
   auxiliary ClinVar-label evaluation, explicitly flagged for circularity
   (predictors trained on ClinVar will look inflated).
3. **Narrow the banner from "non-coding" to "splice region + proximal intronic".**
   That is where functional gold standard genuinely exists (~444 splice +
   ~thousands proximal-intronic variants). **Deep/distal non-coding and promoter
   variants have essentially no functional gold standard** and few ClinVar VUS —
   do not over-claim there.
4. **Circularity caveat specific to BRCA1:** Findlay SGE results were deposited to
   ClinVar, so BRCA1 ClinVar labels are *not* independent of the functional data
   (note the near-100% BRCA1 overlap). The primary evaluation must use the
   continuous functional score as the independent standard; ClinVar labels are
   auxiliary only.
5. Region with functional gold standard to evaluate on: **coding (missense/syn/
   nonsense), splice ±1/2, proximal intronic**, plus partial 3′/5′UTR for BARD1/
   PALB2 (their SGE tiled the UTR). This is the concrete evaluable surface for
   the DNA-LM-vs-protein-LM comparison.

## Deliverables
- `data/processed/functional_scores_master.tsv` — 43,802 functional variants;
  34,582 with GRCh38 chrom/pos/ref/alt, region_class, functional_score, source_urn, assay_type.
- `results/tables/goldstandard_coverage_matrix.tsv` — gene × region with functional
  count + ClinVar overlap (any/≥1★/VUS).
- `results/tables/mavedb_catalog.tsv`, `results/tables/proteingym_coverage.tsv`.

## Open items / caveats
- RAD51C & BARD1 SGE include many codon-level delins (MNV) not mapped as SNVs;
  if needed for a multi-nucleotide analysis they require separate handling.
- Splice vs intronic threshold here is |intron offset| ≤ 2 = splice; refine to
  the splice-region (±1–3 exon / +1..+8, −3..−? intron) in Phase 3 if needed
  (offset is stored in the master).
- 142 Mutalyzer give-ups can be re-run later (cache is incremental).
