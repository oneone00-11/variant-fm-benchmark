# Milestone 3 — Analysis-ready dataset & a-priori precision

Pure data assembly + analytical power priors. **No model scoring was run.** All
precision numbers below are *prior* estimates under *assumed* effect sizes — real
ρ/AUROC come only after Phase-2 scoring.

Final scope: **BRCA1, BRCA2, BARD1, PALB2, RAD51C** (ATM functional-excluded;
ClinVar-only auxiliary). Coordinates consistent with Milestone 2 (GRCh38).

## A. Cleaning & merge

- **ClinVar SNV-only:** kept 54,977 `single_nucleotide_variant`; **excluded 14,286
  non-SNV** → `data/processed/excluded_nonSNV.tsv` (Deletion 8245, Duplication 2978,
  Indel 1159, Insertion 985, Microsatellite 891, Inversion 28).
- **Functional SNV-only:** 34,582 SNV-coord kept; **9,220 codon-level delins/MNV or
  unmapped held out** → `data/processed/functional_nonSNV.tsv`.
- **Inner join on (chrom,pos,ref,alt):** **17,272 variants** with BOTH a functional
  score and a ClinVar record → `data/processed/analysis_ready.tsv`.
  - label_set: **clean_PB 6,447** (P/LP+B/LB) · **VUS 6,171** · conflicting 2,608 · other 2,046.
  - Conflicting and VUS are retained, distinguished only by `label_set`.

**Region rule (unified, documented).** Intron-side variants use the functional
`intron_offset`: |off|≤2 → splice/`splice_core`; 3–8 → splice/`splice_region`;
>8 → intronic (`proximal` ≤50, else `deep`). Exonic/UTR use ClinVar SO mc_terms
(missense/synonymous/nonsense/utr). **Caveat:** ClinVar's MC has no
`splice_region` term, so intron-side splice-region comes from the offset; the
**exon-side last-1-3 splice-region is not detectable** from either source (those
fall under missense/synonymous) — a known, documented gap.

**Join spot-check (passed).** BRCA1 canonical splice P/LP variants land at the
expected coords with strongly damaging functional scores (c.5194-1G>C −4.20,
c.5407-1G>A −3.90, c.4986+1G>A −3.35); deep-intronic B/LB are near-zero. Join is
correct and biologically coherent.

**circularity_flag basis (not assumed — checked).** Per-gene fraction of
functional SNVs already present in ClinVar: **BRCA1 100%** (3893/3893), BRCA2 54%,
BARD1 39%, PALB2 42%, RAD51C 40%. Only BRCA1 shows near-total overlap, consistent
with the documented deposition of the Findlay SGE into ClinVar → **BRCA1 = TRUE**,
others = FALSE. (Partial overlap elsewhere reflects pre-existing clinical ClinVar
entries, not wholesale functional deposition.)

## B. Sample size & precision (key deliverable)

Denominators: **Spearman** (continuous, predictor vs functional_score) needs only
a functional score → N = full functional-SNV set in the region. **AUROC** (binary,
predictor vs ClinVar clean_PB) → n_pos=P/LP, n_neg=B/LB from the merged set.
Thresholds: Spearman CI half-width ≤0.10 and AUROC ≤0.05 = "sufficient".

### Pooled (across 5 genes) — the headline

| subset | N (Spearman) | Spearman ½-width (ρ=0.5) | n_path / n_benign | AUROC ½-width (0.85) | verdict |
|---|--:|--:|--:|--:|---|
| **whole gap (splice+intronic)** | 6,304 | **0.019** | 471 / 926 | **0.023** | both sufficient |
| **splice (core+region)** | 2,430 | **0.030** | 463 / 280 | **0.027** | both sufficient |
| splice_core (±1/2) | 622 | 0.059 | 397 / **0** | — (no benign) | continuous sufficient; binary undefined |
| splice_region (+3..+8) | 1,808 | 0.035 | 66 / 280 | 0.060 | continuous suff.; binary marginal |
| intronic (>8) | 3,874 | 0.024 | 8 / 646 | 0.168 | continuous suff.; binary undefined |
| coding | 14,592 | 0.012 | 1,111 / 3,899 | 0.015 | both sufficient |

### Per-gene splice (banner region)

Continuous: every gene **sufficient** (½-width 0.064–0.075, N 384–536).
Binary: every gene **marginal** (½-width ~0.05–0.076). `splice_core` per gene is
N≈55–148 → continuous **marginal/insufficient** (½-width 0.12–0.15) and binary
**undefined** (zero benign). Full per-gene table in `results/tables/precision_estimate.tsv`.

### Direct answers to "is N enough?"

1. **The whole gap and splice (the banner) are sufficiently powered with the
   current 5 genes — pooled.** Continuous ½-widths 0.019 (gap) / 0.030 (splice)
   and balanced binary ½-widths 0.023 / 0.027 are all well under threshold.
   **No additional genes are required for adequate precision on the gap.**
2. **Per-gene splice is fine for the continuous (functional) standard**, marginal
   for the binary one — so report the gap **pooled**, per-gene as secondary.
3. **`splice_core` per gene cannot be made "sufficient" by adding variants** — a
   gene has only ~100–150 canonical splice positions (biological ceiling). Pooling
   the 5 genes (N=622) already reaches sufficient; more genes only raise pooled N,
   not needed.
4. **The real limitation is class balance, not N:** canonical splice sites are
   essentially never benign (splice_core n_benign=0) and assayed deep-intronic
   SNVs are almost never pathogenic (intronic n_path=8). This is biology, not a
   sampling defect — so:
   - use the **continuous functional score as the primary gold standard** for the
     gap (independent, well-powered);
   - run any **binary ClinVar AUROC pooled**, and report it **with/without BRCA1**
     (circular) as a sensitivity analysis.
5. If a future reviewer demands per-gene binary splice power (≤0.05), that needs a
   balanced P/B splice set that does not exist for these genes; the answer is to
   lean on the continuous standard, not to add ~N more genes.

## Deliverables
- `data/processed/analysis_ready.tsv` — 17,272 rows (gene, coords, region_class,
  splice_class, intronic_depth, intron_offset, functional_score, clnsig_class,
  stars, path_binary, label_set, circularity_flag, source_urn, assay_type, …).
- `data/processed/excluded_nonSNV.tsv` (14,286) · `functional_nonSNV.tsv` (9,220).
- `results/tables/sample_size_matrix.tsv` — gene × region × label_set N + balance.
- `results/tables/precision_estimate.tsv` — Spearman & AUROC CI ½-widths.

## Caveats
- Prior precision only; assumes ρ∈{0.3,0.5,0.7}, AUROC∈{0.75,0.85,0.90}.
- Spearman SE = 1/√(N−3) (Fisher z); AUROC SE = Hanley–McNeil with actual balance.
- proximal/deep intronic split at |offset|≤50 is adjustable (offset stored).
- exon-side splice-region (last 1–3 exonic nt) not separately labeled.
