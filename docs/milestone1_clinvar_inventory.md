# Milestone 1 — ClinVar inventory

**Data:** ClinVar GRCh38 VCF, release **2026-06-15** (`data/raw/clinvar_release.txt`).
**Scope:** 4,436,216 ClinVar records scanned → **69,263** on the 6 target genes;
**66,044** at review stars ≥ 1 (the primary clinical-label reference).

## Clinical significance by gene (review stars ≥ 1)

| gene   | P/LP  | B/LB  | VUS   | Conflicting | Total (≥1★) |
|--------|------:|------:|------:|------------:|------------:|
| BRCA1  | 3926  | 3548  | 2331  | 2909        | 12715       |
| BRCA2  | 5453  | 5680  | 4155  | 5718        | 21006       |
| BARD1  | 608   | 1361  | 2104  | 430         | 4503        |
| RAD51C | 319   | 771   | 735   | 459         | 2284        |
| PALB2  | 1398  | 1650  | 2711  | 723         | 6482        |
| ATM    | 3701  | 5581  | 8277  | 1495        | 19054       |
| **TOTAL** | **15405** | **18591** | **20313** | **11734** | **66044** |

## Consequence classes (review stars ≥ 1)

Missense dominates (≈30k); the **splice + non-coding** subset — the project's
target gap that protein LMs cannot score — totals **15,958** variants
(P/LP 6314, B/LB 6689, VUS 1989).

| gene   | missense | synonymous | nonsense | frameshift | splice | noncoding |
|--------|---------:|-----------:|---------:|-----------:|-------:|----------:|
| BRCA1  | 5148 | 1720 | 863  | 43   | 380 | 4538 |
| BRCA2  | 9539 | 3713 | 1363 | 2944 | 341 | 2779 |
| BARD1  | 2273 | 899  | 190  | 0    | 112 | 1025 |
| RAD51C | 1075 | 405  | 76   | 33   | 81  | 607  |
| PALB2  | 3149 | 1211 | 373  | 565  | 140 | 970  |
| ATM    | 8529 | 3199 | 996  | 1194 | 673 | 4312 |

## Caveats / refinements for next pass
- **`noncoding` bucket is inflated** — it currently absorbs intronic + UTR SNVs
  *and* large deletions/duplications/CNVs annotated with an intron consequence.
  Before drawing gap conclusions we must split by `CLNVC` (keep
  `single_nucleotide_variant` + small indels; set aside structural/CNV events
  that LMs can't score) and separate true near-splice from deep-intronic.
- **Conflicting is large** (esp. BRCA2, BRCA1) — these are excluded from the
  clean P/B label set but worth a sensitivity analysis.
- The high splice/noncoding **P/LP** count needs sanity-checking against known
  pathogenic splice variants; likely includes gross deletions.
- Counts are GRCh38 only. GRCh37 liftover deferred until model scoring needs it.

## Outputs
- `data/processed/clinvar_target_variants.tsv` — variant master seed (all stars, 69,263 rows)
- `results/tables/clinvar_inventory_by_gene.tsv`
- `results/tables/clinvar_inventory_by_class.tsv`

## Next
1. Refine consequence/variant-type split (SNV vs indel vs CNV) — clean the gap subset.
2. **MaveDB / ProteinGym availability check** for BRCA1, BRCA2, BARD1, RAD51C, PALB2, ATM
   (which genes have downloadable SGE/MAVE functional scores — the gold standard).
