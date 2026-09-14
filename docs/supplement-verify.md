# Supplementary tables: verification before rewriting

`phase1/src/build_supp_tables.py` writes every pipeline-derived supplementary table from the report tables. Before it
was allowed to write the frozen-matrix-v2 tables it was run in `--verify` mode against the frozen-matrix-v1 outputs and
the supplement as it stood before the revision: it regenerates each table in memory and compares it with the
document cell by cell, writing nothing. A cell whose printed precision changed is compared as a number. A column the
generator adds is reported beside the count rather than in it, because it has no cell in the older document to
disagree with. The same check is then run against the v2 outputs and the current supplement.

This page records both runs; it is regenerated whenever the generator or the supplement changes.

## 1. v1 outputs against the pre-revision supplement

```
cd phase1
python -m src.build_supp_tables --in "BiB supplementary_material_backup_20260914_preV2.docx" \
    --reports reports/phase1_v1 --verify
```

Supplement sha256 `1d7c485a85daa2f4b80e58abde399637b81bfb9b12798a8277db6d8c6383c333`.

```
[verify] Supplementary Table S2: 10 rows, 0 mismatch(es)
[verify] Supplementary Table S3a: 8 rows, 0 mismatch(es)
[verify] Supplementary Table S3b: 96 rows, 0 mismatch(es)
[verify] Supplementary Table S4: 80 rows, 0 mismatch(es)
[verify] Supplementary Table S5a: 4 rows, 0 mismatch(es)
[verify] Supplementary Table S5b: 13 rows, 0 mismatch(es)
[verify] Supplementary Table S6: 7 rows, 0 mismatch(es)
[verify] Supplementary Table S8a: 4 rows, 0 mismatch(es)
[verify] Supplementary Table S8b: 3 rows, 0 mismatch(es)
[verify] Supplementary Table S10a: 8 rows, 0 mismatch(es)
[verify] Supplementary Table S10b: 5 rows, 0 mismatch(es)
[verify] Supplementary Table S10c: 11 rows, 0 mismatch(es)
[verify] Supplementary Table S10d: 20 rows, 0 mismatch(es)
[verify] Supplementary Table S11: 96 rows, 0 mismatch(es); 5 new column(s) not in the document
[verify] Supplementary Table S12: 9 rows, 0 mismatch(es)
[verify] Supplementary Table S13: 4 rows, 0 mismatch(es)
[verify] Supplementary Table S14: 12 rows, 0 mismatch(es)
[verify] Supplementary Table S15 (a): 2 rows, 0 mismatch(es)
[verify] Supplementary Table S15 (b): 12 rows, 0 mismatch(es)
[verify] Supplementary Table S16 (a): 10 rows, 0 mismatch(es)
[verify] Supplementary Table S16 (b): 21 rows, 0 mismatch(es)
[verify] Supplementary Table S16 (c): 13 rows, 0 mismatch(es)
[verify] Supplementary Table S17: 12 rows, 0 mismatch(es)
[verify] Supplementary Table S19: 24 rows, 0 mismatch(es); 4 new column(s) not in the document
[supp] Supplementary Table S20: not in document
[verify] total mismatches: 0
```

## 2. v2 outputs against the current supplement

```
python -m src.build_supp_tables --in "BiB supplementary_material.docx" --reports reports/phase1 --verify
```

Supplement sha256 `ba8b5bb9471ae659c424bd915db9653822311a2db00129aadf8dcc14f55df310`.

```
[verify] Supplementary Table S2: 10 rows, 0 mismatch(es)
[verify] Supplementary Table S3a: 8 rows, 0 mismatch(es)
[verify] Supplementary Table S3b: 96 rows, 0 mismatch(es)
[verify] Supplementary Table S4: 80 rows, 0 mismatch(es)
[verify] Supplementary Table S5a: 4 rows, 0 mismatch(es)
[verify] Supplementary Table S5b: 13 rows, 0 mismatch(es)
[verify] Supplementary Table S6: 7 rows, 0 mismatch(es)
[verify] Supplementary Table S8a: 4 rows, 0 mismatch(es)
[verify] Supplementary Table S8b: 3 rows, 0 mismatch(es)
[verify] Supplementary Table S10a: 8 rows, 0 mismatch(es)
[verify] Supplementary Table S10b: 5 rows, 0 mismatch(es)
[verify] Supplementary Table S10c: 11 rows, 0 mismatch(es)
[verify] Supplementary Table S10d: 20 rows, 0 mismatch(es)
[verify] Supplementary Table S11: 96 rows, 0 mismatch(es)
[verify] Supplementary Table S12: 9 rows, 0 mismatch(es)
[verify] Supplementary Table S13: 4 rows, 0 mismatch(es)
[verify] Supplementary Table S14: 12 rows, 0 mismatch(es)
[verify] Supplementary Table S15 (a): 2 rows, 0 mismatch(es)
[verify] Supplementary Table S15 (b): 12 rows, 0 mismatch(es)
[verify] Supplementary Table S16 (a): 10 rows, 0 mismatch(es)
[verify] Supplementary Table S16 (b): 21 rows, 0 mismatch(es)
[verify] Supplementary Table S16 (c): 13 rows, 0 mismatch(es)
[verify] Supplementary Table S17: 12 rows, 0 mismatch(es)
[verify] Supplementary Table S19: 24 rows, 0 mismatch(es)
[verify] Supplementary Table S20: 96 rows, 0 mismatch(es)
[verify] total mismatches: 0
```
