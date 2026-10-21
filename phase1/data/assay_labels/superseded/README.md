# Superseded label files

`palb2_mavedb.csv` and `bard1_woo2025.tsv` here are the PALB2 and BARD1 label files
that entered the repository on 2026-07-08 (commit 10d94bb), unchanged.

The label source audit (stage E16c, `docs/label-source-audit.md`) could not reproduce
them from the score sets Supplementary Table S1 cites. Under the investigator-provided
ranges of the cited urn:mavedb:00001259-a-2 and urn:mavedb:00001250-a-2 they differ on
359 and 13 variants.

On 2026-10-21 the author replaced them, under the same names one folder up, with files
built from those ranges (`scripts/stage_E16d/build_labels.py`; the record is
`../PROVENANCE.json`). The decision was made before the refit was seen.

These two are kept for the sensitivity comparison (revision item C3). No stage reads
them.
