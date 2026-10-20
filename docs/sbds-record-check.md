# SBDS: MaveDB record (fifth new external gene)

Written by `phase1/src/evid_newgenes.py` (stage E16b) from the files `scripts/stage_E16/fetch_mavedb.py` saved; do not edit by hand. Each fact is followed by the file it was read from (under `phase1/data/evidence/external/`) and the time that file was fetched.

## Pre-specified label source (预先规定)

Fixed by the author on 2026-10-20, before any predictor was run on SBDS, and not to be changed afterwards:

- SBDS labels come from MaveDB only, as for every external gene so far (DDX3X, TP53 and the four genes of E16). The IGVF copy of the same experiment (`igvf_sbds/IGVFFI9790OCYB.tsv.gz`) is used only for a cross-check (`phase1/reports/evidence/sbds_igvf_crosscheck.csv`) and never to build a label.
- Why the source had to be fixed in advance: in the author's SGE survey (2026-10-20) the 3-10 nt band holds 10 functionally abnormal SNVs by MaveDB and 9 by IGVF, either side of the evaluability threshold (at least 10 damaging), so the choice of source could decide whether the band is evaluated at all.
- Order of class sources: (a) exactly one investigator-provided set of ranges, checked variant by variant against a class column if the score table has one; (b) otherwise a per-variant class column of the score table, copied as it stands; (c) otherwise stop. SBDS has no investigator-provided ranges, so (b) applies: the column `functional_consequence`.
- Class names map by wording only: functionally_abnormal -> 1; functionally_normal -> 0; indeterminate -> NA (row kept); any other value -> unassigned, listed.

## SBDS

- Score set: urn:mavedb:00001284-a-1, "Saturation Genome Editing for SBDS Scores" [`sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]
- Description: We used Saturation Genome Editing to measure the effect of 2,675 SBDS SNVs and 640 3-base pair deletions on SBDS-dependent HAP1 cell survival. [`sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]
- Status: public (private = false), processing state success, mapping state incomplete [`sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]
- Licence: CC0 1.0 (https://creativecommons.org/publicdomain/zero/1.0/) [`sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]
- Publication links: score set not_provided [`sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]; experiment urn:mavedb:00001284-a not_provided [`sbds_raw/experiment_meta.json`, fetched 2026-10-20T10:39:18+00:00]
- Dates: this version created 2026-08-24, published 2026-08-24 [`sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]; first published 2026-08-24, the earliest publication date of its versions [`sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]
- Versions: urn:mavedb:00001284-a-1 (this one), published 2026-08-24, 3,315 variants, mapping incomplete, investigator-provided classes  [`sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]
- Checked against the author's record of 2026-10-20 (URN, licence, publication status, first published 2026-08, the current version dated separately): agrees: urn urn:mavedb:00001284-a-1, license CC0, private False, score_set_publications not_provided, experiment_publications not_provided, first_published 2026-08; first published 2026-08-24, current version published 2026-08-24 [`sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]
- Target: SBDS, ENST00000246868.7 (not_provided) [`sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]
- Variants: 3,315 in the metadata [`sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]; 3,315 rows in the score table [`sbds_raw/urn_mavedb_00001284-a-1.csv`, fetched 2026-10-20T10:39:19+00:00]
- Investigator-provided ranges: none (no scoreCalibrations entry has investigatorProvided true) [`sbds_raw/score_ranges.json`, extracted from `sbds_raw/scoreset_meta.json`, fetched 2026-10-20T10:39:18+00:00]
- Classes: copied from the score-table column `functional_consequence`, mapped by wording only: functionally_abnormal -> abnormal; functionally_normal -> normal; indeterminate -> not_specified (label 1 / 0 / NA); values in the column: functionally_abnormal, functionally_normal, indeterminate [`sbds_raw/urn_mavedb_00001284-a-1.csv`, fetched 2026-10-20T10:39:19+00:00]
