# AlphaGenome by its official definition: the pilot

Written by `scripts/stage_E17/pilot_report.py` from `phase1/data/evidence/alphagenome_official/`; do not edit by hand. Stage E17, commit 1. It computes no likelihood ratio, threshold or tier (E17.0 rule 9): it shows whether the implementation is the notebook's, not how the score does.

## The definition and where it comes from

- Notebook: `colabs/v1/splicing_variant_scoring.ipynb` in https://github.com/google-deepmind/alphagenome, commit `55943112b83439e0c2512110f1ae607d4b5b492a` (the last commit touching the file: `be1dfff647f1`, 2026-10-09T15:18:43Z, "Update alphagenome model imports to use v1 modules.").
- sha256 `c39e6a6a09bc62512107c106cd43cabf7a016cbd073474ff872c669a1ad838df` (13,092 bytes); read 2026-10-22T05:18:48+00:00. The repository's licence: `LICENSE`, "Apache License Version 2.0, January 2004" (Apache-2.0), sha256 `cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30`.
- Printed outputs in the notebook: 0 (none: there is no reference output to reproduce).
- The quoted cells, verbatim, with their numbers: `phase1/data/evidence/alphagenome_official/NOTEBOOK_SOURCE.json` (cells 3, 5, 7, 9, 11, 12, 13, 15).

| item | what the quoted code fixes | cell:line |
|---|---|---|
| client | alphagenome, unpinned in the notebook: the newest release on PyPI on the day of scoring | 3:5 |
| client call | dna_client.create(api key), then dna_model.score_variant(interval, variant, variant_scorers, organism) per variant. The notebook reads the key with colab_utils.get_api_key(); here it comes from phase1/src/evid_score_avi.read_key() (E17.0 rule 1), the only line of the quoted code that is not run as written | 7:1, 11:13, 15:27 |
| scorers | variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_SITES'], ['SPLICE_SITE_USAGE'] and ['SPLICE_JUNCTIONS'] | 9:2, 9:3, 9:4 |
| input window | dna_client.SEQUENCE_LENGTH_1MB nucleotides, an interval of that length centred on the variant (variant.reference_interval.resize) | 11:10, 15:26 |
| organism | dna_client.Organism.HOMO_SAPIENS | 11:17, 15:31 |
| tissue or ontology | none: score_variant is called with no ontology selection, and the maximum runs over every track | 11:16, 13:21 |
| score field | raw_score, not quantile_score | 13:21 |
| aggregation | for each variant and scorer, the maximum of raw_score over all its rows: every gene, every track and, for the junction scorer, every junction. The three scorers are unsigned (is_signed False), so this maximum is the maximum absolute score the notebook's text names | 13:21, 13:22, 9:8, 12:4 |
| combination | SPLICE_SITES + SPLICE_SITE_USAGE + SPLICE_JUNCTIONS / 5.0 | 13:30, 13:31, 13:32 |
| missing values | a NaN raw_score is skipped by the maximum (pandas); a scorer with no row for a variant counts 0.0 (fillna), and one with no row for any variant 0.0 (the .get default) | 13:22, 13:25, 13:30 |
| strand | none chosen: rows of genes on either strand all enter the maximum | 13:21 |
| variant | genome.Variant(chromosome 'chrN', position 1-based, reference_bases, alternate_bases) | 11:2, 11:3, 11:4 |

Read from the installed client (alphagenome 0.10.0): `SEQUENCE_LENGTH_1MB` = 1,048,576; organism `Organism.HOMO_SAPIENS`; the scorers `GeneMaskSplicingScorer(requested_output=SPLICE_SITES, width=None)` (signed: False); `GeneMaskSplicingScorer(requested_output=SPLICE_SITE_USAGE, width=None)` (signed: False); `SpliceJunctionScorer()` (signed: False).

Arithmetic: the client returns raw scores as float32, so cell 13 takes the maxima and adds and divides in float32. The components and the composite are float32 values, written exactly (each number reads back to the same bits), and the composite is (sites + usage) + junctions / 5 in float32 to the last bit (`phase1/tests/test_evid_alphagenome_official.py`). Read them back exactly: pandas' default float parser is not round-trip exact.

## The scoring environment

- Client: alphagenome 0.10.0. The notebook installs it unpinned (cell 3), so it is the newest release on PyPI on the day of scoring; `scripts/stage_E17/requirements.txt` pins it and its 49 packages, each with the sha256 of its file (sha256 of the file: `973a446bfe17437bfac7fd7f35c577c641cf19161f4c24669bd6fc5e51c30895`), and `.venv-ag-official/` (gitignored) is built from it.
- Python: 3.12.13; macOS-26.6.2-arm64-arm-64bit.
- Model version: requested None (dna_client.create(...) as in the notebook's cell 7, with no model_version: the service's default); in the response: not_provided (the response's metadata keys: interval, variant, variant_scorer).
- API key: present True, from `~/.config/alphagenome/key`; its value is recorded nowhere.

## The notebook's own examples

The notebook prints no output (无参照输出), so its three examples have no reference to match. Scored here with the quoted code as written (cell 7's key aside), they show that code running on this client:

| variant | SPLICE_SITES | SPLICE_SITE_USAGE | SPLICE_JUNCTIONS | alphagenome_splicing |
|---|---|---|---|---|
| `chr13:32316462:T>G` | 0.01129150390625 | 0.015625 | 0.09112548828125 | 0.04514160007238388 |
| `chr17:43092919:G>A` | 0.16748046875 | 0.02203369140625 | 0.4111328125 | 0.27174073457717896 |
| `chr7:117559593:T>A` | 0.0078125 | 0.01171875 | 0.04608154296875 | 0.02874755859375 |

## The pilot sample

120 variants of `analysis_set_v1`, an even share in each non-empty gene x stratum cell, drawn with `config.RANDOM_SEED` (20260708) by `scripts/stage_E17/score_official.py`:

| gene | pm12 | s3_10 | s11_50 |
|---|---|---|---|
| BAP1 | 6 | 6 | 6 |
| BARD1 | 6 | 6 | 6 |
| BRCA1 | 6 | 6 | 6 |
| BRCA2 | 6 | 6 | 0 |
| PALB2 | 6 | 6 | 6 |
| RAD51C | 6 | 6 | 6 |
| VHL | 6 | 6 | 6 |

## Determinism

10 of the pilot's variants (drawn with the same seed) were scored again 32 minutes after the pilot ended (2026-10-22T05:24:02+00:00 -> 2026-10-22T05:55:35+00:00). Identical to the last bit in all four values: **10 of 10**.

| variant | SPLICE_SITES | SPLICE_SITE_USAGE | SPLICE_JUNCTIONS | alphagenome_splicing | identical |
|---|---|---|---|---|---|
| `urn:mavedb:00000097-0-2#1593` | 0.9958934783935547 | 0.9175510406494141 | 5.7890625 | 3.0712571144104004 | yes |
| `urn:mavedb:00000097-0-2#1768` | 0.95458984375 | 0.929229736328125 | 7.1015625 | 3.3041319847106934 | yes |
| `urn:mavedb:00000097-0-2#3367` | 0.9915618896484375 | 0.9409332275390625 | 5.109375 | 2.9543700218200684 | yes |
| `urn:mavedb:00000662-0-1#5167` | 0.99237060546875 | 0.9669189453125 | 7.03125 | 3.36553955078125 | yes |
| `urn:mavedb:00000662-0-1#5747` | 0.10546875 | 0.0625 | 0.81787109375 | 0.33154296875 | yes |
| `urn:mavedb:00000662-0-1#928` | 0.03515625 | 0.0078125 | 0.11553955078125 | 0.06607665866613388 | yes |
| `urn:mavedb:00000673-0-1#6462` | 0.01171875 | 0.009765625 | 0.1363525390625 | 0.04875488206744194 | yes |
| `urn:mavedb:00000675-a-1#1088` | 0.1796875 | 0.6982421875 | 3.26953125 | 1.531835913658142 | yes |
| `urn:mavedb:00000675-a-1#1571` | 0.9954566955566406 | 0.9760093688964844 | 6.9296875 | 3.357403516769409 | yes |
| `urn:mavedb:00001250-a-2#3447` | 0.0078125 | 0.01171875 | 0.03656005859375 | 0.02684326097369194 | yes |

## Missing values

| component | variants with no row | NaN raw scores | smallest raw score |
|---|---|---|---|
| SPLICE_SITES | 0 | 0 | 0.00244141 |
| SPLICE_SITE_USAGE | 0 | 0 | 9.53674e-05 |
| SPLICE_JUNCTIONS | 0 | 0 | 1.60933e-06 |

Composite missing: 0 of 120.

## Distribution

| value | min | 25% | median | 75% | max |
|---|---|---|---|---|---|
| SPLICE_SITES | 0.00390625 | 0.0131226 | 0.0944824 | 0.991447 | 0.99979 |
| SPLICE_SITE_USAGE | 0.00390625 | 0.0185547 | 0.0585938 | 0.930346 | 0.983437 |
| SPLICE_JUNCTIONS | 0.0317993 | 0.286743 | 1.45703 | 5.5166 | 9.05469 |
| alphagenome_splicing | 0.0219849 | 0.100391 | 0.472412 | 2.993 | 3.7562 |

| stratum | n | median composite | max composite |
|---|---|---|---|
| pm12 | 42 | 3.07063 | 3.7562 |
| s3_10 | 42 | 0.267188 | 3.49795 |
| s11_50 | 36 | 0.0785619 | 1.26445 |

## Spearman correlation with the other two AlphaGenome columns

Over the pilot's variants, values only, no threshold:

| with | n | Spearman rho |
|---|---|---|
| `alphagenome` (the atlas quantile composite; `alphagenome_atlas_quantile` from commit 3) | 120 | 0.943 |
| `alphagenome_v061` (client 0.6.1) | 57 | 0.994 |

## Throughput

The pilot scored 120 variants in 74 s with 2 workers (1.61 variants a second), with 0 retries and 0 failures. At that rate the whole analysis set, 8,853 variants, takes about 92 minutes and 8,853 calls of `score_variant`, one per variant, before any retry.

## The four Atlas columns: raw scores or quantiles (read-only)

`avi`, `avi_splice_sites`, `avi_splice_site_usage`, `avi_splice_junctions`: **raw: each column is the maximum, over tracks and rows, of X, which client 0.9.0 fills with the service's scores; the calibrated scores (quantiles) sit in the 'quantiles' layer, which the columns do not read**.
- The record: `phase1/data/evidence/avi.parquet.provenance.json` (client 0.9.0, scorers AVI_SCORE, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS; aggregation: "per variant, the maximum over tracks and over rows; SPLICE_JUNCTIONS returns one row per junction").
- Our code: `phase1/src/evid_score_avi.py`, `_aggregate`, reads `ad.X`.
- The client's code: `alphagenome/atlas/atlas.py` in `alphagenome-0.9.0-py3-none-any.whl` (sha256 `a4f35884341ae85b5d2cf088dfe0304961de7ae4a590d4538653069673de32f4`):
  - line 174, the service's scores, as float32: `x = np.frombuffer(score.scores, dtype=np.float32).reshape(score.shape)`
  - line 181, its calibrated scores: `score.calibrated_scores, dtype=np.float32`
  - line 219, the calibrated scores go to a layer named 'quantiles': `layers = {'quantiles': np.concatenate(calibrated_scores)}`
  - line 223, X is the service's scores: `X=scores, obs=obs, var=var, layers=layers`

## Every pilot variant

| variant | gene | stratum | SPLICE_SITES | SPLICE_SITE_USAGE | SPLICE_JUNCTIONS | alphagenome_splicing |
|---|---|---|---|---|---|---|
| `urn:mavedb:00000097-0-2#121` | BRCA1 | s11_50 | 0.0327148 | 0.0117188 | 0.167847 | 0.0780029 |
| `urn:mavedb:00000097-0-2#1335` | BRCA1 | s3_10 | 0.0246582 | 0.0351562 | 0.894531 | 0.238721 |
| `urn:mavedb:00000097-0-2#1452` | BRCA1 | s11_50 | 0.0078125 | 0.0078125 | 0.0463867 | 0.0249023 |
| `urn:mavedb:00000097-0-2#1593` | BRCA1 | pm12 | 0.995893 | 0.917551 | 5.78906 | 3.07126 |
| `urn:mavedb:00000097-0-2#1604` | BRCA1 | s3_10 | 0.0129395 | 0.03125 | 0.290771 | 0.102344 |
| `urn:mavedb:00000097-0-2#1768` | BRCA1 | s3_10 | 0.95459 | 0.92923 | 7.10156 | 3.30413 |
| `urn:mavedb:00000097-0-2#2069` | BRCA1 | s3_10 | 0.00976562 | 0.015625 | 1.14062 | 0.253516 |
| `urn:mavedb:00000097-0-2#2324` | BRCA1 | s3_10 | 0.0146484 | 0.0195312 | 0.107239 | 0.0556274 |
| `urn:mavedb:00000097-0-2#3230` | BRCA1 | pm12 | 0.995548 | 0.949129 | 6.40234 | 3.22515 |
| `urn:mavedb:00000097-0-2#3367` | BRCA1 | pm12 | 0.991562 | 0.940933 | 5.10938 | 2.95437 |
| `urn:mavedb:00000097-0-2#3383` | BRCA1 | s3_10 | 0.00976562 | 0.0117188 | 0.0420227 | 0.0298889 |
| `urn:mavedb:00000097-0-2#3419` | BRCA1 | s11_50 | 0.0078125 | 0.0117188 | 0.0600586 | 0.031543 |
| `urn:mavedb:00000097-0-2#348` | BRCA1 | s11_50 | 0.0859375 | 0.0234375 | 0.357178 | 0.180811 |
| `urn:mavedb:00000097-0-2#627` | BRCA1 | pm12 | 0.995972 | 0.925702 | 7.64062 | 3.4498 |
| `urn:mavedb:00000097-0-2#662` | BRCA1 | s11_50 | 0.0078125 | 0.0078125 | 0.112061 | 0.0380371 |
| `urn:mavedb:00000097-0-2#667` | BRCA1 | s11_50 | 0.0078125 | 0.0078125 | 0.0498657 | 0.0255981 |
| `urn:mavedb:00000097-0-2#964` | BRCA1 | pm12 | 0.913574 | 0.951813 | 5.95703 | 3.05679 |
| `urn:mavedb:00000097-0-2#968` | BRCA1 | pm12 | 0.995338 | 0.951958 | 6.6875 | 3.2848 |
| `urn:mavedb:00000662-0-1#1452` | BAP1 | s11_50 | 0.0493164 | 0.0078125 | 0.225586 | 0.102246 |
| `urn:mavedb:00000662-0-1#2054` | BAP1 | s3_10 | 0.106689 | 0.046875 | 2.80078 | 0.713721 |
| `urn:mavedb:00000662-0-1#2254` | BAP1 | pm12 | 0.996429 | 0.983437 | 8.07812 | 3.59549 |
| `urn:mavedb:00000662-0-1#3363` | BAP1 | pm12 | 0.995988 | 0.952922 | 4.42578 | 2.83407 |
| `urn:mavedb:00000662-0-1#3870` | BAP1 | s3_10 | 0.181641 | 0.367188 | 2.15625 | 0.980078 |
| `urn:mavedb:00000662-0-1#397` | BAP1 | s11_50 | 0.0078125 | 0.0195312 | 0.245972 | 0.0765381 |
| `urn:mavedb:00000662-0-1#4467` | BAP1 | pm12 | 0.99884 | 0.962723 | 7.13672 | 3.38891 |
| `urn:mavedb:00000662-0-1#5157` | BAP1 | s3_10 | 0.0585938 | 0.0117188 | 0.489502 | 0.168213 |
| `urn:mavedb:00000662-0-1#5167` | BAP1 | pm12 | 0.992371 | 0.966919 | 7.03125 | 3.36554 |
| `urn:mavedb:00000662-0-1#547` | BAP1 | pm12 | 0.987305 | 0.963684 | 4.08203 | 2.7674 |
| `urn:mavedb:00000662-0-1#5716` | BAP1 | s11_50 | 0.0234375 | 0.0078125 | 0.0534668 | 0.0419434 |
| `urn:mavedb:00000662-0-1#5747` | BAP1 | s11_50 | 0.105469 | 0.0625 | 0.817871 | 0.331543 |
| `urn:mavedb:00000662-0-1#678` | BAP1 | s3_10 | 0.034668 | 0.0117188 | 0.718262 | 0.190039 |
| `urn:mavedb:00000662-0-1#7324` | BAP1 | s3_10 | 0.0952148 | 0.0429688 | 0.865723 | 0.311328 |
| `urn:mavedb:00000662-0-1#7968` | BAP1 | pm12 | 0.99979 | 0.968641 | 5.50781 | 3.06999 |
| `urn:mavedb:00000662-0-1#8628` | BAP1 | s3_10 | 0.0195312 | 0.03125 | 1.00488 | 0.251758 |
| `urn:mavedb:00000662-0-1#895` | BAP1 | s11_50 | 0.0078125 | 0.00390625 | 0.250732 | 0.0618652 |
| `urn:mavedb:00000662-0-1#928` | BAP1 | s11_50 | 0.0351562 | 0.0078125 | 0.11554 | 0.0660767 |
| `urn:mavedb:00000673-0-1#1166` | RAD51C | pm12 | 0.99599 | 0.959656 | 8.29688 | 3.61502 |
| `urn:mavedb:00000673-0-1#1396` | RAD51C | s3_10 | 0.0131836 | 0.0117188 | 0.0543213 | 0.0357666 |
| `urn:mavedb:00000673-0-1#1779` | RAD51C | s3_10 | 0.09375 | 0.0273438 | 0.515137 | 0.224121 |
| `urn:mavedb:00000673-0-1#1793` | RAD51C | s3_10 | 0.269531 | 0.445312 | 3.20703 | 1.35625 |
| `urn:mavedb:00000673-0-1#1798` | RAD51C | pm12 | 0.995216 | 0.949129 | 5.31641 | 3.00763 |
| `urn:mavedb:00000673-0-1#6409` | RAD51C | pm12 | 0.995224 | 0.944344 | 5.14062 | 2.96769 |
| `urn:mavedb:00000673-0-1#6411` | RAD51C | pm12 | 0.99556 | 0.944283 | 4.62891 | 2.86562 |
| `urn:mavedb:00000673-0-1#6417` | RAD51C | s3_10 | 0.0664062 | 0.140625 | 3.05664 | 0.818359 |
| `urn:mavedb:00000673-0-1#6462` | RAD51C | s11_50 | 0.0117188 | 0.00976562 | 0.136353 | 0.0487549 |
| `urn:mavedb:00000673-0-1#7090` | RAD51C | s3_10 | 0.0078125 | 0.0214844 | 0.429688 | 0.115234 |
| `urn:mavedb:00000673-0-1#7112` | RAD51C | pm12 | 0.996079 | 0.949188 | 9.05469 | 3.7562 |
| `urn:mavedb:00000673-0-1#7294` | RAD51C | pm12 | 0.992152 | 0.949179 | 8.53125 | 3.64758 |
| `urn:mavedb:00000673-0-1#7346` | RAD51C | s11_50 | 0.0078125 | 0.0117188 | 0.294434 | 0.078418 |
| `urn:mavedb:00000673-0-1#7637` | RAD51C | s11_50 | 0.015625 | 0.0214844 | 0.661133 | 0.169336 |
| `urn:mavedb:00000673-0-1#7700` | RAD51C | s3_10 | 0.09375 | 0.660156 | 3.07031 | 1.36797 |
| `urn:mavedb:00000673-0-1#7930` | RAD51C | s11_50 | 0.00952148 | 0.015625 | 0.310791 | 0.0873047 |
| `urn:mavedb:00000673-0-1#7972` | RAD51C | s11_50 | 0.0117188 | 0.0195312 | 0.117004 | 0.0546509 |
| `urn:mavedb:00000673-0-1#8314` | RAD51C | s11_50 | 0.03125 | 0.0117188 | 0.193604 | 0.0816895 |
| `urn:mavedb:00000675-a-1#1088` | VHL | s3_10 | 0.179688 | 0.698242 | 3.26953 | 1.53184 |
| `urn:mavedb:00000675-a-1#1098` | VHL | s3_10 | 0.160156 | 0.658203 | 3.15234 | 1.44883 |
| `urn:mavedb:00000675-a-1#1106` | VHL | s3_10 | 0.960693 | 0.968504 | 7.84375 | 3.49795 |
| `urn:mavedb:00000675-a-1#1110` | VHL | pm12 | 0.996072 | 0.976517 | 7.625 | 3.49759 |
| `urn:mavedb:00000675-a-1#1111` | VHL | pm12 | 0.996062 | 0.976535 | 7.98438 | 3.56947 |
| `urn:mavedb:00000675-a-1#1471` | VHL | pm12 | 0.995906 | 0.956946 | 6.85938 | 3.32473 |
| `urn:mavedb:00000675-a-1#1472` | VHL | pm12 | 0.995893 | 0.956841 | 5.86328 | 3.12539 |
| `urn:mavedb:00000675-a-1#1473` | VHL | pm12 | 0.996047 | 0.956923 | 6.22266 | 3.1975 |
| `urn:mavedb:00000675-a-1#1498` | VHL | s11_50 | 0.601562 | 0.144531 | 2.5918 | 1.26445 |
| `urn:mavedb:00000675-a-1#1509` | VHL | s11_50 | 0.0544434 | 0.0332031 | 0.678711 | 0.223389 |
| `urn:mavedb:00000675-a-1#1523` | VHL | s11_50 | 0.109375 | 0.0195312 | 0.557617 | 0.24043 |
| `urn:mavedb:00000675-a-1#1531` | VHL | s11_50 | 0.191406 | 0.03125 | 0.787598 | 0.380176 |
| `urn:mavedb:00000675-a-1#1535` | VHL | s11_50 | 0.0668945 | 0.0136719 | 0.356201 | 0.151807 |
| `urn:mavedb:00000675-a-1#1561` | VHL | s3_10 | 0.249512 | 0.0546875 | 1.07422 | 0.519043 |
| `urn:mavedb:00000675-a-1#1566` | VHL | s3_10 | 0.570312 | 0.179688 | 2.4375 | 1.2375 |
| `urn:mavedb:00000675-a-1#1571` | VHL | pm12 | 0.995457 | 0.976009 | 6.92969 | 3.3574 |
| `urn:mavedb:00000675-a-1#655` | VHL | s3_10 | 0.157227 | 0.174805 | 1.0498 | 0.541992 |
| `urn:mavedb:00000675-a-1#679` | VHL | s11_50 | 0.0449219 | 0.0332031 | 0.179321 | 0.113989 |
| `urn:mavedb:00001225-a-1#1747` | BRCA2 | pm12 | 0.976438 | 0.878433 | 6.45312 | 3.1455 |
| `urn:mavedb:00001225-a-1#1750` | BRCA2 | pm12 | 0.95459 | 0.877594 | 6.17188 | 3.06656 |
| `urn:mavedb:00001225-a-1#1796` | BRCA2 | s3_10 | 0.00585938 | 0.0429688 | 0.261963 | 0.101221 |
| `urn:mavedb:00001225-a-1#2884` | BRCA2 | s3_10 | 0.00585938 | 0.0195312 | 0.0583801 | 0.0370666 |
| `urn:mavedb:00001225-a-1#3387` | BRCA2 | s3_10 | 0.785156 | 0.873047 | 6.03516 | 2.86523 |
| `urn:mavedb:00001225-a-1#3900` | BRCA2 | s3_10 | 0.0253906 | 0.0273438 | 0.415283 | 0.135791 |
| `urn:mavedb:00001225-a-1#4955` | BRCA2 | s3_10 | 0.0240479 | 0.0390625 | 1.28809 | 0.320728 |
| `urn:mavedb:00001225-a-1#5954` | BRCA2 | pm12 | 0.99453 | 0.89798 | 6.50781 | 3.19407 |
| `urn:mavedb:00001225-a-1#5956` | BRCA2 | pm12 | 0.705078 | 0.856934 | 5.89844 | 2.7417 |
| `urn:mavedb:00001225-a-1#6002` | BRCA2 | s3_10 | 0.00585938 | 0.0195312 | 0.0526428 | 0.0359192 |
| `urn:mavedb:00001225-a-1#6784` | BRCA2 | pm12 | 0.991409 | 0.893593 | 3.41992 | 2.56899 |
| `urn:mavedb:00001225-a-1#6786` | BRCA2 | pm12 | 0.991926 | 0.890268 | 4.09766 | 2.70173 |
| `urn:mavedb:00001250-a-2#1529` | BARD1 | pm12 | 0.995917 | 0.948883 | 6.85547 | 3.31589 |
| `urn:mavedb:00001250-a-2#1557` | BARD1 | s3_10 | 0.0175781 | 0.00976562 | 0.581543 | 0.143652 |
| `urn:mavedb:00001250-a-2#2093` | BARD1 | s11_50 | 0.0751953 | 0.0253906 | 1.62598 | 0.425781 |
| `urn:mavedb:00001250-a-2#2142` | BARD1 | s11_50 | 0.00976562 | 0.00390625 | 0.217529 | 0.0571777 |
| `urn:mavedb:00001250-a-2#2921` | BARD1 | s11_50 | 0.00878906 | 0.0117188 | 0.459961 | 0.1125 |
| `urn:mavedb:00001250-a-2#2963` | BARD1 | s3_10 | 0.0410156 | 0.160156 | 3.4668 | 0.894531 |
| `urn:mavedb:00001250-a-2#3289` | BARD1 | pm12 | 0.987663 | 0.932682 | 4.21484 | 2.76331 |
| `urn:mavedb:00001250-a-2#3362` | BARD1 | s11_50 | 0.0078125 | 0.0195312 | 0.198608 | 0.0670654 |
| `urn:mavedb:00001250-a-2#3447` | BARD1 | s11_50 | 0.0078125 | 0.0117188 | 0.0365601 | 0.0268433 |
| `urn:mavedb:00001250-a-2#3534` | BARD1 | pm12 | 0.994926 | 0.956253 | 7.60156 | 3.47149 |
| `urn:mavedb:00001250-a-2#4254` | BARD1 | s11_50 | 0.0102539 | 0.015625 | 1.13965 | 0.253809 |
| `urn:mavedb:00001250-a-2#7374` | BARD1 | s3_10 | 0.0234375 | 0.0195312 | 0.274658 | 0.0979004 |
| `urn:mavedb:00001250-a-2#7375` | BARD1 | s3_10 | 0.0488281 | 0.0507812 | 0.90625 | 0.280859 |
| `urn:mavedb:00001250-a-2#7904` | BARD1 | pm12 | 0.988081 | 0.905678 | 5.54297 | 3.00235 |
| `urn:mavedb:00001250-a-2#7906` | BARD1 | pm12 | 0.987793 | 0.90522 | 5.48438 | 2.98989 |
| `urn:mavedb:00001250-a-2#7990` | BARD1 | s3_10 | 0.0078125 | 0.0078125 | 0.170532 | 0.0497314 |
| `urn:mavedb:00001250-a-2#8013` | BARD1 | pm12 | 0.876465 | 0.936134 | 6.17578 | 3.04776 |
| `urn:mavedb:00001250-a-2#8320` | BARD1 | s3_10 | 0.0449219 | 0.0195312 | 0.313721 | 0.127197 |
| `urn:mavedb:00001259-a-2#1181` | PALB2 | s3_10 | 0.0571289 | 0.0117188 | 0.616699 | 0.192188 |
| `urn:mavedb:00001259-a-2#1544` | PALB2 | pm12 | 0.996081 | 0.933571 | 6.73047 | 3.27575 |
| `urn:mavedb:00001259-a-2#1941` | PALB2 | s11_50 | 0.0078125 | 0.0078125 | 0.155273 | 0.0466797 |
| `urn:mavedb:00001259-a-2#1956` | PALB2 | s11_50 | 0.00488281 | 0.00976562 | 0.0699463 | 0.0286377 |
| `urn:mavedb:00001259-a-2#2454` | PALB2 | pm12 | 0.991173 | 0.932495 | 4.50391 | 2.82445 |
| `urn:mavedb:00001259-a-2#2502` | PALB2 | s11_50 | 0.0117188 | 0.0078125 | 0.327881 | 0.0851074 |
| `urn:mavedb:00001259-a-2#2825` | PALB2 | s3_10 | 0.0114136 | 0.0273438 | 1.69141 | 0.377039 |
| `urn:mavedb:00001259-a-2#3404` | PALB2 | s11_50 | 0.0283813 | 0.0898438 | 0.674805 | 0.253186 |
| `urn:mavedb:00001259-a-2#3765` | PALB2 | pm12 | 0.947754 | 0.883347 | 3.92383 | 2.61587 |
| `urn:mavedb:00001259-a-2#5222` | PALB2 | pm12 | 0.992118 | 0.929629 | 5.58203 | 3.03815 |
| `urn:mavedb:00001259-a-2#5272` | PALB2 | s11_50 | 0.00390625 | 0.0117188 | 0.0317993 | 0.0219849 |
| `urn:mavedb:00001259-a-2#5367` | PALB2 | s3_10 | 0.00585938 | 0.0351562 | 0.0501709 | 0.0510498 |
| `urn:mavedb:00001259-a-2#5383` | PALB2 | pm12 | 0.995906 | 0.937164 | 4.0625 | 2.74557 |
| `urn:mavedb:00001259-a-2#718` | PALB2 | s3_10 | 0.260864 | 0.167969 | 1.88281 | 0.805395 |
| `urn:mavedb:00001259-a-2#9246` | PALB2 | pm12 | 0.99437 | 0.928093 | 3.49023 | 2.62051 |
| `urn:mavedb:00001259-a-2#9630` | PALB2 | s3_10 | 0.0703125 | 0.5 | 2.25391 | 1.02109 |
| `urn:mavedb:00001259-a-2#9896` | PALB2 | s11_50 | 0.00995588 | 0.0390625 | 0.148438 | 0.0787059 |
| `urn:mavedb:00001259-a-2#9969` | PALB2 | s3_10 | 0.459961 | 0.650391 | 2.53711 | 1.61777 |
