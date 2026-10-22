# Evidence-strength rework: change log by round

A dated log, kept by hand. It records what changed between rounds and why; the
numbers in it are as they stood when each round was written. The current value of
every quantity is in `phase1/reports/evidence/`, and `docs/evidence-delta.md` (written
by `src/evid_delta.py`) is the generated comparison with the published analysis set.

## Round 2 (2026-09-20): what moved relative to the first round's report

A list, no interpretation.

### Decisions that changed which number is quoted

| | first round | now |
|---|---|---|
| the pool every "overall" number is read off | `all_1_50` | `s3_50` (3-50 bp, pm12 excluded) |
| `all_1_50` in the tables | the pooled row | kept, column labelled `out_of_scope_incl_pm12` |
| the 588 ClinVar records with no clinical assertion | inside `recorded_unclassified`, not visible | still inside it; `clinvar_arm_strict` breaks them out beside the three arms, and `clinvar_arm_no_assertion.csv` is the crosstab |
| tools in the panel | 8 + spliceai_walker | + `avi`, `avi_splice_sites`, `avi_splice_site_usage`, `avi_splice_junctions` |

### Numbers that moved because the pool changed

| quantity | `all_1_50` | `s3_50` |
|---|---|---|
| Walker cut point, LR+ (all arms, spliceai_walker) | 9.55 | 8.28 |
| Walker cut point, LR(<=0.1) | 0.126 (Moderate) | 0.252 (Supporting) |
| sensitivity at >= 0.2 | 0.849 | 0.696 |
| ClinVar arms on AUROC: tools where `classified` is highest | 13/13 | 3/13 |
| the same, BRCA1 dropped | 13/13 | 2/13 |
| ClinVar arms on Walker LR+: tools where `classified` is highest | 13/13 | 11/13 |

### New quantities

- `avi` against functional pathogenicity: rho 0.393 (0.424 at 3-10 bp, 0.098 at 11-50 bp);
  `avi_splice_sites` 0.405. Against `alphagenome` (v0.7.0): 0.729 and 0.872.
- Strong under leave-one-gene-out at 3-50 bp, folds whose held-out ratio clears 18.7:
  avi 7/7, alphagenome 6/7, pangolin 6/7, spliceai 5/7, spliceai_walker 5/7.
- Strong under leave-one-gene-out at 11-50 bp: 0/6 for every tool, including the three
  that reach it in-sample.
- In-frame share among the variants a tool calls and the assay does not, 3-10 bp,
  SpliceAI: seven genes 40.1%, DDX3X 26.8%, TP53 9.5%.
- Pangolin on DDX3X under E3's leave-one-gene-out thresholds: Moderate in five of six
  evaluable cells, Strong in one (11-50 bp, FDR labels, LR 29.6).

### Corrections carried from the first round's internal review

Listed in that round's report; the tables here are the corrected ones. The only
correction whose direction is worth restating: the local likelihood ratio's numerator
is no longer floored, so BP4 Strong and Very strong are no longer unreachable by
construction, and BP4 Moderate is reachable at 3-10 bp.

### Correction to the round-2 entry above (2026-09-23)

The row "ClinVar arms on Walker LR+: tools where `classified` is highest" compared
arms by applying SpliceAI's 0.2 cut point to every column's raw score. That cut
point is defined on SpliceAI's delta score only; for AlphaGenome (range about 0.8
to 2.2), CADD or phyloP it selects nearly every variant and the ratio means
nothing. The comparison is withdrawn. It is replaced by `arms_at_tool_threshold.csv`,
which compares arms at each column's own fitted threshold, and at the published cut
point for the two SpliceAI columns only.

## Round 3 (2026-09-23)

- ClinVar arms compared at each column's own fitted threshold
  (`src/evid_arms.py`, `arms_at_tool_threshold.csv`); the direction differs by
  distance band, which the earlier pooled comparison could not show.
- The external gene DDX3X now carries the four AlphaGenome Atlas columns, merged
  from the same Atlas pass as the analysis set; the DDX3X merge and both external
  in-frame attributions are stages of `scripts/reproduce_evidence.py`.
- New stages: training-signal table (E9), fusion coefficients across folds (E10),
  thresholds refitted on every subset of training genes (E11), printed main tables
  (E12), printed supplementary tables (E13), figures (E14).
- The DDX3X assay was earlier described in a docstring as a different kind of
  readout. It is not: like most of the seven, it scores depletion from HAP1 cells.
- Printed thresholds keep as many figures as it takes to select the same variants as
  the fitted value; AlphaGenome's sit just under its ceiling of 2.2.
- The companion atlas is located through `evid_common.ATLAS_REPO`, defaulting to a
  checkout beside this repository.

## Round 4 (2026-09-23)

Changes that alter results, each with the reason:

- **DDX3X labels.** The deposit does publish a per-variant classification, in a
  separate score set (urn:mavedb:00000658-0-1), a random-forest call on the assay's
  fold-changes; it is now the primary label, as the seven genes use their deposits'
  own classifications. The control-anchored label is kept as a sensitivity label,
  ungated and gated: in score set q-1 (the final exon) the nonsense controls are not
  depleted (control AUROC 0.51), and that set supplied 37 of the label's 71 in-scope
  damaging calls. The earlier finding that the Strong thresholds did not transfer to
  DDX3X came from that set.
- **External intervals.** Every external likelihood ratio now carries a
  variant-level bootstrap interval within the gene and the damaging and normal counts
  in its band, and a tier is reported from the lower bound as well as the point
  estimate (`evid_external._band_detail`).
- **TP53 to |offset| 12.** The archived TP53 table stopped at 8 because of the
  earlier study's splice-window constant; `src/evid_tp53_extend.py` scores the 96
  deposited SNVs at offsets 9-12 with the same scorers and carries the 192 archived
  rows unchanged. It runs offline from cached scores as a stage of the entry point.
- **Fusion without the AlphaGenome splice score.** The AlphaGenome terms of service
  bar the use of its outputs to train other models (LICENSE-DATA), so the elastic net
  is refitted on the other seven panel columns (`evid_common.FUSION_FEATURES`). The
  single-column thresholds were rerun with it and came back identical.
- **Atlas combined score provenance.** It is a supervised model trained on gnomAD
  allele frequency, with AlphaMissense and conservation inputs, whose checkpoint was
  selected on four saturation genome editing datasets including the BRCA1, RAD51C and
  DDX3X assays used here. Its provenance row is corrected, and Table 3 counts its
  held-out folds over the other five genes.
- **Walker configuration.** `config/walker2023.yaml` now records the applied weight
  (Supporting, p. 1056) beside the calibrated strength (Moderate, p. 1051), and the
  BRCA1 check at a score of 0.5 (p. 1055); an earlier note told readers not to carry
  "supporting", which was wrong.
- New outputs: `depth_counts.csv` (labelled and damaging variants by intronic depth),
  terms-of-use columns in `predictor_training_provenance.csv` and Table 2, and an
  in-scope split in `clinvar_arm_no_assertion.csv`.

## Round 5 (2026-09-29): reproducing from public sources

A clean-room run from public sources only (a fresh clone, the atlas release archive
from Zenodo, a Python environment built from `requirements-evidence.lock.txt`, the
ClinVar release from NCBI) regenerated 160 of 163 committed outputs byte for byte. The
other three depended on the machine that wrote them, and every result was unchanged.

- **TP53 REF check.** `evid_tp53_extend` checked the 288 REF bases against a 706 MB
  local FASTA that no archive carries, so elsewhere it recorded "not checked" and a
  local path. It now reads TP53's GRCh38 region (MANE exons +/- 5,000 nt, the bounds
  rule used for DDX3X; 29,070 nt from Ensembl REST) from the repository's reference
  cache, where it is registered with its sha256. The fetched sequence is identical to
  the Ensembl r112 FASTA the check used before, and all 288 REF bases match it.
- **Atlas commit.** The analysis-set manifest read the atlas commit with `git`, which
  gave "unknown" for the release archive. An archive now records the commit its
  release was cut from (every file tracked at tag v2.5.0-submission is byte-identical
  in the archive), a checkout records its own HEAD, and an archive unpacked inside
  another checkout is no longer mistaken for it. The manifest also records the sha256
  of the atlas's AlphaGenome v0.6.1 score file, which the set build reads, and a
  missing copy of that file now stops the build instead of dropping the column.
- **Entry point.** `--fetch-inputs` downloads the ClinVar release (NCBI's archive
  first, the weekly directory as fallback) and the atlas release archive, checks
  their checksums and unpacks the archive under `phase1/data/evidence/companion_atlas/`
  (gitignored). Before any stage runs, the entry point checks every atlas file the
  stages read, because the atlas's `results/` files ship in its archive and not in
  its git repository.
- **One command.** `bash reproduce.sh` builds `.venv` from the lock file (and refuses
  a `.venv` that differs from it), then runs the entry point with `--fetch-inputs
  --verify`. `--verify` records every output the checkout carries before the run,
  compares each with what the run writes, and ends with `REPRODUCED` or the list of
  what differs; a manifest may differ only in its wall-clock fields.
- **Machine independence.** The supplementary-table manifest names atlas sources
  `atlas:...` wherever the atlas sits (an archive unpacked inside this repository was
  named by its repository path). `reproduce.sh` runs BLAS single-threaded, so sums
  run in one order on any number of cores (the outputs are identical this way and
  with default threads), and ignores any user matplotlib configuration.
  `.github/workflows/reproduce.yml` runs `bash reproduce.sh` on a GitHub macOS arm64
  runner, where it must pass, and on Linux for information.
- **Downloads that survive the network.** `--fetch-inputs` resumes a download cut
  short or stalled for a minute, asking the server for the remaining bytes (NCBI and
  Zenodo both honour HTTP ranges), and checks the announced size, because a
  connection closed early otherwise ends the read without an error. The checksums
  are checked after every download as before.
- **Quick check.** `bash reproduce.sh --check` (or `scripts/reproduce_evidence.py
  --check`) needs no data download: it checks every checksum the manifests and
  provenance records carry about files in the repository (row-order-free content
  hashes where the record says so; the ClinVar release and the atlas when they are
  here) and runs the test suite. The one test that reads the atlas's
  `predictor_resources.py` now skips when the atlas is absent, as the other
  atlas-dependent tests do, so the suite passes on a bare clone.
- **A third ClinVar source.** The release's Zenodo record (10.5281/zenodo.23043214)
  carries NCBI's `clinvar_20260615.vcf.gz` unchanged (NCBI's ClinVar data are in the
  public domain), and `--fetch-inputs` tries it after NCBI's archive and weekly
  directories, checking it by sha256. On 2026-09-29 the file came from NCBI's
  `archive_2.0/2026/` directory.
- **Review of the one command, before its first CI run.** An adversarial review
  found 32 defects; all are fixed, with a test each where one can hold it:
  - `--check` failed on every fresh clone because it counted the undownloaded ClinVar
    file as missing; it now checks ClinVar and the atlas only when they are here, and
    says so. It also left out whole records: it now covers the provenance records
    under `data/` as well as `phase1/data/`, the TP53 manifests' input hashes, the
    frozen matrix's rescored columns and the supplementary tables' sources.
  - `--verify` took its baseline from the working tree, so a rerun after an
    interrupted run compared the outputs with themselves; it now compares with HEAD
    in a clone and with the tree as first unpacked (`.reproduce_baseline/`) in an
    archive. It passed an output no stage writes any more; such a file now fails as
    "not rewritten" unless `scripts/reproduce_inputs.txt` lists it as an input. A
    manifest could differ in more than its build time (1 for 1.0, reordered keys);
    now only the wall-clock values may differ. A parquet holding missing values
    crashed the comparison.
  - With `--verify` the atlas must be the release file for file
    (`scripts/atlas_release_v2.5.0.sha256`, 46 files) and, for a checkout, at the
    release commit; a local checkout that is not is passed over for the archive. The
    stages, the tests and the entry point now find the atlas by one rule
    (`evid_common.resolve_atlas`), and a relative `EVID_ATLAS_REPO` is made absolute.
  - A ClinVar source whose file failed its checksums ended the fetch; the next
    source is now tried. A 404 was retried for eleven minutes; client errors now end
    a download at once, and a failed atlas download ends with a message, not a
    traceback.
  - The three tests that list files with `git ls-files` fall back to walking the
    tree in a release archive (`tests/_checkout.py`), and the tests that make git
    commits ignore the user's git configuration.
  - `reproduce.yml` skips documentation-only pushes but not the two Markdown files
    the pipeline writes or reads; `reproduce.sh` removes its temporary matplotlib
    folder, and `--verify` its temporary copies.
- **Known, not changed:** `requirements-evidence.lock.txt` pins `git-filter-repo`, a
  development tool the analysis does not use; changing the lock file changes the
  definition of the exact environment, so it stays until a release that needs to.
- **Two verdicts, from the first independent runs.** On a GitHub-hosted M1 (macOS 14)
  every printed table and figure matched byte for byte, while eight intermediate
  tables of the elastic-net combination matched only to 1e-9 in relative terms:
  numpy and scipy call the operating system's Accelerate library on macOS, whose
  kernels differ between chips and releases. Pinning the runner to Python 3.12.13 was
  not possible (the macOS 14 arm64 runners offer 3.12.10 at most). `--verify` now
  ends `REPRODUCED` (byte for byte) or `REPRODUCED NUMERICALLY` (every printed table
  and figure byte for byte, other files within 1e-9, and manifests differing only in
  the hashes of such files); anything else fails. The macOS job must reach one of the
  two. On the Linux runner (x86_64) the TP53 stage first stopped at its exact
  reproduction check, one Nucleotide Transformer score reading back 4.4e-16 apart
  through pandas' default CSV parser, which is not correctly rounded and can differ in
  the last binary digit across CPU architectures; that check now allows a difference
  below 1e-12 and still records the exact counts.
- **What the Linux runner showed.** On x86_64 Linux every number agreed with the
  committed one to a relative 1e-9. Two printed supplementary tables (S6 and S8)
  showed the AlphaGenome splice score's 3–10 bp thresholds to fewer significant
  figures (2.199883 against 2.19988276958): the printer adds figures until the printed
  value selects the same variants as the fitted one, and that fitted threshold lies
  within floating-point precision of an observed score, so a last-digit difference
  decides how many figures it takes. The figures differed in their bytes where Linux
  lacks the Arial font. `--verify` now also reads a text report number by number.

## Release 3.0.0 (2026-09-30): `v3.0.0-submission`

The evidence-strength study with its one-command reproduction, archived at Zenodo as
version DOI 10.5281/zenodo.23043214 (all versions: 10.5281/zenodo.22674887). Nothing
in it changes a number of the study; the rounds above record what was reworked.

- **What changed since 2.0.0:** the evidence-strength analysis (rounds 2 to 5 above)
  and `bash reproduce.sh`, which fetches the two public inputs, rebuilds every output
  and checks it against the committed one, and its five-minute `--check`.
- **Where it reproduces:** byte for byte on the machine the outputs were made on
  (macOS arm64, Python 3.12.13), the analysis-set manifest's build time aside. On a
  GitHub-hosted M1 (macOS 14) every printed table and figure byte for byte, and the
  elastic-net combination's intermediate tables to a relative 1e-9
  (`REPRODUCED NUMERICALLY`). On x86_64 Linux every computed value to 1e-9; the
  figures differ in their bytes without the Arial font, and S6 and S8 print two
  AlphaGenome thresholds to fewer significant figures (see above).
- **The ClinVar mirror:** the Zenodo record carries NCBI's `clinvar_20260615.vcf.gz`
  unchanged, the third source `--fetch-inputs` tries.
- **Title:** `CITATION.cff` now names the release for the evidence-strength study
  ("variant-fm-benchmark: evidence strength of splice-region variant effect predictors
  against saturation genome editing functional assays"); the earlier title described
  the calibration and fusion study of 2.0.0.

## After 3.0.0

- **The guardrail suite on Linux.** With `main` fast-forwarded to the release, the
  pytest workflow ran on Ubuntu (x86_64) and three tests failed that compared a
  rebuild with the committed files byte for byte: the elastic-net refit against E3's
  score grid and against its coefficient table, and the TP53 table rebuilt from its
  cached scores. Both reach the last binary digit through the machine (the operating
  system's maths library; pandas' default CSV parser). The tests now require two
  rebuilds on one machine to be byte-identical, the printed S9 table to match byte
  for byte, and the rest to agree with the committed files to the relative 1e-9 of
  `--verify`'s numerical tier. The full run still checks every file byte for byte on
  the machine that wrote it. Release 3.0.0 is unchanged; on x86_64 Linux its three
  strict tests fail as described here.

## Release 3.0.1 (2026-09-30): `v3.0.1-submission`

Archived at Zenodo as version DOI 10.5281/zenodo.23051009 (all versions:
10.5281/zenodo.22674887).

- **The cross-machine rebuild tests.** The three tests described under "After 3.0.0"
  are archived with the release: two rebuilds on one machine byte-identical, the
  printed S9 table byte for byte, and the rest against the committed files to the
  relative 1e-9 of `--verify`'s numerical tier. On x86_64 Linux the suite passes: 368
  passed, 14 skipped, 0 failed (run 36664902183). In the 3.0.0 archive these three
  tests fail on x86_64 Linux, as described above; that record stays as it is.
- **The ClinVar mirror** now points at this release's record, which holds the same
  file (NCBI's `clinvar_20260615.vcf.gz`, unchanged); the 3.0.0 record keeps its copy.
- **Nothing else:** the pipeline code and every output are those of 3.0.0.

## After 3.0.1: practical metrics (stage E15, 2026-09-30)

- **Why.** A reader of the manuscript pointed out that the likelihood ratio compares
  predictors well but does not tell a user what to expect in classical terms,
  sensitivity and specificity at a threshold, or across thresholds (ROC). The study
  reported AUROC only in the ClinVar-arm comparison and no sensitivity or specificity
  at any fitted threshold.
- **What.** A new stage, `src/evid_practical_metrics.py` ("E15 practical metrics",
  stage 22, after E14 and before `evid_delta`), reads the thresholds the study already
  uses and fits none: Walker's cut points, and E3's thresholds from Supporting to Very
  strong on the PP3 side and Supporting and Moderate on the BP4 side, in-sample, and
  held out on the PP3 side (E3 keeps per-fold thresholds for that side only). At each
  it gives sensitivity and specificity with Wilson 95% intervals, coverage (the share
  of all scored variants in the band that receive the evidence, and the same among
  labelled variants), prevalence, the observed PPV, the band likelihood ratio, the
  posterior at a prior of 0.10, and the tier by E3's interval rule beside the tier by
  the band ratio. Held-out rows apply each fold's threshold to its held-out gene and
  add the counts over the folds that reached the tier. It also writes ROC and
  precision-recall curves with the seven genes pooled, at every distinct score; their
  AUROC and AUPRC beside the gene-pooled AUROC the study reports, with the difference
  split into same-gene pairs, cross-gene pairs and the genes that pool leaves out; the
  same metrics for DDX3X and TP53 at the thresholds E7 carried over; a two-row figure;
  and a one-page summary. The columns are E3's thirteen plus the fusion;
  alphagenome_v061 has no E3 threshold and is left out, as in E2-E14.
- **New files.** Under `phase1/reports/evidence/`: `practical_metrics_by_threshold.csv`,
  `practical_metrics_logo_by_gene.csv`, `practical_metrics_external.csv`,
  `practical_roc_curves.csv`, `practical_pr_curves.csv`, `practical_curve_auc.csv`,
  `figures/figure_practical.pdf` and `.png`; and `docs/practical-metrics-summary.md`,
  which `--verify` now compares as well. Tests: `phase1/tests/test_evid_practical_metrics.py`.
- **Reading thresholds back.** The stage reads every table with pandas'
  `float_precision="round_trip"`: the default parser can land one binary digit away
  from a written threshold, and a threshold is an observed score, so the variants on it
  can fall out of a band that includes its boundary.
- **Unchanged:** every output of E1-E14, byte for byte. The manuscript and the
  supplement do not use these numbers.

## After 3.0.1: thresholds read back exactly, and two printing rules (2026-09-30)

- **Cause.** Nine of the eleven stages that read a threshold back from a CSV used
  pandas' default float parser, which is not correctly rounded: it can round the last
  binary digit the wrong way, and it keeps only sixteen decimal places of a value below
  one, so a small value loses its last significant digits (`float_precision="high"` is
  the same parser). A threshold is an observed score, so a threshold read back high
  leaves the variants that sit exactly on it out of a band that includes its boundary
  (score >= t), and one read back low leaves them out of score <= t. `evid_arms` read
  the combined Atlas score's Moderate threshold at 3-10 bp, written as
  1.0710335969924927, one unit in the last place high, and the damaging PALB2 variant
  urn:mavedb:00001259-a-2#9971, which scores exactly the threshold (ClinVar arm:
  recorded, unclassified), fell out of the band.
- **The four rows of Supplementary Table S10 whose counts change** (combined Atlas
  score, 3-10 bp, Moderate threshold; variants in band, sensitivity and band likelihood
  ratio, 12 cells):

  | rows | before | after |
  |---|---|---|
  | seven genes, all variants | 340, 0.516, 28.1 | 341, 0.518, 28.2 |
  | seven genes, recorded, unclassified | 119, 0.527, 24.2 | 120, 0.532, 24.4 |
  | without BRCA1, all variants | 262, 0.476, 25.2 | 263, 0.479, 25.4 |
  | without BRCA1, recorded, unclassified | 75, 0.486, 18.9 | 76, 0.493, 19.2 |

  The arms' ranking is unchanged (classified highest), so the printed Table 4 (4a-4c)
  and the text that describes it are unchanged. Its supporting file,
  `tables/table4_detail_fitted_threshold_by_arm.csv`, carries the two
  recorded-unclassified ratios above (24.18 to 24.41 and 18.93 to 19.22).
- **Fix.** `evid_common.read_back` reads a table the pipeline wrote with
  `float_precision="round_trip"`, and every stage that reads a threshold back from a CSV
  reads it through it: `evid_arms`, `evid_external`, `evid_territory_metrics`,
  `evid_diagnostics`, `evid_tier_logo`, `evid_fig_data`, `evid_tables`,
  `evid_supp_tables` and `evid_figures`, and `evid_fusion_stability` and
  `evid_practical_metrics`, which already read exactly and now share the function (that
  alone changes nothing they write); the tests that read thresholds back do the same.
  `phase1/tests/test_evid_read_back.py` lists every default-parser read that remains,
  call by call, with the reason none of them uses a threshold, and fails on any other;
  checks that every printed threshold in S6, S7, S8 and S10 selects exactly the
  analysis-set variants its fitted value selects, and every one in S12 exactly the
  external gene's variants its carried value selects (Table 3b rounds its thresholds
  outward by design, and S9's per-fold fusion thresholds are not checked); checks that
  stored thresholds and their copies read back as written; and holds the PALB2 variant
  in its band.
- **Printed thresholds that moved with the fix.** The supplementary tables print a
  threshold with four significant figures, or more where four would move a variant
  across it, and had chosen those digits against the misread values. 39 printed
  thresholds change: S10 now prints this threshold as 1.071 (not 1.071034) in all eight
  of its rows, S6 changes 3 cells in 2 rows, S7 22 cells in 12 rows and S8 6 cells in 6
  rows. Of the 936 thresholds printed in S6, S7, S8 and S10, these 39 and the S8 cell in
  (ii) below were the only ones that selected a different set of variants from their
  fitted values: each left out the one variant that scores exactly its fitted value.
- **Two printing rules, corrected at the same time** (the user's decisions, after review
  found them). (i) A ratio of counts that lies exactly on a decimal tie can come out of
  the arithmetic one unit in the last place below it (243/4 = 60.75 is stored as
  60.74999999999999), and a bootstrap bound that numpy's percentile interpolates between
  two such ratios further below it (exactly 441/20 = 22.05, stored 13 units low; exactly
  259/20 = 12.95, stored 82 units low); half-up rounding of the stored value rounded it
  down. `fmt_lr` now rounds to twelve significant figures first, which puts each such
  value back on its tie; of the 2,955 values the supplementary tables format this way,
  the nearest to a tie that is not one is 2.3e-6 away (relative). Four printed ratios
  change: S6, combined Atlas score, 11-50 bp, Moderate, the held-out maximum and the
  RAD51C fold, 60.7 to 60.8 (243/4); S12, TP53 under the median-split label, 3-10 bp,
  Atlas splice junctions, the Supporting ratio's upper 95% bound, 22.0 to 22.1; and S12,
  TP53 under the mid-band-excluded label, 3-10 bp, SpliceAI (published basis), the upper
  95% bound of the ratio at the published cut point, 12.9 to 13.0. Two ratios of 14.55
  (291/20) in S6 and S8 print 14.6 as before: the default parser had read them back one
  unit up, and the exact read-back alone would have printed 14.5. (ii) `fmt_thr`
  searched four to eleven significant figures and then returned twelve unchecked; it now
  searches to seventeen. One S8 median threshold changes, 0.427076786757 to
  0.4270767867565, which now selects the BAP1 variant that scores it exactly.
- **Copies.** Tables that copy a value from a threshold table now hold it as written:
  345 cells in 9 files move by 1 to 11 units in the last place (at most 1.7e-15
  relative): `arms_at_tool_threshold.csv`, `territory_metrics.csv`,
  `tier_logo_table.csv`, `tables/table4_detail_fitted_threshold_by_arm.csv`,
  `fig_data/fig1_walker_cutpoints.csv`, `fig_data/fig4_external.csv`,
  `external_ddx3x.csv`, `external_tp53.csv` and `practical_metrics_external.csv`. Beyond
  these, the files change only in the S10 rows above (12 cells of
  `arms_at_tool_threshold.csv`: count, sensitivity and ratio in four rows) and the two
  Table 4 detail ratios.
- **Left as they are** (the user's decision: recorded here, not changed). Reads that use
  no threshold still use the default parser, several of them of tables that have a
  threshold column they leave alone (`territory_metrics.csv` for fig3's AUROCs,
  `threshold_dilution_summary.csv` for the shares in Supplementary Figure S2, and others
  read for other columns or a row count), so some copied values keep the default
  parser's digits: in fig1, 27 of the 72 ratios in its twelve pooled seven-gene rows (6
  of 36 at the PP3 cut point, 21 of 36 at the BP4 cut point, where a value below one
  loses its seventeenth significant digit; up to 7 units in the last place) and 7 of the
  18 in its nine rows without BRCA1; 79 of the 560 AUROC and PR-AUC values in
  `fig_data/fig3_evidence_by_territory.csv` (1 or 2 units), and 26 of the 168 that
  `practical_curve_auc.csv` copies from it; and in `fig_data/fig2_interval_lr.csv`
  59,632 of 138,661 scores (up to 8.4e-13 relative; the PALB2 variant's score,
  1.0710335969924927, is held there as 1.071033596992493) and 104,023 of 415,983 ratios.
  No printed table reads the three fig_data files, and the AUROCs that the
  practical-metrics page and its slide figure print from fig3 are the same either way.
- **Scope.** E4.1 asked that only the S10 files change. With the complete fix the user
  chose, 15 outputs changed: S6, S7, S8, S10 and S12, `supplement_tables_manifest.csv`,
  the Table 4 detail file, `arms_at_tool_threshold.csv` and the seven other copies
  above. No other output changed. The manuscript's text cites none of the changed
  values.
- **Review.** An independent review of the first commit found these, now fixed: the
  allow-list test skipped every line that also called read_back, and so missed the
  default read beside one in `evid_diagnostics` (it reads Walker's ratios, not a
  threshold); it now finds reads in the syntax tree, call by call, and three of its
  reasons are corrected (the read in `evid_tp53_extend` keeps the default parser on
  purpose, so that new TP53 scores are stored the way the frozen ones were); the S12
  check compared printed thresholds on the seven genes' scores only; one case of the tie
  test (243/4, exact in binary) passed on the old code as well, and none covered the
  bootstrap bounds; it now asserts on the stored values; rounding to fifteen significant
  figures reaches a tie only from within 5e-14 at these magnitudes (14 to 28 units in
  the last place), so the 12.95 above, 82 units low, still printed as 12.9; and the
  slide figure's overlap check measured the key at 100 dpi, while the PNG is written at
  150, where DejaVu Sans, the fallback without Arial, overlapped the inset's tick labels
  by a point.

## After 3.0.1: the practical metrics for slides (2026-09-30)

- `figures/figure_practical_slides.png`: the practical-metrics curves at 16:9 (2000 x
  1125 pixels) with larger type and the inset kept (shorter than in print, so that the
  key below it leaves its tick labels clear in Arial and in DejaVu Sans; the stage
  measures this at the resolution it saves and stops if they touch); each ROC panel's
  key gives the manuscript's AUROC (per-gene AUROCs combined on the logit scale) instead
  of the pooled curve's, and the caption on the figure says the curves pool the seven
  genes and that each AUROC is the manuscript's, not the area under these curves (the
  two can rank the predictors differently; at 11-50 bp the manuscript's AUROC also
  leaves out BRCA1, with fewer than ten damaging variants there).
- `docs/practical-metrics-summary.md`: coverage is now the share of labelled variants,
  the denominator of the manuscript's shares ("14.8% of the variants" at 0.2), which
  changes five printed cells (SpliceAI 3-10 bp Supporting 0.15 (0.16) to 0.15 (0.17);
  Pangolin 3-10 bp Moderate 0.12 (0.13) to 0.12 (0.14) and Strong 0.07 (0.07) to 0.08
  (0.08); combined Atlas score 3-10 bp Supporting 0.14 (0.17) to 0.15 (0.17) and 11-50
  bp Supporting 0.05 (0.04) to 0.04 (0.04)); a column gives the manuscript's AUROC
  beside the pooled curve's, and its note says which genes it leaves out; a table gives
  the four predictors' sensitivity at a specificity of 0.95 in the three bands; the
  first sentence quotes the manuscript's AUROC and says it is not the area under the
  curves (0.85 against 0.91, where it quoted the pooled 0.81 against 0.90), and the
  third says "held out, ... no more than 8.1% of the labelled variants" (was 8% of all
  variants).
- **Known improvement**, recorded on the summary page: the BP4 side has no held-out
  sensitivity, specificity or coverage, because E3 keeps only the median held-out
  likelihood ratio of its per-fold BP4 thresholds (Supplementary Table S7) and writes
  per-fold thresholds for the PP3 side only; held-out BP4 sensitivity and specificity
  need E3 to store its per-fold BP4 thresholds as well.

## After 3.0.1: four further external genes, not yet scored (stage E16, 2026-10-20)

Why: the external genes go from two (DDX3X, TP53) to six. CTCF, SFPQ, RAD51D and XRCC2
have HAP1 saturation genome editing deposits in MaveDB (CC0), scored for cell fitness
like the seven genes; they are added for a purely prospective test. This round only
records them: nothing is fitted, the seven genes' thresholds as archived with 3.0.1
are untouched, and no predictor is scored on them yet (no score column, likelihood
ratio, ROC or figure).

- Stage 23, `src/evid_newgenes.py`, offline, between E15 and `evid_delta` (now 24):
  each deposit's record (`newgenes_record_check.csv`, `docs/newgenes-record-check.md`,
  every fact with the file it was read from and when that file was fetched), a variant
  table per gene (`<gene>_variants.parquet` and its provenance: GRCh38 position,
  signed offset, side, stratum, score, range class, label), counts by offset and by
  stratum (`newgenes_counts_by_offset.csv`, `newgenes_label_summary.csv`) and three
  check tables (`newgenes_coordinate_check.csv`, `newgenes_offset_sample.csv`,
  `newgenes_controls.csv`, whose synonymous and nonsense counts keep SNVs apart from
  3-bp deletions that MaveDB projects to a stop). For CTCF and SFPQ the record is
  checked against the author's of 2026-10-19 (URN, licence, publication status,
  first-published month); the stage stops on any difference, and on ranges that
  overlap beyond a shared end point (E16.10).
- Inputs, fetched once by `scripts/stage_E16/` (one script per purpose, each writing
  a run record under `scripts/stage_E16/run_log/`): each score set's metadata, its
  superseded versions, its experiment, the score table as served, the
  investigator-provided ranges, and MaveDB's variant-details stream, under
  `phase1/data/evidence/external/<gene>_raw/` with a manifest (path, bytes, sha256,
  source URL, fetch time, API version); and each gene's transcript model and Ensembl
  sequence slice, in the reference cache. RAD51D's and XRCC2's URNs were found with
  MaveDB's search; the requests, candidates and reasons are in the run records.
- Labels: the deposit's own investigator-provided ranges, applied mechanically:
  abnormal 1, normal 0, Indeterminate NA, as `config.ASSAY_LABELS` treats intermediate
  classes; a score outside every range, or on a point two ranges share, is
  `unassigned`. CTCF's Indeterminate and Functionally Normal ranges share one end point;
  the author's rule (2026-10-20) makes a score there `unassigned`. No intron-side SNV
  is scored on it.
- Coordinates: MaveDB served no GRCh38 postMapped record for any intron-side SNV of
  the four score sets (mapping state `incomplete`), and its mapped-variants endpoint
  is gone (HTTP 410; replaced by variant-details). By the author's decision
  (2026-10-20) the coordinates are made as DDX3X's were, through the atlas mapping
  module and the MANE Select transcript. Checks, all 100% on 4,588 SNVs: every
  reference base against the Ensembl release-112 chromosome
  (`scripts/stage_E16/check_reference.py`, recorded in
  `external/newgenes_reference_check.json`) and against the tracked slice; every
  position against MaveDB's own GRCh38 position of the exonic base its c. position
  names, plus the signed offset, which is independent of the mapping; and, as a check
  of the transcript model against itself, every offset against its exon structure,
  on the nearer side.
- Counts, as the literature check found them: 3-10 bp 360 / 396 / 429 / 96 SNVs (42 /
  48 / 55 / 7 damaging), 11-50 bp 1,010 / 623 / 867 / 294 (34 / 52 / 19 / 3), for CTCF /
  SFPQ / RAD51D / XRCC2. The four genes add 108 damaging variants at 11-50 bp, where
  the seven genes have 225. XRCC2 has fewer than ten damaging variants in either band.
- The entry point: the docs page is an output; `phase1/data/evidence/refs/` is a
  fetched-input folder, like `clinvar/`; `--fetch-inputs` also fetches the four
  chromosomes (165 MB; sha256 against each raw manifest), which only the
  reference-base check reads, trying each three times and going on without it if that
  fails, with no wrong or half-unpacked file left behind; `--check` checks every raw
  manifest. `evid_external` gains `intronic_snv_regex`, the DDX3X pattern made
  reusable on RefSeq transcripts (DDX3X's own pattern is unchanged).
- `.gitattributes` keeps `external/*_raw/` byte for byte (no line-end conversion):
  MaveDB serves its score tables with CRLF line ends, and the repository's rules for
  `*.csv` and `* text=auto` would otherwise normalise them to LF when they are
  committed, so that the stored file, and every checkout of it, would differ from the
  sha256 each manifest pins. The DDX3X deposits there are LF already and do not change.
- `tests/_checkout.py` leaves `refs/` out of a release archive's file list; the README
  (stage count, the chromosome download, the build-time sentence), the `reproduce.sh`
  header and the CI workflow (its header, and its path filter, which now includes the
  docs page) describe the branch.
- Run records of `scripts/stage_E16/`: each logs the sha256 of its script and of
  `_runlog.py` as they ran (the commit it logs predates them), hashes an input when it
  is read, and writes a path outside the repository from `~/`. The first records of
  `fetch_reference.py` and `check_reference.py` (2026-10-19, 16:44 and 16:48 UTC) leave
  some inputs out; both scripts were rerun offline from the cache at 18:23 UTC with
  complete records, and `probe_raw.py` once it used the stage's own layers. The
  earlier records are kept.
- E15's stage-order test now asks that E15 run after the figures and before the
  delta, rather than next to it, since E16 sits between.
- Every output of E1-E15 is unchanged.

## After 3.0.1: a fifth new external gene, SBDS, and two checks (stage E16b, 2026-10-20)

Why: the author added SBDS as a fifth further external gene (MaveDB
urn:mavedb:00001284-a-1, CC0, published 2026-08-24; HAP1 saturation genome editing,
cell survival, the platform of the E16 genes), for the same prospective test. As in
E16, nothing is fitted or scored. Two checks come with it, and neither changes a
label, a threshold or a row of the analysis set: the IGVF copy of the SBDS scores
against MaveDB's, and the BARD1 and PALB2 label files of the analysis set against
every archived version of their MaveDB score sets.

- Stage 24, `src/evid_newgenes.py --gene SBDS`, between E16 and `evid_delta` (now 25).
  The stage's gene list and file names became parameters. Run without the argument it
  writes E16's 15 files byte for byte as they were committed at E16 (a35ec6d), the
  provenance build time aside; a test holds them to that commit where the clone has
  it. With it:
  `sbds_{record_check,counts_by_offset,label_summary,coordinate_check,offset_sample,
  controls}.csv`, with E16's columns and strata, SBDS's rows only and no total row;
  `sbds_variants.parquet` and its provenance; `docs/sbds-record-check.md`. The record
  is checked against the author's of 2026-10-20 (URN, licence, publication status,
  first-published month); the current version is -a-1, the only one.
- Labels. The order of class sources was fixed by the author on 2026-10-20, before
  any predictor was run on SBDS, and the record page prints it: (a) exactly one set of
  investigator-provided ranges, checked variant by variant against a class column of
  the score table if there is one, any disagreement stopping the stage; (b) otherwise
  the score table's per-variant class column, copied as it stands; (c) otherwise stop.
  SBDS has no investigator-provided ranges, so (b) applies: the column
  `functional_consequence`, mapped by wording only (functionally_abnormal 1,
  functionally_normal 0, indeterminate NA with the row kept; any other value
  `unassigned` and listed, which stops the stage within 50 nt). No other value occurs.
  The labels come from MaveDB only; the IGVF copy never gives one.
- Counts: 507 intron-side SNVs on chromosome 7, minus strand: 1-2 nt 45 (37 abnormal),
  3-10 nt 189 (10 abnormal, 177 normal, 2 indeterminate), 11-50 nt 258 (all normal),
  beyond 50 nt 15 (the deepest c.460-55). The 3-10 nt band is evaluable exactly at the
  threshold of 10 damaging; the 11-50 nt band has no damaging variant, so it can test
  only the frozen thresholds' false positives among normal variants. The author's
  survey gave the same counts.
- Coordinates, made as E16's; chromosome 7, fetched for XRCC2, is reused. All 507
  reference bases match the Ensembl release-112 chromosome
  (`external/sbds_reference_check.json`), every position matches MaveDB's exonic
  anchor plus the offset, and every offset the exon structure. One check more: SBDS's
  score table carries the authors' own chrom/pos/ref/alt, and all 507 agree with the
  mapped coordinates. MaveDB gives no postMapped allele for any intron-side SNV.
- IGVF (`sbds_igvf_crosscheck.csv`): the score file IGVFFI9790OCYB of analysis set
  IGVFDS8094DYWC, downloaded anonymously by the new `scripts/stage_E16/fetch_igvf.py`
  into `external/igvf_sbds/` with a manifest. Aligned on the c. HGVS, both sides hold
  the same 3,315 variants; the class differs on 60, five of them intron-side SNVs
  within 50 nt, each listed with both sides' scores and classes (the scores differ
  too). At 3-10 nt MaveDB has 10 abnormal and IGVF 9.
- Label versions (`label_version_check.csv`, `label_version_diff.csv`): BARD1
  urn:mavedb:00001250-a-1 to -a-3 and PALB2 urn:mavedb:00001259-a-1 and -a-2, fetched
  by `fetch_mavedb.py --label-versions` into `external/label_versions/`, compared on
  the analysis set's intron-side SNVs at 1-50 nt. A version written in c. notation is
  aligned on the c. HGVS; the two -a-1 versions write only g. notation and are aligned
  on the g. HGVS built from the set's GRCh38 allele. BARD1's file differs from -a-1
  and -a-2 on 3 variants, all at 11-50 nt (0.16% of its 1,931 SNVs at 3-50 nt); the
  current -a-3 has neither ranges nor a class column. PALB2's file differs from the
  current -a-2 on 47 variants (4.51% of 1,041) and from -a-1 on 49 (4.71%), each
  intermediate in the file and normal in the archive. That is over 1% for every
  version, a stop condition of this stage's instructions; by the author's decision of
  2026-10-20 the outputs are committed as they are, nothing is changed, and where the
  two label files came from is the next stage's question (E16c).
- Fetching: `fetch_mavedb.py` gains `--gene`, `--label-versions` and `--a1`;
  `fetch_reference.py` and `check_reference.py` gain `--gene`. The `--a1` run
  (2026-10-20, 10:39 UTC, the same run that fetched the label versions) took the -a-1
  score tables of the four E16 genes for a coordinate check that the author then
  cancelled: those tables name intron-side variants only in genomic notation, and the
  two coordinate checks E16 already has (each reference base against GRCh38, each
  position against MaveDB's exonic anchor, 4,588 of 4,588 each) leave nothing a
  position-based comparison could add. Its files were deleted, not committed; its run
  record stays as the record of the attempt. Every new run record logs the sha256 of
  the script as it ran, and each committed script is that file.
- The entry point: `sbds_raw/` joins the raw manifests that `--check` checks and
  `--fetch-inputs` reads (SBDS names XRCC2's chromosome, which, once fetched, is only
  checked again); `--check` also checks the IGVF and label-version manifests and the
  SBDS reference-base record; the SBDS page is an output; `reproduce_inputs.txt`
  lists the 31 new inputs. `.gitattributes` keeps `external/label_versions/*_raw/` and
  `external/igvf_sbds/` byte for byte, as it keeps `external/*_raw/` (a `*` there does
  not reach a folder below): MaveDB serves the label versions' score tables with CRLF
  line ends. The reference cache registers SBDS's transcript model, Ensembl lookup and
  sequence slice; no earlier entry changes. The README and the CI path filter name the
  new page.
- Tests: the E16 checks now run on SBDS as well, with the column's classes in place of
  ranges; new ones check that no file under `phase1/data` or `phase1/reports` names
  SBDS, its transcripts or its region apart from E16b's own files and the cache
  entries its mapping read, and that none of those has a predictor column; that the
  class source follows the order fixed in advance; that a class outside the table
  stops the stage; the IGVF and label-version tables, recounted from the raw files;
  and the no-argument run against E16's commit. Two existing tests changed without
  changing what they test: E16's stage-order test asks that E16 come after E15 and
  before the delta rather than next to it, and that the raw-manifest list hold SBDS;
  the read-back test lists the new label-file read, which reads text only.
- Every output of E1-E16 is unchanged.

Omitted from the E16 entry above, all in E16's commit (a35ec6d): `reproduce_inputs.txt`
gained the 42 inputs E16 reads (its deposits, cache files and reference-base record);
the reference cache's `MANIFEST.json` registered 12 files, changing no earlier entry;
and `phase1/tests/test_evid_reproducibility.py` gained a test that serves the
chromosome downloads locally and checks that a wrong sha256, a 404 or a stream that
does not unpack leaves no wrong or half-written file. `--check` also checks the
recorded reference-base check, `external/newgenes_reference_check.json` (the sha256
of each variant table it ran on), not only the raw manifests. And XRCC2's
`variant_details.ndjson` holds, inside one GA4GH digest, four letters that read, case
aside, as a journal's abbreviation: a random fragment of the digest, not a mention.
The file is kept byte for byte as served.

## After 3.0.1: the label source audit (stage E16c, 2026-10-21)

Why: E16b found that the PALB2 label file differs from the score set Supplementary Table
S1 cites, and none of the seven label files had ever been compared, variant by
variant, with the source it cites. This stage does that, read-only: no label file,
configuration entry, threshold or analysis-set row changes here.

- Stage 25, `src/evid_label_audit.py`, between E16b and `evid_delta` (now 26). For each
  gene and each candidate rule a row in `label_audit_summary.csv`:
  - R1, the cited version's investigator-provided calibration: its score ranges, or,
    for a class-based calibration, MaveDB's per-variant membership of its classes;
  - R2, a class column of the score table;
  - R3, the ranges with the table's 95% interval;
  - R4, the rule the cited paper states, quoted verbatim, every quote checked against
    the paper's full text;
  - R5, R1 and R2 on the other versions of the chain;
  - R6, a paper supplementary table the author made a source (BRCA2 only).

  Every differing variant goes in `label_audit_diff.csv`. For a gene no rule
  reproduces, `label_audit_impact.csv` gives the counts the closest rule would give
  the analysis set. The cited table's variants the label file does not carry go, by
  stratum, in `label_audit_unclassified.csv`. The record page is
  `docs/label-source-audit.md`.
- Findings:
  - Reproduced:
    - BRCA1, by R1 and by R4 (the paper's Figure 3a cut-offs);
    - RAD51C, by R1, R2 and R4;
    - VHL, by R1, on its 2,200 classed variants. Another 68, all exonic or UTR and
      none in the analysis set, have no class in the archive.
    - BAP1, by R1 as abnormal / normal / intermediate, and by R4, the paper's median
      split, which also gives the sub-class names;
    - BRCA2, by R6 (Supplementary Table 3, 6,959 of 6,959).
  - Not reproduced: PALB2 and BARD1. The closest rule, R1 on the cited a-2, differs
    on 4.51% and 0.16% of the 3-50 nt window. Every one of PALB2's 47 there is
    intermediate in the file and normal in the archive. R3 differs far more (52.9% and
    24.0%).
  - No calibration of a cited version was modified after the label files were
    committed (2026-07-08).
- The author's decisions of 2026-10-21, printed in the record page:
  - relabel PALB2 and BARD1 by R1 and refit (the next commits);
  - BRCA2's label source is the paper's Supplementary Table 3;
  - VHL's 68 unclassified variants stay unclassified;
  - BAP1's sub-class names are kept. They come from the median the paper prints,
    -0.1260642; c.1132del, which scores exactly that median, is strongly depleted in
    the file.
- Sources, in `scripts/stage_E16c/`, each script with a run record:
  - `fetch_sources.py` fetches each cited score set with its chain and its experiment,
    and, for a class-based calibration, MaveDB's per-variant membership
    (`/score-calibrations/{urn}/variants`). BARD1's and PALB2's versions are E16b's,
    reused.
  - `check_paper_quotes.py` checks every quote in each `paper_rule.json` against the
    paper's full text from PMC or Europe PMC, and the IDs against NCBI's converter.
    The audit refuses a rule file the check has not passed.
  - Not committed: three memberships (BAP1's and RAD51C's two, 9-18 MB) and the BRCA2
    xlsx. They sit, as served, in the gitignored
    `phase1/data/evidence/label_audit_cache/`. A compact extract of each (variant ->
    class) is committed, and the audit reads the extracts, so it runs offline.
  - VHL's membership (1.4 MB) is committed as served.
- The entry point:
  - `--fetch-inputs` fetches the cached files, checks each against the sha256 in its
    manifest, and checks every extract row for row against its source. A failure
    stops the run.
  - `--check` checks the audit's manifests and quote checks.
  - The audit's raw folders are kept byte for byte (`.gitattributes`), and their 45
    files are inputs in `reproduce_inputs.txt`.
  - The cache is a fetched-input folder (`NOT_OUTPUTS`, `.gitignore`,
    `tests/_checkout.py`).
  - The record page is an output, and the CI path filter includes it.
- Tests (`phase1/tests/test_evid_label_audit.py`):
  - the stage rerun into a scratch folder leaves the label files, `config.py` and the
    analysis set unchanged and rewrites its own files byte for byte;
  - every gene has a row for every rule, and the counts add up;
  - the diff lists each difference once;
  - the impact table covers exactly the genes no rule reproduces, starting from the
    counts files;
  - the quotes were checked, and every number a rule uses is in a quote;
  - each extract matches its source; a source not in the cache is skipped, with the
    reason printed.

  One existing test changed without changing what it tests: E16b's stage-order test
  asks that E16b come before the delta, rather than right before it.
- Every output of E1-E16b is unchanged. The README still counts 25 stages: README
  changes wait for the manuscript's next revision sheet, by the author's decision of
  2026-10-21.

## After 3.0.1: the seven label files built from their sources (stage E16d, 2026-10-21)

Why: no script in the repository produced the label files in
`phase1/data/assay_labels/`; they entered on 2026-07-08 (commit 10d94bb) with no record
of how. This commit adds the script and runs it into a staging folder. No label file in
use changes here; the next commit replaces two of them.

- `scripts/stage_E16d/build_labels.py` builds each file from the source the label
  source audit pinned, by the rule the author fixed on 2026-10-21:
  - BRCA1, BARD1 and PALB2: the cited version's investigator-provided ranges (R1);
  - BRCA2: the paper's Supplementary Table 3 (R6);
  - RAD51C: the score table's class column (R2);
  - VHL: MaveDB's membership of the class-based calibration, the 2,200 classed
    variants (R1);
  - BAP1: the paper's rule, with the depleted variants split at the printed median.
    The one variant on it is strongly depleted, as the author decided (R4).

  The classes take the wording of the files in use, through an explicit table. The
  columns, delimiter, form of the HGVS and row order follow the files in use too.
  BRCA1's file is ordered by c. position, intron offset and alternative base; the
  others keep the source's order.
- It writes `phase1/data/assay_labels/_generated/` and
  `phase1/data/assay_labels/PROVENANCE.json`. For each gene the record gives:
  - the source and its files with their sha256;
  - the calibration, with its id and dates;
  - the rule and the wording;
  - the script's and the output's sha256;
  - `in_use`: whether the label file in use is, byte for byte, the output.
- Result:
  - BRCA1, BRCA2, RAD51C, VHL and BAP1 are rebuilt byte for byte.
  - PALB2 and BARD1 are the cited a-2's ranges, variant by variant. They differ from
    the files in use on exactly the 359 and 13 variants the audit listed.
  - No score of a ranged source is unassigned or sits on a point two ranges share.
- `scripts/stage_E16d/compare_runs.py`, added here and not run, compares the evidence
  outputs of two commits:
  - for every predictor x stratum x tier x side, the threshold, the interval LR there
    and the held-out folds reaching it, before and after;
  - the count changes;
  - the audit's predicted counts, cell by cell; any difference stops it.

  The next commit runs it.
- `scripts/stage_E16d/check_switch_commit.py` checks that a commit changes no code,
  configuration or workflow file.
- `scripts/reproduce_evidence.py`: `--verify` does not expect a run to rewrite
  `phase1/reports/evidence/label_change_effect.csv`, which that comparison writes once
  (`COMPARISON_REPORTS`, guarded by a test).
- Tests (`phase1/tests/test_label_build.py`):
  - the record matches the generated files and their sources;
  - the script rebuilds its committed outputs. A gene whose source is not in the
    cache is skipped, with the reason.
  - the five files are rebuilt byte for byte;
  - PALB2 and BARD1 are the cited ranges, and differ from the July files on exactly
    the audit's variants, compared with the audit's own list at commit 76ed7c2;
  - after the switch: the switch commit holds no code, and the July files are kept
    in `superseded/`.
- Two existing tests changed so that the next commit, which may change only data,
  passes them. What they test does not change.
  - `test_evid_read_back.py` held Table S10's band LR at the fitted Atlas threshold
    (3-10 nt) to the printed 28.2. It now holds it to that row of
    `arms_at_tool_threshold.csv`, formatted as the table formats it. The band LR
    moves with the number of normal variants.
  - `test_evid_label_audit.py` asked every stratum of an empty impact table. It now
    asks every stratum of each gene listed.

  A dry run of the next commit in a separate worktree, with the two files replaced
  and nothing else, passed the whole suite with these changes.
- Every output is unchanged.

## After 3.0.1: PALB2 and BARD1 relabelled from their cited sources, and refitted (stage E16d, 2026-10-21)

Why: the label source audit (stage E16c) could not reproduce the PALB2 and BARD1 label
files from the score sets Supplementary Table S1 cites. Under the cited a-2 versions'
investigator-provided ranges they differ on 359 and 13 variants; 47 and 3 of those
fall in the analysis window. On 2026-10-21 the author decided to relabel both by those
ranges and refit. The decision was made before any refit was seen, so whatever the
refit shows stands.

- This commit holds no code, configuration or workflow file
  (`scripts/stage_E16d/check_switch_commit.py`, and a test). Every output it changes
  changes because of the two label files.
- `phase1/data/assay_labels/palb2_mavedb.csv` and `bard1_woo2025.tsv` are now the files
  `scripts/stage_E16d/build_labels.py` built from the cited ranges (the previous
  commit's `_generated/`, unchanged).
- The files of 2026-07-08 moved, unchanged, to `phase1/data/assay_labels/superseded/`,
  with a README, for the sensitivity comparison; no stage reads them.
- `PROVENANCE.json` marks both new files as in use.
- The whole pipeline was rerun. The counts moved as the audit predicted, in every cell
  of `label_audit_impact.csv` (44 of 44, `docs/label-change-effect.md`):
  - PALB2: no damaging count changes; normal 324 -> 358 at 3-10 nt and 449 -> 462 at
    11-50 nt.
  - BARD1: only 11-50 nt changes; damaging 21 -> 19, normal 1,414 -> 1,417.
  - The seven genes at 3-50 nt: labelled 7,612 -> 7,660, excluded 385 -> 337, damaging
    795 -> 793. At 3-10 nt damaging stays 570; at 11-50 nt it is 225 -> 223.
- `scripts/stage_E16d/compare_runs.py`, run against the previous commit, wrote
  `phase1/reports/evidence/label_change_effect.csv` and `docs/label-change-effect.md`.
  In sample, 3 of 560 predictor x stratum x tier x side cells changed tier:
  - CADD lost PP3 Supporting at 11-50 nt;
  - the Atlas splice-site score gained BP4 Supporting at 3-50 nt;
  - the Atlas splice-site usage score lost BP4 Moderate at 3-10 nt.

  Held out, one cell changed the number of folds reaching its tier: the Nucleotide
  Transformer, PP3 Supporting at 3-50 nt, 0 -> 1.
- The label source audit now reproduces all seven files.
- The manuscript draft prints four counts that no longer match their named file:
  7,612, 8,453, 224 and 225. Neither the manuscript nor the README is changed here.

## After 3.0.1: the reference image, defined (stage B6, 2026-10-21)

Why: the committed outputs were made on the author's macOS arm64. CI on another Mac
reproduces them numerically; CI on Linux does not (the figures ask for Arial, which
Linux lacks, and numpy and scipy run another BLAS). On 2026-10-20 the author chose one
Docker image as the reference environment, every output to be regenerated in it as the
new committed baseline. This commit defines the image; no output changes.

- `container/Dockerfile`: `python:3.12.13-slim-bookworm` pinned by its linux/amd64
  manifest digest; `git` (which `--verify` calls) pinned to the version apt offers in
  that image; the packages of `requirements-evidence.lock.txt`, from wheels only, each
  checked against its hash in `container/requirements.linux-amd64.txt`, then `pip
  check` and the lock compared with `pip freeze`. No data and no repository code go
  into the image: the repository is mounted when it runs.
- The pins come from command output, not by hand: `scripts/stage_B6/base_image.py`
  read the digests from the registry, and `scripts/stage_B6/pin_versions.py`, in the
  pinned base image, asked apt for git's version and downloaded and hashed every
  wheel. Both write `container/ENVIRONMENT.json`; a test holds the Dockerfile to it.
- `container/run.sh` runs in the image as the host user, with `reproduce.sh`'s thread
  settings, `PYTHONHASHSEED=0`, `TZ=UTC` and `LC_ALL=C.UTF-8`, in two steps:
  `fetch`, online, and `verify`, with no network. The OpenBLAS kernel type and numpy's
  dispatch are left empty, to be set from probes of two CPUs before the rebaseline.
- `scripts/reproduce_evidence.py --fetch-only` fetches and checks what `--fetch-inputs`
  does, the atlas archive included, and runs no stage: the online first step in the
  image. Nothing else in the entry point changes, and the native run here stays
  REPRODUCED.
- `.github/workflows/container.yml`, started by hand only: build and push the image
  (privately), or, given it by digest, probe the runner's CPU or run the whole
  reproduction in it. Its push trigger only registers it, because GitHub dispatches a
  workflow that is not on the default branch only after it has run once; every job is
  skipped on a push.
- A trial in an image built here from this Dockerfile (not the reference image, which
  CI builds): every stage ran with no network, in 34 minutes under amd64 emulation.
  Against 0161aad the comparator said NOT REPRODUCED, listing exactly the files CI's
  native Linux run lists: 21 figures, 21 numerically equal files, and 2 that differ.
  The TP53 manifest's reproduction check counts 23 of 24 re-scored Nucleotide
  Transformer values identical rather than 24, one being 4.4e-16 apart, which the
  stage attributes to pandas' default CSV float parser on x86_64. The
  supplementary-table manifest records that manifest's hash.

## After 3.0.1: the figures' typeface from the repository, and the image's CPU settings (stage B6, 2026-10-21)

Why: the figures asked for Arial, which only some machines have, and matplotlib falls
back to another font without an error where it is missing, so a figure's typeface
depended on the machine. They now load Liberation Sans, metrically compatible with
Arial, from files shipped with the repository. And the reference image, built from the
previous commit, is recorded, and its maths libraries are held to the same instructions
on every CPU, before the outputs are regenerated in it.

- `phase1/assets/fonts/`: Liberation Sans regular, bold, italic and bold italic, and the
  SIL Open Font License 1.1 they are released under, unchanged from the project's 2.1.5
  release (`scripts/stage_B6/fetch_fonts.py`; `SOURCE.json` records the archive's URL
  and sha256 and each file's). The figures open the regular, bold and italic faces.
- The figure modules change in their font settings only:
  - `evid_figures.py` loads these files by path (matplotlib's `addfont`) and drops any
    copy of the family installed on the machine. It names Liberation Sans alone,
    mathtext's faces included, with no fallback.
  - `evid_fig_data.py`'s draft PNGs use the same.
  - `evid_practical_metrics.py`'s note on the slide's label check names the new face.
- `phase1/tests/test_evid_figure_fonts.py` draws every figure again. It holds every font
  file matplotlib opens to that folder, and checks that no glyph is missing.
- Run here (macOS) against the previous commit, the comparator says NOT REPRODUCED with
  "differs (figure): 21; identical: 444; identical (build time aside): 6". Only the 21
  figures changed: every table, `fig_data/` included, is byte for byte the same, and no
  manifest records a figure's hash. The figures in this commit are that run's.
- The reference image was built by `.github/workflows/container.yml` at 9b15c34 (run
  39621528080): `ghcr.io/oneone00-11/variant-fm-benchmark-reference@sha256:0dc03436...`,
  image ID (config digest) `sha256:8ecf15fc...`.
  - Anyone can pull it: the registry answers an anonymous request for its manifest. The
    stage meant to keep it private; on 2026-10-21 the author accepted it as public.
  - `container/ENVIRONMENT.json` records the image as built and as loaded here
    (`scripts/stage_B6/record_image.py`), and what is inside it
    (`record_environment.py`: Python 3.12.13, the 37 packages of both locks, git
    1:2.39.5-0+deb12u3) and the font files.
- CPU settings in `container/run.sh` (B6.4):
  - Probes in the image (`container/cpu_probe/`): both OpenBLAS copies chose Zen on
    GitHub's AMD EPYC 7763 and Nehalem on the Mac's emulated CPU; numpy dispatched up
    to AVX2 on both.
  - Now `OPENBLAS_CORETYPE=Nehalem`, and numpy's dispatch is held to SSE4.2
    (`NPY_DISABLE_CPU_FEATURES`: AVX and every target above it).
  - A probe on the Mac under these settings shows both in force. The runner is probed
    again once this commit is pushed.

## After 3.0.1: every output rebaselined in the reference image (stage B6, 2026-10-22)

Why: the reference image is the reference environment (the author's choice of
2026-10-20), so the committed outputs are now the ones it writes. This commit holds no
code, configuration or workflow file (`scripts/stage_E16d/check_switch_commit.py`).
Every output it changes changes because the stages ran in the image rather than on the
author's Mac.

- The run: the image as built at 9b15c34 (`container/ENVIRONMENT.json`, image), loaded
  from the build's `docker save` archive and checked against the build's record. It ran
  the previous commit's code, with `container/run.sh`'s settings, on the Mac under amd64
  emulation (`container/cpu_probe/mac-rosetta-pinned.json`). The inputs were fetched
  online; every stage and `--verify` then ran with no network, in 40 minutes.
- Against the previous commit, the comparator said NOT REPRODUCED, with "differs: 2;
  differs (figure): 13; differs (numerically equal): 21; identical: 429; identical
  (build time aside): 6" (`docs/rebaseline-report.md`):
  - Every printed table is byte for byte the same.
  - 21 files differ only in the last digits of their numbers; the largest relative
    difference is 1.6e-13.
  - The 13 figures that differ are PNGs whose pixels are identical: only the
    compressed image data differs, as Pillow's zlib writes it on Linux. The 8 PDF
    figures are byte for byte the same.
  - Two files the stage's rules did not allow:
    - The TP53 manifest's reproduction check counts 23 of 24 re-scored Nucleotide
      Transformer values identical rather than 24 (one, 17-7673696-C-A, is 4.4e-16
      apart), which the stage attributes to pandas' default CSV float parser on x86_64.
    - The supplementary-table manifest records that manifest's hash.

    On 2026-10-22 the author accepted both as the image writes them.
- These outputs replace the committed ones, all but the six manifests that differed in
  their build time alone, which keep theirs.
- The runner's CPU, probed again under `run.sh`'s settings, now runs the same OpenBLAS
  kernels (Nehalem) and the same numpy dispatch (up to SSE4.2) as the Mac's emulated
  CPU (`container/cpu_probe/github-linux-pinned.json`).

## After 3.0.1: the required reproduction runs in the reference image (stage B6, 2026-10-22)

- `.github/workflows/reproduce.yml` changes three jobs:
  - The required job now runs in the reference image, pulled by the digest
    `container/ENVIRONMENT.json` records. The inputs are fetched online, then every
    stage and `--verify` run with no network (`container/run.sh`). It must end
    REPRODUCED, byte for byte; REPRODUCED NUMERICALLY does not pass.
  - The native Linux job is gone; the image job takes its place.
  - The native macOS job stays, for the record only (`continue-on-error`, the author's
    decision of 2026-10-22). Against the image's outputs it cannot end REPRODUCED. On
    the rebaseline commit it ended NOT REPRODUCED, for three reasons: Pillow compresses
    the PNG figures into other bytes, with the same pixels; the TP53 reproduction check
    counts 24 of 24 on arm64; and numbers differ in their last digits.
- README: the reproduction section only. It names the image and its digest, gives the
  two steps, says the native path agrees in the numbers only, and lists what is not
  tested.
- Before this, the full run in the image, against the rebaselined outputs, ended
  REPRODUCED on two CPUs: GitHub's runner (AMD EPYC 9V45; `container.yml` run
  39638815249) and the Mac under amd64 emulation.
- The previous commit holds `test_evid_tp53_extend.py` to the TP53 stage's own
  tolerance. That test had failed on the rebaseline commit. The change lies outside
  the stage's scope; the author approved it.
