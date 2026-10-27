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

## After 3.0.1: the PNG figures compared pixel by pixel, and reproduce.sh names the reference image (stage B6, 2026-10-22)

Two decisions of the author that the CI switch left out. No output changes, and the
comparator's categories and criteria are as they were.

- `scripts/stage_B6/compare_png_pixels.py` decodes every PNG under the comparator's
  output folders and compares it, pixel by pixel, with its committed version. Each one
  is one of three: bytes identical; bytes differing with every pixel identical; or
  pixels differing, with how many and where.
  - The macOS job runs it after the comparator and writes its table to the job summary.
  - Run once on the rebaseline (e12220c against 7125246), it finds all 13 PNGs the
    comparator called "differs (figure)" identical in every pixel, differing only in
    their compressed bytes (`docs/rebaseline-png-pixels.md`).
    `docs/rebaseline-report.md` is left as it was written.
- `reproduce.sh` now prints a note at its end, after the run, and its exit status is
  still the run's. The note gives three things:
  - the reference image by digest, with its two commands;
  - the differences a native run should be expected to show: PNG bytes with the same
    pixels, the TP53 manifest's exact count, and last digits;
  - the pixel comparison's command.

## After 3.0.1: AlphaGenome scored by its official definition, not yet read (stage E17, 2026-10-22)

Why: the study's `alphagenome` column is the companion atlas's quantile composite. On
2026-10-19 the author decided to score AlphaGenome by its own definition instead. The
previous commit pinned that definition from AlphaGenome's notebook and scored a pilot;
this one scores the whole analysis set. The pipeline reads none of it yet, and every
output is unchanged.

- The scores: every row of the analysis set, 8,853 variants of the seven training
  genes, with client 0.10.0 (`phase1/data/evidence/alphagenome_official/all_scores.tsv`
  and its provenance). One call per variant, two at a time, from 10:01 to 12:30 UTC;
  the Mac slept for 47 minutes in between and the run went on when it woke. No retry
  and no failure. No variant lacks a row for any of the three scorers, no raw score is
  NaN, and the smallest raw score is 0.0.
- Two checks that the service's model did not change during the run:
  - the pilot's 120 variants, scored again from the full run's own cache, equal the
    pilot's scores bit for bit (120 of 120);
  - the ten of the determinism check, scored afresh after the full run had ended
    (`drift_scores.tsv`), equal them as well (10 of 10).
- `NOTEBOOK_CONTEXT.json` quotes, verbatim, three markdown cells of the same notebook:
  what the score is for (cell 1), why the window is 1 Mb (cell 10), and the score's
  range, with a recommended cutoff still being worked on (cell 16).
  `NOTEBOOK_SOURCE.json` is unchanged.
- `ENVIRONMENT.json` and the new provenance records quote the client's line naming the
  model a request with no version gets (`alphagenome/models/v1/dna_model.py`, line 52,
  `ALL_FOLDS`), with the file's and the wheel's sha256. The request names no version,
  as in the notebook.
- `scripts/stage_E17/score_official.py` gains `--drift`, which runs only after `--all`
  and from an empty cache of its own. The run now stops on an authentication error,
  when three variants in a row have used all their attempts, and on Ctrl+C; then only
  the calls in flight end. Tests check this with a fake client.
- The pilot's files (the examples, the pilot, the repeat, and their records) are
  unchanged.
- Not in this commit: LICENSE-DATA section 2, which is to name the new scores.
  LICENSE-DATA is an input of Supplementary Table S3, whose manifest records its
  sha256, so changing it here would change an output. It goes into the next commit,
  which changes outputs anyway (the author's decision of 2026-10-22).
- A correction. Commit 866ca33 says that the Table S10 cell it holds as a literal "does
  not depend on AlphaGenome". That cell is the AlphaGenome Atlas AVI score (`avi`). It
  does not depend on the `alphagenome` column this stage replaces, and this stage does
  not change the four Atlas columns.

## After 3.0.1: the pipeline reads AlphaGenome by its own definition, and refits its thresholds (stage E17, 2026-10-22)

Why: the switch the previous commit prepared. The definition comes from AlphaGenome's
notebook for scoring variants on their splicing effect
(`colabs/v1/splicing_variant_scoring.ipynb` in google-deepmind/alphagenome at
55943112, sha256 c39e6a6a), quoted item by item in
`phase1/data/evidence/alphagenome_official/NOTEBOOK_SOURCE.json`: the largest raw
score of each of three splicing scorers over genes and tracks, in a 1 Mb window,
combined as sites + usage + junctions / 5; client 0.10.0.

- `src/evid_build_set.py`: the analysis set's `alphagenome` is now that score, read
  back exactly from `alphagenome_official/all_scores.tsv` (equal to it bit for bit). The
  companion atlas's merged-quantile composite (client 0.7.0, a 16-kb window), which the
  column held until now, stays beside it as `alphagenome_atlas_quantile`, outside the
  panel. `alphagenome_v061` is unchanged, and no AlphaGenome column enters the fusion.
  The new column is larger-is-damaging (Spearman rho 0.388 with functional
  pathogenicity), so nothing is flipped.
- A new stage, 27, `src/evid_alphagenome_definitions.py`, writes
  `reports/evidence/alphagenome_definition_comparison.csv`, the one place the old
  composite appears: for each stratum, tier and side, both definitions' in-sample
  threshold, band ratio at it, held-out k/n (PP3, as Table 3 counts it; BP4 has none,
  since E3 writes per-fold thresholds for the PP3 side only, and gives E3's median
  held-out ratio instead) and pooled AUROC, from the functions that report them
  elsewhere. Its new-definition numbers are checked in the stage against E3, E2.2, E4
  and E15's files; its old-definition numbers equal the previous commit's, float for
  float and as Table 3, S6 and S7 print them.
- What moved for AlphaGenome, as the study prints it: Strong is still reached
  in-sample in all three strata. Held out, Strong is cleared in 6/7 (1) genes at
  3–10 bp (was 5/7 (1)), 0/6 (5) at 11–50 bp (was 0/6 (6)) and 5/7 (1) at 3–50 bp
  (was 6/7 (1)); Moderate in 4/6 (2) at 11–50 bp (was 3/6 (3)). Pooled AUROC 0.922,
  0.851 and 0.923 (was 0.915, 0.837 and 0.913). The thresholds sit on the new score's
  own scale: Table S6 prints the Strong threshold at 3–10 bp as 2.673, where it printed
  2.19988276958.
- Supplementary Table S3's AlphaGenome row takes what the column is, its scoring
  configuration and its date from E17's records, not from the atlas's resource row.
  Table S6's note no longer says AlphaGenome's thresholds lie below a ceiling of 2.2.
  The TP53 mismatch (`evid_external.DECLARED_MISMATCH`, and the TP53 manifest's
  definition) names the new analysis-set column; TP53 stays not comparable.
- LICENSE-DATA section 2 is rewritten to name every AlphaGenome value: the four Atlas
  columns, the new `alphagenome` with `alphagenome_official/` and its four score
  columns, `alphagenome_atlas_quantile`, `alphagenome_v061`; section 3's list follows.
  `docs/column-provenance.md` gains a section on the new definition, and its Atlas
  section's comparison now names the old composite. Four comments or docstrings that
  described the old composite as the current column say what it was.
- Outputs were written in the reference image. Against the previous commit,
  `scripts/stage_E17/check_alphagenome_only.py` finds 61 changed and 425 identical,
  every change explained by AlphaGenome: its rows and columns, its curves, the seven
  figures that draw it (PDF and PNG; the slide figure as PNG only) and two draft PNGs,
  the analysis set and the manifests that record it, two pages, and two summaries
  whose counts or ranges take
  AlphaGenome in (Table S11's per-gene summary, three statements of the practical
  metrics page). The fusion's outputs, `alphagenome_v061`'s, and figures 1 and 4 are
  unchanged.
- The manuscript, checked with `scripts/check_evidence_manuscript_numbers.py`: no
  number moves from matched to unmatched (62 matched; the same four unmatched since
  E16d). Three sentences and Supplementary Note 2 describe the old composite or its
  results; the stage report lists them. The manuscript is not changed.

## After the thresholds were frozen: the freeze guard, and the external genes' AlphaGenome scores (stage E18, 2026-10-23)

Why: the seven training genes' thresholds were frozen at the tag `frozen-thresholds-v1`
(tag object c75eb39, commit 1d03019, 2026-10-22T16:34:25Z). The external genes are now
scored, after the tag, for a test against those thresholds. This commit only adds data;
the pipeline reads none of it, and every output is unchanged.

- The freeze guard. `phase1/data/evidence/FROZEN_AT_TAG.json`
  (`scripts/stage_E18/freeze_guard.py`) records the sha256, as tagged, of 96 files:
  - the analysis set and its manifest;
  - every record of `alphagenome_official/`, the seven genes' scores among them;
  - the three threshold tables;
  - the `table3*` tables, and Supplementary Tables S6 and S7;
  - every `interval_lr_*.csv`;
  - the fusion's outputs.

  `test_evid_frozen_guard.py` requires HEAD to hold each of them byte for byte as
  recorded.
- AlphaGenome by the official definition (E17's code and environment, client 0.10.0,
  no model version named) for the external genes' variants:
  - DDX3X, all 1,857 rows E7 reads;
  - TP53, all 288;
  - the five further genes at 3–50 nt, 4,522: CTCF 1,370, SFPQ 1,019, RAD51D 1,296,
    XRCC2 390 and SBDS 447.

  The scores are in `alphagenome_official/external_scores.tsv` with its provenance.
  The run took 2026-10-22 16:58 to 18:13 UTC, with no retry and no failure. No variant
  lacks a row, no raw score is NaN, and none is below zero.
- Two drift checks frame it, each the pilot's ten scored afresh from an empty cache, and
  each equal to the pilot's scores bit for bit (10 of 10):
  - before the external run (`drift_2b_start_scores.tsv`);
  - after it (`drift_2b_end_scores.tsv`, at 01:11 UTC on 2026-10-23).
- Every external score post-dates the tag. The test that allowed only the seven genes
  now allows another gene only with scores after the tag.
- `scripts/stage_E17/score_official.py` gains the three modes. They write
  `ENVIRONMENT_2b.json`, so that E17's `ENVIRONMENT.json` stays as tagged.

## After the thresholds were frozen: the external genes' other columns, by the analysis set's own scorers (stage E18, 2026-10-23)

Why: a threshold fitted on a column applies to an external gene only if the gene's
column is the same variable: the same scorer, version, precision and direction
(E18.0 rule 4). The evidence is either anchors, analysis-set variants re-scored in the
same run and equal to the analysis set's values bit for bit, or a committed
value-for-value record. This commit only adds data; the pipeline reads none of it, and
every output is unchanged.

- The anchors. `phase1/data/evidence/external/e18_anchors.parquet`
  (`scripts/stage_E18/anchors.py`): three analysis-set variants per training gene and
  stratum (s3_10, s11_50), drawn with `config.RANDOM_SEED`, 39 in all (BRCA2 has no
  s11_50 variant). Every run below scored them together with its external variants.
  In every run, every anchor equals the analysis set's stored value bit for bit, in
  every column, the Walker basis's four components included (39/39).
- The runs, all after the tag, each in the column's own pinned environment
  (`scripts/stage_E18/run_*.py`, the scorers unchanged):
  - SpliceAI at distance 50 (the atlas's `score_fullprec.py`) and Walker's basis at
    distance 4,999 (`evid_score_spliceai_walker.py`): 6,706 variants each, the 39
    anchors and all 6,667 external variants;
  - Pangolin (the atlas's `score_fullprec.py`, on one thread, in six processes): the same
    6,706;
  - CADD, phyloP, phastCons and GPN-MSA (the atlas's scorers, run from byte-identical
    copies of the release's files): the anchors, DDX3X's 1,857 and the five further
    genes' 4,522. The sources as served: GPN-MSA at revision cf1718a9 and the pinned
    size, the CADD v1.7 file of 29 January 2024, the UCSC tracks of 8 May 2015;
  - the four Atlas columns (`evid_score_avi.run`, client 0.9.0, the scorer unchanged
    and its client wrapped so that an authentication error or a batch failing every
    attempt stops the run): the anchors, TP53's 288 and the five further genes' 4,522.

  Each scorer first ran the anchors and 50 variants to measure throughput (E18.2);
  every projection of a full run was under three hours. Every variant was scored; no
  value is NaN.
- Not scored again, because a committed record already shows them to be the analysis
  set's variables (`docs/external-column-identity.md`):
  - TP53's CADD, phyloP, phastCons and GPN-MSA: the published study's scorers, value for
    value the atlas's (`docs/column-provenance.md`, section 2); GPN-MSA's raw sign is
    turned to the analysis set's direction when the external test reads it (E18.0
    rule 5);
  - DDX3X's four Atlas columns, from the analysis set's own run (`avi.parquet`).
- Nucleotide Transformer is not in this commit. By the author's decision of 2026-10-23
  it was scored on an RTX A6000 in the original run's recorded environment, and has a
  commit of its own.
- The files: per gene and run `<gene>_<run>.parquet` (DDX3X's full-precision Pangolin
  `ddx3x_pangolin_fullprec.parquet`; a column scored again `<gene>_<run>_rescored.parquet`),
  each with its provenance and environment; the anchors of each run,
  `e18_anchors_<run>.parquet`; the five further genes' Atlas columns,
  `newgenes_avi.parquet`; and the five further genes with every column,
  `newgenes_scored.parquet` with its manifest. `scripts/stage_E18/build_columns.py`
  writes them in the reference image, and its `--offline` rebuild, from these files
  alone, is byte-identical. DDX3X's file names avoid `ddx3x_spliceai_d50.parquet` and
  `ddx3x_nt.parquet`, which E7's merge would read if present.
- Three of Pangolin's six processes finished within the same second, and the run record
  is named to the second, so worker 3's record overwrote those of workers 2 and 5. Their
  caches, environment records and printed anchor checks survive, and `build_columns.py`
  checks every anchor again from the caches.
- LICENSE-DATA covers the new files from the commit that switches the external test
  to them (E18.0 rule 10).
- Tests: `test_evid_e18_columns.py`. The SBDS test also reads the manifest beside an
  assembled table for its scoring time.

## After the thresholds were frozen: NT on the original GPU does not reproduce its anchors, so NT stays out (stage E18, 2026-10-23)

Why: the analysis set's `nt` column is the companion atlas's, scored on a rented RTX
A6000 whose environment the atlas recorded. On this Mac the same scorer reproduces 4 of
the 39 anchors bit for bit on the CPU and 2 on its GPU. By the author's decision of
2026-10-23, NT was scored again on an RTX A6000 in the recorded environment, and the
rule was fixed before the run: NT enters the external test only if all 39 anchors are
identical. This commit only adds data; the pipeline reads none of it, and every output
is unchanged.

- The package (`scripts/stage_E18/nt_gpu_package.py`, `nt_gpu_run.py`, `nt_gpu_README.md`):
  the anchors and the 6,667 external variants, reference slices checked window by
  window against the whole chromosome, the atlas's `models/nt/score.py` verbatim, and a
  runner that starts only on an RTX A6000 (compute capability 8.6) with torch
  2.11.0+cu128, transformers 4.46.3, `pip freeze --all` equal to the atlas's record, and
  the model at revision 06615c16. What the original run left unrecorded comes from
  `score.py`'s defaults: float32, one variant per forward pass, the 6,000-nt window.
- The run (RunPod, 2026-10-23 10:53 to 11:05 UTC): every check passed. The driver,
  which the original run did not record, was 595.91.07 (CUDA 13.2). A first attempt
  stopped before scoring anything, because the pod's template sets
  `HF_HUB_ENABLE_HF_TRANSFER=1` and the recorded environment has no `hf_transfer`; the
  variable was unset rather than the package added. Its log is kept.
- The result: 7 of the 39 anchors identical bit for bit; the largest difference 4.77e-6
  (relative 1.3e-4), 2 to 2,048 float32 units in the last place. A trial run of the
  anchors and 50 variants before it gave the same values to the last bit, so the
  difference is systematic. The analysis set's NT values are float32 values stored as
  such (8,851 of 8,853), so it is not an effect of how they were stored.
- Hence, by the rule fixed beforehand: NT enters the external test for no external
  gene, TP53's published-study NT column included; the external test runs on the other
  twelve columns; no tolerance, and the seven genes are not scored again.
- `phase1/data/evidence/external/nt_gpu_run/`: the pod's environment record, its three
  logs and both manifests, as served; `anchor_check.tsv`, the 39 anchors' stored and
  re-scored values and their difference; `record.json`. The 6,667 external values are
  not the analysis set's variable and are not committed; `record.json` names the two
  result tables by sha256, and they stay in the gitignored cache.
- `docs/external-column-identity.md`: the NT cell of each gene reads "no column of the
  same variable" and cites the anchor check; `newgenes_scored.manifest.json` says the
  same.

## After the thresholds were frozen: DDX3X and TP53 read the analysis set's variables (stage E18, 2026-10-23)

Why: E18.2 made each external gene's columns the analysis set's own variables; E7, the
external test, now reads them. One reason, one commit: the outputs change only for
DDX3X's and TP53's switched columns and what is derived from them
(`scripts/stage_E18/scope_c3.json`, checked by `scripts/stage_E18/check_scope.py`
against the previous commit).

- DDX3X (`evid_external.DDX3X_SCORE_FILES`): SpliceAI at distance 50 and AlphaGenome by
  the official definition are added; Pangolin is full precision, not the CLI's two
  decimals; the Walker basis, CADD, phyloP, phastCons and GPN-MSA are read from the
  files scored again with anchors, value for value the files read before.
- TP53 (`evid_external.TP53_SCORE_FILES`): AlphaGenome is the official definition, scored
  after the tag; the client-0.6.1 column is kept beside it as `alphagenome_v061`, which
  has no frozen threshold and is not tested. SpliceAI, Pangolin, the Walker basis and the
  Atlas columns are read from the files scored again with anchors (SpliceAI differs from
  the table's by at most 1.1e-16, through the TSV the table's values once passed; Pangolin
  by at most 1.3e-9; the Walker basis and the Atlas columns not at all). GPN-MSA, on the
  published study's raw sign, is turned to the analysis set's direction (E18.0 rule 5).
  CADD, phyloP and phastCons are the table's, value for value the analysis set's scorers'
  (`docs/column-provenance.md`, section 2).
- NT has no column of the same variable on either gene (the rule fixed before the GPU
  run; `nt_gpu_run/`): TP53's published-study NT is no longer tested, and DDX3X gets
  none. The column check names NT with its reason.
- `DECLARED_MISMATCH` loses TP53's entry and is empty; the TP53 table's description of
  its 0.6.1 column no longer points to it.
- The column check now reports, for each gene, twelve columns comparable and NT named.
- LICENSE-DATA (E18.0 rule 10): section 1 names stage E18's score files of the predictors
  it covers, section 2 its AlphaGenome and Atlas files; Supplementary Table S3's record
  in the supplement manifest carries LICENSE-DATA's new hash, and that of
  `docs/column-provenance.md`, which gains a section on the external genes' columns.
- What changed with them, all within the declared scope (`check_scope.py`: 14 outputs
  differ from the previous commit, none outside the scope): E7's tables and column checks
  for the two genes (DDX3X gains 40 rows, TP53 loses its 15 NT rows); the TP53 table's
  record; `practical_metrics_external.csv` and the summary's sentence on the external
  genes (DDX3X's largest gain in sensitivity over the held-out seven is now 0.52 and its
  largest loss of specificity 0.06, both on its AlphaGenome rows, which are new; they were
  0.37 and 0.05, both on Pangolin); Supplementary Table S12; Figure 4 and its data; and
  the supplement manifest's records of the sources that changed.
- One change was outside the scope as first declared: the supplement manifest's S1 row
  records the hash of the TP53 table's manifest, which S1 reads, and that manifest's
  description of the 0.6.1 AlphaGenome column changed with the withdrawn entry. S1
  itself is byte for byte the same. `check_scope.py` reported it; by the author's
  decision the scope was amended to that one entry, and the manifest's declaration
  narrowed to the two fields that change (`scope_c3.json`, `amended_after_check`).

## After the thresholds were frozen: TP53's primary label is named for the published ranges it equals (stage E18, 2026-10-23)

Why: TP53's primary label was constructed by this study from the relative fitness
score (RFS > 0), its deposit, urn:mavedb:00001213-a-1, publishing no per-variant class.
The score set urn:mavedb:00001277-a-1 now carries investigator-provided ranges that
classify the same variants. Where they reproduce the label variant by variant, the
label is named for them (E18.4). One reason, one commit: only the label's name changes,
every number byte for byte the same (`scripts/stage_E18/scope_c4.json`, checked by
`check_scope.py`).

- The record (`scripts/stage_E18/fetch_tp53_ranges.py`,
  `phase1/data/evidence/external/tp53_raw/`): the score set's metadata, which carries
  the ranges, and its score table, as served from MaveDB (API 2026.3.2), each under
  25 MB, with a manifest of their sizes, sha256 and fetch times. Nothing else was
  fetched (the author's decision of 2026-10-23).
- The check (`scripts/stage_E18/check_tp53_ranges.py`, `range_check.json` and `.tsv`):
  - every one of the 288 TP53 rows E7 reads names one row of the score table;
  - mapped by E16's route (the atlas mapping module through the Mutalyzer model of
    NM_000546.6, anchored by the MANE Select partner's Ensembl lookup, the reference
    base checked against the cached Ensembl slice), each lands on its existing chrom,
    pos, ref and alt, 288 of 288;
  - classified by the investigator-provided ranges, read from the metadata and applied
    by E16's own rule (normal at or below 0, abnormal above), each gets the class
    `y_control_anchored` gives it, 288 of 288 (64 abnormal, 224 normal).

  E16's cross-check against MaveDB's mapped-variant records was not made: their file is
  40.4 MB, above the 25 MB limit. By the author's decision of 2026-10-23 the 288/288
  identical coordinates stand in for it.
- The name: "Published range (urn:mavedb:00001277-a-1)"
  (`evid_external.TP53_PRIMARY_LABEL_NAME`), in Supplementary Tables S1, S12 and S13 and
  in the TP53 table's manifest; Figures 1 and 4 carry its short form (below). The column
  keeps its name, `y_control_anchored`. DDX3X's control-anchored labels, sensitivity
  labels, keep their names. The median-split and mid-band labels stay as TP53's
  sensitivity labels. `test_evid_read_back.py` reads the printed name back through the
  same gene-aware mapping.
- The names come from one place in the code (`evid_external`): the full name for the
  tables, and for the figures the short form "published range", written as the figures
  write DDX3X's label ("deposit label", "deposit classification"), with the URN for the
  caption (the author's decision of 2026-10-23). The full name in Figure 1's row heading
  ran into panel a's plot; the short form does not.
- A case-insensitive search of every output for "control-anchored" found it left only in
  DDX3X's own labels and in internal column names, after one more place was changed: the
  TP53 table's manifest describes the in-frame subset's rule, and now names the label
  from the same mapping.
- `check_scope.py` first reported two omissions from the declared scope, both inside
  E18.4: the five new files of `tp53_raw/`, and the label's name in lower case in
  Supplementary Table S13's description of TP53's subset; the search above added a third,
  the manifest's rule. By the author's decision all three were added to `scope_c4.json`
  (`amended_after_check`). With the old name replaced by the new, every changed table is
  otherwise byte for byte the previous commit's, and no number changed.

## After the thresholds were frozen: the external test applies the in-sample thresholds of Table 3(b), BP4 included (stage E18, 2026-10-24)

Why: the external test carried over the medians of the leave-one-gene-out folds' PP3
thresholds, which Table 3(b) does not recommend. It now applies the thresholds Table
3(b) does recommend, the in-sample ones, on both sides (E18.5). The medians stay, as a
sensitivity analysis. One reason, one commit (`scripts/stage_E18/scope_c5.json`,
checked by `check_scope.py`).

- The thresholds: read at full precision from the frozen `evidence_thresholds.csv`,
  `pp3_threshold_insample` and `bp4_threshold_insample`, only for the tiers whose
  `*_reachable` is true: PP3 at Supporting, Moderate and Strong, BP4 at Supporting and
  Moderate. A tier not reached in-sample carries no threshold, as before. Nothing is
  fitted.
- The reading: PP3 on the upper side by `_band_detail`, the tier from the point
  estimate, with the tier at the interval's lower bound beside it, as before; BP4 on the
  lower side, the tier read by `evid_common.band_lr_benign_bound`, as the BP4 tier at
  Walker's published cut point is read. Each new cell's bootstrap draws on a seed of its
  own (`insample_<tier>`, `insample_bp4_<tier>`); the medians' cells keep theirs.
- The sensitivity analysis: the medians, computed as before, in columns of their own:
  `pp3_logo_median_<tier>_*` in E7's tables (`e3_<tier>_*` before), "<Tier>,
  sensitivity: …" in Supplementary Table S12. `check_scope.py` compares the renamed
  columns with the previous commit's: every cell the same.
- Walker's published cut points, 0.2 and 0.1: unchanged.
- The outputs that follow, each naming its thresholds' source in its column names:
  - `external_ddx3x.csv`, `external_tp53.csv`: new columns `pp3_insample_<tier>_*`
    and `bp4_insample_<tier>_*`;
  - Supplementary Table S12: "<Tier>: … (in-sample, Table 3(b))" and "BP4 <Tier>: …"
    columns before the sensitivity columns;
  - `practical_metrics_external.csv`: a `side` column; the in-sample PP3 and BP4 rows
    (`threshold_source` in_sample), the BP4 rows with E15's BP4 metrics, each set against
    the seven genes' in-sample row of the same threshold; the medians' rows unchanged;
  - `fig4_external.csv` and Figure 4: the ratio at the in-sample threshold in place of
    the ratio at the median;
  - `docs/practical-metrics-summary.md`, item 5: the sentence now reads the in-sample
    thresholds and names what it compares them with, the seven genes they were fitted
    on (before: the medians against the held-out genes).
- The supplement manifest's S12 row: its file hash, the hashes of the sources that
  changed, and its note, which now states the thresholds as above (it named the
  medians); the note was added to `scope_c5.json` before the scope was checked.
- `check_scope.py` against the previous commit: nine outputs differ, none outside the
  scope. In E7's tables and in S12 the medians' columns, under their new names, hold the
  previous commit's cells unchanged, and of every row only the new columns, and S12's
  Note, differ; in `fig4_external.csv` only the eight columns of the fitted threshold;
  in `practical_metrics_external.csv` the previous rows gain only the new columns, and
  480 rows are new. Figure 4's draft, declared, did not change: it draws the published
  cut point only.
- The sentence on the external genes: on DDX3X the in-sample thresholds are more
  sensitive than in the seven genes they were fitted on by as much as 0.35, and less
  specific by as much as 0.06 (it read 0.52 and 0.06 for the medians against the
  held-out genes).

## After the thresholds were frozen: the external genes' primary labels are named for their sources (stage E18, 2026-10-24)

Why: by the author's decision of 2026-10-24 each external gene's primary label is named
for its source as the source's record calls it, with the URN of the score set it comes
from. TP53's label is the classification by the calibration its record titles
"Investigator-provided functional ranges" (`investigatorProvided: true`); the same record
carries two other calibrations, IGVF's and ExCALIBR's, so "Published range", the name
E18.4 gave it, was imprecise and ambiguous. One reason, one commit, before E18.6: only
text changes, every number byte for byte the same (`scripts/stage_E18/scope_c5b.json`,
declared to the field and checked by `check_scope.py`).

- The names come from one mapping in `evid_external` (`PRIMARY_LABEL_SOURCE`): a label
  read from a score set's investigator-provided ranges is an "Investigator-provided
  range", one copied from a deposit's class column a "Deposit classification", each
  followed by the score set's URN:
  - TP53: "Investigator-provided range (urn:mavedb:00001277-a-1)";
  - DDX3X: "Deposit classification (urn:mavedb:00000658-0-1)", the class column of its
    classification score set.

  Their other labels, sensitivity labels, keep their names.
- Where the names are printed: Supplementary Tables S1, S12 and S13 (S13 names the label
  in lower case in each in-frame subset's description); the TP53 table's manifest, in its
  in-frame rule; and the range check's record, `label_name`. The range check ran again
  in the image: 288 of 288 rows covered, at the same coordinates, with the same label;
  `range_check.tsv` is byte for byte the same.
- Figures 1 and 4 write the short forms, "investigator-provided range" and "deposit
  classification" (Figure 1 wrote "deposit label" for DDX3X); the URNs go into the
  caption. At its old place Figure 1's heading for TP53 ran 62 px into panel a (of a
  2,161-px-wide PNG), and DDX3X's ended 1 px short of it. By the author's decision of
  2026-10-24 Figure 1's group headings now start further left, at axes x -0.70 instead
  of -0.56; no row, axis or panel moves. Figure 4's titles fit as they were.
- `check_scope.py` against the previous commit: ten outputs differ, none outside the
  scope. With each gene's primary-label rename applied to that label's rows
  (`rename_groups`, new in `check_scope.py`), Supplementary Tables S1, S12 and S13 are row
  for row the previous commit's; the supplement manifest changes in the three tables'
  file hashes, and the TP53 table's manifest and the range check's record each in the one
  field that names the label. The drafts of
  Figures 1 and 4, declared, did not change. A case-insensitive search of every output
  finds "published range" nowhere.

## After the thresholds were frozen: five more genes enter the external test (stage E18, 2026-10-24)

Why: the four genes of E16 (CTCF, SFPQ, RAD51D, XRCC2) and SBDS of E16b, scored by E18.2
with the analysis set's variables, are tested at the frozen thresholds like DDX3X and
TP53 (E18.6). One reason, one commit: only the new genes' rows and new tables
(`scripts/stage_E18/scope_c6.json`, checked by `check_scope.py`).

- The labels: the ones each gene's `<gene>_variants.parquet` carries, fixed in E16 and
  E16b before any predictor result. `newgenes_scored.parquet` copies those rows unchanged,
  and a test checks the labels against `<gene>_variants.parquet` value for value.
  - CTCF, SFPQ, RAD51D and XRCC2: the one calibration each record carries, titled
    "Investigator-provided functional classes" (`investigatorProvided: true`), named
    "Investigator-provided range (<URN>)".
  - SBDS: its record carries none; the score table's class column, named "Deposit
    classification (urn:mavedb:00001284-a-1)".

  The names come from the mapping of the previous commit (`PRIMARY_LABEL_SOURCE`).
- The test: E7's new stage (`--apply-new`, stage 12 of `scripts/reproduce_evidence.py`,
  after DDX3X's) runs `apply_thresholds` on each gene, in the three in-scope strata only,
  at least ten labelled on each side; the thresholds, readings and seeds as for DDX3X and
  TP53. Twelve columns are the analysis set's variables on every new gene, by
  `docs/external-column-identity.md` and by E7's column check; NT is named with its
  reason. A column that failed the check would stop the stage; none does.
- The training overlap: `predictor_training_provenance.csv` records none of the five
  assays in any predictor's training or model selection; the check, in the supplement
  stage, finds DDX3X x the combined Atlas score only, as `evid_figures.NOT_EXTERNAL`
  records, and would stop on a new gene.
- New tables:
  - `external_normal_only.csv`: the strata where only the normal variants can be
    evaluated (fewer than ten damaging, at least ten normal), under the primary labels of
    all seven external genes (the author's decision of 2026-10-24): at each frozen PP3
    threshold and, for the SpliceAI columns, at Walker's cut point, the normal variants at
    or above it out of all the normal variants; no ratio, no tier;
  - `external_summary.csv`: per gene set, stratum, column, side and tier, k/n (m) as
    Table 3 writes it: of the n external genes counted, k reach the tier at the frozen
    threshold, as Supplementary Table S12 prints it, and m are not evaluable (the
    author's decision of 2026-10-24, n counting every gene counted); once over the seven
    external genes, once over DDX3X and TP53, in one file. DDX3X x the combined Atlas
    score is not counted, and the note says so.
- Supplementary Table S12 gains the five genes' rows; its title and note in the manifest
  name the genes and their labels. Figure 4 is unchanged (B1).
- `scripts/reproduce_evidence.py` gains the stage; the stages after it are renumbered
  13 to 28 in its list.
- `check_scope.py` against the previous commit: fourteen outputs differ, none outside the
  scope. E7's ten new files hold the five genes' rows only (36 per gene, 13 in each
  column check); S12 gains 180 rows, and its DDX3X and TP53 rows are byte for byte the
  previous commit's; the manifest changes in S12's row alone. The normal-only table has
  97 rows, in four strata: SBDS 11–50, XRCC2 3–10 and 11–50, TP53 11–12. The external
  summary has 195 rows per gene set.

## Supplementary analyses: the five new genes in the practical metrics (stage E19, 2026-10-25)

Why: E18 tested CTCF, SFPQ, RAD51D, XRCC2 and SBDS at the frozen thresholds; E15's
table of practical metrics for the external genes still held DDX3X and TP53 only
(E19.1). This commit adds rows and changes nothing else (`scripts/stage_E19/scope_c1.json`,
checked by `scripts/stage_E19/check_scope.py`, which runs E18's check as it stands).

- `practical_metrics_external.csv`: the five genes' rows, built by E15's own external-gene
  code exactly as DDX3X's and TP53's are. Each gene is under its primary label
  (`evid_external.PRIMARY_LABEL`), in the three in-scope strata, for the twelve columns of
  the analysis set's variable (NT has none on these genes). The thresholds are the ones
  E7 carried over: Table 3(b)'s in-sample PP3 and BP4 thresholds and Walker's published
  cut point; by the author's decision of 2026-10-25 the medians of the leave-one-gene-out
  fits are carried too, as the sensitivity rows they are for DDX3X and TP53.
- `check_scope.py`: one output differs, 1,455 rows added, no existing row changed. The
  summary page `docs/practical-metrics-summary.md` is unchanged (E19.0 rule 1); its
  sentence naming DDX3X and TP53 as the external genes in this table is now incomplete,
  and is listed for the revision.
- `scripts/stage_E19/`: this stage's run records (`_runlog.py`, a copy of E18's), the
  image runner and the pipeline helper (copies of E18's for this folder), the scope check
  and the secrets scan (wrappers that run E18's scripts unchanged, recording here).

## Supplementary analyses, C2: ClinVar's coverage by position (stage E19, 2026-10-25)

Why: revision item C2 asks how ClinVar covers the variants the assays measured,
position by position, with every SNV of the analysis set as the denominator (E19.2). A
new stage, `evid_suppl_c2` (stage 29, after every existing stage), only adds files
under `phase1/reports/evidence/suppl_analyses/`; every existing output is unchanged
(`scripts/stage_E19/scope_c2.json`).

- The denominator: all 8,853 SNVs, the 1-2 nt dinucleotides included, at their signed
  intronic offset (-50 ... -1 acceptor side, +1 ... +50 donor side), read from the c.
  notation by the companion atlas's own offset pattern with the sign kept; its magnitude
  is the set's `intron_offset_abs` for every variant.
- `c2_clinvar_coverage_by_position.csv`: per gene and for the seven genes together, at
  each offset, the variants, those with any ClinVar record, those classified P/LP or
  B/LB, and P/LP, B/LB, VUS, Conflicting and Other (the set's `clnsig_class`), as
  counts and shares.
- `c2_clinvar_coverage_by_stratum.csv`: by band, unrecorded, recorded but unclassified
  and classified, as counts and shares. For the seven genes: 1-2 nt 856 variants,
  31.4% unrecorded; 3-10 nt 3,290, 51.6%; 11-50 nt 4,707, 76.6%.
- `c2_clinvar_coverage_draft.png`: a draft in evid_fig_data's manner; the figure's
  layout is left to the revision.

## Supplementary analyses, C3 items 1, 2 and 4: label sensitivity (stage E19, 2026-10-25)

Why: revision item C3 asks how far the thresholds rest on the labels' handling of
intermediate classes, on the label files the audit replaced, and on BARD1 and PALB2
(E19.5). A new stage, `evid_suppl_c3` (stage 30, after every existing stage), only adds
files under `phase1/reports/evidence/suppl_analyses/`; every existing output is
unchanged (`scripts/stage_E19/scope_c3_labels.json`). Item 3, a second BRCA2 dataset,
follows in commits of its own.

- Every refit is E3's in-sample fit, cell for cell: the same function, bootstrap and
  seed keys, on the bands 3-10, 11-50 and 3-50 nt, for the fourteen columns E3 fits
  (the fusion's out-of-fold column computed as E3 computes it; it is trained on the
  within-gene rank of the functional score, not on these labels). Before any refit, the
  same fit on the analysis set's own labels gives `evidence_thresholds.csv` back
  exactly, or the stage stops.
- `c3_intermediate.csv`: the 337 variants at 3-50 nt that the label files class as
  intermediate, indeterminate or enriched (181 at 3-10, 156 at 11-50 nt), all taken as
  damaging, then all as normal. Of 336 cells each, 8 gain a tier and 26 lose one under
  the first; 2 gain and 4 lose under the second.
- `c3_superseded_labels.csv`: PALB2's and BARD1's labels from
  `data/assay_labels/superseded/`, joined exactly as the current files are; 48 variants
  at 3-50 nt lose their label and 2 change side. 2 cells gain a tier, 1 loses one.
- `c3_without_bard1_palb2.csv`: not refitted; E11's fits of the subsets without BARD1,
  without PALB2 and without both (the PP3 side, Moderate and Strong, its four columns),
  beside every gene's. Of 24 cells each, 2, 1 and 2 lose their tier.

## Supplementary analyses, C4: prior sensitivity (stage E19, 2026-10-25)

Why: revision item C4 asks how the tiers depend on the prior of 0.10 (E19.6). A new
stage, `evid_suppl_c4` (stage 31, after every existing stage), only adds files under
`phase1/reports/evidence/suppl_analyses/`; every existing output is unchanged
(`scripts/stage_E19/scope_c4.json`).

- The priors: 0.05, 0.10 (the frozen fit) and 0.185, fixed in advance (E19.6). Each
  tier's edge keeps the posterior its frozen edge gives at 0.10:
  `LR_p = LR_0.10 x [0.10 / (1 - 0.10)] / [p / (1 - p)]`; the formula and the three
  priors' edges head both files. All four tiers on both sides are converted, Very
  strong included (the author's decision of 2026-10-25).
- The fit: E3's own in-sample function, bootstrap and seed keys, for every cell of
  `evidence_thresholds.csv` (fourteen columns, five strata, both sides, four tiers);
  the run at 0.10 gives that file back exactly, or the stage stops.
- `c4_prior_sensitivity.csv`: of 560 cells, 70 change tier at 0.05 (22 PP3 cells
  lost, 48 BP4 cells gained) and 67 at 0.185 (22 PP3 gained, 45 BP4 lost). Strong is
  reached in 29 column-stratum-side cells at 0.10; at 0.05 it is still reached in 19
  of them, at 0.185 in all 29.
- `c4_band_edges.csv`: the edges and the posteriors they keep.

## Supplementary analyses, C5: held-out intervals, two fold counts, jackknife (stage E19, 2026-10-25)

Why: revision item C5 asks for intervals on the held-out likelihood ratios, Table 3's
fold counts by the point estimate and by the lower bound, and how far the Strong
threshold moves from fold to fold (E19.7). A new stage, `evid_suppl_c5` (stage 32,
after every existing stage), only adds files under
`phase1/reports/evidence/suppl_analyses/`; every existing output, Table 3 among them,
is unchanged (`scripts/stage_E19/scope_c5.json`).

- `c5_heldout_lr_ci.csv`: the 1,120 folds of `evidence_thresholds_logo_folds.csv` at
  3-10, 11-50 and 3-50 nt (fourteen columns, four tiers, the PP3 side); 477 have a
  held-out ratio, each equal to the per-fold file's, and an interval by the resampling
  E7 uses on the external genes (`evid_external._band_detail`: the held-out gene's
  variants resampled within each class, 2,000 draws, seeded from the cell's identity).
  The fusion's rows carry E3's leakage note.
- `c5_fold_counts_point_vs_bound.csv`: Table 3's 78 cells. By the point estimate the
  counts are Table 3's in every cell (the stage stops otherwise); by the lower bound
  24 cells count fewer folds.
- `c5_jackknife_strong.csv`: per column and band, the Strong thresholds of the folds
  (number, smallest, median, largest) and the largest shift from the in-sample Strong
  threshold, read from the per-fold file only; Strong is reached in-sample in 18 of
  the 42 cells.
- The variant-level and gene-level bootstrap thresholds already stand side by side in
  Supplementary Tables S6 and S7; no file repeats them.

## Supplementary analyses, C8, C9, C10: read from existing outputs (stage E19, 2026-10-25)

Why: revision items C8 to C10 (E19.9). A new stage, `evid_suppl_c8_c10` (stage 33,
after every existing stage), reads existing outputs only, fits nothing, and only adds
files under `phase1/reports/evidence/suppl_analyses/`; every existing output is
unchanged (`scripts/stage_E19/scope_c8_c10.json`).

- `c8_without_bap1.csv`: the 11-50 nt band's labelled variants by gene (BAP1 holds 137
  of its 223 damaging variants); from `threshold_dilution.csv`, the fit on every gene
  but BAP1 at 11-50 nt (five genes, BRCA2 having none there) and 3-50 nt (six), with
  its held-out result on BAP1, beside every gene's (E11's four columns, the PP3 side,
  Moderate and Strong). Of 16 cells, 2 change whether the tier is reached; on BAP1, 14
  are evaluable and 10 clear the cut.
- `c9_family_holdout.csv`: (a) from `threshold_dilution.csv`, the fit on the homologous
  recombination genes (the four with variants at 11-50 nt) held out on VHL and BAP1 (48
  rows, 27 evaluable, 23 clearing the cut), and the fit on VHL and BAP1 held out on each
  of those genes (112 rows, 58 evaluable, 58 clearing); (b) the external genes, RAD51C's
  paralogues RAD51D and XRCC2 against the other five, counted from Supplementary Table
  S12 by the supplement's own `external_summary` (Table 3's k/n (m)), beside the seven
  external genes' row of `external_summary.csv`; the two groups add up to it.
- `c10_confusion_mcc_f1.csv`: TP, FP, FN and TN from `practical_metrics_by_threshold.csv`
  (PP3: damaging positive, at or above the threshold called; BP4: normal positive, at or
  below called), precision, recall, F1 and MCC, the formulas in the file's header, for
  its 336 in-sample and 336 held-out rows (133 and 102 with an MCC). For the supplement
  only.

## Supplementary analyses, E19.10: the thresholds refitted on twelve genes (stage E19, 2026-10-25)

Why: E19.10 asks how the in-sample thresholds move when the five further SGE genes of
the same platform join the seven. A new stage, `evid_suppl_twelve_gene` (stage 34, after
every existing stage), only adds a file under `phase1/reports/evidence/suppl_analyses/`;
every existing output, every frozen file among them, is unchanged
(`scripts/stage_E19/scope_twelve_gene.json`).

- The genes: the seven training genes and CTCF, SFPQ, RAD51D, XRCC2 and SBDS (DDX3X and
  TP53 left out: other cell lines and label mechanisms). The five genes' rows are E18's
  `newgenes_scored.parquet`, all at 3-50 nt, read as E7 reads them, under their primary
  labels; they follow the analysis set's rows, in E7's gene order.
- The columns: twelve, the panel without NT (no column of the analysis set's variable
  on the new genes); the fusion is not refitted (a trained model, with NT among its
  inputs). Both have rows that say so.
- The fit: E3's in-sample function, bootstrap and seed keys, the frozen fit's, on 3-10,
  11-50 and 3-50 nt (twelve genes, eleven at 11-50 nt); no leave-one-gene-out fit. The
  same code on the seven genes alone gives `evidence_thresholds.csv` back exactly.
- `twelve_gene_refit.csv`: of 288 cells, 29 change tier (PP3: 5 gained, 8 lost; BP4: 16
  lost); where both fits reach the tier (91 cells) the threshold moves by a median of
  14.0% of the seven-gene value.

## Supplementary analyses, C6 and C7: BP4 at the BP7 positions; the cut points by position (stage E19, 2026-10-26)

Why: revision items C6 and C7 (E19.8). A new stage, `evid_suppl_c6_c7` (stage 35, after
every existing stage), only adds files under `phase1/reports/evidence/suppl_analyses/`;
every existing output is unchanged (`scripts/stage_E19/scope_c6_c7.json`).

- `c6_bp4_bp7.csv`: at Walker 2023's BP7 positions (donor +7 and beyond, acceptor -21
  and beyond), the SGE damaging and normal variants whose published-basis SpliceAI
  score is (a) at or below the published cut, 0.1, or (b) at or below the BP4 Supporting
  threshold of Table 3(b) for the variant's band (+7 ... +10 in 3-10 nt, the rest in
  11-50 nt), per gene and for the seven genes. For (b) the threshold is the fitted
  value at full precision from `evidence_thresholds.csv`, read as E18's external test
  reads Table 3(b) (the author's decision of 2026-10-26); the file's header says that
  with Table 3(b)'s printed values (rounded down to three figures) the two SGE normal
  variants sitting exactly on the fitted thresholds, BAP1 c.255+9G>C and BARD1
  c.1810+43A>G, would not be counted, the damaging counts being the same. For the seven
  genes: (a) 3,658 variants, 57 damaging; (b) 502 variants, 5 damaging.
- `c7_by_position.csv`: per position group (+3...+6, +7...+10, +11...+50; -3...-10,
  -11...-20, -21...-50), the damaging and normal variants and E2's band ratios at the
  published cut points (at or above 0.2, at or below 0.1) with E2's gene-clustered
  intervals; every group has at least ten of each.

## Supplementary analyses, C1: the predictors against ClinVar (stage E19, 2026-10-26)

Why: revision item C1 asks how the frozen thresholds and the published cut points read
against ClinVar's germline classifications, as a concordance check (E19.3). A new stage,
`evid_suppl_c1` (stage 36, after every existing stage), only adds files under
`phase1/reports/evidence/suppl_analyses/`, and one input folder; every existing output
is unchanged (`scripts/stage_E19/scope_c1_main.json`).

- The reference: the analysis set's `clnsig_class`, the project's collapse of the VCF's
  INFO CLNSIG ("Aggregate germline classification for this single variant", ClinVar
  2026-06-15), P/LP positive and B/LB negative; conflicting classifications and the
  rest carry no reference label.
- The stars: from CLNREVSTAT ("review status of germline classification for the
  Variation ID") by ClinVar's own table for germline classifications on aggregate records
  (VCV and RCV), read from its review-status page (`phase1/data/evidence/e19/
  clinvar_review_status/`). The scripted fetch was served a reCAPTCHA page, kept out
  (its run record is kept); the author saved the page with Chrome on 2026-10-26, and its
  manifest records the URL, size, sha256, save time and that. Every review status of the
  set is graded; the conversion is written out (`c1_clinvar_star_conversion.csv`).
- `c1_clinvar_reference_metrics.csv`: at least one star, 74 P/LP and 645 B/LB at 3-10 nt,
  3 and 765 at 11-50 nt (counts only), 77 and 1,410 at 3-50 nt; at least two stars, 59
  and 273, 0 and 304, 59 and 577. Sensitivity and specificity with E15's Wilson
  intervals at Walker's cut points (both SpliceAI columns) and Table 3(b)'s in-sample
  thresholds, by E15's own code and names; AUROC by E4's estimator, each gene's AUROC
  beside the pooled one. Only BRCA1 and BRCA2 have ten of each at 3-10 and 3-50 nt;
  where BRCA2's AUROC is 1 (nine of the fourteen columns) E4's logit pooling is
  undefined and the row says so. No likelihood ratio, no tier (E19.0 rule 3).
- `c1_clinvar_sge_crosstab.csv`: ClinVar (P/LP, B/LB, unrecorded) x SGE (damaging,
  normal) x the side of each threshold, as counts; in P/LP the share of SGE damaging at
  or above each PP3 threshold, in B/LB the share of SGE normal at or below each BP4
  threshold (E15's `bp4_sensitivity`, the author's decision of 2026-10-25), both in the
  unrecorded group.

## Supplementary analyses, C1: the low-circularity subset (stage E19, 2026-10-26)

Why: revision item C1 asks the same ClinVar check on assertions that could not have used
the SGE data (E19.4). A new stage, `evid_suppl_c1_low` (stage 37, after every existing
stage), only adds files under `phase1/reports/evidence/suppl_analyses/`, with its inputs
under `phase1/data/evidence/e19/`; every existing output is unchanged
(`scripts/stage_E19/scope_c1_low.json`).

- The dates of assertions: LastEvaluated in ClinVar's monthly archive
  `variant_summary_2026-06.txt.gz` (439,287,089 bytes, sha256 6da78082...; kept in the
  gitignored `phase1/data/e19_cache/`, not committed and not fetched by CI); committed
  in its place is the extract of its GRCh38 rows for the analysis set's 3,279
  VariationIDs (all found), written in the reference image, with the extraction script.
- The archive against the VCF: classification string and review status, each
  normalised (underscores as spaces, white space collapsed, the VCF's empty field and
  the archive's "-" both absent), must be identical (the author's decision of
  2026-10-25): 3,278 of 3,279 agree (literally, before the normalisation, 82.04%); the
  one that differs is left out.
- The first public dates, each the earliest of its sources, every record saved as
  served: MaveDB's earliest version of the cited score set, the journal's online date
  (Crossref), and a preprint's date, the earlier of Crossref's posted and created (the
  author's review of 2026-10-25). Every gene's study was searched for a preprint in
  Europe PMC's preprint index (`scripts/stage_E19/fetch_preprints.py`), and a preprint
  is the study's when its title names the gene and saturation genome editing is named
  and its authors include the paper's first and last authors (for BARD1 and PALB2, which
  record no paper, a contributor named on the score set): found for BRCA1 (294520),
  BRCA2 (Huang's 10.1101/2023.12.14.571597), VHL (542698) and BARD1 (medRxiv, later than
  MaveDB's date); none for BAP1, PALB2 and RAD51C. BRCA2's date is thereby 2023-12-15,
  no longer the journal's 2025-01-08.
- The subset: P/LP or B/LB, at least one star, the two files agreeing, LastEvaluated
  strictly before the gene's date: 907 variants; 12 P/LP and 305 B/LB at 3-10 nt, 2 and
  348 at 11-50 nt. E19.3's metrics on it (`c1_low_circularity_metrics.csv`, by C1's own
  functions): the 11-50 nt band gives counts only, and no band has two genes with ten
  of each for an AUROC.
- Its limits, as E19.4 states them: some ClinVar classifications may have used predictor
  scores as PP3 or BP4, which no date filter removes; the distal band holds almost no
  usable P/LP.

## Supplementary analyses, C3 item 3: the second BRCA2 set's records, label rule and mapping (stage E19, 2026-10-26)

Why: revision item C3 asks whether BRCA2's labels hold against a second SGE dataset
(E19.5 item 3), and the author specified how to find and label one (2026-10-25). This
commit adds the records, the screening, the rule and the mapping, and no predictor score
is read: the analysis comes in the next commit (the author's order). Every existing
output is unchanged (`scripts/stage_E19/scope_c3_brca2_records.json`).

- The search (`scripts/stage_E19/fetch_brca2_second_set.py`): every one of the 2,938
  published score sets MaveDB served on 2026-10-25, read in URN-prefix parts (its
  paged search repeats and skips score sets that share an experiment title), the rule
  applied to each: the title, a target's name or a target's identifiers (Ensembl,
  RefSeq, HGNC, from saved HGNC, Ensembl and NCBI records) pointing to BRCA2. Ten score
  sets: the catalog's nine and urn:mavedb:00001279-a-1.
- The screening (`screen_brca2_second_set.py`, in the reference image): (a) SGE by the
  catalog's rule, (b) another laboratory by the papers' authors (Crossref), (c) shared
  variants with the set's 551 BRCA2 variants at 3-10 nt, (d) a per-variant experimental
  classification by R1, R2, R4, R6. Only urn:mavedb:00001242-a-1 (Sahu et al. 2025)
  meets all four; (c) 198 shared by its MaveDB table. The others fail as the author
  expected. 00001279-a-1 (the author's item 7, recorded only): no investigator-provided
  calibration; its scores equal 00001225-a-1's on none of the 6,959 variants.
- R4 does not apply (`paper_rule.json`): the paper's classes are cut on a mixture model's
  probability, which the score table does not carry; its quotes are checked against the
  publisher's PDF, which the author obtained and which stays in the gitignored cache
  (`check_brca2_paper_quotes.py`: 8/8 found). R6: Supplementary Table 5's
  `Classification`, never its ACMG column, each class by its side; the table, fetched
  through the article page's own entry, stays in the cache sha256-pinned and
  `--fetch-inputs` fetches it again, its compact extract committed
  (`supplementary_table5_classes.tsv`); grouped as the paper's text counts them, its
  classes are the paper's counts in every class.
- Table 5 checked row by row (`phase1/src/evid_second_set.py`): each row's c., mapped as
  E16 maps coordinates, against the row's g. On 6,519 rows they agree; 28 rows write a c.
  that does not give their g., the g. being GRCh38's (the exon 25 donor site written
  c.9502-c.9505 for c.9501+1 ... +4, so c.9505A>C/G/T stands twice in the table; the exon
  23 acceptor site written c.8953-4 ... c.8953 for c.8954-5 ... -1; one malformed
  c.8633-2-1A>T); 4 rows (c.8610A>C/T, c.8613A>C/T) write a reference base GRCh38 does
  not have. By the author's decisions of 2026-10-26 the rows are joined by g., the 28
  kept and the 4 dropped; any other disagreement stops.
- `coverage_3_10.csv`: the set's 551 variants at 3-10 nt by side and offset: 201 in
  Table 5 by g. (the 198 of MaveDB's table and c.9501+4A>C/G/T), only at donor +3 ... +9
  and acceptor -3 ... -7.
- `evid_label_sources.supplementary_rows` takes an optional g. column (the label audit's
  extracts are unchanged); `reproduce_evidence.py` fetches and checks the second set's
  cached table with the label audit's, and `--check` checks E19's records.

## Supplementary analyses, C3 item 3: BRCA2 against the second set (stage E19, 2026-10-26)

Why: E19.5 item 3, on the records, rule and mapping committed before it. A new stage,
`evid_suppl_c3_brca2` (stage 38, after every existing stage), only adds files under
`phase1/reports/evidence/suppl_analyses/`; every existing output is unchanged
(`scripts/stage_E19/scope_c3_brca2.json`).

- The classes: Sahu et al. 2025's Supplementary Table 5 (R6), joined to the set's BRCA2
  variants by the genomic coordinate after the row-by-row check; BRCA2's own classes are
  its label file's. Shared at 3-50 nt: 201 (all at 3-10 nt).
- `c3_brca2_concordance.csv`: the 3 x 3 table of the two classes, counts and shares; the
  two agree on 151 of the 201 (abnormal 24, normal 127); 25 that BRCA2's file calls
  normal the second set calls abnormal, 13 the other way.
- `c3_brca2_second_set.csv`: two refits, E3's in-sample fit beside the frozen one, in
  `c3_intermediate.csv`'s form: 'only variants the two sets agree on' (a variant they
  class differently loses its label, one side intermediate counting as different) and
  'the second set's labels' (every variant the second set classes takes its class);
  variants it does not cover keep their own. The 11-50 nt band holds no BRCA2 variant and
  cannot change. The first refit changes no tier; the second changes 5 of 336 cells, all
  at 3-10 nt (SpliceAI and its published basis lose PP3 Strong; three AVI BP4 cells).
- Both files' headers name the second set, its range (donor +3 ... +9, acceptor -3 ...
  -7), the rows of Table 5 the author's decisions concern, MaveDB's table lacking the six
  c.9505 rows, and the shared counts (201 by Table 5; the screening's 198 by MaveDB's
  table unchanged).

## Supplementary analyses, C12: the true positives' event records (stage E20, 2026-10-26)

Why: E20.1 (修改单 v2.40 §4, C12). E6 counted how many of the tools' false positives are
predicted to change the transcript in frame, with no true-positive group beside them. This
adds the true positives' event records, the input of the attribution that follows. No stage
is added, and every existing output is unchanged (`scripts/stage_E20/scope_c12_events.json`).

- The subsets: E6's and E2.7's definition with "called damaging" for "called normal", on the
  3-10 and 11-50 nt bands. The seven genes: 553 (SpliceAI, `spliceai_walker` >= 0.2) and 539
  (Pangolin, the same cut on `pangolin`); DDX3X 66 (damaging under any of its four label
  definitions); TP53 15. The same code with "called normal" gives E2.7's stored DDX3X and TP53
  subsets back exactly and contains E6's two; under the current labels E6's would gain 6 and 7
  PALB2 variants (E6 drew them on 2026-09-17, before E16d's relabel).
- The event records: `evid_score_spliceai_events.py` and `evid_score_pangolin_events.py`,
  unchanged, in the atlas's SpliceAI environment and the Pangolin environment E18 names, with
  E6's settings (SpliceAI distance 4,999, no masking; Pangolin -d 50, the CLI's masking, five
  chunks of two threads). Beside each event file: the scorer's provenance and an environment
  record (the command, the interpreter and its pip freeze, the references' sha256).
- Anchor 1: every row's largest score stands in the relation the existing records establish.
  SpliceAI equals `spliceai_walker` bit for bit (E6 609/609, DDX3X 333/333, TP53 28/28; new
  563, 66, 15); Pangolin equals the stored `pangolin` rounded as its CLI prints it (E6 524/524;
  new 549). Anchor 2: 10 SpliceAI and 10 Pangolin variants from E6's event files, scored in the
  same runs, give E6's records back in every field, and Pangolin's CLI lines verbatim.

## Supplementary analyses, C12: the true positives' in-frame attribution (stage E20, 2026-10-26)

Why: E20.2, on the event records of E20.1. A new stage, `evid_suppl_c12` (stage 39, after
every existing stage), only adds a file under `phase1/reports/evidence/suppl_analyses/`;
every existing output is unchanged (`scripts/stage_E20/scope_c12.json`).

- The true positives: each variant of E20.1's subsets, the anchors left out, classified by
  `evid_inframe.classify` (SpliceAI) or `classify_pangolin` (Pangolin), unchanged, on
  `evid_inframe`'s transcript maps; DDX3X and TP53 on SpliceAI's events, each label
  definition counted on its own.
- The false positives: E6's and E2.7's per-variant attribution files as they stand, counted
  as their summaries count them; every count equals the summary's, cell by cell. E2.7's
  DDX3X Pangolin rows have no true-positive counterpart and are not carried.
- `c12_inframe_tp_fp.csv`: 28 rows (gene set x label definition x tool x band x group), the
  in-frame share with E15's Wilson interval and the share among the resolved. Seven genes,
  3-10 nt: in frame 135 of 453 true positives and 155 of 387 false positives (SpliceAI), 145
  of 451 and 131 of 348 (Pangolin).
- The header states that the frame is inferred from the tools' output, that TP53's label is
  its primary label variant by variant, the bands, and that E6's false positives predate
  E16d's relabel (6 and 7 PALB2 variants would be added under the current labels).

## Supplementary analyses, C13: the VCEP specifications' records and rules (stage E20, 2026-10-26)

Why: E20.3 (修改单 v2.40 §4, C13). The VCEP specifications of the training genes, as records,
before any predictor score is read. No stage is added, and every existing output is
unchanged (`scripts/stage_E20/scope_c13_records.json`).

- The records (`phase1/data/evidence/vcep/raw/`, with `manifest.json`): ClinGen's Criteria
  Specification Registry's list of specifications; for each training gene, every
  specification whose rule set names it, with its document and full record; their Zenodo
  deposits' records (licence cc-pddc); ClinGen's terms of use (curated content under CC0 1.0).
  Current and approved: ENIGMA GN092 v1.2 (BRCA1) and GN097 v1.2 (BRCA2), the Hereditary
  Breast, Ovarian and Pancreatic Cancer panel's GN077 v1.3 (PALB2), the VHL panel's GN078
  v1.1 (VHL). None for BAP1 and BARD1; RAD51C's GN096 v1.0 is "Classification Rules In Prep".
- `vcep_rules.json`: 12 rules (PP3, BP4, BP7 of the four specifications), each quote checked
  verbatim against the field of the record it cites; the strengths each allows are the
  record's (Supporting for every PP3 and BP4). SpliceAI cut points: ENIGMA and PALB2 PP3 >=
  0.2, BP4 <= 0.1; VHL PP3 > 0.5 and BP4 <= 0.1, each with VarSeak. PALB2's specification
  names the SVI's thresholds (Walker et al. 2023); the others state no SpliceAI setting.
- `reproduce_evidence.py --check` also reads the VCEP manifest.

## Supplementary analyses, C13: the VCEP rules on this study's data (stage E20, 2026-10-26)

Why: E20.4, on the rules committed before it. A new stage, `evid_suppl_c13` (stage 40, after
every existing stage), only adds a file under `phase1/reports/evidence/suppl_analyses/`;
every existing output is unchanged (`scripts/stage_E20/scope_c13.json`).

- Each PP3 and BP4 rule of `vcep_rules.json`, within its gene, on the 3-10, 11-50 and 3-50
  nt bands (every rule's position scope covers them): E2's band ratio, gates and
  zero-event rule (`evid_walker_thresholds._band_stats`), the external genes' within-gene
  interval (`evid_external._band_detail`, 2,000 draws), the tier at the point estimate and
  at the conservative end, beside the strengths the specification allows (Supporting for
  every PP3 and BP4). The score column is the published basis (`spliceai_walker`) for every
  rule; VHL's "> 0.5" selects the same variants as ">= 0.5" here. VarSeak, PhyloP, RNA
  evidence and PVS1 are not evaluated.
- `c13_vcep_band_lr.csv`: 24 band rows (18 evaluable; BRCA1 and BRCA2 at 11-50 nt and VHL's
  PP3 at 3-10 and 11-50 nt are not), 12 BP7 rows read from `c6_bp4_bp7.csv` beside each
  specification's own BP7 scope, and BAP1, BARD1 and RAD51C without a specification.

## Supplementary analyses, E6 under the current labels (stage E20b, 2026-10-26)

Why: E20b.1 (the author's answer to E20's open items 1 and 2). E6 drew its false-positive
subsets on 2026-09-17, before E16d switched PALB2 and BARD1 to the archived classification;
under the current labels they gain 6 (SpliceAI) and 7 (Pangolin) PALB2 variants at 3-10 nt
and lose none. E6 is drawn again with its own code, its events are re-scored whole, and its
attribution, Supplementary Table S13 and C12's false-positive group follow. No stage is
added; every other output is unchanged (`scripts/stage_E20b/scope_e6.json`).

- The pre-check: every committed input that a one-off script derived from label columns
  (E6's subsets and the Pangolin scorer's variant list, DDX3X's and TP53's in-frame subsets,
  E20.1's true-positive subsets and score inputs, the four further genes' scored table),
  derived again by the same code under the current labels and compared variant by variant.
  Only E6's two subsets and the Pangolin variant list (the same 7 variants) differ.
- The subsets: `python -m src.evid_inframe --subset`, unchanged. SpliceAI 615 (was 609),
  Pangolin 531 (was 524): exactly the committed rows and the 6 and 7 variants.
- The events: re-scored whole by `evid_score_spliceai_events.py` and
  `evid_score_pangolin_events.py` in E20.1's environments and settings (SpliceAI distance
  4,999, no masking; Pangolin -d 50, the CLI's masking, five chunks of two threads). Beside
  each event file, the environment record E6 did not keep (the command, the interpreter and
  its pip freeze, the references' sha256). Anchor 1: the committed variants' records (609
  and 524) are the same in every field, and Pangolin's CLI lines verbatim. Anchor 2: every
  row's largest score stands in E20.1's relation to the stored column (615 and 531).
- `inframe_attribution.csv`: the committed rows unchanged; 13 rows added (SpliceAI: 4 in
  frame, 1 out of frame, 1 undetermined; Pangolin: 3 in frame, 4 out of frame). The summary
  changes only at 3-10 nt: SpliceAI 159 of 393 in frame (was 155 of 387), Pangolin 134 of
  355 (was 131 of 348). Supplementary Table S13's two rows and the supplement's manifest
  follow, and so do C12's false-positive rows of the seven genes at 3-10 nt; C12's header
  says the subsets were drawn again.
- E2.7's DDX3X and TP53 files, C12's true-positive rows and C13 are unchanged. E20.1's
  anchor test reads E6's event files without the added variants, as they stood when the
  anchors were drawn.

## Supplementary analyses, C13: ENIGMA's SpliceAI setting from its specification (stage E20b, 2026-10-26)

Why: E20b.2 (E20's open item 4). ENIGMA's criteria in the registry name SpliceAI and its cut
points and refer to "Specifications Figure 1A and Appendix J", which the registry's records
do not carry. No predictor score is read, no stage is added and every existing output is
unchanged (`scripts/stage_E20b/scope_enigma.json`).

- Where: the registry's page for each specification (GN092, BRCA1; GN097, BRCA2) lists the
  files attached to its rule set, version 1.2 like the records. The two that hold Figure 1A
  and Appendix J (`Specifications_V1.2_2024-11-18.docx`, `Appendices_V1.2_2024-11-18.docx`)
  are the same bytes for both genes. ClinGen's page for the panel links the registry and no
  file. The pages and Crossref's record of the panel's paper are committed in `vcep/raw/`;
  the files themselves are kept in the gitignored `phase1/data/e20b_cache/` and recorded by
  URL and sha256 in `raw/manifest.json`: ClinGen puts its curated content under CC0, but the
  journal's copy of the specification is under the publisher's licences, none of them open.
  PMC served a challenge page, which was not answered.
- `raw/enigma_quote_check.json`: seven sentences of Appendix J, each found verbatim in both
  copies; Figure 1A's caption, its two images (found through the document's relationships)
  and the key inside the first, transcribed (it names SpliceAI and its cut points only);
  every Appendix J paragraph on a source, version, annotation, distance, masking or score,
  and every paragraph of either file that names SpliceAI.
- The setting Appendix J states: the largest of the four delta scores, Max Distance 10 000,
  'unmasked' scores; SpliceAI is cited by its paper, with no version, online tool or
  annotation. `vcep_rules.json`: ENIGMA's six rules record this in `score_definition_stated`
  and carry the sentences in `further_quotes` after the registry's; nothing else changes.
  Neither SpliceAI column of this study has that distance (published basis 4,999,
  `spliceai` 50).

## Supplementary analyses, C13: BP7 as each specification defines it (stage E20b, 2026-10-26)

Why: E20b.3 (E20's open items 5 and 6). E20.4 read every specification's BP7 from
`c6_bp4_bp7.csv`, at Walker 2023's positions; PALB2's BP7 starts at +8 and -22 and carries no
prediction condition, and VHL's names no position and asks for PhyloP <= 0.2. The
definitions are the specifications' own text (E20b.3), not chosen from the counts. Only
`c13_vcep_band_lr.csv` changes (`scripts/stage_E20b/scope_c13_fix.json`).

- BP7 is computed in `evid_suppl_c13`, counted as C6 counts (`evid_suppl_c6_c7.labelled`:
  the SGE labels, the signed offset, `spliceai_walker`). ENIGMA: donor +7 and beyond,
  acceptor -21 and beyond, with BP4 (SpliceAI <= 0.1), per side and both; the same
  definition as C6's condition (a), and the stage stops unless its rows equal C6's cell for
  cell. PALB2: donor +8 and beyond, acceptor -22 and beyond; (i) with no prediction
  condition, as the criterion reads, and (ii) with BP4 (SpliceAI <= 0.1), as the
  instructions to use allow. VHL: the study's 3-10, 11-50 and 3-50 nt, with BP4's SpliceAI
  half (<= 0.1) and PhyloP <= 0.2 (the analysis set's `phylop`, UCSC phyloP100way, hg38; no
  variant lacks it). Each row gives the variants at the positions (or in the band), those
  meeting the condition, the damaging and normal among them and the damaging share.
- PALB2: 5 of 409 at the positions are damaging; all 5 meet (ii), among 368. VHL: none of
  the 43 variants meeting the condition at 3-50 nt is damaging.
- The band rows: every count, ratio, interval and tier as before. A band with no labelled
  variant (BRCA2 at 11-50 nt) now says so in its reason; ENIGMA's rows name the setting its
  Appendix J states (E20b.2), with the published basis kept; the header records the BP7
  definitions, PhyloP's source and ENIGMA's setting.

## The quick check reads a cached file where its manifest caches it (stage E20b, 2026-10-26)

Why: E20b.4 (E20's open item 3). `reproduce_evidence.py --check` read the top-level "file" of
E19's ClinVar archive manifest as a path in the repository; it is the name as served, the
archive being kept only in the gitignored cache, so the check reported it missing. Code and
tests only; no output changes (`scripts/stage_E20b/scope_check_fix.json`).

- A file that a manifest records with a `cached_at` (the cache path) and a sha256, wherever
  the entry sits in E19's and E20's manifests (`cached_entries`), is checked at its cache
  path when the cache holds it, and waited for, like the ClinVar release, when it does not:
  E19's ClinVar archive and E20b's four ENIGMA specification files.
- Tests: the archive read where the cache holds it and waited for where it does not (the
  cache simulated either way, whatever this machine's holds), and entries found wherever a
  manifest records them. `bash reproduce.sh --check` on this machine:
  `scripts/stage_E20b/run_log/run_check_*.json`.

## The manuscript numbers gate: three repairs and a cells mode (stage E5, 2026-10-27)

Why: E5.1 (修改单 v2.44 §3 B2). E5.0, the read-only check of the v5 texts against 72dbcb1
before the revision, found three holes in `scripts/check_evidence_manuscript_numbers.py`, and
that a number found somewhere in the file its sentence names may still not be the number
the sentence means. Code, tests and documents only; no output changes
(`scripts/stage_E5/scope_gate.json`).

- The self-test plants its numbers after the first sentence of the body that names a .csv
  file, after confirming that none of them is a declared constant or a value those files
  hold, and counts each planted number the bound mode reports: "caught k / n", PASSED only
  when k = n. It used to pass whenever the text had any number the bound mode could not
  find, whether or not the planted ones were caught; on v5 it caught five of six, 12.34
  matching a stored proportion of 0.12336 times one hundred. Now six of six, on the main
  text and on the supplementary notes.
- A stored proportion times one hundred is compared only with a number printed as a
  percentage ("%" or "per cent"), in the pool and the bound modes.
- A token is declared by its context first (`CONTEXT_RULES`), then by the plain list. "95"
  is the confidence level only beside "CI", "confidence" or a bound; v5's "in at least 95%
  of evaluable cases" is a measurement. E5.0's thirteen constants with a wrong reason get
  theirs: inline enumerations, the day of a date, the upper end of 3-10 nt, the scoring
  distance of 50, gnomAD's 0.1% allele frequency, the TP53 scale's -1 and +1. The entries
  "(1)" "(2)" "(3)", which no token ever matched, are gone. `--declared` lists each declared
  token with its reason.
- A text without "# Abstract" or "# References" (the supplementary notes) is read from the
  top, or to the end, so the gate runs on it as it is.
- Cells mode, `--cells TABLE --document main|supplement`: a binding table names, for every
  number of the text, either why it is a constant or the file, rows, column and derivation
  (value, sum, min, max, count, count_nonnull, all, text) of the one cell it is read from,
  and the scale (100 for a proportion printed as a percentage). Each measured number is
  compared with that cell at its printed precision; numbers that differ, rows that cannot
  be resolved (no file, row or column; an excerpt that is not the sentence's; no number at
  the row's place) and numbers with no row are listed, and any of the three fails the run.
- `docs/manuscript_bindings_v5.tsv`: the 182 numbers of v5 (main text 120, supplementary
  notes 62), from E5.0's number table, by `scripts/stage_E5/build_bindings_v5.py`; its
  `tests/fixtures/manuscript_bindings_v5_expected.tsv` holds E5.0's verdict on each.
- On v5 at 72dbcb1 the cells mode gives E5.0's verdicts number by number
  (`scripts/stage_E5/check_bindings_v5.py`): main text 53 constants, 48 matching, 19 not
  (E5.0's 14 that do not match and the 5 its bound mode could not find, each with E5.0's
  value); supplementary notes 23 constants, 31 matching, 8 not (E5.0's 3, and the 5 of its
  10 not found that have changed; the other 5 are E5.0 §2.3's four constants and the
  derived 71, now read as a sum). Nothing cannot be bound or is missing from the table.
  At 8e0d609, the release v5 cites, 104 of the 106 measured bindings round to their
  printed number; the other two cannot be read there (a column E18 renamed, a scoring
  record E17 added).
- The bound mode on v5 at 72dbcb1 after the repairs: main text 62 matched, 53 declared,
  5 not found, none unbound (61, 54, 5, 0 before: line 89's 95% is now measured, and its
  file holds a value that rounds to it); supplementary notes 32, 18, 10, 2.
- Tests: `tests/test_evidence_manuscript_numbers.py`, on small texts and tables made in the
  test, on the v5 table against 72dbcb1's files (skipped where that commit is not in the
  clone), and on the v5 texts when `EVID_V5_MAIN` and `EVID_V5_SUPPLEMENT` name them (they
  are not in the repository).

## The repository's documents brought up to date; --check on a cache-only file (stage E5, 2026-10-27)

Why: E5.2 (修改单 v2.44 §3 B2, from E5.0 §5). Documents and code only; no output changes
(`scripts/stage_E5/scope_docs.json`).

- README: 8,501 labelled variants, not 8,453 (`set_counts.csv`, since the labels were rebuilt
  from their cited sources, E16d); the five further external genes are scored and tested
  (E18), not "none of the five is scored yet"; 40 stages, not 25; the VHL assay paper,
  PMID 38969834, is Buckley 2024, as its MaveDB record gives it
  (`phase1/data/evidence/e19/first_public_dates/mavedb/urn_3Amavedb_3A00000675-a-1.json`),
  not 2023. `docs/milestone4_candidate_genes.md`: the same year.
- `reproduce_evidence.py --check`: a file a manifest records with a `cached_at` (E19's
  ClinVar archive, E20b's four ENIGMA specification files) was said to wait for the full
  run when the cache lacks it. A one-off script fetched each, no stage reads it and the
  full run does not fetch it, so the line now says that (`CACHE_ONLY`, `absent_line`;
  E20b's open item 8). An input the full run fetches still waits for it. Tests: both lines.
- Not changed: `supplementary/build_supplementary.py` line 68 calls the BAP1 deposit GRCh38;
  that table is Supplementary Table S1 of Study B, the repository's other analysis, not this one's. `docs/practical-metrics-summary.md`
  line 31 still names only DDX3X and TP53 in `practical_metrics_external.csv`, which holds
  all seven external genes since E19: the file is an output of stage E15, so changing the
  line changes an output, which this commit does not do.
