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
