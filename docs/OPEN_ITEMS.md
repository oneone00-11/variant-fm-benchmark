# OPEN ITEMS (author to verify / decide)

Items below were **deliberately not auto-filled** because they require the author's
own knowledge, an external source, or a re-run. Nothing here was invented. Resolve
before submission to *Bioinformatics*. Companion: `docs/CHANGES.md`.

> **Updated at Milestone 11.** RESOLVED since Milestone 10: §B Supplementary Table S1
> (now `supplementary/Table_S1.tsv`/`.docx`); §C-3 NT scoring script (now archived at
> `scripts/92_score_nt.py`, reproducible on GPU). Still-open parts are kept below and
> marked. New verification items from Table S1 are folded into §A.

---

## A. References to verify or complete (DO NOT rely on as-is)

The revision added only what could be stated responsibly; every uncertain detail is
flagged in the manuscript reference list with "[… to verify]". Confirm each against
the primary source (publisher page / PubMed / DOI):

1. **AlphaGenome (Avsec et al.)** — listed as "[Preprint / DeepMind]". Per the task
   note it reportedly appeared in *Nature* around early 2026; **confirm the year,
   venue, volume/pages and author list** and replace the placeholder. This is the
   single most important citation to pin down (it is a headline comparator).
2. **GPN-MSA (Benegas et al., 2025, Nat. Biotechnol.)** — confirm volume/pages/DOI
   and exact author list (currently "[details to verify]").
3. **Buckley et al. (2024), VHL SGE** — confirm journal, year, volume/pages
   (PMID 38969834 per milestone4; verify the PMID actually corresponds).
4. **Waters et al. (2024), BAP1 SGE** — confirm journal/volume/pages
   (PMID 38969833 per milestone4; verify).
5. **gnomAD version vs citation (see also §C-2):** the text uses **gnomAD v4**, but
   the reference list currently has **Karczewski et al. 2020** (the v2 flagship).
   Either cite the gnomAD **v4** release paper (reportedly Chen et al., 2024, *Nature*
   — **verify before citing; not added automatically to avoid fabrication**) or add a
   sentence clarifying which gnomAD paper is the intended methodological reference.
6. **phyloP (Pollard et al., 2010, Genome Res. 20:110–121)** and **phastCons
   (Siepel et al., 2005, Genome Res. 15:1034–1050)** — added at the standard citation
   level; **verify volume/page numbers** (marked "[bibliographic details to verify]").
7. **Dalla-Torre et al. (Nucleotide Transformer)** — year given as 2024 (Nat. Methods);
   confirm the final published year/volume (some versions are 2024–2025).
8. General: confirm every remaining reference's exact author list, year, volume and
   pages; the draft used best-known forms, not journal-verified records.
9. **Supplementary Table S1 PMIDs (M11):** the per-gene PMIDs were taken from the
   project's milestone records and carry `[verify]` in the table — confirm BRCA1
   (30209399), BRCA2 (39779857), RAD51C (39299233), VHL (38969834), BAP1 (38969833)
   against MaveDB / the source papers. **BARD1 and PALB2 PMIDs are not recorded** in
   the project materials and were left as "PMID not recorded [verify]" (not invented);
   the author should supply them.

## B. Placeholders to fill (left blank on purpose)

- Author list and order; institutional affiliation(s).
- Corresponding-author e-mail (`[contact e-mail]`, appears in title block, Abstract
  Contact, and Data availability).
- Associate Editor (`[TBD]`).
- Funding sources (`Funding` section).
- GitHub URL — replace `[user]` in `https://github.com/[user]/variant-fm-benchmark`
  (title/Abstract Availability + Data availability).
- ~~**Supplementary Table S1** (per-gene MaveDB URN/transcript/licence table)~~
  **RESOLVED (M11):** produced at `supplementary/Table_S1.tsv` and `Table_S1.docx`.
  (PMID verification for it is tracked in §A-9.)
- Decide whether to also deposit the data (e.g. Zenodo DOI) and cite it.

## C. Decisions that may require re-running analysis (NOT executed; no results invented)

1. **Primary set: ClinVar∩functional intersection (21,410) vs full functional SNV
   set (~34,582).** (Task Part 8.) The current main analysis is restricted to the
   intersection so every variant has both a functional score and a ClinVar label. A
   neutral one-line rationale was added to Methods 2.1 and the extension flagged in
   the Discussion. **Decision for the author:** either (a) keep the intersection and
   let the added rationale stand, or (b) re-run scoring/evaluation on the full
   ~34,582 functionally-scored SNV set (higher power; pre-empts selection-bias
   concerns) and report that as primary. Option (b) requires re-running the pipeline
   — **not done here, and no such result was estimated or written.**
2. **gnomAD v4 reference (see §A-5):** if you choose to cite the v4 paper, verify and
   insert it; this is a citation fix, not a re-run.
3. **NT exact configuration (reproducibility gap):** PARTIALLY RESOLVED (M11). The
   original `score_nt2.py` was **not found anywhere** on the machine or in the repo
   (only the git-ignored cache `nt_cache2.tsv` survives), so the scoring flow was
   **rebuilt and committed** as `scripts/92_score_nt.py` (archival reconstruction;
   regenerates the cache on a GPU). The checkpoint
   (InstaDeepAI/nucleotide-transformer-v2-500m-multi-species) and the masked 6-mer
   reference-minus-alternate LLR are now stated in Methods 2.3 and in the script.
   **STILL OPEN — author to confirm:** the **exact input context window** used in the
   original cloud run (unrecorded; the rebuilt script uses a documented default
   `WINDOW_BP=6000` explicitly labelled "to be confirmed", and the manuscript does NOT
   assert a window number). Confirm the window actually used and update the script
   constant + (optionally) Methods 2.3; if the original value is truly unrecoverable,
   state the archived script's window as the reproducible reference and say so.

## D. Optional polish (content, not blocking)

- Consider adding the SpliceAI/Pangolin compute environment (SpliceAI 1.3.1 /
  TensorFlow 2.19.1 CPU; Pangolin / torch 2.12.0 CPU; from milestone6) to the
  Supplementary methods for full reproducibility.
- Consider one sentence in 2.4 stating the bootstrap resample count (2,000) is shared
  across all bootstrap analyses (already implied; currently stated per analysis).
