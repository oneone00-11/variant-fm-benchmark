"""The further external genes (stages E16 and E16b): what their tables may claim.

Every number the stages write traces to a file MaveDB served or a reference the
repository pins: these tests hold the variant tables to the reference genome, the
labels to the deposits' own ranges (or, for SBDS, which has none, to the score
table's own class column) and the seven genes' rule for intermediate classes, the
counts to one another, and every manifest to the files on disk. E16b adds SBDS and
two checks that change nothing: the IGVF copy of the SBDS scores, and the BARD1 and
PALB2 label files against every archived version of their score sets.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
EXT = PHASE1 / "data/evidence/external"
REPORTS = PHASE1 / "reports/evidence"
CACHE = PHASE1 / "data/evidence/reference_cache"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
from src import config as C          # noqa: E402
from src import evid_common as K     # noqa: E402
from src import evid_newgenes as N   # noqa: E402

GENES = N.GENES                                                   # E16
RUN_OF = {g: run for run in N.RUNS.values() for g in run["genes"]}
ALL_GENES = list(RUN_OF)                                          # E16's four, then SBDS (E16b)
E16_COMMIT = "a35ec6d9d7be1574f58191495e71dd1995f7ab3b"           # where E16 was committed


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _manifest(gene: str) -> dict:
    return json.loads((EXT / f"{gene.lower()}_raw/manifest.json").read_text())


def _variants(gene: str) -> pd.DataFrame:
    return pd.read_parquet(EXT / f"{gene.lower()}_variants.parquet")


def _classes(gene: str) -> list[dict]:
    inv = json.loads((EXT / f"{gene.lower()}_raw/score_ranges.json").read_text())
    assert len(inv["investigator_provided"]) == 1
    return inv["investigator_provided"][0]["functionalClassifications"]


def _has_ranges(gene: str) -> bool:
    inv = json.loads((EXT / f"{gene.lower()}_raw/score_ranges.json").read_text())
    return bool(inv["investigator_provided"])


def _class_of(gene: str) -> dict:
    """Class name -> functional classification: the deposit's ranges, or, for a
    deposit without them, the class-name table its score-table column maps through."""
    if _has_ranges(gene):
        return {c["label"]: c["functionalClassification"] for c in _classes(gene)}
    return dict(N.CLASS_NAME_MAP)


def _prefix(gene: str) -> str:
    return RUN_OF[gene]["prefix"]


def _csv(name: str) -> list[dict]:
    with open(REPORTS / name, newline="") as f:
        return list(csv.DictReader(f))


def _raw_table(gene: str) -> list[dict]:
    man = _manifest(gene)
    path = REPO / next(e["path"] for e in man["files"] if e["path"].endswith(".csv"))
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _provenance(gene: str) -> dict:
    return json.loads((EXT / f"{gene.lower()}_variants.parquet.provenance.json").read_text())


# --------------------------------------------------------------------------
# reference bases: 100%, three ways
# --------------------------------------------------------------------------
@pytest.mark.parametrize("gene", ALL_GENES)
def test_the_recorded_reference_check_covers_every_variant_of_the_table(gene):
    """scripts/stage_E16/check_reference.py, against the Ensembl release-112
    chromosome, on exactly the table committed now."""
    rec = json.loads((EXT / f"{_prefix(gene)}_reference_check.json").read_text())
    g = {r["gene"]: r for r in rec["genes"]}[gene]
    table = EXT / f"{gene.lower()}_variants.parquet"
    ref = _manifest(gene)["reference_fasta"]
    assert g["variants_sha256"] == _sha(table)
    assert g["n"] == g["matches"] == len(_variants(gene)) and g["mismatches"] == 0
    assert g["fasta"]["sha256_fa"] == ref["sha256_fa"] and g["fasta"]["url"] == ref["url"]
    assert g["fasta"]["ensembl_release"] == 112


@pytest.mark.parametrize("gene", ALL_GENES)
def test_every_reference_base_matches_the_tracked_sequence_slice(gene):
    """The Ensembl slice the mapping read, tracked and registered in the cache."""
    item = _manifest(gene)["transcript_model"]["cache_files"]["region"]
    path = REPO / item["path"]
    registered = json.loads((CACHE / "MANIFEST.json").read_text())["files"]
    assert _sha(path) == item["sha256"] == registered[str(path.relative_to(CACHE))]["sha256"]
    _, _, chrom, start, end = path.stem.split("_")      # seq_GRCh38_chr<c>_<s>_<e>
    start, end, seq = int(start), int(end), path.read_text().strip().upper()
    assert len(seq) == end - start + 1
    v = _variants(gene)
    assert (v["chrom"].astype(str) == chrom.removeprefix("chr")).all()
    assert [seq[p - start] for p in v["pos"].astype(int)] == list(v["ref"])


@pytest.mark.parametrize("gene", ALL_GENES)
def test_every_reference_base_matches_the_chromosome_when_it_is_here(gene):
    """The whole chromosome, a source the mapping never read; fetched by
    reproduce_evidence.py --fetch-inputs (gitignored), so absent in a bare clone."""
    ref = _manifest(gene)["reference_fasta"]
    fa = REPO / ref["file_fa"]
    if not fa.exists():
        pytest.skip(f"{ref['file_fa']} is not here (reproduce_evidence.py --fetch-inputs)")
    assert _sha(fa) == ref["sha256_fa"], f"{ref['file_fa']} is not the release-112 file"
    pysam = pytest.importorskip("pysam")
    v = _variants(gene)
    with pysam.FastaFile(str(fa)) as f:
        got = [f.fetch(ref["chrom"], p - 1, p).upper() for p in v["pos"].astype(int)]
    assert got == list(v["ref"])


# --------------------------------------------------------------------------
# labels: the deposit's own ranges or class column, the seven genes' rule
# --------------------------------------------------------------------------
@pytest.mark.parametrize("gene", ALL_GENES)
def test_range_classes_are_only_the_deposits_own_and_unassigned(gene):
    allowed = set(_class_of(gene)) | {N.UNASSIGNED}
    assert set(_variants(gene)["range_class"]) <= allowed
    for kind in ("counts_by_offset", "controls"):
        name = f"{_prefix(gene)}_{kind}.csv"
        assert {r["range_class"] for r in _csv(name) if r["gene"] == gene} == allowed


@pytest.mark.parametrize("gene", ALL_GENES)
def test_each_range_class_has_one_label_and_it_follows_the_seven_genes_rule(gene):
    """abnormal -> 1, normal -> 0, any other class NA, like the intermediate classes
    of config.ASSAY_LABELS; unassigned is NA too."""
    fc = _class_of(gene)
    v = _variants(gene)
    for cls, labels in v.groupby("range_class")["label"]:
        values = set(labels.dropna())
        assert len(values) <= 1 and (values == set() or labels.notna().all()), cls
        want = {"abnormal": {1.0}, "normal": {0.0}}.get(fc.get(cls), set())
        assert values == want, (cls, values)


def test_the_label_rule_is_the_one_the_seven_genes_use():
    """Every map in ASSAY_LABELS sends its classes to damaging, normal or None (an
    intermediate class, excluded); E16's map, and SBDS's class-name table, have the
    same three outcomes and no other."""
    for gene, spec in C.ASSAY_LABELS.items():
        outcomes = set(spec["map"].values())
        assert outcomes <= {1.0, 0.0, None} and {1.0, 0.0} <= outcomes, gene
        assert None in outcomes, f"{gene}: no intermediate class mapped to None"
    assert N.LABEL_OF == {"abnormal": 1.0, "normal": 0.0}
    for gene in GENES:
        assert set(_provenance(gene)["label_map"].values()) == {1, 0, None}
    table = _provenance("SBDS")["class_name_map"]
    assert {k: v["class"] for k, v in table.items()} == N.CLASS_NAME_MAP
    assert {v["label"] for v in table.values()} == {1, 0, None}
    assert set(N.CLASS_NAME_MAP.values()) == {"abnormal", "normal", "not_specified"}


@pytest.mark.parametrize("gene", ALL_GENES)
def test_the_classification_is_mechanical(gene):
    """Each score's class recomputed from score_ranges.json alone; for a deposit
    without ranges (SBDS), each class copied from the raw score table's own column,
    a value outside the class-name table being unassigned."""
    v = _variants(gene)
    if _has_ranges(gene):
        classes = _classes(gene)
        for r in v.itertuples():
            hit = [c["label"] for c in classes if N._matches(r.score, c)] if r.score == r.score else []
            assert r.range_class == (hit[0] if len(hit) == 1 else N.UNASSIGNED), r.hgvs_nt
        return
    col = N.CLASS_COLUMNS[0]
    raw = {r["accession"]: r[col] for r in _raw_table(gene)}
    for r in v.itertuples():
        value = raw[r.variant_urn]
        assert r.range_class == (value if value in N.CLASS_NAME_MAP else N.UNASSIGNED), r.hgvs_nt
    within = v[v["offset"].abs() <= N.MAX_DEPTH]
    assert not within["range_class_note"].str.startswith("class value not in").any()
    prov = _provenance(gene)
    assert prov["source"] == "MaveDB" and prov["score_ranges"].startswith("none")
    assert prov["label_source"].startswith(f"the score-table column {col}, not ranges")
    assert "IGVF cross-checked only" in prov["label_source"]
    assert prov["unmapped_values"] == [x for x in prov["column_values"] if x not in N.CLASS_NAME_MAP]


def test_the_shared_point_of_ctcf_is_found_and_none_elsewhere():
    """CTCF's Indeterminate and Functionally Normal ranges both include one point; a
    score there is unassigned (the author's rule of 2026-10-20)."""
    for gene in GENES:
        pts = N.shared_points(_classes(gene))
        assert all("point" in p for p in pts), pts                   # no real overlap
        assert (len(pts) == 1) == (gene == "CTCF"), (gene, pts)
        N.require_mechanical(gene, _classes(gene))                   # passes as it is


def test_ranges_that_overlap_beyond_a_point_stop_the_stage():
    """The author's rule covers a shared end point only; E16.10 stops on the rest."""
    def cls(label, lo, hi, lo_in=False, hi_in=True):
        return {"label": label, "functionalClassification": "abnormal", "range": [lo, hi],
                "inclusiveLowerBound": lo_in, "inclusiveUpperBound": hi_in}
    overlapping = [cls("A", None, -0.5), cls("B", -0.6, -0.2), cls("C", -0.2, None)]
    open_inside = [cls("A", None, -0.5), cls("B", -0.5, None), cls("C", -0.2, None)]
    for classes in (overlapping, open_inside):
        with pytest.raises(SystemExit, match="E16.10"):
            N.require_mechanical("X", classes)
    N.require_mechanical("X", [cls("A", None, -0.5), cls("B", -0.5, None)])   # abutting


def test_no_range_bound_is_written_into_the_code():
    """The ranges are read from the metadata, never typed in (E16.0): neither as the
    text of a bound nor as a number equal to one, in any notation."""
    import ast
    files = [PHASE1 / "src/evid_newgenes.py", *sorted((REPO / "scripts/stage_E16").glob("*.py"))]
    code = "\n".join(p.read_text() for p in files)
    numbers = {abs(n.value) for p in files for n in ast.walk(ast.parse(p.read_text()))
               if isinstance(n, ast.Constant) and type(n.value) in (int, float)}
    for gene in GENES:
        for c in _classes(gene):
            for bound in c["range"]:
                if bound is not None:
                    assert repr(bound) not in code and abs(bound) not in numbers, (gene, bound)


def test_the_author_record_check_stops_on_any_difference(monkeypatch):
    """E16.1, E16b.1: licence, publication status and dates checked against the
    author's record (CTCF and SFPQ of 2026-10-19, SBDS of 2026-10-20); any difference
    stops the stage."""
    import copy
    monkeypatch.chdir(PHASE1)                       # the stage's paths are from phase1/
    assert set(N.AUTHOR_RECORD) == set(N.AUTHOR_RECORD_DATE) == {"CTCF", "SFPQ", "SBDS"}
    for gene in N.AUTHOR_RECORD:
        meta = json.loads((EXT / f"{gene.lower()}_raw/scoreset_meta.json").read_text())
        exp = json.loads((EXT / f"{gene.lower()}_raw/experiment_meta.json").read_text())
        man = _manifest(gene)
        assert N.author_check(gene, meta, exp, man).startswith("agrees")
        paper = {"dbName": "bioRxiv", "identifier": "2026.01.01.000001"}
        for change in (lambda m, e: m.update(private=True),
                       lambda m, e: m["primaryPublicationIdentifiers"].append(paper),
                       lambda m, e: e["secondaryPublicationIdentifiers"].append(paper),
                       lambda m, e: m["license"].update(shortName="CC BY"),
                       lambda m, e: m.update(urn=m["urn"][:-1] + "9")):
            m2, e2 = copy.deepcopy(meta), copy.deepcopy(exp)
            change(m2, e2)
            with pytest.raises(SystemExit, match="E16.10"):
                N.author_check(gene, m2, e2, man)


# --------------------------------------------------------------------------
# E16b: where the SBDS classes may come from, fixed before any predictor result
# --------------------------------------------------------------------------
def _range(label, fc, lo, hi):
    return {"label": label, "functionalClassification": fc, "range": [lo, hi],
            "inclusiveLowerBound": True, "inclusiveUpperBound": False}


RANGES = {"investigator_provided": [{"title": "made up", "functionalClassifications": [
    _range("Functionally Abnormal", "abnormal", None, -0.5),
    _range("Functionally Normal", "normal", -0.5, None)]}]}


def test_the_class_source_follows_the_order_fixed_in_advance():
    """E16b.0: (a) one set of investigator-provided ranges, checked variant by variant
    against a class column; any disagreement stops; (b) else the class column;
    (c) else stop. More than one set of ranges stops too (E16.4)."""
    col = N.CLASS_COLUMNS[0]
    table = [{"accession": "x#1", "score": "-1.0", col: "functionally_abnormal"},
             {"accession": "x#2", "score": "0.2", col: "functionally_normal"}]
    assert N.label_source("X", table, RANGES)["kind"] == "ranges"            # (a), agreeing
    disagree = [dict(table[0]), dict(table[1], **{col: "functionally_abnormal"})]
    with pytest.raises(SystemExit, match="E16b.11"):                       # (a), one differs
        N.label_source("X", disagree, RANGES)
    src = N.label_source("X", table, {"investigator_provided": []})          # (b)
    assert src["kind"] == "column" and src["column"] == col and src["calibration"] is None
    bare = [{"accession": r["accession"], "score": r["score"]} for r in table]
    assert N.label_source("X", bare, RANGES)["kind"] == "ranges"             # (a), no column
    with pytest.raises(SystemExit, match="E16b.11"):                       # (c)
        N.label_source("X", bare, {"investigator_provided": []})
    two = {"investigator_provided": RANGES["investigator_provided"] * 2}
    with pytest.raises(SystemExit, match="E16.4"):
        N.label_source("X", table, two)


def test_a_class_outside_the_class_name_table_is_unassigned_and_stops_the_stage():
    """Only the table's three words map; any other value is unassigned and listed,
    and on an intron-side SNV within 50 nt it stops the stage (E16b.11)."""
    for value, cls in N.CLASS_NAME_MAP.items():
        assert N.column_class(value) == (value, "") and cls in N.CANON
    assert N.column_class("")[0] == N.column_class(None)[0] == N.UNASSIGNED
    cls, why = N.column_class("likely_abnormal")
    assert cls == N.UNASSIGNED and why.endswith("likely_abnormal")
    df = pd.DataFrame({"range_class": [N.UNASSIGNED], "offset": [-12], "range_class_note": [why]})
    with pytest.raises(SystemExit, match="E16b.11"):
        N.require_known_classes("X", df)
    N.require_known_classes("X", df.assign(offset=[-55]))                  # beyond 50 nt
    N.require_known_classes("X", df.assign(range_class_note=["no class in the column"]))


def test_the_prespecified_label_source_is_on_the_record_page():
    """The rule fixed by the author on 2026-10-20 is printed as the code holds it,
    before the SBDS record, and names MaveDB as the only source of labels."""
    text = (REPO / "docs/sbds-record-check.md").read_text()
    block = "\n".join(N.PRESPECIFIED)
    assert block in text and text.index(block) < text.index("\n## SBDS\n")
    assert "SBDS labels come from MaveDB only" in block and "never to build a label" in block
    for value, cls in N.CLASS_NAME_MAP.items():
        assert value in block


def test_the_labels_never_read_the_igvf_copy():
    """IGVF is a cross-check only (E16b.0): nothing that builds a label opens it."""
    import inspect
    for fn in (N.label_source, N.build_gene, N.classify, N.column_class, N.load_raw):
        src = inspect.getsource(fn)
        assert "IGVF" not in src and "igvf" not in src, fn.__name__


# --------------------------------------------------------------------------
# scope, offsets and coordinates
# --------------------------------------------------------------------------
@pytest.mark.parametrize("gene", ALL_GENES)
def test_the_table_holds_every_intron_side_snv_of_the_deposit_as_served(gene):
    """Read from the raw score table itself, with a pattern of this test's own: every
    intron-side SNV is one row, with the deposit's hgvs_nt, its score and the signed
    offset its c. description gives; nothing the stage counted is taken on trust."""
    pattern = re.compile(r":c\.[-*]?\d+([+-])(\d+)[ACGT]>[ACGT]$")
    raw = {r["accession"]: r for r in _raw_table(gene) if pattern.search(r["hgvs_nt"] or "")}
    prov = _provenance(gene)
    v = _variants(gene).set_index("variant_urn")
    assert prov["n_intronic_snv"] == len(raw) and prov["n_duplicates_dropped"] == 0
    assert set(v.index) == set(raw)
    for urn, r in raw.items():
        row, m = v.loc[urn], pattern.search(r["hgvs_nt"])
        assert row["hgvs_nt"] == r["hgvs_nt"], urn
        if r["score"] in ("", "NA"):
            assert row["score"] != row["score"], urn                   # NaN
        else:
            assert row["score"] == float(r["score"]), urn
        assert row["offset"] == (1 if m.group(1) == "+" else -1) * int(m.group(2)), urn


@pytest.mark.parametrize("gene", ALL_GENES)
def test_scope_is_flagged_and_no_row_is_dropped(gene):
    v, prov = _variants(gene), _provenance(gene)
    depth = v["offset"].abs()
    assert (v["in_scope"] == ((depth >= 3) & (depth <= 50))).all()
    assert (depth <= 2).any() and (v.loc[depth <= 2, "stratum"] == "pm12").all()
    assert (v.loc[depth > 50, "stratum"] == N.DEEP).all()
    assert ((v["side"] == "donor") == (v["offset"] > 0)).all()
    assert len(v) == prov["n_rows_written"] == prov["n_intronic_snv"] - prov["n_duplicates_dropped"]
    assert prov["deepest_offset"] == depth.max()


@pytest.mark.parametrize("gene", ALL_GENES)
def test_no_variant_appears_twice(gene):
    v = _variants(gene)
    assert not v.duplicated(["chrom", "pos", "ref", "alt"]).any()
    assert v["variant_urn"].is_unique


@pytest.mark.parametrize("gene", ALL_GENES)
def test_every_coordinate_agrees_with_mavedb_and_the_exon_structure(gene):
    """pos = MaveDB's own GRCh38 position of the exonic base the c. position names,
    plus the strand-signed offset; and |offset| = the distance to that exon."""
    v = _variants(gene)
    assert v["anchor_agrees"].notna().all() and v["anchor_agrees"].all()
    assert v["exon_distance"].notna().all() and (v["exon_distance"] == v["offset"].abs()).all()
    # E16.3 compares with the NEAREST exon boundary: the offset's own side is the nearer
    assert v["exon_distance_other_side"].notna().all()
    assert (v["exon_distance"] <= v["exon_distance_other_side"]).all()
    row = {r["gene"]: r for r in _csv(f"{_prefix(gene)}_coordinate_check.csv")}[gene]
    n = len(v)
    assert int(row["n_mapped"]) == int(row["n_anchor_agree"]) == n
    assert int(row["n_exon_distance_equals_offset"]) == int(row["n_offset_side_is_nearest"]) == n


@pytest.mark.parametrize("gene", ALL_GENES)
def test_every_position_is_mavedbs_own_exonic_anchor_plus_the_offset(gene):
    """Recomputed here from variant_details.ndjson, not read from the stage's columns:
    MaveDB's GRCh38 position of the exonic base c.N (current records, the MANE
    chromosome accession), plus the offset times the strand, is the table's position."""
    tm = _manifest(gene)["transcript_model"]["mane_select"]
    acc, strand = tm["grch38_chr_accession"], (1 if tm["strand"] == "+" else -1)
    exonic = re.compile(r"^(?:N[MR]_|ENST)[\d.]+:(c\.[-*]?\d+)[ACGT]>[ACGT]$")
    genomic = re.compile(rf"^{re.escape(acc)}:g\.(\d+)[ACGT]>[ACGT]$")
    anchor = {}
    with open(EXT / f"{gene.lower()}_raw/variant_details.ndjson") as f:
        for line in f:
            r = json.loads(line) if line.strip() else {}
            m = exonic.match(r.get("targetHgvs") or "")
            if r.get("isCurrent") is not True or not m or not r.get("postMapped"):
                continue
            for e in r["postMapped"].get("expressions", []):
                g = genomic.match(e.get("value", "")) if e.get("syntax") == "hgvs.g" else None
                if g:
                    anchor.setdefault(m.group(1), set()).add(int(g.group(1)))
    v = _variants(gene)
    for r in v.itertuples():
        c = re.match(r"^c\.([-*]?\d+)[+-]", r.hgvs_nt.split(":", 1)[1]).group(1)
        a = anchor.get(f"c.{c}", set())
        assert len(a) == 1 and next(iter(a)) + strand * r.offset == r.pos, r.hgvs_nt


@pytest.mark.parametrize("run", list(N.RUNS), ids=lambda k: k or "E16")
def test_the_offset_sample_draws_ten_per_gene_and_every_one_agrees(run):
    pre = N.RUNS[run]["prefix"]
    rows = _csv(f"{pre}_offset_sample.csv")
    strand = {r["gene"]: r["strand"] for r in _csv(f"{pre}_coordinate_check.csv")}
    assert "-" in strand.values()                       # a minus-strand gene is sampled
    for gene in N.RUNS[run]["genes"]:
        mine = [r for r in rows if r["gene"] == gene]
        assert len(mine) == N.N_SAMPLE
        assert all(int(r["exon_distance"]) == abs(int(r["offset"]))
                   <= int(r["exon_distance_other_side"]) for r in mine)
        assert all(r["anchor_agrees"] == "True" for r in mine)


def test_the_authors_own_coordinates_agree_with_every_sbds_variant():
    """SBDS's score table also carries the authors' chrom/pos/ref/alt (E16b.3); every
    intron-side SNV agrees with the mapped coordinate, recomputed here from the raw
    table. MaveDB gives no postMapped allele for any of them, and the record says so."""
    raw = {r["hgvs_nt"]: r for r in _raw_table("SBDS")}
    v = _variants("SBDS")
    for r in v.itertuples():
        a = raw[r.hgvs_nt]
        assert ((a["chrom"].removeprefix("chr"), int(a["pos"]), a["ref"], a["alt"])
                == (str(r.chrom), int(r.pos), r.ref, r.alt)), r.hgvs_nt
    extra = _provenance("SBDS")["coordinates"]["extra_checks"]
    assert extra["author_columns"]["n_compared"] == extra["author_columns"]["n_agree"] == len(v)
    pm = extra["intronic_postmapped"]
    assert pm["n_compared"] == pm["n_agree"] and not pm["differ"]


# --------------------------------------------------------------------------
# counts
# --------------------------------------------------------------------------
def _counts(prefix: str) -> dict:
    out = {}
    for r in _csv(f"{prefix}_counts_by_offset.csv"):
        out[(r["gene"], r["scheme"], r["bin"], r["range_class"])] = int(r["n"])
    return out


def test_each_stratum_count_is_the_sum_of_its_per_nt_counts():
    spans = {"pm12": range(1, 3), "s3_10": range(3, 11), "s11_50": range(11, 51),
             "s3_50": range(3, 51)}
    for gene in ALL_GENES:
        n = _counts(_prefix(gene))
        for cls in list(_class_of(gene)) + [N.UNASSIGNED]:
            for st, span in spans.items():
                assert n[(gene, "stratum", st, cls)] == sum(
                    n[(gene, "per_nt", str(d), cls)] for d in span), (gene, st, cls)
            assert n[(gene, "stratum", N.DEEP, cls)] == n[(gene, "per_nt", ">50", cls)]


def test_the_stratum_names_are_evid_commons():
    for run in N.RUNS.values():
        used = {r["bin"] for r in _csv(f"{run['prefix']}_counts_by_offset.csv")
                if r["scheme"] == "stratum"}
        assert used - {N.DEEP} <= set(K.STRATA)


@pytest.mark.parametrize("run", list(N.RUNS), ids=lambda k: k or "E16")
def test_the_label_summary_is_the_tables_counted(run):
    spec = N.RUNS[run]
    rows = {(r["gene"], r["stratum"]): r for r in _csv(f"{spec['prefix']}_label_summary.csv")}
    frames = {g: _variants(g) for g in spec["genes"]}
    if spec["total"]:
        frames[spec["total"]] = pd.concat(frames.values(), ignore_index=True)
    assert {g for g, _ in rows} == set(frames)          # SBDS: its own rows, no total row
    for gene, v in frames.items():
        depth = v["offset"].abs()
        for st, keep in {"pm12": depth <= 2, "s3_10": (depth >= 3) & (depth <= 10),
                         "s11_50": (depth >= 11) & (depth <= 50),
                         "s3_50": (depth >= 3) & (depth <= 50), N.DEEP: depth > 50}.items():
            sub, r = v[keep], rows[(gene, st)]
            dmg, nor = int((sub["label"] == 1).sum()), int((sub["label"] == 0).sum())
            una = int((sub["range_class"] == N.UNASSIGNED).sum())
            assert int(r["n_total"]) == len(sub)
            assert (int(r["n_damaging"]), int(r["n_normal"]), int(r["n_unassigned"])) == (dmg, nor, una)
            assert int(r["n_uncertain"]) == len(sub) - dmg - nor - una
            assert (r["evaluable_min10"] == "True") == (dmg >= K.MIN_POS and nor >= K.MIN_NEG)


def test_the_sbds_files_have_e16s_columns_and_strata():
    """E16b writes the same six tables as E16, column for column and stratum for
    stratum, for SBDS's rows only."""
    for kind in ("record_check", "counts_by_offset", "label_summary", "coordinate_check",
                 "offset_sample", "controls"):
        with open(REPORTS / f"newgenes_{kind}.csv", newline="") as a, \
                open(REPORTS / f"sbds_{kind}.csv", newline="") as b:
            assert next(csv.reader(a)) == next(csv.reader(b)), kind
        assert {r["gene"] for r in _csv(f"sbds_{kind}.csv")} == {"SBDS"}, kind
    for kind, col in (("label_summary", "stratum"), ("counts_by_offset", "bin")):
        e16 = [r[col] for r in _csv(f"newgenes_{kind}.csv") if r["gene"] == GENES[0]]
        sbds = [r[col] for r in _csv(f"sbds_{kind}.csv")]
        assert sorted(set(e16)) == sorted(set(sbds)), kind


# --------------------------------------------------------------------------
# E16b: the IGVF cross-check and the label versions, checks that change nothing
# --------------------------------------------------------------------------
def _key(hgvs: str) -> str:
    """This test's own alignment key: the c. description, without the transcript."""
    return hgvs.split(":")[-1].replace(" ", "").lower()


def test_the_igvf_crosscheck_is_recounted_from_the_two_files():
    """Each category recounted here from the MaveDB score table and the IGVF file as
    served; every variant whose class differs is listed with both sides' values."""
    col = N.CLASS_COLUMNS[0]
    man = json.loads((EXT / "igvf_sbds/manifest.json").read_text())
    gz = REPO / next(e["path"] for e in man["files"] if e["path"].endswith(".tsv.gz"))
    for e in man["files"]:
        assert _sha(REPO / e["path"]) == e["sha256"], e["path"]
    with gzip.open(gz, "rt", newline="") as f:
        igvf = {_key(r["hgvs_c"]): r for r in csv.DictReader(f, delimiter="\t")}
    mave = {_key(r["hgvs_nt"]): r for r in _raw_table("SBDS")}
    v = _variants("SBDS")
    intronic = {_key(h) for h, o in zip(v["hgvs_nt"], v["offset"]) if abs(o) <= N.MAX_DEPTH}
    rows = _csv("sbds_igvf_crosscheck.csv")
    for scope, keys in (("all rows", set(mave) | set(igvf)),
                        ("intron-side SNVs 1-50 nt", intronic)):
        want = Counter()
        for k in keys:
            a, b = mave.get(k), igvf.get(k)
            if a and b:
                want["same class" if a[col] == b[col] else "different class"] += 1
            else:
                want["MaveDB only" if a else "IGVF only"] += 1
        got = {r["category"]: int(r["n"]) for r in rows
               if r["row_type"] == "count" and r["scope"] == scope}
        assert got == {c: want[c] for c in ("same class", "different class", "MaveDB only",
                                            "IGVF only")}, scope
        if scope == "all rows":
            listed = [r for r in rows if r["row_type"] == "variant"]
            assert len(listed) == want["different class"]
            for r in listed:
                a, b = mave[_key(r["hgvs_nt"])], igvf[_key(r["hgvs_nt"])]
                assert (r["mavedb_class"], r["mavedb_score"]) == (a[col], a["score"]) != \
                    (r["igvf_class"], r["igvf_score"]) == (b[col], b["score"])


def test_the_listed_igvf_differences_include_every_intron_side_one():
    rows = _csv("sbds_igvf_crosscheck.csv")
    n = {r["category"]: int(r["n"]) for r in rows
         if r["row_type"] == "count" and r["scope"] == "intron-side SNVs 1-50 nt"}
    listed = [r for r in rows if r["row_type"] == "variant" and r["offset"] != ""]
    assert len(listed) == n["different class"]
    assert all(abs(int(r["offset"])) <= N.MAX_DEPTH for r in listed)


def test_the_label_version_check_adds_up_and_lists_every_disagreement():
    """Each archived version of BARD1 and PALB2 against the analysis set's label files,
    over the set's intron-side SNVs at 1-50 nt; every disagreement listed once."""
    check, diff = _csv("label_version_check.csv"), _csv("label_version_diff.csv")
    aset = pd.read_parquet(PHASE1 / "data/evidence/analysis_set_v1.parquet")
    assert {r["gene"] for r in check} == set(N.LABEL_GENES)
    for gene in N.LABEL_GENES:
        folder = EXT / "label_versions" / f"{gene.lower()}_raw"
        man = json.loads((folder / "manifest.json").read_text())
        for e in man["files"]:
            assert _sha(REPO / e["path"]) == e["sha256"], e["path"]
        rows = [r for r in check if r["gene"] == gene]
        assert [r["version"] for r in rows] == [v["urn"] for v in man["versions"]]
        n_set = int(((aset["gene"] == gene) & aset["intron_offset_abs"].between(1, N.MAX_DEPTH)).sum())
        for r, ver in zip(rows, man["versions"]):
            n = {k: int(r[k]) for k in ("n_ours", "n_found", "n_compared", "n_agree",
                                        "n_disagree", "n_found_without_class", "n_only_ours")}
            assert r["current"] == str(ver["superseded_by"] is None)
            assert n["n_ours"] == n_set == n["n_found"] + n["n_only_ours"]
            assert n["n_found"] == n["n_compared"] + n["n_found_without_class"]
            assert n["n_compared"] == n["n_agree"] + n["n_disagree"]
            if r["class_source"] == "no classes":
                assert n["n_compared"] == 0
            mine = [d for d in diff if d["gene"] == gene and d["version"] == r["version"]]
            assert len(mine) == n["n_disagree"] == len({d["hgvs_nt"] for d in mine})
            assert all(d["our_class"] != d["archive_class"] for d in mine)
            if r["align_key"] == "c. HGVS":            # the archived score, as served
                with open(folder / f"urn_mavedb_{ver['urn'].split(':')[-1]}.csv", newline="") as f:
                    score = {_key(x["hgvs_nt"]): x["score"] for x in csv.DictReader(f)}
                assert all(score[_key(d["hgvs_nt"])] == d["archive_score"] for d in mine)


def test_the_label_version_check_changes_no_label():
    """A check only (E16b.6): the two label files are the bytes the analysis set was
    built from, and config.ASSAY_LABELS still names them."""
    srcs = json.loads((PHASE1 / "data/evidence/analysis_set_v1.manifest.json").read_text())
    recorded = {s["file"]: s["sha256"] for s in srcs["sources"]["assay_labels"].values()}
    for gene in N.LABEL_GENES:
        name = C.ASSAY_LABELS[gene]["file"]
        assert _sha(PHASE1 / "data/assay_labels" / name) == recorded[name], name


# --------------------------------------------------------------------------
# E16b.0 rule 3: no predictor score for any SBDS variant anywhere in the repository
# --------------------------------------------------------------------------
PREDICTOR_COLUMNS = sorted(set(K.PANEL + K.OPTIONAL + [K.FUSION] + C.PREDICTOR_FEATURES
                               + C.EVO_FEATURES))


def _repository_files(*folders: str) -> list[Path]:
    """What the next commit would carry under these folders: tracked files and
    untracked ones not ignored, in a clone; the tree, in a release archive."""
    if (REPO / ".git").exists():
        out = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z", "--cached", "--others",
                              "--exclude-standard", "--", *folders],
                             capture_output=True, text=True, check=True).stdout
        return sorted(REPO / p for p in out.split("\0") if p and (REPO / p).is_file())
    sys.path.insert(0, str(REPO / "tests"))
    from _checkout import checkout_files
    return [p for f in folders for p in checkout_files(f)]


def _columns(path: Path) -> list[str]:
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq
        return pq.read_schema(path).names
    opener = gzip.open if path.suffix == ".gz" else open
    stem = path.name.removesuffix(".gz")
    if not stem.endswith((".csv", ".tsv", ".txt")):
        return []
    with opener(path, "rt", newline="") as f:
        line = f.readline().rstrip("\r\n")
    return line.split("," if stem.endswith(".csv") else "\t")


# the files E16b wrote for SBDS, by name (evid_newgenes --gene SBDS); with the files
# its manifests list, these are E16b's own
E16B_OUTPUTS = ["phase1/data/evidence/external/sbds_raw/manifest.json",
                "phase1/data/evidence/external/sbds_variants.parquet",
                "phase1/data/evidence/external/sbds_variants.parquet.provenance.json",
                "phase1/data/evidence/external/sbds_reference_check.json",
                "phase1/data/evidence/external/igvf_sbds/manifest.json"]
# the tagger time of frozen-thresholds-v1, from its tag object
TAG_UTC = "2026-10-22T16:34:25+00:00"


def _scored_after(rel: str) -> str | None:
    """The first scoring time an E18 score file records: the file's own record if it
    is one (a provenance or manifest json with scored_utc_range), else the provenance
    beside it (<file>.provenance.json; for <x>.environment.json, <x>.provenance.json).
    None when there is no such record."""
    path = REPO / rel
    candidates = []
    if path.suffix == ".json" and not rel.endswith(".environment.json"):
        candidates.append(path)
    stem = rel.removesuffix(".environment.json") if rel.endswith(".environment.json") else rel
    candidates.append(REPO / f"{stem}.provenance.json")
    for c in candidates:
        if c.exists():
            rng = json.loads(c.read_text()).get("scored_utc_range")
            if rng:
                return rng[0]
    return None


def test_no_predictor_score_for_any_sbds_variant_is_in_the_repository():
    """The git history is the evidence that the SBDS label source was fixed before any
    predictor result (E16b.0 rule 3). Since the thresholds were frozen (E18), SBDS is
    scored, so the rule now reads: under phase1/data, a file that names SBDS, its
    transcripts, or a position of its region is either one of E16b's own files, which
    still carry no predictor column, or an E18 score file (or its record) whose scoring
    began after the tag frozen-thresholds-v1; and the SBDS labels were built before the
    tag and before any SBDS score. E16b's own files are those its manifests list and the
    files it wrote (E16B_OUTPUTS), not every file named sbds_*: E18's sbds_<column>
    files are score files. phase1/reports is no longer scanned: from E18's external test
    on, its outputs carry SBDS results derived from E18's score files."""
    from datetime import datetime
    man = _manifest("SBDS")
    tm = man["transcript_model"]
    ids = ["sbds", tm["mane_select"]["refseq_nuc"].split(".")[0].lower(),
           tm["mane_select"]["ensembl_nuc"].split(".")[0].lower()]
    chrom, lo, hi = str(tm["region"]["chrom"]), tm["region"]["start"], tm["region"]["end"]
    in_text = re.compile(rf"(?<![0-9A-Za-z])(?:chr)?{chrom}[-:_,\t ](\d{{7,9}})(?!\d)".encode())
    igvf = json.loads((EXT / "igvf_sbds/manifest.json").read_text())
    own_files = ({e["path"] for e in man["files"]}
                 | {item["path"] for item in tm["cache_files"].values()}
                 | {str((CACHE / "MANIFEST.json").relative_to(REPO))}
                 | {e["path"] for e in igvf["files"]} | set(E16B_OUTPUTS))

    def own(rel: str) -> bool:
        return rel in own_files

    def in_region(df: pd.DataFrame) -> int:
        cols = {c.lower(): c for c in df.columns}
        c = next((cols[k] for k in ("chrom", "chr", "chromosome") if k in cols), None)
        p = next((cols[k] for k in ("pos", "position") if k in cols), None)
        if c is None or p is None:
            return 0
        return int(((df[c].astype(str).str.removeprefix("chr") == chrom)
                    & pd.to_numeric(df[p], errors="coerce").between(lo, hi)).sum())

    mentions, scored = {}, []
    for path in _repository_files("phase1/data"):
        rel = path.relative_to(REPO).as_posix()
        why = []
        if path.suffix == ".parquet":
            df = pd.read_parquet(path)
            for c in df.columns:
                if pd.api.types.is_numeric_dtype(df[c]) or pd.api.types.is_bool_dtype(df[c]):
                    continue
                s = df[c].dropna().astype(str).str.lower()
                why += [f"{i} in {c}" for i in ids if s.str.contains(i, regex=False).any()]
                ident = s.str.extract(r"^(?:chr)?(\w+)[-:_](\d+)[-:_]")
                near = (ident[0] == chrom) & pd.to_numeric(ident[1], errors="coerce").between(lo, hi)
                why += [f"{int(near.sum())} ids of the region in {c}"] if near.any() else []
            n = in_region(df)
            why += [f"{n} rows in the region"] if n else []
        else:
            cols = _columns(path)
            if {"chrom", "pos"} <= {c.lower() for c in cols}:
                df = pd.read_csv(path, sep="," if ".csv" in path.name else "\t", dtype=str,
                                 usecols=[c for c in cols if c.lower() in ("chrom", "pos")])
                n = in_region(df)
                why += [f"{n} rows in the region"] if n else []
            data = path.read_bytes()
            if path.suffix == ".gz":
                try:
                    data = gzip.decompress(data)
                except OSError:
                    pass
            low = data.lower()
            why += [i for i in ids if i.encode() in low]
            near = [int(m.group(1)) for m in in_text.finditer(data)]
            why += ["positions of the region"] if any(lo <= p <= hi for p in near) else []
        if why:
            mentions[rel] = why
        if own(rel):
            scored += [f"{rel}: {c}" for c in _columns(path)
                       if any(c.lower() == p or c.lower().startswith(p + "_")
                              for p in PREDICTOR_COLUMNS)]
    built = datetime.fromisoformat(json.loads(
        (EXT / "sbds_variants.parquet.provenance.json").read_text())["built_utc"])
    tagged = datetime.fromisoformat(TAG_UTC)
    assert built < tagged
    after = {rel: _scored_after(rel) for rel in mentions if not own(rel)}
    assert mentions, "the scan found not even E16b's own files"
    assert all(t is not None for t in after.values()), [r for r, t in after.items() if t is None]
    assert all(datetime.fromisoformat(t) > tagged for t in after.values()), after
    assert all(datetime.fromisoformat(t) > built for t in after.values()), after
    assert not scored, scored


# --------------------------------------------------------------------------
# records and manifests
# --------------------------------------------------------------------------
@pytest.mark.parametrize("gene", ALL_GENES)
def test_every_file_a_manifest_lists_is_the_file_on_disk(gene):
    man = _manifest(gene)
    registered = json.loads((CACHE / "MANIFEST.json").read_text())["files"]
    for e in man["files"]:
        path = REPO / e["path"]
        assert path.stat().st_size == e["bytes"] and _sha(path) == e["sha256"], e["path"]
        assert e["bytes"] <= 25_000_000, f"{e['path']} is over the 25 MB limit (E16.2)"
    for item in man["transcript_model"]["cache_files"].values():
        path = REPO / item["path"]
        assert _sha(path) == item["sha256"] == registered[str(path.relative_to(CACHE))]["sha256"]
    prov = _provenance(gene)
    assert prov["inputs"] == {e["path"]: e["sha256"] for e in man["files"]}
    assert prov["sha256"] == _sha(EXT / f"{gene.lower()}_variants.parquet")


def _cross_check_manifests() -> list[Path]:
    return [EXT / "igvf_sbds/manifest.json",
            *(EXT / "label_versions" / f"{g.lower()}_raw/manifest.json" for g in N.LABEL_GENES)]


def test_the_raw_deposits_are_stored_byte_for_byte():
    """MaveDB serves its score tables with CRLF line ends. The repository turns every
    other *.csv to LF on checkout, which would change the bytes each manifest pins."""
    if not (REPO / ".git").exists():
        pytest.skip("a release archive has no .git")
    paths = [e["path"] for g in ALL_GENES for e in _manifest(g)["files"]]
    paths += [e["path"] for m in _cross_check_manifests()
              for e in json.loads(m.read_text())["files"]]
    out = subprocess.run(["git", "check-attr", "text", "--", *paths], cwd=REPO,
                         capture_output=True, text=True, check=True).stdout.splitlines()
    assert len(out) == len(paths) and all(l.endswith(": text: unset") for l in out), out[:2]


@pytest.mark.parametrize("gene", ALL_GENES)
def test_the_record_row_is_read_from_the_metadata(gene):
    meta = json.loads((EXT / f"{gene.lower()}_raw/scoreset_meta.json").read_text())
    row = {r["gene"]: r for r in _csv(f"{_prefix(gene)}_record_check.csv")}[gene]
    assert row["urn"] == meta["urn"] == _manifest(gene)["urn"]
    assert row["license"].startswith(f"{meta['license']['shortName']} {meta['license']['version']}")
    assert meta["license"]["shortName"] == "CC0"
    assert int(row["n_variants_total"]) == meta["numVariants"]
    assert row["target_accession"] == meta["targetGenes"][0]["targetAccession"]["accession"]
    pubs = meta["primaryPublicationIdentifiers"] + meta["secondaryPublicationIdentifiers"]
    assert (row["publication_identifiers"] == "not_provided") == (not pubs)
    assert meta.get("supersedingScoreSet") is None              # the current version
    if gene in N.AUTHOR_RECORD:
        assert row["author_record_check"].startswith("agrees")
        assert row["urn"] == N.AUTHOR_RECORD[gene]["urn"]       # SBDS: -a-1, as decided


@pytest.mark.parametrize("run", list(N.RUNS), ids=lambda k: k or "E16")
def test_every_fact_on_the_record_page_names_its_file(run):
    """A fetched file and its time, or the search record a URN was chosen from, on
    every line of each gene's section (the E16b page also prints the rule fixed in
    advance, a decision rather than a fact read from a file)."""
    spec = N.RUNS[run]
    lines, section = [], None
    for line in (PHASE1 / spec["doc"]).resolve().read_text().splitlines():
        if line.startswith("## "):
            section = line[3:]
        elif line.startswith("- ") and section in spec["genes"]:
            lines.append(line)
    bare = [l for l in lines if "fetched 20" not in l and "scripts/stage_E16/run_log/" not in l]
    assert lines and not bare, bare[:2]


# --------------------------------------------------------------------------
# the entry point knows the stages and their inputs
# --------------------------------------------------------------------------
def _entry_point():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_entry_point_runs_the_stage_after_e15_and_before_the_delta():
    E = _entry_point()
    mods = [m for m, _, _ in E.STAGES]
    i = mods.index("src.evid_newgenes")
    assert mods[i - 1] == "src.evid_practical_metrics" and i < mods.index("src.evid_delta")
    assert [p.name for p in E.NEWGENE_RAW] == [f"{g.lower()}_raw" for g in ALL_GENES]
    assert REPO / "docs" / "newgenes-record-check.md" in E.OUTPUT_ROOTS
    assert PHASE1 / "data" / "evidence" / "refs" in E.NOT_OUTPUTS


def test_the_entry_point_runs_e16b_right_after_e16_and_before_the_delta():
    E = _entry_point()
    rows = [(m, a) for m, _, a in E.STAGES]
    i = rows.index(("src.evid_newgenes", []))
    assert rows[i + 1] == ("src.evid_newgenes", ["--gene", "SBDS"])
    assert [m for m, _ in rows].index("src.evid_delta") > i + 1    # E16c between
    assert REPO / "docs" / "sbds-record-check.md" in E.OUTPUT_ROOTS


def test_every_input_the_stage_reads_is_listed_as_one():
    E = _entry_point()
    inputs = E._run_inputs()
    for gene in ALL_GENES:
        man = _manifest(gene)
        need = [e["path"] for e in man["files"]]
        need += [i["path"] for k, i in man["transcript_model"]["cache_files"].items()]
        need.append(f"phase1/data/evidence/external/{gene.lower()}_raw/manifest.json")
        assert set(need) <= inputs, sorted(set(need) - inputs)
    for m in _cross_check_manifests():
        need = [e["path"] for e in json.loads(m.read_text())["files"]]
        need.append(str(m.relative_to(REPO)))
        assert set(need) <= inputs, sorted(set(need) - inputs)
    assert "phase1/data/evidence/external/newgenes_reference_check.json" in inputs
    assert "phase1/data/evidence/external/sbds_reference_check.json" in inputs


E16_OUTPUTS = ([f"phase1/reports/evidence/newgenes_{k}.csv" for k in (
                   "record_check", "counts_by_offset", "label_summary", "coordinate_check",
                   "offset_sample", "controls")]
               + ["docs/newgenes-record-check.md"]
               + [f"phase1/data/evidence/external/{g.lower()}_variants.parquet{s}"
                  for g in GENES for s in ("", ".provenance.json")])


def test_without_the_argument_the_stage_writes_e16s_files_as_they_were():
    """E16b made the gene list a parameter (E16b.8): run without it, the stage's 15
    outputs are E16's, byte for byte, as committed at E16, the build time of a
    provenance record aside. reproduce_evidence.py --verify reruns the stage; this
    holds the result to E16's commit, where git has it (not in a shallow clone)."""
    assert len(E16_OUTPUTS) == 15
    if not (REPO / ".git").exists():
        pytest.skip("a release archive has no .git")
    if subprocess.run(["git", "-C", str(REPO), "cat-file", "-e", f"{E16_COMMIT}^{{commit}}"],
                      capture_output=True).returncode:
        pytest.skip(f"E16's commit {E16_COMMIT[:7]} is not in this clone (shallow)")
    for rel in E16_OUTPUTS:
        then = subprocess.run(["git", "-C", str(REPO), "show", f"{E16_COMMIT}:{rel}"],
                              capture_output=True, check=True).stdout
        now = (REPO / rel).read_bytes()
        if rel.endswith(".provenance.json"):
            a, b = json.loads(then), json.loads(now)
            a.pop("built_utc"), b.pop("built_utc")
            assert a == b, rel
        else:
            assert then == now, rel
