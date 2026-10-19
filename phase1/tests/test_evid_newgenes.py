"""The four further external genes (stage E16): what their tables may claim.

Every number the stage writes traces to a file MaveDB served or a reference the
repository pins: these tests hold the variant tables to the reference genome, the
labels to the deposits' own ranges and the seven genes' rule for intermediate
classes, the counts to one another, and every manifest to the files on disk.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import re
import subprocess
import sys
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

GENES = N.GENES


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


def _csv(name: str) -> list[dict]:
    with open(REPORTS / name, newline="") as f:
        return list(csv.DictReader(f))


# --------------------------------------------------------------------------
# reference bases: 100%, three ways
# --------------------------------------------------------------------------
@pytest.mark.parametrize("gene", GENES)
def test_the_recorded_reference_check_covers_every_variant_of_the_table(gene):
    """scripts/stage_E16/check_reference.py, against the Ensembl release-112
    chromosome, on exactly the table committed now."""
    rec = json.loads((EXT / "newgenes_reference_check.json").read_text())
    g = {r["gene"]: r for r in rec["genes"]}[gene]
    table = EXT / f"{gene.lower()}_variants.parquet"
    ref = _manifest(gene)["reference_fasta"]
    assert g["variants_sha256"] == _sha(table)
    assert g["n"] == g["matches"] == len(_variants(gene)) and g["mismatches"] == 0
    assert g["fasta"]["sha256_fa"] == ref["sha256_fa"] and g["fasta"]["url"] == ref["url"]
    assert g["fasta"]["ensembl_release"] == 112


@pytest.mark.parametrize("gene", GENES)
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


@pytest.mark.parametrize("gene", GENES)
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
# labels: the deposit's own ranges, the seven genes' rule
# --------------------------------------------------------------------------
@pytest.mark.parametrize("gene", GENES)
def test_range_classes_are_only_the_deposits_own_and_unassigned(gene):
    allowed = {c["label"] for c in _classes(gene)} | {N.UNASSIGNED}
    assert set(_variants(gene)["range_class"]) <= allowed
    for name in ("newgenes_counts_by_offset.csv", "newgenes_controls.csv"):
        assert {r["range_class"] for r in _csv(name) if r["gene"] == gene} == allowed


@pytest.mark.parametrize("gene", GENES)
def test_each_range_class_has_one_label_and_it_follows_the_seven_genes_rule(gene):
    """abnormal -> 1, normal -> 0, any other class NA, like the intermediate classes
    of config.ASSAY_LABELS; unassigned is NA too."""
    fc = {c["label"]: c["functionalClassification"] for c in _classes(gene)}
    v = _variants(gene)
    for cls, labels in v.groupby("range_class")["label"]:
        values = set(labels.dropna())
        assert len(values) <= 1 and (values == set() or labels.notna().all()), cls
        want = {"abnormal": {1.0}, "normal": {0.0}}.get(fc.get(cls), set())
        assert values == want, (cls, values)


def test_the_label_rule_is_the_one_the_seven_genes_use():
    """Every map in ASSAY_LABELS sends its classes to damaging, normal or None (an
    intermediate class, excluded); E16's map has the same three outcomes and no other."""
    for gene, spec in C.ASSAY_LABELS.items():
        outcomes = set(spec["map"].values())
        assert outcomes <= {1.0, 0.0, None} and {1.0, 0.0} <= outcomes, gene
        assert None in outcomes, f"{gene}: no intermediate class mapped to None"
    assert N.LABEL_OF == {"abnormal": 1.0, "normal": 0.0}
    for gene in GENES:
        prov = json.loads((EXT / f"{gene.lower()}_variants.parquet.provenance.json").read_text())
        assert set(prov["label_map"].values()) == {1, 0, None}


@pytest.mark.parametrize("gene", GENES)
def test_the_classification_is_mechanical(gene):
    """Each score's class recomputed from score_ranges.json alone."""
    classes = _classes(gene)
    for r in _variants(gene).itertuples():
        hit = [c["label"] for c in classes if N._matches(r.score, c)] if r.score == r.score else []
        assert r.range_class == (hit[0] if len(hit) == 1 else N.UNASSIGNED), r.hgvs_nt


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
    """E16.1: licence, publication status and dates checked against the author's record
    of 2026-10-19; any difference stops the stage."""
    import copy
    monkeypatch.chdir(PHASE1)                       # the stage's paths are from phase1/
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
# scope, offsets and coordinates
# --------------------------------------------------------------------------
@pytest.mark.parametrize("gene", GENES)
def test_the_table_holds_every_intron_side_snv_of_the_deposit_as_served(gene):
    """Read from the raw score table itself, with a pattern of this test's own: every
    intron-side SNV is one row, with the deposit's hgvs_nt, its score and the signed
    offset its c. description gives; nothing the stage counted is taken on trust."""
    man = _manifest(gene)
    path = REPO / next(e["path"] for e in man["files"] if e["path"].endswith(".csv"))
    pattern = re.compile(r":c\.[-*]?\d+([+-])(\d+)[ACGT]>[ACGT]$")
    with open(path, newline="") as f:
        raw = {r["accession"]: r for r in csv.DictReader(f) if pattern.search(r["hgvs_nt"] or "")}
    prov = json.loads((EXT / f"{gene.lower()}_variants.parquet.provenance.json").read_text())
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


@pytest.mark.parametrize("gene", GENES)
def test_scope_is_flagged_and_no_row_is_dropped(gene):
    v, prov = _variants(gene), json.loads(
        (EXT / f"{gene.lower()}_variants.parquet.provenance.json").read_text())
    depth = v["offset"].abs()
    assert (v["in_scope"] == ((depth >= 3) & (depth <= 50))).all()
    assert (depth <= 2).any() and (v.loc[depth <= 2, "stratum"] == "pm12").all()
    assert (v.loc[depth > 50, "stratum"] == N.DEEP).all()
    assert ((v["side"] == "donor") == (v["offset"] > 0)).all()
    assert len(v) == prov["n_rows_written"] == prov["n_intronic_snv"] - prov["n_duplicates_dropped"]
    assert prov["deepest_offset"] == depth.max()


@pytest.mark.parametrize("gene", GENES)
def test_no_variant_appears_twice(gene):
    v = _variants(gene)
    assert not v.duplicated(["chrom", "pos", "ref", "alt"]).any()
    assert v["variant_urn"].is_unique


@pytest.mark.parametrize("gene", GENES)
def test_every_coordinate_agrees_with_mavedb_and_the_exon_structure(gene):
    """pos = MaveDB's own GRCh38 position of the exonic base the c. position names,
    plus the strand-signed offset; and |offset| = the distance to that exon."""
    v = _variants(gene)
    assert v["anchor_agrees"].notna().all() and v["anchor_agrees"].all()
    assert v["exon_distance"].notna().all() and (v["exon_distance"] == v["offset"].abs()).all()
    # E16.3 compares with the NEAREST exon boundary: the offset's own side is the nearer
    assert v["exon_distance_other_side"].notna().all()
    assert (v["exon_distance"] <= v["exon_distance_other_side"]).all()
    row = {r["gene"]: r for r in _csv("newgenes_coordinate_check.csv")}[gene]
    n = len(v)
    assert int(row["n_mapped"]) == int(row["n_anchor_agree"]) == n
    assert int(row["n_exon_distance_equals_offset"]) == int(row["n_offset_side_is_nearest"]) == n


@pytest.mark.parametrize("gene", GENES)
def test_every_position_is_mavedbs_own_exonic_anchor_plus_the_offset(gene):
    """Recomputed here from variant_details.ndjson, not read from the stage's columns:
    MaveDB's GRCh38 position of the exonic base c.N (current records, the MANE
    chromosome accession), plus the offset times the strand, is the table's position."""
    tm = _manifest(gene)["transcript_model"]["mane_select"]
    acc, strand = tm["grch38_chr_accession"], (1 if tm["strand"] == "+" else -1)
    exonic = re.compile(r"^N[MR]_[\d.]+:(c\.[-*]?\d+)[ACGT]>[ACGT]$")
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


def test_the_offset_sample_draws_ten_per_gene_and_every_one_agrees():
    rows = _csv("newgenes_offset_sample.csv")
    strand = {r["gene"]: r["strand"] for r in _csv("newgenes_coordinate_check.csv")}
    assert "-" in strand.values()                       # a minus-strand gene is sampled
    for gene in GENES:
        mine = [r for r in rows if r["gene"] == gene]
        assert len(mine) == N.N_SAMPLE
        assert all(int(r["exon_distance"]) == abs(int(r["offset"]))
                   <= int(r["exon_distance_other_side"]) for r in mine)
        assert all(r["anchor_agrees"] == "True" for r in mine)


# --------------------------------------------------------------------------
# counts
# --------------------------------------------------------------------------
def _counts() -> dict:
    out = {}
    for r in _csv("newgenes_counts_by_offset.csv"):
        out[(r["gene"], r["scheme"], r["bin"], r["range_class"])] = int(r["n"])
    return out


def test_each_stratum_count_is_the_sum_of_its_per_nt_counts():
    n = _counts()
    spans = {"pm12": range(1, 3), "s3_10": range(3, 11), "s11_50": range(11, 51),
             "s3_50": range(3, 51)}
    for gene in GENES:
        for cls in [c["label"] for c in _classes(gene)] + [N.UNASSIGNED]:
            for st, span in spans.items():
                assert n[(gene, "stratum", st, cls)] == sum(
                    n[(gene, "per_nt", str(d), cls)] for d in span), (gene, st, cls)
            assert n[(gene, "stratum", N.DEEP, cls)] == n[(gene, "per_nt", ">50", cls)]


def test_the_stratum_names_are_evid_commons():
    used = {r["bin"] for r in _csv("newgenes_counts_by_offset.csv") if r["scheme"] == "stratum"}
    assert used - {N.DEEP} <= set(K.STRATA)


def test_the_label_summary_is_the_tables_counted():
    rows = {(r["gene"], r["stratum"]): r for r in _csv("newgenes_label_summary.csv")}
    frames = {g: _variants(g) for g in GENES}
    frames["ALL4"] = pd.concat(frames.values(), ignore_index=True)
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


# --------------------------------------------------------------------------
# records and manifests
# --------------------------------------------------------------------------
@pytest.mark.parametrize("gene", GENES)
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
    prov = json.loads((EXT / f"{gene.lower()}_variants.parquet.provenance.json").read_text())
    assert prov["inputs"] == {e["path"]: e["sha256"] for e in man["files"]}
    assert prov["sha256"] == _sha(EXT / f"{gene.lower()}_variants.parquet")


def test_the_raw_deposits_are_stored_byte_for_byte():
    """MaveDB serves its score tables with CRLF line ends. The repository turns every
    other *.csv to LF on checkout, which would change the bytes each manifest pins."""
    if not (REPO / ".git").exists():
        pytest.skip("a release archive has no .git")
    paths = [e["path"] for g in GENES for e in _manifest(g)["files"]]
    out = subprocess.run(["git", "check-attr", "text", "--", *paths], cwd=REPO,
                         capture_output=True, text=True, check=True).stdout.splitlines()
    assert len(out) == len(paths) and all(l.endswith(": text: unset") for l in out), out[:2]


@pytest.mark.parametrize("gene", GENES)
def test_the_record_row_is_read_from_the_metadata(gene):
    meta = json.loads((EXT / f"{gene.lower()}_raw/scoreset_meta.json").read_text())
    row = {r["gene"]: r for r in _csv("newgenes_record_check.csv")}[gene]
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


def test_every_fact_on_the_record_page_names_its_file():
    """A fetched file and its time, or the search record a URN was chosen from."""
    lines = [l for l in (REPO / "docs/newgenes-record-check.md").read_text().splitlines()
             if l.startswith("- ")]
    bare = [l for l in lines if "fetched 20" not in l and "scripts/stage_E16/run_log/" not in l]
    assert lines and not bare, bare[:2]


# --------------------------------------------------------------------------
# the entry point knows the stage and its inputs
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
    assert mods[i - 1] == "src.evid_practical_metrics" and mods[i + 1] == "src.evid_delta"
    assert [p.name for p in E.NEWGENE_RAW] == [f"{g.lower()}_raw" for g in GENES]
    assert REPO / "docs" / "newgenes-record-check.md" in E.OUTPUT_ROOTS
    assert PHASE1 / "data" / "evidence" / "refs" in E.NOT_OUTPUTS


def test_every_input_the_stage_reads_is_listed_as_one():
    E = _entry_point()
    inputs = E._run_inputs()
    for gene in GENES:
        man = _manifest(gene)
        need = [e["path"] for e in man["files"]]
        need += [i["path"] for k, i in man["transcript_model"]["cache_files"].items()]
        need.append(f"phase1/data/evidence/external/{gene.lower()}_raw/manifest.json")
        assert set(need) <= inputs, sorted(set(need) - inputs)
    assert "phase1/data/evidence/external/newgenes_reference_check.json" in inputs
