"""E19 C1, low circularity -- the ClinVar check on assertions that could not have used the SGE data.

The subset keeps a variant's ClinVar assertion only if ClinVar last evaluated it
before its gene's SGE results were first public (E19.4); such an assertion cannot have
used those data. It removes that circularity and no other: some ClinVar classifications
may have used predictor scores as PP3 or BP4, which no date removes, and the distal band
holds almost no usable P/LP.

The first public date of each training gene is the earliest of three dates, each read
offline from a record saved as served (scripts/stage_E19/fetch_first_public_dates.py;
phase1/data/evidence/e19/first_public_dates/):
  * MaveDB: the publishedDate of the earliest version of the cited score set
    (Supplementary Table S1), walking supersededScoreSet back; where the chain names an
    earlier version MaveDB does not serve, the earliest version's date is not known and
    this source gives none;
  * the journal: Crossref's published-online date of a journal article whose DOI MaveDB
    records; where Crossref gives no published-online date, its `created` date (the
    author's decision of 2026-10-25; in the records that carry both, `created` is the
    online date or the day before it);
  * a preprint: the earlier of Crossref's `posted` and `created` dates (the author's
    review of 2026-10-25, item 1), for a preprint DOI that MaveDB records, that a
    "has-preprint" relation of a Crossref record lists, or that Europe PMC's preprint
    index gives for the study: the searches of scripts/stage_E19/fetch_preprints.py
    (first_public_dates/preprint_search.json), every gene searched, not left to the
    relation fields.
A preprint those searches return is the study's when its title names the gene (as a
word), its title or abstract names saturation genome editing, and its authors include
both the first and the last author of the paper MaveDB records for the score set
(family names, case and accents folded) or, for a score set that records no paper, a
contributor a version of it names. Copies of one study share its title (markup, case
and punctuation folded); preprints of more than one title for a gene stop the stage for
the author (E19.12), as does a gene without any date.

The date of an assertion is LastEvaluated in ClinVar's monthly archive of June 2026
(the committed extract, scripts/stage_E19/extract_clinvar_variant_summary.py), joined on
VariationID. The archive (2026-06-04) and the VCF the set was built from (2026-06-15)
are compared variant by variant on two fields, each normalised and then required to be
identical (the author's decision of 2026-10-25): the classification string (the
archive's ClinicalSignificance against the VCF's CLNSIG) and the review status
(ReviewStatus against CLNREVSTAT). The normalisation: underscores read as spaces, runs of
white space as one, and an absent value -- the VCF's empty field, the archive's "-" -- as
the same absent value. A variant that differs on either is left out of the subset and
counted; under 98% agreement the stage stops (E19.4).

The subset: P/LP or B/LB with at least one star, agreeing in both files, with a
LastEvaluated date strictly before the gene's first public date. On it, E19.3's metrics
are computed again by evid_suppl_c1's own functions, at the main level only.

Outputs, in reports/evidence/suppl_analyses/:
  c1_preprint_search.csv             every preprint the searches returned for each gene,
                                     the three conditions and whether it is the study's;
                                     a row per gene says which was found, or none
  c1_first_public_dates.csv          each training gene's dates by source and the
                                     earliest, with the record each is read from
  c1_low_circularity_variants.csv    every variant with a ClinVar record: the two files'
                                     fields, whether they agree, the dates, whether in the
                                     subset and if not why
  c1_low_circularity_composition.csv P/LP and B/LB in the subset by gene and band, beside
                                     the main analysis's
  c1_low_circularity_metrics.csv     E19.3's metrics on the subset

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_c1_low
"""
from __future__ import annotations

import csv
import html
import json
import re
import unicodedata
import urllib.parse
from datetime import date, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import evid_common as K
from . import evid_suppl_c1 as C1
from . import evid_territory_metrics as T

E19 = Path("data/evidence/e19")
DATES = E19 / "first_public_dates"
ARCHIVE = E19 / "clinvar_variant_summary"
EXTRACT = ARCHIVE / "variant_summary_2026-06_analysis_set.tsv"
OUT_DIR = C1.OUT_DIR
OUT_DATES = OUT_DIR / "c1_first_public_dates.csv"
OUT_PREPRINTS = OUT_DIR / "c1_preprint_search.csv"
PREPRINT_SEARCH = DATES / "preprint_search.json"
SGE_WORDS = "saturation genome editing"
OUT_VARIANTS = OUT_DIR / "c1_low_circularity_variants.csv"
OUT_COMPOSITION = OUT_DIR / "c1_low_circularity_composition.csv"
OUT_METRICS = OUT_DIR / "c1_low_circularity_metrics.csv"
MIN_AGREEMENT = 0.98
MAIN = [lvl for lvl in C1.STAR_LEVELS if lvl[2] == "main"]


def _repo_path(p: str) -> Path:
    """A path the manifests record from the repository root, read from phase1/."""
    return Path("..") / p


def _ymd(parts: list | None, what: str) -> date | None:
    """A Crossref date given to the day; absent, None. A date given to the month or the
    year only is not a day, and is not filled in: it stops the stage (E19.12)."""
    if not parts or not parts[0] or parts[0][0] is None:
        return None
    if len(parts[0]) < 3:
        raise SystemExit(f"[C1 low] {what} is given as {parts[0]}, not to the day: stopped")
    y, m, d = parts[0][:3]
    return date(int(y), int(m), int(d))


def _fold(s: str) -> str:
    s = unicodedata.normalize("NFKD", s or "")
    return "".join(c for c in s if not unicodedata.combining(c)).lower().strip()


def _plain(s: str) -> str:
    """Text without its markup and entities, white space collapsed."""
    return " ".join(html.unescape(re.sub(r"<[^>]+>", " ", s or "")).split())


def _title_key(s: str) -> str:
    return " ".join(re.sub(r"[^0-9a-z]+", " ", _fold(_plain(s))).split())


def _preprint_date(m: dict, doi: str) -> tuple[str, str]:
    """The earlier of a preprint record's posted and created dates, and which it is."""
    got = [(d, f) for f, d in (("posted", _ymd((m.get("posted") or {}).get("date-parts"),
                                                f"{doi} posted")),
                               ("created", _ymd((m.get("created") or {}).get("date-parts"),
                                                f"{doi} created"))) if d is not None]
    if not got:
        return None, "neither posted nor created"
    d, f = min(got)
    return d.isoformat(), f"{f}, the earlier of posted and created"


def _records() -> tuple[dict, dict, dict]:
    """The two fetch records, and every saved file by the URL it was fetched from."""
    man = json.loads((DATES / "manifest.json").read_text())
    ps = json.loads(PREPRINT_SEARCH.read_text())
    saved = {f["url"]: f.get("file") for f in man["files"] + ps["files"]}
    return man, ps, saved


def _crossref(man: dict, saved: dict, doi: str) -> tuple[dict | None, str | None]:
    path = saved.get(man["crossref_api"] + urllib.parse.quote(doi, safe="/"))
    if path is None:
        return None, None
    return json.loads(_repo_path(path).read_text())["message"], path


def preprint_matches() -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Every preprint the Europe PMC searches returned, per gene, against the rule of the
    module's docstring; and each gene's matched preprint DOIs."""
    man, ps, saved = _records()
    rows, found = [], {}
    for gene, searches in ps["searches"].items():
        g = man["genes"][gene]
        papers = [d for d in g["dois"] if not d.startswith("10.1101/")]
        need, contributors = [], set()
        for doi in papers:
            m, _ = _crossref(man, saved, doi)
            people = [_fold(a["family"]) for a in m.get("author") or [] if a.get("family")]
            need.append({people[0], people[-1]})
        if not papers:
            for urn in g["versions"]:
                meta = json.loads(_repo_path(saved[f"{man['mavedb_api']}/score-sets/{urn}"])
                                  .read_text())
                contributors |= {_fold(c["familyName"]) for c in meta.get("contributors") or []
                                 if c.get("familyName")}
        hits = {}
        for f in ps["files"]:
            if f.get("file") and f["url"].startswith(ps["europepmc_api"]) and \
                    Path(f["file"]).name.startswith(f"{gene}_"):
                for h in json.loads(_repo_path(f["file"]).read_text())["resultList"]["result"]:
                    hits.setdefault(h["id"], h)
        word = re.compile(rf"(?<![A-Za-z0-9]){re.escape(gene)}(?![A-Za-z0-9])")
        for pid, h in sorted(hits.items()):
            title = _plain(h.get("title"))
            text = f"{title} {_plain(h.get('abstractText'))}".lower()
            authors = {_fold(a.get("lastName", "")) for a in
                       (h.get("authorList") or {}).get("author") or [] if a.get("lastName")}
            c_gene = bool(word.search(title))
            c_sge = SGE_WORDS in text
            c_auth = (any(n <= authors for n in need) if papers else bool(authors & contributors))
            rows.append({"gene": gene, "preprint": pid, "doi": (h.get("doi") or "").lower(),
                         "server": (h.get("bookOrReportDetails") or {}).get("publisher", ""),
                         "europepmc_first_publication_date": h.get("firstPublicationDate"),
                         "title": title, "title_names_gene": c_gene,
                         "names_saturation_genome_editing": c_sge,
                         "authors_condition": ("first and last author of " + ", ".join(papers)
                                               if papers else "a contributor named on the score set"),
                         "authors_met": c_auth, "is_the_studys": c_gene and c_sge and c_auth})
        mine = [r for r in rows if r["gene"] == gene and r["is_the_studys"]]
        titles = {_title_key(r["title"]) for r in mine}
        if len(titles) > 1:
            raise SystemExit(f"[C1 low] {gene}: preprints of {len(titles)} titles meet the rule "
                             f"({[r['doi'] for r in mine]}): stopped for the author (E19.12)")
        found[gene] = sorted({r["doi"] for r in mine if r["doi"]})
        rows.append({"gene": gene, "preprint": "", "doi": ";".join(found[gene]),
                     "title": ("found: " + "; ".join(found[gene])) if found[gene] else "not found",
                     "is_the_studys": bool(found[gene])})
    return pd.DataFrame(rows), found


def first_public_dates(found: dict[str, list[str]]) -> pd.DataFrame:
    man, ps, saved = _records()
    rows = []
    for gene, g in man["genes"].items():
        dated = []
        # MaveDB: the earliest version, if MaveDB serves the whole chain
        if g["versions_not_served"]:
            rows.append({"gene": gene, "source": "mavedb_earliest_version",
                         "identifier": g["versions_not_served"][-1], "date": None,
                         "field": "the chain's earliest version is not served by MaveDB "
                                  "(HTTP 404)", "record": ""})
        else:
            oldest = g["versions"][-1]
            path = saved[f"{man['mavedb_api']}/score-sets/{oldest}"]
            rec = json.loads(_repo_path(path).read_text())
            dated.append(("mavedb_earliest_version", oldest, rec["publishedDate"],
                          "publishedDate", path))
        # every DOI MaveDB records, the relations' preprints, and the searches' preprints
        routes = {d: ["MaveDB record or Crossref relation"] for d in g["dois"]}
        for d in found.get(gene, []):
            routes.setdefault(d, []).append("Europe PMC preprint search")
        for doi, how in routes.items():
            m, path = _crossref(man, saved, doi)
            if m is None:
                rows.append({"gene": gene, "source": "crossref", "identifier": doi,
                             "date": None, "field": "no Crossref record (HTTP 404)",
                             "record": "", "found_by": "; ".join(how)})
                continue
            if m.get("type") == "posted-content":
                d, field = _preprint_date(m, doi)
                dated.append(("preprint", doi, d, field, path, "; ".join(how)))
            else:
                d = _ymd((m.get("published-online") or {}).get("date-parts"),
                         f"{doi} published-online")
                field = "published-online"
                if d is None:
                    d = _ymd((m.get("created") or {}).get("date-parts"), f"{doi} created")
                    field = "created (no published-online; the author's decision of 2026-10-25)"
                dated.append(("journal_online", doi, d and d.isoformat(), field, path,
                              "; ".join(how)))
        for x in dated:
            rows.append({"gene": gene, "source": x[0], "identifier": x[1], "date": x[2],
                         "field": x[3], "record": x[4], "found_by": x[5] if len(x) > 5 else ""})
        have = [r for r in rows if r["gene"] == gene and r["date"]]
        if not have:
            raise SystemExit(f"[C1 low] {gene}: no record gives a first public date: "
                             "stopped (E19.12)")
        first = min(have, key=lambda r: r["date"])
        rows.append({"gene": gene, "source": "first_public_date",
                     "identifier": first["identifier"], "date": first["date"],
                     "field": f"the earliest: {first['source']}", "record": first["record"],
                     "found_by": first.get("found_by", "")})
    return pd.DataFrame(rows)


def _norm(v) -> str:
    """The normalisation of E19.4's comparison, as the docstring states it."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        v = ""
    v = " ".join(str(v).replace("_", " ").split())
    return "" if v in ("", "-") else v


def _evaluated(v: str) -> date | None:
    v = (v or "").strip()
    if v in ("", "-"):
        return None
    return datetime.strptime(v, "%b %d, %Y").date()


def low_circularity_variants(df: pd.DataFrame, dates: pd.DataFrame) -> pd.DataFrame:
    x = K.read_back(EXTRACT, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)
    if x["VariationID"].duplicated().any():
        raise SystemExit("[C1 low] a VariationID has two GRCh38 rows in the extract")
    rec = df[df["variation_id"].notna()].copy()
    rec = rec.merge(x[["VariationID", "ClinicalSignificance", "ReviewStatus", "LastEvaluated"]],
                    left_on="variation_id", right_on="VariationID", how="left")
    if rec["VariationID"].isna().any():
        raise SystemExit("[C1 low] a ClinVar record of the set is not in the archive extract")
    first = dates[dates.source == "first_public_date"].set_index("gene")["date"]
    rec["same_classification"] = rec["clnsig_raw"].map(_norm) == rec["ClinicalSignificance"].map(_norm)
    rec["same_review_status"] = rec["clnrevstat"].map(_norm) == rec["ReviewStatus"].map(_norm)
    rec["agree"] = rec["same_classification"] & rec["same_review_status"]
    rec["last_evaluated"] = rec["LastEvaluated"].map(_evaluated)
    rec["first_public_date"] = rec["gene"].map(first)
    reason = []
    for r in rec.itertuples():
        if r.clnsig_class not in ("P/LP", "B/LB"):
            reason.append("not P/LP or B/LB")
        elif not (r.clinvar_stars >= 1):
            reason.append("fewer than one star")
        elif not r.agree:
            reason.append("archive and VCF differ")
        elif r.last_evaluated is None:
            reason.append("no LastEvaluated date")
        elif not (r.last_evaluated < date.fromisoformat(r.first_public_date)):
            reason.append("evaluated on or after the first public date")
        else:
            reason.append("")
    rec["in_subset"] = [not w for w in reason]
    rec["why_not"] = reason
    cols = ["variant_id", "gene", "stratum", "variation_id", "clnsig_raw", "ClinicalSignificance",
            "same_classification", "clnrevstat", "ReviewStatus", "same_review_status", "agree",
            "clnsig_class", "clinvar_stars", "LastEvaluated", "last_evaluated",
            "first_public_date", "in_subset", "why_not"]
    return rec[cols]


def composition(df: pd.DataFrame, var: pd.DataFrame) -> pd.DataFrame:
    rows = []
    main = df[df["y_clinvar"].notna() & (df["clinvar_stars"] >= 1)]
    sub = var[var["in_subset"]]
    for gene in sorted(df["gene"].unique()) + ["Seven genes"]:
        for stratum in C1.BANDS:
            m = K.stratum_frame(main if gene == "Seven genes" else main[main.gene == gene], stratum)
            s = K.stratum_frame(sub if gene == "Seven genes" else sub[sub.gene == gene], stratum)
            rows.append({"gene": gene, "stratum": stratum,
                         "main_n_plp": int((m.clnsig_class == "P/LP").sum()),
                         "main_n_blb": int((m.clnsig_class == "B/LB").sum()),
                         "subset_n_plp": int((s.clnsig_class == "P/LP").sum()),
                         "subset_n_blb": int((s.clnsig_class == "B/LB").sum())})
    return pd.DataFrame(rows)


def main() -> None:
    cfg = yaml.safe_load(C1.CONFIG_PATH.read_text())
    df = K.load_set()
    tools = K.panel_of(df)
    df[K.FUSION] = K.logo_fusion(df, [t for t in K.FUSION_FEATURES if t in df.columns])
    tools = tools + [K.FUSION]
    df["clinvar_stars"] = C1.stars_of(df, C1.star_table())
    df["y_clinvar"] = df["clnsig_class"].map({"P/LP": 1.0, "B/LB": 0.0})
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    search, found = preprint_matches()
    search.to_csv(OUT_PREPRINTS, index=False)
    dates = first_public_dates(found)
    dates.to_csv(OUT_DATES, index=False)
    var = low_circularity_variants(df, dates)
    rate = float(var["agree"].mean())
    var.to_csv(OUT_VARIANTS, index=False)
    print(f"[C1 low] archive and VCF agree on {int(var.agree.sum())} of {len(var)} records "
          f"({rate:.4f}); {int((~var.agree).sum())} left out")
    if rate < MIN_AGREEMENT:
        raise SystemExit(f"[C1 low] agreement {rate:.4f} is under {MIN_AGREEMENT}: stopped (E19.4)")
    keep = set(var.loc[var["in_subset"], "variant_id"])
    sub = df.assign(y_clinvar=df["y_clinvar"].where(df["variant_id"].isin(keep)))
    composition(df, var).to_csv(OUT_COMPOSITION, index=False)
    ev = K.read_back(C1.REPORT_DIR / "evidence_thresholds.csv")
    m = C1.reference_metrics(sub, tools, ev, cfg, T._atlas_clinical_evidence(), levels=MAIN)
    m.to_csv(OUT_METRICS, index=False)
    print(f"[C1 low] wrote {OUT_DATES}, {OUT_VARIANTS}, {OUT_COMPOSITION} and {OUT_METRICS} "
          f"({len(m)} rows); {len(keep)} variants in the subset")
    print(dates[dates.source == "first_public_date"][["gene", "date", "field"]].to_string(index=False))


if __name__ == "__main__":
    main()
