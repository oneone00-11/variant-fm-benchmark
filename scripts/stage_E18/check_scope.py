"""E18.0 rule 9: a commit that changes outputs changes them only within its declared scope.

Every output file (scripts/reproduce_evidence.py's OUTPUT_ROOTS, NOT_OUTPUTS left out)
as committed at --base is compared with the same file at --head (a commit; by default
the working tree, so that a commit is checked before it is made). A file byte for byte
the same passes. A file that differs, appears or disappears must be declared in the
scope file (scripts/stage_E18/scope_<commit>.json), in one of three ways:

  "rows"    a table keyed by its declared key columns; every row added, removed or
            changed must belong to a declared (gene, tool) pair of the scope, and a
            changed row may change only the declared fields ("*": any); a table of one
            gene with no gene column declares the gene ("gene_const");
  "fields"  a table keyed the same way, of which only the declared rows' declared
            fields may change (the supplement manifest's hash of an input, say);
  "lines"   a text, of which only the lines matching a declared pattern may change,
            appear or go (a sentence derived from the external rows);
  "json"    a JSON record, of which only the declared fields may change, each a path of
            keys from the top (a description, the hash of a source the record lists);
  "file"    any change, with the reason declared (a figure that draws the external
            genes, a record whose wording changed).

The scope's pairs are written with the pipeline's tool names; a table that prints a
tool's label instead (Supplementary Table S12) is mapped through the supplement's own
label table. Cells are compared as the files print them, text for text. A file that
is declared but did not change is reported, so that the declaration says what
happened. The result, file by file, goes to scripts/stage_E18/run_log/ through the run
record; exit 1 on anything outside the scope.

Run from the repository root:
    .venv/bin/python scripts/stage_E18/check_scope.py --scope scripts/stage_E18/scope_c3.json --base <commit>
"""
from __future__ import annotations

import argparse
import io
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

sys.path.insert(0, str(REPO / "scripts"))
import reproduce_evidence as R  # noqa: E402

ROOTS = [str(r.relative_to(REPO)) for r in R.OUTPUT_ROOTS]
NOT = [str(r.relative_to(REPO)) for r in R.NOT_OUTPUTS]


def _git(*args: str) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), *args], check=True, capture_output=True).stdout


def _in_scope_paths(paths: list[str]) -> list[str]:
    return sorted(p for p in paths if p and not any(p == n or p.startswith(n + "/") for n in NOT))


def files_at(commit: str | None) -> list[str]:
    if commit:
        out = _git("ls-tree", "-r", "-z", "--name-only", commit, "--", *ROOTS)
    else:
        out = _git("ls-files", "-z", "--cached", "--others", "--exclude-standard", "--", *ROOTS)
    return _in_scope_paths([p for p in out.decode().split("\0") if p and
                            (commit or (REPO / p).is_file())])


def content(commit: str | None, path: str) -> bytes | None:
    if commit:
        try:
            return _git("show", f"{commit}:{path}")
        except subprocess.CalledProcessError:
            return None
    p = REPO / path
    return p.read_bytes() if p.is_file() else None


def _table(data: bytes | None) -> pd.DataFrame | None:
    if data is None:
        return None
    return pd.read_csv(io.BytesIO(data), dtype=str, keep_default_na=False)


def tool_labels() -> dict[str, str]:
    """label -> tool, from the supplement's own table of tool labels."""
    sys.path.insert(0, str(REPO / "phase1"))
    from src import evid_supp_tables as S
    return {v: k for k, v in S.tool_labels().items()}


def check_rows(path: str, decl: dict, old: pd.DataFrame | None, new: pd.DataFrame | None,
               pairs: dict, labels: dict) -> dict:
    key = decl["key"]
    fields = decl.get("fields", "*")
    gene_col, tool_col = decl.get("gene"), decl["tool"]
    gene_const = decl.get("gene_const")          # a table of one gene, without a gene column

    def keyed(df):
        if df is None:
            return {}
        dup = df.duplicated(key, keep=False)
        if dup.any():
            raise SystemExit(f"{path}: the declared key {key} is not unique")
        return {tuple(r[k] for k in key): r for _, r in df.iterrows()}

    def allowed(row) -> bool:
        tool = row[tool_col]
        tool = labels.get(tool, tool) if decl.get("tool_labels") else tool
        return tool in pairs.get(gene_const or row[gene_col], [])

    a, b = keyed(old), keyed(new)
    out = {"added": 0, "removed": 0, "changed": 0, "outside": []}
    for k in sorted(set(a) | set(b), key=str):
        ra, rb = a.get(k), b.get(k)
        if ra is not None and rb is not None:
            cols = sorted(set(ra.index) | set(rb.index))
            diff = [c for c in cols if ra.get(c, None) != rb.get(c, None)]
            if not diff:
                continue
            out["changed"] += 1
            bad_fields = [c for c in diff if fields != "*" and c not in fields]
            if not allowed(rb) or bad_fields:
                out["outside"].append({"key": list(k), "fields": diff})
        else:
            row = rb if rb is not None else ra
            out["added" if rb is not None else "removed"] += 1
            if not allowed(row):
                out["outside"].append({"key": list(k), "row": "added" if rb is not None else "removed"})
    if old is not None and new is not None and list(old.columns) != list(new.columns):
        added = [c for c in new.columns if c not in old.columns]
        removed = [c for c in old.columns if c not in new.columns]
        if (added or removed) and not decl.get("columns_may_change"):
            out["outside"].append({"columns_added": added, "columns_removed": removed})
    return out


def check_fields(path: str, decl: dict, old: pd.DataFrame | None, new: pd.DataFrame | None) -> dict:
    key = decl["key"]
    ka = {tuple(r[k] for k in key): r for _, r in old.iterrows()}
    kb = {tuple(r[k] for k in key): r for _, r in new.iterrows()}
    out = {"changed": 0, "outside": []}
    if set(ka) != set(kb) or list(old.columns) != list(new.columns):
        out["outside"].append({"rows_or_columns_differ": True})
        return out
    allowed = {tuple([r]) if isinstance(r, str) else tuple(r): f for r, f in decl["rows"].items()}
    parts_from = decl.get("parts_from", "Source files")
    for k in ka:
        diff = [c for c in old.columns if ka[k][c] != kb[k][c]]
        if not diff:
            continue
        out["changed"] += 1
        spec = allowed.get(k, [])
        for c in diff:
            if isinstance(spec, dict) and c in spec:
                # a "; "-joined list, one entry per source named in parts_from: only the
                # declared sources' entries may change
                names = [x.strip() for x in str(kb[k][parts_from]).split(";")]
                pa, pb = str(ka[k][c]).split("; "), str(kb[k][c]).split("; ")
                if not (len(pa) == len(pb) == len(names)) or ka[k][parts_from] != kb[k][parts_from]:
                    out["outside"].append({"key": list(k), "field": c, "parts": "list changed"})
                    continue
                moved = [n for n, x, y in zip(names, pa, pb) if x != y]
                if not set(moved) <= set(spec[c]):
                    out["outside"].append({"key": list(k), "field": c, "parts": moved})
            elif c not in spec:
                out["outside"].append({"key": list(k), "fields": [c]})
    return out


def check_json(decl: dict, old: bytes | None, new: bytes | None) -> dict:
    allowed = [tuple(p) for p in decl["paths"]]
    out = {"changed": 0, "outside": []}

    def walk(a, b, path=()):
        if isinstance(a, dict) and isinstance(b, dict):
            for k in sorted(set(a) | set(b)):
                yield from walk(a.get(k), b.get(k), path + (k,))
        elif a != b:
            yield path

    if old is None or new is None:
        out["outside"].append({"record": "added" if old is None else "removed"})
        return out
    for path in walk(json.loads(old), json.loads(new)):
        out["changed"] += 1
        if path not in allowed:
            out["outside"].append({"path": list(path)})
    return out


def check_lines(decl: dict, old: bytes | None, new: bytes | None) -> dict:
    import difflib
    import re
    pats = [re.compile(p) for p in decl["patterns"]]
    a = (old or b"").decode().splitlines()
    b = (new or b"").decode().splitlines()
    out = {"changed": 0, "outside": []}
    for line in difflib.unified_diff(a, b, lineterm="", n=0):
        if line.startswith(("---", "+++", "@@")):
            continue
        out["changed"] += 1
        if not any(p.search(line[1:]) for p in pats):
            out["outside"].append({"line": line[:200]})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scope", required=True)
    ap.add_argument("--base", required=True)
    ap.add_argument("--head", default=None, help="a commit; the working tree if left out")
    args = ap.parse_args()
    scope = json.loads(Path(args.scope).read_text())
    pairs = scope["pairs"]
    with RunLog(__file__) as log:
        log.input(Path(args.scope))
        labels = tool_labels()
        paths = sorted(set(files_at(args.base)) | set(files_at(args.head)))
        report, failed = {}, []
        for p in paths:
            a, b = content(args.base, p), content(args.head, p)
            if a == b:
                continue
            decl = scope["files"].get(p)
            if decl is None:
                report[p] = {"status": "OUTSIDE: not declared",
                             "change": "added" if a is None else "removed" if b is None else "differs"}
                failed.append(p)
                continue
            kind = decl["kind"]
            if kind == "file":
                report[p] = {"status": "declared", "why": decl["why"],
                             "change": "added" if a is None else "removed" if b is None else "differs"}
                continue
            if kind == "rows":
                r = check_rows(p, decl, _table(a), _table(b), pairs, labels)
            elif kind == "lines":
                r = check_lines(decl, a, b)
            elif kind == "json":
                r = check_json(decl, a, b)
            else:
                r = check_fields(p, decl, _table(a), _table(b))
            r["status"] = "OUTSIDE" if r["outside"] else "within the scope"
            report[p] = r
            if r["outside"]:
                failed.append(p)
        unchanged = [p for p in scope["files"] if p not in report]
        log.note("scope", rel(Path(args.scope)))
        log.note("base", args.base)
        log.note("head", args.head or "working tree")
        log.note("files", report)
        log.note("declared_but_unchanged", unchanged)
        log.note("outside", failed)
        for p, r in report.items():
            print(f"  {r['status']:24s} {p}"
                  + (f"  (+{r.get('added', 0)} -{r.get('removed', 0)} ~{r.get('changed', 0)})"
                     if "changed" in r else ""))
            for o in r.get("outside", [])[:5]:
                print(f"      {o}")
        for p in unchanged:
            print(f"  {'declared, unchanged':24s} {p}")
        print(f"\n{len(report)} output file(s) differ from {args.base}; "
              f"{len(failed)} outside the scope of {scope['commit']}")
        if failed:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
