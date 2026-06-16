"""Classify a functional variant into a region_class from its transcript HGVS,
using the SAME label vocabulary as Milestone 1 plus an explicit `intronic` split.

Labels: missense | synonymous | nonsense | splice | intronic | noncoding
        | frameshift | inframe_indel | startloss | coding_other | unknown

Splice vs intronic: an intronic offset |n| <= 2 (the essential GT..AG donor/
acceptor dinucleotides) is `splice`; deeper intronic is `intronic`. The signed
offset is also returned so the threshold can be refined later.
"""
import re

# c. position token: capture the coordinate part between "c." and the change op
_CDOT = re.compile(r":c\.([\*\-]?\d+(?:[+\-]\d+)?(?:_[\*\-]?\d+(?:[+\-]\d+)?)?)")
_INTRON_OFFSET = re.compile(r"\d+([+\-]\d+)")


def _intron_offset(cpos: str):
    """Return signed intron offset if cpos is intronic (e.g. '5074+1' -> 1,
    '5075-2' -> -2), else None. A leading '-' (5'UTR) or '*' (3'UTR) is NOT
    an intron offset."""
    # strip a possible range; take the first endpoint
    first = cpos.split("_")[0]
    m = _INTRON_OFFSET.search(first)
    if m:
        return int(m.group(1))
    return None


def classify_region(hgvs_nt: str, hgvs_pro: str = "", hgvs_splice: str = ""):
    """Return (region_class, intron_offset_or_None)."""
    nt = (hgvs_nt or "").strip()
    pro = (hgvs_pro or "").strip()

    # Prefer the splice-coordinate HGVS when the nt one is plain (some MaveDB
    # rows carry the intronic coordinate only in hgvs_splice).
    src = nt
    m = _CDOT.search(nt)
    if not m and hgvs_splice and hgvs_splice not in ("NA", ""):
        src = hgvs_splice
        m = _CDOT.search(hgvs_splice)

    if not m:
        # n. (non-coding transcript) or unparseable
        if ":n." in nt or ":n." in (hgvs_splice or ""):
            return "noncoding", None
        return "unknown", None

    cpos = m.group(1)

    # UTRs
    if cpos.startswith("*"):
        return "noncoding", None          # 3' UTR
    if cpos.startswith("-"):
        return "noncoding", None          # 5' UTR

    # Intronic?
    off = _intron_offset(cpos)
    if off is not None:
        return ("splice" if abs(off) <= 2 else "intronic"), off

    # Exonic coding -> use protein consequence
    if pro and pro not in ("NA", ""):
        p = pro
        if p.endswith("=") or "(=)" in p or p == "p.=":
            return "synonymous", None
        if "Ter" in p or p.endswith("*"):
            # start-loss sometimes written p.Met1? ; stop-gain ends with Ter
            return "nonsense", None
        if "Met1?" in p or p.endswith("Met1?"):
            return "startloss", None
        if re.search(r"p\.[A-Za-z]{3}\d+[A-Za-z]{3}", p):
            return "missense", None
        if "del" in p or "ins" in p or "dup" in p:
            return "inframe_indel", None
        return "coding_other", None

    # No protein annotation: infer from nt change shape
    if "del" in nt or "dup" in nt or "ins" in nt:
        # rough: frameshift vs inframe needs length; mark generic
        return "coding_other", None
    return "coding_other", None
