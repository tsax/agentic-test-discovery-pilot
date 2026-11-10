"""GRIM baseline (Brown & Heathers 2017), implemented from the published definition.

A mean M reported to d decimals, computed from n observations of an integer-valued
measure, is GRIM-consistent iff there exists an integer total T with
round(T / n, d) == M.  GRIM only has power when n < 10**d.

For a k-item scale whose score is the mean of k integer items, the granularity is
1 / (n * k), so the effective n is n * k.

Rounding: we accept a match under round-half-up, round-half-even or round-half-down,
so that a reporting convention difference never produces a flag.

Inputs are found in tables: a column (or row) labelled Mean / M / "Mean (SD)" /
"M ± SD", with n available from (in priority order) the row label "(n = ..)",
an N/n column, an "(n = ..)" in the column header, a row labelled N/n, or a single
unambiguous "N = .." in the caption/footnotes.

Conservatism: an item is tested only if there is positive evidence that the
measure is integer-valued (Likert / k-point / rating / count / number of / items /
total or sum score / age in years) and no evidence it is continuous (%, ms, kg,
BMI, time, adjusted/estimated means, ...).  A multi-item scale with unknown k is
skipped.  Every candidate gets a recorded status and reason (see `analyze`).

Interface (docs/SCHEMA.md §3): check(paper) -> list[flag]; applicable(paper) -> bool.
Stdlib only, deterministic, offline.
"""
from __future__ import annotations

import math
import re
from decimal import Decimal, ROUND_HALF_DOWN, ROUND_HALF_EVEN, ROUND_HALF_UP, InvalidOperation

# ---------------------------------------------------------------- core GRIM test

_ROUNDINGS = (ROUND_HALF_UP, ROUND_HALF_EVEN, ROUND_HALF_DOWN)


def decimals_of(mean_str: str) -> int:
    s = mean_str.strip()
    return len(s.split(".", 1)[1]) if "." in s else 0


def grim_has_power(n: int, d: int) -> bool:
    return n >= 1 and d >= 1 and n < 10 ** d


def grim_consistent(mean_str: str, n: int) -> bool | None:
    """True/False for consistency; None if GRIM is not applicable (no power / bad input)."""
    try:
        m = Decimal(mean_str.strip().replace("−", "-"))
    except InvalidOperation:
        return None
    d = decimals_of(mean_str)
    if not grim_has_power(n, d):
        return None
    q = Decimal(1).scaleb(-d)
    approx = m * n
    lo, hi = math.floor(approx) - 1, math.ceil(approx) + 1
    for total in range(lo, hi + 1):
        x = Decimal(total) / Decimal(n)
        for mode in _ROUNDINGS:
            if x.quantize(q, rounding=mode) == m:
                return True
    return False


# ---------------------------------------------------------------- table parsing

_NUM = r"[−\-]?\d+(?:\.\d+)?"
_LEAD_NUM = re.compile(r"^\s*(" + _NUM + r")(?=$|[\s(±\[,;/*†‡§a-zA-Z])")
_N_EQ = re.compile(r"(?<![A-Za-z])[nN]\s*[=:]\s*(\d{1,3}(?:,\d{3})*|\d+)(?!\.\d)")

_MEAN_HDR_EXCLUDE = re.compile(
    r"median|mean\s*diff|difference|mean\s*rank|adjust|ls\s*mean|least|geometric|estimated|"
    r"marginal|predicted|change|\bse\b|s\.e\.|sem\b|log",
    re.I,
)
_MEAN_HDR = re.compile(
    r"^\s*(?:mean|average|avg\.?|m)\s*(?:$|[(\[±+/,]|\(sd\)|sd\b)|\bmeans?\b|\bM\s*\(\s*SD\s*\)|\bM\s*±\s*SD",
    re.I,
)

_INT_CUES = re.compile(
    r"likert|\b\d+\s*-?\s*point\b|rating|number\s+of|\bno\.\s+of|\bcount|\bitems?\b|"
    r"total\s+score|sum\s+score|sum\s+of|age\s*\(\s*y(?:ea)?rs?\s*\)|age\s+in\s+years|"
    r"age,?\s+years|\bscore\b",
    re.I,
)
_NONINT_CUES = re.compile(
    r"%|percent|\bms\b|millisec|\bsec\b|seconds|\bmin\b|minutes|hours|\bkg\b|\bcm\b|\bmm\b|"
    r"\bm2\b|bmi|height|weight|time|duration|\brt\b|reaction|latency|income|salary|ratio|"
    r"\brate\b|proportion|\blog|z-?score|t-?score|standardi[sz]|adjust|estimated|marginal|"
    r"predicted|residual|\bindex\b|mmhg|\bmg\b|\bml\b|\bhz\b|°|temperature|concentration|"
    r"\bdose|\bse\b|s\.e\.|accuracy|d['′]|\bpercentile|\$|€|£|usd|wage|distance|speed|volume|"
    r"\bfactor\b|loading|\bbeta\b|coefficient|\bor\b\s*\(|odds",
    re.I,
)
_MULTI_ITEM_CUES = re.compile(
    r"subscale|sub-scale|composite|cronbach|\balpha\b|α|averaged|average\s+of|mean\s+of\s+(?:the\s+)?items|"
    r"\bscales\b|questionnaire|inventory|dimension",
    re.I,
)
_K_ITEMS = re.compile(r"\b(\d{1,3})\s*-?\s*items?\b", re.I)
_AGE = re.compile(r"\bage\b", re.I)
_AGE_YEARS = re.compile(r"age\s*\(\s*y(?:ea)?rs?\s*\)|age\s+in\s+years|age,?\s+years|\byears\b", re.I)
_N_LABEL = re.compile(r"^\s*(?:n|N|no\.?|sample\s*size|participants|subjects)\s*$", re.I)
_K_LABEL = re.compile(r"^\s*(?:k|items|no\.?\s*(?:of\s*)?items|number\s+of\s+items|#\s*items)\s*$", re.I)


def _parse_int(s: str) -> int | None:
    s = s.strip().replace(",", "")
    s = re.sub(r"[*†‡§a-z]+$", "", s)
    return int(s) if re.fullmatch(r"\d+", s) else None


def _lead_mean(cell: str) -> str | None:
    m = _LEAD_NUM.match(cell.replace("−", "-"))
    if not m:
        return None
    v = m.group(1)
    return v if "." in v else None  # d == 0 never has power


def _single_n(text: str) -> int | None:
    vals = {int(v.replace(",", "")) for v in _N_EQ.findall(text or "")}
    return vals.pop() if len(vals) == 1 else None


def _col_header(table: dict, c: int) -> str:
    parts = []
    for hr in table.get("header_rows") or []:
        if c < len(hr) and hr[c] and hr[c] not in parts:
            parts.append(hr[c])
    return " / ".join(parts)


def _is_mean_label(text: str) -> bool:
    if not text or _MEAN_HDR_EXCLUDE.search(text):
        return False
    return bool(_MEAN_HDR.search(text))


def _integer_judgement(label: str, header: str, caption: str) -> tuple[bool, int | None, str]:
    """Return (ok, k_items, reason). k_items=None means 1 (single integer item / sum)."""
    local = f"{label} || {header}"
    if _NONINT_CUES.search(local):
        return False, None, f"continuous cue in label/header: {_NONINT_CUES.search(local).group(0)!r}"
    if _AGE.search(label) and not _AGE_YEARS.search(local + " " + caption):
        return False, None, "age without 'years' unit"
    k = None
    for src in (label, header, caption):
        km = _K_ITEMS.search(src)
        if km:
            k = int(km.group(1))
            break
    has_local_int = bool(_INT_CUES.search(local)) or bool(_AGE.search(label))
    has_cap_int = bool(_INT_CUES.search(caption))
    if not (has_local_int or has_cap_int):
        return False, None, "no evidence measure is integer-valued"
    if _NONINT_CUES.search(caption) and not has_local_int:
        return False, None, "caption has continuous cue and label has no integer cue"
    if k is None and not _AGE.search(label):
        if _MULTI_ITEM_CUES.search(local + " " + caption) or re.search(r"\bitems\b", local + " " + caption, re.I):
            if not re.search(r"total|sum", local, re.I):
                return False, None, "multi-item scale with unknown number of items"
    return True, k, "integer-valued cue found"


def _find_candidates(table: dict):
    """Yield dicts: row, col, mean_str, n, n_source, label, header."""
    rows = table.get("rows") or []
    if not rows:
        return
    ncol = max(len(r) for r in rows)
    headers = [_col_header(table, c) for c in range(ncol)]
    caption = " ".join(x for x in (table.get("caption") or "", table.get("footnotes") or "") if x)
    global_n = _single_n(caption)

    n_cols = [c for c in range(ncol) if _N_LABEL.match(headers[c] or "")]
    k_cols = [c for c in range(ncol) if _K_LABEL.match(headers[c] or "")]
    n_row_idx = [i for i, r in enumerate(rows) if r and _N_LABEL.match(r[0] or "")]
    header_n = [_single_n(h) for h in headers]

    def row_n(i, c):
        r = rows[i]
        lab = r[0] if r else ""
        v = _single_n(lab)
        if v:
            return v, "row label"
        for nc in n_cols:
            if nc < len(r):
                v = _parse_int(r[nc])
                if v:
                    return v, "N column"
        if header_n[c]:
            return header_n[c], "column header"
        for ni in n_row_idx:
            if c < len(rows[ni]):
                v = _parse_int(rows[ni][c])
                if v:
                    return v, "N row"
        if global_n:
            return global_n, "caption/footnote"
        return None, None

    def row_k(i):
        for kc in k_cols:
            if kc < len(rows[i]):
                v = _parse_int(rows[i][kc])
                if v:
                    return v
        return None

    mean_cols = [c for c in range(1, ncol) if _is_mean_label(headers[c])]
    for c in mean_cols:
        for i, r in enumerate(rows):
            if c >= len(r) or i in n_row_idx:
                continue
            ms = _lead_mean(r[c])
            if ms is None:
                continue
            n, src = row_n(i, c)
            yield dict(row=i, col=c, mean_str=ms, n=n, n_source=src, label=r[0] if r else "",
                       header=headers[c], caption=caption, k_col=row_k(i))
    if mean_cols:
        return
    # Row-wise: a row whose first cell is Mean/M; columns are groups/measures.
    for i, r in enumerate(rows):
        if not r or not _is_mean_label(r[0]):
            continue
        for c in range(1, len(r)):
            if _N_LABEL.match(headers[c] or "") or re.fullmatch(r"\s*(%|sd|s\.d\.|se)\s*", headers[c] or "", re.I):
                continue
            ms = _lead_mean(r[c])
            if ms is None:
                continue
            n, src = row_n(i, c)
            yield dict(row=i, col=c, mean_str=ms, n=n, n_source=src, label=headers[c],
                       header=r[0], caption=caption, k_col=None)


# ---------------------------------------------------------------- public API

def analyze(paper: dict) -> list[dict]:
    """Every candidate mean with status in {consistent, inconsistent, skipped} and a reason."""
    out = []
    for t in paper.get("tables") or []:
        for cand in _find_candidates(t):
            rec = {"table_id": t.get("table_id"), "row": cand["row"], "col": cand["col"],
                   "mean": cand["mean_str"], "n": cand["n"], "n_source": cand["n_source"],
                   "label": cand["label"], "header": cand["header"]}
            if not cand["n"]:
                rec.update(status="skipped", reason="no n found")
                out.append(rec)
                continue
            ok, k, why = _integer_judgement(cand["label"], cand["header"], cand["caption"])
            k = cand["k_col"] or k
            if not ok:
                rec.update(status="skipped", reason=why)
                out.append(rec)
                continue
            n_eff = cand["n"] * (k or 1)
            d = decimals_of(cand["mean_str"])
            rec["k"], rec["n_eff"], rec["decimals"] = k or 1, n_eff, d
            if not grim_has_power(n_eff, d):
                rec.update(status="skipped", reason=f"no GRIM power (n_eff={n_eff} >= 10^{d})")
                out.append(rec)
                continue
            res = grim_consistent(cand["mean_str"], n_eff)
            rec.update(status="consistent" if res else "inconsistent", reason=why)
            out.append(rec)
    return out


def applicable(paper: dict) -> bool:
    return any(r["status"] in ("consistent", "inconsistent") for r in analyze(paper))


def check(paper: dict) -> list[dict]:
    flags = []
    for r in analyze(paper):
        if r["status"] != "inconsistent":
            continue
        flags.append({
            "location": {"kind": "cell", "table_id": r["table_id"], "row": r["row"], "col": r["col"]},
            "message": (f"GRIM-inconsistent mean {r['mean']} with n={r['n']}"
                        f"{' x k=' + str(r['k']) + ' items' if r.get('k', 1) > 1 else ''}"
                        f" (n from {r['n_source']}): no integer total T gives round(T/{r['n_eff']}, "
                        f"{r['decimals']}) = {r['mean']}"),
        })
    return flags
