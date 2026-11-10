"""Python port of statcheck (Nuijten et al. 2016; Epskamp & Nuijten, R package `statcheck`).

FALLBACK ONLY: the reference baseline is baselines/statcheck_check.py, which calls the
real R package.  This port exists because R cannot be installed right now.  It is
written from the published method description, not from the R source, so small
behavioural differences are expected (listed under "Known deviations" below).

Method
------
1. Extraction of APA-style results from paragraph text:
     t(df) = x, p ...          F(df1, df2) = x, p ...      r(df) = x, p ...
     chi2(df [, N = n]) = x, p ...  (χ2, χ², chi2, chi-square(d), X2)
     z = x, p ...              Q(df) = x, p ...  (also Qw / Qb)
   The test statistic may be compared with =, < or >; the p-value with =, <, >
   (≤ / ≥ accepted as non-strict) or be reported as "ns".  Unicode minus / dashes,
   unicode spaces, thousands commas in the statistic, and leading-zero-less decimals
   (.05) are handled.  Statistic and p must be separated by a comma or semicolon.
2. Recompute the two-tailed p (scipy.stats): t; F; r via t = r*sqrt(df/(1-r^2));
   chi2; z; Q via the chi2 distribution.
3. Rounding: a statistic reported with d decimals stands for the interval
   [x - 0.5*10^-d, x + 0.5*10^-d]; this gives an interval of computable p-values
   [lowP, upP].  A result is an ERROR iff no p in that interval is consistent with
   the reported p:
     p = v   (v has k decimals): [lowP, upP] must overlap [v - .5*10^-k, v + .5*10^-k]
     p < v : lowP < v            p ≤ v : lowP <= v
     p > v : upP > v             p ≥ v : upP >= v
     ns    : treated as p > .05
   p = 0 / p = .000 is always an error (statcheck's pZeroError = TRUE default).
4. One-tailed: no correction by default.  If the paper text mentions "one-tailed",
   "one-sided" or "directional" (statcheck's OneTailedTxt search), a t / r / z result
   that is an error two-tailed but consistent with half the p-value is NOT an error.
5. Decision error (alpha = .05, p = .05 counts as significant): an ERROR where
   reported significance differs from the significance of the p computed from the
   reported statistic:  "= v": (v <= .05) != (p <= .05);  "< v" with v <= .05 and
   p > .05;  "> v" / ns with v >= .05 and p <= .05.

Interface (docs/SCHEMA.md §3): extract(text) -> list[dict]; check(paper) -> flags;
applicable(paper) -> bool (>= 1 result extracted).  Flag location is the matched span
{"kind": "text", "para_id", "start", "end"} in the original paragraph text.

Known deviations from the R package
-----------------------------------
- Regexes are re-implemented, so the set of extracted results can differ at the
  margins (e.g. we require a non-letter before the test letter; we accept ";" as
  separator and ≤/≥ comparisons; we do not parse "p = 1.2 × 10^-5").
- Exact-p rounding uses interval overlap (accepts round-half-up and half-even),
  slightly more lenient than R's round() on exact ties.
- One-tailed text search runs over the whole paper (paragraphs) and covers t, r, z.
- Only paragraph text is scanned (tables are not), matching the R wrapper here.
"""
from __future__ import annotations

import json
import math
import re
from decimal import Decimal, InvalidOperation

from scipy import stats

ALPHA = 0.05

# ---------------------------------------------------------------- normalisation
# Same-length, char-for-char normalisation so regex offsets index the original text.
_CHAR_MAP = {}
for _c in "−–—‐‑‒﹣－":
    _CHAR_MAP[_c] = "-"
for _c in "              　\t\n\r":
    _CHAR_MAP[_c] = " "
_CHAR_MAP.update({"²": "2", "＝": "=", "＜": "<", "＞": ">", "Χ": "χ",
                  "≤": "≤", "≥": "≥", " ": " ", "，": ","})


def _norm(text: str) -> str:
    return "".join(_CHAR_MAP.get(c, c) for c in text)


# ---------------------------------------------------------------- regexes
_S = r" ?"                                   # optional single space (statcheck's \s?)
_NUM = r"\d*(?:,\d{3})*\.?\d+"              # 1,234.5 / .30 / 2
_DF = r"\d*\.?\d+"                          # df may be fractional (Welch)
_DF_T = r"\d{1,3}(?:,\d{3})+|\d*\.?\d+"     # t/r df may have thousands commas
_CMP = r"[<>=≤≥]"
_STAT = rf"(?P<scmp>{_CMP}){_S}(?P<sign>-)?{_S}(?P<stat>{_NUM})"
_P = (rf"(?:(?P<ns>(?<![A-Za-z])ns\b)|"
      rf"(?<![A-Za-z])p{_S}(?P<pcmp>{_CMP}){_S}(?P<p>\d?\.\d+(?:[eE]-?\d+)?))")
_SEP = rf"{_S}[,;]{_S}"
_LB = r"(?<![A-Za-z0-9_])"

_PATTERNS = [
    ("t", re.compile(rf"{_LB}t{_S}\({_S}(?P<df1>{_DF_T}){_S}\){_S}{_STAT}{_SEP}{_P}")),
    ("F", re.compile(rf"{_LB}F{_S}\({_S}(?P<df1>{_DF}){_S},{_S}(?P<df2>{_DF}){_S}\){_S}{_STAT}{_SEP}{_P}")),
    ("r", re.compile(rf"{_LB}r{_S}\({_S}(?P<df1>{_DF_T}){_S}\){_S}{_STAT}{_SEP}{_P}")),
    ("chi2", re.compile(
        rf"{_LB}(?:χ{_S}2?|[Cc]hi{_S}2|[Cc]hi-squared?|X2|x2){_S}"
        rf"\({_S}(?P<df1>{_DF}){_S}(?:,{_S}[Nn]{_S}={_S}\d{{1,3}}(?:,\d{{3}})*\d*{_S})?\){_S}{_STAT}{_SEP}{_P}")),
    ("z", re.compile(rf"{_LB}[zZ]{_S}{_STAT}{_SEP}{_P}")),
    ("Q", re.compile(rf"{_LB}Q(?:w|b|within|between)?{_S}\({_S}(?P<df1>{_DF}){_S}\){_S}{_STAT}{_SEP}{_P}")),
]

_ONE_TAILED = re.compile(r"one[- ]?tailed|one[- ]?sided|directional", re.I)


def _num(s: str) -> float:
    return float(s.replace(",", ""))


def _decimals(s: str) -> int:
    """Number of decimals implied by a number string (handles e-notation)."""
    try:
        e = Decimal(s.replace(",", "")).as_tuple().exponent
    except InvalidOperation:
        return 0
    return max(0, -int(e))


# ---------------------------------------------------------------- p computation
def compute_p(test: str, stat: float, df1: float | None, df2: float | None) -> float:
    """Two-tailed p for the given test statistic (F / chi2 / Q are upper-tail)."""
    if test == "t":
        return float(2 * stats.t.sf(abs(stat), df1))
    if test == "z":
        return float(2 * stats.norm.sf(abs(stat)))
    if test == "r":
        r = abs(stat)
        if r >= 1:
            return 0.0
        t = r * math.sqrt(df1 / (1 - r * r))
        return float(2 * stats.t.sf(t, df1))
    if test == "F":
        return float(stats.f.sf(max(stat, 0.0), df1, df2))
    if test in ("chi2", "Q"):
        return float(stats.chi2.sf(max(stat, 0.0), df1))
    raise ValueError(test)


def p_bounds(test, stat, stat_str, scmp, df1, df2):
    """[lowP, upP] of p-values compatible with the reported (rounded) statistic."""
    a = abs(stat)
    half = 0.5 * 10 ** -_decimals(stat_str)
    lo_s, hi_s = max(0.0, a - half), a + half
    if test == "r":
        hi_s = min(hi_s, 1.0)
    # '<' on the statistic: true |stat| in [0, hi]; '>' : [lo, inf)
    if scmp in ("<", "≤"):
        lo_s = 0.0
    elif scmp in (">", "≥"):
        hi_s = math.inf
    low_p = 0.0 if hi_s == math.inf else compute_p(test, hi_s, df1, df2)
    up_p = compute_p(test, lo_s, df1, df2)
    return low_p, up_p


def _consistent(low_p, up_p, pcmp, rep, p_dec):
    if pcmp == "=":
        h = 0.5 * 10 ** -p_dec
        tol = 1e-12
        return low_p <= rep + h + tol and up_p >= rep - h - tol
    if pcmp == "<":
        return low_p < rep
    if pcmp == "≤":
        return low_p <= rep
    if pcmp == ">":
        return up_p > rep
    if pcmp == "≥":
        return up_p >= rep
    raise ValueError(pcmp)


def _decision_error(pcmp, rep, computed):
    if pcmp == "=":
        return (rep <= ALPHA) != (computed <= ALPHA)
    if pcmp in ("<", "≤"):
        return rep <= ALPHA and computed > ALPHA
    if pcmp in (">", "≥"):
        return rep >= ALPHA and computed <= ALPHA
    return False


# ---------------------------------------------------------------- extraction
def extract(text: str, one_tailed_txt: bool | None = None) -> list[dict]:
    """Extract and evaluate APA results from text.

    one_tailed_txt: whether one-tailed wording is present in the surrounding paper;
    None = search `text` itself.
    """
    if not text:
        return []
    norm = _norm(text)
    if one_tailed_txt is None:
        one_tailed_txt = bool(_ONE_TAILED.search(norm))
    out = []
    taken: list[tuple[int, int]] = []
    for test, rx in _PATTERNS:
        for m in rx.finditer(norm):
            s, e = m.span()
            if any(s < te and ts < e for ts, te in taken):
                continue
            g = m.groupdict()
            df1 = _num(g["df1"]) if g.get("df1") else None
            df2 = _num(g["df2"]) if g.get("df2") else None
            if (df1 is not None and df1 <= 0) or (df2 is not None and df2 <= 0):
                continue
            stat = _num(g["stat"]) * (-1 if g.get("sign") else 1)
            if test == "r" and abs(stat) > 1:
                continue  # statcheck drops impossible correlations
            if g.get("ns"):
                pcmp, rep, p_dec = ">", ALPHA, 2
            else:
                pcmp, rep, p_dec = g["pcmp"], float(g["p"]), _decimals(g["p"])
            if rep > 1:
                continue
            taken.append((s, e))
            computed = compute_p(test, stat, df1, df2)
            low_p, up_p = p_bounds(test, stat, g["stat"], g["scmp"], df1, df2)
            p_zero = pcmp == "=" and rep == 0
            consistent = (not p_zero) and _consistent(low_p, up_p, pcmp, rep, p_dec)
            one_tail = False
            if not consistent and not p_zero and one_tailed_txt and test in ("t", "r", "z"):
                if _consistent(low_p / 2, up_p / 2, pcmp, rep, p_dec):
                    consistent, one_tail = True, True
            error = not consistent
            dec_err = error and not p_zero and _decision_error(pcmp, rep, computed)
            if p_zero:
                dec_err = computed > ALPHA
            out.append({
                "test": test, "df1": df1, "df2": df2,
                "stat_comparison": g["scmp"], "stat": stat,
                "p_comparison": "ns" if g.get("ns") else pcmp,
                "reported_p": rep, "p_decimals": p_dec,
                "computed_p": computed, "low_p": low_p, "up_p": up_p,
                "error": error, "decision_error": bool(dec_err),
                "one_tailed_accepted": one_tail,
                "raw": text[s:e], "start": s, "end": e,
            })
    out.sort(key=lambda r: r["start"])
    return out


# ---------------------------------------------------------------- paper interface
def _paper_one_tailed(paper: dict) -> bool:
    return any(_ONE_TAILED.search(_norm(p.get("text") or "")) for p in paper.get("paragraphs") or [])


def analyze(paper: dict) -> list[dict]:
    ot = _paper_one_tailed(paper)
    res = []
    for para in paper.get("paragraphs") or []:
        for r in extract(para.get("text") or "", one_tailed_txt=ot):
            r["para_id"] = para["para_id"]
            res.append(r)
    return res


def _fmt_p(p: float) -> str:
    return f"{p:.4f}" if p >= 1e-4 else f"{p:.2e}"


def flags_from(results: list[dict]) -> list[dict]:
    flags = []
    for r in results:
        if not r["error"]:
            continue
        de = " (decision error)" if r["decision_error"] else ""
        flags.append({
            "location": {"kind": "text", "para_id": r["para_id"], "start": r["start"], "end": r["end"]},
            "message": (f"statcheck_py: reported '{r['raw']}' but recomputed p = {_fmt_p(r['computed_p'])} "
                        f"(range from statistic rounding {_fmt_p(r['low_p'])}-{_fmt_p(r['up_p'])}){de}"),
        })
    return flags


def check(paper: dict) -> list[dict]:
    return flags_from(analyze(paper))


def applicable(paper: dict) -> bool:
    return len(analyze(paper)) > 0


# ---------------------------------------------------------------- corpus run
def main():
    import glob
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    files = sorted(glob.glob(os.path.join(root, "corpus", "parsed", "*.json")))
    out = {}
    n_res = n_err = n_de = 0
    for f in files:
        paper = json.load(open(f, encoding="utf-8"))
        if "__" in paper["paper_id"]:
            continue
        res = analyze(paper)
        fl = flags_from(res)
        out[paper["paper_id"]] = {"applicable": len(res) > 0, "n_extracted": len(res), "flags": fl}
        n_res += len(res)
        n_err += sum(r["error"] for r in res)
        n_de += sum(r["decision_error"] for r in res)
    with open(os.path.join(here, "results_statcheck_py.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    app = sum(v["applicable"] for v in out.values())
    print(f"papers applicable: {app} / {len(out)}")
    print(f"results extracted: {n_res}")
    print(f"errors: {n_err}  (papers with >=1 error: {sum(bool(v['flags']) for v in out.values())})")
    print(f"decision errors: {n_de}")


if __name__ == "__main__":
    main()
