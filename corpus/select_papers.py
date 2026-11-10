"""Select the final corpus from parsed candidates.

Usage: python3 -I corpus/select_papers.py [--target 40] [--seed 20261110]

A table is "results-like" if it is
  - REGRESSION: mostly numeric, and (>=2 starred numbers AND (>=2 parenthesised numbers
    OR an SE/coefficient header OR a regression-ish caption)) OR (>=3 parenthesised numbers
    AND regression-ish caption/header) OR (coefficient header AND SE/t/z/p header column);
  - DESCRIPTIVE: mostly numeric and caption/header/stub mention mean/SD/N/%/descriptive etc.
A paper is eligible if it has >=1 results-like table. Papers whose title/abstract-ish text
signals a formal model / simulation and that have no regression table are deprioritised.
Ranking: regression-table papers first; then score = 10*[has a classic stars+(SE) table]
+ min(in-text APA-style stats, 10) (useful for statcheck); ties broken by seeded random. Top --target are kept.
Writes corpus/selection.csv; moves non-kept JSON (and .txt) to corpus/parsed_dropped/.
"""
import argparse
import csv
import json
import random
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARSED, DROPPED, TEXT = HERE / "parsed", HERE / "parsed_dropped", HERE / "text"
TEXT_DROPPED = HERE / "text_dropped"

MINUS = "-−–"
NUM = re.compile(rf"^[{MINUS}+]?\s*(\d+(,\d{{3}})*)?\.?\d+\s*%?\s*[*†‡§a-d]*$")
STAR = re.compile(rf"^[{MINUS}+]?\s*\d*\.?\d+\s*(\*+|†+)$")
PAREN = re.compile(rf"^[\(\[]\s*[{MINUS}+]?\s*\d*\.?\d+\s*[\)\]]\s*[*†]*$")
CI = re.compile(rf"^[\(\[]\s*[{MINUS}]?\d*\.?\d+\s*[,;––-]\s*[{MINUS}]?\d*\.?\d+\s*[\)\]]$")
NPCT = re.compile(rf"^[{MINUS}]?\d*\.?\d+\s*%?\s*[\(\[]\s*[{MINUS}]?\d*\.?\d+\s*%?\s*[\)\]]\s*[*†]*$")
PM = re.compile(rf"^[{MINUS}]?\d*\.?\d+\s*±\s*\d*\.?\d+\s*[*†]*$")
REG_WORD = re.compile(r"regress|logit|probit|ols\b|fixed[- ]effect|odds ratio|hazard|marginal effect|gmm|"
                      r"multilevel|predict", re.I)
REG_CAP = re.compile(r"regress|logit|probit|ols\b|coefficient|estimat|model|multilevel|fixed[- ]effect|"
                     r"odds ratio|hazard|marginal effect|predictors? of|determinants", re.I)
COEF_HDR = re.compile(r"\b(coef\w*|b|β|beta|estimate|or|odds ratio|ame)\b", re.I)
SE_HDR = re.compile(r"\b(s\.?e\.?|std\.? ?err\w*|standard error|t|z|p|p[- ]?value|t[- ]?stat\w*)\b", re.I)
DESC_KW = re.compile(r"\bmean\b|\bmedian\b|\bs\.?d\.?\b|std\.? ?dev|standard deviation|\bn\b|%|percent|"
                     r"descriptive|summary statistic|frequenc|proportion|share", re.I)
FORMAL = re.compile(r"\b(simulat\w*|agent[- ]based|game[- ]theor\w*|formal model|equilibri\w*|"
                    r"evolutionary game|differential equation|a simple political model)\b", re.I)
APA = re.compile(r"\b(t|F|z|r|χ2|χ ?2|chi2|chi-square[d]?)\s*\(\s*\d+(\.\d+)?(\s*,\s*\d+(\.\d+)?)?\s*\)\s*[=<>]|"
                 r"\bp\s*[<=>]\s*0?\.\d+")


def classify(t):
    cells = [c for r in t["rows"] for c in r]
    nonempty = [c for c in cells if c.strip()]
    if len(nonempty) < 4:
        return None, "tiny"
    numeric = [c for c in nonempty if NUM.match(c) or PAREN.match(c) or CI.match(c) or NPCT.match(c) or PM.match(c)]
    frac = len(numeric) / len(nonempty)
    if frac < 0.3:
        return None, f"numfrac={frac:.2f}"
    n_star = sum(bool(STAR.match(c)) for c in nonempty)
    n_paren = sum(bool(PAREN.match(c)) for c in nonempty)
    hdr = " | ".join(c for r in t["header_rows"] for c in r)
    cap = t["caption"] + " " + t["label"]
    stub = " ".join(r[0] for r in t["rows"] if r)
    regcap = bool(REG_CAP.search(cap))
    coef_h, se_h = bool(COEF_HDR.search(hdr)), bool(SE_HDR.search(hdr))
    if re.search(r"correlat", cap, re.I) and not REG_WORD.search(cap):
        return "descriptive", f"correlation table stars={n_star} numfrac={frac:.2f}"
    if (n_star >= 2 and (n_paren >= 2 or coef_h or se_h or regcap)) or \
       (n_paren >= 3 and (regcap or coef_h)) or (coef_h and se_h):
        kind = "regression"
        return kind, f"classic={int(n_star >= 2 and n_paren >= 2)} stars={n_star} parens={n_paren} numfrac={frac:.2f}"
    if DESC_KW.search(cap + " " + hdr + " " + stub):
        return "descriptive", f"numfrac={frac:.2f}"
    return None, f"numeric but no reg/desc cue numfrac={frac:.2f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20261110)
    args = ap.parse_args()
    # re-include anything previously dropped so the script is idempotent
    for src, dst in ((DROPPED, PARSED), (TEXT_DROPPED, TEXT)):
        if src.exists():
            for f in src.iterdir():
                shutil.move(str(f), dst / f.name)

    rng = random.Random(args.seed)
    rows = []
    for f in sorted(PARSED.glob("*.json")):
        p = json.loads(f.read_text())
        cls = [classify(t) for t in p["tables"]]
        kinds = [k for k, _ in cls]
        n_classic = sum(k == "regression" and "classic=1" in why for k, why in cls)
        n_reg, n_desc = kinds.count("regression"), kinds.count("descriptive")
        text = " ".join(x["text"] for x in p["paragraphs"])
        n_apa = len(APA.findall(text))
        formal = bool(FORMAL.search(p["title"])) or len(FORMAL.findall(text)) >= 8
        rows.append(dict(paper_id=p["paper_id"], n_tables=len(p["tables"]), n_reg=n_reg, n_desc=n_desc,
                         n_apa=n_apa, retracted=p["title"].upper().startswith(("RETRACTED", "EXPRESSION OF CONCERN")), n_classic=n_classic, formal=formal, tie=rng.random()))

    def eligible(r):
        return (r["n_reg"] + r["n_desc"] > 0) and not r["retracted"] and not (r["formal"] and r["n_reg"] == 0)

    ranked = sorted([r for r in rows if eligible(r)],
                    key=lambda r: (r["n_reg"] == 0, r["formal"],
                                   -(10 * (r["n_classic"] > 0) + min(r["n_apa"], 10)), r["tie"]))
    keep = {r["paper_id"] for r in ranked[: args.target]}

    DROPPED.mkdir(exist_ok=True)
    TEXT_DROPPED.mkdir(exist_ok=True)
    with open(HERE / "selection.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["paper_id", "n_tables", "n_results_tables", "kept", "reason"])
        for r in rows:
            kept = r["paper_id"] in keep
            detail = f"reg={r['n_reg']} (stars+SE-parens={r['n_classic']}) desc={r['n_desc']} apa_stats={r['n_apa']}" + (" formal/sim" if r["formal"] else "")
            if kept:
                reason = "kept: " + ("regression table" if r["n_reg"] else "descriptive table") + f" ({detail})"
            elif r["retracted"]:
                reason = f"dropped: retracted / expression of concern ({detail})"
            elif r["n_tables"] == 0:
                reason = "dropped: no parseable tables"
            elif r["n_reg"] + r["n_desc"] == 0:
                reason = f"dropped: no results-like table ({detail})"
            elif r["formal"] and r["n_reg"] == 0:
                reason = f"dropped: formal model/simulation without regression table ({detail})"
            else:
                reason = f"dropped: eligible but below target cutoff ({detail})"
            w.writerow([r["paper_id"], r["n_tables"], r["n_reg"] + r["n_desc"], int(kept), reason])
            if not kept:
                shutil.move(str(PARSED / f"{r['paper_id']}.json"), DROPPED / f"{r['paper_id']}.json")
                t = TEXT / f"{r['paper_id']}.txt"
                if t.exists():
                    shutil.move(str(t), TEXT_DROPPED / t.name)
    kr = [r for r in rows if r["paper_id"] in keep]
    print(f"candidates={len(rows)} eligible={len(ranked)} kept={len(kr)} "
          f"tables={sum(r['n_tables'] for r in kr)} reg_tables={sum(r['n_reg'] for r in kr)} "
          f"desc_tables={sum(r['n_desc'] for r in kr)} papers_with_reg={sum(r['n_reg']>0 for r in kr)} "
          f"papers_with_classic={sum(r['n_classic']>0 for r in kr)} papers_desc_only={sum(r['n_reg']==0 for r in kr)}")


if __name__ == "__main__":
    main()
