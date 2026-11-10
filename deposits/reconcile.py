"""Reconcile numbers reported in a paper with numbers recomputed from its deposited data.

Run from repo root:
  uv run --with pandas --with openpyxl --with pyreadstat --with xlrd python -I deposits/reconcile.py

Inputs : corpus/parsed/<id>.json, deposits/availability.csv, deposits/data/<id>/*
Output : deposits/reconcile_results.json

Deterministic; no network. Two kinds of target, both deliberately simple:

1. Total sample size (one comparison per paper). Reported Ns = every integer >= 10 given as
   "N = x" / "n = x" / "sample of x" in captions, table footnotes, header cells or body text, or in a
   table row labelled N / Obs / Observations / Sample size ("34.449" in a count cell is read as
   34449). Recomputed = row count of each loaded data table (file or sheet, all-empty rows dropped).
   Match if some data table's row count equals a reported N, or does so after subtracting a single
   exclusion stated in the text as "<k> <participants|respondents|...> were excluded/removed/dropped".
   The record also gives the "headline" N (most frequent table-level N) for context. This is a
   consistency check ("the deposit's N is reported in the paper"), not proof it is the analysis N.

2. Descriptive-statistics tables. A table qualifies if its (last-level) column headers name
   Mean / SD / Min / Max / N / % (or a combined "M (SD)" / "M±SD" column), with each statistic in at
   most one column (tables with the same statistic for several groups are skipped: the subgroup
   cannot be reconstructed safely). A row is used only if exactly one of its label cells (cells left
   of the first statistic column) matches exactly one data column:
     tier 1 exact column name; tier 2 same name after lower-casing and dropping non-alphanumerics,
     against column names or SPSS/Stata variable labels; tier 3 difflib ratio >= 0.92 on the
     normalised strings (both >= 8 chars) with the runner-up < 0.85.
   Anything else -> abstain (logged with a reason).
   Panels: a row with text but no statistics starts a new panel; listwise samples are per panel.
   Sample: if the row reports N, the recomputed sample is either all non-missing values of the column
   or listwise-complete rows across all columns matched in that table (same data table) -- whichever
   reproduces N. If neither does, the N comparison is recorded and the row's other statistics are
   skipped (subsample unknown).
   Statistic tolerance: half a unit of the last reported decimal (rounding at reported precision).
   SD matches with ddof=1 or ddof=0. "%" columns are compared with 100*mean for 0/1 columns only;
   a "Mean (%)" column accepts mean or 100*mean for 0/1 columns.
"""
import csv
import difflib
import json
import os
import re
import warnings

import pandas as pd
import pyreadstat

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEP = os.path.join(ROOT, "deposits")
DATA = os.path.join(DEP, "data")

# ---------------------------------------------------------------- data loading

def _good_header(cols):
    cols = [str(c) for c in cols]
    bad = sum(c.startswith("Unnamed") or c.strip() == "" or c == "nan" for c in cols)
    return len(cols) >= 2 and bad / len(cols) <= 0.3


def load_tables(pid):
    """Return list of dicts: name, df, labels(col->label). Skips sheets that are not rectangular data."""
    d = os.path.join(DATA, pid)
    out, skipped = [], []
    if not os.path.isdir(d):
        return out, skipped
    for fn in sorted(os.listdir(d)):
        p = os.path.join(d, fn)
        low = fn.lower()
        try:
            frames = []
            if low.endswith((".csv", ".tsv", ".tab", ".txt")):
                sep = "\t" if low.endswith((".tsv", ".tab")) else None
                frames.append((fn, pd.read_csv(p, sep=sep, engine="python"), {}))
            elif low.endswith((".xlsx", ".xls")):
                for sh, df in pd.read_excel(p, sheet_name=None).items():
                    frames.append((f"{fn}[{sh}]", df, {}))
            elif low.endswith(".dta"):
                df, m = pyreadstat.read_dta(p)
                frames.append((fn, df, {k: v for k, v in m.column_names_to_labels.items() if v}))
            elif low.endswith(".sav"):
                df, m = pyreadstat.read_sav(p)
                frames.append((fn, df, {k: v for k, v in m.column_names_to_labels.items() if v}))
            for name, df, lab in frames:
                df = df.dropna(how="all")
                if df.shape[0] < 10 or not _good_header(df.columns):
                    skipped.append(f"{name}: not a rectangular data table (rows={df.shape[0]})")
                    continue
                df.columns = [str(c).strip() for c in df.columns]
                out.append(dict(name=name, df=df, labels=lab))
        except Exception as e:
            skipped.append(f"{fn}: load error {type(e).__name__}: {e}")
    return out, skipped

# ---------------------------------------------------------------- number parsing

NUM_RE = re.compile(r"[-+]?\d[\d,]*\.?\d*|[-+]?\.\d+")


def norm_minus(s):
    return s.replace("−", "-").replace("–", "-").replace("—", "-")


def parse_num(cell):
    """Single number in a cell (stars/daggers/% allowed). Returns (value, decimals) or None."""
    s = norm_minus(str(cell)).strip()
    s = re.sub(r"[*†‡§#a-zA-Z%\s]+$", "", s)  # trailing stars, footnote letters, %
    s = s.strip()
    m = re.fullmatch(r"\(?\s*([-+]?(?:\d{1,3}(?:,\d{3})+|\d+)?(?:\.\d+)?)\s*\)?", s)
    if not m or not m.group(1) or m.group(1) in "+-":
        return None
    t = m.group(1).replace(",", "")
    try:
        v = float(t)
    except ValueError:
        return None
    dec = len(t.split(".")[1]) if "." in t else 0
    return v, dec


def parse_pair(cell):
    """'3.45 ± 0.67' or '3.45 (0.67)' -> ((m,dec),(sd,dec))"""
    s = norm_minus(str(cell))
    m = re.fullmatch(r"\s*([-+]?\d[\d,]*\.?\d*)\s*(?:±|\+/-|\()\s*([-+]?\d[\d,]*\.?\d*)\s*\)?\s*", s)
    if not m:
        return None
    a, b = parse_num(m.group(1)), parse_num(m.group(2))
    return (a, b) if a and b else None


def parse_count(cell):
    """Integer count; '34.449' / '1.243' (dot as thousands separator) is read as 34449 / 1243,
    since a count cannot be fractional."""
    s = norm_minus(str(cell)).strip()
    s = re.sub(r"[*†‡§a-zA-Z\s]+$", "", s)
    if re.fullmatch(r"\d{1,3}(?:\.\d{3})+", s):
        return int(s.replace(".", ""))
    pv = parse_num(s)
    return int(pv[0]) if pv and pv[1] == 0 else None


def within(rep, dec, val):
    return abs(rep - val) <= 0.5 * 10 ** (-dec) + 1e-9 * max(1.0, abs(val))

# ---------------------------------------------------------------- reported N

N_PATS = [re.compile(r"(?<![A-Za-z])[Nn]\s*=\s*(\d{1,3}(?:,\d{3})+|\d+)(?![\d.])"),
          re.compile(r"\b(?:sample|total) of (\d{1,3}(?:,\d{3})+|\d+)(?![\d.])")]
N_LABELS = {"n", "obs", "obs.", "observations", "number of observations", "sample size", "no. of observations"}
EXCL_RE = re.compile(r"(\d{1,3}(?:,\d{3})+|\d+)\s+(?:participants|respondents|subjects|cases|observations|"
                     r"individuals|responses|questionnaires)\s+(?:were|was)\s+(?:excluded|removed|dropped|eliminated|discarded)")


def reported_ns(paper):
    found = []  # (value, where)
    def scan(text, where):
        for pat in N_PATS:
            for m in pat.finditer(norm_minus(text)):
                v = int(m.group(1).replace(",", ""))
                if v >= 10:
                    found.append((v, where))
    for p in paper["paragraphs"]:
        scan(p["text"], f"text {p['para_id']}")
    for t in paper["tables"]:
        tid = t["table_id"]
        scan(t["caption"], f"{tid} caption")
        scan(t.get("footnotes", ""), f"{tid} footnote")
        for hr in t["header_rows"]:
            for c in hr:
                scan(c, f"{tid} header")
        for i, r in enumerate(t["rows"]):
            if r and r[0].strip().lower() in N_LABELS:
                for c in r[1:]:
                    v = parse_count(c)
                    if v and v >= 10:
                        found.append((v, f"{tid} row {i} '{r[0]}'"))
    excl = sorted({int(m.group(1).replace(",", "")) for p in paper["paragraphs"]
                   for m in EXCL_RE.finditer(p["text"])})
    return found, excl

# ---------------------------------------------------------------- descriptive tables

STAT_NAMES = {
    "mean": {"mean", "means", "average", "average value", "m", "avg", "avg.", "mean value"},
    "mean_pct": {"mean (%)", "mean(%)", "mean/%", "mean / %"},
    "sd": {"sd", "s.d.", "s.d", "std", "std.", "std. dev.", "std. dev", "std dev", "std. deviation",
           "standard deviation", "the standard deviation", "st. dev.", "stdev", "std.dev."},
    "min": {"min", "min.", "minimum", "minimum value", "the minimum", "minimum values"},
    "max": {"max", "max.", "maximum", "maximum value", "the maximum", "maximum values"},
    "n": {"n", "obs", "obs.", "observations", "number of observations", "sample size", "no. of obs."},
    "msd": {"m (sd)", "m(sd)", "mean (sd)", "mean(sd)", "m±sd", "m ± sd", "mean ± sd", "mean±sd"},
    "pct": {"%", "(%)", "percent", "percentage"},
}


def stat_of(h):
    h = " ".join(norm_minus(h).lower().split())
    for k, names in STAT_NAMES.items():
        if h in names:
            return k
    if any(h.startswith(n) for n in STAT_NAMES["msd"]):  # e.g. "M±SD, n (%)"
        return "msd"
    return None


def col_headers(t):
    """Last-level header text per column; if no stat word in header rows, try first body rows."""
    ncol = len(t["rows"][0]) if t["rows"] else 0
    hdr = [""] * ncol
    for hr in t["header_rows"]:
        for j, c in enumerate(hr[:ncol]):
            if c.strip():
                hdr[j] = c.strip()
    body_start = 0
    if not any(stat_of(h) for h in hdr):
        for i, r in enumerate(t["rows"][:2]):
            if sum(bool(stat_of(c)) for c in r) >= 2:
                hdr, body_start = [c.strip() for c in r], i + 1
                break
    return hdr, body_start


def norm(s):
    return re.sub(r"[^0-9a-z\u0080-￿]", "", str(s).lower())


def match_label(label, tables):
    """Return (table_index, column, tier) or (None, reason)."""
    lab = label.strip()
    if not lab or parse_num(lab):
        return None, "empty/numeric label"
    hits = [(i, c) for i, tb in enumerate(tables) for c in tb["df"].columns if c == lab]
    tier = 1
    if not hits:
        nl = norm(lab)
        if len(nl) < 2:
            return None, "label too short"
        hits = [(i, c) for i, tb in enumerate(tables) for c in tb["df"].columns
                if norm(c) == nl or norm(tb["labels"].get(c, "")) == nl]
        tier = 2
    if not hits:
        nl = norm(lab)
        if len(nl) >= 8:
            sc = sorted(((difflib.SequenceMatcher(None, nl, norm(c)).ratio(), i, c)
                         for i, tb in enumerate(tables) for c in tb["df"].columns if len(norm(c)) >= 8),
                        reverse=True)
            if sc and sc[0][0] >= 0.92 and (len(sc) == 1 or sc[1][0] < 0.85):
                hits, tier = [(sc[0][1], sc[0][2])], 3
    if not hits:
        return None, "no data column matches label"
    if len(hits) > 1:
        # identical columns in several tables are not ambiguous
        first = tables[hits[0][0]]["df"][hits[0][1]]
        if all(tables[i]["df"][c].reset_index(drop=True).equals(first.reset_index(drop=True)) for i, c in hits[1:]):
            return (hits[0][0], hits[0][1], tier), None
        return None, f"ambiguous: {len(hits)} data columns match"
    return (hits[0][0], hits[0][1], tier), None


def numeric(series):
    s = series.dropna()
    s = s[~s.astype(str).str.strip().isin(["", "."])]
    x = pd.to_numeric(s, errors="coerce")
    if len(s) and x.isna().mean() > 0.05:
        return None
    return x.dropna()


WINSOR_SUFFIXES = ("_w", "_win", "_wins", "_winsor", "_winsorized")


def winsor_sibling(df, c):
    for suf in WINSOR_SUFFIXES:
        if c + suf in df.columns:
            return c + suf
    return None


def describe_tables(paper, tables):
    comps, abst = [], []
    for t in paper["tables"]:
        if not t["rows"]:
            continue
        hdr, b0 = col_headers(t)
        stats = {}
        for j, h in enumerate(hdr):
            k = stat_of(h)
            if k:
                stats.setdefault(k, []).append(j)
        useful = set(stats) & {"mean", "mean_pct", "sd", "min", "max", "msd", "pct"}
        if not useful:
            continue
        if any(len(v) > 1 for v in stats.values()):
            abst.append(dict(table=t["table_id"], reason="same statistic in several columns (grouped table)"))
            continue
        first_stat = min(j for v in stats.values() for j in v)
        rows = t["rows"][b0:]
        # pass 1: match rows; a row with text but no statistics starts a new panel segment
        matched, seg_of, seg = [], {}, 0
        for ri, r in enumerate(rows):
            if not any(r[j].strip() for v in stats.values() for j in v):
                if any(c.strip() for c in r):
                    seg += 1
                continue
            seg_of[ri + b0] = seg
            cands = []
            reasons = []
            for j in range(first_stat):
                if r[j].strip():
                    m, why = match_label(r[j], tables)
                    if m:
                        cands.append((j, m))
                    else:
                        reasons.append(f"'{r[j][:40]}': {why}")
            if len({c[1][:2] for c in cands}) == 1:
                matched.append((ri + b0, r, cands[0][1], r[cands[0][0]]))
            elif len(cands) > 1:
                abst.append(dict(table=t["table_id"], row=ri + b0, reason="label cells match different columns"))
            else:
                abst.append(dict(table=t["table_id"], row=ri + b0, reason="; ".join(reasons)[:160]))
        # listwise sample per (panel segment, data table) across all matched columns
        bytab = {}
        for ri, _, (ti, c, _), _ in matched:
            bytab.setdefault((seg_of[ri], ti), []).append(c)
        listwise = {}
        for key, cols in bytab.items():
            num = {c: pd.to_numeric(tables[key[1]]["df"][c], errors="coerce") for c in set(cols)}
            listwise[key] = pd.DataFrame(num).dropna().index
        # pass 2: compare
        for ri, r, (ti, c, tier), matched_cell in matched:
            tb = tables[ti]
            x_all = numeric(tb["df"][c])
            base = dict(table=t["table_id"], row=ri, row_label=matched_cell,
                        data_file=tb["name"], data_column=c, match_tier=tier)
            if x_all is None:
                abst.append(dict(table=t["table_id"], row=ri, reason=f"data column {c} not numeric"))
                continue
            x = x_all
            if "n" in stats:
                rep = parse_count(r[stats["n"][0]])
                if rep is not None:
                    x_lw = pd.to_numeric(tb["df"].loc[listwise[(seg_of[ri], ti)], c], errors="coerce").dropna()
                    if rep == len(x_all):
                        sample = "all non-missing"
                    elif rep == len(x_lw):
                        x, sample = x_lw, "listwise across matched table columns"
                    else:
                        sample = None
                    comps.append(dict(base, stat="n", col=stats["n"][0], reported=r[stats["n"][0]],
                                      recomputed=len(x_all), alt_recomputed=len(x_lw),
                                      result="match" if sample else "mismatch"))
                    if not sample:
                        abst.append(dict(table=t["table_id"], row=ri,
                                         reason="row N not reproducible; other statistics skipped"))
                        continue
            binary = len(x) > 0 and set(x.unique()) <= {0, 1}
            vals = {"mean": [x.mean()], "sd": [x.std(ddof=1), x.std(ddof=0)],
                    "min": [x.min()], "max": [x.max()]}
            todo = []
            for k in ("mean", "sd", "min", "max"):
                if k in stats:
                    todo.append((k, stats[k][0], r[stats[k][0]], vals[k]))
            if "mean_pct" in stats:
                j = stats["mean_pct"][0]
                todo.append(("mean", j, r[j], vals["mean"] + ([100 * x.mean()] if binary else [])))
            if "pct" in stats and binary:
                j = stats["pct"][0]
                todo.append(("pct", j, r[j], [100 * x.mean()]))
            if "msd" in stats:
                j = stats["msd"][0]
                pp = parse_pair(r[j])
                if pp:
                    for (v, d), k in zip(pp, ("mean", "sd")):
                        res = "match" if any(within(v, d, w) for w in vals[k]) else "mismatch"
                        comps.append(dict(base, stat=k, col=j, reported=r[j], recomputed=round(vals[k][0], 6),
                                          result=res))
            sib = winsor_sibling(tb["df"], c)
            xs = numeric(tb["df"][sib]) if sib else None
            if xs is not None and len(xs) != len(x):
                xs = None
            for k, j, cell, cands in todo:
                pv = parse_num(cell)
                if not pv:
                    continue
                v, d = pv
                res = "match" if any(within(v, d, w) for w in cands) else "mismatch"
                comp = dict(base, stat=k, col=j, reported=cell, recomputed=round(float(cands[0]), 6), result=res)
                if res == "mismatch" and xs is not None and k in ("mean", "sd", "min", "max"):
                    alt = {"mean": [xs.mean()], "sd": [xs.std(ddof=1), xs.std(ddof=0)],
                           "min": [xs.min()], "max": [xs.max()]}[k]
                    # annotation only: the result stays "mismatch" (exact-name rule); reviewed by hand
                    comp["winsorized_sibling"] = dict(column=sib, value=round(float(alt[0]), 6),
                                                      reproduces=any(within(v, d, w) for w in alt))
                comps.append(comp)
    return comps, abst


def total_n(paper, tables):
    """Is the deposit's row count reported as an N in the paper?
    Reported Ns = every "N = x" / "n = x" / "sample of x" in captions, footnotes, headers and body text,
    plus N/Observations rows. Match if some data table's row count (or row count minus one stated
    exclusion) equals one of them. For context the record also gives the "headline" N = most frequent
    table-level N (ties -> larger; text only if tables give none); a mismatch shows that value."""
    found, excl = reported_ns(paper)
    tab = [(v, w) for v, w in found if not w.startswith("text")]
    pool = tab or found
    if not pool or not tables:
        return None
    freq = {}
    for v, _ in pool:
        freq[v] = freq.get(v, 0) + 1
    rep = max(freq, key=lambda v: (freq[v], v))
    where = sorted({w for v, w in pool if v == rep})
    counts = {tb["name"]: len(tb["df"]) for tb in tables}
    allv = {v for v, _ in found}
    ok = [(n, None) for n in counts.values() if n in allv] + \
         [(n, k) for n in counts.values() for k in excl if n - k in allv]
    elsewhere = sorted({w for v, w in found for n in counts.values() if v == n})
    return dict(table="; ".join(where)[:200], stat="total_n", headline_n_equals_rows=rep in counts.values(),
                reported=rep, reported_source="tables" if tab else "text",
                recomputed=max(counts.values()), data_row_counts=counts, stated_exclusions=excl,
                all_reported_ns=sorted({v for v, _ in found}),
                result="match" if ok else "mismatch",
                note=(f"after stated exclusion of {ok[0][1]}" if ok and ok[0][1] else ""),
                row_count_reported_elsewhere=elsewhere[:10])


def main():
    av = list(csv.DictReader(open(os.path.join(DEP, "availability.csv"))))
    out = {}
    for row in av:
        pid = row["paper_id"]
        paper = json.load(open(os.path.join(ROOT, "corpus", "parsed", pid + ".json")))
        tables, skipped = load_tables(pid) if row["category"] == "fetchable" else ([], [])
        rec = dict(category=row["category"], data_found=bool(tables),
                   data_tables=[f"{tb['name']} ({tb['df'].shape[0]}x{tb['df'].shape[1]})" for tb in tables],
                   skipped_files=skipped)
        comps, abst = [], []
        if tables:
            tn = total_n(paper, tables)
            if tn:
                comps.append(tn)
            else:
                abst.append(dict(reason="no reported N found"))
            c2, a2 = describe_tables(paper, tables)
            comps += c2
            abst += a2
        rec.update(n_comparisons=len(comps), n_match=sum(c["result"] == "match" for c in comps),
                   n_mismatch=sum(c["result"] == "mismatch" for c in comps),
                   mismatches=[c for c in comps if c["result"] == "mismatch"],
                   matches=[c for c in comps if c["result"] == "match"],
                   abstentions=abst)
        out[pid] = rec
    tot = dict(papers=len(out), data_found=sum(r["data_found"] for r in out.values()),
               papers_with_comparisons=sum(r["n_comparisons"] > 0 for r in out.values()),
               comparisons=sum(r["n_comparisons"] for r in out.values()),
               match=sum(r["n_match"] for r in out.values()),
               mismatch=sum(r["n_mismatch"] for r in out.values()))
    json.dump(dict(summary=tot, papers=out), open(os.path.join(DEP, "reconcile_results.json"), "w"),
              indent=1, ensure_ascii=False, default=str)
    print(json.dumps(tot))
    for pid, r in out.items():
        if r["data_found"]:
            print(pid[-7:], r["n_comparisons"], r["n_match"], r["n_mismatch"], len(r["abstentions"]))


if __name__ == "__main__":
    main()
