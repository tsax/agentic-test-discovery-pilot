#!/usr/bin/env python3
"""Frozen error planter (answer key). Stdlib only.

Usage:
    python3 -I planted/plant.py [--corpus corpus/parsed] [--out planted]

Reads parsed papers (docs/SCHEMA.md section 1), writes planted copies with exactly one
error each plus manifest.jsonl per split (SCHEMA section 2) and stats.json.
Error-type catalog: planted/ERROR_TYPES.md. Split: planted/split.json.

Conventions: cell locations use 0-based indices into `rows` (header rows excluded) and
0-based column indices. Text locations are character offsets into the NEW paragraph text,
spanning the whole inconsistent expression (number, or estimate+CI, or statistic+p).
"""
import argparse
import copy
import json
import math
import os
import random
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MAX_COPIES = 3

# --------------------------------------------------------------------------- numbers

NUMU = r'(?:\d+(?:\.\d+)?|\.\d+)'                    # unsigned
NUMS = r'[−\-–]?(?:\d+(?:\.\d+)?|\.\d+)'             # signed
INT = r'(?:\d{1,3}(?:,\d{3})+|\d+)'
MARKS = r'[*†‡+§#]*'
TOK_RE = re.compile(
    r'(?<![\w.,\-−–])(?P<sign>[−\-–])?(?P<core>\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+\.\d+|\.\d+|\d+)(?![\w]|[.,]\d)')
# Like TOK_RE but sign attaches even after letters/brackets (for scanning text numbers).
TXT_TOK_RE = re.compile(
    r'(?<![\d.,])(?P<sign>[−\-–](?=[\d.]))?(?P<core>\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+\.\d+|\.\d+|\d+)(?![\w]|[.,]\d)')


def to_float(s):
    s = s.strip().replace('−', '-').replace('–', '-').replace(',', '')
    return float(s)


def decimals(s):
    s = s.strip()
    return len(s.split('.', 1)[1]) if '.' in s else 0


def half_unit(s):
    return 0.5 * 10 ** (-decimals(s))


def _split_sign(s):
    if s and s[0] in '−-–':
        return s[0], s[1:]
    return '', s


def _commafy(intpart):
    out = ''
    while len(intpart) > 3:
        out = ',' + intpart[-3:] + out
        intpart = intpart[:-3]
    return intpart + out


def _well_formed(core, leading_dot):
    """Normalise an unsigned core with no commas; reject malformed ones."""
    if core.count('.') > 1 or core in ('', '.'):
        return None
    if '.' in core:
        ip, fp = core.split('.')
        ip = ip.lstrip('0') or '0'
        if leading_dot and ip == '0':
            ip = ''
        if fp == '':
            return ip
        return ip + '.' + fp
    if len(core) > 1 and core[0] == '0':
        return None
    return core


def num_variants(s, kinds=('transpose', 'digit', 'shift')):
    """Realistic typo variants of a printed number. Returns list of (kind, new_string)."""
    sign, core = _split_sign(s)
    commas = ',' in core
    core = core.replace(',', '')
    leading_dot = core.startswith('.')
    out = []
    seen = {s}

    def add(kind, newcore):
        nc = _well_formed(newcore, leading_dot)
        if nc is None:
            return
        if commas and '.' not in nc:
            nc = _commafy(nc)
        elif commas:
            ip, fp = nc.split('.')
            nc = _commafy(ip) + '.' + fp
        ns = sign + nc
        if ns not in seen:
            seen.add(ns)
            out.append((kind, ns))

    intlen = core.index('.') if '.' in core else len(core)
    if 'transpose' in kinds:
        for i in range(len(core) - 1):
            a, b = core[i], core[i + 1]
            if a.isdigit() and b.isdigit() and a != b:
                if i == 0 and b == '0' and intlen > 1:
                    continue
                if i == 0 and a == '0' and intlen == 1:
                    continue  # don't move the leading "0." zero
                add('transpose', core[:i] + b + a + core[i + 2:])
    if 'digit' in kinds:
        for i, ch in enumerate(core):
            if not ch.isdigit():
                continue
            if i == 0 and ch == '0' and intlen == 1 and len(core) > 1:
                continue
            for d in (int(ch) - 1, int(ch) + 1):
                if 0 <= d <= 9:
                    if i == 0 and d == 0 and intlen > 1:
                        continue
                    add('digit', core[:i] + str(d) + core[i + 1:])
    if 'shift' in kinds:
        digits = core.replace('.', '')
        p = intlen
        for np_ in (p + 1, p - 1):
            if np_ < 0:
                continue
            if np_ > len(digits):
                nc = digits + '0' * (np_ - len(digits))
            else:
                nc = digits[:np_] + '.' + digits[np_:]
                if nc.endswith('.'):
                    nc = nc[:-1]
            if nc.startswith('.'):
                nc = '0' + nc
            add('shift_right' if np_ > p else 'shift_left', nc)
    return out


# --------------------------------------------------------------------------- distributions

def norm_p2(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def _betacf(a, b, x):
    fpmin = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > fpmin else fpmin)
    h = d
    for m in range(1, 500):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > fpmin else fpmin)
        c = 1.0 + aa / c
        c = c if abs(c) > fpmin else fpmin
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > fpmin else fpmin)
        c = 1.0 + aa / c
        c = c if abs(c) > fpmin else fpmin
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-14:
            break
    return h


def betai(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbt = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x)
    bt = math.exp(lbt)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def gammq(a, x):
    """Upper regularised incomplete gamma Q(a, x)."""
    if x <= 0:
        return 1.0
    lpre = -x + a * math.log(x) - math.lgamma(a)
    if x < a + 1.0:
        ap, s = a, 1.0 / a
        d = s
        for _ in range(2000):
            ap += 1.0
            d *= x / ap
            s += d
            if abs(d) < abs(s) * 1e-15:
                break
        return max(0.0, 1.0 - s * math.exp(lpre))
    fpmin = 1e-300
    b = x + 1.0 - a
    c = 1.0 / fpmin
    d = 1.0 / b
    h = d
    for i in range(1, 2000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = 1.0 / (d if abs(d) > fpmin else fpmin)
        c = b + an / c
        c = c if abs(c) > fpmin else fpmin
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-15:
            break
    return math.exp(lpre) * h


def t_p2(t, df):
    return betai(df / 2.0, 0.5, df / (df + t * t))


def f_p(f, d1, d2):
    if f <= 0:
        return 1.0
    return betai(d2 / 2.0, d1 / 2.0, d2 / (d2 + d1 * f))


def chi2_p(x, df):
    return gammq(df / 2.0, x / 2.0)


def stat_p(kind, val, df1, df2):
    v = abs(val)
    if kind == 't':
        return t_p2(v, df1)
    if kind == 'F':
        return f_p(v, df1, df2)
    if kind == 'chi2':
        return chi2_p(v, df1)
    if kind == 'z':
        return norm_p2(v)
    if kind == 'r':
        if v >= 1:
            return 0.0
        return t_p2(v * math.sqrt(df1 / (1 - v * v)), df1)
    raise ValueError(kind)


# --------------------------------------------------------------------------- helpers

def table_text(t):
    parts = [t.get('caption') or '', t.get('footnotes') or '', t.get('label') or '']
    for hr in t.get('header_rows') or []:
        parts.extend(hr)
    return ' '.join(parts)


def col_header(t, c):
    return ' '.join(hr[c] for hr in (t.get('header_rows') or []) if c < len(hr))


def all_table_cores(paper):
    cores = set()
    for t in paper.get('tables', []):
        texts = [t.get('caption') or '', t.get('footnotes') or '']
        for hr in t.get('header_rows') or []:
            texts.extend(hr)
        for row in t.get('rows', []):
            texts.extend(row)
        for s in texts:
            for m in TOK_RE.finditer(s):
                cores.add(m.group('core').replace(',', ''))
    return cores


def set_cell(paper, tid, r, c, new):
    for t in paper['tables']:
        if t['table_id'] == tid:
            old = t['rows'][r][c]
            t['rows'][r][c] = new
            return old
    raise KeyError(tid)


def get_para(paper, pid):
    for p in paper['paragraphs']:
        if p['para_id'] == pid:
            return p
    raise KeyError(pid)


def splice(s, start, end, new):
    return s[:start] + new + s[end:]


def cell_record(tid, r, c, before, after):
    return {"kind": "cell", "table_id": tid, "row": r, "col": c}, before, after


def text_record(paper, pid, start, end, new_sub, span_start, span_end_old):
    """Replace text[start:end] by new_sub; manifest span is [span_start, span_end_old) in OLD
    coordinates, recomputed for the NEW text."""
    p = get_para(paper, pid)
    old_text = p['text']
    before = old_text[span_start:span_end_old]
    new_text = splice(old_text, start, end, new_sub)
    delta = len(new_sub) - (end - start)
    ns, ne = span_start, span_end_old + delta
    p['text'] = new_text
    return {"kind": "text", "para_id": pid, "start": ns, "end": ne}, before, new_text[ns:ne]


# --------------------------------------------------------------------------- star legend / coef-SE pairs

LEG_SYM = re.compile(
    r'(\*{1,4}|†{1,2}|‡|\+|§|#)\s*[:=]?\s*\(?\s*[pP]\s*(?:-?\s*values?\s*)?([<≤])\s*(0?\.\d+)')
LEG_PCT = re.compile(r'(\*{1,3}|†)\s*[:=]?\s*(?:significant\s+at\s+(?:the\s+)?)?(\d+(?:\.\d+)?)\s*%')


STARG = re.compile(r'\*(?:\s?\*){0,2}')
LIST_SEP = re.compile(r'^[\s,;/()]*(?:and|or)?[\s,;/()]*$')
NUMLIST = re.compile(r'(?<![\w.])(\d*\.?\d+)\s*(%|percent)?')


def _list_legend(text):
    """'***, ** and * denote significance at the 1%, 5% and 10% levels' (any order, any count)."""
    groups = list(STARG.finditer(text))
    i = 0
    while i < len(groups):
        seq = [groups[i]]
        j = i + 1
        while j < len(groups):
            gap = text[seq[-1].end():groups[j].start()]
            if len(gap) > 8 or not LIST_SEP.match(gap):
                break
            seq.append(groups[j])
            j += 1
        if len(seq) >= 2:
            syms = [g.group().replace(' ', '') for g in seq]
            tail = text[seq[-1].end():seq[-1].end() + 140]
            if len(set(syms)) == len(syms) and re.search(r'signific|level|denote|indicate|represent', tail, re.I):
                nums = []
                for m in NUMLIST.finditer(tail):
                    v = float(m.group(1))
                    nums.append(v / 100 if (m.group(2) or v >= 1) else v)
                    if len(nums) == len(syms):
                        return list(zip(syms, nums))
        i = j
    return None


def parse_legend(text):
    """Return list [(symbol, threshold)] sorted by threshold descending, or None."""
    found = {}

    def put(sym, thr):
        if not (0 < thr <= 0.1):
            return False
        if sym in found and abs(found[sym] - thr) > 1e-12:
            found[sym] = None
            return False
        found[sym] = thr
        return True

    for m in LEG_SYM.finditer(text):
        put(m.group(1), float(m.group(3)))
    if not found:
        for sym, v in _list_legend(text) or []:
            put(sym, v)
    if not found:
        for m in LEG_PCT.finditer(text):
            v = float(m.group(2))
            if v <= 10:
                put(m.group(1), v / 100)
    if not found or any(v is None for v in found.values()):
        return None
    leg = sorted(found.items(), key=lambda kv: -kv[1])
    thrs = [v for _, v in leg]
    if len(set(thrs)) != len(thrs):
        return None
    stars = [(s, v) for s, v in leg if set(s) == {'*'}]
    for (s1, v1), (s2, v2) in zip(stars, stars[1:]):
        if len(s2) <= len(s1):
            return None
    return leg


def level_of(marks, legend):
    if marks == '':
        return 0
    for i, (sym, _) in enumerate(legend):
        if sym == marks:
            return i + 1
    return None


def level_for_p(p, legend):
    lvl = 0
    for i, (_, thr) in enumerate(legend):
        if p < thr:
            lvl = i + 1
    return lvl


def robust_levels(b_str, se_str, legend):
    """Set of star levels compatible with printed b and se (rounding; normal and t(40))."""
    b = abs(to_float(b_str))
    se = to_float(se_str)
    hb, hs = half_unit(b_str), half_unit(se_str)
    if se - hs <= 0:
        return None
    zlo = max(0.0, (b - hb) / (se + hs))
    zhi = (b + hb) / (se - hs)
    ps = [norm_p2(zlo), norm_p2(zhi), t_p2(zlo, 40), t_p2(zhi, 40)]
    return {level_for_p(p, legend) for p in ps}


COEF_RE = re.compile(r'^\s*(?P<num>' + NUMS + r')\s*(?P<marks>' + MARKS + r')\s*$')
SE_RE = re.compile(r'^\s*\(\s*(?P<num>' + NUMU + r')\s*\)\s*(?P<marks>' + MARKS + r')\s*$')
PLAIN_RE = re.compile(r'^\s*(?P<num>' + NUMS + r')\s*$')
SAME_RE = re.compile(r'^\s*(?P<num>' + NUMS + r')\s*(?P<m1>' + MARKS + r')\s*\(\s*(?P<se>' + NUMU +
                     r')\s*\)\s*(?P<m2>' + MARKS + r')\s*$')
SE_HDR = re.compile(r'^\s*(?:robust\s+)?(?:s\.?\s?e\.?|std\.?\s*err(?:or|\.)?|standard\s+error)s?\s*$', re.I)
RATIO_TABLE = re.compile(
    r'odds ratio|hazard ratio|incidence rate ratio|relative risk|risk ratio|exponentiated|exp\s*\(|'
    r'\bIRR\b|\bOR\b|\bHR\b|\bRRR\b|\bAOR\b|(?<![tz][- ])\bratios?\b', re.I)
MSD_HDR = re.compile(r'\bM\s*\(\s*SD\s*\)|\bmeans?\b|\bSD\b|\bS\.D\.|standard deviation|std\.?\s*dev|average', re.I)
NOT_SE_PAREN = re.compile(
    r'\b[tz][- ]?(?:stat|statistic|value|score|ratio)s?\b[^.]{0,40}(?:parenthes|bracket)|'
    r'(?:parenthes|bracket)[^.]{0,40}\b[tz][- ]?(?:stat|statistic|value|score|ratio)|'
    r'p[- ]?values?[^.]{0,30}parenthes|parenthes[^.]{0,40}p[- ]?values?|'
    r'confidence interval[^.]{0,30}parenthes|parenthes[^.]{0,40}confidence', re.I)
NEG_PAREN = re.compile(r'^\s*\(\s*[−\-–]\s*\d')
SD_PAREN = re.compile(r'standard deviations?[^.]{0,30}parenthes|parenthes[^.]{0,40}standard deviation', re.I)
SE_EXPLICIT = re.compile(r'standard error|\bS\.?E\.?s?\b|std\.?\s*err', re.I)


def coef_pairs(t):
    """All coefficient / standard-error pairs in a table (layout-aware)."""
    rows = t.get('rows') or []
    pairs = []
    used = set()

    def mk(layout, crc, src, cm, sm, cgrp, sgrp, cmk, smk):
        b_str, se_str = cm.group(cgrp), sm.group(sgrp)
        m_c, m_s = cm.group(cmk), sm.group(smk)
        if m_c and m_s:
            return None
        if '.' not in b_str or '.' not in se_str:
            return None
        if MSD_HDR.search(col_header(t, crc[1])) or MSD_HDR.search(col_header(t, src[1])):
            return None  # descriptive mean (sd) cells, not coefficient / SE
        if m_s:
            marks, mrc, mspan = m_s, src, sm.span(smk)
        else:
            marks, mrc, mspan = m_c, crc, cm.span(cmk)
        return dict(layout=layout, coef_rc=crc, se_rc=src, b_str=b_str, se_str=se_str,
                    b_span=cm.span(cgrp), se_span=sm.span(sgrp), marks=marks, marks_rc=mrc,
                    marks_span=mspan, marks_on_se=bool(m_s))

    for r, row in enumerate(rows):
        for c in range(1, len(row)):
            m = SAME_RE.match(row[c])
            if m:
                p = mk('same', (r, c), (r, c), m, m, 'num', 'se', 'm1', 'm2')
                if p:
                    pairs.append(p)
                    used.add((r, c))
    for r in range(len(rows) - 1):
        for c in range(1, min(len(rows[r]), len(rows[r + 1]))):
            if (r, c) in used:
                continue
            mc, ms = COEF_RE.match(rows[r][c]), SE_RE.match(rows[r + 1][c])
            if mc and ms and (rows[r + 1][0].strip() in ('', rows[r][0].strip())):
                p = mk('stacked', (r, c), (r + 1, c), mc, ms, 'num', 'num', 'marks', 'marks')
                if p:
                    pairs.append(p)
    hdr_n = max((len(h) for h in t.get('header_rows') or []), default=0)
    for c in range(2, hdr_n):
        if SE_HDR.match(col_header(t, c).strip() or 'x') and not SE_HDR.match(col_header(t, c - 1).strip() or 'x'):
            for r, row in enumerate(rows):
                if c >= len(row):
                    continue
                mc = COEF_RE.match(row[c - 1])
                ms = PLAIN_RE.match(row[c]) or SE_RE.match(row[c])
                if mc and ms and not to_float(ms.group('num')) < 0:
                    if 'marks' not in ms.groupdict():
                        ms = re.match(r'^\s*(?P<num>' + NUMU + r')\s*(?P<marks>)$', row[c])
                        if not ms:
                            continue
                    p = mk('side', (r, c - 1), (r, c), mc, ms, 'num', 'num', 'marks', 'marks')
                    if p:
                        pairs.append(p)
    return pairs


def star_sites(paper):
    """Coef/SE pairs in tables with a valid legend whose printed stars are robustly consistent."""
    sites = []
    legends = {t['table_id']: parse_legend((t.get('footnotes') or '') + ' ' + (t.get('caption') or ''))
               for t in paper.get('tables', [])}
    distinct = {json.dumps(v) for v in legends.values() if v}
    paper_legend = json.loads(distinct.pop()) if len(distinct) == 1 else None
    if any(NOT_SE_PAREN.search(table_text(t)) for t in paper.get('tables', [])):
        paper_legend = None  # some tables print t/z/p in parentheses: don't guess for unlabelled ones
    for t in paper.get('tables', []):
        tt = table_text(t)
        if RATIO_TABLE.search(tt) or NOT_SE_PAREN.search(tt):
            continue
        if any(NEG_PAREN.match(x) for row in t.get('rows', []) for x in row):
            continue  # negative numbers in parentheses cannot be standard errors
        # A table without its own legend inherits the paper's legend if all tables agree on one.
        legend = legends[t['table_id']] or paper_legend
        if not legend:
            continue
        legend = [tuple(x) for x in legend]
        pairs = coef_pairs(t)
        good, n = [], 0
        for p in pairs:
            lvl = level_of(p['marks'], legend)
            if lvl is None:
                continue
            n += 1
            rl = robust_levels(p['b_str'], p['se_str'], legend)
            if rl is not None and rl == {lvl}:
                p = dict(p, table_id=t['table_id'], level=lvl, legend=legend,
                         stars_on_se=any(q['marks_on_se'] for q in pairs))
                good.append(p)
        need = 2 if SE_EXPLICIT.search(tt) and not SD_PAREN.search(tt) else 3
        if n >= need and len(good) >= 0.7 * n:
            sites.extend(good)
    return sites


# --------------------------------------------------------------------------- error types

def sites_star_flip(paper):
    return star_sites(paper)


def apply_star_flip(paper, s, rng):
    L = len(s['legend'])
    k = s['level']
    opts = [x for x in (k - 1, k + 1) if 0 <= x <= L]
    newk = rng.choice(opts)
    newmarks = '' if newk == 0 else s['legend'][newk - 1][0]
    if s['marks'] == '' and s['stars_on_se'] and s['layout'] == 'stacked':
        r, c = s['se_rc']
        t = _table(paper, s['table_id'])
        cell = t['rows'][r][c]
        m = SE_RE.match(cell)
        a, b = m.span('marks')
    else:
        r, c = s['marks_rc']
        t = _table(paper, s['table_id'])
        cell = t['rows'][r][c]
        a, b = s['marks_span']
    new = splice(cell, a, b, newmarks)
    set_cell(paper, s['table_id'], r, c, new)
    return cell_record(s['table_id'], r, c, cell, new)


def _table(paper, tid):
    for t in paper['tables']:
        if t['table_id'] == tid:
            return t
    raise KeyError(tid)


def sites_se_misreport(paper):
    return star_sites(paper)


def apply_se_misreport(paper, s, rng):
    k = s['level']
    b = abs(to_float(s['b_str']))
    se = to_float(s['se_str'])
    z0 = b / se
    cands = []
    for kind, ns in num_variants(s['se_str']):
        if '.' not in ns or to_float(ns) <= 0:
            continue
        rl = robust_levels(s['b_str'], ns, s['legend'])
        if rl is None or k in rl or len(rl) != 1:
            continue
        z1 = b / to_float(ns)
        if not (z1 / z0 >= 2.0 or z1 / z0 <= 0.5):
            continue
        cands.append((kind, ns))
    if not cands:
        return None
    shifts = [c for c in cands if c[0].startswith('shift')]
    kind, ns = rng.choice(shifts if shifts and rng.random() < 0.5 else cands)
    r, c = s['se_rc']
    t = _table(paper, s['table_id'])
    cell = t['rows'][r][c]
    a, bb = s['se_span']
    new = splice(cell, a, bb, ns)
    set_cell(paper, s['table_id'], r, c, new)
    return cell_record(s['table_id'], r, c, cell, new)


def _norm_core(core):
    return core.replace(',', '')


def _sign_of(sign):
    return '-' if sign else '+'


def sites_sign_flip(paper):
    cell_count = {}
    coefs = []
    for t in paper.get('tables', []):
        if RATIO_TABLE.search(table_text(t)):
            continue
        for p in coef_pairs(t):
            coefs.append((t['table_id'], p))
        for row in t.get('rows', []):
            for cell in row:
                for m in TOK_RE.finditer(cell):
                    k = _norm_core(m.group('core'))
                    cell_count[k] = cell_count.get(k, 0) + 1
    sites = []
    for tid, p in coefs:
        sign, core = _split_sign(p['b_str'])
        if to_float(core) == 0 or len(core.replace('.', '').lstrip('0')) < 2:
            continue
        if cell_count.get(core, 0) != 1:
            continue
        for para in paper.get('paragraphs', []):
            hit = False
            for m in TXT_TOK_RE.finditer(para['text']):
                if m.group('core') == core and _sign_of(m.group('sign')) == _sign_of(sign):
                    hit = True
                    break
            if hit:
                sites.append(dict(table_id=tid, pair=p, para_id=para['para_id']))
                break
    return sites


def apply_sign_flip(paper, s, rng):
    p = s['pair']
    r, c = p['coef_rc']
    t = _table(paper, s['table_id'])
    cell = t['rows'][r][c]
    a, b = p['b_span']
    sign, core = _split_sign(p['b_str'])
    if sign:
        nb = core
    else:
        uses_true_minus = any('−' in x for tt in paper['tables'] for row in tt['rows'] for x in row)
        nb = ('−' if uses_true_minus else '-') + core
    new = splice(cell, a, b, nb)
    set_cell(paper, s['table_id'], r, c, new)
    return cell_record(s['table_id'], r, c, cell, new)


N_ROW = re.compile(r'^\s*(?:n|n\.|obs\.?|observations|no\.? of (?:obs\.?|observations)|number of (?:obs\.?|observations|'
                   r'respondents|cases|participants)|num\.? obs\.?|sample size|respondents|participants|cases)\s*$', re.I)
N_EQ = re.compile(r'\b[Nn]\s*=\s*(' + INT + r')\b')
N_CTX = re.compile(r'\b[Nn]\s*=|respondent|participant|observation|individual|subject|sample|\bcases\b|voter|'
                   r'citizen|people|persons|users|legislator|countries|surveyed|interviewee|households|students|'
                   r'members|completed', re.I)


def table_ns(paper):
    ns = set()
    for t in paper.get('tables', []):
        for row in t.get('rows', []):
            if row and N_ROW.match(row[0]):
                for cell in row[1:]:
                    m = re.match(r'^\s*(' + INT + r')\s*$', cell)
                    if m:
                        ns.add(int(m.group(1).replace(',', '')))
        for s in [t.get('caption') or '', t.get('footnotes') or ''] + [x for hr in t.get('header_rows') or [] for x in hr]:
            for m in N_EQ.finditer(s):
                ns.add(int(m.group(1).replace(',', '')))
    return {n for n in ns if n >= 50}


def sites_n_text_mismatch(paper):
    ns = table_ns(paper)
    if not ns:
        return []
    sites = []
    for para in paper.get('paragraphs', []):
        txt = para['text']
        for m in TXT_TOK_RE.finditer(txt):
            core = m.group('core')
            if '.' in core or m.group('sign'):
                continue
            v = int(core.replace(',', ''))
            if v not in ns:
                continue
            if 1900 <= v <= 2035 and ',' not in core and not re.search(r'[Nn]\s*=\s*$', txt[max(0, m.start() - 6):m.start()]):
                continue
            ctx = txt[max(0, m.start() - 60):m.end() + 60]
            if not N_CTX.search(ctx):
                continue
            sites.append(dict(para_id=para['para_id'], start=m.start('core'), end=m.end('core'), core=core))
    return sites


def apply_n_text_mismatch(paper, s, rng):
    cores = all_table_cores(paper)
    ns = table_ns(paper)
    cands = []
    for kind, nv in num_variants(s['core'], kinds=('transpose', 'digit')):
        v = int(nv.replace(',', ''))
        if str(v) in cores or v in ns or v < 10:
            continue
        cands.append(nv)
    if not cands:
        return None
    nv = rng.choice(cands)
    return text_record(paper, s['para_id'], s['start'], s['end'], nv, s['start'], s['end'])


PCT_CELL = re.compile(r'^\s*(?P<num>' + NUMU + r')\s*%\s*$')
PCT_PLAIN = re.compile(r'^\s*(?P<num>' + NUMU + r')\s*$')
PCT_CTX = re.compile(r'%|percent', re.I)


def _pct_val(cell, ctx):
    m = PCT_CELL.match(cell) or (PCT_PLAIN.match(cell) if ctx else None)
    if not m:
        return None
    v = to_float(m.group('num'))
    if 0 <= v <= 100:
        return m.group('num'), m.span('num'), v
    return None


def _segments(run):
    """run: list of (pos, numstr, span, val). Partition into contiguous ~100 segments of size>=3."""
    maxdec = max(decimals(x[1]) for x in run)
    segs = []
    i = 0
    while i < len(run):
        s = 0.0
        found = None
        for j in range(i, len(run)):
            s += run[j][3]
            k = j - i + 1
            tol = min(1.5, max(0.5 * k * 10 ** (-maxdec), 0.2))
            if abs(s - 100) <= tol:
                nxt = run[j + 1][3] if j + 1 < len(run) else None
                if nxt is None or nxt > 0:
                    found = (j, tol)
                    break
            if s > 100 + 1.5:
                break
        if found is None:
            return segs if segs else []
        j, tol = found
        if j - i + 1 >= 3:
            segs.append((run[i:j + 1], tol))
        i = j + 1
    return segs


def sites_pct_sum_break(paper):
    sites = []
    for t in paper.get('tables', []):
        rows = t.get('rows') or []
        tctx = bool(PCT_CTX.search((t.get('caption') or '') + ' ' + (t.get('footnotes') or '')))
        ncol = max((len(r) for r in rows), default=0)
        # columns
        for c in range(1, ncol):
            ctx = tctx or bool(PCT_CTX.search(col_header(t, c)))
            run = []
            for r in range(len(rows) + 1):
                pv = _pct_val(rows[r][c], ctx) if r < len(rows) and c < len(rows[r]) else None
                if pv:
                    run.append(((r, c),) + pv)
                else:
                    if len(run) >= 3:
                        for seg, tol in _segments(run):
                            sites.append(dict(table_id=t['table_id'], seg=seg, tol=tol))
                    run = []
        for r, row in enumerate(rows):
            ctx = tctx or bool(PCT_CTX.search(row[0] if row else ''))
            run = []
            for c in range(1, len(row) + 1):
                pv = _pct_val(row[c], ctx) if c < len(row) else None
                if pv:
                    run.append(((r, c),) + pv)
                else:
                    if len(run) >= 3:
                        for seg, tol in _segments(run):
                            sites.append(dict(table_id=t['table_id'], seg=seg, tol=tol))
                    run = []
    return sites


def apply_pct_sum_break(paper, s, rng):
    seg, tol = s['seg'], s['tol']
    total = sum(x[3] for x in seg)
    cands = []
    for idx, (rc, numstr, span, v) in enumerate(seg):
        for kind, nv in num_variants(numstr):
            if nv.startswith(('-', '−')):
                continue
            fv = to_float(nv)
            if not (0 <= fv <= 100):
                continue
            dev = abs(total - v + fv - 100)
            if dev >= max(3.0, 2 * tol) and dev <= 40:
                cands.append((idx, nv))
    if not cands:
        return None
    idx, nv = rng.choice(cands)
    (r, c), numstr, (a, b), v = seg[idx]
    t = _table(paper, s['table_id'])
    cell = t['rows'][r][c]
    new = splice(cell, a, b, nv)
    set_cell(paper, s['table_id'], r, c, new)
    return cell_record(s['table_id'], r, c, cell, new)


MEAN_H = re.compile(r'^\s*(?:mean|m|average|avg\.?)\b', re.I)
MIN_H = re.compile(r'^\s*(?:min\.?|minimum)\b', re.I)
MAX_H = re.compile(r'^\s*(?:max\.?|maximum)\b', re.I)


def sites_mean_outside_range(paper):
    sites = []
    for t in paper.get('tables', []):
        rows = t.get('rows') or []
        hdrs = t.get('header_rows') or []
        # variables in rows, stats in columns
        for hr in hdrs:
            cm = [i for i, h in enumerate(hr) if MEAN_H.match(h)]
            cn = [i for i, h in enumerate(hr) if MIN_H.match(h)]
            cx = [i for i, h in enumerate(hr) if MAX_H.match(h)]
            if len(cm) == 1 and len(cn) == 1 and len(cx) == 1:
                for r, row in enumerate(rows):
                    _mean_site(sites, t, row, r, cm[0], cn[0], cx[0], orient='row')
        # stats in rows, variables in columns
        lab = {}
        for r, row in enumerate(rows):
            if not row:
                continue
            for name, rx in (('mean', MEAN_H), ('min', MIN_H), ('max', MAX_H)):
                if rx.match(row[0]):
                    lab.setdefault(name, []).append(r)
        if all(len(lab.get(k, [])) == 1 for k in ('mean', 'min', 'max')):
            rm, rn, rx_ = lab['mean'][0], lab['min'][0], lab['max'][0]
            for c in range(1, min(len(rows[rm]), len(rows[rn]), len(rows[rx_]))):
                cells = (rows[rm][c], rows[rn][c], rows[rx_][c])
                ms = [PLAIN_RE.match(x) for x in cells]
                if all(ms):
                    _mean_check(sites, t, (rm, c), ms)
    return sites


def _mean_site(sites, t, row, r, im, imin, imax, orient):
    if max(im, imin, imax) >= len(row):
        return
    ms = [PLAIN_RE.match(row[i]) for i in (im, imin, imax)]
    if all(ms):
        _mean_check(sites, t, (r, im), ms)


def _mean_check(sites, t, rc, ms):
    mean, mn, mx = (to_float(m.group('num')) for m in ms)
    if mn < mx and mn <= mean <= mx:
        sites.append(dict(table_id=t['table_id'], rc=rc, mean_str=ms[0].group('num'),
                          span=ms[0].span('num'), mn=mn, mx=mx))


def apply_mean_outside_range(paper, s, rng):
    rng_w = s['mx'] - s['mn']
    lim = 10 * max(abs(s['mn']), abs(s['mx']), 1.0)
    cands = []
    for kind, nv in num_variants(s['mean_str']):
        fv = to_float(nv)
        margin = max(0.05 * rng_w, 2 * half_unit(nv))
        if (fv < s['mn'] - margin or fv > s['mx'] + margin) and abs(fv) <= lim:
            cands.append(nv)
    if not cands:
        return None
    nv = rng.choice(cands)
    r, c = s['rc']
    t = _table(paper, s['table_id'])
    cell = t['rows'][r][c]
    a, b = s['span']
    new = splice(cell, a, b, nv)
    set_cell(paper, s['table_id'], r, c, new)
    return cell_record(s['table_id'], r, c, cell, new)


SEP = r'\s*(?:,|;|–|—|\bto\b|-)\s*'
CI_CELL = re.compile(r'^\s*(?P<est>' + NUMS + r')\s*(?P<m1>' + MARKS + r')\s*[\(\[]\s*(?P<lo>' + NUMS + r')' + SEP +
                     r'(?P<hi>' + NUMS + r')\s*[\)\]]\s*(?P<m2>' + MARKS + r')\s*$')
CI_ONLY = re.compile(r'^\s*[\(\[]?\s*(?P<lo>' + NUMS + r')' + SEP + r'(?P<hi>' + NUMS + r')\s*[\)\]]?\s*$')
CI_HDR = re.compile(r'\bCI\b|confidence|credible|interval', re.I)
LO_HDR = re.compile(r'^\s*(?:lower|LL|LB|low)\b', re.I)
HI_HDR = re.compile(r'^\s*(?:upper|UL|UB|high)\b', re.I)
NONEST_HDR = re.compile(r'^\s*(?:s\.?e\.?|std|standard|p\b|p-|t\b|z\b|sig|wald|df|n\b)', re.I)
CI_TXT = re.compile(
    r'(?<![\w.])(?P<est>' + NUMS + r')(?P<mid>\s*[,;]?\s*\(?\s*(?:95\s?%\s*)?(?:CI|C\.I\.|CrI|confidence interval|'
    r'credible interval)s?\s*[:=,]?\s*|\s*)(?P<open>[\(\[])?\s*(?P<lo>' + NUMS + r')' + SEP + r'(?P<hi>' + NUMS +
    r')(?![\d.])\s*(?P<close>[\)\]])?')


def _ci_ok(est, lo, hi):
    e, l, h = to_float(est), to_float(lo), to_float(hi)
    return l < h and l <= e <= h


def sites_ci_excludes_estimate(paper):
    sites = []
    for t in paper.get('tables', []):
        tt = table_text(t)
        rows = t.get('rows') or []
        if CI_HDR.search(tt):
            for r, row in enumerate(rows):
                for c in range(1, len(row)):
                    m = CI_CELL.match(row[c])
                    if m and _ci_ok(m['est'], m['lo'], m['hi']):
                        sites.append(dict(kind='cell', table_id=t['table_id'], parts={
                            k: ((r, c), m.group(k), m.span(k)) for k in ('est', 'lo', 'hi')}))
        ncol = max((len(r) for r in rows), default=0)
        for c in range(2, ncol):
            h = col_header(t, c)
            two = None
            if LO_HDR.search(h.split()[-1] if h.split() else '') or LO_HDR.match(h):
                two = 'lohi'
            elif CI_HDR.search(h) and c + 1 < ncol and col_header(t, c + 1) == h:
                two = 'span'
            elif not CI_HDR.search(h):
                continue
            for r, row in enumerate(rows):
                if c >= len(row):
                    continue
                ec = None
                for cc in range(c - 1, 0, -1):
                    if NONEST_HDR.match((t['header_rows'][-1][cc] if t.get('header_rows') and cc < len(t['header_rows'][-1]) else '')):
                        continue
                    if COEF_RE.match(row[cc]):
                        ec = cc
                    break
                if ec is None:
                    continue
                me = COEF_RE.match(row[ec])
                if two and c + 1 < len(row):
                    ml, mh = PLAIN_RE.match(row[c]), PLAIN_RE.match(row[c + 1])
                    if two == 'lohi' and not HI_HDR.match(col_header(t, c + 1).split()[-1] if col_header(t, c + 1).split() else ''):
                        continue
                    if ml and mh and _ci_ok(me['num'], ml['num'], mh['num']):
                        sites.append(dict(kind='cell', table_id=t['table_id'], parts={
                            'est': ((r, ec), me['num'], me.span('num')),
                            'lo': ((r, c), ml['num'], ml.span('num')),
                            'hi': ((r, c + 1), mh['num'], mh.span('num'))}))
                elif not two:
                    mc = CI_ONLY.match(row[c])
                    if mc and _ci_ok(me['num'], mc['lo'], mc['hi']):
                        sites.append(dict(kind='cell', table_id=t['table_id'], parts={
                            'est': ((r, ec), me['num'], me.span('num')),
                            'lo': ((r, c), mc['lo'], mc.span('lo')),
                            'hi': ((r, c), mc['hi'], mc.span('hi'))}))
    for para in paper.get('paragraphs', []):
        for m in ci_text_matches(para['text']):
            sites.append(dict(kind='text', para_id=para['para_id'], m=m))
    return sites


def ci_text_matches(txt):
    out = []
    for m in CI_TXT.finditer(txt):
        has_ci = bool(re.search(r'CI|C\.I\.|CrI|interval', m.group('mid')))
        bracketed = bool(m.group('open')) and bool(m.group('close'))
        if not (has_ci or (bracketed and m.group('mid').strip() == '')):
            continue
        if not has_ci and not all('.' in m.group(k) for k in ('est', 'lo', 'hi')):
            continue
        if bracketed and (m.group('open') + m.group('close')) not in ('()', '[]'):
            continue
        if _ci_ok(m['est'], m['lo'], m['hi']):
            out.append({k: (m.group(k), m.span(k)) for k in ('est', 'lo', 'hi')} | {'span': m.span()})
    return out


def _ci_mut_ok(est, lo, hi):
    e, l, h = to_float(est), to_float(lo), to_float(hi)
    if not l < h:
        return False
    margin = 2 * max(half_unit(est), half_unit(lo), half_unit(hi))
    return e < l - margin or e > h + margin


def apply_ci_excludes_estimate(paper, s, rng):
    if s['kind'] == 'cell':
        vals = {k: v[1] for k, v in s['parts'].items()}
    else:
        vals = {k: s['m'][k][0] for k in ('est', 'lo', 'hi')}
    cands = []
    for k in ('est', 'lo', 'hi'):
        for kind, nv in num_variants(vals[k]):
            trial = dict(vals, **{k: nv})
            if _ci_mut_ok(trial['est'], trial['lo'], trial['hi']):
                cands.append((k, nv))
    if not cands:
        return None
    k, nv = rng.choice(cands)
    if s['kind'] == 'cell':
        (r, c), old, (a, b) = s['parts'][k]
        t = _table(paper, s['table_id'])
        cell = t['rows'][r][c]
        new = splice(cell, a, b, nv)
        set_cell(paper, s['table_id'], r, c, new)
        return cell_record(s['table_id'], r, c, cell, new)
    a, b = s['m'][k][1]
    ss, se = s['m']['span']
    return text_record(paper, s['para_id'], a, b, nv, ss, se)


STAT_TXT = re.compile(
    r'(?<![\w])(?P<stat>t|F|χ\s?2|χ²|chi2|chi-squared?|X2|z|Z|r)\s*(?:\(\s*(?P<df1>\d+(?:\.\d+)?)\s*'
    r'(?:,\s*(?:N\s*=\s*)?(?P<df2>' + INT + r'(?:\.\d+)?))?\s*\))?\s*=\s*(?P<val>' + NUMS + r')\s*[,;]\s*'
    r'(?P<pl>[pP])\s*(?P<op>[<=>≤≥])\s*(?P<p>\d?\.\d+|1\.0+|1)(?![\d.])')


def _stat_kind(s):
    s = s.replace(' ', '')
    if s == 't':
        return 't'
    if s == 'F':
        return 'F'
    if s in ('z', 'Z'):
        return 'z'
    if s == 'r':
        return 'r'
    return 'chi2'


def stat_prange(m):
    kind = _stat_kind(m.group('stat'))
    df1 = float(m.group('df1')) if m.group('df1') else None
    df2 = float(m.group('df2').replace(',', '')) if m.group('df2') else None
    if kind in ('t', 'r', 'chi2') and df1 is None:
        return None
    if kind == 'F' and (df1 is None or df2 is None):
        return None
    if kind == 'chi2' and df2 is not None and m.group('df2') and 'N' not in m.group(0)[:m.start('df2') - m.start()]:
        return None
    if kind == 'z' and df1 is not None:
        return None
    return _prange(kind, m.group('val'), df1, df2)


def _stat_args(m):
    kind = _stat_kind(m.group('stat'))
    df1 = float(m.group('df1')) if m.group('df1') else None
    df2 = float(m.group('df2').replace(',', '')) if m.group('df2') else None
    return kind, df1, df2


def _prange(kind, vs, df1, df2):
    v = abs(to_float(vs))
    h = half_unit(vs)
    if kind == 'r' and v + h >= 1:
        return None
    lo_v, hi_v = max(0.0, v - h), v + h
    return stat_p(kind, hi_v, df1, df2), stat_p(kind, lo_v, df1, df2)


def p_consistent(op, pstr, pmin, pmax):
    p = to_float(pstr)
    hp = half_unit(pstr)
    if op == '=':
        return pmax >= p - hp and pmin <= p + hp
    if op in '<≤':
        return pmin < p or (op == '≤' and pmin <= p)
    return pmax > p or (op == '≥' and pmax >= p)


def p_gross(op, pstr, pmin, pmax):
    p = to_float(pstr)
    hp = half_unit(pstr)
    if not (0 < p <= 1):
        return False
    if op == '=':
        return p + hp < pmin / 2 or p - hp > min(1.0, pmax * 2)
    if op in '<≤':
        return p <= pmin / 2
    return p >= pmax * 2


def sites_statcheck_p_mismatch(paper):
    sites = []
    for para in paper.get('paragraphs', []):
        for m in STAT_TXT.finditer(para['text']):
            pr = stat_prange(m)
            if pr is None:
                continue
            if p_consistent(m.group('op'), m.group('p'), *pr):
                kind, df1, df2 = _stat_args(m)
                sites.append(dict(para_id=para['para_id'], span=m.span(), p=(m.group('p'), m.span('p')),
                                  op=(m.group('op'), m.span('op')), prange=pr,
                                  val=(m.group('val'), m.span('val')), kind=kind, df1=df1, df2=df2))
    return sites


def apply_statcheck_p_mismatch(paper, s, rng):
    pmin, pmax = s['prange']
    op, (oa, ob) = s['op']
    pstr, (pa, pb) = s['p']
    vstr, (va, vb) = s['val']
    cands = {'p': [], 'op': [], 'stat': []}
    for kind, nv in num_variants(pstr):
        if p_gross(op, nv, pmin, pmax):
            cands['p'].append(nv)
    flip = {'<': '>', '>': '<', '≤': '≥', '≥': '≤'}
    if op in flip and p_gross(flip[op], pstr, pmin, pmax):
        cands['op'].append(flip[op])
    for kind, nv in num_variants(vstr):
        pr = _prange(s['kind'], nv, s['df1'], s['df2'])
        if pr is not None and p_gross(op, pstr, *pr):
            cands['stat'].append(nv)
    cats = [k for k in ('p', 'op', 'stat') if cands[k]]
    if not cats:
        return None
    what = rng.choice(cats)
    nv = rng.choice(cands[what])
    a, b = {'p': (pa, pb), 'op': (oa, ob), 'stat': (va, vb)}[what]
    ss, se = s['span']
    return text_record(paper, s['para_id'], a, b, nv, ss, se)


P_BEFORE = re.compile(r'[pP]\s*[<=>≤≥]\s*$')


def sites_text_table_transposition(paper):
    cores = all_table_cores(paper)
    sites = []
    for para in paper.get('paragraphs', []):
        txt = para['text']
        excl = [x['span'] for x in ci_text_matches(txt)] + [m.span() for m in STAT_TXT.finditer(txt)]
        for m in TXT_TOK_RE.finditer(txt):
            core = m.group('core')
            if '.' not in core or ',' in core:
                continue
            if sum(ch.isdigit() for ch in core.lstrip('0.')) < 2:
                continue
            if core not in cores:
                continue
            a, b = m.span('core')
            if any(x <= a < y for x, y in excl):
                continue
            if P_BEFORE.search(txt[max(0, m.start() - 8):m.start()]):
                continue
            sites.append(dict(para_id=para['para_id'], start=a, end=b, core=core))
    return sites


def apply_text_table_transposition(paper, s, rng):
    cores = all_table_cores(paper)
    cands = [nv for kind, nv in num_variants(s['core'], kinds=('transpose',)) if nv not in cores]
    if not cands:
        cands = [nv for kind, nv in num_variants(s['core'], kinds=('digit',)) if nv not in cores]
    if not cands:
        return None
    nv = rng.choice(cands)
    return text_record(paper, s['para_id'], s['start'], s['end'], nv, s['start'], s['end'])


CP1 = re.compile(r'^\s*(?P<cnt>' + INT + r')\s*\(\s*(?P<pct>' + NUMU + r')\s*(?P<pc>%?)\s*\)\s*$')
CP2 = re.compile(r'^\s*(?P<pct>' + NUMU + r')\s*%\s*\(\s*(?P<cnt>' + INT + r')\s*\)\s*$')


def _cp_parse(cell, ctx):
    m = CP2.match(cell)
    if not m:
        m = CP1.match(cell)
        if m and not (m.group('pc') or ctx):
            return None
    if not m:
        return None
    cnt = int(m.group('cnt').replace(',', ''))
    pct = to_float(m.group('pct'))
    if pct > 100:
        return None
    return dict(cnt=cnt, cnt_s=m.group('cnt'), pct_s=m.group('pct'), cnt_span=m.span('cnt'),
                pct_span=m.span('pct'))


def _n_interval(cnt, pct_s):
    p = to_float(pct_s)
    h = half_unit(pct_s)
    if cnt == 0:
        return None
    if p + h <= 0:
        return (math.inf, math.inf)
    lo = 100.0 * cnt / (p + h)
    hi = 100.0 * cnt / (p - h) if p - h > 0 else math.inf
    return (lo, hi)


def _intersect(ivs):
    lo = max(i[0] for i in ivs)
    hi = min(i[1] for i in ivs)
    return lo, hi


def _has_int(lo, hi):
    return math.ceil(lo - 1e-9) <= hi + 1e-9


def sites_count_pct_mismatch(paper):
    sites = []
    for t in paper.get('tables', []):
        rows = t.get('rows') or []
        tctx = bool(re.search(r'%|percent|n\s*\(\s*%\s*\)', table_text(t), re.I))
        ncol = max((len(r) for r in rows), default=0)
        groups = []
        for c in range(1, ncol):
            run = []
            for r in range(len(rows) + 1):
                cp = _cp_parse(rows[r][c], tctx) if r < len(rows) and c < len(rows[r]) else None
                if cp:
                    run.append(((r, c), cp))
                else:
                    groups.append(run)
                    run = []
        for r, row in enumerate(rows):
            run = []
            for c in range(1, len(row) + 1):
                cp = _cp_parse(row[c], tctx) if c < len(row) else None
                if cp:
                    run.append(((r, c), cp))
                else:
                    groups.append(run)
                    run = []
        for g in groups:
            ivs = [_n_interval(cp['cnt'], cp['pct_s']) for _, cp in g]
            ivs_nz = [iv for iv in ivs if iv is not None]
            if len(ivs_nz) < 3:
                continue
            lo, hi = _intersect(ivs_nz)
            if lo <= hi and _has_int(lo, hi):
                for i, (rc, cp) in enumerate(g):
                    if cp['cnt'] == 0:
                        continue
                    others = [iv for j, iv in enumerate(ivs) if j != i and iv is not None]
                    olo, ohi = _intersect(others)
                    sites.append(dict(table_id=t['table_id'], rc=rc, cp=cp, olo=olo, ohi=ohi))
    return sites


def apply_count_pct_mismatch(paper, s, rng):
    cp = s['cp']
    olo, ohi = s['olo'] * 0.97, s['ohi'] * 1.03
    cands = []
    for kind, nv in num_variants(cp['pct_s'], kinds=('transpose', 'digit')):
        if nv.startswith(('-', '−')) or to_float(nv) > 100 or to_float(nv) <= 0:
            continue
        iv = _n_interval(cp['cnt'], nv)
        if iv and (iv[1] < olo or iv[0] > ohi):
            cands.append(('pct', nv))
    for kind, nv in num_variants(cp['cnt_s'], kinds=('transpose', 'digit')):
        if '.' in nv or nv.startswith(('-', '−')):
            continue
        c2 = int(nv.replace(',', ''))
        if c2 <= 0 or c2 > s['ohi']:
            continue
        iv = _n_interval(c2, cp['pct_s'])
        if iv and (iv[1] < olo or iv[0] > ohi):
            cands.append(('cnt', nv))
    if not cands:
        return None
    what, nv = rng.choice(cands)
    r, c = s['rc']
    t = _table(paper, s['table_id'])
    cell = t['rows'][r][c]
    a, b = cp['pct_span'] if what == 'pct' else cp['cnt_span']
    new = splice(cell, a, b, nv)
    set_cell(paper, s['table_id'], r, c, new)
    return cell_record(s['table_id'], r, c, cell, new)


ERROR_TYPES = {
    'star_flip': (sites_star_flip, apply_star_flip),
    'se_misreport': (sites_se_misreport, apply_se_misreport),
    'sign_flip': (sites_sign_flip, apply_sign_flip),
    'n_text_mismatch': (sites_n_text_mismatch, apply_n_text_mismatch),
    'pct_sum_break': (sites_pct_sum_break, apply_pct_sum_break),
    'mean_outside_range': (sites_mean_outside_range, apply_mean_outside_range),
    'ci_excludes_estimate': (sites_ci_excludes_estimate, apply_ci_excludes_estimate),
    'statcheck_p_mismatch': (sites_statcheck_p_mismatch, apply_statcheck_p_mismatch),
    'count_pct_mismatch': (sites_count_pct_mismatch, apply_count_pct_mismatch),
    'text_table_transposition': (sites_text_table_transposition, apply_text_table_transposition),
}


# --------------------------------------------------------------------------- driver

def load_split(path=os.path.join(HERE, 'split.json')):
    with open(path) as f:
        sp = json.load(f)
    assert not set(sp['visible']) & set(sp['heldout'])
    assert set(sp['visible']) | set(sp['heldout']) == set(ERROR_TYPES)
    return sp


def plant_paper(paper, etype, seed, max_copies=MAX_COPIES):
    """Return list of (planted_paper, location, before, after) for one paper and type."""
    finder, applier = ERROR_TYPES[etype]
    try:
        sites = finder(paper)
    except Exception as e:  # never let one odd paper kill the run
        print(f"[warn] {paper.get('paper_id')} {etype} finder: {e!r}", file=sys.stderr)
        return [], 0
    rng = random.Random(f"{seed}:{paper['paper_id']}:{etype}")
    order = list(range(len(sites)))
    rng.shuffle(order)
    out = []
    for i in order:
        if len(out) >= max_copies:
            break
        cp = copy.deepcopy(paper)
        try:
            res = applier(cp, sites[i], rng)
        except Exception as e:
            print(f"[warn] {paper.get('paper_id')} {etype} apply: {e!r}", file=sys.stderr)
            res = None
        if res is None:
            continue
        loc, before, after = res
        if before == after or cp == paper:
            continue
        out.append((cp, loc, before, after))
    return out, len(sites)


def run(corpus_dir, out_dir, split=None, verbose=True):
    split = split or load_split()
    seed = split['seed']
    files = sorted(f for f in os.listdir(corpus_dir) if f.endswith('.json'))
    type_split = {t: 'visible' for t in split['visible']} | {t: 'heldout' for t in split['heldout']}
    for sp in ('visible', 'heldout'):
        d = os.path.join(out_dir, sp)
        os.makedirs(d, exist_ok=True)
        for f in os.listdir(d):
            if f.endswith('.json') or f == 'manifest.jsonl':
                os.remove(os.path.join(d, f))
    manifests = {'visible': [], 'heldout': []}
    counters = {}
    stats = {'n_papers': len(files), 'seed': seed, 'per_split': {'visible': 0, 'heldout': 0},
             'per_type': {t: {'split': type_split[t], 'copies': 0, 'papers_applicable': 0, 'sites': 0}
                          for t in ERROR_TYPES},
             'papers_with_any_applicable': 0,
             'papers_with_any_applicable_by_split': {'visible': 0, 'heldout': 0}}
    for fn in files:
        with open(os.path.join(corpus_dir, fn)) as f:
            paper = json.load(f)
        orig_id = paper['paper_id']
        snapshot = json.dumps(paper, sort_keys=True)
        any_app = {'visible': False, 'heldout': False}
        for etype in sorted(ERROR_TYPES):
            sp = type_split[etype]
            copies, nsites = plant_paper(paper, etype, seed)
            stats['per_type'][etype]['sites'] += nsites
            if copies:
                stats['per_type'][etype]['papers_applicable'] += 1
                any_app[sp] = True
            for cp, loc, before, after in copies:
                n = counters.get((sp, orig_id), 0) + 1
                counters[(sp, orig_id)] = n
                pid = f"{orig_id}__{n}"
                cp['paper_id'] = pid
                with open(os.path.join(out_dir, sp, pid + '.json'), 'w') as f:
                    json.dump(cp, f, ensure_ascii=False, indent=1)
                manifests[sp].append({"planted_id": pid, "orig_id": orig_id, "error_type": etype,
                                      "location": loc, "before": before, "after": after})
                stats['per_type'][etype]['copies'] += 1
                stats['per_split'][sp] += 1
        assert json.dumps(paper, sort_keys=True) == snapshot, 'original mutated'
        if any_app['visible'] or any_app['heldout']:
            stats['papers_with_any_applicable'] += 1
        for sp in any_app:
            stats['papers_with_any_applicable_by_split'][sp] += any_app[sp]
    for sp, lines in manifests.items():
        with open(os.path.join(out_dir, sp, 'manifest.jsonl'), 'w') as f:
            for rec in lines:
                f.write(json.dumps(rec, ensure_ascii=False) + '\n')
    with open(os.path.join(out_dir, 'stats.json'), 'w') as f:
        json.dump(stats, f, indent=2)
    if verbose:
        print(json.dumps(stats, indent=2))
    return stats


def main():
    ap = argparse.ArgumentParser()
    root = os.path.dirname(HERE)
    ap.add_argument('--corpus', default=os.path.join(root, 'corpus', 'parsed'))
    ap.add_argument('--out', default=HERE)
    a = ap.parse_args()
    run(a.corpus, a.out)


if __name__ == '__main__':
    main()
