import re
from collections import Counter

NUM = r'(\d{1,3}(?:,\d{3})+|\d+)'
INT_RE = re.compile(r'^\s*' + NUM + r'\s*[*†‡§]*\s*$')
PCT_RE = re.compile(r'^\s*(\d+(?:\.\d+)?)\s*%?\s*[*†‡§]*\s*$')
PCTSIGN_RE = re.compile(r'^\s*(\d+(?:\.\d+)?)\s*%\s*[*†‡§]*\s*$')
COMB_RE = re.compile(r'^\s*' + NUM + r'\s*[\(\[]\s*(\d+(?:\.\d+)?)\s*%?\s*[\)\]]\s*[*†‡§]*\s*$')
NEQ_RE = re.compile(r'(?<![A-Za-z])[Nn]\s*[=:]\s*' + NUM)
SKIP_LABEL = re.compile(r'(mean|median|\bsd\b|s\.d\.|iqr|±|range|std)', re.I)
OBS_LABEL = re.compile(r'^\s*(observations|obs\.?|n|number of \w+|sample size|total|no\. of \w+)\s*$', re.I)

TXT_A = re.compile(r'(?<![\d.,])' + NUM + r'\s*(?:/|of|out of)\s*' + NUM + r'\s*'
                   r'(?:[A-Za-z][\w\-]*\s+){0,6}?[\(\[]\s*(\d+(?:\.\d+)?)\s*%')
TXT_B = re.compile(r'(\d+(?:\.\d+)?)\s*%\s*[\(\[]\s*(?:n\s*=\s*)?' + NUM + r'\s*(?:/|of|out of)\s*'
                   + NUM + r'\s*[\)\]]')
TP1 = re.compile(r'\b[nN]\s*=\s*' + NUM + r'\s*[,;]?\s*\(?\s*(\d+(?:\.\d+)?)\s*%')
TP2 = re.compile(r'(?<![\d.,])' + NUM + r'\s*\(\s*(\d+(?:\.\d+)?)\s*%\s*\)')

POP = (r'respondents|participants|subjects|teachers|students|patients|individuals|observations|samples|'
       r'firms|companies|households|children|adults|people|persons|cases|women|men|employees|workers|'
       r'nurses|physicians|residents|users|villagers|farmers|enterprises|schools|parents|employees|'
       r'interviewees|volunteers|families|members|citizens|voters|consumers|customers|managers')
MENT1 = re.compile(r'\b[Nn]\s*=\s*' + NUM)
MENT2 = re.compile(r'(?<![\d.,])' + NUM + r'\s+(?:[A-Za-z\-]+\s+){0,2}?(?:' + POP + r')\b', re.I)
MENT3 = re.compile(r'(?:total of|sample of|sample size of|sample size was|sample consisted of|'
                   r'sample comprised)\s+' + NUM, re.I)


def _int(s):
    return int(s.replace(',', ''))


def _dec(s):
    return len(s.split('.')[1]) if '.' in s else 0


def _consistent(n, pstr, N):
    if N <= 0:
        return False
    p = float(pstr)
    d = _dec(pstr)
    true = 100.0 * n / N
    tol = 0.5 * 10 ** (-d) + 1e-6
    if abs(true - p) <= tol:
        return True
    f = 10 ** d
    if abs(int(true * f + 1e-9) / f - p) < 1e-6:
        return True
    return False


def _mindev(n, pstr, cands):
    devs = [abs(100.0 * n / N - float(pstr)) for N in cands if N > 0]
    return min(devs) if devs else 999.0


def _infer(pairs, explicit, min_support, frac):
    support = {}
    nconstr = 0
    for n, ps in pairs:
        p = float(ps)
        tol = 0.5 * 10 ** (-_dec(ps))
        if n == 0 or p - tol <= 0:
            continue
        lo = int(100.0 * n / (p + tol)) - 1
        hi = int(100.0 * n / (p - tol)) + 1
        if hi - lo > 5000:
            continue
        nconstr += 1
        for N in range(max(1, lo), hi + 1):
            if _consistent(n, ps, N):
                support[N] = support.get(N, 0) + 1
    for N in explicit:
        if N > 0 and N not in support:
            s = sum(1 for n, ps in pairs if _consistent(n, ps, N))
            if s:
                support[N] = s
    if not support:
        return None
    best = max(support.values())
    if best < min_support or best < frac * max(nconstr, 1):
        return None
    return [N for N, s in support.items() if s == best]


def _hdr(table, j):
    return ' '.join(hr[j] for hr in table.get('header_rows', []) if j < len(hr))


def _pairs(table):
    rows = table.get('rows', [])
    if not rows:
        return []
    ncol = max(len(r) for r in rows)
    out = []
    for j in range(ncol):
        h = _hdr(table, j)
        if re.search(r'mean|\bsd\b|m\s*\(sd\)|median', h, re.I):
            continue
        comb = sum(1 for r in rows if j < len(r) and COMB_RE.match(r[j]))
        if comb >= 3:
            out.append((j, 'comb', j))
            continue
        if '%' in h and not re.search(r'valid|cum', h, re.I):
            for k in (j - 1, j + 1):
                if k < 0 or k >= ncol:
                    continue
                hk = _hdr(table, k)
                if '%' in hk or re.search(r'mean|\bsd\b|median', hk, re.I):
                    continue
                cnt = sum(1 for r in rows if j < len(r) and k < len(r)
                          and INT_RE.match(r[k]) and PCT_RE.match(r[j]))
                if cnt >= 3:
                    out.append((j, 'sep', k))
                    break
    return out


def _get(table, r, pair):
    j, kind, k = pair
    row = table['rows'][r]
    if j >= len(row) or k >= len(row):
        return None
    if kind == 'comb':
        m = COMB_RE.match(row[j])
        if not m:
            return None
        return _int(m.group(1)), m.group(2)
    mi = INT_RE.match(row[k])
    mp = PCT_RE.match(row[j])
    if not mi or not mp:
        return None
    return _int(mi.group(1)), mp.group(1)


def _blocks(rows, rlist, two_level):
    blocks, cur, prev = [], [], None
    for r in rlist:
        brk = prev is None or r != prev + 1 or (two_level and rows[r][0].strip() != '')
        if brk and cur:
            blocks.append(cur)
            cur = []
        cur.append(r)
        prev = r
    if cur:
        blocks.append(cur)
    return blocks


def _is_weighted(table):
    ctx = (table.get('caption', '') or '') + ' ' + (table.get('footnotes', '') or '') + ' ' + \
        ' '.join(' '.join(h) for h in table.get('header_rows', []))
    return bool(re.search(r'weight', ctx, re.I))


# ---------------- count vs percentage in tables ----------------
def _count_pct_table(table, inferred_sets):
    flags = []
    rows = table.get('rows', [])
    if _is_weighted(table):
        return []
    ctx = (table.get('caption', '') or '') + ' ' + (table.get('footnotes', '') or '')
    pairs = _pairs(table)
    if not pairs:
        return []
    fpc = min(min(p[0], p[2]) for p in pairs)

    def label(r):
        return ' '.join(rows[r][:max(1, fpc)])

    vals = {}
    for pi, pair in enumerate(pairs):
        for r in range(len(rows)):
            if SKIP_LABEL.search(label(r)):
                continue
            v = _get(table, r, pair)
            if v is None or float(v[1]) > 100.0:
                continue
            vals[(pi, r)] = v

    if len(pairs) >= 2:
        nrow, nrp = 0, 0
        for r in range(len(rows)):
            ps = [float(vals[(pi, r)][1]) for pi in range(len(pairs)) if (pi, r) in vals]
            if len(ps) >= 2:
                nrow += 1
                s = sum(p for p in ps if abs(p - 100) > 1e-9)
                if abs(s - 100) <= 1.5:
                    nrp += 1
        if nrow and nrp >= 0.5 * nrow:
            return []

    gctx_ns = set(_int(x) for x in NEQ_RE.findall(ctx))
    for pi, pair in enumerate(pairs):
        j, kind, k = pair
        rlist = [r for r in range(len(rows)) if (pi, r) in vals]
        if len(rlist) < 3:
            continue
        explicit = set(gctx_ns)
        explicit |= set(_int(x) for x in NEQ_RE.findall(_hdr(table, j) + ' ' + _hdr(table, k)))
        for r in rlist:
            if re.search(r'\btotal\b', label(r), re.I):
                explicit.add(vals[(pi, r)][0])
        bestNs = _infer([vals[(pi, r)] for r in rlist], explicit, 3, 0.6)
        if not bestNs:
            continue
        two_level = fpc >= 2 and any(rows[r][0].strip() == '' and rows[r][1].strip() for r in rlist)
        bsum = {}
        sums = Counter()
        for b in _blocks(rows, rlist, two_level):
            if len(b) >= 2:
                s = sum(vals[(pi, r)][0] for r in b if not re.search(r'\btotal\b', label(r), re.I))
                sums[s] += 1
                for r in b:
                    bsum[r] = s
        if len(bestNs) <= 5:
            inferred_sets.append(set(bestNs))
        for s, c in sums.items():
            if c >= 2:
                inferred_sets.append({s})
        for r in rlist:
            n, ps = vals[(pi, r)]
            cands = list(bestNs) + list(explicit)
            if r in bsum:
                cands.append(bsum[r])
            if any(_consistent(n, ps, N) for N in cands):
                continue
            if _mindev(n, ps, cands) <= 0.3:
                continue
            Nshow = min(bestNs, key=lambda N: abs(100.0 * n / N - float(ps)))
            flags.append({
                "location": {"kind": "cell", "table_id": table.get('table_id'), "row": r, "col": j},
                "message": "Percentage %s%% inconsistent with count %d: with the table's denominator "
                           "N=%d it should be %.2f%%." % (ps, n, Nshow, 100.0 * n / Nshow)})
    return flags


# ---------------- percentages that should sum to 100 ----------------
def _pct_sum_table(table):
    flags = []
    rows = table.get('rows', [])
    if not rows:
        return []
    ncol = max(len(r) for r in rows)
    paired_cols = set()
    for j, kind, k in _pairs(table):
        paired_cols.add(j)
        paired_cols.add(k)
    hdrs = [_hdr(table, j) for j in range(ncol)]

    def cellpct(r, j):
        row = rows[r]
        if j >= len(row) or j == 0:
            return None
        if re.search(r'total', hdrs[j], re.I):
            return None
        c = row[j]
        m = PCTSIGN_RE.match(c)
        if not m and '%' in hdrs[j] and j not in paired_cols and \
                not re.search(r'mean|\bsd\b|ci|change|valid|cum', hdrs[j], re.I):
            m = PCT_RE.match(c)
        if not m:
            return None
        s = m.group(1)
        v = float(s)
        if v > 100 or abs(v - 100) < 1e-9:
            return None
        return s

    # row sums
    rowinfo = []
    for r in range(len(rows)):
        if SKIP_LABEL.search(rows[r][0] if rows[r] else ''):
            continue
        cells = [(j, cellpct(r, j)) for j in range(ncol)]
        cells = [(j, s) for j, s in cells if s is not None]
        if len(cells) >= 2:
            d = max(_dec(s) for _, s in cells)
            tol = len(cells) * 0.5 * 10 ** (-d) + 0.05
            tot = sum(float(s) for _, s in cells)
            rowinfo.append((r, cells, tot, tol))
    if len(rowinfo) >= 3:
        good = [x for x in rowinfo if abs(x[2] - 100) <= x[3]]
        if len(good) >= 2 and len(good) >= 0.6 * len(rowinfo):
            for r, cells, tot, tol in rowinfo:
                if abs(tot - 100) > max(tol, 1.0):
                    flags.append({"location": {"kind": "cell", "table_id": table.get('table_id'),
                                               "row": r, "col": cells[0][0]},
                                  "message": "Percentages in this row sum to %.2f, whereas the other rows "
                                             "of the table sum to 100." % tot})
    # column-block sums (pure percentage columns only)
    two_level = ncol > 2 and any(r and r[0].strip() == '' and len(r) > 1 and r[1].strip() for r in rows)
    for j in range(1, ncol):
        if j in paired_cols:
            continue
        rlist = [r for r in range(len(rows)) if cellpct(r, j) is not None
                 and not SKIP_LABEL.search(' '.join(rows[r][:j]))]
        if len(rlist) < 4:
            continue
        info = []
        for b in _blocks(rows, rlist, two_level):
            if len(b) < 2:
                continue
            ss = [cellpct(r, j) for r in b]
            d = max(_dec(s) for s in ss)
            tol = len(ss) * 0.5 * 10 ** (-d) + 0.05
            info.append((b, sum(float(s) for s in ss), tol))
        if len(info) < 3:
            continue
        good = [x for x in info if abs(x[1] - 100) <= x[2]]
        if len(good) >= 2 and len(good) >= 0.6 * len(info):
            for b, tot, tol in info:
                if abs(tot - 100) > max(tol, 1.0):
                    for r in b:
                        flags.append({"location": {"kind": "cell", "table_id": table.get('table_id'),
                                                   "row": r, "col": j},
                                      "message": "Category percentages in this block sum to %.2f, whereas "
                                                 "other blocks in the column sum to 100." % tot})
    return flags


# ---------------- sample size mentions ----------------
def _text_mentions(paper):
    out = []
    for para in paper.get('paragraphs', []) or []:
        t = para.get('text', '') or ''
        seen = set()
        for rx, strict in ((MENT1, False), (MENT2, True), (MENT3, False)):
            for m in rx.finditer(t):
                s, e = m.span(1)
                if (s, e) in seen:
                    continue
                raw = m.group(1)
                v = _int(raw)
                if v < 100:
                    continue
                if strict and ',' not in raw and 1900 <= v <= 2035:
                    continue
                seen.add((s, e))
                out.append((v, para.get('para_id'), s, e))
    return out


def _table_mentions(paper):
    vals = []
    for t in paper.get('tables', []) or []:
        ctx = (t.get('caption', '') or '') + ' ' + (t.get('footnotes', '') or '') + ' ' + \
            ' '.join(' '.join(h) for h in t.get('header_rows', []))
        for x in NEQ_RE.findall(ctx):
            v = _int(x)
            if v >= 100:
                vals.append(v)
        for row in t.get('rows', []):
            if row and OBS_LABEL.match(row[0] or ''):
                for c in row[1:]:
                    m = INT_RE.match(c) or COMB_RE.match(c)
                    if m:
                        v = _int(m.group(1))
                        if v >= 100:
                            vals.append(v)
    return vals


def _close(a, b):
    sa, sb = str(a), str(b)
    if a == b or len(sa) != len(sb) or len(sa) < 3:
        return False
    diff = [i for i in range(len(sa)) if sa[i] != sb[i]]
    if len(diff) == 1:
        return True
    if len(diff) == 2 and diff[1] == diff[0] + 1 and sa[diff[0]] == sb[diff[1]] and sa[diff[1]] == sb[diff[0]]:
        return True
    return False


def _sample_size_flags(paper, inferred_sets):
    flags = []
    tm = _text_mentions(paper)
    pool = [m[0] for m in tm] + _table_mentions(paper)
    cnt = Counter(pool)
    for v, pid, s, e in tm:
        if cnt[v] != 1:
            continue
        partner = None
        for y in cnt:
            if _close(v, y):
                partner = y
                break
        if partner is None:
            for S in inferred_sets:
                if v not in S:
                    for y in S:
                        if _close(v, y):
                            partner = y
                            break
                if partner is not None:
                    break
        if partner is not None:
            flags.append({"location": {"kind": "text", "para_id": pid, "start": s, "end": e},
                          "message": "Sample size %d is reported only here, while the paper elsewhere reports "
                                     "or implies N=%d (differs by one digit)." % (v, partner)})
    return flags


# ---------------- percentages in text ----------------
def _text_pct_flags(paper, inferred_sets, pool_vals):
    flags = []
    pair_flags = []
    for para in paper.get('paragraphs', []) or []:
        t = para.get('text', '') or ''
        pid = para.get('para_id')
        spans = []
        for m in TXT_A.finditer(t):
            x, y, ps = _int(m.group(1)), _int(m.group(2)), m.group(3)
            spans.append(m.span())
            if y == 0 or x > y or float(ps) > 100:
                continue
            if not _consistent(x, ps, y) and _mindev(x, ps, [y]) > 0.3:
                flags.append({"location": {"kind": "text", "para_id": pid, "start": m.start(), "end": m.end()},
                              "message": "%d of %d is %.2f%%, but text reports %s%%." % (x, y, 100.0 * x / y, ps)})
        for m in TXT_B.finditer(t):
            if any(not (m.end() <= s or m.start() >= e) for s, e in spans):
                continue
            spans.append(m.span())
            ps, x, y = m.group(1), _int(m.group(2)), _int(m.group(3))
            if y == 0 or x > y or float(ps) > 100:
                continue
            if not _consistent(x, ps, y) and _mindev(x, ps, [y]) > 0.3:
                flags.append({"location": {"kind": "text", "para_id": pid, "start": m.start(), "end": m.end()},
                              "message": "%d/%d is %.2f%%, but text reports %s%%." % (x, y, 100.0 * x / y, ps)})
        # n (p%) pairs with implicit denominator
        prs = []
        for rx in (TP1, TP2):
            for m in rx.finditer(t):
                if any(not (m.end() <= s or m.start() >= e) for s, e in spans):
                    continue
                if any(not (m.end() <= s or m.start() >= e) for s, e, _, _ in prs):
                    continue
                n, ps = _int(m.group(1)), m.group(2)
                if float(ps) > 100 or float(ps) <= 0:
                    continue
                prs.append((m.start(), m.end(), n, ps))
        if len(prs) >= 3:
            bestNs = _infer([(n, ps) for _, _, n, ps in prs], set(), max(2, len(prs) - 1), 0.75)
            if bestNs:
                if len(bestNs) <= 5:
                    inferred_sets.append(set(bestNs))
                cands = list(bestNs) + list(pool_vals)
                for s, e, n, ps in prs:
                    if any(_consistent(n, ps, N) for N in cands):
                        continue
                    if _mindev(n, ps, bestNs) <= 0.3:
                        continue
                    Nshow = min(bestNs, key=lambda N: abs(100.0 * n / N - float(ps)))
                    pair_flags.append({"location": {"kind": "text", "para_id": pid, "start": s, "end": e},
                                       "message": "n=%d with %s%% is inconsistent with the denominator N=%d "
                                                  "implied by the other figures in this paragraph (%.2f%%)."
                                                  % (n, ps, Nshow, 100.0 * n / Nshow)})
    if len(pair_flags) <= 2:
        flags.extend(pair_flags)
    return flags


def applicable(paper):
    for t in paper.get('tables', []) or []:
        try:
            if _pairs(t):
                return True
            for row in t.get('rows', []):
                if any(PCTSIGN_RE.match(c or '') for c in row):
                    return True
        except Exception:
            pass
    for p in paper.get('paragraphs', []) or []:
        tx = p.get('text', '') or ''
        if TXT_A.search(tx) or TXT_B.search(tx) or TP1.search(tx) or MENT1.search(tx) or MENT2.search(tx):
            return True
    return False


def check(paper):
    flags = []
    inferred_sets = []
    cp = []
    for t in paper.get('tables', []) or []:
        try:
            cp.extend(_count_pct_table(t, inferred_sets))
        except Exception:
            pass
    if len(cp) <= 2:
        flags.extend(cp)
    for t in paper.get('tables', []) or []:
        try:
            flags.extend(_pct_sum_table(t))
        except Exception:
            pass
    try:
        pool_vals = set(m[0] for m in _text_mentions(paper)) | set(_table_mentions(paper))
    except Exception:
        pool_vals = set()
    try:
        flags.extend(_text_pct_flags(paper, inferred_sets, pool_vals))
    except Exception:
        pass
    try:
        flags.extend(_sample_size_flags(paper, inferred_sets))
    except Exception:
        pass
    return flags
