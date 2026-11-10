import re

INT_RE = re.compile(r'^\s*(\d{1,3}(?:,\d{3})+|\d+)\s*[*†‡§]*\s*$')
PCT_RE = re.compile(r'^\s*(\d+(?:\.\d+)?)\s*%?\s*[*†‡§]*\s*$')
COMB_RE = re.compile(r'^\s*(\d{1,3}(?:,\d{3})+|\d+)\s*[\(\[]\s*(\d+(?:\.\d+)?)\s*%?\s*[\)\]]\s*[*†‡§]*\s*$')
NEQ_RE = re.compile(r'(?<![A-Za-z])[Nn]\s*[=:]\s*(\d{1,3}(?:,\d{3})+|\d+)')
SKIP_LABEL = re.compile(r'(mean|median|\bsd\b|s\.d\.|iqr|±|range|std)', re.I)
COUNT_HDR = re.compile(r'^\s*(n|no\.?|number|count|freq\.?|frequency|cases|participants)\b', re.I)

TXT_A = re.compile(r'\b(\d{1,3}(?:,\d{3})+|\d+)\s*(?:/|of|out of)\s*(\d{1,3}(?:,\d{3})+|\d+)\s*'
                   r'(?:[A-Za-z][\w\-]*\s+){0,6}?[\(\[]\s*(\d+(?:\.\d+)?)\s*%')
TXT_B = re.compile(r'(\d+(?:\.\d+)?)\s*%\s*[\(\[]\s*(?:n\s*=\s*)?(\d{1,3}(?:,\d{3})+|\d+)\s*(?:/|of|out of)\s*'
                   r'(\d{1,3}(?:,\d{3})+|\d+)\s*[\)\]]')


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
    # truncation instead of rounding
    f = 10 ** d
    if abs(int(true * f + 1e-9) / f - p) < 1e-6:
        return True
    return False


def _hdr(table, j):
    return ' '.join(hr[j] for hr in table.get('header_rows', []) if j < len(hr))


def _pairs(table):
    """Return list of (pct_col, kind, count_col)."""
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
                if cnt >= 3 and (COUNT_HDR.search(hk) or hk.strip() == h.replace('%', '').strip()
                                 or True):
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


def _check_table(table):
    flags = []
    rows = table.get('rows', [])
    ctx = (table.get('caption', '') or '') + ' ' + (table.get('footnotes', '') or '')
    if re.search(r'weight', ctx + ' ' + ' '.join(' '.join(h) for h in table.get('header_rows', [])), re.I):
        return []
    pairs = _pairs(table)
    if not pairs:
        return []
    first_pair_col = min(min(p[0], p[2]) for p in pairs)

    def label(r):
        return ' '.join(rows[r][:max(1, first_pair_col)])

    vals = {}
    for pi, pair in enumerate(pairs):
        for r in range(len(rows)):
            if SKIP_LABEL.search(label(r)):
                continue
            v = _get(table, r, pair)
            if v is None:
                continue
            if float(v[1]) > 100.0:
                continue
            vals[(pi, r)] = v

    # row-percentage cross-tab detection
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
        # support counting
        support = {}
        nconstr = 0
        for r in rlist:
            n, ps = vals[(pi, r)]
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
            if N not in support:
                support[N] = sum(1 for r in rlist if _consistent(vals[(pi, r)][0], vals[(pi, r)][1], N))
        if not support:
            continue
        best = max(support.values())
        if best < 3 or best < 0.6 * max(nconstr, 1):
            continue
        bestNs = [N for N, s in support.items() if s == best]
        # blocks
        two_level = any(rows[r][0].strip() == '' and len(rows[r]) > 1 and rows[r][1].strip()
                        for r in rlist) and first_pair_col >= 2
        blocks = []
        cur = []
        prev = None
        for r in rlist:
            brk = prev is None or r != prev + 1 or (two_level and rows[r][0].strip() != '')
            if brk and cur:
                blocks.append(cur)
                cur = []
            cur.append(r)
            prev = r
        if cur:
            blocks.append(cur)
        bsum = {}
        for b in blocks:
            if len(b) >= 2:
                s = sum(vals[(pi, r)][0] for r in b if not re.search(r'\btotal\b', label(r), re.I))
                for r in b:
                    bsum[r] = s
        for r in rlist:
            n, ps = vals[(pi, r)]
            cands = list(bestNs) + list(explicit)
            if r in bsum:
                cands.append(bsum[r])
            if any(_consistent(n, ps, N) for N in cands):
                continue
            if n > max(bestNs) and float(ps) < 100:
                pass
            Nshow = min(bestNs, key=lambda N: abs(100.0 * n / N - float(ps)))
            flags.append({
                "location": {"kind": "cell", "table_id": table.get('table_id'), "row": r, "col": j},
                "message": "Percentage %s%% inconsistent with count %d: with the table's denominator "
                           "N=%d it should be %.2f%% (no plausible denominator fits)."
                           % (ps, n, Nshow, 100.0 * n / Nshow)})
    return flags


def _check_text(paper):
    flags = []
    for para in paper.get('paragraphs', []):
        t = para.get('text', '') or ''
        seen = []
        for m in TXT_A.finditer(t):
            x, y, ps = _int(m.group(1)), _int(m.group(2)), m.group(3)
            seen.append((m.start(), m.end()))
            if y == 0 or x > y or float(ps) > 100:
                continue
            if not _consistent(x, ps, y):
                flags.append({"location": {"kind": "text", "para_id": para.get('para_id'),
                                           "start": m.start(), "end": m.end()},
                              "message": "%d of %d is %.2f%%, but text reports %s%%."
                                         % (x, y, 100.0 * x / y, ps)})
        for m in TXT_B.finditer(t):
            if any(not (m.end() <= s or m.start() >= e) for s, e in seen):
                continue
            ps, x, y = m.group(1), _int(m.group(2)), _int(m.group(3))
            if y == 0 or x > y or float(ps) > 100:
                continue
            if not _consistent(x, ps, y):
                flags.append({"location": {"kind": "text", "para_id": para.get('para_id'),
                                           "start": m.start(), "end": m.end()},
                              "message": "%d/%d is %.2f%%, but text reports %s%%."
                                         % (x, y, 100.0 * x / y, ps)})
    return flags


def applicable(paper):
    for t in paper.get('tables', []) or []:
        try:
            if _pairs(t):
                return True
        except Exception:
            pass
    for p in paper.get('paragraphs', []) or []:
        tx = p.get('text', '') or ''
        if TXT_A.search(tx) or TXT_B.search(tx):
            return True
    return False


def check(paper):
    flags = []
    for t in paper.get('tables', []) or []:
        try:
            flags.extend(_check_table(t))
        except Exception:
            continue
    try:
        flags.extend(_check_text(paper))
    except Exception:
        pass
    return flags
