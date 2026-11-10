import re
import math

NUM = r'\d{1,3}(?:,\d{3})+|\d+'
TAIL = r'\s*[\*†‡§#a-d]*\s*$'
COMB = re.compile(r'^\s*(' + NUM + r')\s*(?:/\s*(' + NUM + r')\s*)?[\(\[]\s*(\d{1,3}(?:\.\d+)?)\s*(%?)\s*[\)\]]' + TAIL)
REV = re.compile(r'^\s*(\d{1,3}(?:\.\d+)?)\s*%?\s*[\(\[]\s*(' + NUM + r')\s*[\)\]]' + TAIL)
INTC = re.compile(r'^\s*(' + NUM + r')' + TAIL)
PCTC = re.compile(r'^\s*(\d{1,3}(?:\.\d+)?)\s*%?' + TAIL)
NEQ = re.compile(r'\bn\s*=\s*(' + NUM + r')', re.I)
LABEL_SKIP = re.compile(r'(\bmean\b|\bsd\b|s\.d\.|\bmedian\b|\biqr\b|\brange\b|±|\baverage\b|\bm\s*\()', re.I)
NCOL = re.compile(r'^\(?\s*(n|no\.?|number|count|counts|frequency|frequencies|freq\.?|f|cases)\s*\)?$', re.I)
TEXT_PAT = re.compile(r'\b(' + NUM + r')\s*(?:/|out of|of)\s*(' + NUM + r')\b([^\d%;]{0,40}?)\(?\s*(\d{1,3}(?:\.\d+)?)\s*%')


def _num(s):
    return int(s.replace(',', ''))


def _dec(s):
    return len(s.split('.')[1]) if '.' in s else 0


def _tol(d):
    return 0.5 * (10 ** -d) * 1.1 + 1e-9


def _consistent(n, pct, d, N):
    if N <= 0:
        return False
    return abs(n * 100.0 / N - pct) <= _tol(d)


def _interval(n, pct, d):
    t = _tol(d)
    if n == 0:
        if pct <= t:
            return (0.0, math.inf)
        return None
    if pct <= t:
        return (n * 100.0 / t, math.inf)
    lo = n * 100.0 / (pct + t)
    hi = n * 100.0 / (pct - t)
    return (lo, hi)


def _has_int(lo, hi):
    if hi == math.inf:
        return True
    return math.floor(hi) >= math.ceil(lo)


def _cell(row, c):
    return row[c] if c < len(row) else ''


def _is_numericish(s):
    return bool(re.search(r'\d', s or ''))


def _col_header(table, c):
    hr = table.get('header_rows') or []
    full = ' '.join(_cell(h, c) for h in hr).strip()
    last = _cell(hr[-1], c).strip() if hr else ''
    return full, last


def _row_label(row, upto):
    parts = []
    for c in range(min(upto, len(row))):
        s = row[c]
        if not re.fullmatch(r'[\d\s.,%()\-–*]+', s or ''):
            parts.append(s)
    return ' '.join(parts)


def _extract_columns(table):
    """Return list of column specs: dict(entries={row:(n,pct,d,explN,flagcol)}, cand=[N...])"""
    rows = table.get('rows') or []
    if not rows:
        return []
    ncols = max(len(r) for r in rows)
    cap = ((table.get('caption') or '') + ' ' + (table.get('footnotes') or '')).lower()
    cap_pct = bool(re.search(r'n\s*\(%\)|no\.?\s*\(%\)|number\s*\(%\)|\(%\)|percent|frequenc|%', cap))
    cap_N = [_num(m) for m in NEQ.findall(table.get('caption') or '')]
    specs = []
    used = set()
    # separate n / % columns
    for c in range(1, ncols):
        full, last = _col_header(table, c)
        pfull, plast = _col_header(table, c - 1)
        lastl = last.lower()
        if ('%' in lastl or 'percent' in lastl) and not re.search(r'\bci\b|confidence', lastl) \
                and NCOL.match(re.sub(r'[\s:]+$', '', plast.strip())):
            entries = {}
            for r, row in enumerate(rows):
                a, b = _cell(row, c - 1), _cell(row, c)
                ma, mb = INTC.match(a or ''), PCTC.match(b or '')
                if not (ma and mb):
                    continue
                if LABEL_SKIP.search(_row_label(row, c - 1)):
                    continue
                pct = float(mb.group(1))
                if pct > 100.0001:
                    continue
                entries[r] = (_num(ma.group(1)), pct, _dec(mb.group(1)), None, c)
            cand = [_num(m) for m in NEQ.findall(full + ' ' + pfull)] + cap_N
            if entries:
                specs.append({'entries': entries, 'cand': cand, 'firstcol': c - 1})
                used.add(c)
                used.add(c - 1)
    # combined columns
    for c in range(1, ncols):
        if c in used:
            continue
        full, last = _col_header(table, c)
        fl = full.lower()
        reversed_ = bool(re.search(r'%\s*\(\s*n\s*\)', fl))
        if re.search(r'mean|\bsd\b|median|±|\bci\b', fl) and '%' not in fl:
            continue
        cells = [(r, _cell(row, c)) for r, row in enumerate(rows)]
        nonempty = [(r, s) for r, s in cells if s.strip() and _is_numericish(s)]
        if len(nonempty) < 2:
            continue
        pat = REV if reversed_ else COMB
        matches = [(r, s, pat.match(s)) for r, s in nonempty]
        good = [(r, s, m) for r, s, m in matches if m]
        if len(good) < 2 or len(good) < 0.5 * len(nonempty):
            continue
        any_pct_cell = any('%' in s for _, s, _ in good)
        hdr_pct = ('%' in fl) or ('percent' in fl)
        if not (any_pct_cell or hdr_pct or cap_pct):
            continue
        entries = {}
        for r, s, m in good:
            if '±' in s:
                continue
            if LABEL_SKIP.search(_row_label(rows[r], c)):
                continue
            if reversed_:
                pct_s, n_s, N_s = m.group(1), m.group(2), None
            else:
                n_s, N_s, pct_s = m.group(1), m.group(2), m.group(3)
            pct = float(pct_s)
            if pct > 100.0001:
                continue
            n = _num(n_s)
            explN = _num(N_s) if N_s else None
            if explN is not None and (explN < n or explN == 0):
                continue
            entries[r] = (n, pct, _dec(pct_s), explN, c)
        cand = [_num(m) for m in NEQ.findall(full)] + cap_N
        if entries:
            specs.append({'entries': entries, 'cand': cand, 'firstcol': c})
    return specs


def _groups(table, spec):
    rows = table.get('rows') or []
    two_label = False
    for row in rows:
        if len(row) >= 3 and spec['firstcol'] >= 2:
            if not (row[0] or '').strip() and (row[1] or '').strip() and not re.fullmatch(r'[\d\s.,%()\-–*]+', row[1]):
                two_label = True
                break
    gid = {}
    block = 0
    key = ''
    for r, row in enumerate(rows):
        if r not in spec['entries']:
            block += 1
        if two_label and (row[0] or '').strip():
            key = row[0].strip() + '#' + str(r)
        gid[r] = (key, block)
    return gid


def _check_spec(table, spec, flags):
    entries = spec['entries']
    tid = table.get('table_id')
    # explicit n/N (pct)
    for r, (n, pct, d, explN, fc) in entries.items():
        if explN is not None and not _consistent(n, pct, d, explN):
            flags.append({"location": {"kind": "cell", "table_id": tid, "row": r, "col": fc},
                          "message": f"Percentage {pct} does not match {n}/{explN} = {n*100.0/explN:.2f}%."})
    S = [(r, e) for r, e in entries.items() if e[3] is None and (e[0] > 0 or e[1] > 0)]
    if len(S) < 3:
        return
    ivs = {}
    for r, e in S:
        ivs[r] = _interval(e[0], e[1], e[2])
    cands = set(x for x in spec['cand'] if x > 0)
    for r, e in S:
        iv = ivs[r]
        if iv is None or iv[1] == math.inf:
            continue
        cands.add(int(math.ceil(iv[0])))
        cands.add(int(math.floor(iv[1])))
    cands = [N for N in cands if N > 0]
    if not cands:
        return
    scores = {}
    for N in cands:
        scores[N] = sum(1 for r, e in S if _consistent(e[0], e[1], e[2], N))
    best_score = max(scores.values())
    if best_score < 3 or best_score < 0.5 * len(S):
        return
    best = sorted(N for N, s in scores.items() if s == best_score)
    if best[-1] > 1.5 * best[0] + 2:
        return  # ambiguous denominator
    gid = _groups(table, spec)

    def ok_best(e):
        return any(_consistent(e[0], e[1], e[2], N) for N in best)

    for r, e in S:
        if ok_best(e):
            continue
        others = [(r2, e2) for r2, e2 in S if r2 != r and gid.get(r2) == gid.get(r)]
        iv = ivs[r]
        flag = False
        if others:
            if not any(ok_best(e2) for _, e2 in others):
                continue
            L, H = 0.0, math.inf
            bad = False
            for r2, e2 in others:
                iv2 = ivs[r2]
                if iv2 is None:
                    bad = True
                    break
                L, H = max(L, iv2[0]), min(H, iv2[1])
            if bad:
                continue
            if iv is None:
                flag = True
            else:
                lo, hi = max(L, iv[0]), min(H, iv[1])
                flag = not (lo <= hi and _has_int(lo, hi))
        else:
            if best_score < 5 or best_score < 0.7 * len(S):
                continue
            if iv is None:
                flag = True
            else:
                lo, hi = max(0.75 * best[0], iv[0]), min(best[-1] + 0.5, iv[1])
                flag = not (lo <= hi and _has_int(lo, hi))
        if flag:
            Nb = best[0]
            flags.append({"location": {"kind": "cell", "table_id": table.get('table_id'), "row": r, "col": e[4]},
                          "message": f"Count {e[0]} with {e[1]}% is inconsistent with the denominator N={Nb} "
                                     f"shared by other rows ({e[0]}/{Nb} = {e[0]*100.0/Nb:.2f}%)."})


def _text_flags(paper):
    out = []
    for p in paper.get('paragraphs') or []:
        txt = p.get('text') or ''
        for m in TEXT_PAT.finditer(txt):
            a, b = _num(m.group(1)), _num(m.group(2))
            if b == 0 or a > b:
                continue
            if 1800 <= a <= 2100 and 1800 <= b <= 2100:
                continue
            pct_s = m.group(4)
            pct = float(pct_s)
            if pct > 100:
                continue
            if not _consistent(a, pct, _dec(pct_s), b):
                out.append({"location": {"kind": "text", "para_id": p.get('para_id'),
                                         "start": m.start(), "end": m.end()},
                            "message": f"{pct}% does not match {a}/{b} = {a*100.0/b:.2f}%."})
    return out


def applicable(paper):
    try:
        for t in paper.get('tables') or []:
            if _extract_columns(t):
                return True
        for p in paper.get('paragraphs') or []:
            if TEXT_PAT.search(p.get('text') or ''):
                return True
    except Exception:
        return False
    return False


def check(paper):
    flags = []
    try:
        for t in paper.get('tables') or []:
            try:
                for spec in _extract_columns(t):
                    _check_spec(t, spec, flags)
            except Exception:
                continue
        flags.extend(_text_flags(paper))
    except Exception:
        return []
    seen = set()
    out = []
    for f in flags:
        k = repr(sorted(f['location'].items()))
        if k in seen:
            continue
        seen.add(k)
        out.append(f)
    return out
