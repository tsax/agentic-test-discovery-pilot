import re
import math
from scipy.stats import norm, t as tdist

MINUS = str.maketrans({'\u2212': '-', '\u2012': '-', '\ufe63': '-', '\uff0d': '-'})


def nt(s):
    return (s or '').translate(MINUS).replace('\u2009', ' ').replace('\xa0', ' ').strip()


NUMP = r'[-+]?(?:\d+\.\d*|\.\d+|\d+)'
SEPC = r'\s*(,|;|–|—|to|-|‑|‐)\s*'
SINGLE = re.compile(r'^(' + NUMP + r')\s*%?\s*(\*{0,4})\s*([†‡§#]?)\s*[a-d]?$')
PAREN = re.compile(r'^[\(\[]\s*(' + NUMP + r')\s*[\)\]]\s*(\*{0,4})\s*([†‡§#]?)$')
SAMECELL = re.compile(r'^(' + NUMP + r')\s*(\*{0,4})\s*([†‡§#]?)\s*[\(\[]\s*(' + NUMP + r')\s*[\)\]]\s*(\*{0,4})\s*([†‡§#]?)$')
RANGE = re.compile(r'^[\(\[]?\s*(' + NUMP + r')\s*%?' + SEPC + r'(' + NUMP + r')\s*%?\s*[\)\]]?$')
ESTCI = re.compile(r'^(' + NUMP + r')\s*%?\s*(\*{0,4})\s*([†‡§#]?)\s*[\(\[]\s*(' + NUMP + r')\s*%?' + SEPC + r'(' + NUMP + r')\s*%?\s*[\)\]]\s*(\*{0,4})\s*([†‡§#]?)$')

EST_H = re.compile(r'^(b|coef|coeff|coefficient|coefficients|estimate|estimates|est|unstandardized b|unstandardized coefficient|b coefficient|log odds|b \(unstandardized\))$')
BETA_H = re.compile(r'^(β|beta|std β|standardized β|standardized coefficient|standardized beta)$')
SE_H = re.compile(r'^(se|s e|std err|std error|standard error|se b|se \(b\)|se\(b\)|robust se|s e b)$')
T_H = re.compile(r'^(t|z|t value|z value|t-value|z-value|t stat|t-stat|t statistic|t-statistic|t ratio|z score|z-score)$')
P_H = re.compile(r'^(p|p value|p-value|p-values|sig|significance|pr\(>\|t\|\)|pr\(>\|z\|\)|p>\|t\||p>\|z\||p ?> ?\|t\||p ?> ?\|z\|)$')
RATIO_COL = re.compile(r'^(exp ?\( ?b ?\)|e\^?b|a?or|odds ratio|adjusted odds ratio|a?hr|hazard ratio|irr|rr)$')
RATIO_H = re.compile(r'\b(a?or|a?hr|a?rr|irr|a?pr|odds|hazard|risk ratio|rate ratio|ratio|exp ?\( ?b ?\))\b')
LIN_H = re.compile(r'^(b|β|beta|coef|coefficient|estimate|difference|mean difference|md|smd|effect|ame|marginal effect)')

L1 = {'*': .05, '**': .01, '***': .001}
L2 = {'*': .1, '**': .05, '***': .01}


def dec(s):
    s = s.lstrip('+-')
    return len(s.split('.')[1]) if '.' in s else 0


def half(d):
    return 0.5 * 10 ** (-d) + 1e-12


def hnorm(s):
    s = nt(s).lower()
    s = re.sub(r'[\s\.]+', ' ', s).strip()
    return s


def parse_legend(txt):
    txt = nt(txt)
    res = {}
    for m in re.finditer(r'(\*{1,4}|†|‡)\s*[:=]?\s*\(?\s*[pP]\s*(?:-?value)?\s*(?:<=|<|≤)\s*(0?\.\d+|\d+(?:\.\d+)?\s*%)', txt):
        v = m.group(2)
        thr = float(v.replace('%', '').strip()) / 100 if '%' in v else float(v)
        res.setdefault(m.group(1), thr)
    if not res:
        m = re.search(r'(\*+)\s*,\s*(\*+)\s*,?\s*(?:and|&)?\s*(\*+)[^.]{0,150}?(\d+(?:\.\d+)?)\s*%\s*,\s*(\d+(?:\.\d+)?)\s*%\s*,?\s*(?:and|&|or)?\s*(\d+(?:\.\d+)?)\s*%', txt)
        if m:
            for i in range(3):
                res[m.group(1 + i)] = float(m.group(4 + i)) / 100
        else:
            m = re.search(r'(\*+)\s*,\s*(\*+)\s*,?\s*(?:and|&)?\s*(\*+)[^.]{0,150}?(0?\.\d+)\s*,\s*(0?\.\d+)\s*,?\s*(?:and|&|or)?\s*(0?\.\d+)', txt)
            if m:
                for i in range(3):
                    res[m.group(1 + i)] = float(m.group(4 + i))
    if not res:
        for m in re.finditer(r'(\d+(?:\.\d+)?)\s*%\s*(?:level)?\s*[\(:]?\s*(\*{1,3})(?!\*)', txt):
            res.setdefault(m.group(2), float(m.group(1)) / 100)
    if not res:
        for m in re.finditer(r'[pP]\s*(?:<|≤)\s*(0?\.\d+)\s*[:=]?\s*(\*{1,4}|†)', txt):
            res.setdefault(m.group(2), float(m.group(1)))
    stars = sorted([k for k in res if k.startswith('*')], key=len)
    for a, b in zip(stars, stars[1:]):
        if not res[b] < res[a]:
            return {}
    for k, v in res.items():
        if not (0 < v < 0.5):
            return {}
    return res


def claim(stars, dag, legend):
    if stars:
        if stars not in legend:
            return None
        return legend[stars], legend.get(stars + '*', 0.0)
    if dag:
        if dag not in legend:
            return None
        return legend[dag], legend.get('*', 0.0)
    return 1.0, max(legend.values())


def mismatch(pmin, pmax, cl):
    up, lo = cl
    if pmin > up * 1.25:
        return 'more stars than the implied p-value supports'
    if lo > 0 and pmax < lo / 1.5:
        return 'fewer stars than the implied p-value warrants'
    return None


def zint(b, bd, v, vd, kind):
    b = abs(b)
    hb, hv = half(bd), half(vd)
    if kind == 'se':
        if v <= 0:
            return None
        zmin = max(0.0, b - hb) / (v + hv)
        zmax = (b + hb) / (v - hv) if v - hv > 0 else math.inf
    else:
        a = abs(v)
        zmin, zmax = max(0.0, a - hv), a + hv
    return zmin, zmax


def pint(zmin, zmax):
    pmin = 0.0 if zmax == math.inf else 2 * norm.sf(zmax)
    pmax = 2 * norm.sf(zmin)
    return pmin, pmax


def parse_p(s):
    s = nt(s).replace(' ', '')
    s = re.sub(r'[\*†‡§#]+$', '', s)
    m = re.match(r'^(?:p)?(<|≤|<=)(0?\.\d+)$', s, re.I)
    if m:
        return 0.0, float(m.group(2))
    m = re.match(r'^(?:p=)?(0?\.\d+|1\.0*|0|1)$', s, re.I)
    if m:
        v = float(m.group(1))
        h = half(dec(m.group(1)))
        return max(0.0, v - h), min(1.0, v + h)
    return None


def cflag(tid, r, c, msg):
    return {"location": {"kind": "cell", "table_id": tid, "row": r, "col": c}, "message": msg}


def notes_of(table):
    return nt((table.get('caption') or '') + ' ' + (table.get('footnotes') or ''))


def table_has_stars(table):
    for row in table.get('rows') or []:
        for cell in row[1:]:
            if re.search(r'\d\s*\*', cell or ''):
                return True
    return False


def headers(table, ncols):
    hr = table.get('header_rows') or []
    if hr:
        last = [hnorm(hr[-1][c]) if c < len(hr[-1]) else '' for c in range(ncols)]
    else:
        last = [''] * ncols
    full = [hnorm(' '.join(h[c] for h in hr if c < len(h))) for c in range(ncols)]
    return last, full


def table_legend(table):
    return parse_legend(table.get('footnotes') or '') or parse_legend(table.get('caption') or '')


def paper_legend(paper):
    counts = {}
    for t in paper.get('tables') or []:
        lg = table_legend(t)
        if lg:
            k = tuple(sorted(lg.items()))
            counts[k] = counts.get(k, 0) + 1
    if not counts:
        return {}
    if len(counts) > 1:
        return {}
    return dict(next(iter(counts)))


# ---------------- A: coefficient with parenthetical SE / t ----------------
def check_paren(table, flags, plegend):
    rows = table.get('rows') or []
    if not rows:
        return
    ncols = max(len(r) for r in rows)
    last, full = headers(table, ncols)
    notes = notes_of(table).lower()
    pairs = []
    for r, row in enumerate(rows):
        for c in range(1, len(row)):
            if re.search(r'mean|\bsd\b|median|\bm \(|range|iqr|%', full[c] if c < len(full) else ''):
                continue
            cell = nt(row[c])
            m = SAMECELL.match(cell)
            if m:
                pairs.append(dict(r=r, c=c, b=float(m.group(1)), bd=dec(m.group(1)), stars=m.group(2) or m.group(5),
                                  dag=m.group(3) or m.group(6), v=float(m.group(4)), vd=dec(m.group(4))))
                continue
            m1 = SINGLE.match(cell)
            if m1 and '%' not in cell and r + 1 < len(rows) and c < len(rows[r + 1]):
                m2 = PAREN.match(nt(rows[r + 1][c]))
                if m2:
                    pairs.append(dict(r=r, c=c, b=float(m1.group(1)), bd=dec(m1.group(1)), stars=m1.group(2) or m2.group(2),
                                      dag=m1.group(3) or m2.group(3), v=float(m2.group(1)), vd=dec(m2.group(1))))
    if len(pairs) < 3:
        return
    if ('standard deviation' in notes or re.search(r'\bsd\b', notes)) and 'standard error' not in notes:
        return
    tid = table.get('table_id')
    anyneg = any(p['v'] < 0 for p in pairs)
    has_t = bool(re.search(r'\b[tz][- ]?(stat|statistic|statistics|value|values|ratio|ratios|score|scores)\b', notes))
    has_se = bool(re.search(r'standard error|\bs\.\s?e\.|\bses?\b|std\.? ?err', notes))
    if anyneg:
        kinds = ['t']
    elif has_t and not has_se:
        kinds = ['t']
    elif has_se and not has_t:
        kinds = ['se']
    else:
        kinds = ['se', 't']
    if anyneg:
        evald = [p for p in pairs if p['b'] != 0 and p['v'] != 0]
        bad = [p for p in evald if (p['b'] > 0) != (p['v'] > 0)]
        if len(evald) >= 3 and len(bad) <= max(1, 0.2 * len(evald)):
            for p in bad:
                flags.append(cflag(tid, p['r'], p['c'],
                                   f"Coefficient {p['b']} and its t/z statistic ({p['v']}) have opposite signs."))
    legend = table_legend(table)
    if legend:
        legends = [legend]
    elif table_has_stars(table):
        legends = [plegend] if plegend else [L1, L2]
    else:
        return
    best = None
    for kind in kinds:
        for lg in legends:
            ne, bad = 0, []
            for p in pairs:
                cl = claim(p['stars'], p['dag'], lg)
                if cl is None:
                    continue
                zi = zint(p['b'], p['bd'], p['v'], p['vd'], kind)
                if zi is None:
                    continue
                pmin, pmax = pint(*zi)
                ne += 1
                mm = mismatch(pmin, pmax, cl)
                if mm:
                    bad.append((p, mm, pmin, pmax))
            if ne == 0:
                continue
            key = (len(bad) / ne, -ne)
            if best is None or key < best[0]:
                best = (key, ne, bad, kind)
    if best is None:
        return
    _, ne, bad, kind = best
    ncomb = len(kinds) * len(legends)
    if ne < (5 if ncomb > 1 else 3):
        return
    if len(bad) > max(1, 0.2 * ne):
        return
    for p, mm, pmin, pmax in bad:
        flags.append(cflag(tid, p['r'], p['c'],
                           f"Estimate {p['b']} with {'SE' if kind == 'se' else 't'} {p['v']} implies p in "
                           f"[{pmin:.4f}, {pmax:.4f}], but the cell shows {mm}."))


# ---------------- B: separate B / SE / t / p / β / exp(B) columns ----------------
def check_columns(table, flags, plegend):
    rows = table.get('rows') or []
    if not rows:
        return
    ncols = max(len(r) for r in rows)
    last, full = headers(table, ncols)
    tid = table.get('table_id')
    legend = table_legend(table) or plegend
    for s in range(1, ncols):
        if not SE_H.match(last[s]):
            continue
        est, est_lin = None, False
        for c in range(s - 1, -1, -1):
            if EST_H.match(last[c]):
                est, est_lin = c, True
                break
            if BETA_H.match(last[c]):
                est = c
                for c2 in range(c - 1, -1, -1):
                    if EST_H.match(last[c2]):
                        est, est_lin = c2, True
                        break
                break
        if est is None:
            continue
        tc = next((c for c in range(s + 1, min(ncols, s + 4)) if T_H.match(last[c])), None)
        pc = next((c for c in range(s + 1, min(ncols, s + 6)) if P_H.match(last[c])), None)
        if tc is None and pc is None:
            continue
        ne, bad = 0, []
        for r, row in enumerate(rows):
            def g(c):
                return nt(row[c]) if c is not None and c < len(row) else ''
            mb, ms = SINGLE.match(g(est)), SINGLE.match(g(s))
            if not (mb and ms):
                continue
            b, bd, se, sd_ = float(mb.group(1)), dec(mb.group(1)), float(ms.group(1)), dec(ms.group(1))
            if se <= 0:
                continue
            mt = SINGLE.match(g(tc)) if tc is not None else None
            tv = float(mt.group(1)) if mt else None
            td = dec(mt.group(1)) if mt else None
            problems = []
            evaluated = False
            zi_bse = zint(b, bd, se, sd_, 'se')
            if tv is not None:
                evaluated = True
                if b != 0 and tv != 0 and (b > 0) != (tv > 0):
                    problems.append(f"estimate {b} and t/z {tv} have opposite signs")
                elif est_lin:
                    zmin, zmax = zi_bse
                    a, ht = abs(tv), half(td)
                    if a + ht < zmin * 0.97 or (zmax != math.inf and a - ht > zmax * 1.03):
                        problems.append(f"t/z {tv} differs from estimate/SE = {b / se:.3f}")
            if est_lin or tv is not None:
                if tv is not None:
                    zi = (max(0.0, abs(tv) - half(td)), abs(tv) + half(td))
                else:
                    zi = zi_bse
                pmin = 0.0 if zi[1] == math.inf else 2 * norm.sf(zi[1])
                pmax = 2 * tdist.sf(zi[0], 20)
                pp = parse_p(g(pc)) if pc is not None else None
                if pp is not None:
                    evaluated = True
                    if pp[1] < pmin / 1.25 - 1e-9:
                        problems.append(f"reported p ({g(pc)}) is smaller than any p implied by the estimate (≥{pmin:.4f})")
                    elif pp[0] > pmax * 1.25 + 0.002:
                        problems.append(f"reported p ({g(pc)}) is larger than implied by the estimate (≤{pmax:.4f})")
                if legend and (mb.group(2) or mb.group(3)):
                    cl = claim(mb.group(2), mb.group(3), legend)
                    if cl:
                        evaluated = True
                        mm = mismatch(pmin, 2 * norm.sf(zi[0]), cl)
                        if mm:
                            problems.append(mm)
            if evaluated:
                ne += 1
                if problems:
                    bad.append((r, problems))
        if ne >= 2 and len(bad) <= max(1, 0.2 * ne):
            for r, pr in bad:
                flags.append(cflag(tid, r, est, "Inconsistent estimate/SE/test statistic: " + "; ".join(pr) + "."))
    for e in range(ncols):
        if not EST_H.match(last[e]):
            continue
        bc = next((c for c in range(e + 1, min(ncols, e + 5)) if BETA_H.match(last[c])), None)
        oc = next((c for c in range(e + 1, min(ncols, e + 7)) if RATIO_COL.match(last[c])), None)
        if bc is None and oc is None:
            continue
        ne, bad = 0, []
        for r, row in enumerate(rows):
            mb = SINGLE.match(nt(row[e])) if e < len(row) else None
            if not mb:
                continue
            b, hb = float(mb.group(1)), half(dec(mb.group(1)))
            probs = []
            ev = False
            if bc is not None and bc < len(row):
                mbe = SINGLE.match(nt(row[bc]))
                if mbe:
                    beta = float(mbe.group(1))
                    ev = True
                    if b != 0 and beta != 0 and (b > 0) != (beta > 0):
                        probs.append(f"B {b} and standardized β {beta} have opposite signs")
            if oc is not None and oc < len(row):
                mo = SINGLE.match(nt(row[oc]))
                if mo and abs(b) < 20:
                    orv, ho = float(mo.group(1)), half(dec(mo.group(1)))
                    if orv > 0:
                        ev = True
                        lo, hi = math.exp(b - hb), math.exp(b + hb)
                        if hi * 1.01 < orv - ho or lo / 1.01 > orv + ho:
                            probs.append(f"exp(B) = {math.exp(b):.3f} but the ratio column shows {orv}")
            if ev:
                ne += 1
                if probs:
                    bad.append((r, probs))
        if ne >= 2 and len(bad) <= max(1, 0.2 * ne):
            for r, pr in bad:
                flags.append(cflag(tid, r, e, "Inconsistent coefficient companions: " + "; ".join(pr) + "."))


# ---------------- C: confidence intervals in tables ----------------
def parse_range(a, sep, b):
    lo, hi = float(a), float(b)
    if lo > hi and sep in ('-', '‑', '‐') and hi > 0:
        hi = -hi
        if lo > hi:
            return None
    if lo > hi:
        return None
    return lo, hi


def null_for(hdr, caption):
    h = hdr or ''
    if RATIO_H.search(h):
        return 1.0
    if LIN_H.match(h):
        return 0.0
    if not h and RATIO_H.search(caption):
        return 1.0
    return None


CI_WORD = r'\bci\b|\bc i\b|confidence|credible|\bcri\b'


def check_ci(table, flags, plegend):
    rows = table.get('rows') or []
    if not rows:
        return
    ncols = max(len(r) for r in rows)
    last, full = headers(table, ncols)
    tid = table.get('table_id')
    notes = hnorm(notes_of(table))
    ci_ctx = bool(re.search(CI_WORD, notes + ' ' + ' '.join(full)))
    if not ci_ctx:
        return
    other_level = bool(re.search(r'\b(90|99)\s*%', notes + ' ' + ' '.join(full)))
    legend = table_legend(table) or plegend
    pcols = [c for c in range(ncols) if P_H.match(last[c])]
    triples = []
    caption = hnorm(table.get('caption') or '')
    for c in range(1, ncols):
        hc = full[c]
        if re.search(r'iqr|range|min|max|\bsd\b|quartile|percentile', hc):
            continue
        if not re.search(CI_WORD, hc) and re.search(r'iqr|range|quartile|\bsd\b|min|max', notes):
            continue
        nv = null_for(hc, caption)
        for r, row in enumerate(rows):
            if c >= len(row):
                continue
            m = ESTCI.match(nt(row[c]))
            if not m:
                continue
            rg = parse_range(m.group(4), m.group(5), m.group(6))
            if rg is None:
                continue
            d = min(dec(m.group(1)), dec(m.group(4)), dec(m.group(6)))
            triples.append((r, c, float(m.group(1)), rg[0], rg[1], d, m.group(2) or m.group(7), m.group(3) or m.group(8), nv, 'comb'))

    # separate CI columns: collect bounds per row first, then choose estimate column by data
    for c in range(1, ncols):
        hc = full[c]
        kind = None
        if re.search(CI_WORD, hc) and not re.search(r'lower|upper|\bll\b|\bul\b', last[c]):
            kind = 'range'
        elif re.search(r'\blower\b|\bll\b|\blcl\b|\blb\b', last[c]) and c + 1 < ncols and \
                re.search(r'\bupper\b|\bul\b|\bucl\b|\bub\b', last[c + 1]):
            kind = 'lu'
        if kind is None:
            continue
        bounds = {}
        for r, row in enumerate(rows):
            if c >= len(row):
                continue
            if kind == 'range':
                mr = RANGE.match(nt(row[c]))
                if not mr:
                    continue
                rg = parse_range(mr.group(1), mr.group(2), mr.group(3))
                if rg is None:
                    continue
                bounds[r] = (rg[0], rg[1], min(dec(mr.group(1)), dec(mr.group(3))))
            else:
                if c + 1 >= len(row):
                    continue
                ml, mu = SINGLE.match(nt(row[c])), SINGLE.match(nt(row[c + 1]))
                if not (ml and mu):
                    continue
                lo, hi = float(ml.group(1)), float(mu.group(1))
                if lo > hi:
                    continue
                bounds[r] = (lo, hi, min(dec(ml.group(1)), dec(mu.group(1))))
        if not bounds:
            continue
        best = None
        for e in range(c - 1, max(-1, c - 5), -1):
            h = last[e]
            if SE_H.match(h) or T_H.match(h) or P_H.match(h) or re.search(r'^(n|df|%|wald|χ2|chi)', h):
                continue
            nin, ntot = 0, 0
            for r, (lo, hi, bd) in bounds.items():
                row = rows[r]
                me = SINGLE.match(nt(row[e])) if e < len(row) else None
                if not me:
                    continue
                ntot += 1
                est = float(me.group(1))
                tol = half(min(dec(me.group(1)), bd))
                if lo - tol <= est <= hi + tol:
                    nin += 1
            if ntot >= 1:
                frac = nin / ntot
                if best is None or frac > best[0] + 1e-9:
                    best = (frac, e, ntot)
        if best is None or best[0] < 0.8:
            continue
        e = best[1]
        nv = null_for(full[e] if e < len(full) else '', caption)
        for r, (lo, hi, bd) in bounds.items():
            row = rows[r]
            me = SINGLE.match(nt(row[e])) if e < len(row) else None
            if not me:
                continue
            d = min(dec(me.group(1)), bd)
            triples.append((r, e, float(me.group(1)), lo, hi, d, me.group(2), me.group(3), nv, 'sep'))
    if not triples:
        return
    results = []
    for (r, c, est, lo, hi, d, stars, dag, nv, src) in triples:
        tol = half(d)
        probs = []
        if est < lo - tol or est > hi + tol:
            probs.append(f"estimate {est} lies outside its confidence interval [{lo}, {hi}]")
        if nv is not None and not (nv == 1.0 and lo <= 0) and not other_level:
            excl = lo > nv + tol or hi < nv - tol
            incl = lo < nv - tol and hi > nv + tol
            row = rows[r]
            for pc in pcols:
                if pc < len(row):
                    pp = parse_p(row[pc])
                    if pp is None:
                        continue
                    if excl and pp[0] >= 0.06:
                        probs.append(f"95% CI [{lo}, {hi}] excludes {nv:g} but p = {nt(row[pc])}")
                    if incl and pp[1] < 0.045:
                        probs.append(f"95% CI [{lo}, {hi}] includes {nv:g} but p = {nt(row[pc])}")
                    break
            if legend and (stars or dag):
                cl = claim(stars, dag, legend)
                if cl and cl[0] <= 0.05 + 1e-9 and incl:
                    probs.append(f"estimate is starred (p<{cl[0]}) but 95% CI [{lo}, {hi}] includes {nv:g}")
        results.append((r, c, probs, src))
    n = len(results)
    nbad = sum(1 for x in results if x[2])
    nsep = sum(1 for x in results if x[3] == 'sep')
    if nbad > max(1, 0.2 * n):
        return
    for r, c, probs, src in results:
        if probs:
            if src == 'sep' and nsep < 3:
                continue
            flags.append(cflag(tid, r, c, "; ".join(probs) + "."))


# ---------------- D: text ----------------
TXT_CI = re.compile(r'(?:9[05]\s*%\s*)?(?:CI|C\.I\.|confidence interval|CrI|credible interval)\s*[:=,]?\s*(?:of|was|is|from)?\s*[\(\[]?\s*(' + NUMP + r')\s*%?' + SEPC + r'(' + NUMP + r')')
TXT_RATIO_CI = re.compile(r'\b(a?OR|a?HR|a?RR|IRR|a?PR)\s*[=:]?\s*(' + NUMP + r')\s*[\(\[]\s*(' + NUMP + r')' + SEPC + r'(' + NUMP + r')\s*[\)\]]')
TXT_BRACK = re.compile(r'(?<![A-Za-z])([A-Za-zβγ]{1,12})\s*=\s*(' + NUMP + r')\s*,?\s*\[\s*(' + NUMP + r')' + SEPC + r'(' + NUMP + r')\s*\]')
E_EQ = re.compile(r'([A-Za-zβχ²]*)\s*(?:=|:)\s*(' + NUMP + r')\s*%?\s*[,;]?\s*[\(\[]?\s*$')
E_KW = re.compile(r'\b(a?OR|a?HR|a?RR|IRR|a?PR|β|coefficient|estimate|odds ratio|hazard ratio|risk ratio|mean difference|MD|SMD)\s*(?:of|was|is)?\s*(' + NUMP + r')\s*%?\s*[,;]?\s*[\(\[]?\s*$')
E_PCT = re.compile(r'(?<![\d.])(' + NUMP + r')\s*%\s*[,;]?\s*[\(\[]?\s*$')
SKIP_PREF = re.compile(r'^(p|n|df|t|f|z|se|sd|chi2|χ2|χ²|k|r2|i2)$', re.I)
P_AFTER = re.compile(r'^[^.;]{0,40}?(?<![A-Za-z])p\s*(=|<|≤|>)\s*(' + NUMP + r')')
EST_T = re.compile(r'(?<![A-Za-z])(β|b|B|γ)\s*=\s*(' + NUMP + r')')
RATIO_KW = re.compile(r'^(a?or|a?hr|a?rr|irr|a?pr|odds ratio|hazard ratio|risk ratio)$', re.I)
LIN_KW = re.compile(r'^(β|b|coefficient|estimate|mean difference|md|smd|d|g|difference)$', re.I)
EQNUM = re.compile(r'[=:]\s*(' + NUMP + r')')


def tflag(pid, s, e, msg):
    return {"location": {"kind": "text", "para_id": pid, "start": s, "end": e}, "message": msg}


def text_ci_p(txt, end, est_label, lo, hi, tol):
    nv = None
    if est_label:
        if RATIO_KW.match(est_label):
            nv = 1.0
        elif LIN_KW.match(est_label):
            nv = 0.0
    if nv is None or (nv == 1.0 and lo <= 0):
        return None
    m = P_AFTER.match(txt[end:end + 60])
    if not m:
        return None
    op, v = m.group(1), float(m.group(2))
    excl = lo > nv + tol or hi < nv - tol
    incl = lo < nv - tol and hi > nv + tol
    if op == '=':
        h = half(dec(m.group(2)))
        if excl and v - h >= 0.06:
            return f"CI [{lo}, {hi}] excludes {nv:g} but p = {m.group(2)}"
        if incl and v + h < 0.045:
            return f"CI [{lo}, {hi}] includes {nv:g} but p = {m.group(2)}"
    elif op in '<≤':
        if incl and v <= 0.045:
            return f"CI [{lo}, {hi}] includes {nv:g} but p < {m.group(2)}"
    elif op == '>':
        if excl and v >= 0.06:
            return f"CI [{lo}, {hi}] excludes {nv:g} but p > {m.group(2)}"
    return None


def other_eq_inside(txt, s, lo, hi, tol):
    pre = txt[max(0, s - 150):s]
    for mn in EQNUM.finditer(pre):
        try:
            v = float(mn.group(1))
        except ValueError:
            continue
        if lo - tol <= v <= hi + tol:
            return True
    return False


def check_text(paper, flags):
    for p in paper.get('paragraphs') or []:
        raw = p.get('text') or ''
        txt = raw.translate(MINUS)
        pid = p.get('para_id')
        covered = []
        for rx, lab_g, est_g, a_g in ((TXT_RATIO_CI, 1, 2, 3), (TXT_BRACK, 1, 2, 3)):
            for m in rx.finditer(txt):
                if any(a <= m.start() < b for a, b in covered):
                    continue
                label = m.group(lab_g)
                if rx is TXT_BRACK and SKIP_PREF.match(label):
                    continue
                rg = parse_range(m.group(a_g), m.group(a_g + 1), m.group(a_g + 2))
                if rg is None:
                    continue
                es = m.group(est_g)
                est = float(es)
                d = min(dec(es), dec(m.group(a_g)), dec(m.group(a_g + 2)))
                tol = half(d)
                lo, hi = rg
                covered.append((m.start(), m.end()))
                msgs = []
                if (est < lo - tol or est > hi + tol) and not other_eq_inside(txt, m.start(est_g), lo, hi, tol):
                    msgs.append(f"Estimate {est} lies outside its reported interval [{lo}, {hi}]")
                pm = text_ci_p(txt, m.end(), label, lo, hi, tol)
                if pm:
                    msgs.append(pm)
                if msgs:
                    flags.append(tflag(pid, m.start(est_g), m.end(), "; ".join(msgs) + "."))
        for m in TXT_CI.finditer(txt):
            if any(a <= m.start() < b for a, b in covered):
                continue
            rg = parse_range(m.group(1), m.group(2), m.group(3))
            if rg is None:
                continue
            win_start = max(0, m.start() - 80)
            win = txt[win_start:m.start()]
            est_m, gi, label = None, None, ''
            me = E_EQ.search(win)
            if me and not SKIP_PREF.match(me.group(1) or ''):
                est_m, gi, label = me, 2, me.group(1) or ''
            else:
                mk = E_KW.search(win)
                if mk:
                    est_m, gi, label = mk, 2, mk.group(1)
                else:
                    mp = E_PCT.search(win)
                    if mp and not re.search(r'9[05]\s*%\s*$', win[:mp.end()].rstrip('([ ,;')):
                        est_m, gi = mp, 1
            if est_m is None:
                continue
            es = est_m.group(gi)
            est = float(es)
            d = min(dec(es), dec(m.group(1)), dec(m.group(3)))
            tol = half(d)
            lo, hi = rg
            s = win_start + est_m.start(gi)
            msgs = []
            if (est < lo - tol or est > hi + tol) and not other_eq_inside(txt, s, lo, hi, tol):
                msgs.append(f"Estimate {est} lies outside its reported CI [{lo}, {hi}]")
            pm = text_ci_p(txt, m.end(), label, lo, hi, tol)
            if pm:
                msgs.append(pm)
            if msgs:
                flags.append(tflag(pid, s, m.end(), "; ".join(msgs) + "."))
        ests = list(EST_T.finditer(txt))
        for i, m in enumerate(ests):
            seg_end = min(len(txt), m.end() + 150)
            if i + 1 < len(ests):
                seg_end = min(seg_end, ests[i + 1].start())
            seg = txt[m.end():seg_end]
            sym = m.group(1)
            b = float(m.group(2))
            bd = dec(m.group(2))
            se_m = re.search(r'\bSE\s*=\s*(' + NUMP + r')', seg)
            t_m = re.search(r'(?<![A-Za-z])([tz])\s*(?:\(\s*(\d+(?:\.\d+)?)\s*\))?\s*=\s*(' + NUMP + r')', seg)
            p_m = re.search(r'(?<![A-Za-z])p\s*(=|<|≤|>)\s*(' + NUMP + r')', seg)
            msgs = []
            endpos = m.end()
            if t_m:
                tv = float(t_m.group(3))
                if b != 0 and tv != 0 and (b > 0) != (tv > 0):
                    msgs.append(f"estimate {b} and test statistic {tv} have opposite signs")
                    endpos = m.end() + t_m.end()
                elif se_m and sym in ('b', 'B', 'γ'):
                    se = float(se_m.group(1))
                    if se > 0:
                        zmin, zmax = zint(b, bd, se, dec(se_m.group(1)), 'se')
                        a, ht = abs(tv), half(dec(t_m.group(3)))
                        if a + ht < zmin * 0.97 or (zmax != math.inf and a - ht > zmax * 1.03):
                            msgs.append(f"test statistic {tv} differs from estimate/SE = {b / se:.3f}")
                            endpos = m.end() + t_m.end()
            elif se_m and p_m and sym in ('b', 'B', 'γ'):
                se = float(se_m.group(1))
                if se > 0:
                    zmin, zmax = zint(b, bd, se, dec(se_m.group(1)), 'se')
                    pmin = 0.0 if zmax == math.inf else 2 * norm.sf(zmax)
                    pmax = 2 * tdist.sf(zmin, 20)
                    op, v = p_m.group(1), float(p_m.group(2))
                    bad = False
                    if op == '=':
                        h = half(dec(p_m.group(2)))
                        bad = (v + h < pmin / 1.25 - 1e-9) or (v - h > pmax * 1.25 + 0.002)
                    elif op in '<≤':
                        bad = pmin > v * 1.25
                    elif op == '>':
                        bad = pmax < v / 1.25
                    if bad:
                        msgs.append(f"p {op} {p_m.group(2)} is inconsistent with estimate/SE = {b / se:.3f}")
                        endpos = m.end() + p_m.end()
            if msgs:
                flags.append(tflag(pid, m.start(), endpos, "; ".join(msgs) + "."))


def applicable(paper):
    try:
        for t in paper.get('tables') or []:
            txt = ' '.join(' '.join(r) for r in (t.get('rows') or []))
            hdr = ' '.join(' '.join(h) for h in (t.get('header_rows') or []))
            if re.search(r'\(\s*-?\d*\.\d+\s*\)', txt) or re.search(r'\bCI\b|C\.I\.|SE|confidence|Exp|β|\bB\b|Lower|Upper', hdr + ' ' + (t.get('caption') or '')):
                return True
        for p in paper.get('paragraphs') or []:
            tx = (p.get('text') or '').translate(MINUS)
            if TXT_CI.search(tx) or TXT_RATIO_CI.search(tx) or EST_T.search(tx) or TXT_BRACK.search(tx):
                return True
    except Exception:
        return False
    return False


def check(paper):
    flags = []
    try:
        plegend = paper_legend(paper)
    except Exception:
        plegend = {}
    for t in paper.get('tables') or []:
        for fn in (check_paren, check_columns, check_ci):
            try:
                fn(t, flags, plegend)
            except Exception:
                continue
    try:
        check_text(paper, flags)
    except Exception:
        pass
    seen, out = set(), []
    for f in flags:
        k = repr(sorted(f['location'].items()))
        if k in seen:
            continue
        seen.add(k)
        out.append(f)
    return out
