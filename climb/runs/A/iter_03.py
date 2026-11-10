import re, math
from scipy import stats as _st

NUM = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)'
MK = r'[*\u2020\u2021\u00a7#]+'
SEP = r'(?:,|;|\u2013|\u2014|to|-)'
INF = float('inf')

RE_EM = re.compile(rf'^({NUM})\s*%?\s*({MK})?$')
RE_PO = re.compile(rf'^[\(\[]\s*({NUM})\s*[\)\]]\s*({MK})?$')
RE_EP = re.compile(rf'^({NUM})\s*({MK})?\s*[\(\[]\s*({NUM})\s*[\)\]]\s*({MK})?$')
RE_EC = re.compile(rf'^({NUM})\s*%?\s*({MK})?\s*[\(\[]\s*({NUM})\s*%?\s*{SEP}\s*({NUM})\s*%?\s*[\)\]]\s*({MK})?$')
RE_CI = re.compile(rf'^[\(\[]?\s*({NUM})\s*%?\s*{SEP}\s*({NUM})\s*%?\s*[\)\]]?$')
RE_PV = re.compile(rf'^(?:[pP]\s*)?([<>=\u2264\u2265]?)\s*({NUM})\s*({MK})?$')


def nt(s):
    s = (s or '').replace('\u2212', '-').replace('\u2012', '-').replace('\ufe63', '-').replace('\uff0d', '-')
    s = re.sub(r'(?:^|(?<=[\(\[,;:\s=]))[\u2013\u2014](?=\.?\d)', '-', s)
    return s


def fnum(s):
    s = s.strip()
    if not re.fullmatch(NUM, s):
        return None
    v = float(s)
    d = len(s.split('.')[1]) if '.' in s else 0
    return (v, 0.5 * 10 ** (-d))


def nsf(z):
    return float(_st.norm.sf(z))


def prange_z(zlo, zhi, df_lo, tails):
    plo = 0.0 if zhi == INF else tails * nsf(zhi)
    phi = tails * (float(_st.t.sf(zlo, df_lo)) if df_lo else nsf(zlo))
    return min(plo, 1.0), min(phi, 1.0)


def parse_legend(txt):
    t = nt(txt).replace('\u2264', '<').replace('\u2a7d', '<')
    leg = {}
    # list style: "***, **, and * denote significance at 1%, 5%, and 10%"
    for m in re.finditer(r'((?:[*\u2020\u2021]+\s*(?:,\s*and|,\s*or|,|and|or|&|/|;)?\s*){2,})', t):
        syms = re.findall(r'[*\u2020\u2021]+', m.group(1))
        if len(set(syms)) != len(syms):
            continue

        def vals_of(seg):
            vals = []
            for a, b in re.findall(r'(\d+(?:\.\d+)?)\s*(?:%|percent)|(?<![\d.])(0?\.\d+)', seg):
                v = float(a) / 100 if a else float(b)
                if 0 < v <= 0.2:
                    vals.append(v)
            return vals
        after = t[m.end(): m.end() + 160]
        vals = vals_of(after) if re.search(r'signific|p\s*<|level', after, re.I) else []
        if len(vals) < len(syms):
            before = t[max(0, m.start() - 160): m.start()]
            if re.search(r'signific|level', before, re.I):
                vb = vals_of(before)
                vals = vb[-len(syms):] if len(vb) >= len(syms) else []
        if len(vals) >= len(syms):
            for s, v in zip(syms, vals):
                leg[s] = v
            return leg
    for m in re.finditer(rf'({MK})[^*\u2020\u2021\u00a7#\d]{{0,50}}?(?:[pP]\s*(?:-?values?)?\s*(?:<|=|less than)\s*|\bat(?: the)?\s*(?:\u03b1\s*=\s*)?)({NUM})\s*(%)?', t):
        try:
            v = float(m.group(2))
        except ValueError:
            continue
        if m.group(3):
            v /= 100
        if 0 < v <= 0.2 and m.group(1) not in leg:
            leg[m.group(1)] = v
    return leg


def mark_interval(leg, mark):
    thr = sorted(set(leg.values()))
    if not thr:
        return None
    if not mark:
        return (max(thr), 1.0)
    if mark not in leg:
        return None
    T = leg[mark]
    return (max([x for x in thr if x < T], default=0.0), T)


def p_interval(cell):
    m = RE_PV.match(nt(cell).strip())
    if not m:
        return None
    op, f = m.group(1), fnum(m.group(2))
    if f is None:
        return None
    v, h = f
    if not (0 <= v <= 1):
        return None
    if op in ('<', '\u2264'):
        return (0.0, v)
    if op in ('>', '\u2265'):
        return (v, 1.0)
    if v == 0:
        return (0.0, h)
    return (max(0, v - h), min(1, v + h))


def ovl(a, b):
    return a[0] <= b[1] and b[0] <= a[1]


def hdrs(table):
    hr = table.get('header_rows') or []
    ncol = max([len(r) for r in hr] + [len(r) for r in table.get('rows', [])] + [0])
    out = []
    for c in range(ncol):
        parts = []
        for r in hr:
            if c < len(r):
                t = r[c].strip()
                if t and (not parts or parts[-1] != t):
                    parts.append(t)
        out.append(nt(' '.join(parts)))
    return out


def cell(rows, r, c):
    if c is not None and 0 <= r < len(rows) and c < len(rows[r]):
        return nt(rows[r][c]).strip()
    return ''


def z_from(H, e, s):
    (b, hb), (sv, hs) = e, s
    if H == 't':
        return (max(0, abs(sv) - hs), abs(sv) + hs)
    if sv <= 0:
        return None
    slo = sv - hs
    shi = sv + hs
    if H == 'se':
        lo = max(0, abs(b) - hb) / shi
        hi = (abs(b) + hb) / slo if slo > 0 else INF
        return (lo, hi)
    if b - hb <= 0:
        return None
    pts = [b - hb + k * (2 * hb) / 6 for k in range(7)]
    if H == 'se_log':
        f = [abs(math.log(x)) for x in pts]
    else:
        f = [abs(math.log(x)) * x for x in pts]
    if (b - hb) <= 1 <= (b + hb):
        f.append(0)
    return (min(f) / shi, max(f) / slo if slo > 0 else INF)


def get_N(txt):
    Ns = []
    for m in re.finditer(r'\b[Nn]\s*=\s*(\d[\d,]*)(?:\s*(?:\u2013|-|to)\s*(\d[\d,]*))?', nt(txt)):
        for g in m.groups():
            if g:
                try:
                    v = int(g.replace(',', ''))
                    if v >= 20:
                        Ns.append(v)
                except ValueError:
                    pass
    return Ns


def ci_z(e, lo, hi, ratio, crit):
    (v, h), (l, hl), (u, hu) = e, lo, hi
    if ratio:
        if l - hl <= 0 or v <= 0:
            return None
        f = math.log
    else:
        f = lambda x: x
    try:
        V, L, U = f(v), f(l), f(u)
        Lw, Uw = f(l - hl), f(u + hu)
        Ln = f(l + hl) if l + hl < u - hu else L
        Un = f(u - hu) if u - hu > l + hl else U
    except ValueError:
        return None
    w = U - L
    if w <= 0:
        return None
    if abs((V - L) - (U - V)) > 0.3 * w + (f(v + h) - f(max(v - h, 1e-9)) if ratio else 2 * h) + (Uw - Un):
        return None
    se_hi = (Uw - Lw) / (2 * crit)
    se_lo = max((Un - Ln) / (2 * crit), 1e-12)
    mid = (L + U) / 2
    zl = max(0, abs(mid) - (Uw - Un)) / se_hi
    zh = (abs(mid) + (Uw - Un)) / se_lo
    return (zl * 0.8, zh * 1.25)


def crit_of(txt):
    m = re.search(r'(9\d(?:\.\d)?)\s*%', txt)
    if m:
        lv = float(m.group(1)) / 100
        if 0.5 < lv < 1:
            return float(_st.norm.ppf(1 - (1 - lv) / 2))
    return 1.959964


RATIO_RE = re.compile(r'odds ratio|hazard ratio|risk ratio|rate ratio|relative risk|prevalence ratio|\ba?ors?\b|\ba?hrs?\b|\birrs?\b|\brrs?\b|exp\s*\(\s*[b\u03b2]\s*\)|exponentiated', re.I)


class Acc:
    def __init__(self):
        self.groups = {}
        self.n = 0

    def add(self, key, ok, loc=None, msg=None):
        g = self.groups.setdefault(key, [0, []])
        g[0] += 1
        self.n += 1
        if not ok:
            g[1].append((loc, msg))

    def flags(self):
        out = []
        for key, (n, bad) in self.groups.items():
            if bad and len(bad) <= max(1, int(0.2 * n)) and len(bad) < max(n, 2):
                for loc, msg in bad:
                    out.append({"location": loc, "message": msg})
        return out


def cloc(t, r, c):
    return {"kind": "cell", "table_id": t['table_id'], "row": r, "col": c}


def check_pairs(t, leg, ctx, df_lo, tails, acc):
    rows = t.get('rows') or []
    hd = hdrs(t)
    items = []
    for c in range(1, len(hd)):
        if re.search(r'\bmean\b|\bsd\b|\bm\b|median|%|\bn\b|range|iqr|freq|count|\bmin|\bmax', hd[c].lower()):
            continue
        ci = []
        for r in range(len(rows)):
            x = cell(rows, r, c)
            m = RE_EP.match(x)
            if m:
                if '.' in m.group(1):
                    ci.append((r, c, fnum(m.group(1)), m.group(2) or m.group(4) or '', fnum(m.group(3))))
                continue
            m = RE_EM.match(x)
            if m and '.' in m.group(1) and '%' not in x:
                m2 = RE_PO.match(cell(rows, r + 1, c))
                if m2:
                    ci.append((r, c, fnum(m.group(1)), m.group(2) or m2.group(2) or '', fnum(m2.group(1))))
        if any(it[3] for it in ci):
            items.extend(ci)
    if len(items) < 2:
        return
    has_se = re.search(r'standard error|\bs\.?e\.?s?\b|std\.? ?err|robust', ctx)
    has_t = re.search(r'\b[tz][- ]?(?:stat|value|ratio|score)|absolute value of [tz]|\b[tz]-values?\b|parenthes[ie]s\s*(?:are|is|represent|denote|show|report|indicate)?\s*(?:the\s*)?(?:robust\s*)?[tz]\b|\b[tz][- ]statistics?\b', ctx)
    has_p = re.search(r'p[- ]?values? (?:are )?(?:in|within) (?:parenthes|bracket)', ctx)
    allpos = all(it[2][0] > 0 for it in items)
    anyneg_par = any(it[4][0] < 0 for it in items)
    ratio = bool(RATIO_RE.search(ctx + ' ' + ' '.join(hd))) and allpos
    if has_p and not has_se and not has_t:
        hyps = ['p']
    elif has_t and not has_se:
        hyps = ['t']
    elif has_se and not has_t:
        hyps = ['se'] + (['se_delta', 'se_log'] if ratio else [])
    else:
        hyps = ['se', 't'] + (['se_delta', 'se_log'] if ratio else [])
    if anyneg_par:
        hyps = [h for h in hyps if h == 't']
    best = None
    for H in hyps:
        res = []
        for (r, c, e, mk, s) in items:
            iv = mark_interval(leg, mk)
            if iv is None:
                continue
            if H == 'p':
                if not (0 <= s[0] <= 1):
                    continue
                pr = (max(0, s[0] - s[1]), s[0] + s[1])
            else:
                zr = z_from(H, e, s)
                if zr is None:
                    continue
                pr = prange_z(zr[0], zr[1], df_lo, tails)
            ok = pr[0] < iv[1] and pr[1] >= iv[0]
            res.append((ok, r, c, e, mk, s, pr))
        if not res:
            continue
        nbad = sum(1 for x in res if not x[0])
        if best is None or nbad < best[0]:
            best = (nbad, H, res)
    if not best:
        return
    nbad, H, res = best
    for ok, r, c, e, mk, s, pr in res:
        acc.add((t['table_id'], 'pairs'), ok, cloc(t, r, c),
                f"Significance marker '{mk or 'none'}' on estimate {e[0]} inconsistent with parenthesized value {s[0]} "
                f"(read as {H}; implied p in [{pr[0]:.3g}, {pr[1]:.3g}]) under the table's star legend.")
    # signed t-statistics must share the estimate's sign
    if H == 't' and anyneg_par:
        for (r, c, e, mk, s) in items:
            if abs(e[0]) > e[1] and abs(s[0]) > s[1]:
                ok = (e[0] > 0) == (s[0] > 0)
                acc.add((t['table_id'], 'sign_t'), ok, cloc(t, r, c),
                        f"Sign of estimate {e[0]} disagrees with sign of its t-statistic ({s[0]}).")


def check_est_ci_cells(t, leg, ctx, tails, acc):
    rows = t.get('rows') or []
    hd = hdrs(t)
    crit = crit_of(ctx + ' ' + ' '.join(hd))
    ratio_ctx = bool(RATIO_RE.search(ctx + ' ' + ' '.join(hd)))
    for c in range(len(hd)):
        col = []
        for r in range(len(rows)):
            m = RE_EC.match(cell(rows, r, c))
            if not m:
                continue
            e, l, u = fnum(m.group(1)), fnum(m.group(3)), fnum(m.group(4))
            if not (e and l and u) or l[0] > u[0]:
                continue
            mk = m.group(2) or m.group(5) or ''
            col.append((r, e, l, u, mk))
            ok = l[0] - l[1] - 1e-9 <= e[0] <= u[0] + u[1] + 1e-9
            acc.add((t['table_id'], 'contain'), ok, cloc(t, r, c),
                    f"Point estimate {e[0]} lies outside its own interval [{l[0]}, {u[0]}].")
        if not leg or not any(x[4] for x in col):
            continue
        for (r, e, l, u, mk) in col:
            iv = mark_interval(leg, mk)
            if iv is None or not (l[0] <= e[0] <= u[0]):
                continue
            if l[0] > 0:
                if not ratio_ctx:
                    continue
                ratio = True
            else:
                ratio = False
            zr = ci_z(e, l, u, ratio, crit)
            if zr is None:
                continue
            pr = prange_z(zr[0], zr[1], None, tails)
            ok = pr[0] < iv[1] and pr[1] >= iv[0]
            acc.add((t['table_id'], 'ci_star'), ok, cloc(t, r, c),
                    f"Significance marker '{mk or 'none'}' inconsistent with CI [{l[0]}, {u[0]}] around {e[0]} (implied p in [{pr[0]:.3g}, {pr[1]:.3g}]).")


def role_of(h):
    hl = h.lower().strip()
    if not hl:
        return None
    if len(hl) < 40 and re.search(r'(^|\s)(lower|ll|lb)(\s|$|\b)', hl) and not re.search(r'secondary|school|class|income', hl):
        return 'lo'
    if len(hl) < 40 and re.search(r'(^|\s)(upper|ul|ub)(\s|$|\b)', hl) and not re.search(r'secondary|school|class|income', hl):
        return 'hi'
    if re.search(r'\bci\b|confidence interval|conf\. ?int', hl):
        return 'ci'
    if re.search(r'(^|\s)(p|p[- ]?values?|sig\.?|significance|pr\s*\(>\s*\|?[tz]\|?\))$', hl):
        return 'p'
    if re.search(r'(^|\s|\()(s\.?\s?e\.?|std\.? ?err(or)?\.?|standard error|se\s*b|seb|robust se)\)?$', hl):
        return 'se'
    if re.search(r'(^|\s)(sd|s\.d\.|std\.? deviation|standard deviation)$', hl):
        return 'sd'
    if re.search(r'\bwald\b', hl):
        return 'wald'
    if re.search(r'(^|\s)([tz]|[tz][- ]?values?|[tz][- ]?stat(istic)?s?|z[- ]?score|t ratio)$', hl):
        return 'tz'
    if hl in ('df', 'd.f.'):
        return 'df'
    if RATIO_RE.search(hl):
        return 'ratio'
    last = hl.split()[-1]
    if last in ('b', 'coef', 'coef.', 'coefficient', 'coefficients', 'estimate', 'estimates', 'est.', 'est', '\u03b3'):
        return 'est'
    if last in ('\u03b2', 'beta'):
        return 'beta'
    return None


SKIP_ROLES = ('se', 'tz', 'wald', 'p', 'df', 'sd')


def ci_specs(t):
    """Return dict est_col -> (spec_type, c1, c2) for CI columns and the estimate col to their left."""
    rows = t.get('rows') or []
    hd = hdrs(t)
    roles = [role_of(h) for h in hd]

    def single_frac(c):
        ns = nr = 0
        for r in range(len(rows)):
            x = cell(rows, r, c)
            if RE_EM.match(x):
                ns += 1
            elif RE_CI.match(x):
                nr += 1
        return ns, nr

    specs = []
    c = 0
    while c < len(hd):
        ro = roles[c]
        if ro == 'ci':
            if c + 1 < len(hd) and roles[c + 1] == 'ci':
                a, b = single_frac(c), single_frac(c + 1)
                if a[0] > a[1] and b[0] > b[1] and a[0] >= 1:
                    specs.append(('pair', c, c + 1))
                    c += 2
                    continue
            specs.append(('one', c, None))
        elif ro == 'lo':
            for d in range(c + 1, min(c + 3, len(hd))):
                if roles[d] == 'hi':
                    specs.append(('pair', c, d))
                    break
        c += 1
    out = {}
    for sp in specs:
        d = sp[1] - 1
        while d >= 1 and roles[d] in SKIP_ROLES:
            d -= 1
        if d < 1 or roles[d] in ('ci', 'lo', 'hi'):
            continue
        nn = num = 0
        for r in range(len(rows)):
            x = cell(rows, r, d)
            if not x:
                continue
            nn += 1
            if RE_EM.match(x) or RE_EC.match(x):
                num += 1
        if nn and num >= max(1, 0.5 * nn):
            out[d] = sp
    return out


def ci_at(rows, r, sp):
    if sp[0] == 'pair':
        a = RE_EM.match(cell(rows, r, sp[1]))
        b = RE_EM.match(cell(rows, r, sp[2]))
        if a and b:
            return fnum(a.group(1)), fnum(b.group(1))
        return None, None
    m = RE_CI.match(cell(rows, r, sp[1]))
    if m:
        return fnum(m.group(1)), fnum(m.group(2))
    return None, None


def est_at(rows, r, c):
    x = cell(rows, r, c)
    m = RE_EM.match(x)
    if m:
        return fnum(m.group(1)), (m.group(2) or '')
    m = RE_EC.match(x)
    if m:
        return fnum(m.group(1)), (m.group(2) or m.group(5) or '')
    return None, ''


def check_ci_columns(t, acc):
    rows = t.get('rows') or []
    for ec, sp in ci_specs(t).items():
        for r in range(len(rows)):
            E, _ = est_at(rows, r, ec)
            L, U = ci_at(rows, r, sp)
            if not (E and L and U) or L[0] > U[0]:
                continue
            ok = L[0] - L[1] - 1e-9 <= E[0] <= U[0] + U[1] + 1e-9
            acc.add((t['table_id'], 'contain'), ok, cloc(t, r, ec),
                    f"Estimate {E[0]} lies outside its CI [{L[0]}, {U[0]}].")


def check_columns(t, leg, ctx, df_lo, tails, acc):
    rows = t.get('rows') or []
    hd = hdrs(t)
    specs = ci_specs(t)
    groups, cur = [], {}
    for c, h in enumerate(hd):
        ro = role_of(h)
        if ro is None or ro == 'sd':
            continue
        if ro in ('est', 'ratio', 'beta') and (ro in cur or ('p' in cur and ro != 'beta' and ('est' in cur or 'ratio' in cur))):
            groups.append(cur)
            cur = {}
        if ro not in cur:
            cur[ro] = c
    groups.append(cur)
    tid = t['table_id']
    crit = crit_of(ctx + ' ' + ' '.join(hd))
    tz_neg = any(cell(rows, r, g.get('tz')).startswith('-') for g in groups if 'tz' in g for r in range(len(rows)))
    for g in groups:
        if 'est' not in g and 'beta' in g:
            g['est'] = g['beta']
        ec = g.get('est')
        rc = g.get('ratio')
        if ec is None and rc is None:
            continue
        mcol = ec if ec is not None else rc
        col_has_mk = any(re.search(MK, cell(rows, rr, mcol)) for rr in range(len(rows)))
        for r in range(len(rows)):
            def num_at(c):
                if c is None:
                    return None
                m = RE_EM.match(cell(rows, r, c))
                return fnum(m.group(1)) if m else None
            E, emk = est_at(rows, r, ec) if ec is not None else (None, '')
            Q, qmk = est_at(rows, r, rc) if rc is not None else (None, '')
            S = num_at(g.get('se'))
            T = num_at(g.get('tz'))
            W = num_at(g.get('wald'))
            D = num_at(g.get('df'))
            pc = g.get('p')
            P = p_interval(cell(rows, r, pc)) if pc is not None else None
            mk = emk or qmk
            zr = None
            if E and S and S[0] > 0:
                zr = z_from('se', E, S)
            if zr and T:
                tr = (max(0, abs(T[0]) - T[1]), abs(T[0]) + T[1])
                ok = ovl((zr[0] * 0.92 - 0.05, zr[1] * 1.08 + 0.05), tr)
                acc.add((tid, 'bse_t'), ok, cloc(t, r, g['tz']),
                        f"Test statistic {T[0]} inconsistent with estimate {E[0]} / SE {S[0]} (implied {zr[0]:.2f}-{zr[1]:.2f}).")
            if E and T and tz_neg and abs(E[0]) > E[1] and abs(T[0]) > T[1]:
                acc.add((tid, 'sign_t'), (E[0] > 0) == (T[0] > 0), cloc(t, r, ec),
                        f"Sign of estimate {E[0]} disagrees with sign of its test statistic {T[0]}.")
            if zr and W and W[0] >= 0:
                wz = (math.sqrt(max(0, W[0] - W[1])), math.sqrt(W[0] + W[1]))
                ok = ovl((zr[0] * 0.88 - 0.1, zr[1] * 1.12 + 0.1), wz)
                acc.add((tid, 'bse_w'), ok, cloc(t, r, g['wald']),
                        f"Wald statistic {W[0]} inconsistent with (B/SE)^2 from B={E[0]}, SE={S[0]}.")
            pr = None
            if T:
                pr = prange_z(max(0, abs(T[0]) - T[1]), abs(T[0]) + T[1], df_lo, tails)
            elif W and W[0] >= 0:
                dfw = D[0] if D and D[0] >= 1 and D[0] == int(D[0]) else 1
                pr = (float(_st.chi2.sf(W[0] + W[1], dfw)), float(_st.chi2.sf(max(0, W[0] - W[1]), dfw)))
            elif zr:
                pr = prange_z(zr[0], zr[1], df_lo, tails)
            if pr and P is not None:
                ok = ovl((pr[0] * 0.9, min(1, pr[1] * 1.1 + 1e-4)), P)
                acc.add((tid, 'stat_p'), ok, cloc(t, r, pc),
                        f"Reported p {cell(rows, r, pc)} inconsistent with test statistic / estimate-SE (implied p in [{pr[0]:.3g}, {pr[1]:.3g}]).")
            if E and Q and Q[0] > 0 and ec is not None and rc is not None:
                ex = (math.exp(E[0] - E[1]), math.exp(E[0] + E[1]))
                ok = ovl((ex[0] * 0.98 - Q[1], ex[1] * 1.02 + Q[1]), (Q[0] - Q[1], Q[0] + Q[1]))
                acc.add((tid, 'expb'), ok, cloc(t, r, rc), f"Exp(B)={Q[0]} inconsistent with B={E[0]}.")
            # CI vs p
            for pcol, pt in ((rc, Q), (ec, E)):
                if pcol is None or pcol not in specs or not pt or P is None:
                    continue
                L, U = ci_at(rows, r, specs[pcol])
                if not (L and U) or L[0] > U[0] or not (L[0] <= pt[0] <= U[0]):
                    continue
                isr = pcol == rc
                if L[0] > 0 and not isr:
                    continue
                zc = ci_z(pt, L, U, isr, crit)
                if zc:
                    prc = prange_z(zc[0], zc[1], None, tails)
                    acc.add((tid, 'ci_p'), ovl(prc, P), cloc(t, r, pc),
                            f"Reported p {cell(rows, r, pc)} inconsistent with CI [{L[0]}, {U[0]}] around {pt[0]} (implied p in [{prc[0]:.3g}, {prc[1]:.3g}]).")
                break
            if leg and col_has_mk and (E or Q):
                iv = mark_interval(leg, mk)
                if iv is not None:
                    if P is not None:
                        ok = P[0] < iv[1] and P[1] >= iv[0]
                        acc.add((tid, 'star_p'), ok, cloc(t, r, mcol),
                                f"Significance marker '{mk or 'none'}' inconsistent with reported p {cell(rows, r, pc)}.")
                    elif pr:
                        ok = pr[0] < iv[1] and pr[1] >= iv[0]
                        acc.add((tid, 'star_stat'), ok, cloc(t, r, mcol),
                                f"Significance marker '{mk or 'none'}' inconsistent with implied p in [{pr[0]:.3g}, {pr[1]:.3g}].")


def check_corr(t, leg, ctx, tails, acc):
    rows = t.get('rows') or []
    hr = t.get('header_rows') or []
    hd = hdrs(t)
    if not leg or not re.search(r'correlat|\br\b', ctx + ' ' + ' '.join(hd), re.I):
        return
    if not re.search(r'correlat', ctx + ' ' + ' '.join(hd), re.I):
        return
    last = hr[-1] if hr else []
    vc = {}
    for c in range(1, len(hd)):
        h = (last[c] if c < len(last) else '').strip()
        if re.search(r'\bmean\b|\bsd\b|\bm\b|alpha|\u03b1|range|%|\bn\b', h.lower()):
            continue
        m = re.fullmatch(r'\(?(\d{1,2})\)?\.?', h)
        vc[c] = int(m.group(1)) if m else None
    thr = sorted(set(leg.values()), reverse=True)
    lvl = lambda mk: 0 if not mk else (thr.index(leg[mk]) + 1 if mk in leg else None)
    cells = []
    for r, row in enumerate(rows):
        m = re.match(r'\s*(\d{1,2})\s*[.)]', row[0] if row else '')
        ri = int(m.group(1)) if m else None
        for c, j in vc.items():
            x = cell(rows, r, c)
            m2 = RE_EM.match(x)
            if not m2 or '.' not in m2.group(1) or '%' in x:
                continue
            f = fnum(m2.group(1))
            if abs(f[0]) > 1:
                continue
            L = lvl(m2.group(2) or '')
            if L is None:
                continue
            if ri is None or j is None:
                tri = 'X'
            else:
                tri = 'L' if ri > j else ('U' if ri < j else 'D')
            if tri == 'D':
                continue
            cells.append((r, c, f, L, tri))
    if len(cells) < 6 or not any(x[3] for x in cells):
        return
    tid = t['table_id']
    viol = {}
    groups = {}
    for x in cells:
        groups.setdefault((x[1], x[4]), []).append(x)
    for key, g in groups.items():
        for a in g:
            for b in g:
                if abs(a[2][0]) - a[2][1] >= abs(b[2][0]) + b[2][1] + 0.02 and a[3] < b[3]:
                    viol[(a[0], a[1])] = viol.get((a[0], a[1]), 0) + 1
                    viol[(b[0], b[1])] = viol.get((b[0], b[1]), 0) + 1
    for _ in range(len(cells)):
        acc.add((tid, 'corr_mono'), True)
    if viol:
        mx = max(viol.values())
        for (r, c), v in viol.items():
            if v == mx:
                acc.add((tid, 'corr_mono'), False, cloc(t, r, c),
                        "Correlation significance markers are not monotone in |r| within this column (larger |r| marked less significant than smaller |r|).")
    Ns = get_N(ctx)
    if Ns and not re.search(r'within|multilevel|daily|level 1|level-1|person-level|bootstrap', ctx, re.I):
        Nmax, Nmin = max(Ns), min(Ns)
        for (r, c, f, L, tri) in cells:
            rv = min(abs(f[0]) + f[1], 0.9999)
            if L >= 1:
                T = thr[L - 1]
                tt = rv * math.sqrt((Nmax - 2) / (1 - rv * rv))
                p = tails * float(_st.t.sf(tt, Nmax - 2))
                ok = p < T
            else:
                rl = max(0, abs(f[0]) - f[1])
                n3 = max(Nmin // 3, 10)
                tt = rl * math.sqrt((n3 - 2) / (1 - rl * rl))
                p = tails * float(_st.t.sf(tt, n3 - 2))
                ok = p >= thr[0]
            acc.add((tid, 'corr_N'), ok, cloc(t, r, c),
                    f"Correlation {f[0]} with marker level {L} inconsistent with N={Nmin}-{Nmax} (p\u2248{p:.3g}).")


def check_text(paper, acc):
    ciRE = re.compile(rf'(?:(\d{{2}}(?:\.\d+)?)\s*%\s*(?:CI|C\.I\.|confidence intervals?|CrI|credible intervals?)|\bCI\s*(\d{{2}})\s*%)\s*(?:[:=,]|of|was|from)?\s*[\(\[]?\s*({NUM})\s*%?\s*{SEP}\s*({NUM})\s*%?\s*[\)\]]?')
    seRE = re.compile(rf'\b(b|B|\u03b2|beta|Beta|coefficient|estimate|\u03b3)\s*=\s*({NUM})\s*[,;]\s*(?:SE|S\.E\.|s\.e\.|SE\s*B|SEB|standard error)\s*=\s*({NUM})((?:\s*[,;]\s*(?:t|z|p)\b[^,;)]{{0,25}}){{0,3}})')
    for para in paper.get('paragraphs') or []:
        txt = nt(para.get('text') or '')
        pid = para.get('para_id')
        for m in ciRE.finditer(txt):
            lv = m.group(1) or m.group(2)
            try:
                lvv = float(lv) / 100
            except (TypeError, ValueError):
                continue
            if not 0.5 < lvv < 1:
                continue
            crit = float(_st.norm.ppf(1 - (1 - lvv) / 2))
            L, U = fnum(m.group(3)), fnum(m.group(4))
            if not L or not U or L[0] > U[0]:
                continue
            w0 = max(0, m.start() - 60)
            win = txt[w0:m.start()]
            me = re.search(rf'({NUM})\s*%?\s*(?:[\(\[,;]\s*)*$', win)
            if not me:
                continue
            pre = win[:me.start()]
            if re.search(r'(?:\bp|\bn|\bN|SE|SD|df|\bt|\bF|\bz|\u03c72|\bM|\bI2|I\u00b2)\s*[=<>:]\s*$', pre) or re.search(r'[\d.]$', pre):
                continue
            E = fnum(me.group(1))
            if not E:
                continue
            es = w0 + me.start(1)
            keyw = pre[-35:]
            isr_key = bool(re.search(r'\b(a?ORs?|a?HRs?|RRs?|IRRs?|a?PRs?|SIR|SMR)\b', keyw) or re.search(r'odds ratio|hazard ratio|risk ratio|rate ratio|relative risk|prevalence ratio|exp\(b\)', keyw, re.I))
            is_r = bool(re.search(r'\br\s*[=:]?\s*$|correlation', keyw))
            loc = {"kind": "text", "para_id": pid, "start": es, "end": m.end()}
            ok = L[0] - L[1] - 1e-9 <= E[0] <= U[0] + U[1] + 1e-9
            acc.add(('text', 'contain'), ok, loc, f"Estimate {E[0]} lies outside its reported CI [{L[0]}, {U[0]}].")
            if not ok or is_r:
                continue
            mp = re.match(rf'\s*[\)\]]?\s*[,;]?\s*(?:and\s*)?[pP]\s*(?:-?value)?\s*([=<>\u2264\u2265])\s*({NUM})', txt[m.end():m.end() + 40])
            if not mp:
                continue
            P = p_interval(mp.group(1) + mp.group(2))
            if P is None:
                continue
            if isr_key:
                if L[0] <= 0:
                    continue
                ratio = True
            elif L[0] > 0:
                continue
            else:
                ratio = False
            zc = ci_z(E, L, U, ratio, crit)
            if not zc:
                continue
            prc = prange_z(zc[0], zc[1], None, 2)
            loc2 = {"kind": "text", "para_id": pid, "start": es, "end": m.end() + mp.end()}
            acc.add(('text', 'ci_p'), ovl(prc, P), loc2,
                    f"Reported p inconsistent with CI [{L[0]}, {U[0]}] around {E[0]} (implied p in [{prc[0]:.3g}, {prc[1]:.3g}]).")
        for m in seRE.finditer(txt):
            E, S = fnum(m.group(2)), fnum(m.group(3))
            if not E or not S or S[0] <= 0:
                continue
            zr = z_from('se', E, S)
            tail = m.group(4) or ''
            loc = {"kind": "text", "para_id": pid, "start": m.start(), "end": m.end()}
            mt = re.search(rf'\b[tz]\s*(?:\(\s*([\d.,]+)\s*\))?\s*=\s*({NUM})', tail)
            df = None
            if mt:
                T = fnum(mt.group(2))
                if mt.group(1):
                    try:
                        df = float(mt.group(1).split(',')[-1])
                    except ValueError:
                        df = None
                if T:
                    tr = (max(0, abs(T[0]) - T[1]), abs(T[0]) + T[1])
                    acc.add(('text', 'se_t'), ovl((zr[0] * 0.92 - 0.05, zr[1] * 1.08 + 0.05), tr), loc,
                            f"t/z = {T[0]} inconsistent with estimate {E[0]} / SE {S[0]}.")
                    if T[0] != 0 and E[0] != 0 and (T[0] > 0) != (E[0] > 0) and abs(T[0]) > T[1] and abs(E[0]) > E[1]:
                        acc.add(('text', 'sign'), False, loc,
                                f"Sign of t/z ({T[0]}) disagrees with sign of estimate {E[0]}.")
            mp = re.search(rf'\bp\s*([=<>\u2264\u2265])\s*({NUM})', tail)
            if mp:
                P = p_interval(mp.group(1) + mp.group(2))
                if P is not None:
                    pr = prange_z(zr[0], zr[1], df if df and df > 0 else 20, 2)
                    acc.add(('text', 'se_p'), ovl((pr[0] * 0.9, min(1, pr[1] * 1.1 + 1e-4)), P), loc,
                            f"p inconsistent with estimate {E[0]} / SE {S[0]} (implied p in [{pr[0]:.3g}, {pr[1]:.3g}]).")


def _run(paper):
    acc = Acc()
    for t in paper.get('tables') or []:
        try:
            if 'table_id' not in t:
                continue
            ctx = nt((t.get('caption') or '') + ' ' + (t.get('footnotes') or ''))
            leg = parse_legend(ctx)
            tails = 1 if re.search(r'one[- ](?:tailed|sided)', ctx, re.I) else 2
            Ns = get_N(ctx)
            df_lo = max(5, min(Ns) - 20) if Ns else 20
            if leg:
                try:
                    check_pairs(t, leg, ctx.lower(), df_lo, tails, acc)
                except Exception:
                    pass
                try:
                    check_corr(t, leg, ctx, tails, acc)
                except Exception:
                    pass
            for fn in (lambda: check_est_ci_cells(t, leg, ctx, tails, acc),
                       lambda: check_ci_columns(t, acc),
                       lambda: check_columns(t, leg, ctx, df_lo, tails, acc)):
                try:
                    fn()
                except Exception:
                    pass
        except Exception:
            continue
    try:
        check_text(paper, acc)
    except Exception:
        pass
    return acc


def applicable(paper: dict) -> bool:
    try:
        return _run(paper).n > 0
    except Exception:
        return False


def check(paper: dict) -> list:
    try:
        flags = _run(paper).flags()
    except Exception:
        return []
    seen, out = set(), []
    for f in flags:
        loc = f['location']
        if loc is None:
            continue
        k = tuple(sorted(loc.items()))
        if k in seen:
            continue
        seen.add(k)
        out.append(f)
    return out
