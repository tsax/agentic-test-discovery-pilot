#!/usr/bin/env python3
"""Tests for the frozen planter on a synthetic paper.  Run: python3 -I planted/test_plant.py"""
import copy
import hashlib
import json
import math
import os
import re
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plant  # noqa: E402

SYN = {
    "paper_id": "synthetic_0001", "doi": "synthetic/0001", "title": "Synthetic turnout study", "year": 2022,
    "data_availability": "All data are available.",
    "paragraphs": [
        {"para_id": "p0", "section": "Methods",
         "text": "Our sample includes 1,842 respondents interviewed after the 2016 election."},
        {"para_id": "p1", "section": "Results",
         "text": "Income is positively associated with turnout (b = 0.234, SE = 0.041). Education has a "
                 "negative coefficient of −0.112 in Model 1."},
        {"para_id": "p2", "section": "Results",
         "text": "The odds of voting were higher among union members (OR = 1.45, 95% CI [1.12, 1.88])."},
        {"para_id": "p3", "section": "Results",
         "text": "A difference-in-means test confirms this, t(120) = 2.50, p = .014, for the main sample."},
    ],
    "tables": [
        {"table_id": "t1", "label": "Table 1", "caption": "OLS models of turnout.",
         "footnotes": "Standard errors in parentheses. * p<0.05, ** p<0.01, *** p<0.001",
         "header_rows": [["", "Model 1", "Model 2"]],
         "rows": [["Income", "0.234***", "0.198**"], ["", "(0.041)", "(0.065)"],
                  ["Education", "−0.112*", "−0.087"], ["", "(0.050)", "(0.061)"],
                  ["Age", "0.015", "0.021*"], ["", "(0.012)", "(0.009)"],
                  ["Constant", "1.502***", "1.344***"], ["", "(0.210)", "(0.198)"],
                  ["Observations", "1,842", "1,842"], ["R²", "0.12", "0.15"]]},
        {"table_id": "t2", "label": "Table 2", "caption": "Descriptive statistics.", "footnotes": "",
         "header_rows": [["Variable", "Mean", "SD", "Min", "Max"]],
         "rows": [["Age", "45.3", "12.1", "18", "90"], ["Income (1-10)", "5.42", "2.10", "1", "10"],
                  ["Female", "0.52", "0.50", "0", "1"]]},
        {"table_id": "t3", "label": "Table 3", "caption": "Sample characteristics, n (%).", "footnotes": "",
         "header_rows": [["", "Treatment (N = 200)", "Control (N = 200)"]],
         "rows": [["Gender", "", ""], ["Male", "96 (48.0)", "102 (51.0)"], ["Female", "100 (50.0)", "94 (47.0)"],
                  ["Other", "4 (2.0)", "4 (2.0)"]]},
        {"table_id": "t4", "label": "Table 4", "caption": "Vote choice by party identification (%).",
         "footnotes": "",
         "header_rows": [["Party", "Democrat", "Republican", "Other", "Total"]],
         "rows": [["Democrats", "81.2", "12.5", "6.3", "100.0"], ["Republicans", "9.8", "84.1", "6.1", "100.0"],
                  ["Independents", "38.4", "35.2", "26.4", "100.0"]]},
        {"table_id": "t5", "label": "Table 5", "caption": "Logistic regression of turnout (odds ratios).",
         "footnotes": "", "header_rows": [["", "OR", "95% CI"]],
         "rows": [["Member", "1.45", "1.12–1.88"], ["Union", "0.87", "0.70–1.08"]]},
    ],
}

ALL_TYPES = sorted(plant.ERROR_TYPES)


def f(s):
    return plant.to_float(s)


def cell(paper, loc):
    t = next(t for t in paper['tables'] if t['table_id'] == loc['table_id'])
    return t['rows'][loc['row']][loc['col']]


def table(paper, tid):
    return next(t for t in paper['tables'] if t['table_id'] == tid)


def diff_locations(a, b):
    """Every cell / paragraph that differs between two papers."""
    out = []
    for ta, tb in zip(a['tables'], b['tables']):
        for r, (ra, rb) in enumerate(zip(ta['rows'], tb['rows'])):
            for c, (x, y) in enumerate(zip(ra, rb)):
                if x != y:
                    out.append(('cell', ta['table_id'], r, c))
        assert ta['header_rows'] == tb['header_rows'] and ta['caption'] == tb['caption']
    for pa, pb in zip(a['paragraphs'], b['paragraphs']):
        if pa['text'] != pb['text']:
            out.append(('text', pa['para_id']))
    return out


def stars_expected(b, se, thresholds=(0.05, 0.01, 0.001)):
    z = abs(b) / se
    p = math.erfc(z / math.sqrt(2))
    return sum(p < t for t in thresholds)


def plant_all(etype, paper=SYN):
    copies, _ = plant.plant_paper(paper, etype, 12345)
    return copies


class TestDistributions(unittest.TestCase):
    def test_critical_values(self):
        self.assertAlmostEqual(plant.t_p2(2.228, 10), 0.05, places=3)
        self.assertAlmostEqual(plant.chi2_p(3.841, 1), 0.05, places=3)
        self.assertAlmostEqual(plant.chi2_p(18.307, 10), 0.05, places=3)
        self.assertAlmostEqual(plant.f_p(4.965, 1, 10), 0.05, places=3)
        self.assertAlmostEqual(plant.norm_p2(1.95996), 0.05, places=4)
        self.assertAlmostEqual(plant.t_p2(2.50, 120), 0.0138, places=3)

    def test_num_variants(self):
        v = dict((n, k) for k, n in plant.num_variants('0.021'))
        self.assertIn('0.21', v)
        self.assertIn('0.0021', v)
        self.assertIn('0.012', v)
        self.assertTrue(all(not n.startswith('00') for n in v))
        v = [n for _, n in plant.num_variants('1,842', kinds=('transpose',))]
        self.assertIn('8,142', v)
        self.assertIn('1,824', v)
        v = [n for _, n in plant.num_variants('−0.112')]
        self.assertTrue(all(n.startswith('−') for n in v))

    def test_legend(self):
        self.assertEqual(plant.parse_legend('* p<0.05, ** p<0.01, *** p<0.001'),
                         [('*', 0.05), ('**', 0.01), ('***', 0.001)])
        self.assertEqual(plant.parse_legend('***p < .01; **p < .05; *p < .1'),
                         [('*', 0.1), ('**', 0.05), ('***', 0.01)])
        self.assertEqual(plant.parse_legend('*, **, and *** denote significance at the 10%, 5%, and 1% levels.'),
                         [('*', 0.1), ('**', 0.05), ('***', 0.01)])
        self.assertIsNone(plant.parse_legend('Standard errors in parentheses.'))
        self.assertEqual(plant.parse_legend('*** (**) ((*)) indicate significance at the 1% (5%) ((10%)) level.'),
                         [('*', 0.1), ('**', 0.05), ('***', 0.01)])
        self.assertEqual(plant.parse_legend('*,**, and *** represent significance at 10, 5, and 1 percent.'),
                         [('*', 0.1), ('**', 0.05), ('***', 0.01)])
        self.assertEqual(plant.parse_legend('*, * *, * * * indicate the significance levels of 10%, 5%, and 1%'),
                         [('*', 0.1), ('**', 0.05), ('***', 0.01)])
        self.assertIsNone(plant.parse_legend('*** and * indicate significance at the 10% and 1% levels.'))


class TestPreconditions(unittest.TestCase):
    """Each type applies to the synthetic paper and does NOT apply once its precondition is removed."""

    def strip(self, mutate):
        p = copy.deepcopy(SYN)
        mutate(p)
        return p

    def assert_pre(self, etype, mutate):
        finder = plant.ERROR_TYPES[etype][0]
        self.assertTrue(finder(SYN), f'{etype} should apply')
        self.assertEqual(finder(self.strip(mutate)), [], f'{etype} should not apply after stripping')

    def test_star_flip(self):
        self.assert_pre('star_flip', lambda p: p['tables'][0].update(footnotes='Standard errors in parentheses.'))

    def test_star_flip_tstats(self):
        self.assert_pre('star_flip', lambda p: p['tables'][0].update(
            footnotes='t statistics in parentheses. * p<0.05, ** p<0.01, *** p<0.001'))

    def test_star_flip_negative_paren(self):
        self.assert_pre('star_flip', lambda p: p['tables'][0]['rows'][5].__setitem__(1, '(−0.012)'))

    def test_star_legend_inherited(self):
        p = copy.deepcopy(SYN)
        t = copy.deepcopy(p['tables'][0])
        t.update(table_id='t9', footnotes='See note to Table 1.')
        p['tables'].append(t)
        self.assertIn('t9', {s['table_id'] for s in plant.sites_star_flip(p)})

    def test_se_misreport(self):
        self.assert_pre('se_misreport', lambda p: p['tables'][0].update(footnotes=''))

    def test_sign_flip(self):
        def m(p):
            p['paragraphs'][1]['text'] = 'Income and education matter.'
        self.assert_pre('sign_flip', m)

    def test_n_text(self):
        def m(p):
            p['tables'][0]['rows'][8][0] = 'Fixed effects'
        sites = plant.sites_n_text_mismatch(SYN)
        self.assertEqual([(s['para_id'], s['core']) for s in sites], [('p0', '1,842')])  # 2016 (year) excluded
        self.assert_pre('n_text_mismatch', m)

    def test_pct_sum(self):
        def m(p):
            p['tables'][3]['rows'][0][1] = '71.2'  # rows no longer sum to 100 -> not a partition
            p['tables'][3]['rows'][1][1] = '19.8'
            p['tables'][3]['rows'][2][1] = '48.4'
        sites = plant.sites_pct_sum_break(SYN)
        self.assertEqual(len(sites), 3)
        self.assert_pre('pct_sum_break', m)

    def test_mean_range(self):
        self.assert_pre('mean_outside_range', lambda p: p['tables'][1]['header_rows'][0].__setitem__(3, 'Median'))

    def test_ci(self):
        def m(p):
            p['tables'][4]['header_rows'][0][2] = 'Range'
            p['paragraphs'][2]['text'] = 'The odds of voting were higher among union members (OR = 1.45).'
        self.assertEqual(len(plant.sites_ci_excludes_estimate(SYN)), 3)
        self.assert_pre('ci_excludes_estimate', m)

    def test_statcheck(self):
        def m(p):
            p['paragraphs'][3]['text'] = 'A difference-in-means test confirms this (p = .014).'
        self.assert_pre('statcheck_p_mismatch', m)

    def test_statcheck_rejects_inconsistent_original(self):
        p = copy.deepcopy(SYN)
        p['paragraphs'][3]['text'] = 'Test: t(120) = 2.50, p = .14.'
        self.assertEqual(plant.sites_statcheck_p_mismatch(p), [])

    def test_count_pct(self):
        def m(p):
            p['tables'][2]['rows'][1][1] = '96 (38.0)'
            p['tables'][2]['rows'][1][2] = '102 (61.0)'
        self.assertEqual(len(plant.sites_count_pct_mismatch(SYN)), 6)
        self.assert_pre('count_pct_mismatch', m)

    def test_text_table(self):
        def m(p):
            for para in p['paragraphs']:
                para['text'] = 'Nothing numeric here.'
        cores = sorted(s['core'] for s in plant.sites_text_table_transposition(SYN))
        self.assertEqual(cores, ['0.041', '0.112', '0.234'])
        self.assert_pre('text_table_transposition', m)


class TestMutations(unittest.TestCase):
    """Each planted copy has exactly one change, at the manifest location, and it breaks the invariant."""

    def check_common(self, etype):
        copies = plant_all(etype)
        self.assertTrue(copies, etype)
        self.assertLessEqual(len(copies), plant.MAX_COPIES)
        locs = set()
        for cp, loc, before, after in copies:
            self.assertNotEqual(before, after)
            d = diff_locations(SYN, cp)
            self.assertEqual(len(d), 1, (etype, d))
            if loc['kind'] == 'cell':
                self.assertEqual(d[0], ('cell', loc['table_id'], loc['row'], loc['col']))
                self.assertEqual(cell(SYN, loc), before)
                self.assertEqual(cell(cp, loc), after)
            else:
                self.assertEqual(d[0], ('text', loc['para_id']))
                new = plant.get_para(cp, loc['para_id'])['text']
                old = plant.get_para(SYN, loc['para_id'])['text']
                self.assertEqual(new[loc['start']:loc['end']], after)
                self.assertEqual(new[:loc['start']], old[:loc['start']])
                self.assertEqual(old.replace(before, after, 1) if old.count(before) == 1 else new, new)
            locs.add(json.dumps(loc, sort_keys=True))
        self.assertEqual(len(locs), len(copies), 'copies must be at distinct locations')
        return copies

    def test_star_flip(self):
        for cp, loc, before, after in self.check_common('star_flip'):
            t = table(cp, loc['table_id'])
            r = loc['row']
            coef = t['rows'][r][loc['col']]
            se = t['rows'][r + 1][loc['col']]
            b = f(re.match(r'[−\-\d.]+', coef).group())
            s = f(se.strip('()'))
            printed = len(re.search(r'\**$', coef).group())
            self.assertNotEqual(printed, stars_expected(b, s))
            self.assertEqual(abs(printed - stars_expected(b, s)), 1)

    def test_se_misreport(self):
        for cp, loc, before, after in self.check_common('se_misreport'):
            t = table(cp, loc['table_id'])
            r = loc['row']
            self.assertTrue(after.startswith('('))
            coef = t['rows'][r - 1][loc['col']]
            b = f(re.match(r'[−\-\d.]+', coef).group())
            printed = len(re.search(r'\**$', coef).group())
            self.assertNotEqual(printed, stars_expected(b, f(after.strip('()'))))
            self.assertEqual(printed, stars_expected(b, f(before.strip('()'))))

    def test_sign_flip(self):
        for cp, loc, before, after in self.check_common('sign_flip'):
            num_b = re.match(r'[−\-]?[\d.]+', before).group()
            num_a = re.match(r'[−\-]?[\d.]+', after).group()
            self.assertEqual(f(num_a), -f(num_b))
            self.assertIn(num_b, ' '.join(p['text'] for p in cp['paragraphs']))  # text still states old sign

    def test_n_text(self):
        for cp, loc, before, after in self.check_common('n_text_mismatch'):
            self.assertEqual(before, '1,842')
            self.assertNotIn(after.replace(',', ''), plant.all_table_cores(cp))
            self.assertEqual(len(after.replace(',', '')), 4)

    def test_pct_sum(self):
        for cp, loc, before, after in self.check_common('pct_sum_break'):
            row = table(cp, loc['table_id'])['rows'][loc['row']]
            s = sum(f(x) for x in row[1:4])
            self.assertGreaterEqual(abs(s - 100), 3.0)
            self.assertLessEqual(abs(s - 100), 40)

    def test_mean_range(self):
        for cp, loc, before, after in self.check_common('mean_outside_range'):
            row = table(cp, loc['table_id'])['rows'][loc['row']]
            self.assertEqual(loc['col'], 1)
            self.assertTrue(f(row[1]) < f(row[3]) or f(row[1]) > f(row[4]))

    def test_ci(self):
        for cp, loc, before, after in self.check_common('ci_excludes_estimate'):
            if loc['kind'] == 'text':
                m = re.search(r'([\d.]+), 95% CI \[([\d.]+), ([\d.]+)\]', after)
                self.assertTrue(after.startswith(m.group(1)) and after.endswith(']'))
                e, lo, hi = map(f, m.groups())
            else:
                row = table(cp, loc['table_id'])['rows'][loc['row']]
                e = f(row[1])
                lo, hi = map(f, row[2].split('–'))
            self.assertLess(lo, hi)
            self.assertTrue(e < lo or e > hi)

    def test_statcheck(self):
        for cp, loc, before, after in self.check_common('statcheck_p_mismatch'):
            m = re.fullmatch(r't\(120\) = ([\d.]+), p ([<=>]) (\d?\.\d+)', after)
            self.assertTrue(m, after)
            self.assertEqual(sum(x != y for x, y in zip(m.groups(), ('2.50', '=', '.014'))), 1)
            true_p = plant.t_p2(f(m.group(1)), 120)
            m = re.search(r'p ([<=>]) (\d?\.\d+)', after)
            op, p = m.group(1), f(m.group(2))
            if op == '=':
                self.assertTrue(p > 2 * true_p or p < true_p / 2)
            elif op == '<':
                self.assertLessEqual(p, true_p)
            else:
                self.assertGreaterEqual(p, true_p)

    def test_count_pct(self):
        for cp, loc, before, after in self.check_common('count_pct_mismatch'):
            m = re.match(r'(\d+) \(([\d.]+)\)', after)
            n_implied = int(m.group(1)) / f(m.group(2)) * 100
            self.assertGreater(abs(n_implied - 200) / 200, 0.03)

    def test_text_table(self):
        for cp, loc, before, after in self.check_common('text_table_transposition'):
            cores = plant.all_table_cores(cp)
            self.assertIn(before, cores)
            self.assertNotIn(after, cores)
            self.assertEqual(sorted(before), sorted(after))  # transposition keeps the digits


class TestRun(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.corpus = os.path.join(self.tmp, 'corpus')
        os.makedirs(self.corpus)
        p2 = copy.deepcopy(SYN)
        p2['paper_id'] = 'synthetic_0002'
        p2['paragraphs'] = p2['paragraphs'][:1]
        for pid, p in (('synthetic_0001', SYN), ('synthetic_0002', p2)):
            with open(os.path.join(self.corpus, pid + '.json'), 'w') as fh:
                json.dump(p, fh, ensure_ascii=False)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def digest(self, d):
        h = hashlib.sha256()
        for root, _, files in sorted(os.walk(d)):
            for fn in sorted(files):
                h.update(fn.encode())
                with open(os.path.join(root, fn), 'rb') as fh:
                    h.update(fh.read())
        return h.hexdigest()

    def test_deterministic_and_schema(self):
        before = self.digest(self.corpus)
        o1, o2 = os.path.join(self.tmp, 'o1'), os.path.join(self.tmp, 'o2')
        s1 = plant.run(self.corpus, o1, verbose=False)
        s2 = plant.run(self.corpus, o2, verbose=False)
        self.assertEqual(self.digest(o1), self.digest(o2))
        self.assertEqual(s1, s2)
        self.assertEqual(self.digest(self.corpus), before, 'originals must be untouched')
        split = plant.load_split()
        self.assertFalse(set(split['visible']) & set(split['heldout']))
        for sp in ('visible', 'heldout'):
            with open(os.path.join(o1, sp, 'manifest.jsonl')) as fh:
                recs = [json.loads(line) for line in fh]
            self.assertTrue(recs)
            ids = [r['planted_id'] for r in recs]
            self.assertEqual(len(ids), len(set(ids)))
            for r in recs:
                self.assertEqual(set(r), {'planted_id', 'orig_id', 'error_type', 'location', 'before', 'after'})
                self.assertIn(r['error_type'], split[sp])
                self.assertTrue(r['planted_id'].startswith(r['orig_id'] + '__'))
                with open(os.path.join(o1, sp, r['planted_id'] + '.json')) as fh:
                    pp = json.load(fh)
                self.assertEqual(pp['paper_id'], r['planted_id'])
                loc = r['location']
                if loc['kind'] == 'cell':
                    self.assertEqual(set(loc), {'kind', 'table_id', 'row', 'col'})
                    self.assertEqual(cell(pp, loc), r['after'])
                else:
                    self.assertEqual(set(loc), {'kind', 'para_id', 'start', 'end'})
                    self.assertEqual(plant.get_para(pp, loc['para_id'])['text'][loc['start']:loc['end']], r['after'])
        self.assertEqual(s1['papers_with_any_applicable'], 2)
        self.assertEqual(s1['per_split']['visible'] + s1['per_split']['heldout'],
                         sum(v['copies'] for v in s1['per_type'].values()))
        for t in ALL_TYPES:
            self.assertGreater(s1['per_type'][t]['copies'], 0, t)


if __name__ == '__main__':
    unittest.main(verbosity=1)
