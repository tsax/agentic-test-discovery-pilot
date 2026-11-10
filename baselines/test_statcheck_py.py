"""Tests for the Python statcheck port. Run:
uv run --quiet --with scipy --with pytest python -I -m pytest baselines/test_statcheck_py.py
(or: python -m unittest baselines.test_statcheck_py)"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import statcheck_py as S  # noqa: E402


def one(text, **kw):
    res = S.extract(text, **kw)
    assert len(res) == 1, res
    return res[0]


class TestCoreCases(unittest.TestCase):
    def test_t_error(self):
        r = one("t(28) = 2.20, p = .01")
        self.assertTrue(r["error"])
        self.assertFalse(r["decision_error"])
        self.assertAlmostEqual(r["computed_p"], 0.036, places=3)

    def test_t_consistent(self):
        self.assertFalse(one("t(28) = 2.20, p = .036")["error"])

    def test_F_less_than(self):
        r = one("F(1, 98) = 4.12, p < .05")
        self.assertFalse(r["error"])
        self.assertLess(r["computed_p"], .05)

    def test_r_consistent(self):
        r = one("r(48) = .30, p = .03")
        self.assertFalse(r["error"])
        self.assertAlmostEqual(r["computed_p"], 0.034, places=3)

    def test_chi2_error(self):
        r = one("χ2(2, N = 120) = 3.10, p = .02")
        self.assertTrue(r["error"])
        self.assertTrue(r["decision_error"])  # .02 reported, ~.21 computed

    def test_decision_error_t(self):
        r = one("t(28) = 1.50, p = .04")
        self.assertTrue(r["error"])
        self.assertTrue(r["decision_error"])

    def test_decision_error_ns(self):
        r = one("t(40) = 3.10, ns")
        self.assertTrue(r["error"])
        self.assertTrue(r["decision_error"])

    def test_rounding_edge_accepted(self):
        # p(t=2.0, 28) = .0553 -> .06, but t in [1.95, 2.05] gives p down to .0498 -> .05.
        r = one("t(28) = 2.0, p = .05")
        self.assertGreater(round(r["computed_p"], 2), .05)
        self.assertFalse(r["error"])
        # a statistic whose whole rounding interval gives p > .055 is still an error
        self.assertTrue(one("t(28) = 1.9, p = .05")["error"])


class TestExtraction(unittest.TestCase):
    def test_variants(self):
        txt = ("t (28) = −2.20, p = 0.036; F(2, 120) = 3.5, p < .05; "
               "chi-square(1) = 3.84, p = .05, X2(1) = 10.83, p < .001, "
               "Chi2(3, N = 1,200) = 7.81, p = .05; z = 1.96, p = .05; Q(10) = 18.31, p = .05; "
               "t(1,234) = 1,234.5, p < .001; r(100) = –.20, p = .046")
        res = S.extract(txt)
        self.assertEqual([r["test"] for r in res], ["t", "F", "chi2", "chi2", "chi2", "z", "Q", "t", "r"])
        self.assertFalse(any(r["error"] for r in res), [r["raw"] for r in res if r["error"]])
        self.assertEqual(res[0]["stat"], -2.2)
        self.assertEqual(res[7]["df1"], 1234)
        self.assertEqual(res[7]["stat"], 1234.5)

    def test_span_maps_to_original(self):
        txt = "As predicted, t(28) = 2.20, p = .01 in the sample."
        r = one(txt)
        self.assertEqual(txt[r["start"]:r["end"]], "t(28) = 2.20, p = .01")

    def test_not_extracted(self):
        self.assertEqual(S.extract("β = .21, p < .01; r = .40, p < .01; (p = 0.004)"), [])

    def test_p_zero_error(self):
        self.assertTrue(one("F(1, 50) = 40.2, p = .000")["error"])

    def test_one_tailed(self):
        # p(t=1.80, 30) two-tailed = .082; one-tailed .041
        self.assertTrue(one("t(30) = 1.80, p = .041")["error"])
        self.assertFalse(one("t(30) = 1.80, p = .041", one_tailed_txt=True)["error"])
        self.assertTrue(one("t(30) = 1.80, p = .041, one-tailed")["one_tailed_accepted"])

    def test_greater_than(self):
        self.assertFalse(one("F(1, 30) = 0.50, p > .05")["error"])
        self.assertTrue(one("F(1, 30) = 9.50, p > .05")["error"])


class TestPaperInterface(unittest.TestCase):
    def test_check_and_applicable(self):
        paper = {"paper_id": "x", "paragraphs": [
            {"para_id": "p0", "text": "No stats here."},
            {"para_id": "p1", "text": "We found t(28) = 2.20, p = .01 and r(48) = .30, p = .03."}]}
        self.assertTrue(S.applicable(paper))
        fl = S.check(paper)
        self.assertEqual(len(fl), 1)
        loc = fl[0]["location"]
        self.assertEqual(loc["para_id"], "p1")
        self.assertEqual(paper["paragraphs"][1]["text"][loc["start"]:loc["end"]], "t(28) = 2.20, p = .01")
        self.assertIn("0.03", fl[0]["message"])
        self.assertFalse(S.applicable({"paper_id": "y", "paragraphs": []}))


if __name__ == "__main__":
    unittest.main()
