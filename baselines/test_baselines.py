"""Unit tests for the baselines. Run: python3 -m unittest baselines.test_baselines  (or pytest)."""
import os
import shutil
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import grim_check as G  # noqa: E402

HAVE_R = bool(shutil.which("Rscript") or os.path.exists("/opt/homebrew/bin/Rscript"))


def _paper(tables=(), paragraphs=()):
    return {"paper_id": "test_paper", "doi": "x", "title": "t", "year": 2020,
            "data_availability": "", "paragraphs": list(paragraphs), "tables": list(tables)}


class TestGrimCore(unittest.TestCase):
    def test_known_inconsistent(self):
        self.assertIs(G.grim_consistent("5.19", 28), False)

    def test_known_consistent(self):
        self.assertIs(G.grim_consistent("5.18", 28), True)   # 145/28 = 5.1786
        self.assertIs(G.grim_consistent("5.21", 28), True)   # 146/28 = 5.2143

    def test_rounding_conventions(self):
        # 1/8 = 0.125: half-up -> 0.13, half-even -> 0.12; both must be accepted.
        self.assertIs(G.grim_consistent("0.13", 8), True)
        self.assertIs(G.grim_consistent("0.12", 8), True)

    def test_no_power(self):
        self.assertIsNone(G.grim_consistent("5.19", 100))
        self.assertIsNone(G.grim_consistent("5.1", 28))

    def test_effective_n_items(self):
        # n=20, 3 items -> n_eff=60; 3.48*60 = 208.8 -> 209/60 = 3.4833 -> 3.48 consistent
        self.assertIs(G.grim_consistent("3.48", 60), True)
        self.assertIs(G.grim_consistent("3.48", 20), False)


class TestGrimTables(unittest.TestCase):
    def _table(self, caption, rows, header=(("Measure", "M", "SD"),), footnotes=""):
        return {"table_id": "t1", "label": "Table 1", "caption": caption, "footnotes": footnotes,
                "header_rows": [list(h) for h in header], "rows": [list(r) for r in rows]}

    def test_flags_inconsistent_cell(self):
        t = self._table("Ratings on a 7-point Likert item (N = 28)",
                        [["Liking", "5.19", "1.20"], ["Trust", "5.18", "1.10"]])
        flags = G.check(_paper([t]))
        self.assertEqual(len(flags), 1)
        self.assertEqual(flags[0]["location"], {"kind": "cell", "table_id": "t1", "row": 0, "col": 1})
        self.assertTrue(G.applicable(_paper([t])))

    def test_mean_sd_combined_and_header_n(self):
        t = self._table("Single-item ratings (1-7 Likert).",
                        [["Liking", "5.19 (1.20)", "5.18 ± 1.1"]],
                        header=(("Item", "Control (n = 28) M (SD)", "Treatment (n = 28) M (SD)"),))
        flags = G.check(_paper([t]))
        self.assertEqual([f["location"]["col"] for f in flags], [1])

    def test_n_column(self):
        t = self._table("Number of errors per condition",
                        [["A", "28", "5.19"], ["B", "28", "5.18"]],
                        header=(("Condition", "n", "Mean"),))
        flags = G.check(_paper([t]))
        self.assertEqual([(f["location"]["row"], f["location"]["col"]) for f in flags], [(0, 2)])

    def test_skips_continuous(self):
        t = self._table("Reaction time (ms) by condition, N = 28", [["RT", "5.19", "1.2"]])
        self.assertEqual(G.check(_paper([t])), [])
        self.assertFalse(G.applicable(_paper([t])))
        self.assertIn("continuous", G.analyze(_paper([t]))[0]["reason"])

    def test_skips_without_integer_evidence(self):
        t = self._table("Descriptives, N = 28", [["Thing", "5.19", "1.2"]])
        self.assertEqual(G.check(_paper([t])), [])
        self.assertIn("integer", G.analyze(_paper([t]))[0]["reason"])

    def test_skips_multi_item_unknown_k(self):
        t = self._table("Means of the subscales (Likert 1-5), N = 28", [["Anxiety", "3.19", "1.2"]])
        self.assertEqual(G.check(_paper([t])), [])

    def test_multi_item_known_k(self):
        # 3.48 with n=20 alone is inconsistent, but 3-item scale -> n_eff 60 -> consistent
        t = self._table("Mean of 3 items on a 5-point Likert scale, N = 20", [["Anxiety", "3.48", "1.2"]])
        self.assertEqual(G.check(_paper([t])), [])
        self.assertTrue(G.applicable(_paper([t])))

    def test_no_power_large_n(self):
        t = self._table("Likert item ratings, N = 280", [["Liking", "5.19", "1.20"]])
        self.assertEqual(G.check(_paper([t])), [])
        self.assertFalse(G.applicable(_paper([t])))

    def test_row_wise_mean(self):
        t = {"table_id": "t2", "label": "Table 2", "caption": "Number of correct answers", "footnotes": "",
             "header_rows": [["", "Group A", "Group B"]],
             "rows": [["N", "28", "28"], ["Mean", "5.19", "5.18"], ["SD", "1.2", "1.3"]]}
        flags = G.check(_paper([t]))
        self.assertEqual([f["location"] for f in flags], [{"kind": "cell", "table_id": "t2", "row": 1, "col": 1}])

    def test_no_tables(self):
        self.assertEqual(G.check(_paper()), [])
        self.assertFalse(G.applicable(_paper()))


@unittest.skipUnless(HAVE_R, "Rscript not installed")
class TestStatcheck(unittest.TestCase):
    def setUp(self):
        import statcheck_check as S
        self.S = S

    def test_inconsistent_vs_consistent(self):
        bad = "Participants differed, t(28) = 2.20, p = .01, as expected."
        good = "Participants differed, t(28) = 2.20, p = .036, as expected."
        none = "No statistics in this paragraph."
        p1 = _paper(paragraphs=[{"para_id": "p0", "section": "Results", "text": bad}])
        p1["paper_id"] = "bad"
        p2 = _paper(paragraphs=[{"para_id": "p0", "section": "Results", "text": good}])
        p2["paper_id"] = "good"
        p3 = _paper(paragraphs=[{"para_id": "p0", "section": "Results", "text": none}])
        p3["paper_id"] = "none"
        res = self.S.run_batch_full([p1, p2, p3])
        self.assertEqual(len(res["bad"]["flags"]), 1)
        loc = res["bad"]["flags"][0]["location"]
        self.assertEqual(loc["kind"], "text")
        self.assertEqual(loc["para_id"], "p0")
        self.assertIn("2.20", bad[loc["start"]:loc["end"]])
        self.assertLess(loc["end"] - loc["start"], len(bad))
        self.assertEqual(res["good"]["flags"], [])
        self.assertTrue(res["good"]["applicable"])
        self.assertFalse(res["none"]["applicable"])
        self.assertEqual(self.S.check(p1), res["bad"]["flags"])

    def test_locate_tolerant(self):
        self.assertEqual(self.S.locate("t(28) = 2.20, p = .01", "x t (28)=2.20, p=.01 y"), (2, 20))


if __name__ == "__main__":
    unittest.main()
