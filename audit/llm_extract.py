"""Audit: are statcheck/GRIM applicability numbers an artifact of our regex extraction?

An LLM independently reads each paper (text AND tables) and lists:
  (a) every reported null-hypothesis test with statistic, df and p (statcheck's inputs)
  (b) every reported mean with the n it is based on and whether the measure is
      integer-valued (GRIM's inputs)
The SAME deterministic code then checks them: statcheck_py.extract() on a normalized
APA string, and grim_check.grim_has_power()/grim_consistent().

Usage: uv run --with anthropic --with scipy python -I audit/llm_extract.py
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import anthropic

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "baselines"))
sys.path.insert(0, str(ROOT / "verify"))
import grim_check as G  # noqa: E402
import statcheck_py as SC  # noqa: E402
from verify import render_paper  # noqa: E402

MODEL = "claude-opus-5-5"
PRICE = {"input": 4.00, "output": 20.00}
OUT = ROOT / "audit"

PROMPT = """Read this published paper (body text and every table, with captions and footnotes) and extract two lists exactly as printed. Be exhaustive; include items in tables as well as text. Do not compute anything.

(a) tests: every reported null-hypothesis significance test where a test statistic of type t, F, chi-square, r (correlation with df), z, or Q is printed together with a p-value (exact or inequality) for that same test. Include its degrees of freedom if printed (for chi-square also N if printed). Skip regression coefficient tables that give only estimates/SEs/stars without a test statistic AND p for the same row. Skip goodness-of-fit indices like CFI, RMSEA, chi2/df ratios.

(b) means: every reported mean (average) of a variable that is printed together with, or clearly attributable to, a specific sample size n for that mean (n of that group/cell, not the whole study if the cell n differs). Record whether the underlying per-person measure is integer-valued (a single Likert item, a count, a yes/no coded 0/1, age in whole years) and how many items are averaged if it is a multi-item scale.

PAPER:
{paper}"""

SCHEMA = {
    "type": "object", "additionalProperties": False, "required": ["tests", "means"],
    "properties": {
        "tests": {"type": "array", "items": {"type": "object", "additionalProperties": False,
            "required": ["where", "stat_type", "df1", "df2", "stat_value", "p_comparator", "p_value", "as_printed"],
            "properties": {
                "where": {"type": "string", "description": "para_id like p12 or table_id like t3"},
                "stat_type": {"type": "string", "enum": ["t", "F", "chi2", "r", "z", "Q"]},
                "df1": {"type": "string", "description": "as printed, or empty"},
                "df2": {"type": "string", "description": "as printed (F only), or empty"},
                "stat_value": {"type": "string", "description": "as printed, e.g. 2.20 or -1.96"},
                "p_comparator": {"type": "string", "enum": ["=", "<", ">", "ns"]},
                "p_value": {"type": "string", "description": "as printed, e.g. .036 or 0.001"},
                "as_printed": {"type": "string"}}}},
        "means": {"type": "array", "items": {"type": "object", "additionalProperties": False,
            "required": ["where", "variable", "mean_as_printed", "n", "n_evidence", "integer_valued", "n_items"],
            "properties": {
                "where": {"type": "string"},
                "variable": {"type": "string"},
                "mean_as_printed": {"type": "string"},
                "n": {"type": "integer"},
                "n_evidence": {"type": "string", "description": "short quote showing where n comes from"},
                "integer_valued": {"type": "string", "enum": ["yes", "no", "unclear"]},
                "n_items": {"type": "integer", "description": "items averaged; 1 if single item; 0 if unknown"}}}},
    },
}


def apa(t: dict) -> str:
    st = {"t": "t", "F": "F", "chi2": "χ2", "r": "r", "z": "z", "Q": "Q"}[t["stat_type"]]
    df = ", ".join(x for x in (t["df1"], t["df2"]) if x.strip())
    head = f"{st}({df})" if df and t["stat_type"] != "z" else st
    p = "ns" if t["p_comparator"] == "ns" else f"p {t['p_comparator']} {t['p_value']}"
    return f"{head} = {t['stat_value']}, {p}"


def one(client, path: Path) -> dict:
    paper = json.loads(path.read_text())
    with client.messages.stream(
        model=MODEL, max_tokens=32000, messages=[{"role": "user", "content": PROMPT.format(paper=render_paper(paper, 200_000))}],
        thinking={"type": "adaptive"}, output_config={"effort": "medium", "format": {"type": "json_schema", "schema": SCHEMA}},
    ) as stream:
        resp = stream.get_final_message()
    cost = (resp.usage.input_tokens * PRICE["input"] + resp.usage.output_tokens * PRICE["output"]) / 1e6
    if resp.stop_reason != "end_turn":
        return {"paper_id": paper["paper_id"], "error": resp.stop_reason, "cost_usd": cost}
    data = json.loads(next(b.text for b in resp.content if b.type == "text"))

    # statcheck: same deterministic checker, fed a normalized APA string
    sc_checked, sc_errors, sc_unparsed = [], 0, 0
    for t in data["tests"]:
        s = apa(t)
        res = SC.extract(s)
        if not res:
            sc_unparsed += 1
            continue
        r = res[0]
        sc_checked.append({**t, "apa": s, "error": bool(r.get("error")), "computed_p": r.get("computed_p")})
        sc_errors += bool(r.get("error"))

    # GRIM: same deterministic functions
    grim_rows = []
    for m in data["means"]:
        n_eff = m["n"] * (m["n_items"] if m["n_items"] and m["n_items"] > 1 else 1)
        d = G.decimals_of(m["mean_as_printed"]) if any(ch.isdigit() for ch in m["mean_as_printed"]) else 0
        powered = m["integer_valued"] == "yes" and n_eff > 0 and d > 0 and G.grim_has_power(n_eff, d)
        cons = G.grim_consistent(m["mean_as_printed"], n_eff) if powered else None
        grim_rows.append({**m, "n_eff": n_eff, "decimals": d, "grim_testable": powered, "grim_consistent": cons})

    unparsed = [{**t, "apa": apa(t)} for t in data["tests"] if not SC.extract(apa(t))]
    return {"paper_id": paper["paper_id"], "cost_usd": round(cost, 4), "unparsed_tests": unparsed,
            "n_tests_found": len(data["tests"]), "n_tests_checkable": len(sc_checked), "n_tests_unparsed": sc_unparsed,
            "statcheck_errors": sc_errors, "tests": sc_checked,
            "n_means_found": len(data["means"]),
            "n_means_integer": sum(m["integer_valued"] == "yes" for m in data["means"]),
            "n_grim_testable": sum(r["grim_testable"] for r in grim_rows),
            "grim_inconsistent": sum(r["grim_consistent"] is False for r in grim_rows), "means": grim_rows}


def main() -> None:
    client = anthropic.Anthropic()
    paths = sorted((ROOT / "corpus" / "parsed").glob("*.json"))
    if len(sys.argv) > 1:  # optional: re-run only listed paper ids, writing to a separate file
        keep = set(sys.argv[1].split(","))
        paths = [p for p in paths if p.stem[-7:] in keep]
    with ThreadPoolExecutor(8) as ex:
        results = list(ex.map(lambda p: one(client, p), paths))
    suffix = "_rerun" if len(sys.argv) > 1 else ""
    (OUT / f"llm_extraction{suffix}.json").write_text(json.dumps(results, indent=1, default=str))
    ok = [r for r in results if "error" not in r]
    summary = {
        "papers": len(results), "failed": len(results) - len(ok),
        "papers_with_any_test_found": sum(r["n_tests_found"] > 0 for r in ok),
        "papers_statcheck_applicable": sum(r["n_tests_checkable"] > 0 for r in ok),
        "tests_found": sum(r["n_tests_found"] for r in ok), "tests_checkable": sum(r["n_tests_checkable"] for r in ok),
        "statcheck_errors": sum(r["statcheck_errors"] for r in ok),
        "papers_with_statcheck_error": sum(r["statcheck_errors"] > 0 for r in ok),
        "papers_with_any_mean_and_n": sum(r["n_means_found"] > 0 for r in ok),
        "means_found": sum(r["n_means_found"] for r in ok), "means_integer_valued": sum(r["n_means_integer"] for r in ok),
        "papers_grim_applicable": sum(r["n_grim_testable"] > 0 for r in ok),
        "means_grim_testable": sum(r["n_grim_testable"] for r in ok), "grim_inconsistent": sum(r["grim_inconsistent"] for r in ok),
        "cost_usd": round(sum(r["cost_usd"] for r in results), 2),
    }
    (OUT / f"summary{suffix}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
