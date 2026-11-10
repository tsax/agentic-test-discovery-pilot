"""Hill-climbing test discovery: an LLM agent proposes a consistency test, the
black-box scorer (climb/score.py) scores it on VISIBLE planted errors over DEV
papers, the agent reads the feedback and revises. Every step is logged.

The agent never sees: the planter code, held-out error types, or test papers.

Usage:
  uv run --with anthropic python -I climb/climb.py --name A --family A --iters 15 --budget 8
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import anthropic

ROOT = Path(__file__).resolve().parent.parent
MODEL = "claude-opus-5-5"
# USD per million tokens, Claude Opus 5.5 first-party API (cache write = 1.25x input, 5-min TTL)
PRICE = {"input": 4.00, "output": 20.00, "cache_read": 0.20, "cache_write": 5.00}

FAMILIES = {
    "A": "Start from INFERENTIAL consistency: whether reported estimates, standard errors, "
         "test statistics, confidence intervals, p-values and significance markers agree with "
         "each other and with the table's own stated conventions.",
    "B": "Start from ARITHMETIC and ACCOUNTING consistency: whether counts, sample sizes, "
         "percentages, totals, means and ranges reported in tables and text agree with each "
         "other within the paper.",
    "C": "Start wherever you think published quantitative social-science papers are most likely "
         "to contain checkable internal inconsistencies that existing tools (statcheck, GRIM) miss. "
         "Choose your own principle.",
}

SYSTEM = """You are a research agent discovering NEW automated checks for errors in published quantitative social-science papers.

A "test" is a small deterministic Python module that reads one parsed paper and flags specific locations where the paper's own reported numbers are mutually inconsistent. Existing hand-written checks (statcheck: recomputes p-values from APA-style test statistics in text; GRIM: tests whether means are possible given N for integer data) cover only a narrow slice. Your job is to find a different, general consistency principle and implement it well.

## Paper format (input to check)
{schema_paper}

## Interface you must implement
```python
def applicable(paper: dict) -> bool:
    \"\"\"True if the paper contains at least one input this test can evaluate.\"\"\"

def check(paper: dict) -> list[dict]:
    \"\"\"Return flags: {{"location": {{"kind": "cell", "table_id": str, "row": int, "col": int}} | {{"kind": "text", "para_id": str, "start": int, "end": int}}, "message": str}}.
    Deterministic, offline, < 2 s per paper, Python stdlib (+ numpy/scipy allowed). Abstain (return []) when inputs are missing or ambiguous.\"\"\"
```
`row` indexes table["rows"] (header rows excluded); `start`/`end` are character offsets into that paragraph's text.

## How you are scored (a black box you can only query by submitting code)
- Errors have been planted, one per copy, into copies of real papers. Recall = fraction of planted errors where your test raises a NEW flag at the planted location (same table and row ±1, or overlapping text span in the same paragraph) that it did not raise on the unmodified original.
- False alarms: unmodified original papers are mostly, but not certainly, correct. Each original paper you flag costs you. score = recall − 0.5 × (fraction of original papers flagged).
- You see feedback only for a VISIBLE set of error types on DEVELOPMENT papers. Your final test will be scored once on DIFFERENT, unseen error types and on unseen papers. Detectors overfitted to the exact visible mutations will generalize poorly; a sound, general consistency principle will generalize well.
- Never hard-code paper ids, table ids, or specific values. No network access, no file access.

## Direction for this run
{family}

## Each turn
Reply with exactly:
HYPOTHESIS: <one sentence: the consistency principle this version checks>
CHANGES: <one or two sentences: what you changed from the last version and why>
```python
<the complete module>
```
Keep ONE coherent test (one principle, possibly with several sub-cases), not a grab bag of unrelated checks."""


def extract(text: str) -> tuple[str | None, str, str]:
    m = re.findall(r"```python\n(.*?)```", text, flags=re.S)
    code = m[-1] if m else None
    hyp = re.search(r"HYPOTHESIS:\s*(.+)", text)
    chg = re.search(r"CHANGES:\s*(.+)", text)
    return code, hyp.group(1).strip() if hyp else "", chg.group(1).strip() if chg else ""


def cost(u) -> float:
    cr = getattr(u, "cache_read_input_tokens", 0) or 0
    cw = getattr(u, "cache_creation_input_tokens", 0) or 0
    return (u.input_tokens * PRICE["input"] + u.output_tokens * PRICE["output"]
            + cr * PRICE["cache_read"] + cw * PRICE["cache_write"]) / 1e6


def sample_papers(ids: list[str], k: int = 2) -> str:
    out = []
    for pid in ids[:k]:
        p = json.loads((ROOT / "corpus" / "parsed" / f"{pid}.json").read_text())
        out.append({
            "paper_id": p["paper_id"], "title": p.get("title"),
            "paragraphs": [{**q, "text": q["text"][:600]} for q in p["paragraphs"][:6]],
            "tables": [{**t, "rows": t["rows"][:14]} for t in p["tables"][:2]],
        })
    return json.dumps(out, indent=1)[:30000]


def run_score(path: Path, out: Path) -> dict:
    proc = subprocess.run(
        [sys.executable, "-I", str(ROOT / "climb" / "score.py"), str(path),
         "--split", "visible", "--papers", "dev", "--feedback", "6", "--out", str(out)],
        capture_output=True, text=True, timeout=900,
    )
    if proc.returncode != 0:
        return {"score": -1.0, "load_error": proc.stderr[-2000:]}
    res = json.loads(out.read_text())
    res.pop("_orig_flags", None)
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--family", required=True, choices=sorted(FAMILIES))
    ap.add_argument("--iters", type=int, default=15)
    ap.add_argument("--budget", type=float, default=8.0, help="USD cap for this climber")
    a = ap.parse_args()

    run_dir = ROOT / "climb" / "runs" / a.name
    run_dir.mkdir(parents=True, exist_ok=True)
    log = (run_dir / "trajectory.jsonl").open("a")

    schema = (ROOT / "docs" / "SCHEMA.md").read_text()
    schema_paper = schema.split("## 2.")[0].split("## 1.")[1]
    system = SYSTEM.format(schema_paper=schema_paper, family=FAMILIES[a.family])
    dev_ids = json.loads((ROOT / "climb" / "paper_split.json").read_text())["dev"]

    client = anthropic.Anthropic()
    messages: list[dict] = [{"role": "user", "content":
        "Here are two development papers as parsed (truncated) so you can see real data:\n"
        + sample_papers(dev_ids) + "\n\nWrite version 1 of your test."}]

    spent, best = 0.0, {"score": float("-inf"), "iter": None}
    for it in range(1, a.iters + 1):
        if spent >= a.budget:
            print(f"[{a.name}] budget reached (${spent:.2f})")
            break
        t0 = time.time()
        with client.beta.messages.stream(
            model=MODEL, max_tokens=32000, system=system, messages=messages,
            thinking={"type": "adaptive"}, output_config={"effort": "high"},
            cache_control={"type": "ephemeral"},
            betas=["server-side-fallback-2026-07-01"], fallbacks="default",
        ) as stream:
            resp = stream.get_final_message()
        c = cost(resp.usage)
        spent += c
        messages.append({"role": "assistant", "content": resp.content})
        if resp.stop_reason == "refusal":
            log.write(json.dumps({"iter": it, "event": "refusal", "cost_usd": round(c, 4)}) + "\n"); log.flush()
            messages.append({"role": "user", "content": "That reply was declined. Continue with the next version of the test."})
            continue

        text = "".join(b.text for b in resp.content if b.type == "text")
        code, hyp, chg = extract(text)
        if not code:
            messages.append({"role": "user", "content": "No ```python block found. Reply in the required format."})
            log.write(json.dumps({"iter": it, "event": "no_code", "cost_usd": round(c, 4)}) + "\n"); log.flush()
            continue

        path = run_dir / f"iter_{it:02d}.py"
        path.write_text(code)
        res = run_score(path, run_dir / f"iter_{it:02d}.metrics.json")
        kept = res.get("score", -1) > best["score"]
        if kept:
            best = {"score": res["score"], "iter": it}
            (run_dir / "best.py").write_text(code)
        entry = {
            "iter": it, "hypothesis": hyp, "changes": chg,
            "code_sha256": hashlib.sha256(code.encode()).hexdigest()[:16],
            "score": res.get("score"), "recall": res.get("recall"), "fp_papers": res.get("fp_papers"),
            "applicable_papers": res.get("applicable_papers"), "recall_by_type": res.get("recall_by_type"),
            "n_errors": res.get("n_errors"), "decision": "keep" if kept else "reject",
            "best_so_far": best, "tokens": {"in": resp.usage.input_tokens, "out": resp.usage.output_tokens,
            "cache_read": getattr(resp.usage, "cache_read_input_tokens", 0),
            "cache_write": getattr(resp.usage, "cache_creation_input_tokens", 0)},
            "cost_usd": round(c, 4), "spent_usd": round(spent, 4), "seconds": round(time.time() - t0, 1),
        }
        log.write(json.dumps(entry) + "\n"); log.flush()
        print(f"[{a.name}] iter {it}: score={res.get('score')} recall={res.get('recall')} "
              f"fp={res.get('fp_papers')} {entry['decision']} ${spent:.2f}", flush=True)

        feedback = {k: v for k, v in res.items() if k not in ("test",)}
        messages.append({"role": "user", "content":
            f"Scorer feedback for version {it} (best so far: version {best['iter']}, score {best['score']}):\n"
            + json.dumps(feedback, indent=1, default=str)[:40000]
            + "\n\nRevise the test. If a direction has stopped improving, you may change approach within the same principle."})

    (run_dir / "summary.json").write_text(json.dumps({"name": a.name, "family": a.family, "best": best,
                                                       "spent_usd": round(spent, 4)}, indent=1))
    print(f"[{a.name}] done. best={best} spent=${spent:.2f}")


if __name__ == "__main__":
    main()
