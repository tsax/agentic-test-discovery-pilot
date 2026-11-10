"""Blind verification of flags raised on REAL (unmodified) papers.

Each flag is shown to a fresh model call with the paper's tables and text, the
location, and the claimed inconsistency. The verifier is not told which test
raised the flag, and is instructed to be skeptical. Output: one label per flag:
  confirmed  - the paper's own numbers are genuinely inconsistent as claimed
  false_alarm - the numbers are consistent once conventions/context are read correctly
  unclear    - cannot be decided from the paper alone

Usage:
  uv run --with anthropic python -I verify/verify.py flags.jsonl --out verify/verdicts.jsonl
Input lines: {"flag_id", "paper_id", "location", "message"} (other keys are ignored and never sent).
"""
from __future__ import annotations

import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import anthropic

ROOT = Path(__file__).resolve().parent.parent
MODEL = "claude-opus-5-5"
PRICE = {"input": 4.00, "output": 20.00, "cache_read": 0.20, "cache_write": 5.00}

PROMPT = """You are checking one claimed numerical inconsistency in a published paper. An automated tool raised it; such tools are often wrong because they misread table conventions (e.g. which rows are standard errors vs coefficients, what the star legend is, rounding, subsamples, reported vs derived quantities). Be skeptical and precise.

Decide, using ONLY the paper content below:
- "confirmed": the paper's own reported numbers are genuinely inconsistent in the way claimed (or in a closely related way at the same location), beyond rounding.
- "false_alarm": once the table/text is read correctly, the numbers are consistent, or the tool misparsed something.
- "unclear": it cannot be decided from the paper alone (e.g. needs the data, or the convention is ambiguous).

Claimed problem: {message}
Location: {location}
{located}

Paper (title, table captions/footnotes/rows, and body paragraphs; tables indexed with row numbers matching the location):
{paper}
"""

SCHEMA = {
    "type": "object",
    "properties": {
        "label": {"type": "string", "enum": ["confirmed", "false_alarm", "unclear"]},
        "explanation": {"type": "string"},
        "evidence": {"type": "string", "description": "the exact numbers/cells/text you relied on"},
    },
    "required": ["label", "explanation", "evidence"],
    "additionalProperties": False,
}


def render_paper(p: dict, limit: int = 120_000) -> str:
    parts = [f"TITLE: {p.get('title')}"]
    for t in p["tables"]:
        parts.append(f"\n[{t['table_id']}] {t.get('label')}: {t.get('caption','')}")
        for h in t.get("header_rows", []):
            parts.append("   header | " + " | ".join(h))
        for i, r in enumerate(t["rows"]):
            parts.append(f"   row {i} | " + " | ".join(r))
        if t.get("footnotes"):
            parts.append(f"   footnotes: {t['footnotes']}")
    for q in p["paragraphs"]:
        parts.append(f"\n[{q['para_id']}] ({q.get('section','')}) {q['text']}")
    return "\n".join(parts)[:limit]


def locate(p: dict, loc: dict) -> str:
    if loc.get("kind") == "cell":
        t = next((t for t in p["tables"] if t["table_id"] == loc.get("table_id")), None)
        if t and 0 <= loc.get("row", -1) < len(t["rows"]):
            return f"Flagged row: {t['rows'][loc['row']]}"
    if loc.get("kind") == "text":
        q = next((q for q in p["paragraphs"] if q["para_id"] == loc.get("para_id")), None)
        if q:
            s, e = loc.get("start", 0), loc.get("end", len(q["text"]))
            return f"Flagged text: ...{q['text'][max(0, s-200):e+200]}..."
    return ""


def verify_one(client: anthropic.Anthropic, item: dict) -> dict:
    p = json.loads((ROOT / "corpus" / "parsed" / f"{item['paper_id']}.json").read_text())
    content = PROMPT.format(message=item["message"], location=json.dumps(item["location"]),
                            located=locate(p, item["location"]), paper=render_paper(p))
    resp = client.messages.create(
        model=MODEL, max_tokens=16000, messages=[{"role": "user", "content": content}],
        thinking={"type": "adaptive"}, output_config={"effort": "high", "format": {"type": "json_schema", "schema": SCHEMA}},
    )
    u = resp.usage
    c = (u.input_tokens * PRICE["input"] + u.output_tokens * PRICE["output"]
         + (u.cache_read_input_tokens or 0) * PRICE["cache_read"]
         + (u.cache_creation_input_tokens or 0) * PRICE["cache_write"]) / 1e6
    if resp.stop_reason == "refusal":
        return {"flag_id": item["flag_id"], "label": "unclear", "explanation": "verifier declined", "evidence": "", "cost_usd": c}
    data = json.loads(next(b.text for b in resp.content if b.type == "text"))
    return {"flag_id": item["flag_id"], **data, "cost_usd": round(c, 4)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("flags", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    items = [json.loads(l) for l in a.flags.read_text().splitlines() if l.strip()]
    blind = [{k: it[k] for k in ("flag_id", "paper_id", "location", "message")} for it in items]
    client = anthropic.Anthropic()
    with ThreadPoolExecutor(a.workers) as ex:
        results = list(ex.map(lambda it: verify_one(client, it), blind))
    with a.out.open("w") as f:
        for r in results:
            f.write(json.dumps(r) + "\n")
    from collections import Counter
    print(Counter(r["label"] for r in results), f"cost=${sum(r['cost_usd'] for r in results):.2f}")


if __name__ == "__main__":
    main()
