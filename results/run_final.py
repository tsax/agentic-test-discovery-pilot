"""Final, one-shot scoring. Pre-registered rule (committed before climbers finished):

1. The agent's test = the single climber `best.py` with the highest VISIBLE-type score
   on DEV papers (as recorded in climb/runs/*/summary.json). Ties -> fewer fp_papers.
2. Every test (GRIM, statcheck port, all three climbers' best) is scored ONCE on TEST
   papers for both visible and held-out error types. Held-out results are not used for
   any selection.
3. Every test runs on all 40 ORIGINAL papers; their flags go to blind verification.
Recall is reported pooled and macro-averaged over error types (one held-out type is
over half of the held-out copies).

Usage: python -I results/run_final.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "climb"))
import score as S  # noqa: E402

OUT = ROOT / "results"


def macro(rbt: dict) -> float:
    vals = [int(h) / int(n) for h, n in (v.split("/") for v in rbt.values()) if int(n)]
    return round(sum(vals) / len(vals), 4) if vals else 0.0


def main() -> None:
    summaries = {p.parent.name: json.loads(p.read_text()) for p in (ROOT / "climb" / "runs").glob("*/summary.json")}
    dev_scores = {}
    for name in summaries:
        m = S.score(ROOT / "climb" / "runs" / name / "best.py", "visible", "dev", 0)
        dev_scores[name] = (m["score"], -m["fp_papers"])
    winner = max(dev_scores, key=lambda k: dev_scores[k])

    tests = {"grim": ROOT / "baselines" / "grim_check.py",
             "statcheck_py": ROOT / "baselines" / "statcheck_py.py"}
    tests.update({f"climber_{n}": ROOT / "climb" / "runs" / n / "best.py" for n in sorted(summaries)})

    all_ids = sorted(p.stem for p in S.PARSED.glob("*.json"))
    table, flags_out = {}, []
    for name, path in tests.items():
        row = {"path": str(path.relative_to(ROOT))}
        for split in ("visible", "heldout"):
            m = S.score(path, split, "test", 0)
            row[split] = {"recall": m["recall"], "macro_recall": macro(m["recall_by_type"]),
                          "recall_by_type": m["recall_by_type"], "n_planted": m["n_planted"]}
        mod = S.load_module(path)
        n_app, n_flag_papers = 0, 0
        for pid in all_ids:
            paper = json.loads((S.PARSED / f"{pid}.json").read_text())
            fl, _ = S.run_check(mod, paper)
            a = S.run_applicable(mod, paper)
            n_app += bool(a or (a is None and fl))
            n_flag_papers += bool(fl)
            for i, f in enumerate(fl):
                flags_out.append({"flag_id": f"{name}|{pid}|{i}", "test": name, "paper_id": pid,
                                  "location": f["location"], "message": str(f.get("message", ""))[:500]})
        row["original_40"] = {"applicable_papers": n_app, "papers_flagged": n_flag_papers,
                              "flags": sum(1 for x in flags_out if x["test"] == name)}
        table[name] = row
        print(name, json.dumps({k: v for k, v in row.items() if k != "path"}))

    (OUT / "final_scores.json").write_text(json.dumps(
        {"winner": winner, "dev_visible_scores": dev_scores, "tests": table}, indent=1))
    with (OUT / "real_flags.jsonl").open("w") as f:
        for x in flags_out:
            f.write(json.dumps(x) + "\n")
    print("winner:", winner, "| real flags to verify:", len(flags_out))


if __name__ == "__main__":
    main()
