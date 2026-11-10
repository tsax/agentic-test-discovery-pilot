"""Symmetric 'applicability': on how many of the 40 papers does each test catch >=1 planted error
(any split)? Same rule for every test, so self-declared applicable() functions don't matter."""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "climb")); import score as S
T = {"grim": ROOT/"baselines/grim_check.py", "statcheck": ROOT/"baselines/statcheck_py.py",
     "A": ROOT/"climb/runs/A/best.py", "B": ROOT/"climb/runs/B/best.py", "C": ROOT/"climb/runs/C/best.py"}
ids = sorted(p.stem for p in S.PARSED.glob("*.json"))
caught = {}
for name, path in T.items():
    mod = S.load_module(path)
    orig = {i: S.run_check(mod, json.loads((S.PARSED/f"{i}.json").read_text()))[0] for i in ids}
    papers = set()
    for split in ("visible", "heldout"):
        for l in (S.PLANTED/split/"manifest.jsonl").read_text().splitlines():
            m = json.loads(l)
            fl, _ = S.run_check(mod, json.loads((S.PLANTED/split/f"{m['planted_id']}.json").read_text()))
            if any(S.loc_matches(f["location"], m["location"]) for f in fl) and not any(S.loc_matches(f["location"], m["location"]) for f in orig[m["orig_id"]]):
                papers.add(m["orig_id"])
    caught[name] = sorted(papers)
planted_papers = {json.loads(l)["orig_id"] for s in ("visible","heldout") for l in (S.PLANTED/s/"manifest.jsonl").read_text().splitlines()}
res = {k: len(v) for k, v in caught.items()}
res["A+B+C"] = len(set(caught["A"]) | set(caught["B"]) | set(caught["C"]))
res["papers_with_any_planted_error"] = len(planted_papers)
print(json.dumps(res, indent=1)); (ROOT/"audit/engagement.json").write_text(json.dumps({"counts": res, "papers": caught}, indent=1))
