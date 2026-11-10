"""EXPLORATORY (not pre-registered): which planted errors does a union of tests catch?"""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "climb")); import score as S

def hits(path, split):
    mod = S.load_module(path); ids = set(S.paper_split()["test"])
    orig = {i: S.run_check(mod, json.loads((S.PARSED/f"{i}.json").read_text()))[0] for i in ids}
    out = set()
    for l in (S.PLANTED/split/"manifest.jsonl").read_text().splitlines():
        m = json.loads(l)
        if m["orig_id"] not in ids: continue
        fl, _ = S.run_check(mod, json.loads((S.PLANTED/split/f"{m['planted_id']}.json").read_text()))
        if any(S.loc_matches(f["location"], m["location"]) for f in fl) and not any(S.loc_matches(f["location"], m["location"]) for f in orig[m["orig_id"]]):
            out.add((m["planted_id"], m["error_type"]))
    return out

T = {"statcheck": ROOT/"baselines/statcheck_py.py", "A": ROOT/"climb/runs/A/best.py", "B": ROOT/"climb/runs/B/best.py", "C": ROOT/"climb/runs/C/best.py"}
res = {}
for split in ("visible", "heldout"):
    H = {k: hits(p, split) for k, p in T.items()}
    n = sum(1 for l in (S.PLANTED/split/"manifest.jsonl").read_text().splitlines() if json.loads(l)["orig_id"] in S.paper_split()["test"])
    for combo in (["B","A"], ["B","A","C"], ["statcheck","B","A","C"]):
        u = set().union(*(H[c] for c in combo))
        res[f"{split}:{'+'.join(combo)}"] = f"{len(u)}/{n} = {len(u)/n:.2f}"
print(json.dumps(res, indent=1)); (ROOT/"results/union_exploratory.json").write_text(json.dumps(res, indent=1))
