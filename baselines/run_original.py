"""Run both baselines on every original paper in corpus/parsed/ and write
baselines/results_original.json. Usage: python3 baselines/run_original.py [--no-statcheck]"""
import collections
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import grim_check  # noqa: E402


def main():
    files = sorted(glob.glob(os.path.join(ROOT, "corpus", "parsed", "*.json")))
    papers = [json.load(open(f, encoding="utf-8")) for f in files]
    papers = [p for p in papers if "__" not in p["paper_id"]]  # originals only
    out = {"grim": {}, "statcheck": {}}

    skip_reasons = collections.Counter()
    for p in papers:
        an = grim_check.analyze(p)
        for r in an:
            if r["status"] == "skipped":
                skip_reasons[r["reason"].split(":")[0].split("(")[0].strip()] += 1
        out["grim"][p["paper_id"]] = {
            "applicable": any(r["status"] in ("consistent", "inconsistent") for r in an),
            "flags": grim_check.check(p),
            "n_candidates": len(an),
            "n_tested": sum(r["status"] != "skipped" for r in an),
        }

    if "--no-statcheck" not in sys.argv:
        import statcheck_check
        res = statcheck_check.run_batch_full(papers)
        for pid, e in res.items():
            out["statcheck"][pid] = {"applicable": e["applicable"], "flags": e["flags"],
                                     "n_extracted": e["n_extracted"]}

    with open(os.path.join(HERE, "results_original.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    print(f"papers: {len(papers)}")
    for name, d in out.items():
        if not d:
            print(f"{name}: not run")
            continue
        app = sum(v["applicable"] for v in d.values())
        nfl = sum(len(v["flags"]) for v in d.values())
        fpp = sum(bool(v["flags"]) for v in d.values())
        print(f"{name}: applicable {app}/{len(d)}, total flags {nfl}, papers with >=1 flag {fpp}")
    print("grim candidates:", sum(v["n_candidates"] for v in out["grim"].values()),
          "tested:", sum(v["n_tested"] for v in out["grim"].values()))
    print("grim skip reasons:", dict(skip_reasons.most_common()))


if __name__ == "__main__":
    main()
