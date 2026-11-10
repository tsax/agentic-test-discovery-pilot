"""Black-box scorer for candidate tests (see docs/SCHEMA.md §3-4).

Usage:
  python -I climb/score.py <test_module.py> --split visible --papers dev [--feedback N] [--out metrics.json]

Runs the candidate's check() on (a) every ORIGINAL paper in the paper set and
(b) every planted copy of those papers in the chosen split. Each paper runs
under a per-paper time limit; exceptions count as "no flags" and are reported.

A planted error is HIT only if the test raises a flag at the planted location
on the planted copy that it did NOT raise at that location on the original
(so pre-existing flags cannot be credited as detections).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import random
import signal
import sys
import traceback
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PARSED = ROOT / "corpus" / "parsed"
PLANTED = ROOT / "planted"
SPLIT_FILE = ROOT / "climb" / "paper_split.json"
PER_PAPER_SECONDS = 2


class _Timeout(Exception):
    pass


def _alarm(signum, frame):  # noqa: ARG001
    raise _Timeout()


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location(f"candidate_{abs(hash(str(path)))}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    if not hasattr(mod, "check"):
        raise AttributeError("module has no check(paper) function")
    return mod


def paper_split() -> dict[str, list[str]]:
    """Deterministic 60/40 dev/test split of original papers, created once."""
    if SPLIT_FILE.exists():
        return json.loads(SPLIT_FILE.read_text())
    ids = sorted(p.stem for p in PARSED.glob("*.json"))
    rng = random.Random(20261110)
    rng.shuffle(ids)
    k = round(len(ids) * 0.6)
    split = {"dev": sorted(ids[:k]), "test": sorted(ids[k:])}
    SPLIT_FILE.write_text(json.dumps(split, indent=1))
    return split


def run_check(mod, paper: dict) -> tuple[list[dict], str | None]:
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(PER_PAPER_SECONDS)
    try:
        flags = mod.check(paper) or []
        if not isinstance(flags, list):
            return [], "check() did not return a list"
        return [f for f in flags if isinstance(f, dict) and isinstance(f.get("location"), dict)], None
    except _Timeout:
        return [], f"timeout > {PER_PAPER_SECONDS}s"
    except Exception:  # noqa: BLE001
        return [], traceback.format_exc(limit=3)
    finally:
        signal.alarm(0)


def run_applicable(mod, paper: dict) -> bool | None:
    fn = getattr(mod, "applicable", None)
    if fn is None:
        return None
    try:
        return bool(fn(paper))
    except Exception:  # noqa: BLE001
        return False


def loc_matches(flag_loc: dict, planted_loc: dict) -> bool:
    if flag_loc.get("kind") != planted_loc.get("kind"):
        return False
    if planted_loc["kind"] == "cell":
        return (
            flag_loc.get("table_id") == planted_loc["table_id"]
            and isinstance(flag_loc.get("row"), int)
            and abs(flag_loc["row"] - planted_loc["row"]) <= 1
        )
    if planted_loc["kind"] == "text":
        if flag_loc.get("para_id") != planted_loc["para_id"]:
            return False
        s, e = flag_loc.get("start"), flag_loc.get("end")
        if not isinstance(s, int) or not isinstance(e, int):
            return True  # paragraph-level flag counts for that paragraph
        return s < planted_loc["end"] and planted_loc["start"] < e
    return False


def row_context(paper: dict, loc: dict) -> dict:
    """Small excerpt around a location, for feedback to the searching agent."""
    if loc.get("kind") == "cell":
        t = next((t for t in paper["tables"] if t["table_id"] == loc.get("table_id")), None)
        if not t:
            return {}
        r = loc.get("row", 0)
        rows = t["rows"][max(0, r - 1): r + 2]
        return {"table_label": t.get("label"), "caption": (t.get("caption") or "")[:200],
                "header": t.get("header_rows", [])[-1:] , "rows_around": rows,
                "footnotes": (t.get("footnotes") or "")[:300]}
    if loc.get("kind") == "text":
        p = next((p for p in paper["paragraphs"] if p["para_id"] == loc.get("para_id")), None)
        if not p:
            return {}
        s = max(0, loc.get("start", 0) - 150)
        return {"text_excerpt": p["text"][s: loc.get("end", s) + 150]}
    return {}


def score(test_path: Path, split_name: str, paper_set: str, n_feedback: int = 6) -> dict:
    mod = load_module(test_path)
    ids = paper_split()[paper_set]
    originals = {i: json.loads((PARSED / f"{i}.json").read_text()) for i in ids}

    errors: list[str] = []
    orig_flags: dict[str, list[dict]] = {}
    applicable = 0
    for pid, paper in originals.items():
        flags, err = run_check(mod, paper)
        if err:
            errors.append(f"{pid}: {err}")
        orig_flags[pid] = flags
        a = run_applicable(mod, paper)
        if a or (a is None and flags):
            applicable += 1

    manifest = [json.loads(l) for l in (PLANTED / split_name / "manifest.jsonl").read_text().splitlines() if l.strip()]
    manifest = [m for m in manifest if m["orig_id"] in originals]

    hits_by_type: Counter = Counter()
    total_by_type: Counter = Counter()
    missed, hit_examples = [], []
    for m in manifest:
        copy = json.loads((PLANTED / split_name / f"{m['planted_id']}.json").read_text())
        flags, err = run_check(mod, copy)
        if err:
            errors.append(f"{m['planted_id']}: {err}")
        pre = [f for f in orig_flags[m["orig_id"]] if loc_matches(f["location"], m["location"])]
        new = [f for f in flags if loc_matches(f["location"], m["location"])]
        hit = bool(new) and not pre
        total_by_type[m["error_type"]] += 1
        if hit:
            hits_by_type[m["error_type"]] += 1
            hit_examples.append(m)
        else:
            missed.append((m, copy))

    n_planted = len(manifest)
    recall = sum(hits_by_type.values()) / n_planted if n_planted else 0.0
    fp_papers = sum(1 for f in orig_flags.values() if f)
    n_orig = len(originals)
    result = {
        "test": str(test_path.relative_to(ROOT)) if test_path.is_relative_to(ROOT) else str(test_path),
        "split": split_name, "paper_set": paper_set,
        "n_original_papers": n_orig, "n_planted": n_planted,
        "recall": round(recall, 4),
        "fp_papers": fp_papers, "fp_flags_total": sum(len(f) for f in orig_flags.values()),
        "applicable_papers": applicable,
        "score": round(recall - 0.5 * (fp_papers / n_orig if n_orig else 0), 4),
        "recall_by_type": {t: f"{hits_by_type[t]}/{total_by_type[t]}" for t in sorted(total_by_type)},
        "n_errors": len(errors), "errors_sample": errors[:3],
    }
    if n_feedback:
        rng = random.Random(0)
        sample = rng.sample(missed, min(n_feedback, len(missed)))
        result["missed_examples"] = [
            {"error_type": m["error_type"], "location": m["location"], "before": m["before"],
             "after": m["after"], "context_in_copy": row_context(c, m["location"])}
            for m, c in sample
        ]
        fps = [(pid, f) for pid, fl in orig_flags.items() for f in fl]
        rng.shuffle(fps)
        result["flags_on_original_papers"] = [
            {"paper_id": pid, "flag": f, "context": row_context(originals[pid], f["location"])}
            for pid, f in fps[:n_feedback]
        ]
    result["_orig_flags"] = orig_flags  # stripped before showing to the agent
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("test", type=Path)
    ap.add_argument("--split", default="visible", choices=["visible", "heldout"])
    ap.add_argument("--papers", default="dev", choices=["dev", "test"])
    ap.add_argument("--feedback", type=int, default=6)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    res = score(a.test.resolve(), a.split, a.papers, a.feedback)
    if a.out:
        a.out.write_text(json.dumps(res, indent=1, default=str))
    res.pop("_orig_flags", None)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
