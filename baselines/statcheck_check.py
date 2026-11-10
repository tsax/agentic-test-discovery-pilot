"""statcheck baseline: wraps the real statcheck R package (Nuijten et al. 2016).

statcheck extracts APA-style NHST results (t, F, r, chi2, z, Q) from text, recomputes
the p-value from the test statistic and df, and marks `error` when the reported p is
inconsistent with the recomputed one.  A flag here = statcheck error == TRUE.

R is called ONCE per batch (`run_batch`) via baselines/run_statcheck.R.  `check` and
`applicable` call `run_batch` on a single paper; results are cached by paper content
so check+applicable on the same paper cost one R call.

Location mapping: statcheck gives no offsets, so we search for the matched raw string in
the paragraph text (exact, then a whitespace/dash/minus-tolerant regex); if not found,
the flag spans the whole paragraph.

Note: this baseline needs Rscript + statcheck, so it is not stdlib-only as SCHEMA §3
asks of candidate tests; it is a reference baseline.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
R_SCRIPT = os.path.join(HERE, "run_statcheck.R")
_CACHE: dict[str, dict] = {}


def _rscript() -> str:
    exe = shutil.which("Rscript") or "/opt/homebrew/bin/Rscript"
    if not os.path.exists(exe):
        raise RuntimeError("Rscript not found; install R and the statcheck package")
    return exe


def _clean(s: str) -> str:
    # Same length as input so offsets would still be valid; TSV-safe.
    return re.sub(r"[\t\r\n]", " ", s or "")


def _paper_key(paper: dict) -> str:
    h = hashlib.sha1()
    h.update(str(paper.get("paper_id")).encode())
    for p in paper.get("paragraphs") or []:
        h.update(b"\x1e" + str(p.get("para_id")).encode() + b"\x1f" + (p.get("text") or "").encode())
    return h.hexdigest()


_EQUIV = {
    "-": "[-−–—]", "−": "[-−–—]", "–": "[-−–—]",
    "=": "[=＝]", "<": "[<≤]", ">": "[>≥]", ".": "[.·]",
}


def locate(raw: str, text: str) -> tuple[int, int] | None:
    if not raw:
        return None
    i = text.find(raw)
    if i >= 0:
        return i, i + len(raw)
    parts = []
    for ch in raw:
        if ch.isspace():
            continue
        parts.append(_EQUIV.get(ch, re.escape(ch)))
    if not parts:
        return None
    pat = r"\s*".join(parts)
    m = re.search(pat, text, flags=re.I)
    if m:
        return m.start(), m.end()
    return None


def _truthy(v: str) -> bool:
    return str(v).strip().upper() in ("TRUE", "T", "1")


def run_batch_full(papers: list[dict]) -> dict[str, dict]:
    """paper_id -> {"applicable": bool, "n_extracted": int, "flags": [...], "results": [...]}."""
    out: dict[str, dict] = {}
    todo = []
    for p in papers:
        k = _paper_key(p)
        if k in _CACHE:
            out[p["paper_id"]] = _CACHE[k]
        else:
            todo.append((k, p))
    if not todo:
        return out
    with tempfile.TemporaryDirectory() as td:
        tsv = os.path.join(td, "in.tsv")
        csv_out = os.path.join(td, "out.csv")
        with open(tsv, "w", encoding="utf-8", newline="") as f:
            f.write("paper_id\tpara_id\ttext\n")
            for _, p in todo:
                for para in p.get("paragraphs") or []:
                    f.write(f"{_clean(p['paper_id'])}\t{_clean(para['para_id'])}\t{_clean(para.get('text'))}\n")
        proc = subprocess.run([_rscript(), R_SCRIPT, tsv, csv_out], capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"run_statcheck.R failed: {proc.stderr[-2000:]}")
        with open(csv_out, encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))
    by_paper: dict[str, list[dict]] = {}
    for r in rows:
        by_paper.setdefault(r["paper_id"], []).append(r)
    for k, p in todo:
        texts = {para["para_id"]: para.get("text") or "" for para in p.get("paragraphs") or []}
        res = by_paper.get(p["paper_id"], [])
        flags = []
        for r in res:
            if not _truthy(r.get("error")):
                continue
            text = texts.get(r["para_id"], "")
            span = locate(r.get("raw", ""), text) or (0, len(text))
            de = " (decision error)" if _truthy(r.get("decision_error")) else ""
            flags.append({
                "location": {"kind": "text", "para_id": r["para_id"], "start": span[0], "end": span[1]},
                "message": (f"statcheck: reported '{r.get('raw')}' but recomputed p = "
                            f"{r.get('computed_p')}{de}"),
            })
        entry = {"applicable": len(res) > 0, "n_extracted": len(res), "flags": flags, "results": res}
        _CACHE[k] = entry
        out[p["paper_id"]] = entry
    return out


def run_batch(papers: list[dict]) -> dict[str, list[dict]]:
    return {pid: e["flags"] for pid, e in run_batch_full(papers).items()}


def check(paper: dict) -> list[dict]:
    return run_batch_full([paper])[paper["paper_id"]]["flags"]


def applicable(paper: dict) -> bool:
    return run_batch_full([paper])[paper["paper_id"]]["applicable"]


if __name__ == "__main__":
    import sys
    papers = [json.load(open(fp)) for fp in sys.argv[1:]]
    print(json.dumps(run_batch_full(papers), indent=1, default=str)[:5000])
