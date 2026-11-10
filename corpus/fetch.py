"""Fetch candidate PLOS ONE political-science papers (JATS XML).

Usage: python3 -I corpus/fetch.py [--n 250] [--seed 20261110]

1. Queries the PLOS Search API (Solr) for all matching DOIs (paged).
2. Sorts DOIs lexicographically, samples N with a fixed seed.
3. Writes corpus/candidates.tsv and downloads JATS XML to corpus/raw/<paper_id>.xml
   (<= 4 concurrent requests, skips existing files, retries once).
"""
import argparse
import json
import random
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
API = "https://api.plos.org/search"
PARAMS = [
    ("q", 'subject:"Political science"'),
    ("fq", "doc_type:full"),
    ("fq", "journal_key:PLoSONE"),
    ("fq", "publication_date:[2015-01-01T00:00:00Z TO *]"),
    ("fl", "id,title_display,publication_date"),
    ("sort", "id asc"),
    ("wt", "json"),
]
XML_URL = "https://journals.plos.org/plosone/article/file?id={doi}&type=manuscript"
UA = {"User-Agent": "amazon-pilot-research/0.1 (academic corpus build)"}


def get(url, timeout=60):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def query_all(page=1000):
    docs, start = [], 0
    while True:
        qs = urllib.parse.urlencode(PARAMS + [("rows", page), ("start", start)])
        data = json.loads(get(f"{API}?{qs}"))["response"]
        docs.extend(data["docs"])
        start += page
        if start >= data["numFound"]:
            return docs, data["numFound"]
        time.sleep(1)


def paper_id(doi):
    return doi.replace("/", "_")


def download(doi):
    out = RAW / f"{paper_id(doi)}.xml"
    if out.exists() and out.stat().st_size > 0:
        return doi, "skip"
    for attempt in range(2):
        try:
            data = get(XML_URL.format(doi=urllib.parse.quote(doi, safe="/")))
            if b"<article" not in data[:5000]:
                raise ValueError("not JATS")
            out.write_bytes(data)
            return doi, "ok"
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(2)
    return doi, f"fail: {err}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=250)
    ap.add_argument("--seed", type=int, default=20261110)
    args = ap.parse_args()

    docs, found = query_all()
    # only research articles' DOIs (exclude e.g. corrections with suffixes)
    docs = sorted({d["id"]: d for d in docs if d["id"].count("/") == 1}.values(), key=lambda d: d["id"])
    print(f"numFound={found} unique={len(docs)}")
    sample = sorted(random.Random(args.seed).sample(docs, min(args.n, len(docs))), key=lambda d: d["id"])
    with open(HERE / "candidates.tsv", "w") as f:
        f.write("paper_id\tdoi\tpublication_date\ttitle\n")
        for d in sample:
            title = (d.get("title_display") or "").replace("\t", " ").replace("\n", " ")
            f.write(f"{paper_id(d['id'])}\t{d['id']}\t{d['publication_date'][:10]}\t{title}\n")

    RAW.mkdir(exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as ex:
        results = list(ex.map(download, [d["id"] for d in sample]))
    fails = [r for r in results if r[1].startswith("fail")]
    print(f"downloaded={sum(r[1]=='ok' for r in results)} skipped={sum(r[1]=='skip' for r in results)} failed={len(fails)}")
    for f in fails:
        print("  ", *f)


if __name__ == "__main__":
    main()
