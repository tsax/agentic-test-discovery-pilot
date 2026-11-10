"""Side-by-side check of parsed tables vs raw JATS for N random kept papers.

Usage: python3 -I corpus/sanity_check.py [--n 3] [--seed 7]
Prints, for one starred/parenthesised table per paper, the first rows of the parsed grid
next to the raw <tr> cell texts, and checks that every '*', '(', U+2212 and '†'
in the raw table text survives into the parsed cells.
"""
import argparse
import json
import random
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from parse_jats import flatten, local, norm_minus  # noqa: E402

MARKS = ["*", "(", ")", "−", "†", "‡", "-"]


def raw_table_text(xml_path, label):
    root = ET.parse(xml_path).getroot()
    for tw in root.find("body").iter("table-wrap"):
        if (tw.findtext("label") or "").strip() == label:
            tbl = next(tw.iter("table"))
            trs = [tr for tr in tbl.iter("tr")]
            return [[flatten(c) for c in tr if local(c.tag) in ("td", "th")] for tr in trs], tbl
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    papers = sorted((HERE / "parsed").glob("*.json"))
    ok_all = True
    for f in random.Random(args.seed).sample(papers, args.n):
        p = json.loads(f.read_text())
        tabs = sorted(p["tables"], key=lambda t: -sum(("*" in c) + ("(" in c) + ("−" in c)
                                                       for r in t["rows"] for c in r))
        t = tabs[0]
        raw_rows, tbl = raw_table_text(HERE / "raw" / f"{p['paper_id']}.xml", t["label"])
        print("=" * 100)
        print(p["paper_id"], t["table_id"], t["label"], "|", t["caption"][:80])
        print("footnotes:", t["footnotes"][:160])
        parsed_rows = t["header_rows"] + t["rows"]
        for i, (a, b) in enumerate(zip(parsed_rows[:8], raw_rows[:8])):
            print(f" parsed[{i}]: {a}")
            print(f"    raw[{i}]: {b}")
        raw_txt = "".join(tbl.itertext())
        parsed_txt = "".join(c for r in parsed_rows for c in r)
        for m in MARKS:
            rc, pc = raw_txt.count(m), parsed_txt.count(m)
            # parsed may exceed raw because colspan cells are repeated
            status = "OK" if pc >= rc else "LOST"
            ok_all &= status == "OK"
            print(f"   {m!r:10} raw={rc:4d} parsed={pc:4d} {status}")
        stored_u2212 = sum("−" in c for r in parsed_rows for c in r)
        normed = sum("−" in norm_minus(c) for r in parsed_rows for c in r)
        print(f"   cells with U+2212 stored={stored_u2212} after norm_minus()={normed}")
        if t["rows"] and len(t["rows"][0]) != len(parsed_rows[0]):
            ok_all = False
    print("ALL OK" if ok_all else "PROBLEMS FOUND")


if __name__ == "__main__":
    main()
