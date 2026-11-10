"""Parse PLOS ONE JATS XML into the corpus schema (docs/SCHEMA.md section 1). Stdlib only.

Usage: python3 -I corpus/parse_jats.py [raw_dir] [out_dir] [text_dir]
Defaults: corpus/raw -> corpus/parsed (+ corpus/text/<paper_id>.txt).
Dropped (non-rectangular / unparseable) tables are logged to corpus/parse_log.tsv.
"""
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
MML = "{http://www.w3.org/1998/Math/MathML}"
WS = re.compile(r"\s+")

# Elements whose text must not leak into flattened paragraph / cell text.
SKIP_TAGS = {"table-wrap", "fig", "tex-math", "graphic", "inline-graphic",
             "object-id", "disp-formula-group", "supplementary-material"}
# Block-ish elements: put a space around them so words don't glue together.
BREAK_TAGS = {"break", "p", "list-item", "title", "fn", "label"}


def local(tag):
    return tag.split("}", 1)[-1] if isinstance(tag, str) else ""


def flatten(el, skip=SKIP_TAGS):
    """Visible text of el with inline markup (bold/italic/sup/sub/xref/MathML) flattened."""
    parts = []

    def walk(e):
        name = local(e.tag)
        if name in skip:
            return
        if name == "alternatives":
            # prefer MathML / table over graphic/tex alternatives
            kids = [k for k in e if local(k.tag) not in skip]
            if kids:
                walk(kids[0])
            return
        if name in BREAK_TAGS:
            parts.append(" ")
        if e.text:
            parts.append(e.text)
        for k in e:
            walk(k)
            if k.tail:
                parts.append(k.tail)
        if name in BREAK_TAGS:
            parts.append(" ")

    walk(el)
    return WS.sub(" ", "".join(parts)).strip()


def norm_minus(s):
    """Helper for downstream numeric use ONLY. Stored text keeps U+2212 as printed."""
    return s.replace("−", "-").replace("–", "-")


# ---------------------------------------------------------------- metadata

def get_meta(root):
    am = root.find("./front/article-meta")
    doi = am.findtext("./article-id[@pub-id-type='doi']").strip()
    tg = am.find("./title-group/article-title")
    title = flatten(tg) if tg is not None else ""
    year = None
    for pd in am.findall("./pub-date"):
        y = pd.findtext("year")
        if y and (pd.get("pub-type") in ("epub", None) or pd.get("date-type") == "pub"):
            year = int(y)
            break
    if year is None and am.find("./pub-date/year") is not None:
        year = int(am.findtext("./pub-date/year"))
    da = ""
    for cm in am.iter("custom-meta"):
        if (cm.findtext("meta-name") or "").strip().lower() == "data availability":
            mv = cm.find("meta-value")
            da = flatten(mv) if mv is not None else ""
            break
    return doi, title, year, da


# ---------------------------------------------------------------- paragraphs

def collect_paragraphs(body):
    out = []

    def walk(e, section):
        for k in e:
            name = local(k.tag)
            if name == "sec":
                sec_title = section
                if section is None:  # top-level section title is the label we keep
                    t = k.find("title")
                    sec_title = flatten(t) if t is not None else ""
                walk(k, sec_title)
            elif name == "p":
                txt = flatten(k)
                if txt:
                    out.append({"para_id": f"p{len(out)}", "section": section or "", "text": txt})
                # paragraphs nested inside paragraphs (e.g. list items) are already flattened into txt
            elif name in ("list", "boxed-text", "disp-quote"):
                walk(k, section)
            elif name == "list-item":
                walk(k, section)
            # skip table-wrap, fig, disp-formula, etc.

    walk(body, None)
    return out


# ---------------------------------------------------------------- tables

def expand_rows(trs):
    """Expand colspan (repeat text) and rowspan (fill continued rows with '') into a grid.
    Returns list of rows (lists of str)."""
    grid = []
    pending = {}  # (row_idx, col_idx) -> text placeholder for rowspan continuation
    for r, tr in enumerate(trs):
        row = []
        col = 0
        cells = [c for c in tr if local(c.tag) in ("td", "th")]
        ci = 0
        while ci < len(cells) or (r, col) in pending:
            if (r, col) in pending:
                row.append(pending.pop((r, col)))
                col += 1
                continue
            c = cells[ci]
            ci += 1
            txt = flatten(c)
            cs = max(1, int(c.get("colspan") or 1))
            rs = max(1, int(c.get("rowspan") or 1))
            for j in range(cs):
                row.append(txt)
                for i in range(1, rs):
                    pending[(r + i, col + j)] = ""
            col += cs
        stray = [k for k in pending if k[0] == r]
        if stray:  # rowspan continuation beyond the end of this row's cells
            return grid, True
        grid.append(row)
    leftover = [k for k in pending if k[0] >= len(trs)]
    return grid, bool(leftover)


def parse_table(tw, idx):
    """Return (table_dict, None) or (None, reason)."""
    label = (tw.findtext("label") or "").strip()
    cap = tw.find("caption")
    caption = flatten(cap) if cap is not None else ""
    foot = tw.find("table-wrap-foot")
    footnotes = flatten(foot) if foot is not None else ""
    tables = [t for t in tw.iter("table")]
    if not tables:
        return None, "no <table> (image only)"
    if len(tables) > 1:
        return None, f"{len(tables)} <table> elements in one table-wrap"
    tbl = tables[0]
    head_trs = [tr for th in tbl.findall("thead") for tr in th.findall("tr")]
    body_trs = [tr for tb in tbl.findall("tbody") for tr in tb.findall("tr")]
    if not body_trs:
        body_trs = tbl.findall("tr")
    if not body_trs:
        return None, "no body rows"
    grid, overflow = expand_rows(head_trs + body_trs)
    widths = {len(r) for r in grid}
    if len(widths) != 1:
        return None, f"non-rectangular widths={sorted(widths)}"
    if overflow:
        return None, "rowspan overflows table"
    nh = len(head_trs)
    return {
        "table_id": f"t{idx}",
        "label": label,
        "caption": caption,
        "footnotes": footnotes,
        "header_rows": grid[:nh],
        "rows": grid[nh:],
    }, None


# ---------------------------------------------------------------- main

def parse_file(path):
    root = ET.parse(path).getroot()
    doi, title, year, da = get_meta(root)
    body = root.find("body")
    paragraphs = collect_paragraphs(body) if body is not None else []
    tables, dropped = [], []
    tws = list(body.iter("table-wrap")) if body is not None else []
    for n, tw in enumerate(tws, 1):
        t, why = parse_table(tw, n)
        if t is None:
            dropped.append((f"t{n}", (tw.findtext("label") or "").strip(), why))
        else:
            tables.append(t)
    paper = {
        "paper_id": doi.replace("/", "_"),
        "doi": doi,
        "title": title,
        "year": year,
        "data_availability": da,
        "paragraphs": paragraphs,
        "tables": tables,
    }
    return paper, dropped


def main():
    raw = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "raw"
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "parsed"
    txt = Path(sys.argv[3]) if len(sys.argv) > 3 else HERE / "text"
    out.mkdir(parents=True, exist_ok=True)
    txt.mkdir(parents=True, exist_ok=True)
    log = ["paper_id\ttable_id\tlabel\treason"]
    n_ok = n_tab = n_drop = 0
    for p in sorted(raw.glob("*.xml")):
        try:
            paper, dropped = parse_file(p)
        except Exception as e:  # noqa: BLE001
            log.append(f"{p.stem}\t-\t-\tPAPER PARSE ERROR: {e}")
            continue
        pid = paper["paper_id"]
        (out / f"{pid}.json").write_text(json.dumps(paper, ensure_ascii=False, indent=1))
        (txt / f"{pid}.txt").write_text("\n\n".join(x["text"] for x in paper["paragraphs"]) + "\n")
        for d in dropped:
            log.append(f"{pid}\t" + "\t".join(d))
        n_ok += 1
        n_tab += len(paper["tables"])
        n_drop += len(dropped)
    (HERE / "parse_log.tsv").write_text("\n".join(log) + "\n")
    print(f"papers={n_ok} tables_kept={n_tab} tables_dropped={n_drop}")


if __name__ == "__main__":
    main()
