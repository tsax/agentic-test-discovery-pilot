# Shared contract (all tracks build against this)

Clean-room rule: this repo contains no code from, and was not written by consulting,
`pdf-table-extractor` or any DNH project code. Sources are public JATS XML from PLOS ONE.

## 1. Parsed paper — `corpus/parsed/<paper_id>.json`

`paper_id` = DOI with `/` replaced by `_`, e.g. `10.1371_journal.pone.0273951`.

```json
{
  "paper_id": "10.1371_journal.pone.0273951",
  "doi": "10.1371/journal.pone.0273951",
  "title": "...",
  "year": 2022,
  "data_availability": "free text of the data availability statement",
  "paragraphs": [
    {"para_id": "p0", "section": "Results", "text": "..."}
  ],
  "tables": [
    {
      "table_id": "t1",
      "label": "Table 1",
      "caption": "...",
      "footnotes": "free text of table-wrap-foot, incl. star legend like '* p<0.05'",
      "header_rows": [["", "Model 1", "Model 2"]],
      "rows": [["Income", "0.123***", "0.101**"], ["", "(0.021)", "(0.030)"]]
    }
  ]
}
```

Rules:
- Cell text is the visible text, whitespace-collapsed, with significance stars and
  daggers kept exactly as printed (`***`, `†`). Superscripts are inlined as text.
- `rows` excludes header rows. Cells spanning multiple columns are repeated
  across the spanned columns so every row has the same length as the header.
- A table that cannot be parsed into a rectangular grid is dropped and logged.

## 2. Planted copy — `planted/<split>/<paper_id>__<n>.json` + `planted/<split>/manifest.jsonl`

A planted copy is a parsed-paper JSON with exactly ONE error inserted. Its
`paper_id` is `<orig_id>__<n>`. Each manifest line:

```json
{"planted_id": "10.1371_journal.pone.0273951__3", "orig_id": "10.1371_journal.pone.0273951",
 "error_type": "star_flip", "location": {"kind": "cell", "table_id": "t2", "row": 4, "col": 2},
 "before": "0.123***", "after": "0.123"}
```

`location` is either `{"kind":"cell","table_id","row","col"}` or
`{"kind":"text","para_id","start","end"}` (character offsets into the NEW paragraph text).

`<split>` is `visible` or `heldout`. The two splits use DISJOINT error types.

## 3. A candidate test — a Python module exposing

```python
def check(paper: dict) -> list[dict]:
    """Return flags. Must be deterministic, offline, < 2 s per paper, stdlib + numpy/scipy only."""
```

Each flag:
```json
{"location": {"kind": "cell", "table_id": "t2", "row": 4, "col": 2}, "message": "why"}
```
A test may also return `[]` (no flag) and should abstain (return `[]`) when its inputs are missing.

## 4. Scoring (`climb/score.py`)

A flag HITS a planted error if it is on the same paper and:
- cell location: same `table_id` and same `row` (column may differ, because a
  coefficient/SE pair spans two rows the flag on either row `row` or `row±1` counts), or
- text location: same `para_id` and overlapping character span.

Per test, on a split:
- `recall` = planted errors hit / planted errors in split
- `fp_papers` = number of ORIGINAL (unmodified) papers with ≥1 flag
- `applicable` = original papers where the test returned anything or declared applicability
- `score` = recall − 0.5 × (fp_papers / n_original_papers)

Unmodified originals are NOT guaranteed error-free; flags on them are reported, then
blind-verified, never assumed false.
