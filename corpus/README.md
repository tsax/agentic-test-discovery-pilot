# Corpus: PLOS ONE political-science papers (pilot)

Final corpus: **40 papers** in `parsed/` (schema: `../docs/SCHEMA.md` section 1), with
plain body text in `text/<paper_id>.txt` (paragraphs joined by blank lines, for statcheck).
Publication years 2016-2026. All 40 have a data availability statement.

| | count |
|---|---|
| candidates queried / sampled / downloaded | 3,126 / 250 / 250 |
| candidates parsed | 250 (0 paper-level failures) |
| eligible (>= 1 results-like table, not retracted, not formal-only) | 182 |
| kept | 40 |
| tables in kept papers | 262 |
|  - regression-like | 130 (23 papers have a classic "coef*** / (SE)" table) |
|  - descriptive (means/SD/N/%, correlation matrices) | 81 |
|  - other (lists, fit indices, loadings, text tables) | 51 |
| paragraphs in kept papers | 1,950 |
| table cells containing U+2212 (kept as printed) | 461 |

## Pipeline (all stdlib, `python3 -I`)

1. `fetch.py` - query (see `QUERY.md`), sample 250, download JATS to `raw/`.
2. `parse_jats.py` - `raw/*.xml` -> `parsed/*.json` + `text/*.txt`; dropped tables -> `parse_log.tsv`.
3. `select_papers.py` - classify tables, write `selection.csv`, move non-kept papers to
   `parsed_dropped/` (text to `text_dropped/`). Idempotent: it first moves dropped files back.
4. `sanity_check.py` - prints parsed vs raw rows for 3 random kept papers and checks that
   `*`, `(`, `)`, U+2212, `-`, dagger counts survive. Passed for seeds 1-5 and 7.

## Selection rule (`select_papers.py`)

A table is **regression** if >= 30% of non-empty body cells are numeric and any of:
>= 2 starred numbers plus (>= 2 parenthesised numbers, or a coefficient/SE/t/p header, or a
regression-ish caption); >= 3 parenthesised numbers plus regression caption/coef header; or a
coefficient header together with an SE/t/z/p header. Correlation tables (caption says
"correlat" without regression words) are counted as **descriptive**. A table is
**descriptive** if mostly numeric (counting "n (%)" and "mean ± sd" cells) and caption /
header / stub mention mean, SD, N, %, frequency, descriptive, etc.

A paper is eligible if it has >= 1 results-like table, its title is not "RETRACTED"/"Expression
of concern", and it is not a formal-model/simulation paper without a regression table
(title or >= 8 body mentions of simulation/agent-based/game-theory/equilibrium).
Ranking: papers with a regression table first, then score = 10 x [has classic stars+(SE)
table] + min(#APA-style in-text stats like `t(23) =`, `p < .05`, 10), ties by seeded random.
Top 40 kept. Per-paper reasons are in `selection.csv`.

## Parsing conventions and known limits

- Cell/paragraph text = `itertext` of the element with inline markup flattened and
  whitespace collapsed; stars, daggers, superscripts are inlined as printed
  (`0.123<sup>***</sup>` -> `0.123***`). U+2212 is NOT normalised in stored text; use
  `parse_jats.norm_minus()` (U+2212 and U+2013 -> `-`) when converting to numbers.
- `colspan` cells are repeated across spanned columns (per schema). **`rowspan`
  continuation cells are filled with `""`** (not repeated) so numbers are never duplicated;
  a row label spanning two rows therefore appears only on the first row.
- `header_rows` = `<thead>` rows only. 13 kept tables have no `<thead>`, so their header
  appears as the first entries of `rows`.
- `table_id` = position of the `table-wrap` in the body (`t1`, `t2`, ...); ids of dropped
  tables are skipped, not renumbered (e.g. 10.1371_journal.pone.0284692 has no `t2`).
- Dropped tables (whole candidate set, 2 of 1,099): one table-wrap with 2 `<table>`
  elements (pone.0284692 t2, a kept paper), one non-rectangular grid (pone.0355484 t4).
  Image-only tables would also be dropped (none seen). Supplementary tables are not parsed.
- Inline math: MathML is flattened to its text (e.g. `R2`, `χ2`); `<tex-math>` and graphics
  are skipped, so some formula symbols render oddly. Display formulas (`disp-formula`
  outside `<p>`) are excluded from paragraphs.
- `section` = top-level `<sec>` title (e.g. "Results"), not the subsection. Paragraphs in
  lists/boxed-text are included; captions, table footnotes and figure text are not.
- Citations (`<xref>`) stay inline as printed (e.g. `[12]`).
- The PLOS "Political science" subject tag is broad: the corpus includes public-health,
  management, and economics papers that carry the tag, not only classic political science.
