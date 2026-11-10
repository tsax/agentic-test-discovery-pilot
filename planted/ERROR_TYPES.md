# Planted error catalog (frozen)

Each planted copy contains exactly ONE error. Every error is detectable in principle from the
paper's own tables and text, and is produced by a realistic slip (dropped/added star, adjacent-digit
transposition, off-by-one digit, decimal-point shift, sign drop, `<`/`>` swap). Originals are
required to be internally consistent at the planted site before mutation, and the mutation must
break that consistency by a clear margin (robust to the rounding of printed numbers).

Shared number mutations (`num_variants` in plant.py): `transpose` (swap two adjacent different
digits), `digit` (one digit +/-1), `shift` (move the decimal point one place). Leading zeros are
never created; sign, comma grouping and `.5`-style formatting are preserved.

Locations: cell = `rows` index (0-based, header rows excluded) and column (0-based). Text = offsets
into the NEW paragraph text, covering the inconsistent expression.

| # | name | split | where | difficulty |
|---|------|-------|-------|------------|
| 1 | star_flip | visible | cell | easy |
| 2 | se_misreport | heldout | cell | medium |
| 3 | sign_flip | visible | cell | hard |
| 4 | n_text_mismatch | visible | text | medium |
| 5 | pct_sum_break | visible | cell | easy |
| 6 | mean_outside_range | heldout | cell | easy |
| 7 | ci_excludes_estimate | visible | cell or text | easy |
| 8 | statcheck_p_mismatch | heldout | text | medium |
| 9 | count_pct_mismatch | visible | cell | medium |
| 10 | text_table_transposition | heldout | text | hard |

## 1. star_flip
- **Change:** significance marks on one coefficient move one level (add or drop one level, e.g.
  `0.198**` -> `0.198***`, `-0.46` -> `-0.46*`), using the table's own legend symbols.
- **Where:** the cell carrying the marks (coefficient cell, or the SE cell if the table prints
  stars there).
- **Precondition:** a coefficient (with a decimal point) paired with a standard error, in one of
  three layouts: stacked (`coef` row then `(se)` row), same cell (`coef (se)`), or side-by-side
  (adjacent column headed SE / Std. Err.). The table's caption/footnotes contain a parsable star
  legend (`* p<0.05 ...`, `***p < .01`, or `*, **, *** ... 10%, 5%, 1%`). The table is not an
  odds/hazard/risk-ratio table and the footnotes do not say the parentheses hold t/z statistics,
  p-values, CIs or SDs. If SEs are not named explicitly, >= 3 pairs are needed. >= 70% of the
  table's pairs must be star-consistent, and the chosen pair's star level must be unambiguous across
  the rounding interval of coef and SE under both normal and t(40) reference distributions.

## 2. se_misreport
- **Change:** the standard error is mistyped (decimal shift, digit transposition or one-digit
  change) so that |coef/SE| changes by a factor >= 2 and the printed stars are now unambiguously
  wrong (e.g. `(0.120)` -> `(1.20)` under a `***` coefficient).
- **Where:** the SE cell (the row below the coefficient in stacked layout).
- **Precondition:** same as star_flip.

## 3. sign_flip
- **Change:** the minus sign of a regression coefficient is dropped or added in the table
  (`0.234***` -> `−0.234***`, using the paper's own minus character). Stars stay unchanged.
- **Where:** coefficient cell.
- **Precondition:** a coefficient/SE pair (any of the three layouts; not a ratio table, not a
  mean(SD) column) whose absolute value has >= 2 significant digits, appears in exactly one table
  cell of the paper, and is quoted in a paragraph with the same sign as the table. After the
  flip the text states the opposite direction from the table.

## 4. n_text_mismatch
- **Change:** a sample size quoted in the text gets an adjacent-digit transposition or a one-digit
  change (`1,842 respondents` -> `1,824 respondents`).
- **Where:** text, the number.
- **Precondition:** a table reports N (a row labelled N / Observations / Obs. / Sample size /
  Respondents ..., or `N = 1234` in a header/caption/footnote), N >= 50, and the same integer
  appears in a paragraph within 60 characters of a sample-size word (respondents, participants,
  observations, sample, N = ...). Years (1900-2035 without comma) are excluded unless written
  `N = `. The new value does not appear anywhere in the paper's tables.

## 5. pct_sum_break
- **Change:** one percentage in a set of shares that sums to 100 is mistyped so the set now sums
  to 100 +/- 3..40 (`38.4` -> `48.4`, `6.3` -> `0.63`).
- **Where:** cell.
- **Precondition:** >= 3 contiguous cells in a column or row that are percentages (cell contains
  `%`, or the column header / row label / caption / footnotes mention `%` or percent) and sum to
  100 within rounding tolerance (<= 1.5). Contiguous runs may be partitioned into several
  100-sum blocks; a trailing `100` total cell is ignored.

## 6. mean_outside_range
- **Change:** a reported mean is mistyped (usually a decimal shift) so it falls outside its own
  reported [min, max] by >= 5% of the range (`45.3` -> `4.53` with min 18, max 90).
- **Where:** the mean cell.
- **Precondition:** a descriptive-statistics table with exactly one Mean, one Min and one Max
  column (or row), with plain numeric cells satisfying min <= mean <= max, min < max.

## 7. ci_excludes_estimate
- **Change:** the point estimate or one CI bound is mistyped so the estimate lies outside its own
  confidence interval (`1.45, 95% CI [1.12, 1.88]` -> `0.45, 95% CI [1.12, 1.88]`), bounds still
  ordered.
- **Where:** cell (estimate or bound cell) or text (span = estimate through upper bound).
- **Precondition:** an estimate with CI that originally contains it: `est (lo–hi)` / `est [lo, hi]`
  cells in a table that mentions CI/interval; a CI column (`[lo, hi]`, `lo–hi`) or lower/upper
  column pair next to an estimate column; or text `est, 95% CI [lo, hi]`, `est (95% CI lo–hi)`,
  `est [lo, hi]` (all decimals if no CI keyword).

## 8. statcheck_p_mismatch
- **Change:** in a reported test (`t(df) = x, p = y`; also F(df1,df2), χ2(df), z, r(df)) either
  the p-value is mistyped, the inequality is flipped (`<` -> `>`), or the test statistic is
  mistyped, such that the reported p is grossly incompatible (factor >= 2) with the p implied by
  the statistic and df.
- **Where:** text, span = the whole `stat(df) = x, p op y` expression.
- **Precondition:** such an expression whose original p is consistent with the statistic over the
  rounding interval of the statistic.

## 9. count_pct_mismatch
- **Change:** in an `n (%)` cell the count or the percentage is mistyped (transposition / one
  digit) so the implied denominator differs by > 3% from the one implied by its neighbours
  (`94 (47.0)` -> `94 (74.0)` when the column's N is 200).
- **Where:** cell.
- **Precondition:** >= 3 contiguous `count (pct)` or `pct% (count)` cells in a column or row (with
  `%` in the cell or a %/percent/n (%) cue in the table) whose implied denominators 100*n/pct share
  a common integer N within rounding.

## 10. text_table_transposition
- **Change:** a decimal number quoted in the text that matches a table value gets two adjacent
  digits transposed (fallback: one digit +/-1) (`b = 0.234` -> `b = 0.324`).
- **Where:** text, the number.
- **Precondition:** a decimal number in a paragraph with >= 2 significant digits that appears
  verbatim among the numbers of the paper's tables, is not a p-value, and is not inside a CI or
  test-statistic expression. The mutated number does not appear in any table.
