# Pilot results (aggregate)

**Setup.** 40 PLOS ONE papers tagged "Political science" (2015+; the tag also covers public health and economics) with 262 results tables. 24 development papers, 16 test papers. Ten planted error kinds: six guided the search (175 altered copies), four were held back (187 copies). Three search runs with different starting directions; 44 rounds in total.

## Scores on the 16 unseen test papers

| | GRIM | statcheck (Python port) | Best single test (pre-registered) | Three tests together |
|---|---|---|---|---|
| Papers it can check (of 40) | 0 | 8 | 40 | 40 |
| Guiding error kinds caught | 0 of 71 | 0 of 71 | 28 of 71 (39%) | 47 of 71 (66%; 95% CI 55–76%) |
| Held-back error kinds caught | 0 of 78 | 9 of 78 (12%) | 12 of 78 (15%) | 23 of 78 (29%; 95% CI 21–40%) |

- "Three tests together" was chosen as a framing after seeing the results; the pre-registered rule picked one test.
- All 9 statcheck catches were of the held-back kind that mimics the errors statcheck was built for. On the other three held-back kinds it caught 0 of 69.
- No test caught the held-back kind "a number in the text doesn't match its table" (0 of 39).
- **Overfitting:** the best single test caught 58% of guiding-kind errors on development papers and 39% on unseen papers.

## What the discovered tests check

- **Run A (inferential):** recomputes significance from each coefficient and standard error under the table's own star legend; checks that estimates fall inside their confidence intervals.
- **Run B (arithmetic):** checks that counts, percentages, totals and sample sizes reconcile across tables and text.
- **Run C (open direction):** converged independently on the same principle as A: an estimate must agree with every other report of itself (SE, stars, interval, sign in the text).

## On the 40 original, unmodified papers

| | statcheck | Test A | Test B | Test C | A + B + C |
|---|---|---|---|---|---|
| Papers flagged | 4 | 6 | 4 | 6 | 11 |
| Papers with ≥1 model-confirmed inconsistency | 3 | 6 | 2 | 4 | 9 |
| Flags judged false alarms | 1 of 5 | 0 of 38 | 8 of 14 | 2 of 20 | 10 of 72 |

12 of the 40 papers have at least one model-confirmed inconsistency. The most common kind was significance stars that contradict the table's own stated legend; others were percentages that do not add up and an interval that does not contain its estimate. These look like typos or convention slips. A human spot check is in progress. Paper-level details are not published here (see README).

## Cost

- Search: $22.84 for 44 rounds, about **$0.52 per round** (Claude Opus 5.5, adaptive thinking, high effort, prompt caching).
- Verification: $7.50 for 77 flags.
- Running a discovered test: deterministic Python, effectively free per paper.

## Limitations

- The answer key and the search use the same model family, so planted-error recall measures generalization within shared assumptions about what errors look like.
- Confirmation of real flags is by a model; human review is pending.
- The pilot uses text and tables only, not deposited data and code.
- Small samples: 40 papers, 9–39 altered copies per held-back kind.
