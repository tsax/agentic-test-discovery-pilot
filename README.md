# Agents that discover tests of empirical research: pilot

Can an AI agent invent new automatic checks for errors in published research papers, and can we tell whether those checks really work?

This repository is a one-day pilot (November 10, 2026) of that idea, run on 40 open-access social-science papers.

## What we did, in plain terms

1. **Collected papers.** 40 PLOS ONE papers with results tables, downloaded as publisher XML and converted into a standard format (`corpus/`, format in `docs/SCHEMA.md`).
2. **Built an answer key.** A separate agent wrote a program that copies a paper and inserts one known error: for example, significance stars that don't match the numbers, or percentages that no longer add up. It covers ten kinds of errors. Six were used to guide the search; four were **held back** and only used at the very end. The program was frozen (hashes in `planted/FROZEN.txt`) before any search began (`planted/`).
3. **Ran existing checks.** GRIM and statcheck (`baselines/`).
4. **Let the agent invent checks.** Three independent runs of Claude each wrote a checking program, got a score from a black-box judge (`climb/score.py`), read its mistakes, and rewrote the program, for about 15 rounds each. Every attempt, score and keep/reject decision is logged (`climb/runs/*/trajectory.jsonl`).
5. **Final exam.** The resulting checks were scored once, on 16 papers and 4 error kinds the agent never saw (`results/run_final.py`; the selection rule was committed before the search finished).
6. **Audit of the baselines.** Our first statcheck/GRIM counts came from our own text extraction, so Claude independently re-read all 40 papers and the same deterministic code re-ran on what it found (`audit/`).
7. **Deposited data (stretch).** For papers whose authors posted data, a reconciliation check compares reported statistics with values recomputed from the deposit (`deposits/`).
8. **Second opinion on real flags.** Each flag on the original, unmodified papers was judged by a separate model that was not told which check raised it (`verify/verify.py`).

The git history shows the order: search setup → frozen answer key → baselines → pre-registered scoring rule → search results → held-out scores.

## Results

See [`RESULTS.md`](RESULTS.md). In short:

| | GRIM | statcheck | Agent-discovered tests (3 together) |
|---|---|---|---|
| Papers with inputs the check can use (of 40, audited) | at most 5 | 9 (14 with AI extraction of tables) | n/a |
| Papers where it caught a planted error (of 37) | 0 | 6 | 32 |
| Held-back error kinds caught, unseen papers | 0 of 78 | 9 of 78 (all of the one kind it was built for) | 23 of 78 |
| Original papers with a model-confirmed inconsistency | 0 | 3 | 9 |

The pilot also shows the main weakness the full project must address: the agent's checks overfit, and recall fell from 58% on development papers to 39% on unseen papers for the best single check.

## What is not here, and why

Flags that name specific papers are kept private. We contact authors and journals before publishing any paper-specific flag. Aggregate results are public.

## How to run

Python 3.11+ and [uv](https://github.com/astral-sh/uv). Model calls need `ANTHROPIC_API_KEY`.

```bash
python3 -I climb/score.py climb/runs/B/best.py --split visible --papers dev     # score a test
uv run --with anthropic --with numpy --with scipy python -I climb/climb.py --name X --family A --iters 15 --budget 8
uv run --with numpy --with scipy python -I results/run_final.py                # one-shot final scoring
```

## Notes

- Everything here was written from scratch for this pilot, with AI coding assistance (Claude).
- statcheck results come from `baselines/statcheck_py.py`, a Python port of the published method (Nuijten et al., 2016). The R package wrapper is in `baselines/statcheck_check.py`.
- `planted/visible/` and `planted/heldout/` contain **deliberately altered** copies of published papers. They are test fixtures, not the papers as published.
- Paper content is from PLOS ONE under CC BY 4.0; each file keeps its DOI. Code is MIT-licensed (`LICENSE`).
