# Amendment 1: label mode was not passed into the container (27 November 2026)

Written after attempt 1 and before attempt 2. The evaluators, inputs, records and reported quantities are unchanged.

## What happened

- **Attempt 1:** job `44884353`, COMPLETED 0:0 in 11:46. It ran all 26 evaluations, and every one exited 0.
- **Problem:** `singularity exec --cleanenv` drops the environment of the calling shell. The per-command `LABELS=... RELABEL=...` therefore never reached Python, so every evaluation ran with the default `LABELS=old`. Each output says so in its first line, `labels: old; entrapment precursors relabelled as human: 0`, including the files named `*_corrected.txt`.
- **How it was caught:** by reading the outputs. No corrected result was reported.

## What attempt 1 still shows

It is a valid run of `LABELS=old`, and so it is the reproduction check of the plan.
- The per-seed rows and pooled lines are identical to the published EVAL tables: `EVAL_ODP_D/H`, `EVAL_PLAN_SIB_D/H`, `EVAL_PLAN_GBT_D/H`, `EVAL_PLAN_V7M6_D` (pooled) and `EVAL_DIANN_R5` (D and H).
- RTR D matches the RESULT_RTR_D numbers (14,895 / 15,369).
- The outputs are kept unchanged under [attempt1/](attempt1/).

## Change for attempt 2

`relabel.sbatch` now:
1. sets `LABELS` as a shell variable and passes it, with `RELABEL`, by `singularity --env`;
2. stops with exit 4 if an output's first line does not state the requested mode;
3. keeps any earlier `out/` as `out_attempt1/` and writes a fresh `out/`.

A local check in the evaluator container with `--cleanenv` and `--env` gives `labels: corrected` and relabels 10 precursors on the synthetic library, against 0 for `old`.

`bundle.sha256` is regenerated for the changed script. Only the script changed; the evaluators, list and container are the same files.
