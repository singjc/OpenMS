# AMEND_CALIBRATION (PLAN_FASTA), 28 November 2026

**When it was written.** After the v1 arms at the frozen code (`v1_frozen_code/INTERIM_V1.md`), and before any run of the changed code. Results already seen:

- P: N = 15,942, completed.
- F: failed at calibration with 24 anchors.
- Pe, H: their logs showed intermediate scoring passes before I cancelled them. Pe passed 13,556 targets at q ≤ 0.01 and H 37,654. These are not final results, and I read no entrapment counts.
- Fe: still predicting MS2 when cancelled.

**Who decided.** I adopted this without the user's explicit choice. I offered A (this amendment) and B (close with v1). The user did not answer, so I went ahead with A, the minimal variant A1. The user may still reject it; v1 then stands as the result.

## What changes
PLAN_FASTA step 3 says the panel pairs "are the only calibration pairs". At a445c1892e the code applied that to the survey only. The final search calibrated on every pair of the final library. The calibrator takes 40,000 of these at an even RT stride, and in F most of them are absent peptides.

The fix is 7fb7a6fb8d:
- The final library writes its kept panel pairs first. Its footer key `opendia.fasta.panel_precursors` gives their count.
- After reading a `-fasta` library, OpenDIA keeps only those as calibration pairs.
- **Variant A1:** the kept set is unchanged, so the prefilter alone decides. A panel pair below the top K is not in P's library, and P calibrates on the kept panel pairs only. F keeps every pair, so it has all 38,372.
- The final-library cache key changes (`final=2`).

## What does not change
Criteria, arms, K* = 800,000, the entrapment FASTA, r0, the evaluator (`eval_fasta.py`), seeds and thresholds.

The rerun repeats every arm from a cold cache, on the same codon layout, with a new SIF. v1 is reported next to v2 in RESULT_FASTA:
- P_v1 = 15,942 is kept in the report.
- If v2 changes P, both numbers are given.
- Criterion 1 is rechecked at the new commit, since the empirical path should be untouched.
