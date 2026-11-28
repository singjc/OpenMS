# Interim: arms at the frozen code (a445c1892e, SIF 2d1178da…065a), 2026-11-28

Written after these runs, and before any decision on a fix. This is not RESULT_FASTA.

| arm | status | N (q ≤ 0.01) | wall | peak RSS |
|---|---|---|---|---|
| P (D, seed 0, K* = 800k) | completed | **15,942** | 1315 s | 16.4 GB |
| F (D, K = 0) | **failed** after its 3.03M-pair library (540 s MS2): `Calibration found only 24 whole-run stage-1 anchors at q <= 0.05 of 40000 panel targets` | - | - | - |
| Pe, Fe, H | cancelled by me mid-run (see below) | - | - | - |
| P1, P2 | cancelled before start | - | - | - |

- P's survey calibrated on the 38,372 panel pairs, with 608 anchors. P's final search calibrated with 453 anchors and rt_tolerance 39.3 s. Before the prefilter, 721,568 of the 800,000 kept pairs had ≥ 3 fragments.
- Goal N(P) ≥ 24,669: **not met** (15,942, which is 0.58 × 27,410). This matches design recall 0.581 at K*.
- **Cause of the F failure (diagnosis, not yet fixed).** PLAN_FASTA step 3 says the panel pairs are the only calibration pairs. The code applies this to the survey only (`erase_if` in FastaWorkflow). The final search is calibrated on every pair of the final library, and the calibrator picks 40,000 at an even RT stride. In F, these are mostly absent peptides, so calibration collapses (24 anchors). In P, the prefiltered library still gave 453 anchors.
- Criteria 2-4 cannot be evaluated from v1. Criterion 1 passed (see `identity/`).
- I cancelled Pe/Fe/H on noticing the deviation. I read no entrapment counts or final results from them. Their logs did show intermediate scoring passes: Pe 13,556 and H 37,654 targets at q ≤ 0.01. The logs are kept here.
- Full outputs are on codon in `fasta-20261128/v1_frozen_code/`.
