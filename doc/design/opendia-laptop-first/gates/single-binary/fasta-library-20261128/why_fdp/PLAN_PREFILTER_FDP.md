# PLAN_PREFILTER_FDP: why do prefiltered searches show higher entrapment FDP? (28 November 2026)

Written after the diagnostics below (from existing outputs only), and **before** any of the causal arms run. The
causal arms, hypotheses, predictions and decision rules are fixed here.

## Facts from existing outputs (`d1_exchange.py`, `d2_survey_symmetry.py`, `d3_null_flip.py`, `s1_selection.py`)
1. **TDC itself is honest.** In every arm FDR_tdc = (D + 1) / T = 1.00 % at the cut. No decoy loser carries q ≤ 0.01.
2. **The surplus is confined to the very top of the null.** Over all entrapment (null) pairs, the target beats its
   own decoy in 0.500-0.501 of pairs in every arm. That holds for all pairs and for the 10-100 % ranks. Only the
   top 0.1 % (where the q ≤ 0.01 cut falls) departs:

   | arm | top 0.1 % null target-win share | GBM rounds run |
   |---|---|---|
   | FbF (K = 0) | 0.509 ± 0.020 | 1 (r2 stopped: zero decoy spread) |
   | Fe (K = 0, custom) | 0.528 ± 0.026 | 1 (r2 stopped) |
   | FbP | 0.590 ± 0.052 | 3 |
   | FbV1 | 0.615 ± 0.050 | 3 |
   | FbPp | 0.586 ± 0.058 | 3 |
   | Pe | 0.563 ± 0.060 | 3 |

   Accepted null pairs (q ≤ 0.01): T_e/D_e = 1.29 (FbF), 1.67 (FbP), 2.08 (FbV1), 1.49 (FbPp). Under exchangeable
   decoys this ratio should be 1.
3. **The survey selection is symmetric.** Among kept entrapment pairs, the survey favours the target in 0.502-0.506
   (FbP, FbV1, FbPp), against 0.519-0.527 for kept human pairs. Selecting on max(T, D) does not tilt null pairs
   toward targets.
4. **The tilt is not a property of the selected pairs.** Scored by the K = 0 arm (FbF), entrapment pairs that the
   prefilter keeps are not over-won by targets beyond the K = 0 baseline (FbPp kept: T_e/D_e 1.33 against 1.31 for
   not kept). The FbP/FbV1 kept subsets have too few FbF-accepted null pairs (12 and 7) to say.
5. **Confound: the K = 0 arms never ran GBM rounds 2-3.** In FbF and Fe, GBM round 2 throws "zero decoy spread"
   and the run keeps round 1. Every prefilter arm runs all 3 rounds. The library size is what differs: K = 0 has
   about 40M windows, so a fold's hard-negative set is mostly tied. So "prefilter against K = 0" also compares
   "3 GBM rounds against 1".
6. **Precedent:** RESULT_GBMNULL (Astral) found the GBM can make q optimistic only near the threshold, invisible in
   the null bulk. That is the same signature as fact 2.

## Hypotheses
- **H1 (GBM rounds):** the extra hard-negative GBM rounds learn features that rank top-null targets above their
  decoys (targets and decoys are not exchangeable at the top under the trained score). The prefilter raises FDP
  only because it lets those rounds run.
- **H2 (selection):** restricting the search to survey-positive pairs makes top null targets more target-like than
  their decoys. This is weakened by facts 3 and 4 but not excluded.
- **H3 (small-library scoring, any round):** in a 0.8M-pair library the first-pass scorer and round-1 GBM are
  already tilted, independent of extra rounds.

## Causal arms (D, seed 0, SIF v4, FDRBench database, cold caches, `arms_pf.sbatch`)
| arm | options | isolates |
|---|---|---|
| FbV1_g1 | `-survey_predicted -gbm_rounds 1` | H1: FbV1 with the K = 0 arm's effective GBM depth |
| FbP_g1 | default + `-gbm_rounds 1` | H1 replicate on V0 |
| FbV1_s1 | `-survey_predicted -seed 1` | seed noise of FbV1 at 3 rounds |
| FbV1_s2 | `-survey_predicted -seed 2` | seed noise of FbV1 at 3 rounds |

Each arm also writes `-features_out` (network_score = the pre-GBM score), so H3 can be read without another run:
the FDP and top-null share under the network score of the same windows.

## Predictions and decision rules (FDP_c(r0) via `eval_fb.py`, top-null share via `d3_null_flip.py`)
- **H1 supported** if FbV1_g1 and FbP_g1 both have FDP_c within 0.003 of FbF (1.21 %), each at or below
  FbF + 0.003 = 1.51 %. Their N must also drop by no more than 15 % from FbV1 / FbP (so the GBM is the cause, not
  a loss of everything).
- **H1 rejected** if FbV1_g1 FDP_c ≥ 1.60 % (still above FbF's upper CI 1.55 %).
- **H3 read-out:** in each causal arm's features table, FDP_c at network q ≤ 0.01 (pair winners on network_score,
  the pre-GBM score) against FDP_c at the final q ≤ 0.01. If the network score is already ≥ 1.60 % in FbV1_s1/s2,
  the tilt predates the GBM (H3), whatever the -gbm_rounds 1 arms show.
- **Seed noise:** if FbV1_s1 and FbV1_s2 span below 1.51 %, the FbV1 elevation is not robust to the seed, and H1
  is judged on the mean.
- No result here changes a default. A fix (for example a GBM guard) needs its own frozen plan.
