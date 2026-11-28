# RESULT_ENTRAP_FDRBENCH: custom entrapment against FDRBench (28 November 2026)

Plan: PLAN_ENTRAP_FDRBENCH (10a8e7586c), frozen before any search. Jobs 45301910-14 (D, seed 0, SIF v4, code
143d3ac7c6). Database: FDRBench v1.1.1 protein-level shuffled entrapment of the same human FASTA (4d6272d4...f20f).
Evaluator: `eval_fb.py` (label `_p_target`). FDRBench's own combined estimator agrees to within 1e-6 (ESTIMATOR_CHECK
and `eval/*_fdp.csv`).

## Results (N = target precursors at run q ≤ 0.01; 95 % Poisson CI on N_e)
| arm | N | N_e | r0 | FDP_c(r0) | r_lib | FDP_c(r_lib) | N@FDP_c 1 % | wall | peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| FbF (K = 0) | 8,074 | 49 | 1.000 | 1.21 % (0.87-1.55) | 0.999 | 1.21 % | 7,212 | 10,449 s | 71.5 GB |
| FbP (V0, K*) | 11,120 | 85 | 1.000 | 1.53 % (1.20-1.85) | 0.970 | 1.55 % | 10,440 | 1,641 s | 35.0 GB |
| FbV1 (V1, K*) | 14,994 | 129 | 0.999 | 1.72 % (1.42-2.02) | 0.959 | 1.76 % | 14,210 | 2,714 s | 34.0 GB |
| FbPp (V1+pool, K*) | 20,258 | 149 | 0.999 | 1.47 % (1.24-1.71) | 0.804 | 1.65 % | 19,502 | 2,419 s | 33.9 GB |

FDRBench combined_fdp at q ≤ 0.01: FbF 0.01214, FbP 0.01529, FbV1 0.01721, FbPp 0.01471. Lower bounds: 0.0061,
0.0076, 0.0086, 0.0074.

Custom-entrapment counterparts (RESULT_FASTA, RESULT_V4): Fe 1.63 % (1.20-2.06); Pe 1.34 %; Pe' 0.94 % (r0) /
2.54 % (r_lib).

## Comparisons as frozen
| | custom | FDRBench | outcome |
|---|---|---|---|
| **C1** base FDP_c, K = 0 | Fe 1.63 % (1.20-2.06) | FbF 1.21 % (0.87-1.55) | **agree** (CIs overlap) |
| **C2** V0 prefilter minus K = 0 | -0.29 (-0.82, +0.24) | +0.32 (-0.16, +0.79) | signs differ; both CIs contain 0. **Not the same sign** (fails as worded). Neither detects an effect. |
| **C2** V1+pool minus K = 0, r0 | -0.69 (-1.16, -0.22) | +0.26 (-0.16, +0.67) | **disagree** (opposite signs, CIs do not overlap) |
| **C2** V1+pool minus K = 0, r_lib | +0.93 (+0.24, +1.62) | +0.44 (+0.01, +0.87) | **agree** (same sign, CIs overlap) |
| **C3** r_lib / r0 under V1+pool | 0.27 | 0.80 | prediction met (between 0.27 and 1, below 1) |
| **C4** N ratio to K = 0 | P/F 0.63, P'/F 0.99 (pure human) | FbP 1.38, FbV1 1.86, FbPp 2.51 | not comparable: FbF loses power in a 2× search space |

## What this says
1. **The custom entrapment gives the same base picture.** At K = 0, both put FDP_c at or slightly above 1 %
   (1.21 % and 1.63 %, overlapping). FbF's CI contains 1 %, Fe's does not. So the 1.6× seen with the custom database
   is partly a property of that database or noise. FDRBench puts the base search closer to nominal.
2. **The V1+pool problem reproduces in its conservative form.** FDRBench, whose entrapment proteins are full-length
   and can collect protein evidence, shows V1+pool raising FDP above K = 0. That agrees with the custom r_lib reading
   and contradicts the custom r0 reading (-0.69). The optimistic custom r0 estimate was the misleading one, as
   AMEND_SURVEY_POOL predicted. The asymmetry is much milder with FDRBench (0.80 against 0.27), so its r0 and r_lib
   estimates differ by only 0.18 points.
3. **New: the prefilter itself raises FDP, and V1 alone the most.** FbV1 - FbF = +0.51 (+0.06, +0.96), the only
   single-variant difference with a CI excluding 0. Every prefilter arm is at 1.5-1.8 %, CIs excluding 1 %. Pooling
   does not add to this (FbPp r0 1.47 % ≤ FbV1 1.72 %). The custom database had no V1-only arm, so this was not seen
   before. It contradicts option 1 of RESULT_V4 ("V1 without pooling keeps r0 valid"): r0 does stay valid, but the
   FDP of V1 alone is not lower.
4. Applying PLAN_FASTA's criterion-2 tolerance (reference + 0.003) with FbF as the reference (1.51 %), for
   information only, since this plan does not decide defaults: FbP 1.53 % and FbV1 1.72 % exceed it. FbPp passes with
   r0 (1.47 %) and exceeds it with r_lib (1.65 %).

## Caveats
- FDRBench uses r = 1 and converts I to L in the human half. So N is not comparable with the pure-human or custom
  arms, only ratios within the FDRBench arms.
- One run (D), one seed. Differences of about 0.3-0.5 points are at the edge of what 49-149 entrapment hits resolve.
- A hypothesis, not tested: the prefilter selects pairs on max(T, D) of the survey. Target-decoy competition stays
  symmetric, but the search is now restricted to pairs with survey signal, where false targets may resemble real
  signal more than OpenDIA's decoys do.
