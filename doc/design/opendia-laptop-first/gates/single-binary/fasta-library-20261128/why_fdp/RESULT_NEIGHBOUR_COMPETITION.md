# RESULT_NEIGHBOUR_COMPETITION (28 November 2026)

Plan: PLAN_NEIGHBOUR_COMPETITION (frozen aac76d23cc). The evaluator `nc_eval.py` was committed before the results
(939cc3bcf2). The G0 dry run, job 45461949, passes: in all 11 arms the recomputed N equals the tool's N exactly, and
no fragments are missing. Main job 45474615 took 3 min and 14 GB. Outputs: `nc_main.json`, `nc_main_eval.json`,
`nc_eval.txt`.

## Verdict (primary S = 0.5): **rejected in this form (E1 fails)**
| criterion | before | after NC, S = 0.5 | result |
|---|---|---|---|
| E1 excess (mean of 7 prefilter arms − FbF) | +0.37 points | +0.24 points (limit ≤ +0.185) | **fail** |
| E1 FbV1 − FbF, 95 % CI | +0.51 (+0.06, +0.96) | +0.49 (+0.05, +0.93) | **fail** (excludes 0) |
| E2 pooled T_e / D_e (CI) against FbF | 1.77 (1.58-1.99) against 1.29 | 1.49 (1.33-1.67) against 1.21 | fail (supporting only) |
| E3 human IDs after / before, prefilter arms | | 0.996-1.013 | pass |
| E4 custom ΔFDP_c; human ratio (Fe / Pe / Pe') | | −0.15 / −0.23 / −0.15 points; 0.984 / 0.989 / 0.997 | pass |

Per arm FDP_c, before → after (S = 0.5):

| arm | before | after |
|---|---|---|
| FbF | 1.21 % | 1.16 % |
| FbV1 | 1.72 % | 1.64 % |
| FbV1_s1 | 1.60 % | 1.52 % |
| FbV1_s2 | 1.61 % | 1.38 % |
| FbV1_g1 | 1.68 % | 1.60 % |
| FbP | 1.53 % | 1.23 % |
| FbP_g1 | 1.47 % | 1.17 % |
| FbPp | 1.47 % | 1.21 % |
| Fe | 1.63 % | 1.48 % |
| Pe | 1.34 % | 1.11 % |
| Pe' | 0.94 % | 0.79 % |

## Sensitivity (descriptive only; the plan fixed S = 0.5 as primary, so these do not change the verdict)
- **S = 1/3:** excess +0.18, FbV1 − FbF +0.35 (−0.08, +0.79), pooled T_e/D_e 1.36 (1.21-1.52) against 1.22. E1-E4 would all pass.
- **S = 2/3:** excess +0.24, FbV1 − FbF +0.38 (−0.05, +0.82). E1 would fail on the excess.
- **Not monotone in S.** Removing a winner promotes its partner and moves the TDC cut, so nearby levels can differ.
  The counts are small (about 100 entrapment hits per arm). The differences between S levels are within noise.

## Reading
- **Direction is right everywhere:**
  - FDP_c falls in all 11 arms at every S.
  - Human IDs are essentially unchanged.
  - Custom-entrapment FDP falls too. That entrapment has no designed neighbours, but it does have real co-eluting
    human chimeras.
- **The rule works for V0 (FbP, FbP_g1) and V1+pool (FbPp).** Their FDP falls to 1.17-1.23 %, at the FbF level.
- **It barely moves V1-only arms (FbV1 s0/g1).** These keep most of the excess, and they drive the E1 failure. This
  fits RESULT_PREFILTER_FDP: about 40 % of the excess was not parent co-elution, and V1's survey favours exactly the
  pairs with predicted-spectrum agreement. That remainder is not explained by an identified co-eluting neighbour, so
  this rule cannot reach it.
- **Validity:**
  - The rule is label-blind and the same for targets and decoys. It only swaps which pair member competes, so TDC's
    counting argument is unchanged.
  - It is a post hoc filter on final scores, not integrated into OpenDIA.
- Per the plan, no further variants are tried without a new frozen plan.

## Next options (each needs its own plan)
1. **Integrated neighbour competition in OpenDIA scoring** (fragment-claim sharing before the GBM, not after). Worth
   it for V0 and V1+pool, where this post hoc form already closes the gap.
2. **Diagnose the V1-only residual first:**
   - Characterise accepted null target winners in FbV1 that NC does not remove (unidentified neighbour, partial
     fragment overlap, RT or IM).
   - Do this before designing a fix for it.
