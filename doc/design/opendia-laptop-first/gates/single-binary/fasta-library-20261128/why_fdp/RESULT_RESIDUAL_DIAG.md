# RESULT_RESIDUAL_DIAG (28 November 2026)

Plan: PLAN_RESIDUAL_DIAG (frozen b86d8de81d). Dry-run gate, job 45490520: passes in 8 of 8 arms (100 % parent,
survey and features joins, 0 missing fragments). Main job 45494485 → `rd_main.json`.

**Disclosed change after the main run.** Rerun 45499697 → `rd_main3.json`.
- What was missing: the plan's D rule needs Q3 and Q2 for **D-only** winners, but the first script reported only
  pooled C+D. I added C-only and D-only groups, the SMD of every feature, and a per-pair dump.
- What did not change: the category counts are identical between the two runs (checked). No rule or threshold
  changed.

## Where the excess null targets sit (excess = T − D among accepted null pair winners)
| group | T / D | A | B_lib | B_cut | C | D |
|---|---|---|---|---|---|---|
| **V1** (FbV1 s0/s1/s2/g1, pooled) | 490 / 250 | 77 / 0 → **32 %** | 0 / 0 | 0 / 0 | 195 / 143 → **22 %** | 218 / 107 → **46 %** |
| V0 (FbP, FbP_g1) | 168 / 102 | 33 / 0 → 50 % | 0 | 0 | 70 / 52 → 27 % | 65 / 50 → 23 % |
| V1+pool (FbPp) | 149 / 98 | 21 / 2 → 37 % | 0 | 0 | 68 / 50 → 35 % | 60 / 46 → 27 % |
| K = 0 (FbF) | 49 / 37 | 2 / 0 → 17 % | 0 | 0 | 13 / 13 → 0 % | 34 / 24 → 83 % |

**Categories:**
- A: the parent is accepted and co-elutes.
- B_lib, B_cut: the parent is present at K = 0 but lost in this arm, or was cut by the prefilter.
- C: another accepted target co-elutes and shares a fragment.
- D: none of the above.

The per-arm split is stable across seeds and GBM rounds (`rd_main3.json`).

## Frozen decision rules applied (V1 group)
Both branches reach 40 %: A + C = 54 % and D = 46 %. The plan did not say what happens when two branches qualify. I
report both and do not choose between them.

**A + C ≥ 40 %:**
- Q2 asks whether the interference resolver sees the neighbour. A+B has no decoy winners at all, so the comparison
  uses C. The SMD of i_explained (target − decoy) is −0.19 (s1), −0.45 (s2) and +0.06 (g1). All are below 0.5.
- **The rule's verdict: the resolver does not see the neighbour. Fix: integrated claim competition.**
- Supporting evidence: among accepted null winners, i_nclaim > 0 for 100 % of targets **and** 100 % of decoys, and
  mean i_explained is 0.92-0.96 for both. The resolver links nearly every window, so in its current form the feature
  is saturated and carries no target/decoy contrast.

**D ≥ 40 %:**
- Q3 for D target winners is 0.90-0.93 (≥ 0.60), and for D decoy winners 0.12-0.26 (≤ 0.50). The rule's verdict is
  "the V1 survey favours shuffled targets through a non-neighbour route; next, a prediction-quality diagnostic".
- **Caveat: this condition turned out not to be diagnostic.** The same pattern holds in category C (0.89-0.92
  against 0.13-0.16) and in V0 (0.86-0.90 against 0.30-0.35), and V0 uses a different survey. Among accepted pairs,
  the survey simply agrees with whichever member won the final score. That is selection on the outcome.
  RESULT_PREFILTER_FDP had already found the survey symmetric among kept pairs (0.502-0.506). So this branch's
  survey attribution is weak.

**B_lib = B_cut = 0 in every arm.** The prefilter never cuts, or loses, the parent of an accepted null target. A
neighbour-aware survey selection would not help.

**Predictions:**

| prediction | outcome |
|---|---|
| V1: A + C largest share but < 60 %, D ≥ 25 % | met (54 %, 46 %) |
| V0 and V1+pool: A + C ≥ 60 % | met (77 %, 73 %) |
| K = 0: no dominant stage | **not met** (D is 83 % of a small excess of 12) |

## Q4: why post hoc NC missed category A (V1)
- In category A, the parent outscores the null winner in 73 of 77 pairs (95 %). The shared count is ≥ 3 of 6 in
  71 of 77 (mostly 4 of 6). So the S = 0.5 rule did fire on most of them.
- The rule also fired on null **decoy** winners at the same rate as on null target winners:
  - FbV1: 43 of 129 null targets and 19 of 62 null decoys (about 1/3 each);
  - FbP: 33 of 85 and 18 of 51.
- Each removal promotes the pair partner and loosens the TDC cut. The net change in N_e was therefore small (FbV1
  129 → 124).
- **Conclusion:** a rule that asks "does any co-eluting accepted target in the 25 Th window share k of 6 fragment
  m/z" is not specific. Random fragment coincidences among many co-eluting anchors hit decoys as often as targets.
  This is the same saturation the resolver features show.

## Post hoc (not in the plan; descriptive only)
The plan's shares are within-arm, so they include the excess already present at K = 0. Here I instead split the
**increase** over K = 0, per identification: (T_c − D_c) / N, arm mean, minus FbF.

| | A | C | D | total |
|---|---|---|---|---|
| V1 increase over K = 0 | +0.105 (41 %) | +0.088 (34 %) | +0.063 (25 %) | +0.256 points |
| V0 | +0.122 | +0.080 | −0.057 | +0.146 points |
| V1+pool | +0.069 | +0.089 | −0.055 | +0.103 points |

About 75 % of what V1 adds is neighbour-related (A + C). D is mostly a baseline that K = 0 also has.

## What this means for a built-in fix (needs its own frozen plan)
1. **The neighbour route is the main prefilter-driven component.** It is not a survey-selection problem
   (B_cut = 0).
2. **The existing InterferenceResolver is the right place, but it is not specific enough.**
   - Its features are saturated: every window has claimants, for targets and decoys alike.
   - Fragment-m/z overlap with "any co-eluting accepted precursor" is equally non-specific.
3. **What a specific claim competition needs, before the GBM:**
   - **Shared signal:** a fragment's signal counts as shared only when the claimant's *intensity-predicted*
     contribution on that fragment (with apex and peak-shape agreement) accounts for it.
   - **Residual scoring:** the candidate is scored on its residual. This matches the resolver's design intent.
   - **Validation before wiring:** on the existing features, the residual signal must separate category A and C
     targets from decoys without firing equally on decoys. This is the property the post hoc rule lacked.
4. **D (25 % of the V1 increase) is not explained** at this resolution. The Q3 survey contrast is not diagnostic.
   The plan's next step for D, a check of prediction quality for shuffled targets against decoys, stays open.

Caveats: one run; small counts (about 50-130 null winners per arm); differences below 15 points of share are not
distinguished.
