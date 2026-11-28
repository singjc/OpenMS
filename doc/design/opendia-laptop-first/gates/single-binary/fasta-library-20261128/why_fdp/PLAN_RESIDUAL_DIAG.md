# PLAN_RESIDUAL_DIAG: where in the pipeline do null targets get their edge over decoys? (28 November 2026)

Written before any target/decoy split of the categories below has been seen. The only run so far is the label-blind
dry run, job 45490520 (`DRY=1`). It reports joins and mapping rates only.

## Why
- The user wants a fix built into OpenDIA from the start, not a post hoc filter.
- RESULT_NEIGHBOUR_COMPETITION showed that post hoc NC closes the gap for V0 and V1+pool but not for V1 alone. About
  40 % of the excess is not explained by an identified, co-eluting parent.
- Before designing an integrated fix, I need to know **at which stage** each excess null target gains its edge:
  - survey or prefilter;
  - interference resolution (the existing `InterferenceResolver`, features i_*);
  - identification competition.

## Categories (`diag_residual.py`)
Each accepted null pair winner (q ≤ 0.01, FDRBench entrapment pair) gets exactly one category. The same rule applies
to target and decoy winners, using the winner's own apex and fragments.

| category | definition | stage it points to |
|---|---|---|
| **A** | FDRBench parent accepted in this arm, span contains the winner apex | Interference resolver failed to explain the borrowed signal, or competition never ran |
| **B_lib** | Not A. Parent accepted in K = 0 (FbF) with a co-eluting span, scored in this arm's library but not accepted | Parent present but lost here; the entrapment target won the shared signal |
| **B_cut** | As B_lib, but the parent is not in this arm's searched library | Prefilter cut the parent, so its signal is unclaimed |
| **C** | Not A/B. Another accepted target co-elutes (window, IM), sharing ≥ 1 of 6 fragments | Non-parent neighbour; same stage as A |
| **D** | None of the above | Not neighbour-related: survey or scoring behaviour |

**Per arm and pooled** (groups: V1 = FbV1 s0/s1/s2/g1, V0 = FbP/FbP_g1, V1+pool = FbPp, K = 0 = FbF): T_c, D_c,
excess_c = T_c − D_c, and share_c = excess_c / (T − D).

**Secondary questions:**
- **Q2 (arms with `-features_out`: s1, s2, g1, FbP_g1):** for the winner window, compare target winners with decoy
  winners (standardised mean difference of all features) in groups A+B and C+D. Do i_explained and i_nclaim see the
  neighbour?
- **Q3 (survey arms):** among accepted null winners, the share whose survey target co-elution exceeded the decoy's,
  by group.
- **Q4 (category A targets):** whether the parent outscores the winner, and how many of 6 fragments are shared. This
  explains why post hoc NC at S = 0.5 missed them.

## Decision rules (fixed now; each points to one integrated fix, built and tested under its own plan)
For the **V1 group**, the stage carrying ≥ 40 % of the excess is named the primary target:
- **A + C ≥ 40 %:**
  - If the target winners' i_explained mean exceeds the decoy winners' by < 0.5 SD (Q2), the resolver does not see
    the neighbour. Fix: integrated claim competition. Shared-fragment area is assigned to the better-supported
    precursor before scoring, and the GBM sees the residual-only features.
  - Otherwise the resolver sees it but the GBM does not use it. Fix: make the shared-fragment evidence a hard input,
    or give it a monotone constraint.
- **B_cut ≥ 40 %:** the prefilter removes the parent and leaves its signal unclaimed. Fix: neighbour-aware survey
  selection. Selecting a pair also selects isobaric or shared-fragment neighbours in the same window, for targets and
  decoys symmetrically. The alternative is to compute claims against the full candidate set.
- **B_lib ≥ 40 %:** the parent is scored but loses the shared signal to the entrapment target. Same fix as A, with
  the competition done at window level before identification.
- **D ≥ 40 %:**
  - **Survey cause:** Q3 for D target winners is ≥ 0.60 while D decoy winners are ≤ 0.50. The V1 survey's predicted
    spectrum favours shuffled targets through a non-neighbour route. Next step: a diagnostic of decoy versus shuffled
    target prediction quality, not a fix yet.
  - **Otherwise:** the cause is not located at this resolution. I report that and stop.
- **No stage reaches 40 %:** I report the split and propose no single fix.

**Predictions** (written now, for honesty):
- **V1:** A + C is the largest share but < 60 %, and D is ≥ 25 %.
- **V0 and V1+pool:** A + C ≥ 60 % (where post hoc NC worked).
- **K = 0:** small excess with no dominant stage.

## Plumbing gate (dry run)
- Every arm needs a parent mapping for ≥ 99 % of accepted null pairs, and missing_fragments = 0.
- The features join must be ≥ 99 % where features exist, and the survey join ≥ 99 % where a survey exists.
- If any check fails, I fix the plumbing only, then rerun the dry run.

## Caveats (state in the result)
- These are descriptive diagnostics on existing outputs, with one run and small counts (about 50-130 nulls per
  arm). Categories with a share below 15 points are not distinguished from each other.
- Category A uses this arm's identifications, B uses FbF's, and C uses any accepted target. A neighbour that no arm
  identified falls into D.
- Nothing here changes OpenDIA. Any fix it points to gets a frozen plan, an implementation inside the pipeline, and
  fresh runs.
