# PLAN_NEIGHBOUR_COMPETITION: does a label-blind neighbour-competition rule remove the prefilter FDP excess? (28 November 2026)

Written before any removal result is computed. The only run that has happened so far is the label-blind dry run
(job 45461949, `NC_DRY=1`). It reports plumbing only (tool N against recomputed N, candidate, anchor and
missing-fragment counts). It does not remove anything and never reads entrapment labels.

## Why
RESULT_PREFILTER_FDP (764cdff496) found that accepted null target winners co-elute with an identified neighbour that
shares their fragments (the FDRBench parent) far more often than null decoy winners do (15.6 % against 0.2 %). The
prefilter raises the share of such pairs among accepted nulls from 2.3 % to 10.1 %. Proposed remedy: neighbour
competition. A candidate whose top fragments are already explained by a co-eluting, better-scoring accepted
precursor loses.

## Rule (`nc_remedy.py`, post hoc on existing outputs; no new search)
Label-blind and identical for target and decoy rows:
- **Pair winners:** per pair (precursor_id // 2), the higher discriminant_score wins, with ties to the decoy. Plain
  TDC is used: q = cummin((D + 1) / T).
- **Anchors:** target pair winners with q ≤ 0.01. **Candidates:** pair winners with q ≤ 0.05.
- **Explained:** candidate c is explained by anchor a if all of these hold:
  - score(a) > score(c);
  - a is not a charge sibling of c's target sequence;
  - rt_start(a) ≤ rt_apex(c) ≤ rt_end(a);
  - |Δ precursor m/z| ≤ 25 Th (the window width);
  - |Δ library 1/K0| ≤ 0.05;
  - share ≥ S, where share is the fraction of c's top-6 search-library fragments within 10 ppm of any of a's
    library fragments.
- **Removal:** an explained candidate is removed. Its pair partner then becomes the pair winner with its own score.
  TDC is recomputed.
- **Thresholds:** **S = 0.5 is primary** (3 of 6 fragments). S = 1/3 and 2/3 are reported as sensitivity only. All
  constants are fixed here (MZ 25, IM 0.05, 10 ppm, top 6), and none are tuned afterwards.
- **Synthetic test:** `nc_remedy_test.py` plants one entrapment target sharing 4/6 fragments with a co-eluting human
  anchor. The test passes locally and on codon.

## Arms
- **FDRBench:** FbF (K = 0 reference), and the prefilter arms FbV1 s0/s1/s2, FbV1_g1, FbP, FbP_g1 and FbPp.
- **Custom entrapment:** Fe, Pe and Pe'.
- FDP_c = N_e (1 + 1/r0) / N with each arm's r0 (as eval_fb / eval_v4). The CI is Poisson on N_e. A difference
  between arms combines their two variances (independent runs).

## Plumbing gate (dry run; must pass before the criteria are read)
- **G0:** for every arm, recomputed N equals the tool's N within 0.2 %, and missing_fragments = 0.
  - On FbV1, a local check of the same TDC already gives an exact match (14,994 = 14,994, max |Δq| 7.6e-5).
- If G0 fails, I fix the plumbing only (never the rule) and redo the dry run before any removal result.

## Criteria (primary S = 0.5)
- **E1 (excess removed):**
  - The 7 FDRBench prefilter arms are pooled as the mean of their FDP_c.
  - Before NC, excess_before = mean FDP_c(prefilter) − FDP_c(FbF). After NC, excess_after is the same quantity,
    with both sides after NC.
  - **E1 passes if excess_after ≤ 0.5 × excess_before, and FbV1 (s0) after minus FbF after has a 95 % CI that
    includes 0.**
- **E2 (null balance):** pooled T_e / D_e over the 7 prefilter arms after NC has a 95 % CI (binomial on
  T_e / (T_e + D_e)) that includes FbF's T_e / D_e after NC.
- **E3 (cost):** in every prefilter arm, human target IDs after NC are ≥ 0.95 × before. A rule that wins only by
  discarding real IDs is not a fix.
- **E4 (no harm on the custom entrapment):** FDP_c after NC for Fe, Pe and Pe' is each ≤ its value before + 0.25
  points. Human IDs are ≥ 0.95 × before.
- **Verdict:**
  - "NC works" requires E1, E3 and E4.
  - E2 is supporting evidence.
  - If E1 holds but E3 fails, I report "removes the excess at an unacceptable ID cost".
  - If E1 fails, the remedy is rejected for this form, and I make no further variants without a new plan.

## Validity notes (to state in the result whatever the outcome)
- The rule never reads labels and is the same for targets and decoys. It changes only which pair member competes,
  so TDC's counting argument is unchanged. It removes exactly the non-exchangeable component (a null target with a
  borrowed signal).
- It is a post hoc filter on final scores, not integrated scoring. A real implementation would act inside OpenDIA
  and needs its own plan and runs.
- It also removes genuinely co-eluting chimeric true IDs. E3 measures that cost.
- One run (D), one data set. The S sensitivity levels are descriptive.
