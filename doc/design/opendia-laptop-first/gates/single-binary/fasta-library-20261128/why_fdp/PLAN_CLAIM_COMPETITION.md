# PLAN_CLAIM_COMPETITION: built-in, label-blind claim competition before scoring (28 November 2026)

Written before any run with the flag. Code: OpenMS 86f3199513, `-claim_competition` (off by default). All OpenDIA
tests pass (29/29), including a new resolver unit test.

## Why (RESULT_RESIDUAL_DIAG, RESULT_NEIGHBOUR_COMPETITION)
- **Where the excess comes from:** about 75 % of the null-target excess that the V1 prefilter adds is co-eluting
  neighbour signal. Groups A and C are false targets sharing fragments with an accepted target.
- **Why the existing resolver misses it:** InterferenceResolver lets *every* co-eluting candidate window claim. There
  are about 30 co-eluting anchors per null winner, plus all non-accepted windows. So i_explained and i_nclaim saturate
  (100 % of target and decoy null winners have a claimant), and the GBM learns nothing from them.
- **Why the post hoc rule failed:** it removed decoys as often as targets.
- **The user's requirement:** do it correctly at the start, not as a post hoc filter.

## Mechanism (code 86f3199513)
- **Passes:** OpenDIA already runs two search-extract-score passes (the RT refinement). With `-claim_competition`,
  the first pass's results define the claimants for the second pass.
- **Claimants:**
  - every first-pass pooled paired TDC winner at q ≤ `claim_q` (0.01), **targets and decoys by the same rule**;
  - claimant rank = first-pass score;
  - only the window within `anchor_apex_s` (5 s) of its first-pass apex claims.
- **Claim rule** (`InterferenceResolver::ClaimFilter`): b claims from window a only when
  - rank(b) > own rank(a) (a non-claimant has rank −∞), and
  - b is of a different peptide family (so the pair partner, a charge sibling or a modified form never claims).
  - All other link tests are unchanged: same isolation window, RT overlap, shared fragment within the calibrated
    ppm, and 1/K0.
  - The explained signal uses the existing intensity model: the claimant's scale from its unshared fragments, times
    its library intensity on the shared fragment, times the RT-profile overlap.
- **Four new second-pass features**, the same quantities as i_* but with the filter:
  - cc_explained;
  - cc_libcos_resid (residual cosine to the library);
  - cc_dlibcos;
  - cc_nclaim.
  - They are appended to the feature table, so the network and the GBM see them. Nothing is removed or filtered after
    scoring.
- **Label-blind:**
  - The filter reads no labels.
  - A decoy pair winner at q ≤ 0.01 is a claimant just like a target.
  - A decoy and its target have the same family, so they never claim each other.
  - The symmetry comes from the rule, not from per-label tuning.

## Arms (D, FDRBench database, SIF built from 86f3199513, seed 0; each a fresh full run with cold caches)
| arm | flags | compare to |
|---|---|---|
| CcV1 | `-survey_predicted -claim_competition` | FbV1 (1.72 %; seeds 1.60 / 1.61 %) |
| CcP | `-claim_competition` | FbP (1.53 %) |
| CcPp | `-survey_predicted -survey_pool -claim_competition` | FbPp (1.47 %) |
| CcF | `-prefilter_pairs 0 -claim_competition` | FbF (1.21 %) |
| Ce | custom entrapment, `-prefilter_pairs 0 -claim_competition` | Fe (1.63 %) |
| CcPe' | custom entrapment, `-survey_predicted -survey_pool -claim_competition` | Pe' (0.94 %) |

- **CcV1_s1:** CcV1 with seed 1, for seed spread.
- **Libraries:** the FDRBench and custom FASTAs are already on codon, and the library caches are rebuilt per arm (the
  survey depends only on the library and data, so it is unchanged).
- **Estimators:** FDP_c(r0) by `eval_fb.py` (FDRBench arms) and `eval_v4.py` / `eval_fasta.py` conventions (custom
  arms), each arm with its own r0. Ratios come from `ratio_fb.sh`-style counts on each arm's cache.

## Gates before the criteria are read
- **G1, default path unchanged:**
  - Without the flag, the empirical-library identity run on D (`identity_D_cc.sh`) gives the 3 reference files
    byte-identical and the 15 v1 files the same.
  - Any difference means I stop and fix it.
- **G2, the mechanism fires as designed:**
  - In each Cc run's log, claimant counts are reported for targets and decoys.
  - cc_nclaim > 0 must be **less** saturated than i_nclaim: among all windows, the share with cc_nclaim > 0 is
    ≤ 0.5 × the share with i_nclaim > 0.
  - Checked from `-features_out` on CcV1.
- **G3:** every run completes (no OOM, no exception), and N at q ≤ 0.01 is reported.

## Criteria
- **K1, excess:**
  - mean FDP_c(r0) over the prefilter arms (CcV1, CcV1_s1, CcP, CcPp), minus FDP_c(CcF), is ≤ 0.5 × the same
    quantity before (FbV1 s0/s1, FbP, FbPp minus FbF);
  - **and** CcV1 − CcF has a 95 % CI (independent Poisson on N_e) that includes 0.
- **K2, no FDP harm:**
  - CcF ≤ FbF + 0.25 points;
  - Ce ≤ Fe + 0.25 points;
  - CcPe' ≤ Pe' + 0.25 points (r0 estimate, as reported for Pe').
- **K3, sensitivity cost:** human target IDs at run q ≤ 0.01 are ≥ 0.97 × the matching Fb/F arm in every arm.
- **K4, mechanism (supporting):**
  - With RESULT_RESIDUAL_DIAG's categories rerun on CcV1, the combined A + C excess (T − D) per 1,000 IDs falls by
    ≥ 50 % against FbV1.
  - If K1 passes but K4 does not, the gain is not from the intended mechanism; I report that.

**Verdict:**
- **"Claim competition works":** K1, K2 and K3 pass. It stays opt-in until the user decides on defaults.
- **"Fails":** K1 fails. I report it, with no retuning of claim_q, anchor_apex_s or the rule.
- **"Costs IDs":** K1 passes but K3 fails.

## Caveats (to state in the result)
- One raw file (D) and two entrapment designs; the counts are small (about 100 entrapment hits per arm).
- Claimants come from the first pass, which already contains the neighbour-driven false targets. A false target
  accepted in the first pass becomes a claimant itself, but it can only claim from windows it outranks in another
  family, never its own decoy.
- The first pass is not changed, and neither is the GBM. Only the second pass gets four extra columns.
- Runtime roughly doubles the interference step (two resolver calls per group). Wall time and peak RSS are reported.
