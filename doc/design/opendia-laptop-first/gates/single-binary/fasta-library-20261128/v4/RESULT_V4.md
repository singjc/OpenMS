# RESULT_V4: P', Pe', H' with the V1+pool survey (28 November 2026)

Code 143d3ac7c6 (`-survey_predicted -survey_pool`, default K* = 800k). SIF v4 e1915832...047b. Jobs 45271258-60.
References: F and Fe from v2 (K = 0, so the survey options do not apply), FH for H'. Evaluator: `eval_v4.py`. It
runs the frozen `eval_fasta.py` and adds the AMEND_SURVEY_POOL estimator. The script was committed before results
(072113f3fc) and smoke-tested on v2 data, where it reproduced RESULT_FASTA exactly.

## Results (D unless noted; N = target precursors at run q ≤ 0.01)
| arm | N | N_e | FDP_c(r0) | FDP_c(r_lib) | N@FDP_c 1 % (r0 / r_lib) | wall | peak RSS |
|---|---|---|---|---|---|---|---|
| P' | **24,924** | | | | | 1,688 s | 15.9 GB |
| Pe' | 23,546 | 82 | 0.0094 | **0.0254** | 23,639 / 21,371 | 2,169 s | 25.4 GB |
| Fe (ref) | 9,270 | 56 | 0.0163 | 0.0160 | 8,543 / 8,563 | 7,453 s | 57.4 GB |
| F (ref) | 25,121 | | | | | 3,981 s | 36 GB |
| H' (run H) | 47,513 | | | | | 1,920 s | 23.7 GB |

- Ratios: Pe' r0 = 0.5902 (its candidates; V1 drops pairs without survey fragments, so r0 differs slightly from Fe's
  0.5894). Pe' r_lib = 101,674 / 638,605 = **0.159**. Fe r_lib = 0.606.
- Compared with v2 (V0, no pool): P 15,917 → 24,924; H 37,830 → 47,513 (0.899 of FH's 52,847).

## Criteria
| criterion | value | result |
|---|---|---|
| 1 no regression (defaults) | identity_v4_compare.txt: 3 byte-identical, 15 same as v1 | pass |
| 2 FDR, frozen estimator r0 | 0.0094 ≤ 0.0163 + 0.003 | pass |
| 2 FDR, AMEND_SURVEY_POOL r_lib | 0.0254 > 0.0193 | **fail** |
| 3 IDs | P'/F = 0.992 ≥ 0.97; N@FDP_c 23,639 (r_lib 21,371) ≥ 0.97 × 8,543 | pass |
| 4 cost | 1,688 / 3,981 = 0.42 ≤ 0.5 | pass |
| goal N_P ≥ 24,669 | 24,924 | met |

**By the amended rule, written before these runs, V1+pool does not qualify as a default.** Criterion 2 must pass
under both estimators, and it fails under r_lib.

## Interpretation (not a rule)
- The prediction in AMEND_SURVEY_POOL was correct. Pooling removed entrapment from the library much faster than
  human peptides: r_lib 0.159 against r0 0.590. With r0, each entrapment hit stands for about 2.7 false hits; with
  r_lib, about 7.3. The truth lies between the two:
  - r0 is optimistic, because absent human peptides are boosted and entrapment peptides are not;
  - r_lib is pessimistic, because the library also contains the true human peptides, which inflates the denominator.
- So the FDP of V1+pool at q = 0.01 lies between 0.9 % and 2.5 %. The run's own q-values (target-decoy) are
  unaffected by this asymmetry, because pooling is symmetric between targets and decoys. What is unresolved is
  whether that symmetry is enough, or whether pooling lets wrong human peptides through faster than decoys do.
- That question is open: decoys inherit the boost through the pair statistic (max(T, D)), so a decoy competes for
  its target's pooled rank. Entrapment cannot show whether the FDR is controlled here, because pooling breaks the
  equivalence between entrapment and absent human peptides that the entrapment estimate assumes.
- One way to resolve it without retuning is an entrapment design that pools like human peptides do. For example,
  shuffled entrapment proteins that carry a fraction of real human peptides, so that entrapment proteins also get
  protein evidence. That would be a new plan, frozen before it is run.

## Options for the user (no default change made)
1. Ship `-survey_predicted` (V1) without pooling as the default prefilter, after its own P/Pe run. V1 alone is
   target-decoy symmetric and pooling-free, so r0 remains valid. Survey recall is 0.706 against 0.596.
2. Keep V1+pool opt-in (`-survey_predicted -survey_pool`) and document the FDR caveat.
3. Default K = 0 (no prefilter): 25,121 IDs, 3,981 s, 36 GB.
4. Run a new frozen plan with pooling-aware entrapment to settle criterion 2 for V1+pool.
