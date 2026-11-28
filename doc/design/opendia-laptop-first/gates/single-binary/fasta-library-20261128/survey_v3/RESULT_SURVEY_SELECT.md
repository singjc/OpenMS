# RESULT_SURVEY_SELECT and AMEND_SURVEY_POOL (28 November 2026)

## Survey-only results (D, codon; truth = F's IDs, 25,121; recall at K = 800k)
| variant | base | +pool | +learn | +learn+pool | charge 3-4 (base / +pool) | E_c (base) | wall |
|---|---|---|---|---|---|---|---|
| V0 current | 0.596 | 0.757 | 0.613 | 0.766 | 0.39 / 0.72 | 0.581 | 6:32 |
| V1 predicted fragments | **0.706** | **0.878** | 0.712 | 0.879 | 0.54 / 0.87 | 0.659 | 14:27 |
| V2 normalized co-elution | 0.598 | 0.756 | 0.617 | 0.773 | 0.39 / 0.73 | 0.581 | 5:04 |
| V3 both | 0.708 | 0.880 | 0.709 | 0.883 | 0.54 / 0.88 | 0.659 | 12:43 |

Held out H (truth = FH's IDs, 52,847): V0 base 0.704, V0+pool 0.759, V1 base 0.802, **V1+pool 0.851**.

V0 on the variant code (ac2107c592) reproduces the P survey byte for byte.

## Rule application (PLAN_SURVEY, as frozen)
1. Highest base recall is V3 at 0.7079; V1 has 0.7056. The difference of 0.0023 is within 0.005, so the tie goes to the cheaper variant. That is **V1**: V3 is V1 plus normalization, so it does strictly more work. The wall times are noisy across nodes; V2 also ran faster than V0.
2. **+pool** on V1 adds 0.172 (≥ 0.02), so it is adopted. **+learn** adds 0.006 (< 0.02), so it is not adopted. +learn+pool over +pool adds 0.001.
3. **On H**, V1+pool reaches 0.851 against V0's 0.704, a gain of 0.147 (≥ 0.02). Confirmed.
4. Next: P' (human, V1+pool, K* = 800k) and Pe' (human + entrapment, V1+pool) on D. F and Fe stay as the references, because the survey options do not apply at K = 0.

Findings, not rules:
- Normalized co-elution (V2) changes almost nothing.
- The learned target/decoy score helps only slightly. On H it lowers recall (0.673 against 0.704), which looks like a train/test shift.

## Amendment written before P'/Pe' run: the entrapment estimate under pooling
Pooling is target-decoy symmetric, so the q-values stay valid. A null pair's statistic depends only on max(T, D) of itself and of label-free neighbours, so swapping T and D does not change the selection.

Pooling is not symmetric between **absent human peptides and entrapment peptides**, though. An absent human peptide gets a boost from its present charge siblings and from its present proteins. An entrapment peptide never gets one, because its proteins are absent. So among false candidates, the prefilter keeps human ones at a higher rate than entrapment ones. The entrapment share of false IDs then falls below r0/(1+r0), and **FDP_c with r0 from the candidates underestimates FDP**.

Therefore, for Pe', both of these are reported:
- FDP_c(r0), the frozen estimator;
- FDP_c(r_lib), with r_lib = entrapment/human target precursors in Pe's final library. r_lib is at most the entrapment/human ratio among kept false candidates, because the final library also holds the true human peptides. So FDP_c(r_lib) is the conservative one.

Criterion 2 for a default change must pass with **both** estimators, against FDP_c(Fe). For Fe there is no prefilter, so r_lib ≈ r0. This makes criterion 2 stricter, not looser.
