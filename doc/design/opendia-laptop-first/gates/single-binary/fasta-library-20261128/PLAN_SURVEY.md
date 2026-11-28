# PLAN_SURVEY: a better survey signal for the FASTA prefilter (28 November 2026)

**When it was written.** After RESULT_FASTA and its diagnostics (`why/`). Before any run of the variants below. I wrote it because the user asked to try the ideas from the diagnostic: better fragments, interference reduction, pooling and a learned score.

## Why
At K* = 800k, the prefilter keeps 0.60 of F's IDs on D, the IDs the full library finds. The IDs it cuts are seen in the survey, but at noise level. 70 % of them are charge 3-4, and the survey fragments are singly charged y ions without intensities.

## Variants (all opt-in, target-decoy symmetric, defaults unchanged)
Each survey variant is a function applied the same way to the target and the decoy of a pair. Selection stays the top K of a pair statistic that is symmetric in (target, decoy).

| id | survey fragments | co-elution | pair statistic |
|---|---|---|---|
| V0 | theoretical y1 (current) | sum of positive r against the raw consensus (current) | max(T, D) |
| V1 | **predicted** top 6 PeptDeep fragments for every candidate (b/y, fragment charge up to min(2, z-1)) | current | max(T, D) |
| V2 | theoretical y1 | **normalized**: each fragment scaled to its own maximum before the consensus and the correlations, so one intense interferer cannot dominate | max(T, D) |
| V3 | predicted | normalized | max(T, D) |
| +pool | on the best of V0-V3 | | max over charge siblings (same modified sequence) of max(T, D), plus 0.5 × the protein's top-3 mean (R5/R2 offline rules) |
| +learn | on the best of V0-V3 | | 3-fold cross-fitted logistic score on the survey features of the best candidate (co-elution, support, \|apex RT error\|, peaks in window, log apex intensity, FWHM, charge). Folds are by peptide family. Pair statistic = max(f(T), f(D)), with f fit on the other folds. |

Multi-run pooling is not tested: there is one run per arm.

## Metric (D, survey only, codon, same SIF recipe)
- **Primary:** recall at K = 800k of F's target IDs at q ≤ 0.01 (`runs_v2/F`, 25,121). This is a proxy truth, about 1.6 % false.
- **Secondary:**
  - recall of the empirical IDs E_c (26,185);
  - recall of charge 3-4;
  - the label-free excess of target wins over decoy wins among the kept pairs;
  - survey wall time.
- **Held out (H):** the same metrics against FH, F's library searched on H (job 45224807). H is never used to choose.

## Rule (fixed now)
1. Pick the variant with the highest primary recall on D. Ties within 0.005 go to the cheaper variant.
2. Then try +pool and +learn on top of it. Each is adopted only if it raises primary recall on D by ≥ 0.02.
3. The chosen variant is then confirmed on H. It must raise H recall over V0 by ≥ 0.02, or nothing ships.
4. If confirmed, the arms P', Pe' and Fe' run with it on D. For a default change, it must pass PLAN_FASTA criterion 2 (FDP_c(Pe') ≤ FDP_c(Fe') + 0.003) and criterion 4.
5. Criterion 3 is reported against F. Nothing is tuned after the H confirmation.

## Caveats stated now
- F's IDs carry the full search's own biases.
- The learned score uses labels. Cross-fitting keeps a null pair's T and D exchangeable only if the folds hold. Its FDR must be checked by Pe'/Fe' before any default change.
- Data for single-cell or sparse use is not available on codon. Whether these variants help there remains untested.
