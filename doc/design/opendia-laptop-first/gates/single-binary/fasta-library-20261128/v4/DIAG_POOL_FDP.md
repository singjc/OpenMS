# DIAG_POOL_FDP: where the V1+pool FDP estimate comes from (diagnostic, not a rule; after RESULT_V4)

Scripts: `diag_pool.py` (codon jobs 45292731, 45293190), `overlap.py`, `tdc_check.py`. Output: `diag_pool.json`.

1. **The tool's q-values are plain TDC.** With π0 = 1, (D+1)/T gives the same IDs (Pe' 23,530 against 23,546; Fe
   identical), so no π0 leniency is involved.
2. **The base search is already about 1.6× anti-conservative without any prefilter.** Fe has FDP_c 1.63 % (95 %
   Poisson CI 1.20-2.06 %) at q = 1 %. This is independent of pooling.
3. **The r_lib estimate is dominated by noise.** Stratify the library by each pair's own (unpooled) survey rank into
   20 equal-count bins. In the boost-admitted ranks (about 220k-575k), entrapment/human r_b is only 0.03-0.04, so
   each entrapment hit there stands for about 30 false IDs. Stratified FDP_c = **2.13 % (CI 1.50-2.76 %)**. That
   overlaps Fe's 1.63 % (1.20-2.06 %), so the excess due to pooling is not resolved by these 82 hits.
4. **The protein term causes the entrapment asymmetry; the sibling term does not.** At K = 800k on Pe''s survey:

   | statistic | recall of F | r_sel / r0 |
   |---|---|---|
   | base | 0.623 | 0.98 |
   | sibling only | 0.752 | 1.02 |
   | sibling + 0.1 × protein | 0.785 | 0.77 |
   | sibling + 0.25 × protein | 0.824 | 0.50 |
   | sibling + 0.5 × protein (shipped opt-in) | 0.856 | 0.26 |

5. **Overlap with F.** Of P''s IDs, 13 % of those admitted by their own rank are absent from F's IDs, against 18 % of
   those admitted by the pooling boost. Of F's IDs, 3,826 are absent from P''s.
