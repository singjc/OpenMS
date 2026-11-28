# ESTIMATOR_CHECK: FDRBench's own estimator equals eval_fasta's FDP_c (before any FDRBench-arm result)

FDRBench v1.1.1 `-level precursor -r 0.5894 -score score:1` on Fe (custom entrapment, v2), input from `to_fdrbench.py`
(targets only, ENTRAP_ relabelled to the _p_target suffix):

| | FDRBench combined_fdp | eval_fasta FDP_c(r0) |
|---|---|---|
| Fe at q ≤ 0.01 | 0.016290 (n_t 9,214, n_p 56) | 0.016291 (N_human 9,214, N_e 56) |

FDRBench lower bound: 0.0060. A q ≤ 0.05 input cap (`to_fdrbench.py ... 0.05`) gives identical values at q ≤ 0.01
(13 s against 5 min for the full table), so the cap is used for all arms.
