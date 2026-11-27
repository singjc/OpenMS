# Full-library structural round: identity PASS, 2026-11-27

## Frozen verdict

**PASS: 15/15 fixed cells**, three Astral raws times five seeds. Every task in array 44880255 completed with exit 0:0. Evaluator 44880256 completed with exit 0:0 in 52 seconds. No replacement runs or exclusions were made.

For every cell, `precursors.tsv`, the run-specific precursor TSV and calibration JSON are byte-identical to the frozen original baseline. The unchanged evaluator also gives 15/15 against the compact-Parquet and RT2 references. Thus all 45 output-file comparisons pass against each reference. The denominator remained 15.

The source is `a251c23b40cc511781710909759e79b0bd156151`. [PLAN_STRUCT.md](PLAN_STRUCT.md), frozen in `e9ae73e920`, and [PREFLIGHT_APPROVED.md](PREFLIGHT_APPROVED.md), pushed in `b0b89bcca3`, preceded submission. Inputs, binary hashes, container checks and submission evidence are retained with the plan.

## Descriptive resources

GNU time values, all 15 cells, 16 allocated CPUs. GiB means KiB / 1,048,576. The speedup column is the median of paired baseline/candidate wall ratios, not a ratio of medians.

| Version | Median wall, min | Median peak RSS, GiB | Median paired speedup to current | Median paired RSS reduction |
|---|---:|---:|---:|---:|
| Original baseline | 94.60 | 87.91 | 4.60x | 44.4% |
| Compact Parquet | 37.79 | 64.47 | 1.90x | 24.2% |
| RT2 | 35.09 | 52.84 | 1.73x | 7.2% |
| Current structural round | **20.77** | **48.76** | 1.00x | 0.0% |

Current GNU wall times range from 18.70 to 25.99 minutes. The matched cells used different Slurm nodes and uncontrolled filesystem caches. These are descriptive comparisons, not controlled causal speedups. No FDR or quantification claim follows from the resource result.

The structural round makes a substantial runtime improvement, but **does not achieve 2x lower whole-run memory** on this full library. Relative to RT2, the median peak reduction is only about 7.7% when computed from the two medians, or 7.2% as the median paired reduction. The full-library scorer now sets the peak, so reduced extraction/calibration state does not translate into a proportional process-peak reduction. The smaller-library and HeLa memory results must not be extrapolated to this workload.

## Current per-stage timing

Medians across the same 15 runs, seconds. Pass 2 here means the RT-refined search/extraction/scoring repeat, not calibration's internal pass 2.

| Stage | First pass | RT-refined pass |
|---|---:|---:|
| Search | 87.670 | 81.948 |
| Grouped extraction, including group-local interference | 216.397 | 201.654 |
| Scoring | 159.802 | 148.401 |

Calibration total median is 74.908 seconds. Stage medians are not additive run totals. Full per-stage RSS, peak RSS, fold training and node provenance remain in [runs/](runs/). Scoring time varied across runs: raw01 seed4 took 505.727 seconds for first-pass scoring, with no stall or restart.

The cache's recorded size is 1,625-1,640 MB (median 1,629 MB, using the program's `cache_mb` units). A post-completion remote scan found **zero `opendia_cache_*` directories** under the 15 run directories. Cache cleanup succeeded.

## Evidence

- [Primary comparison](results/comparison.json): original frozen baseline, 15/15.
- [Compact comparison](results/comparison_compact.json): compact Parquet, 15/15.
- [RT2 comparison](results/comparison_rt2.json): RT2, 15/15.
- [Evaluator exit codes](results/evaluation_exit_codes.txt): all three comparisons return 0.
- [Slurm accounting](results/slurm_accounting.tsv), [job logs](logs/), and [per-run evidence](runs/).

Remote work directory: `/nfs/research/juan/singj/opendia/astral-struct-20261127`. Reference directories were read-only throughout. The separate corrected-label HeLa/DIA-NN evaluation is recorded in [RESULT_RELABEL](../relabel-20261127/RESULT_RELABEL.md).

## Next performance target

The remaining full-library peak is scorer state, not raw replay or active grids. A next output-identical change should bound or compact scorer working sets, with a fresh frozen identity/resource gate. Float32 training or changed scoring semantics must remain opt-in and are not justified by this identity result.
