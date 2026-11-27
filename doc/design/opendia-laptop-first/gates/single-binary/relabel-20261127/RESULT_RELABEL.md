# Corrected HeLa and DIA-NN evaluation, 2026-11-27

## Outcome

The separate evaluation job 44889904 completed with exit 0:0 in 12:56. All 26 evaluations succeeded: 13 old-label reproductions and 13 corrected-label evaluations. The scope is 100 existing OpenDIA runs and two DIA-NN reports. No search, scoring, ranking, q-value calculation or primary output was changed.

Old-label evaluations reproduce the historical per-seed and pooled records. The old-mode outputs also match the first attempt byte for byte. The corrected table is in [SUMMARY.md](SUMMARY.md), with complete evaluator output and hash inventory in [out/](out/), and execution evidence in [logs/](logs/).

## Correction

The corrected rule uses exact `instr(accession, '_p_target') > 0`, not the case-insensitive SQL wildcard `LIKE '%_p_target%'`. The frozen list of 2,912 precursor entries matching human substrings, with I=L, is labelled human. The library entrapment ratio changes from 7.0539 to 6.9176. Labels are evaluation-only and never enter training or ranking. No new similarity exclusion was introduced.

## Current default and DIA-NN

OpenDIA values are medians over five seeds. DIA-NN has one existing report per set.

| Set | Tool | N at combined FDP 1%, old | Corrected | Own-q lower FDP, corrected | Own-q combined FDP, corrected |
|---|---|---:|---:|---:|---:|
| D | OpenDIA default GBT | 21,396 | 21,805 | 0.94% | 1.07% |
| D | DIA-NN | 26,006 | 26,391 | 0.79% | 0.90% |
| H | OpenDIA default GBT | 37,342 | 37,897 | 1.12% | 1.28% |
| H | DIA-NN | 41,272 | 41,835 | 0.82% | 0.93% |

OpenDIA/DIA-NN N at combined FDP 1% changes from 82.3% to 82.6% on D and from 90.5% to 90.6% on H. Absolute calibration improves for both tools, but the relative identification gap changes little. H OpenDIA's lower estimate remains above 1%. These evaluation corrections are not evidence of a species-independent FDR fix.

ODP, RTR, SIB, GBT, V7M6 D and BLEND are all re-derived in the full table. Frozen historical verdicts are retained, not retroactively re-decided using corrected labels. The older pre-ODP BOUNDS_R5 library and r~1.56 records are outside this frozen scope.

## Reproducibility and failed first attempt

[PLAN_RELABEL.md](PLAN_RELABEL.md) was frozen before submission. Attempt 1, job 44884353, completed but `singularity --cleanenv` removed label-selection environment variables. Its supposed corrected outputs were therefore old-mode outputs. This was detected from the mode header, and those outputs were not accepted as corrected results. They remain under [attempt1/](attempt1/).

[AMENDMENT_1.md](AMENDMENT_1.md) was committed and pushed before attempt 2. Explicit container environment arguments and output-mode guards fixed the wrapper without changing evaluation rules. Attempt 2 has 13 verified old-mode headers and 13 corrected-mode headers. Retrieved output checksums match the remote manifest. The evaluator SIF SHA256 is `717d3282ee6e557b081d11b8633ace217f6650e0a8475d59e1a82a86fb283f61`.

## Separate structural gate and cleanup

The structural full-library gate is independent: array 44880255, 15 fixed raw/seed cells, followed by evaluation 44880256. At 17:41 UTC, three cells had completed successfully, three were running, and nine were pending. No acceptance claim is made before the fixed-denominator identity evaluation finishes.

The authorized old-container cleanup removed 50 SIF files, about 17.55 GB. References, DIA-NN tools, library-generation tools, retained RT2 and today's active containers were preserved. The exact deletion inventory is [deleted_sif.tsv](deleted_sif.tsv).
