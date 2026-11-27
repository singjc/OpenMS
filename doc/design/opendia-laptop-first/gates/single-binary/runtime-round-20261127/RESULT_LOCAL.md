# Local result: runtime round 3 (cache locality), 27 November 2026

**Verdict: every arm byte-identical to the base on Astral raw01 seed0 (3/3 files each, 6 arms). 19/19 OpenDIA CTest tests pass after each commit.** This is a local laptop diagnostic with the 603k-precursor library and `-scoring_lib library12_full.pqp`, 16 threads. The laptop was under other load (rustc and others) during several arms, and arms ran one after another, not interleaved. Timings are descriptive only. Byte identity is the result.

## Source (OpenMS `feature/opendia-clean-rewrite`, fork singjc/OpenMS)

| Commit | Change |
|---|---|
| `ecabd5d0ce` | Removes the failed FDR opt-ins (`-gbm_blend`, `-evidence_share`, `-gbm_exclude`, `-gbm_monotone`, `-gbm_fold_local`, `-bulk_null_q`). Their records stay in `fdr-20261125/`. Default path unchanged |
| `ba108d3d18` | Scorer: `readout` decoy covariance loop interchange (inner loop over contiguous j), `__restrict` on the backward pass, `target_clones("avx2","default")` of `embed`, `readout` and `objective` (GCC, x86_64, Linux). With `-ffp-contract=off` there is no FMA, so every sum keeps its order. New test: `readout` identical for 1 and 16 threads |
| `4ca5022f27` | Dead code: `TraceCollector::releaseWindowIds()` and `evidenceEnabled()` (no callers) |
| `17f9a98c81` | `FragmentGridCollector` keeps apex RT per window-index entry and precursor 1/K0 in flat arrays (was a load from a 224-byte window struct and a 384-byte precursor struct per step). `CalibrationCollector` keeps each member's last RT in a flat array (was a pointer chase into `member.rt.back()`) |
| `042386aa68` | `TraceCollector` RT range searches run on sorted flat apex arrays parallel to `by_rt_` and to each m/z bucket (was a 128-byte `Sample` load per binary-search step). MS1 hit lists and MS2 isolation windows reuse member scratch vectors |

## Profiles that drove the changes

- **CPU (`perf record`, clean binary):** the scorer kernels took about 65 % of CPU: `WhitenedNetwork::embed` 24 %, `readout` 20 %, `objective` 20 %. Those loops were scalar SSE2. The AVX2 clones give 4-wide doubles with the same operation order.
- **LLC load misses by symbol (`perf record -e LLC-load-misses`, sc2 binary, `-no_rt_refine`, [miss_top.txt](evidence/miss_top.txt)):** `FragmentGridCollector::consumeSpectrum` (serial part) 28 %, `CalibrationCollector::consumeSpectrum` 24 %, `TraceCollector::consumeMS1_` 8 %. The scorer is below 1 % each.
- **Annotated hot instructions:**
  - In the grid collector, the misses fell on the `windows_[window_index_[k]].rt_apex` load (the open test and the bin test) and on the first grid-cell add.
  - In calibration, 57 % of the function's misses were on the `member.rt.back()` load: one random heap line per active member per isolation window.

These are the latency patterns from "It's the memory latency, stupid": one useful double per cache line, behind an index. The fix keeps the same values in a dense array that the loop walks in order. The grid-cell add misses remain: each window owns a 12 × 32 double grid, so that store is inherently scattered.

## Local measurements (raw01 seed0, 16 threads, GNU time and `perf stat`)

| Arm | Binary sha256 (12) | Wall | User s | Sys s | Cal consume p1 / p2 s | Extraction 1 / 2 s | Scoring 1 / 2 s | LLC load misses (core, 1e9) | Peak RSS GB |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| base | `8ca771eeda87` | 21:38 | 4781 | 955 | 46 / 59 | 97 / 195 | 261 / 413 | 6.8 | 12.8 |
| sc2 (`ba108d3d18`) | `4c1e48a24209` | 16:29 | 3458 | 814 | 60 / 72 | 135 / 111 | 221 / 191 | 6.7 | 12.8 |
| loc1 (`17f9a98c81`) | `051e78fa0fdc` | 13:58 | 3234 | 784 | 48 / 49 | 101 / 93 | 177 / 191 | 5.5 | 12.7 |
| loc2 (flat trace apex) | `44d76c976bbf` | 14:01 | 3288 | 786 | 38 / 52 | 99 / 109 | 178 / 195 | 5.0 | 12.9 |
| loc3 (`042386aa68`) | `e4dd72785710` | 13:29 | 3221 | 774 | 36 / 48 | 98 / 90 | 190 / 184 | 5.1 | 12.7 |

Against the base, loc3 uses 33 % less user CPU and 19 % less system time, and has 25 % fewer core LLC load misses. Peak RSS is unchanged: it is set by the 600k-member pass-2 calibration panel, which this round did not change.

## Rejected: `OMP_WAIT_POLICY=active` (environment only, no code)

The arm `loc1act` (loc1 binary, active wait) was byte-identical. System time fell from 784 s to 46 s and voluntary context switches from 40 M to 0.2 M. But user CPU rose from 3234 s to 17097 s, and wall time rose from 13:58 to 18:12, because spinning threads compete with the reader and with other processes. The default passive policy stays.

## Not done in this round

- The fragment grids stay double. A float grid would halve the 6 GB grid memory at HeLa r ≈ 7, but it changes sums, so it needs its own opt-in and validation.
- The calibration panel layout (about 300-byte members with four heap trace vectors each) sets the Astral peak. A structure-of-arrays layout is the next memory candidate.
- The 40 M voluntary context switches come from OpenMP barriers in per-spectrum parallel loops (about 300k spectra × several regions × 16 threads). Merging those regions into fewer, larger regions is the next runtime candidate.
- Codon speed and identity on HeLa D r ≈ 7 is not measured for this round.

## Evidence

The [`evidence/`](evidence/) folder contains:

- `run_arm.sh`, the exact run script;
- `base_outputs.sha256`, the sha256 of the base outputs;
- `miss_top.txt`, the LLC-miss profile;
- per arm: the binary sha256, GNU time output, `perf stat` CSV, the gzipped OpenDIA log, and `identity.txt` (`cmp` against base, 3/3 same).
