# Local result: structural round (fold parallel, spectrum cache, grouped extraction), 27 November 2026

**Verdict: byte-identical to the round-3 base on Astral raw01 seed0 for the cache and grouped arms (3/3 files each). Grouped extraction at 1 group (serial) and at 64 groups (8 concurrent) is byte-identical too. 20/20 OpenDIA CTest tests pass.**

- **Setup:** laptop, 603k-precursor library, `-scoring_lib library12_full.pqp`, 16 threads.
- **Timings:** descriptive only. The laptop ran other work, and the arms ran one after another.
- **HeLa diaPASEF:** checked on codon under [PRECHECK_HELA.md](PRECHECK_HELA.md), frozen before submission.

## Source (OpenMS `feature/opendia-clean-rewrite`, fork singjc/OpenMS)

| Commit | Change |
|---|---|
| `6977925ec3` | Stage2 module renamed to `OpenDIAPeakFeatures`, and code identifiers lose the `stage2` prefix. The calibration JSON key `"stage2"` is unchanged |
| `819205cf35` | Cross-fit scorer: the folds train concurrently, with threads split evenly between them. Each fold has its own `z` and `round_score`, its log lines are printed in fold order, and a fold's error is rethrown |
| `eeca94ea64` | Run-local spectrum cache. The first raw read writes MS1 and MS2 to a flat binary file under `<out_dir>/opendia_cache_<unique>`, and later passes replay it through a prefetch thread. It is removed when `run()` exits, including on exceptions. `-no_spectrum_cache` turns it off and `-spectrum_cache_dir` sets its location |
| `ff4fdc067f` | Extraction in precursor m/z groups (details below) |

Details of `ff4fdc067f`:
- **Grouping:** windows are split into groups contiguous in precursor m/z, with about equal window counts. The default is 16 groups, with threads/4 running concurrently.
- **Replay:** each group replays all MS1 spectra, but only the MS2 spectra whose isolation window overlaps its m/z range.
- **Results:** per-window results are copied out. Interference isolation indices are rebuilt from the cache index.
- **Options:** `-extraction_groups`, `-extraction_parallel`.

## Why grouped extraction is byte-identical

Every per-window extraction value depends only on the window's own precursor and on the spectra that isolate its precursor m/z:
- the 6-fragment trace tensor, MS1 evidence, b/y counts and the 12-fragment grid;
- this holds for every value except the MS1 sums, which use all MS1 spectra, and each group gets all of them.

A group therefore sees every spectrum that any of its windows would have used, in the same RT order, and computes the same sums in the same order.

Each collector keeps its own close order and slot pools. These affect memory only, not values.

The steps that couple windows run after the groups merge, on the full window set as before:
- the noise median;
- interference, which links windows only within one isolation window;
- sibling features.

The isolation-window list is the set of valid isolation windows of non-empty, finite-RT MS2 spectra, the same set the single collector's lanes held. It is rebuilt from the cache index.

## Local measurements (raw01 seed0, 16 threads, GNU time)

| Arm | Binary sha256 (12) | Wall | User s | Sys s | Extraction 1 / 2 s | Scoring 1 / 2 s | Peak held tensors 1 / 2 | Peak RSS GB |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| loc3 (round 3) | `e4dd72785710` | 13:29 | 3221 | 774 | 98 / 90 | 190 / 184 | 118k / 96k | 12.7 |
| fold (`819205cf35`) | `6b99f3ff3aef` | 12:07 | 3298 | 547 | 102 / 97 | 139 / 141 | 118k / 96k | 13.2 |
| cache (`eeca94ea64`) | `8cd8c5269cf8` | 11:56 | 3259 | 585 | 112 / 85 | 147 / 154 | 118k / 96k | 10.9 |
| part (`ff4fdc067f`) | `e581501325df` | **9:00** | 3061 | 480 | **40 / 32** | 141 / 135 | **15k / 11k per group** | 10.6 |

Against the round-3 base (21:38, 12.8 GB), wall time is now 2.4× shorter. On these arms:

- **Grouped extraction:**
  - Extraction is 2.4-2.8× faster, because 4 groups use the 16 threads better than one collector whose serial parts limit it.
  - Peak held tensors per group fall 8×. With 4 groups in flight, about 60k tensors are held against 118k.
- **Peak RSS:** set by calibration pass 2 (600k-member panel, about 10.4 GB), not by extraction. Extraction ends at 6.2 GB RSS on Astral. The memory gain from grouping shows where extraction set the peak, as on HeLa D r≈7 (34.7 GB): see the codon check.
- **Cache:** read time per later pass falls from 15-22 s to 1.2-2.1 s. The cache is 1.6 GB.

Extremes (`-no_rt_refine`, 16 threads; the two outputs are byte-identical to each other):

| Groups / parallel | Wall | Extraction s | Peak held tensors per group | Peak RSS GB |
|---|---:|---:|---:|---:|
| 1 / 1 (serial, single collector) | 5:25 | 71 | 118k | 10.6 |
| 64 / 8 | 4:21 | 30 | 6.9k | 10.8 |

## Evidence

`evidence/{fold,cache,part,nr_g1,nr_g64}/` holds the OpenDIA log (gzip), the GNU time output and the binary SHA256 for each arm, and `identity.txt` (cmp against `out_base`) where one was taken. The fold arm's identity was checked at the time, and its outputs were cleaned before this record. The cache and part arms include fold.
