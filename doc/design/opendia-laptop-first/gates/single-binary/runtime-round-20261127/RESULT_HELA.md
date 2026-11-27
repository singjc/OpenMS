# Result: HeLa diaPASEF identity and memory, structural round (27 November 2026)

**Verdict: PASS on both gates for D and H.**

- **Identity:** every output file of `new_r1`, `new_r2` and `base_r2` is `cmp`-identical to `base_r1`, and the file sets are the same.
- **Cache cleanup:** no `opendia_cache_*` directory is left behind.
- **Scope:** this is the first `.d` test of the spectrum cache. The mobility arrays, scalar drift times and per-frame spectra round-trip exactly.

The plan is [PRECHECK_HELA.md](PRECHECK_HELA.md) (frozen in `3eb30fd255` before submission). Codon jobs 44786934 (D, node cn-087) and 44786935 (H, node cn-090) both ran on AMD EPYC 9555 nodes. Arms alternated on one node, with 8 threads and seed 0. The logs are in [hela-evidence/codon_logs.tar.gz](hela-evidence/codon_logs.tar.gz).

## Results (GNU time; wall and peak RSS are same-node A/B, descriptive)

| Set | Arm | Wall r1 / r2 | Peak RSS GiB r1 / r2 | Extraction pass 1 / 2 s (r1) | Peak held tensors pass 1 |
|---|---|---|---|---|---|
| D r≈7 | base `8ca771ee` | 14:39 / 14:32 | 29.9 / 29.7 | 208 / 186 | 1.98 M |
| D r≈7 | new `e3533d93` | **10:57 / 11:08** | **17.0 / 17.1** | **111 / 98** | 151k per group (2 in flight) |
| H | base | 19:36 / 20:35 | 21.6 / 22.0 | 304 / 301 | 674k |
| H | new | **14:00 / 16:04** | 21.4 / 21.5 | **154 / 141** | 71k per group (2 in flight) |

- **D:**
  - Peak RSS is 43 % lower (29.8 to 17.1 GiB on the means).
  - Wall time is 25 % shorter (14:36 to 11:03).
  - Extraction is 1.9× faster.
  - The peak used to be the extraction tensors and grids (1.98 M windows open at once). Now it is 17.4 GB at extraction and features, and the per-window arrays of 10.1 M windows dominate.
- **H:**
  - Wall time is 25 % shorter (20:06 to 15:02), and extraction is 2× faster.
  - Peak RSS is about unchanged (21.8 to 21.5 GiB). H held only 674k tensors at once. Its peak (about 21.9 GB at pass-1 extraction in both arms) comes from the arrays of 12.1 M windows: search candidates, `RunWindow`, the interference windows (384 B each) and the feature table. Calibration on H peaks at 14.7 GB.

## Follow-up: interference per extraction group (OpenMS `738f88d0e3`)

Extraction groups now hold whole isolation windows. The interference features, which link windows only within one isolation window, are computed inside each group. Only the 12-fragment central areas (96 B per window) outlive the group, so the run-wide `InterferenceResolver::Window` array (384 B per window) is gone.

The check ran in codon jobs 44834711 (D, cn-072) and 44834712 (H, cn-091), with the script in [ifg-evidence/](ifg-evidence/). Each job alternated `prev` (`e3533d93`) and `ifg` (`8cd6186d`) on one node.

Result: every `ifg` and `prev` output is `cmp`-identical to `base_r1` of the A/B above, on both D and H. Astral raw01 seed0 is also byte-identical, both with the cache (16 groups) and with `-no_spectrum_cache` (a single group); see `evidence/ifg*`.

| Set | prev wall r3 / r4 | ifg wall r1 / r2 | prev peak RSS GiB | ifg peak RSS GiB |
|---|---|---|---|---|
| D | 11:03 / 11:11 | 11:03 / 10:54 | 16.9 / 16.6 | 16.1 / 15.9 |
| H | 13:29 / 13:10 | 13:08 / 13:12 | 21.4 / 19.4 | 21.8 / 18.6 |

- **D:** peak RSS is about 0.8 GiB lower. Wall time is unchanged.
- **H:** peak RSS varies by about 2 GiB between repeats of the same binary, which is larger than the change.
- **Astral (laptop):** wall time fell from 9:00 to 7:18. Features plus interference now take 0.5 s plus 32 s summed over the groups, against 12.7 s single-threaded before, and the groups run in parallel. This timing is descriptive.
- **Where the peak is now:**
  - On both D and H it is the per-window state: the feature table, `RunWindow` and the search candidates, at extraction and in the first scorer fold (H r1 reached 21.4 GB there).
  - On H, calibration peaks at 14.7 GB.

## Follow-up: streaming windowed calibration (OpenMS `a251c23b40`)

The change applies to the windowed (second) calibration pass, which Astral uses with a panel of 300k targets and 300k decoys:
- **Search on close:** a member's stage-3 search runs as soon as the run passes its collection window, batched and in parallel, instead of after the whole pass.
- **Trace cut:** the member's traces are then cut to the ±64-spectrum slice around its best apex, which is all that stage 4 reads.
- **Deferred reserve:** the trace capacity is reserved at the member's first spectrum instead of at panel setup, so only open members hold it.

Ranges, traces and candidates are the same as before.

| Check | Result |
|---|---|
| Astral raw01 seed0 (laptop), byte identity | 3/3 SAME ([evidence/cal2/](evidence/cal2/)) |
| Astral calibration pass-2 RSS | 10.9 to 7.8 GB |
| Astral run peak RSS | **10.6 to 7.9 GB** (now set by calibration pass 1); wall 7:11 (descriptive) |
| HeLa D, codon 44865425 (cn-005) | 3/3 SAME against `D_base_r1`, file set same; 11:26, 16.4 GiB |
| HeLa H, codon 44865426 (cn-007) | 3/3 SAME against `H_base_r1`, file set same; 14:14, 21.2 GiB |

HeLa uses only the whole-run pass, so its memory is unchanged. The run checks that the refactored stage-3 path still gives identical output there. The scripts and hashes are in [cal-evidence/](cal-evidence/).

## What is left for memory (no change made)

- Per-window state that scales with the window count: the feature table (45 floats per window), `RunWindow` with its `SearchCandidate` (about 100 B) and the search candidates before the population filter. Shrinking these needs a lighter `RunWindow` or features built per group. Both change internal layouts, but not outputs.
- On Astral, the peak is now the whole-run calibration pass 1. Its 80k-member panel holds traces over the full run (about 7.2 GB at the end of the pass), because stages 1 and 2 search the whole run. On H, pass 1 peaks at 14.7 GB for the same reason.
