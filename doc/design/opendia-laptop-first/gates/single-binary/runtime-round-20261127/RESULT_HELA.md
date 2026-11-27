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

## What is left for memory (no change made)

- Per-window state that scales with the window count: `InterferenceResolver::Window` (core area and 32-bin profile), the feature table, `RunWindow` and the search candidates. Shrinking these would mean computing interference per isolation-window group and not materialising every window at once. That is a larger restructuring.
- On Astral, the peak is calibration pass 2, with 600k members' traces held over the run (about 10.4 GB).
