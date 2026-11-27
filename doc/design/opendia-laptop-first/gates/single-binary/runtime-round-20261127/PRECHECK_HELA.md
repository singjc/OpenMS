# Precheck: HeLa diaPASEF identity and memory, structural round (27 November 2026)

This was frozen before submission. There are no tuning knobs. The gate is byte identity, and the timings are descriptive.

## What is tested

- **New binary:** OpenMS `ff4fdc067f` with default options, OpenDIA sha256 `e3533d93709c…`. It includes three changes:
  - fold-parallel cross-fit scorer training (`819205cf35`);
  - the run-local spectrum cache (`eeca94ea64`);
  - extraction in 16 precursor m/z groups, 2 running concurrently at 8 threads (`ff4fdc067f`).
- **Base binary:** OpenDIA sha256 `8ca771eeda87…`, the round-3 base (`ecabd5d0ce` lineage).
- **Container:** both binaries are in one container, `opendia-part-20261127.sif` (sha256 `e94a20c3a30b…`), built on the `opendia-evshare-20261126:tools` image with the current `libOpenMS*.so`. The hashes are in [hela-evidence/](hela-evidence/).

## Run

- **Script:** `part_ab.sbatch` (sha256 `12f59382…`). It runs one sbatch job per set on one node: base r1, new r1, base r2, new r2.
- **Settings:** `-threads 8 -seed 0`, library `hela_r5.predicted_decoys.pqp`.
- **Sets:** D first (r≈7, `…_2737.d`), then H (`…_22474.d`) as a second job.

## Gates

1. Every output file of `new_r1`, `new_r2` and `base_r2` is `cmp`-identical to `base_r1`, and the file sets are identical.
2. No `opendia_cache_*` directory is left in the output directory.

Reported in all cases:
- wall time and peak RSS per arm, from GNU time;
- the Extraction log lines (groups, peak open grids and tensors per group).

This is the first `.d` test of the cache. It checks that the ion-mobility array name and unit, the scalar drift time and the per-frame spectra all round-trip.
