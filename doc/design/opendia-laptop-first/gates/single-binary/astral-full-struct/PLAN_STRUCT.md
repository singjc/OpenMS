# Frozen plan: Astral full-library structural round identity and resource experiment

Date: 2026-11-27. Status: FROZEN PROTOCOL. Submission requires the pushed preflight approval.
Source commit: `a251c23b40cc511781710909759e79b0bd156151` (`singjc/OpenMS`, `feature/opendia-clean-rewrite`).

## What changed since the RT2 round (`bdb0159ff4`)

The source adds these commits, each byte-identical on its local and codon checks:

| Commit | Change |
|---|---|
| `ecabd5d0ce` | Removes the failed FDR opt-ins |
| `ba108d3d18` | AVX2 scorer clones |
| `4ca5022f27`, `17f9a98c81`, `042386aa68` | Dead accessors and flat hot-loop arrays |
| `6977925ec3` | Module rename |
| `819205cf35` | Fold-parallel scorer training |
| `eeca94ea64` | Run-local spectrum cache, on by default and removed after the run |
| `ff4fdc067f` | Extraction in 16 precursor m/z groups replaying filtered cache spectra |
| `738f88d0e3` | Interference per extraction group |
| `a251c23b40` | Streaming windowed calibration |

The RT2 records give the evidence for the earlier rounds, and `runtime-round-20261127/` gives it for this one:
- Astral raw01 seed0 is byte-identical against the round-3 base for every commit.
- HeLa diaPASEF D and H are byte-identical on codon.
- 20/20 OpenDIA CTest tests pass.

The local library has 603k precursors. That predicts neither the full-library peak nor the time.

## Scope and immutable inputs

| Name | Path | Note |
|---|---|---|
| WORK | `/nfs/research/juan/singj/opendia/astral-struct-20261127` | Fresh directory; must not exist before submission |
| B | `/nfs/research/juan/singj/opendia/astral-full-20261123` | Frozen attempt-3 baseline, PQP input; primary identity reference |
| C | `/nfs/research/juan/singj/opendia/astral-parquet-20261124` | Compact Parquet result, 15/15 PASS; resource reference |
| RT2 | `/nfs/research/juan/singj/opendia/astral-rt2-20261125` | RT2 result, 15/15 PASS; resource reference |
| CONTAINER | `/nfs/research/juan/singj/opendia/part-20261127/opendia-struct-20261127.sif` | |

**Container contents:** a single OpenDIA binary, SHA256 `6e8f7d02ef9fbc8a18f2809950f73307b9cd823c598a934ef68f10d98371361f`, built from the clean, pushed source.

**Library:** the same two Parquet caches as RT2, symlinked unchanged into `WORK/library` and read directly:
- `C/library/search6.opendia.parquet`, SHA256 `1943673979d97b8da9b87cd6b9df0a48ba5519d72af1537bf64f582f05dafaeb`;
- `C/library/scoring12.opendia.parquet`, SHA256 `425c9e02a3e9fe0a21ccfce44d63b542e1a4664a6d5a9607b0e349a0310603fa`.

**Fixed:** the same raws, seeds, options, 16 cores and 256G as B, C and RT2. Options stay at their defaults, so the new defaults apply: spectrum cache on, 16 extraction groups, parallel = threads / 4. The cache goes to a unique directory inside each run's `out_dir` and is removed when the run ends.

**Not allowed:** tuning, filtering, prediction or library changes. Never modify B, C or RT2 inputs, results, plans or jobs.

## Root release gates

1. The implementation commit is pushed, with local byte identity and all 20 tests passing.
2. `board_sync` runs, and this PLAN is frozen and pushed before any remote action for this experiment.
3. The container holds exactly the committed binary: local, Docker and SIF SHA256 are equal. Container-native `OpenDIAWorkflow_test`, `OpenDIALibraryParquet_test`, `OpenDIAScorer_test`, `OpenDIASpectrumCache_test` and `OpenDIACalibration_test` pass.
4. A preflight amendment `PREFLIGHT_APPROVED.md` is pushed, recording:
   - the commits and evidence;
   - `preflight.sha256`, with absolute codon paths for the container, caches and raws;
   - `runtime.sha256`, listing exactly `run.sbatch`, `evaluate.sbatch`, `submit.sh` and `compare_resources.py`.
5. The frozen bundle is transferred into the fresh WORK. No jobs run before the amendment is pushed.

## Fixed job sequence

The scripts are those of RT2, with only the job names, WORK and CONTAINER changed. `evaluate.sbatch` also runs the unchanged evaluator against RT2.

`submit.sh`:
- verifies the preflight once;
- links and hash-checks the two caches;
- submits the primary array `0-14%3` (raw = floor(i/5) + 1, seed = i mod 5) and an `afterany` evaluation.

Each task verifies `runtime.sha256` and `library.sha256` before it runs.

A partial submission is inspected, never cleaned or resubmitted. Other jobs are never cancelled or touched.

## Evaluation and interpretation

**Primary verdict:** the unchanged `compare_resources.py`, byte-identical to the RT2 and compact evaluator, run against B on the same three files. Each cell PASSes only if all of these hold:
- all three SHA256 match;
- candidate exit is 0;
- baseline exit is 0;
- Slurm accounting is exactly COMPLETED 0:0.

The denominator is fixed at 15. There are no replacement seeds and no exclusions.

**Agreement check:** the same evaluator against C and RT2 must give the same identity verdict as the primary, because both are byte-identical to B.

**Reported:** wall seconds and peak RSS from GNU time, matched per cell against B, C and RT2, with node provenance, plus per-stage seconds and RSS from the OpenDIA logs.

**Expected direction:**
- peak RSS lower than RT2's 52.8 GiB median;
- shorter calibration, extraction and scoring stages.

These expectations are not acceptance criteria. Resources are descriptive only: Slurm assigns nodes, and filesystem caches are not controlled. The spectrum cache adds about 12 bytes per peak of temporary disk under each `out_dir` while a run is active, and its peak size is reported. No FDR or quantification claim is made.

## Local bundle checks

`python3 -m unittest discover -s <this directory> -p 'test_*.py' -v` gives 16/16, and `bash -n` passes on all three shell files. These are local synthetic checks only.
