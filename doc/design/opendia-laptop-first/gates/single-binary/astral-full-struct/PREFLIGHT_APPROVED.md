# Structural round preflight approval, 27 November 2026

Status: APPROVED FOR THE FROZEN SUBMISSION ONLY AFTER THIS AMENDMENT IS PUSHED.
This is not a result. No codon workload for this experiment has run at approval preparation.

## Release gates met

1. **Source.** `a251c23b40cc511781710909759e79b0bd156151` is pushed to `singjc/OpenMS`, branch `feature/opendia-clean-rewrite`.
2. **Frozen plan.** Commit `e9ae73e920` is pushed. The 15 primary cells, unchanged options and acceptance rules remain in force.
3. **Local tests.** All 20 local OpenDIA CTest tests pass. The raw01 seed0 run of this binary matches the three retained base outputs byte for byte; see [evidence/cal2](../runtime-round-20261127/evidence/cal2/). HeLa D and H on codon are byte-identical for this binary (jobs 44865425 and 44865426, [RESULT_HELA](../runtime-round-20261127/RESULT_HELA.md)).
4. **Container.**
   - It is built from the clean, pushed source tree.
   - Local, Docker and SIF OpenDIA SHA256 all equal `6e8f7d02ef9fbc8a18f2809950f73307b9cd823c598a934ef68f10d98371361f`, and the image holds no other OpenDIA binary.
   - Docker and SIF `OpenDIA --help` exit 0.
   - Container-native `OpenDIAWorkflow_test`, `OpenDIALibraryParquet_test`, `OpenDIAScorer_test`, `OpenDIASpectrumCache_test` and `OpenDIACalibration_test` pass ([container-evidence/](container-evidence/)).
   - It was built before the plan push, as [source_provenance.txt](source_provenance.txt) records.
5. **SIF.** SHA256 `7b2ceba837e75f67d3d7cfd53ee864e7f0c5c974fe340b02de12815d21be3033`.
   - Uploaded to the new path `/nfs/research/juan/singj/opendia/part-20261127/opendia-struct-20261127.sif`, and the hash was verified on codon.
   - No existing SIF was overwritten.
6. **Inputs.** The two Parquet caches and three raws are the RT2 preflight's, at the same paths and with the same hashes. They are used unchanged and only read.
7. **Bundle checks.** All 16 evaluator and shell-contract unit tests pass, and all three shell files pass `bash -n`. `compare_resources.py` is byte-identical to the RT2 and compact evaluator.
8. **Board sync.** `board_sync.py` ran before the plan commit ([board_sync_freeze.log](board_sync_freeze.log)) and runs again before this amendment's commit.

## Exact immutable inputs and scripts

**`preflight.sha256`** lists:
- this approval, the frozen PLAN, the source provenance and the evaluator tests;
- `runtime.sha256` and the four runtime scripts;
- the SIF, both Parquet caches and all three raws, at their absolute codon paths.

`submit.sh` verifies it once, before the first `sbatch`.

**`runtime.sha256`** lists only `run.sbatch`, `evaluate.sbatch`, `submit.sh` and `compare_resources.py`. Each primary task verifies it, together with the two-cache `library.sha256` that `submit.sh` writes and checks.

**Rules from here on:**
- Only the fresh `/nfs/research/juan/singj/opendia/astral-struct-20261127` work directory may be created.
- The baseline, compact and RT2 plans, inputs, outputs and jobs stay untouched.
- A partial submission is inspected, not cleaned or repeated.
- Mandatory for every cell: actual executable exit, exact three-file baseline bytes, and exact COMPLETED 0:0 task accounting.
- A re-evaluation for accounting lag is read-only and never a replacement seed.

Record the pushed amendment commit and the job IDs in a later submission record. Do not rewrite this approval or either manifest after submission. No FDR or quantification claim is authorised.
