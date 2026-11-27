# OpenDIA-owned library preparation, 27 November 2026

Implementation commit: `86bb3c6654`, branch `feature/opendia-clean-rewrite` on `singjc/OpenMS`.

## Scope

The private OpenDIA implementation now owns the parts of predicted-decoy preparation used by `-decoys predict`. Headers and implementations mirror one another under `include/library/preparation/` and `source/library/preparation/`:

- `Assay`: library-building compounds, transitions, proteins and modifications. These records do not depend on the compact analysis library or on OpenSWATH records.
- `DecoySequences`: fixed pseudo-reversal with terminal K/R switching, modification relocation, collision reservation and the historical annotation eligibility exclusions.
- `AssayPrediction`: direct low-level PeptDeep RT, CCS and MS2 inference, including model resolution and fragment mapping.
- `Preparation`: composition of sequence generation, instrument/NCE estimation, RT and mobility mapping, fragment selection, intensity normalization and the historical 80% compound/protein coverage check.

The CLI retains `TransitionPQPFile` as its sole OpenSWATH include, at the PQP conversion edge. The prior target-only PQP write/read round trip is preserved because it canonicalizes identifiers and numeric values. Cache naming, cache keys, temporary-file cleanup, options and output schemas are unchanged. Shared OpenMS/OpenSWATH implementations are not changed or removed. Legacy classes remain in the preparation tests as independent oracles, not as core dependencies.

This is the minimal fixed policy already used by OpenDIA, not a general decoy-generation option surface. Future FASTA-to-library prediction and assay generation can compose the preparation records and predictor. Those future entry points are not implemented here.

## Validation

- OpenDIA and the new preparation test executable build successfully. Sequence eligibility, metadata and order match the legacy generator. Real-model predictor fields and a newly generated target/decoy PQP match exactly.
- A separate actual CLI gate used the saved pre-change executable and the new executable with distinct, initially empty output/cache directories. Both generated 24 decoys. Their PQPs are byte-identical, SHA256 `ffc95e7192ab7cb742641ff663c62c10de4773fbdbc656ca5421f2808fd3e17c`. Both then exit 8 at the same expected calibration precondition because the fixture has only 24 target-decoy pairs and three target fragments. This validates preparation and CLI conversion, not successful analysis on that tiny fixture. Prediction temporary files were removed.
- All three new implementation files compile with an isolated ONNX-disabled config overlay. The predictor constructor and preparation entry point both throw `Precondition` as intended. This is a targeted compile/runtime smoke check, not a second complete no-ONNX OpenMS build.
- The recursive architectural guard passes on the actual source. Isolated negative fixtures reject an extra OpenSWATH CLI include and a forbidden dependency nested under `source/library/preparation/`. Symbol inspection finds no undefined `OpenSwathLibraryPreparation`, `MRMDecoy` or `PeptDeepLibraryPredictor` references in the OpenDIA executable.

- Representative Astral raw01 seed0 full analysis completed successfully. All three files are byte-identical to the read-only pre-refactor reference, the file set is identical, and the spectrum cache was removed. The loaded-laptop run took 9:51.36 with peak RSS 8,104,264 KiB. This is an identity check, not a performance comparison.

- Real-model oracle tests cover CCS enabled and disabled, exact fragment/compound/protein fields, internal and terminal modifications, missing-mobility sentinels, retained short decoy assays, and independent compound/protein coverage failures. Fresh generated PQPs are byte-identical in both mobility modes. The oracle explicitly uses the old CLI's `predict_keep_short=true` and four inference threads.

All 22 OpenDIA CTests passed in the final integrated run, including the two new preparation tests, architecture guard, CLI rejection checks and TOPP INI/CTD checks. Logs, binary/output hashes and validation harnesses are preserved in [evidence/](evidence/).
