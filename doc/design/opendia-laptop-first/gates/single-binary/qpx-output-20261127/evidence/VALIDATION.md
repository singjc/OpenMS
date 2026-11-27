# Reproducing the output acceptance gate

The source change is `ba3125b20f` on the fork's `feature/opendia-clean-rewrite`. The QPX specification and validator checkout is pinned to `2dadc0e120cf8b868aab0c4b6489f047efc19561`.

The scripts in this directory are the actual local gate artifacts, not a second implementation of QPX validation. Paths in the frozen runner and validator reflect this workstation. Use new output directories when adapting them. Do not modify the reference files.

## Native build and tests

```bash
cmake --build "$BUILD" --target OpenDIA OpenDIAQuantification_test OpenDIAQpxExport_test OpenDIAReportExport_test OpenDIAPeakMapExport_test OpenDIAAnnotations_test OpenDIASpectrumCache_test -j4
ctest --test-dir "$BUILD" -R OpenDIA --output-on-failure -j4
```

The five new CTests cover quantification, optional annotations, canonical/diagnostic QPX, reports and native XIPM. Existing spectrum-cache tests include exact precursor mobility metadata roundtrips. Prediction tests retain their real-model gate.

## Representative analysis

`run_acceptance.sh` runs Astral raw01, seed 0, 16 threads, unchanged search/scoring libraries, both high-level formats, legacy TSVs and optional XIPM. It refuses an existing output directory, records the executable and frozen plan SHA256, compares three reference files byte-for-byte and verifies cache/partial cleanup.

## Independent validation

Install the unmodified pinned QPX checkout in an isolated Python environment. The actual invocation is:

```bash
python validate_opendia_qpx.py OUTPUT/qpx --json-output conformance.json
```

This calls upstream `Dataset.validate(strict=True, semantic=True, check_identity=True)` and fails on errors or warnings. Supplied opaque feature IDs produce an expected upstream informational finding, not a warning. The script additionally re-derives every producer candidate ID with upstream hashing and reports stale IDs explicitly.

A separate environment containing DuckDB, PyArrow and NumPy runs:

```bash
python check_outputs.py OUTPUT --refdir BASELINE
python check_outputs_selftest.py
```

The independent harness checks TSV/Parquet equality, complete candidate lineage and exact numeric values, legacy comparisons, target-only matrix membership and sums, independent peptide/protein cuts, null versus zero, producer and protein-group hashes, XIPM binary arrays and candidate links, and output cleanup. Its positive/typed fixtures and 13 corruption or unavailable-evidence cases test the checks themselves. This is not a validation of empirical true FDR or a replacement for native XIPM decoding.

Finally, compile `read_native_xipm.cpp` with the source/build OpenMS include paths, link the current `libOpenMS`, and run it on the actual exported `.xipm`. It loads every map through `XIPMParquetFile`, checks metadata, paired array lengths, finite measured values and monotone RT, and reports map/precursor/point counts. The native unit test also verifies exact values and raw/cache mobility-window equivalence on controlled fixtures.

The archived evidence is intentionally compact. Large spectra, cached spectra and generated result datasets remain scratch artifacts rather than repository fixtures. Timings are descriptive measurements on a loaded laptop, not a controlled performance comparison.
