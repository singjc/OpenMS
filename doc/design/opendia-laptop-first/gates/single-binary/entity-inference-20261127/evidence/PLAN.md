# Frozen acceptance plan: native entity inference

Date: 2026-11-27. Baseline OpenMS ba3125b20f0505445946049a739ce36ce4006197.

## Method and scope

Reuse OpenMS LevelContextInference with its default ErrorEstimationConfig, adapted from PyProphet and Rosenberger et al. 2017, DOI 10.1038/nmeth.4398 (PMC5581544). Python oracle pinned to PyProphet 4344e2b5cb30871085eb33b953210791b920624d. No decoy-pair requirement and no fallback estimator or parameter retuning.

Canonical modified peptide sequences collapse charge states, not modifications. Keep target and decoy classes separate and retain the actual decoy sequences. Protein identities are supplied accessions plus class. Only peptides uniquely mapped to one supplied accession contribute protein scores and areas. Protein groups are singleton accessions, not invented indistinguishable groups.

Compute best entity scores from retained finite scorer winners before reporting cuts. Compute run-specific, experiment-wide and global statistics. Default report eligibility combines run-specific and global cuts. Experiment-wide selection is configurable and still requires global cuts. Unestimable contexts produce null statistics and explicit warnings, without substituting another method. Precursor scoring and run/global q values remain unchanged.

## Fixed gates

1. Build OpenDIA and all unit tests. Run all registered OpenDIA CTests.
2. Multi-run synthetic independent oracle: use actual pinned PyProphet entity SQL readers and error_statistics/lookup routines. Exact entity, context and score rows. Absolute numerical tolerances fixed before observing results: p-values 1e-12, q-values 1e-10, PEP 1e-2, reflecting existing native regression precision. Include modified sequences, charge collapse, ties, shared proteins, sparse runs and actual decoys.
3. Synthetic output gates: independent peptide/protein thresholds, selected-context eligibility, null versus zero, exact inference companion statistics and TSV/Parquet agreement.
4. Full Astral raw01 seed0, 16 threads, both report formats, legacy outputs and native XIPM in a new unique directory. Compare all three legacy outputs byte-for-byte against perf-20261127/out_base. No identity changes accepted in calibration JSON, per-run precursor TSV or pooled precursor TSV.
5. Validate actual QPX with unmodified upstream Python validation pinned to 2dadc0e120cf8b868aab0c4b6489f047efc19561 and independently check producer IDs and report lineage. Native XIPM roundtrip check. Entity statistics and corresponding report eligibility intentionally differ from the baseline.
6. Record wall/RSS descriptively on loaded laptop, no causal runtime claim. Astral protein accessions are peptide-derived placeholders, not biological protein inference validation. Method reuse is not an empirical FDR calibration guarantee for OpenDIA data.

All previous output/reference directories are read-only. No remote jobs needed. No retuning after seeing results.

7. Added 2026-11-27 22:57 at user request, before any run: diaPASEF Evosep HeLa 200ng .d (2737) with 200623_hela_evosep_fractions_library_ptypic_decoy.pqp (146,788 precursors, 72,866 decoys, real protein accessions). Run the saved pre-change binary OpenDIA.before and the new binary, seed 0, 16 threads, legacy outputs, into new folders. Legacy precursor TSVs and calibration JSON must be byte-identical. Then run check_outputs.py on the new output. This library has real protein accessions, so peptide/protein counts are reported descriptively. No accuracy claim: there is no ground truth for this file. Caveat: the before binary links today's libOpenMS, so it runs the fixed pi0 too. That doesn't matter for precursor identity, since precursor q-values never use MultipleTesting.

Environment note: installed SciPy 1.18.1 returns floating rank storage for method=min, which the pinned PyProphet uses directly as NumPy indices. The external oracle shim casts only those ranks to integer after exact value-equality assertion. No rank values, tie rules, statistical parameters or pinned upstream files change. Record this compatibility adjustment in results.
