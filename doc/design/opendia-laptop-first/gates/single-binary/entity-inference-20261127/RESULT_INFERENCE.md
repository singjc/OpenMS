# OpenDIA: PyProphet peptide/protein inference (OpenMS f0a62dfbf5)

Date: 2026-11-27. Baseline ba3125b20f. Plan frozen before results: [PLAN.md](evidence/PLAN.md).

## What changed

- The experimental picked peptide/protein estimator is gone. OpenDIA now calls OpenMS
  `LevelContextInference`, the same code OpenSwathInfer uses. It is adapted from PyProphet and
  Rosenberger et al. 2017 (10.1038/nmeth.4398) and runs with its default settings.
- Run-specific, experiment-wide and global p, q and PEP are exported for peptides and proteins.
  `-inference_context run-specific|experiment-wide` picks the reporting gate (default
  run-specific). The global cut always applies as well. Full-precision values are in
  `qpx/opendia/inference.parquet`.
- Peptide = modified sequence, with charge states combined. Protein = supplied accession, scored
  and quantified only by peptides that map to exactly one accession. Target and decoy classes
  stay separate, and no target/decoy pairing is required.
- A context that can't be estimated is NA with a warning, as in PyProphet, which rejects pi0 <= 0.
  There is no fallback estimator.
- **Shared code fix:** in `MultipleTesting::pi0Est`, bootstrap minpi0 took a nearest-rank 10%
  quantile. PyProphet uses `np.percentile(pi0, 0.1)`, the 0.1th percentile with linear
  interpolation. It now matches. This changes bootstrap pi0, so q and PEP, for **all callers,
  including OpenSwathInfer**. A new regression test uses a value taken from pinned PyProphet.
- Precursor scoring and precursor q-values are unchanged.

## Gates

| Gate | Result |
| --- | --- |
| Build plus affected CTests (OpenDIA, MultipleTesting, LevelContextInference, OpenSwathInfer, Peptidoform) | 39/39 pass |
| Pinned PyProphet 4344e2b oracle (actual SQL readers + `error_statistics`), 480+480 peptides x 3 runs, 10 level/context batches | Pass. Rows and scores exact; max p 2.2e-16, q 1.7e-16, PEP 3.7e-15 (limits 1e-12, 1e-10, 1e-2). Before the pi0 fix, 8/10 batches failed (q up to 0.026, PEP up to 0.13). All-tie input: upstream raises pi0 <= 0; native gives NA + 6 warnings |
| Astral raw01, seed 0, 16 threads, both formats, legacy, XIPM | 3/3 legacy files byte-identical to perf-20261127/out_base; upstream QPX validator VALID, 1,675,700 producer IDs, 0 stale; independent checker ACCEPTANCE_PASS; native XIPM 15,799 maps |
| diaPASEF Evosep HeLa 200 ng .d (2737) + ptypic decoy PQP, pre-change vs new binary | 3/3 legacy files byte-identical; QPX VALID; independent checker ACCEPTANCE_PASS (89 checks) |

## Descriptive counts (target entities passing 1%, run-specific + global)

These are not accuracy claims. There is no ground truth here.

| Data | Level | Before | After |
| --- | --- | --- | --- |
| Astral raw01 | peptides | 1,652 (unmodified sequence) | 1,705 (modified sequence) |
| Astral raw01 | proteins | 1,652 | 1,675 |
| diaPASEF | peptides | 22,300 of 54,366 | 24,560 of 59,659 |
| diaPASEF | proteins | 4,768 of 7,033 | 4,576 of 7,033 |

The before and after peptide counts use different identities, so they aren't directly
comparable. The Astral protein accessions are placeholders made from peptides, so only diaPASEF
says anything about real protein inference. No inference warnings occurred on either dataset.

Resources, descriptive only (loaded laptop): Astral 14:33 wall / 7.7 GiB peak vs 11:53 /
7.7 GiB before. diaPASEF after 3:47 / 4.1 GiB vs before 4:39 / 4.1 GiB.

## Caveats

- Reusing the published method doesn't guarantee its FDR is calibrated on OpenDIA scores.
- The oracle needed a SciPy 1.18.1 workaround: `rankdata(method="min")` ranks are cast to
  integers, with equality asserted. No values change.
- The checker's fixed UniMod table gained `UniMod:1 -> Acetyl` for the diaPASEF library. Unknown
  spellings still fail.
- The pre-change diaPASEF binary links today's libOpenMS, so it already has the pi0 fix. That
  doesn't affect precursor identity, because precursor q-values don't use `MultipleTesting`.

Script hashes: check_outputs.py eb48fee5, parity.py b8c57b6d, scipy_compat.py f0f59c6c.
