# OpenDIA

Identification-free DIA analysis on a predicted spectral library. Targets and decoys are treated identically; no labels
or external identifications are used.

## One command: `OpenDIA`

```
OpenDIA -in run1.raw run2.raw ... -lib library.pqp [-scoring_lib library12.pqp] -out_dir results -threads 8
```

Per run:
1. **Calibration**, identification-free: library RT to run RT map and tolerance, fragment m/z offset and fitted ppm
   tolerance, 1/K0 line, peak-boundary cap. A windowed second pass runs when the whole-run pass finds too few anchors.
   The search, grid and interference use 0.6 × the calibrated 1/K0 half width (`-im_tolerance_scale`; PLAN_IMNARROW,
   +4.3 % IDs at 1 % true FDP on the r ≈ 7 HeLa entrapment, calibration unchanged). `-im_tolerance_scale 1` restores
   the calibrated width.
2. **Candidate search** with the calibrated library: native-scan XIC peak picking, up to 5 peaks per precursor.
3. **Population**: every candidate of each target-decoy pair whose best candidate has at least 3 supporting fragments.
4. **One more pass over the spectra** collects the scoring evidence: 6-fragment traces (each fragment within ±10 ppm,
   `-trace_tolerance_ppm`) and apex b/y ions, and the 12-fragment grid of every window. A grid is held only while its
   window is open.
5. **49 scorer features**:
   - search 4, boundary 5, b/y ions 2, 12-fragment evidence 6, fragment deficit 3, shared-signal interference 4;
   - evidence block 12 (core fragment and MS1 mass error, MS1 co-elution, isotope and area, MS2 area, ion-mobility
     centroid, co-elution and share);
   - matched filter 9 (`mt_`; per-fragment robust z, leave-one-out correlation, cosine, dominance, LLR);
   - charge-state siblings 4 (`s_n`, `s_co`, `s_best`, `s_best_d`, from the first scoring pass; see step 6).
6. **Scoring**:
   - `OpenDIA::CrossFitScorer`, three folds by peptide family. Seeds are the calibration anchors' windows and decoy
     windows. The scorer is a small network with a decoy-whitened readout and two self-training rounds, then one window
     per precursor, decoy-tail calibrated scores and pooled paired target-decoy q-values.
   - Sibling rescoring: the sibling features are computed from those out-of-fold scores and the network scores again.
   - Hard-negative GBM rescoring: three cross-fitted histogram gradient boosting rounds on the final features, trained
     on q ≤ 0.01 targets against all windows of the hardest 2 % of decoy precursors. Its scores and q-values are final.
7. **RT refine**: a ridge RT correction on sequence composition is fitted from the run's own q ≤ 0.01 IDs, then steps
   2-6 run again with corrected RTs and a narrower RT window (`-no_rt_refine` skips it, about half the wall time).

Decoys come from the library, or with `-decoys predict` (`-decoy_model_dir`, `-decoy_cache`) pseudo-reversed decoys are
predicted like the targets.

**Library format.** `-lib` and `-scoring_lib` take a PQP or an OpenDIA Parquet library (`.parquet`). A PQP of at least
100,000 fragments is converted once, without any option: `<stem>.opendia.parquet` is written next to it (temporary name,
then rename) and later runs read it instead while the PQP keeps its size and modification time. If the directory is
not writable the PQP is read each run. Both formats give the same in-memory library, so results are byte-identical.
The Parquet file (schema version 1, zstd) has one row per fragment in library order: precursor index, ID, TRAML_ID,
decoy, m/z, charge, library RT and 1/K0, sequence, proteins, transition ID, product m/z, intensity, ion type, ordinal,
charge and intensity rank within the precursor; precursors without fragments are one row without a fragment.
`-decoys predict` still needs a PQP.

**Library memory.** Parquet inputs and valid automatic caches are decoded one row group at a time into a compact
analysis table. Search fragments retain only m/z, intensity and the numeric transition ID. Precursor IDs and protein
sets are interned, and folds, pair hashes, decoy partners and sibling groups are computed once. A separate scoring
library is streamed into the twelve lowest transition IDs per precursor, not retained as another library. Calibration
and trace collectors read indexed coordinates. Search applies calibration and RT corrections on the fly, without
copying the library for either pass. RT refinement keeps feature rows only for its training IDs and predicts the
remaining rows individually. The first PQP read still materializes the source library to create the automatic cache.

Across runs, a precursor's global q-value comes from its best run score, by paired target-decoy competition.

## Results and evidence

QPX is the canonical low-level evidence dataset. It is a directory of linked Parquet views, not a single `.qpx`
file. High-level reports are separate convenience views. Scoring, candidate selection, RT boundaries and raw MS2
areas are unchanged by the export layer.

| Output in `-out_dir` | Contents |
| --- | --- |
| `qpx/opendia.feature.parquet` | Every candidate peak, including decoys, with QPX identities, peptide metadata, boundaries, scores and q-values |
| `qpx/opendia.pg.parquet` | Run-level target and decoy protein-group evidence, quantities and independent protein-group q-values |
| `qpx/opendia.run.parquet`, `qpx/opendia.sample.parquet`, `qpx/opendia.provenance.parquet` | Run/channel bindings, unannotated injection sample identifiers, calibration and processing/threshold metadata |
| `qpx/opendia/<run>.candidates.parquet` | All candidate scorer features and full-precision diagnostic values, linked by `feature_id` |
| `precursors.tsv` | Every retained precursor winner across all runs, including decoys and rows failing reporting cuts |
| `precursor_matrix.tsv` | Target precursor raw MS2 areas passing both run and global peak-group cuts |
| `peptides.tsv`, `peptide_matrix.tsv` | Unmodified-sequence evidence and target peptide quantities |
| `protein_groups.tsv`, `protein_group_matrix.tsv` | Protein-group evidence and target protein-group quantities |
| `<run>.calibration.json` | Full calibration sidecar |
| `parameters.ini` | Effective CLI parameters, including libraries, seed, threads and reporting cuts |
| `<run>.xipm`, with `-peakmaps_out` | Optional measured MS1/MS2 peak maps for visualization |
| `<run>.features.tsv`, with `-features_out` | Optional legacy-style scorer-feature dump |
| `legacy/`, with `-legacy_outputs` | Original per-run and filtered combined precursor TSVs, for exact comparisons |

`-output_format tsv|parquet|both` chooses the high-level table/matrix formats. The default is TSV. QPX and XIPM
remain Parquet regardless of this option. The precursor table includes original precursor ID, modified and
unmodified sequence, charge, proteins, available genes, RT apex/start/end in seconds, raw MS2 area, discriminant
score, decoy status, run/global precursor and peptide/protein-group q-values, and reporting pass flags. Missing
annotations and unestimable level q-values are `NA` in TSV or null in Parquet, not invented values.

### Configurable reporting filters

`-q 0.01`, `-peptide_q 0.01` and `-protein_q 0.01` are the defaults. Both run-level and global q-values must pass
the corresponding cut. These options filter reports and quantities, not the underlying candidate evidence.

- Precursor quantities require the peak-group cuts.
- Peptide quantities additionally require peptide cuts. A peptide is an **unmodified sequence**, summed across its
  accepted charge states and modifications.
- Protein-group quantities require the peak-group and protein-group cuts, independently of the peptide reporting cut.
  Indistinguishable groups are defined by identical
  target-peptide evidence in the input library. Peptides shared across distinct groups do not contribute to protein
  scoring or area, avoiding double counting.
- Peptide and protein-group q-values are separate picked target-decoy reporting estimates, computed from the best
  eligible entity score before report filtering. They are not copied from precursor q-values. Missing observed
  target/decoy counterparts have no level estimate. These new level estimators are experimental, not an empirical
  guarantee of the stated true FDR on every dataset.
- Quantities are **raw integrated MS2 areas summed without normalization, imputation or MaxLFQ**. Missing matrix
  cells are `NA`/null. A genuinely measured zero stays zero. Matrices contain targets only.

OpenDIA supplies opaque QPX feature IDs generated from run name, original precursor ID and candidate row, avoiding
rounded-RT collisions. The producer composite is recorded in the footer and checked independently of QPX's default
identity recipe. Full-precision numeric `opendia_*` columns remain in the linked producer candidate diagnostics.
Canonical QPX uses standards-compliant nullable string `x_opendia_*` extensions for producer-specific values, with
round-trip numeric representations where needed, alongside its standard float32 fields. Nonwinning candidates retain
their features and diagnostics, but no measured intensity or unmasked final
score is fabricated. Unavailable observed m/z, PEP, scan IDs and fragment measurements remain absent/null. The
library's predicted fragment intensities are never presented as measured fragment evidence.

### Optional visualization evidence

`-peakmaps_out` writes native OpenMS `.xipm` peak maps from a separate spectrum replay after final scoring, while
the run-local cache is still available. `-peakmaps_q` selects run-level precursor winners (default 0.01),
`-peakmaps_rt_padding` adds seconds around each winning RT boundary (default 15), and `-peakmaps_decoys` includes
passing decoys. Maps contain actual m/z, RT, available mobility and intensity samples for precursor MS1 and search
fragment MS2 channels. This is optional additional I/O, not another scoring pass. Without the flag there is no
peak-map extraction overhead. XIPM is separate from QPX for now, with candidate identifiers linking the evidence.

**Raw files** need a .NET 8 runtime (set `DOTNET_ROOT` if it is not in a standard location); mzML works without it.

**Cost** on one Astral 15 min raw file (8 threads): about 8-10 min and 9.6 GB peak.
- Search: 30-40 s.
- Evidence pass: 5-6 min. It holds about 90,000 open 12-fragment grids and the 6-fragment trace tensor of 1.6 M windows.
- Scorer: 1-2 min.
- Calibration adds 1-2 passes; the RT refine second pass repeats search, evidence and scoring.

## Interference features

Co-eluting precursors in one isolation window often share fragment m/z. A candidate that borrows a strong ion from its
neighbour then looks better than it is, and a real peptide whose fragments overlap a neighbour looks worse.

`OpenDIA::InterferenceResolver` links candidate peaks of different precursors that:
- are in the same isolation window;
- overlap in RT;
- share a fragment within the fitted tolerance.

For every link, the claimant's abundance is estimated from its own unshared fragments. Its elution profile gives how
much of the shared signal it explains in the candidate's window, and that part is removed. Four features describe the
result:
- `i_explained`: the explained share of the central signal;
- `i_libcos_resid`: the library cosine of the residual;
- `i_dlibcos`: the change in that cosine;
- `i_nclaim`: log(1 + claimants).

**Validation** (`gap-20261114`, PLAN_GAP arm I; Astral entrapment raws 01/02, 12-fragment DIALibGen library):
- **Python prototype:** added to the prototype's scorer, it raised N at entrapment FDP 1 % from 1694 / 1628 to 1879 / 1921 (frozen seed).
- **C++ implementation:**
  - +143..+185 (raw01) and +260..+367 (raw02) at 1 % over five training seeds, against +118..+185 and +199..+351 for the Python features.
  - At 2 % the raw01 gain depends on the seed: +14..+131.
  - Entrapment FDP_c at q ≤ 0.01 stayed ≤ 0.0071.
- **Benchmark-specific share:** about a third of the raw02 gain comes from same-composition entrapment/parent pairs that only exist in the benchmark.

## Module layout

`OpenDIA.cpp` is the command-line entry point. Core headers in `include/` and implementations in `source/` use
the same module directories, with no `OpenDIA` prefix on header or implementation filenames:

| Module | Files |
| --- | --- |
| `workflow/` | `Workflow` |
| `library/` | `Library`, `LibraryParquet` |
| `library/preparation/` | `Assay`, `DecoySequences`, `AssayPrediction`, `Preparation` |
| `io/` | `SpectrumCache`, `PeakMapExport` |
| `results/` | `Annotations`, `Quantification`, `QpxExport`, `ReportExport` |
| `calibration/` | `Calibration`, `RTMap`, `RTRefine` |
| `extraction/` | `Search`, `Traces`, `FragmentGrid`, `PeakPicker`, `MzIndex` |
| `features/` | `PeakFeatures`, `Interference` |
| `scoring/` | `Scorer`, `GBM` |

Includes are qualified relative to `include/`, for example `#include <extraction/Search.h>` and
`#include <library/Library.h>`. Header-only components and implementations without a separate header remain in
their module directory, so the mirrored layout does not require artificial header/source pairs.

Library preparation owns separate assay records, pseudo-reverse sequence generation and direct PeptDeep inference.
`Preparation` composes these into the existing predicted-decoy workflow. It has no analysis-state or file-format
ownership. The CLI alone adapts PQP records through `TransitionPQPFile`. Future FASTA prediction and assay generation
can compose the same preparation components without coupling them to the compact analysis library. These future
entry points are not implemented here.

Tests remain flat in `test/`, with filenames such as `Workflow_test.cpp`, `FragmentGrid_cases.inc` and
`no_openswath.cmake`. Build targets and CTest names retain their `OpenDIA` prefixes, independently of source filenames.

## Tests

Run `ctest -R 'OpenDIA'` to include the TOPP INI/CTD checks as well as the core tests.
- **Output tests:** `OpenDIAQuantification_test` checks independent hand-counted level competition, shared peptides,
  missing counterparts and raw-area aggregation. `OpenDIAQpxExport_test` and `OpenDIAReportExport_test` inspect
  linked IDs, schemas, precision, null/zero handling, decoys, matrices and atomic publication.
  `OpenDIAPeakMapExport_test` round-trips real sampled points through the native XIPM reader.
  `OpenDIAAnnotations_test` checks optional PQP/Parquet gene metadata.
- **`OpenDIAPreparation_test`:** sequence and annotation eligibility against the legacy decoy generator.
  With ONNX, **`OpenDIAPreparation_prediction_test`** also compares real-model assay predictions and fresh generated
  PQP bytes against the legacy preparation path. Legacy dependencies occur only in these test oracles.
- **`OpenDIAWorkflow_test`:** compact-table grouping, top-twelve scoring, identities, moves and partial-cache fallback
  against an independent materialized library oracle, including shuffled/interleaved fragments.
- **`OpenDIACalibration_test` and `OpenDIATraces_test`:** indexed and Library-based collector equivalence, including
  fragment ties, repeated precursor windows, mobility and b/y-ion callbacks.
- **`OpenDIARTRefine_test`:** byte-identical ridge fits, out-of-fold predictions and final predictions from full-library
  and training-only feature matrices, with an identical all-library vocabulary.
- **`OpenDIAScorer_test`:**
  - gradients against finite differences;
  - the readout in closed form;
  - paired competition on a hand example;
  - tail calibration monotonicity;
  - thread independence and recovery of a synthetic population.
- **`OpenDIAFragmentGrid_test`:** the 12-fragment and deficit formulas against numpy reference values, and the streaming grid collector against a brute-force fill.
- **`OpenDIAInterference_test`:** the explained share against an all-pairs brute force.
