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

**Outputs** (in `-out_dir`):
- `<run>.precursors.tsv`: one row per precursor, with its window, score and run q-value;
- `<run>.calibration.json`;
- `precursors.tsv`: target precursors at run and global q ≤ `-q`, all runs;
- with `-features_out`, `<run>.features.tsv`: every scored window with its features.

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

## Tests

Run `ctest -R 'OpenDIA'` to include the TOPP INI/CTD checks as well as the core tests.
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
