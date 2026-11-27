# OpenDIA QPX evidence and quantification outputs, 27 November 2026

Implementation commit: `ba3125b20f`, branch `feature/opendia-clean-rewrite` on `singjc/OpenMS`.

## Scope and compatibility

The output layer is separated from calibration, extraction and scoring. Canonical evidence is a QPX directory of linked Parquet views. High-level precursor tables and quantification matrices are convenience reports, selectable as TSV, Parquet or both. Native OpenMS XIPM is an optional measured-evidence export, not a new QPX XIC implementation.

- `results/Annotations` reads optional gene metadata without changing the compact analysis library or existing library caches.
- `results/Quantification` retains the existing global precursor competition ordering and adds independent peptide/protein-group reporting estimates.
- `results/QpxExport` streams canonical evidence and exact producer candidate diagnostics before releasing scorer feature matrices.
- `results/ReportExport` emits all retained precursor winners, including decoys and failed reporting cuts, and target-only matrices.
- `io/PeakMapExport` replays spectra once after final scoring, while the temporary cache remains alive, and writes native XIPM rows with candidate links.

The main precursor table schema deliberately changes. `-legacy_outputs` writes the original per-run and filtered combined TSVs under `legacy/`, preserving an explicit byte-identity gate. Calibration JSON stays a sidecar. `parameters.ini` records effective inputs and options. New behavior does not alter analysis kernels, candidate winners or precursor identification estimates.

## Reporting policy

Default cuts are `-q 0.01`, `-peptide_q 0.01` and `-protein_q 0.01`. Both run and global estimates must pass. Precursor matrices apply peak-group cuts. Peptide and protein matrices additionally apply only their corresponding level cuts, independently of one another.

A peptide is an unmodified sequence, combining accepted modifications and charge states. Protein groups contain accessions with identical target-peptide evidence in the complete input library, independent of identified hits. A peptide shared across resolved groups contributes to neither protein scoring nor protein area. Quantities are sums of raw integrated MS2 precursor area, without normalization, imputation, MaxLFQ or transfer between runs. Null/NA is distinct from a genuinely measured zero.

Peptide and protein estimates use picked target-decoy competition on the best entity score, with ties to decoys, conservative equal-score blocks and reverse cumulative minima. They are computed before reporting filters. Missing observed counterparts produce no estimate. These are experimental reporting estimators, not a claim of empirically proven true FDR for all data. Entrapment species labels do not enter inference or filtering.

## Evidence format

The reference QPX checkout is pinned to `2dadc0e120cf8b868aab0c4b6489f047efc19561`. The actual upstream Python reader/validator is used, without altering its schema or suppressing findings. OpenDIA supplies opaque candidate identities generated from run, original precursor ID and candidate row, rather than float32-rounded RT. Canonical footer metadata declares `identity_source=supplied`, and records the actual producer composite separately. The upstream validator consequently reports an informational identity-not-checked message. A separate full producer-hash check recomputes identities independently.

Canonical producer extensions use nullable `x_opendia_*` strings, including round-trip representations of exact numeric values. Linked producer candidate diagnostics retain typed values, original numeric precision and scorer features. A technical sample identifier binds each run's LFQ channel without inventing biological sample metadata. Canonical nullable fields do not invent PEP, observed precursor m/z, scan IDs, measured fragment intensities or peptide positions. Predicted library intensities are not exported as observations.

Optional XIPM maps contain actual m/z, RT, available mobility and intensity points, selected by run q-value and padded final winning boundaries. MS1 precursor and search-fragment MS2 maps use calibrated extraction coordinates and actual isolation/mobility rectangles. Only RT-overlapping windows retain points, completed maps are streamed in bounded row groups, and one spectrum replay serves all selected precursors. Partial files are cleaned on exceptions. The ephemeral spectrum cache additionally preserves precursor mobility rectangles for raw/cache equivalence.

## Validation

Final integrated build and all 27 OpenDIA CTests passed. The rebuilt positive and empty canonical fixtures pass the pinned upstream QPX validator without errors or warnings. The independent acceptance harness also passes canonical/typed fixtures and 13 corruption or unavailable-evidence gates.

Astral raw01, seed 0, 16 threads, unchanged search/scoring libraries, with `-output_format both -legacy_outputs -peakmaps_out`:

| Gate | Observed result |
| --- | --- |
| Legacy analysis compatibility | Calibration JSON, per-run precursor TSV and combined precursor TSV all byte-identical to `out_base`, 3/3 |
| Retained evidence | 1,675,700 candidate rows, 336,609 retained winners including 168,192 decoys |
| Upstream QPX conformance | Feature, protein-group, run, sample and provenance views valid, zero warnings/errors, expected supplied-ID informational finding only |
| Independent producer identity | All 1,675,700 feature IDs re-derived, zero stale; separate independent hash implementation also passes |
| Numeric/report lineage | Exact candidate/report/legacy values, TSV/Parquet equality, independent level cuts, matrix eligibility/sums, null/zero and cleanup checks pass |
| Target matrices | 2,221 precursors, 1,652 unmodified peptides, 1,652 supplied-accession protein groups |
| Native XIPM | All 15,799 maps decode through `XIPMParquetFile`: 2,257 precursor MS1 maps and 13,542 MS2 maps, 151,462 actual points; every map links to an exported winner |
| Cleanup | Run-local spectrum cache removed, no orphan partial files |

This Astral fixture uses peptide-derived placeholder protein accessions, so its 1,652 groups are not a biological protein-identification claim. Real multi-protein grouping is exercised by controlled unit tests. Biological protein-level calibration and a realistic annotated multi-run acceptance dataset remain follow-up work. Missing Astral mobility remains NaN/unknown in the XIPM arrays, not a fabricated mobility coordinate.

## Output cost and remaining limits

Local wall time was 11:53.43 and GNU time peak RSS was 8,125,220 KiB, or 7.75 GiB. Inside this same run, final analysis finished after 394.75 s, optional XIPM replay took 5.46 s, candidate diagnostics took 41.46 s, and QPX plus high-level reports took 262.14 s. Output serialization is therefore a material additional cost, not a performance improvement. Avoiding repeated per-candidate peptide/ProForma and metadata construction is a candidate follow-up, not a measured optimization in this gate.

Canonical feature and typed candidate files occupy 181,767,499 and 305,994,388 bytes respectively. Optional XIPM occupies 2,807,075 bytes. Both high-level formats plus legacy files and all evidence occupy about 747 MiB. These deliberately include all candidates and unfiltered winners, unlike the old small filtered combined TSV.

The full validation commands, executable/plan hashes, logs, unmodified-upstream validation reports, independently tested acceptance scripts and native reader check are archived in [`evidence/`](evidence/VALIDATION.md). No runtime or memory improvement is claimed for this output refactor. Local laptop timing is descriptive only. No remote jobs or upstream pushes are part of this gate.
