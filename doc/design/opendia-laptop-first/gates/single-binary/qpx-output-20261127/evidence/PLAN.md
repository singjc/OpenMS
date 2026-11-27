# OpenDIA QPX output acceptance plan

Date: 2026-11-27
Source baseline: 86bb3c6654 (owned library preparation)
QPX specification snapshot: 2dadc0e120cf8b868aab0c4b6489f047efc19561
Saved reference executable: OpenDIA.ref in this directory

## Contract
- No change to extraction, scores, candidate selection, run/global precursor q-values, boundaries or MS2 areas.
- Canonical QPX feature, protein-group, run and provenance Parquet views, plus exact candidate feature diagnostics.
- Unfiltered retained-winner precursor table and target-only precursor, unmodified-peptide and protein-group matrices.
- Default run/global level cuts 0.01, independently configurable.
- Optional native OpenMS XIPM with actual MS1/MS2 points, no invented fragment observations.
- Legacy output flag allows exact baseline comparison in a separate legacy directory.
- Missing quantities/annotations/q-values are null/NA, not zero or fabricated estimates.

## Gates
1. Integrated OpenDIA and all output test targets compile.
2. Full OpenDIA CTest suite, including independent level FDR examples and native XIPM round-trip.
3. Actual upstream QPX Python schema, semantic, identity and conformance checks on synthetic and representative outputs. Any upstream limitations are recorded, not hidden.
4. Astral raw01, seed 0, 16 threads, unchanged library/scoring library, final CLI exports both TSV/Parquet, legacy output and optional XIPM.
5. Legacy per-run and combined TSVs and calibration JSON byte-identical to out_base.
6. Independent Python/Arrow checks for complete candidate IDs, exact retained numeric values, TSV/Parquet consistency, matrix eligibility, sums and null/zero behavior.
7. Cache directory removed, no orphan partial files, XIPM links to exported candidate IDs, one replay regardless of precursor count.

Local timing and RSS are descriptive only on this loaded laptop. No controlled speedup claim.
New peptide/protein reporting estimators remain experimental, not a claim of empirically established true FDR.
References are read-only. No remote jobs or upstream pushes are part of this task.
