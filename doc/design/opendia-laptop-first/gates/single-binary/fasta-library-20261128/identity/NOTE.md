# Criterion 1: empirical byte identity (D, seed 0, code a445c1892e)

The `-in D -tr <empirical library>` run matches `entity-inference-20261127/out_pasef_owned` (read-only):

- `legacy/*.precursors.tsv`, `legacy/precursors.tsv` and `calibration.json` are **byte-identical**. These are the files criterion 1 names.
- The top-level `precursors.tsv/.parquet` and `peptides.tsv/.parquet` differ only in the six entity columns (peptide run/experiment/global q and PEP). The largest absolute difference is 2.22e-16, and the row count and every other column are identical. The reason: the reference was written at 23:41 on 2026-11-27, before 3e71a141b9 (00:08, exact PyProphet lambdas). That commit documented this same ≤2.2e-16 change. The FASTA work does not cause it.
