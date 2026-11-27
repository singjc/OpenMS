# Relabel re-derivation: HeLa r ≈ 7 records and DIA-NN with corrected entrapment labels (27 November 2026)

Frozen before the job runs. This is evaluation only: no new search, no tool change, no rule chosen after seeing the result.

## Why

[DIAG_STEPA](../fdr-20261125/DIAG_STEPA.md) found two labelling errors in the evaluators of every HeLa r ≈ 7 record:

1. **SQL wildcard match.** The SQL `LIKE '%_p_target%'` treats `_` as a wildcard, so it labels human proteins whose accession ends in `P_target` as entrapment. This affects 19 proteins on HeLa.
2. **Human sub-peptides.** The library build re-digests foreign entrapment entries with missed cleavages. 2,648 resulting entrapment peptides (2,912 precursors) are exact human substrings (I = L), so they can truly be present. Step A lists them in [relabel_hela.tsv](relabel_hela.tsv), a copy of `stepA/label_audit_hela.tsv`, sha256 `66411e6d…`.

DIAG_STEPA re-evaluated only the network, GBT and blend arms. The other HeLa records still report old-label bounds, and so does DIA-NN. This job re-derives all of them in one place.

## Evaluators

Only two changes are made, both switched by `LABELS`. Everything else is the frozen code.

| Script | Base | Change |
|---|---|---|
| `eval_pdec_relabel.py` | frozen `pdec/eval_pdec.py`, sha256 `bb03387d…` | adds `LABELS` |
| `eval_diann_relabel.py` | frozen `pdec/eval_diann_r5.py`, sha256 `3fd7b103…` | adds `LABELS` |

- **`LABELS=old`** is the frozen rule. It reproduces the published `EVAL_*.txt`.
- **`LABELS=corrected`** uses `instr(accession, '_p_target')` and counts the RELABEL peptides as human. These are the corrections DIAG_STEPA applied in `stepA/eval_pdec_relabel.py`.

Both change the labels only. Rankings, q-values, r and the bound formulas are unchanged, and r is recomputed from the corrected labels, as the frozen formula defines it.

[test_relabel.py](test_relabel.py) runs 3 synthetic tests, which pass both locally and inside the evaluator container:
- `old` gives output identical to the frozen evaluator;
- `corrected` excludes `P_target` and relabels the listed peptides;
- both hold for DIA-NN too.

## Inputs (read only; sha256 checked on codon 27 Nov)

| Input | Path under `/nfs/research/juan/singj/opendia` | sha256 |
|---|---|---|
| Predicted-decoy library of every r ≈ 7 OpenDIA run | `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp` | `76f2a6c2…` |
| Targets-only library DIA-NN searched | `entlib-r5-20261118/targets.pqp` | `b29d8b43…` |
| DIA-NN D, H reports | `entlib-diann5-20261119/{D,H}/report.parquet` | `1a8f89b4…`, `ab5464cb…` |
| OpenDIA run tables | `entlib-odp5-20261118`, `entlib-sib-20261121`, `entlib-rtr-20261121`, `entlib-gbt-20261122`, `entlib-v7m6-20261123`, `blend-20261125` (B5 arm) | as written by the runs |

**Records covered:** ODP D/H, SIB D/H, RTR D/H, GBT D/H, V7M6 D, BLEND D/H (B5), and DIA-NN D/H. Their G0 and CUR arms are the previous record's default and are already present.

**Not covered:**
- BOUNDS_R5, whose single seed-0 r5 run used the pre-ODP library;
- the r ≈ 1.56 records.

These are listed and not re-derived.

## Job

- **Script:** [relabel.sbatch](relabel.sbatch). It is one separate Slurm job, writes only under `relabel-20261127/` (symlinked layouts and outputs), and verifies `bundle.sha256` first.
- **Container:** `relabel-eval-20261127.sif`, built on `python:3.12-slim` with pandas 2.2.3, numpy 2.2.6, scipy 1.15.3, pyarrow 20.0.0 and sqlite 3.46.1.

## Reported, with no pass or fail rule

1. **Reproduction.** For each record, does `LABELS=old` reproduce the published EVAL table: N, E, bounds and median N@FDP_c 1 %, per seed and pooled? Any mismatch is reported and explained, never silently replaced.
2. **Corrected values.** The corrected lower and combined bounds and median N@FDP_c 1 % per record, with the change against `old`.
3. **DIA-NN.** The corrected DIA-NN bounds and N@FDP_c 1 %, and OpenDIA as a fraction of DIA-NN on corrected labels for the current default (GBT arm).

These are evaluation corrections of existing results. They do not change any verdict made under a frozen rule. Any verdict that would read differently under corrected labels is flagged, not re-decided.
