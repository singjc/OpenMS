# PLAN_STRUCT submission record, 27 November 2026

| Item | Value |
|---|---|
| Frozen plan | docs `e9ae73e920` |
| Preflight approval | docs `b0b89bcca3`, pushed before any remote action for this experiment |
| WORK | `/nfs/research/juan/singj/opendia/astral-struct-20261127`, created fresh at transfer |
| Primary array `0-14%3` | Slurm `44880255` |
| Evaluation, `afterany` | Slurm `44880256` |

At transfer `submit.sh` verified the whole preflight once:
- `runtime.sha256`;
- every bundle file;
- the SIF, both Parquet caches and all three raws, each at its absolute codon path.

It then linked both caches and checked `library.sha256` before the first `sbatch`. All checks were OK; see [submission-evidence/submit.log](submission-evidence/submit.log). Both caches and the three raw files are read only.

No other job was cancelled or touched. The verdict follows the frozen rule in [PLAN_STRUCT.md](PLAN_STRUCT.md).
