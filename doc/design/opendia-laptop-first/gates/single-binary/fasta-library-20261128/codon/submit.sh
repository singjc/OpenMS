#!/bin/bash
# Submit PLAN_FASTA arms. P1/P2 depend on P (reuse its library cache).
set -euo pipefail
cd /nfs/research/juan/singj/opendia/fasta-20261128
ln -sfn ../../entlib-r5-20261118/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.fasta fasta/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.fasta
sha256sum -c fasta.sha256
mkdir -p logs runs cache
p=$(sbatch --parsable -J fasta-P -o logs/%x.%j.out arms.sbatch P)
for a in F Pe Fe H; do sbatch --parsable -J fasta-$a -o logs/%x.%j.out arms.sbatch $a; done
for a in P1 P2; do sbatch --parsable -d afterok:$p -J fasta-$a -o logs/%x.%j.out arms.sbatch $a; done
squeue -u pst_prd -n fasta-P,fasta-F,fasta-Pe,fasta-Fe,fasta-H,fasta-P1,fasta-P2
