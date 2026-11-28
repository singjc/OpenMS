#!/bin/bash
# Submit PLAN_FASTA arms, v2 (AMEND_CALIBRATION). P1/P2 depend on P (reuse its library cache).
set -euo pipefail
cd /nfs/research/juan/singj/opendia/fasta-20261128
sha256sum -c fasta.sha256
printf '8ce130185fb69cb3500e2f61201519ec08817e4327deb36b17c79479d92f3e5e  opendia-fasta-v2-20261128.sif\n' | sha256sum -c -
test ! -e runs_v2 && test ! -e cache_v2
mkdir -p logs_v2 runs_v2 cache_v2
p=$(sbatch --parsable -J fastav2-P -o logs_v2/%x.%j.out arms_v2.sbatch P)
for a in F Pe Fe H; do sbatch --parsable -J fastav2-$a -o logs_v2/%x.%j.out arms_v2.sbatch $a; done
for a in P1 P2; do sbatch --parsable -d afterok:$p -J fastav2-$a -o logs_v2/%x.%j.out arms_v2.sbatch $a; done
squeue -u pst_prd -h -o '%i %j %T'
