#!/bin/bash
# Submit PLAN_SURVEY follow-up arms P', Pe', H' (V1+pool).
set -euo pipefail
cd /nfs/research/juan/singj/opendia/fasta-20261128
sha256sum -c fasta.sha256
printf 'e1915832a272391e409448a08c813a600058b91230d4d85e84900b92b871047b  opendia-fasta-v4-20261128.sif\n' | sha256sum -c -
test ! -e runs_v4 && test ! -e cache_v4
mkdir -p logs_v4 runs_v4 cache_v4
for a in Pp Pep Hp; do sbatch --parsable -J fastav4-$a -o logs_v4/%x.%j.out arms_v4.sbatch $a; done
squeue -u pst_prd -h -o '%i %j %T'
