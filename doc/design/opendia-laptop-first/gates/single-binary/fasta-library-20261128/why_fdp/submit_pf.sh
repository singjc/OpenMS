#!/bin/bash
set -euo pipefail
cd /nfs/research/juan/singj/opendia/fasta-20261128
sha256sum -c fasta_fb.sha256.codon
printf "e1915832a272391e409448a08c813a600058b91230d4d85e84900b92b871047b  opendia-fasta-v4-20261128.sif\n" | sha256sum -c -
test ! -e runs_pf && test ! -e cache_pf
mkdir -p logs_pf runs_pf cache_pf
for a in FbV1_g1 FbP_g1 FbV1_s1 FbV1_s2; do sbatch --parsable -J fastapf-$a -o logs_pf/%x.%j.out why_fdp/arms_pf.sbatch $a; done
