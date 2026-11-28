#!/bin/bash
set -euo pipefail
cd /nfs/research/juan/singj/opendia/fasta-20261128
sha256sum -c fasta_fb.sha256.codon
printf 'e1915832a272391e409448a08c813a600058b91230d4d85e84900b92b871047b  opendia-fasta-v4-20261128.sif\n' | sha256sum -c -
test ! -e runs_fb && test ! -e cache_fb
mkdir -p logs_fb runs_fb cache_fb
for a in FbF FbP FbV1 FbPp; do sbatch --parsable -J fastafb-$a -o logs_fb/%x.%j.out arms_fb.sbatch $a; done
