#!/bin/bash
# Run on codon (pst_prd): slim copies of FINISHED arms (those with a .time file) for local evaluation. Columns are the
# ones eval_fasta.py reads; sha256 of the full files. Usage: slim_arms.sh ARM...
set -euo pipefail
cd /nfs/research/juan/singj/opendia/fasta-20261128
out=${TMPDIR:-/tmp}/fastav2_slim_$$; mkdir -p $out/runs_v2
for a in "$@"; do
  test -s runs_v2/$a.time || { echo "skip unfinished $a" >&2; continue; }
  f=runs_v2/$a/precursors.tsv
  sha256sum $f >> $out/full.sha256
  mkdir -p $out/runs_v2/$a
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} {print $c["precursor_id"]"\t"$c["proteins"]"\t"$c["discriminant_score"]"\t"$c["run_q_value"]"\t"$c["decoy"]}' $f | gzip -1 > $out/runs_v2/$a/precursors.slim.tsv.gz
  cp runs_v2/$a.log runs_v2/$a.time $out/runs_v2/
done
cp logs_v2/*.out $out/
tar -C $out -cf - .
rm -rf $out
