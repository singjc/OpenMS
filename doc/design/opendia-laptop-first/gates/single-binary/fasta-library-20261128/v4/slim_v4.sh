#!/bin/bash
# Run on codon (pst_prd): slim copies of FINISHED arms (those with a .time file) for local evaluation. Columns are the
# ones eval_fasta.py reads; sha256 of the full files. Usage: slim_arms.sh ARM...
set -euo pipefail
cd /nfs/research/juan/singj/opendia/fasta-20261128
out=${TMPDIR:-/tmp}/fastav4_slim_$$; mkdir -p $out/runs_v4
for a in "$@"; do
  test -s runs_v4/$a.time || { echo "skip unfinished $a" >&2; continue; }
  f=runs_v4/$a/precursors.tsv
  sha256sum $f >> $out/full.sha256
  mkdir -p $out/runs_v4/$a
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} {print $c["precursor_id"]"\t"$c["proteins"]"\t"$c["discriminant_score"]"\t"$c["run_q_value"]"\t"$c["decoy"]}' $f | gzip -1 > $out/runs_v4/$a/precursors.slim.tsv.gz
  cp runs_v4/$a.log runs_v4/$a.time $out/runs_v4/
done
cp logs_v4/*.out $out/
tar -C $out -cf - .
rm -rf $out
