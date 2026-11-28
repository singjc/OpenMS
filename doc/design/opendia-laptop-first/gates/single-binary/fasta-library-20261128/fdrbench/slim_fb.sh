#!/bin/bash
# Run on codon: slim tables of finished FDRBench arms (eval_fb.py columns plus sequence/charge for FDRBench's own FDP).
set -euo pipefail
cd /nfs/research/juan/singj/opendia/fasta-20261128
out=${TMPDIR:-/tmp}/fastafb_slim_$$; mkdir -p $out/runs_fb
for a in "$@"; do
  test -s runs_fb/$a.time || { echo "skip unfinished $a" >&2; continue; }
  f=runs_fb/$a/precursors.tsv
  sha256sum $f >> $out/full.sha256
  mkdir -p $out/runs_fb/$a
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} {print $c["precursor_id"]"\t"$c["proteins"]"\t"$c["discriminant_score"]"\t"$c["run_q_value"]"\t"$c["decoy"]"\t"$c["sequence"]"\t"$c["unmodified_sequence"]"\t"$c["charge"]"\t"$c["run"]}' $f | gzip -1 > $out/runs_fb/$a/precursors.slim.tsv.gz
  cp runs_fb/$a.log runs_fb/$a.time $out/runs_fb/
done
cp logs_fb/*.out $out/ 2>/dev/null || true
tar -C $out -cf - .
rm -rf $out
