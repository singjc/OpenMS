#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 2:00:00
cd /nfs/research/juan/singj/opendia/fasta-20261128
R=20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737
for a in FbV1_g1 FbP_g1 FbV1_s1 FbV1_s2; do
  r=0.9994; [ $a = FbP_g1 ] && r=0.9996
  s=${TMPDIR:-/tmp}/slim_$a.tsv.gz
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} {print $c["precursor_id"]"\t"$c["proteins"]"\t"$c["decoy"]}' runs_pf/$a/precursors.tsv | gzip -1 > $s
  echo "== $a"; ls -la runs_pf/$a/$R.features.tsv | awk '{print $5}'
  uv run -q --with numpy python why_fdp/h3_network.py runs_pf/$a/$R.features.tsv $s $r
  rm -f $s
done
