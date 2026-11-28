#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=96G
#SBATCH -t 3:00:00
cd /nfs/research/juan/singj/opendia/fasta-20261128
python3 why_fdp/h4_parent.py why_fdp/h4.json runs_fb:FbF runs_fb:FbP runs_fb:FbV1 runs_fb:FbPp runs_pf:FbV1_s1 runs_pf:FbV1_s2 runs_pf:FbV1_g1 runs_pf:FbP_g1
