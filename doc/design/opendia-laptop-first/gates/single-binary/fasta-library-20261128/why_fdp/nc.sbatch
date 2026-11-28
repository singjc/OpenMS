#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=48G
#SBATCH -t 6:00:00
# PLAN_NEIGHBOUR_COMPETITION. Usage: sbatch -J nc-<tag> -o why_fdp/nc-%x.%j.out nc.sbatch OUT.json [NC_DRY=1 via env]
set -euo pipefail
W=/nfs/research/juan/singj/opendia/fasta-20261128; cd $W
OUT=$1; shift
NC_DRY=${NC_DRY:-0} uv run -q --with pyarrow --with numpy python why_fdp/nc_remedy.py $OUT \
  fdrbench:FbF:runs_fb/FbF/precursors.tsv:cache_fb/FbF/human_fdrbench.library.search.parquet:0.9996135357569843 \
  fdrbench:FbV1:runs_fb/FbV1/precursors.tsv:cache_fb/FbV1/human_fdrbench.library.search.parquet:0.9994012382639613 \
  fdrbench:FbV1_s1:runs_pf/FbV1_s1/precursors.tsv:cache_pf/FbV1_s1/human_fdrbench.library.search.parquet:0.9994012382639613 \
  fdrbench:FbV1_s2:runs_pf/FbV1_s2/precursors.tsv:cache_pf/FbV1_s2/human_fdrbench.library.search.parquet:0.9994012382639613 \
  fdrbench:FbV1_g1:runs_pf/FbV1_g1/precursors.tsv:cache_pf/FbV1_g1/human_fdrbench.library.search.parquet:0.9994012382639613 \
  fdrbench:FbP:runs_fb/FbP/precursors.tsv:cache_fb/FbP/human_fdrbench.library.search.parquet:0.9996135357569843 \
  fdrbench:FbP_g1:runs_pf/FbP_g1/precursors.tsv:cache_pf/FbP_g1/human_fdrbench.library.search.parquet:0.9996135357569843 \
  fdrbench:FbPp:runs_fb/FbPp/precursors.tsv:cache_fb/FbPp/human_fdrbench.library.search.parquet:0.9994012382639613 \
  entrap:Fe:runs_v2/Fe/precursors.tsv:cache_v2/Fe/human_entrapment.library.search.parquet:0.5893838777524728 \
  entrap:Pe:runs_v2/Pe/precursors.tsv:cache_v2/Pe/human_entrapment.library.search.parquet:0.5893838777524728 \
  entrap:Pep:runs_v4/Pep/precursors.tsv:cache_v4/Pep/human_entrapment.library.search.parquet:0.5902413896340998
