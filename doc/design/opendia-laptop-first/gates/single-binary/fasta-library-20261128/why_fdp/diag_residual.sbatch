#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 6:00:00
# PLAN_RESIDUAL_DIAG. Usage: sbatch --export=ALL,DRY=0|1 -J rd-<tag> -o why_fdp/rd-%x.%j.out diag_residual.sbatch OUT.json
set -euo pipefail
W=/nfs/research/juan/singj/opendia/fasta-20261128; cd $W
sv() { test -s $1/human_fdrbench.survey.tsv && echo $1/human_fdrbench.survey.tsv || echo -; }
ft() { f=$(ls $1/*.features.tsv 2>/dev/null | head -1); test -n "$f" && echo $f || echo -; }
spec() { echo "$1:$2/precursors.tsv:$3/human_fdrbench.library.search.parquet:$(sv $3):$(ft $2)"; }
DRY=${DRY:-0} $HOME/.local/bin/uv run -q --with pyarrow --with numpy python why_fdp/diag_residual.py $1 \
  $(spec FbF runs_fb/FbF cache_fb/FbF) $(spec FbV1 runs_fb/FbV1 cache_fb/FbV1) \
  $(spec FbV1_s1 runs_pf/FbV1_s1 cache_pf/FbV1_s1) $(spec FbV1_s2 runs_pf/FbV1_s2 cache_pf/FbV1_s2) \
  $(spec FbV1_g1 runs_pf/FbV1_g1 cache_pf/FbV1_g1) $(spec FbP runs_fb/FbP cache_fb/FbP) \
  $(spec FbP_g1 runs_pf/FbP_g1 cache_pf/FbP_g1) $(spec FbPp runs_fb/FbPp cache_fb/FbPp)
