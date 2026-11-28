#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=240G
#SBATCH -t 24:00:00
# PLAN_PREFILTER_FDP causal arms on D, FDRBench database, SIF v4. Usage: sbatch -J fastapf-<ARM> -o logs_pf/%x.%j.out arms_pf.sbatch ARM
set -uo pipefail
C=/nfs/research/juan/singj/opendia; W=$C/fasta-20261128; V=$C/entlib-hela-20261118/data
DD=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
FASTA=$W/fasta/human_fdrbench.fasta
SIF=$W/opendia-fasta-v4-20261128.sif
ARM=$1
SEED=0
case $ARM in
  FbV1_g1) X="-survey_predicted -gbm_rounds 1" ;;
  FbP_g1)  X="-gbm_rounds 1" ;;
  FbV1_s1) X="-survey_predicted"; SEED=1 ;;
  FbV1_s2) X="-survey_predicted"; SEED=2 ;;
  *) echo "unknown arm $ARM"; exit 2 ;;
esac
O=$W/runs_pf/$ARM; CACHE=$W/cache_pf/$ARM
test ! -e $O || { echo "exists $O"; exit 1; }
test ! -e $CACHE || { echo "cache exists $CACHE"; exit 1; }
mkdir -p $O $CACHE
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//') arm $ARM"
sha256sum $SIF $FASTA
SC=${TMPDIR:-/tmp}/fasta-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v -o $O.time singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $SIF OpenDIA -in $DD -fasta $FASTA $X -library_cache $CACHE -out_dir $O -threads 16 -seed $SEED -features_out \
  -legacy_outputs -trace_tensor_gb 200 -spectrum_cache_dir /tmp > $O.log 2>&1
rc=$?
rm -rf $SC
echo "$ARM exit=$rc $(grep -E 'Elapsed' $O.time | awk '{print $NF}') maxrss_kb=$(grep 'Maximum resident' $O.time | awk '{print $NF}')"
grep -E 'FASTA|Combined|Calibration: passes' $O.log | grep -v 'pairs predicted' | cut -c1-260
