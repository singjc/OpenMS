#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=240G
#SBATCH -t 24:00:00
# PLAN_CLAIM_COMPETITION arms on D, SIF cc (OpenMS 86f3199513), cold caches. Usage:
#   sbatch -J fastacc-<ARM> -o logs_cc/%x.%j.out arms_cc.sbatch ARM
#   CcV1    FDRBench, -survey_predicted -claim_competition            <-> FbV1
#   CcV1_s1 as CcV1, seed 1                                           <-> FbV1_s1
#   CcP     FDRBench, -claim_competition (V0 survey)                  <-> FbP
#   CcPp    FDRBench, -survey_predicted -survey_pool -claim_competition <-> FbPp
#   CcF     FDRBench, -prefilter_pairs 0 -claim_competition           <-> FbF
#   Ce      custom entrapment, -prefilter_pairs 0 -claim_competition  <-> Fe
#   CcPep   custom entrapment, -survey_predicted -survey_pool -claim_competition <-> Pe'
# CcV1 also writes -features_out (gate G2).
set -uo pipefail
C=/nfs/research/juan/singj/opendia; W=$C/fasta-20261128; V=$C/entlib-hela-20261118/data
DD=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
FB=$W/fasta/human_fdrbench.fasta; EN=$W/fasta/human_entrapment.fasta
SIF=$W/opendia-fasta-cc-20261128.sif
ARM=$1
SEED=0; FASTA=$FB; FO=""
case $ARM in
  CcV1)    X="-survey_predicted"; FO="-features_out" ;;
  CcV1_s1) X="-survey_predicted"; SEED=1 ;;
  CcP)     X="" ;;
  CcPp)    X="-survey_predicted -survey_pool" ;;
  CcF)     X="-prefilter_pairs 0" ;;
  Ce)      X="-prefilter_pairs 0"; FASTA=$EN ;;
  CcPep)   X="-survey_predicted -survey_pool"; FASTA=$EN ;;
  *) echo "unknown arm $ARM"; exit 2 ;;
esac
O=$W/runs_cc/$ARM; CACHE=$W/cache_cc/$ARM
test ! -e $O || { echo "exists $O"; exit 1; }
test ! -e $CACHE || { echo "cache exists $CACHE"; exit 1; }
mkdir -p $O $CACHE
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//') arm $ARM"
sha256sum $SIF $FASTA
SC=${TMPDIR:-/tmp}/fastacc-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v -o $O.time singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $SIF OpenDIA -in $DD -fasta $FASTA $X -claim_competition -library_cache $CACHE -out_dir $O -threads 16 -seed $SEED $FO \
  -legacy_outputs -trace_tensor_gb 200 -spectrum_cache_dir /tmp > $O.log 2>&1
rc=$?
rm -rf $SC
echo "$ARM exit=$rc $(grep -E 'Elapsed' $O.time | awk '{print $NF}') maxrss_kb=$(grep 'Maximum resident' $O.time | awk '{print $NF}')"
grep -E 'FASTA|Combined|Calibration: passes|Claim competition' $O.log | grep -v 'pairs predicted' | cut -c1-260
