#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=240G
#SBATCH -t 24:00:00
# PLAN_ENTRAP_FDRBENCH arms on D, FASTA = FDRBench v1.1.1 protein-level shuffled entrapment of the human FASTA
# (human I->L + one shuffled entrapment protein per human protein, r = 1). SIF v4 (code 143d3ac7c6), seed 0, cold caches.
#   FbF    K = 0 (no prefilter)              <-> Fe  (custom entrapment)
#   FbP    K* default, V0 survey             <-> Pe
#   FbV1   K* default, -survey_predicted     (no custom counterpart; candidate default)
#   FbPp   K* default, -survey_predicted -survey_pool   <-> Pe'
set -uo pipefail
C=/nfs/research/juan/singj/opendia; W=$C/fasta-20261128; V=$C/entlib-hela-20261118/data
DD=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
FASTA=$W/fasta/human_fdrbench.fasta
SIF=$W/opendia-fasta-v4-20261128.sif
ARM=$1
case $ARM in
  FbF)  X="-prefilter_pairs 0" ;;
  FbP)  X="" ;;
  FbV1) X="-survey_predicted" ;;
  FbPp) X="-survey_predicted -survey_pool" ;;
  *) echo "unknown arm $ARM"; exit 2 ;;
esac
O=$W/runs_fb/$ARM; CACHE=$W/cache_fb/$ARM
test ! -e $O || { echo "exists $O"; exit 1; }
test ! -e $CACHE || { echo "cache exists $CACHE"; exit 1; }
mkdir -p $O $CACHE
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//') arm $ARM"
sha256sum $SIF $FASTA
SC=${TMPDIR:-/tmp}/fasta-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v -o $O.time singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $SIF OpenDIA -in $DD -fasta $FASTA $X -library_cache $CACHE -out_dir $O -threads 16 -seed 0 \
  -legacy_outputs -trace_tensor_gb 200 -spectrum_cache_dir /tmp > $O.log 2>&1
rc=$?
rm -rf $SC
echo "$ARM exit=$rc $(grep -E 'Elapsed' $O.time | awk '{print $NF}') maxrss_kb=$(grep 'Maximum resident' $O.time | awk '{print $NF}')"
grep -E 'FASTA|Combined|Calibration: passes' $O.log | grep -v 'pairs predicted' | cut -c1-260
