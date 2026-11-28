#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=240G
#SBATCH -t 24:00:00
# PLAN_FASTA arms, v2 (AMEND_CALIBRATION, code 7fb7a6fb8d). Usage: sbatch -J fastav2-<ARM> -o logs_v2/%x.%j.out arms_v2.sbatch ARM
#   P   human, K* (default), D seed 0, cold library cache      (criterion 3, 4)
#   F   human, K = 0, D seed 0, cold library cache             (criterion 3, 4)
#   Pe  human + entrapment, K*, D seed 0, cold                 (criterion 2, 3)
#   Fe  human + entrapment, K = 0, D seed 0, cold              (criterion 2, 3)
#   P1, P2  human, K*, D seeds 1/2, reusing P's library (after P)
#   H   human, K*, held-out run H, own cold cache (the survey uses H, so its library differs from D's)
# -trace_tensor_gb 200 for every arm: a memory guard only (outputs are unchanged below it); the K = 0 arms may exceed 64.
set -uo pipefail
C=/nfs/research/juan/singj/opendia; W=$C/fasta-20261128; V=$C/entlib-hela-20261118/data
DD=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
HD=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
HUMAN=$W/fasta/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.fasta
ENTRAP=$W/fasta/human_entrapment.fasta
SIF=$W/opendia-fasta-v2-20261128.sif
BIN=OpenDIA
ARM=$1
IN=$DD; FASTA=$HUMAN; K=""; SEED=0; CACHE=$W/cache_v2/$ARM
case $ARM in
  P)  ;;
  F)  K="-prefilter_pairs 0" ;;
  Pe) FASTA=$ENTRAP ;;
  Fe) FASTA=$ENTRAP; K="-prefilter_pairs 0" ;;
  P1) SEED=1; CACHE=$W/cache_v2/P ;;
  P2) SEED=2; CACHE=$W/cache_v2/P ;;
  H)  IN=$HD ;;
  *) echo "unknown arm $ARM"; exit 2 ;;
esac
O=$W/runs_v2/$ARM
test ! -e $O || { echo "exists $O"; exit 1; }
case $ARM in P1|P2) test -e $W/cache_v2/P/$(basename ${FASTA%.fasta}).library.parquet.key || { echo "P library missing"; exit 1; } ;; *) test ! -e $CACHE || { echo "cache exists $CACHE"; exit 1; } ;; esac
mkdir -p $O $CACHE
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//') arm $ARM"
sha256sum $SIF $FASTA
SC=${TMPDIR:-/tmp}/fasta-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v -o $O.time singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $SIF $BIN -in $IN -fasta $FASTA $K -library_cache $CACHE -out_dir $O -threads 16 -seed $SEED -legacy_outputs \
  -trace_tensor_gb 200 -spectrum_cache_dir /tmp > $O.log 2>&1
rc=$?
rm -rf $SC
echo "$ARM exit=$rc $(grep -E 'Elapsed' $O.time | awk '{print $NF}') maxrss_kb=$(grep 'Maximum resident' $O.time | awk '{print $NF}')"
grep -E 'FASTA|Combined|Calibration: passes' $O.log | grep -v 'pairs predicted' | cut -c1-260
