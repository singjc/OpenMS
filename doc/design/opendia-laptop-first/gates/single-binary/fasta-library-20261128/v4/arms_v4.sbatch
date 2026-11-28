#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=240G
#SBATCH -t 24:00:00
# PLAN_SURVEY follow-up arms (V1+pool, code 143d3ac7c6, SIF v4). Usage: sbatch -J fastav4-<ARM> -o logs_v4/%x.%j.out arms_v4.sbatch ARM
#   Pp   P'  human, K* (default 800k), -survey_predicted -survey_pool, D seed 0, cold  (criteria 3, 4)
#   Pep  Pe' human + entrapment, same options, D seed 0, cold                         (criterion 2, both r0 and r_lib)
#   Hp   H'  human, same options, held-out run H, own cold cache                       (optional)
# References unchanged: F, Fe (runs_v2; survey options do not apply at K = 0), FH for H'.
set -uo pipefail
C=/nfs/research/juan/singj/opendia; W=$C/fasta-20261128; V=$C/entlib-hela-20261118/data
DD=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
HD=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
HUMAN=$W/fasta/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.fasta
ENTRAP=$W/fasta/human_entrapment.fasta
SIF=$W/opendia-fasta-v4-20261128.sif
ARM=$1
IN=$DD; FASTA=$HUMAN
case $ARM in
  Pp) ;;
  Pep) FASTA=$ENTRAP ;;
  Hp) IN=$HD ;;
  *) echo "unknown arm $ARM"; exit 2 ;;
esac
O=$W/runs_v4/$ARM; CACHE=$W/cache_v4/$ARM
test ! -e $O || { echo "exists $O"; exit 1; }
test ! -e $CACHE || { echo "cache exists $CACHE"; exit 1; }
mkdir -p $O $CACHE
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//') arm $ARM"
sha256sum $SIF $FASTA
SC=${TMPDIR:-/tmp}/fasta-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v -o $O.time singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $SIF OpenDIA -in $IN -fasta $FASTA -survey_predicted -survey_pool -library_cache $CACHE -out_dir $O -threads 16 -seed 0 \
  -legacy_outputs -trace_tensor_gb 200 -spectrum_cache_dir /tmp > $O.log 2>&1
rc=$?
rm -rf $SC
echo "$ARM exit=$rc $(grep -E 'Elapsed' $O.time | awk '{print $NF}') maxrss_kb=$(grep 'Maximum resident' $O.time | awk '{print $NF}')"
grep -E 'FASTA|Combined|Calibration: passes' $O.log | grep -v 'pairs predicted' | cut -c1-260
