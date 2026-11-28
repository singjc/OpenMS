#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
# PLAN_SURVEY: survey-only variants. Usage: sbatch -J survey-<V>-<RUN> -o logs_v3/%x.%j.out survey_v3.sbatch V RUN
#   V in V0 (current), V1 (predicted fragments), V2 (normalized co-elution), V3 (both); RUN in D, H.
# Every variant writes the survey features too (-survey_features) for +learn.
set -uo pipefail
C=/nfs/research/juan/singj/opendia; W=$C/fasta-20261128; V=$C/entlib-hela-20261118/data
DD=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
HD=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
HUMAN=$W/fasta/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.fasta
SIF=$W/opendia-fasta-v3-20261128.sif
VAR=$1; RUN=$2
case $VAR in
  V0) X="" ;;
  V1) X="-survey_predicted" ;;
  V2) X="-survey_normalized" ;;
  V3) X="-survey_predicted -survey_normalized" ;;
  *) echo "unknown variant $VAR"; exit 2 ;;
esac
case $RUN in D) IN=$DD ;; H) IN=$HD ;; *) echo "unknown run $RUN"; exit 2 ;; esac
O=$W/survey_v3/$VAR-$RUN
test ! -e $O || { echo "exists $O"; exit 1; }
mkdir -p $O
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//') variant $VAR run $RUN"
sha256sum $SIF $HUMAN
SC=${TMPDIR:-/tmp}/survey-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v -o $O.time singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $SIF OpenDIA -in $IN -fasta $HUMAN -library_cache $O -out_dir $O -threads 16 -survey_only -survey_features $X \
  -spectrum_cache_dir /tmp > $O.log 2>&1
rc=$?
rm -rf $SC
echo "$VAR-$RUN exit=$rc $(grep -E 'Elapsed' $O.time | awk '{print $NF}') maxrss_kb=$(grep 'Maximum resident' $O.time | awk '{print $NF}')"
grep -E 'FASTA survey|Calibration: passes|Search:' $O.log | grep -v 'predicted [0-9]' | cut -c1-200
gzip -1 $O/*.survey.tsv $O/*.survey.features.tsv
