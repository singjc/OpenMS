#!/usr/bin/env bash
# PLAN_FASTA design stage: the survey alone on HeLa D with the human FASTA (defaults), then recall(K) against the
# empirical IDs. The empirical outputs are read only.
set -euo pipefail
S=/home/singjc/.jcode/scratch/fasta-library-20261128
D=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published
IN=$D/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
FASTA=$D/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.fasta
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA
OUT=$S/design_D
THREADS=${THREADS:-16}
if test -e "$OUT"; then echo "Refusing to reuse $OUT" >&2; exit 1; fi
mkdir -p "$OUT"
sha256sum "$BIN" "$FASTA" > "$S/design.sha256"
echo 'JCODE_PROGRESS {"message":"design survey started"}'
/usr/bin/time -v -o "$S/design.time.txt" "$BIN" -in "$IN" -fasta "$FASTA" -out_dir "$OUT" -threads "$THREADS" -survey_only \
  > "$S/design.log" 2>&1
grep -E 'Elapsed|Maximum resident' "$S/design.time.txt"
echo 'JCODE_PROGRESS {"message":"design survey finished"}'
