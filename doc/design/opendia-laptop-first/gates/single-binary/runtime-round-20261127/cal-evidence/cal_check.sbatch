#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 4:00:00
# Streaming windowed calibration check (OpenMS a251c23b40). One node: prev (e3533d93, ff4fdc067f) then ifg, then prev, ifg,
# D or H seed 0, -threads 8. Gate: every ifg output cmp-identical to runs/<set>_base_r1 of the frozen part A/B.
# Usage: sbatch -J ifg-D -o $W/logs/%x.%j.out ifg_ab.sbatch D|H
set -uo pipefail
C=/nfs/research/juan/singj/opendia; W=$C/part-20261127; V=$C/entlib-hela-20261118/data
DD=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
HD=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
LIB=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp
SIF=$W/opendia-cal-20261127.sif
P=${1:-D}
case $P in D) IN=$DD;; H) IN=$HD;; *) echo "unknown set $P"; exit 2;; esac
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
echo "node $(hostname) set $P"; sha256sum $SIF
run() {
  local arm=$1 exe=$2 rep=$3
  local O=$W/runs/${P}_${arm}_r$rep
  test ! -e $O || { echo "exists $O"; return 1; }
  mkdir -p $O
  local SC=${TMPDIR:-/tmp}/cal-$SLURM_JOB_ID-$arm-$rep; mkdir -p $SC
  /usr/bin/time -v -o $O.time singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
    $SIF $BIN/$exe -in $IN -lib $LIB -out_dir $O -threads 8 -seed 0 > $O.log 2>&1
  local rc=$?
  rm -rf $SC
  echo "$arm r$rep exit=$rc $(grep -E 'Elapsed' $O.time | awk '{print $NF}') maxrss_kb=$(grep 'Maximum resident' $O.time | awk '{print $NF}')"
}
run cal OpenDIA.cal 1
REF=$W/runs/${P}_base_r1
for O in $W/runs/${P}_cal_r1; do
  ls -d $O/opendia_cache_* 2>/dev/null && echo "GATE CACHE_LEFT $O"
  for f in $(cd $REF && ls); do cmp -s $O/$f $REF/$f && echo "GATE SAME $(basename $O) $f" || echo "GATE DIFF $(basename $O) $f"; done
  [ "$(cd $REF && ls)" = "$(cd $O && ls)" ] && echo "GATE FILESET_SAME $(basename $O)" || echo "GATE FILESET_DIFF $(basename $O)"
done
