#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=48G
#SBATCH -t 6:00:00
#SBATCH --job-name=relabel-eval
# Relabel re-derivation (relabel-20261127). Read-only on every run, report and library; writes only under $W.
# Each record is evaluated twice with the same evaluator logic:
#   LABELS=old        the frozen rule, a reproduction check against the published EVAL_*.txt;
#   LABELS=corrected  instr('_p_target') plus the step A list of exact human substrings counted as human.
# Usage: sbatch -o $W/logs/%x.%j.out relabel.sbatch
set -uo pipefail
C=/nfs/research/juan/singj/opendia
W=$C/relabel-20261127
SIF=$W/relabel-eval-20261127.sif
LIB=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp   # predicted-decoy library of every r~7 OpenDIA run
TGT=$C/entlib-r5-20261118/targets.pqp                                     # targets-only library DIA-NN searched
cd "$W"
sha256sum -c bundle.sha256 || exit 3
mkdir -p lay out
py() { singularity exec --cleanenv --no-home --bind "$C" "$SIF" python "$@"; }
# Layouts in the evaluator's <P>_<arm>_s<seed> form, symlinks only (no copy, no write to run directories).
# record:runs_dir:set:arms (arm names as published). Only run tables are linked.
records=(
  "odp:entlib-odp5-20261118:D:cur,new"
  "odp:entlib-odp5-20261118:H:cur,new"
  "sib:entlib-sib-20261121:D:cur,sib"
  "sib:entlib-sib-20261121:H:cur,sib"
  "rtr:entlib-rtr-20261121:D:cur,rtr"
  "rtr:entlib-rtr-20261121:H:cur,rtr"
  "gbt:entlib-gbt-20261122:D:cur,gbt"
  "gbt:entlib-gbt-20261122:H:cur,gbt"
  "v7m6:entlib-v7m6-20261123:D:cur,v7m6"
  "blend:blend-20261125:D:B5"
  "blend:blend-20261125:H:B5"
)
for r in "${records[@]}"; do
  IFS=: read -r name dir set arms <<< "$r"
  for arm in ${arms//,/ }; do
    for s in 0 1 2 3 4; do
      src=$C/$dir/runs/${set}_${arm}_s$s
      dst=lay/$name/${set}_${arm}_s$s
      [ -d "$src" ] || { echo "MISSING $src"; continue; }
      mkdir -p "$dst"
      for f in "$src"/*.precursors.tsv*; do
        case "$(basename "$f")" in precursors.tsv|precursors.tsv.gz) ;; *) ln -sfn "$f" "$dst/";; esac
      done
    done
  done
done
for name in odp sib rtr gbt v7m6 blend; do
  for set in D H; do
    ls -d lay/$name/${set}_* > /dev/null 2>&1 || continue
    for labels in old corrected; do
      LABELS=$labels RELABEL=$W/relabel_hela.tsv py eval_pdec_relabel.py "$LIB" "lay/$name" "$set" > out/${name}_${set}_${labels}.txt 2>&1
      echo "$name $set $labels exit=$?"
    done
  done
done
for set in D H; do
  for labels in old corrected; do
    LABELS=$labels RELABEL=$W/relabel_hela.tsv py eval_diann_relabel.py "$TGT" "$C/entlib-diann5-20261119/$set/report.parquet" > out/diann_${set}_${labels}.txt 2>&1
    echo "diann $set $labels exit=$?"
  done
done
sha256sum out/*.txt > out/outputs.sha256
echo done
