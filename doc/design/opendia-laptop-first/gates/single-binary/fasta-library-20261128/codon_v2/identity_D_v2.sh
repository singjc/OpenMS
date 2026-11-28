#!/usr/bin/env bash
# PLAN_FASTA criterion 1: the empirical-library path on D must give legacy outputs and calibration byte-identical to
# entity-inference-20261127/out_pasef_owned (read only), with the same command as run_pasef.sh.
set -euo pipefail
S=/home/singjc/.jcode/scratch/fasta-library-20261128
R=/home/singjc/.jcode/scratch/entity-inference-20261127/out_pasef_owned
D=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published
IN=$D/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
LIB=$D/200623_hela_evosep_fractions_library_ptypic_decoy.pqp
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA
OUT=$S/identity_D
test ! -e "$OUT"
mkdir -p "$OUT"
sha256sum "$BIN" > "$S/identity.sha256"
/usr/bin/time -v -o "$S/identity.time.txt" "$BIN" -in "$IN" -lib "$LIB" -out_dir "$OUT" -threads 16 -seed 0 -legacy_outputs \
  -output_format both > "$S/identity.log" 2>&1
for file in "$R"/legacy/* "$R"/*.calibration.json; do
  name=$(basename "$file")
  case "$name" in *.calibration.json) actual=$OUT/$name ;; *) actual=$OUT/legacy/$name ;; esac
  cmp "$file" "$actual" && echo "BYTE_IDENTICAL $name"
done
for f in "$R"/*.tsv "$R"/*.parquet; do cmp -s "$f" "$OUT/$(basename "$f")" && echo "SAME $(basename "$f")" || echo "DIFF $(basename "$f")"; done
