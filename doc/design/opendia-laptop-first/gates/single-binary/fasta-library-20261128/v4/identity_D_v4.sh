#!/usr/bin/env bash
# PLAN_FASTA criterion 1 recheck at the PLAN_SURVEY +pool code commit (143d3ac7c6): empirical-library path on D vs
# entity-inference-20261127/out_pasef_owned (read only). Same command as identity_D.sh, new output paths.
set -euo pipefail
S=/home/singjc/.jcode/scratch/fasta-library-20261128
R=/home/singjc/.jcode/scratch/entity-inference-20261127/out_pasef_owned
D=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published
IN=$D/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
LIB=$D/200623_hela_evosep_fractions_library_ptypic_decoy.pqp
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA
OUT=$S/identity_D_v4
test ! -e "$OUT"
mkdir -p "$OUT"
sha256sum "$BIN" > "$S/identity_v4.sha256"
/usr/bin/time -v -o "$S/identity_v4.time.txt" "$BIN" -in "$IN" -lib "$LIB" -out_dir "$OUT" -threads 16 -seed 0 -legacy_outputs \
  -output_format both > "$S/identity_v4.log" 2>&1
for file in "$R"/legacy/* "$R"/*.calibration.json; do
  name=$(basename "$file")
  case "$name" in *.calibration.json) actual=$OUT/$name ;; *) actual=$OUT/legacy/$name ;; esac
  cmp "$file" "$actual" && echo "BYTE_IDENTICAL $name"
done
# v1 identity run (same command, a445c1892e) against v2: every output file must be byte-identical.
for file in "$S"/identity_D/*.tsv "$S"/identity_D/*.parquet "$S"/identity_D/legacy/* "$S"/identity_D/*.calibration.json; do
  rel=${file#"$S"/identity_D/}
  cmp -s "$file" "$OUT/$rel" && echo "V1_V2_SAME $rel" || echo "V1_V2_DIFF $rel"
done
