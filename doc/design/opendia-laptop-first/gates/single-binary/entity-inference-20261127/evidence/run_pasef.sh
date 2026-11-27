#!/usr/bin/env bash
# diaPASEF acceptance: same .d and PQP through the pre-change and new OpenDIA binaries.
# Precursor outputs must be byte-identical. Entity inference is expected to differ.
set -euo pipefail
S=/home/singjc/.jcode/scratch/entity-inference-20261127
D=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published
IN=$D/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
LIB=$D/200623_hela_evosep_fractions_library_ptypic_decoy.pqp
NEW=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA
OLD=$S/OpenDIA.before
THREADS=${THREADS:-16}
for out in out_pasef_before out_pasef_after; do
  if test -e "$S/$out"; then echo "Refusing to reuse existing output $S/$out" >&2; exit 1; fi
done
sha256sum "$OLD" "$NEW" "$LIB" > "$S/pasef.sha256"
run() {
  local bin=$1 out=$2; shift 2
  mkdir -p "$S/$out"
  echo "JCODE_PROGRESS {\"message\":\"diaPASEF $out started\"}"
  /usr/bin/time -v -o "$S/$out.time.txt" "$bin" -in "$IN" -lib "$LIB" -out_dir "$S/$out" \
    -threads "$THREADS" -seed 0 -legacy_outputs "$@" > "$S/$out.log" 2>&1
  grep -E 'Elapsed|Maximum resident' "$S/$out.time.txt"
}
run "$OLD" out_pasef_before
run "$NEW" out_pasef_after -output_format both
for file in "$S"/out_pasef_before/legacy/* "$S"/out_pasef_before/*.calibration.json; do
  name=$(basename "$file")
  case "$name" in *.calibration.json) actual=$S/out_pasef_after/$name ;; *) actual=$S/out_pasef_after/legacy/$name ;; esac
  cmp "$file" "$actual"
  echo "BYTE_IDENTICAL $name"
done
if find "$S/out_pasef_after" -name '*.partial*' -o -type d -name 'opendia_cache_*' | grep -q .; then echo 'FAIL leftover partial/cache'; exit 1; fi
echo 'JCODE_PROGRESS {"message":"diaPASEF identity gate complete"}'
