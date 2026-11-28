#!/usr/bin/env bash
# Rerun Astral and diaPASEF after the exact-lambda fix, into new folders. Precursor outputs must stay
# byte-identical. Entity statistics may change and are compared descriptively.
set -euo pipefail
S=/home/singjc/.jcode/scratch/entity-inference-20261127
P=/home/singjc/.jcode/scratch/perf-20261127
Q=/home/singjc/.jcode/scratch/qpx-output-20261127
D=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8
for out in out_pasef_lambda out_astral_lambda; do
  if test -e "$S/$out"; then echo "Refusing to reuse $S/$out" >&2; exit 1; fi
  mkdir -p "$S/$out"
done
echo 'JCODE_PROGRESS {"message":"diaPASEF and Astral started"}'
/usr/bin/time -v -o "$S/out_pasef_lambda.time.txt" "$B" -in "$D/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d" \
  -lib "$D/200623_hela_evosep_fractions_library_ptypic_decoy.pqp" -out_dir "$S/out_pasef_lambda" -threads 8 -seed 0 \
  -legacy_outputs -output_format both > "$S/out_pasef_lambda.log" 2>&1 &
pasef=$!
/usr/bin/time -v -o "$S/out_astral_lambda.time.txt" "$B" \
  -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw \
  -lib "$P/lib/library.pqp" -scoring_lib "$P/lib/library12_full.pqp" -out_dir "$S/out_astral_lambda" \
  -threads 12 -seed 0 -output_format both -legacy_outputs -peakmaps_out > "$S/out_astral_lambda.log" 2>&1 &
astral=$!
wait $pasef; echo 'JCODE_PROGRESS {"message":"diaPASEF done"}'
wait $astral; echo 'JCODE_PROGRESS {"message":"Astral done"}'
identical() { # reference_dir output_dir
  for f in "$1"/*; do
    n=$(basename "$f")
    case "$n" in *.calibration.json) a=$2/$n ;; *.precursors.tsv|precursors.tsv) a=$2/legacy/$n ;; *) continue ;; esac
    cmp "$f" "$a"; echo "BYTE_IDENTICAL $n"
  done
}
echo '== diaPASEF'; identical "$S/pasef_ref" "$S/out_pasef_lambda"
"$Q/check-env/bin/python" "$S/check_outputs.py" "$S/out_pasef_lambda" --refdir "$S/pasef_ref" > "$S/pasef_lambda_check.log" 2>&1
tail -1 "$S/pasef_lambda_check.log"
echo '== Astral'; identical "$P/out_base" "$S/out_astral_lambda"
"$Q/check-env/bin/python" "$S/check_outputs.py" "$S/out_astral_lambda" --refdir "$P/out_base" --require-xipm > "$S/astral_lambda_check.log" 2>&1
tail -1 "$S/astral_lambda_check.log"
/home/singjc/.jcode/scratch/qpx-conformance-env/bin/python "$Q/validate_opendia_qpx.py" "$S/out_astral_lambda/qpx" > "$S/astral_lambda_qpx.log" 2>&1
grep Overall "$S/astral_lambda_qpx.log"
/home/singjc/.jcode/scratch/qpx-conformance-env/bin/python "$Q/validate_opendia_qpx.py" "$S/out_pasef_lambda/qpx" > "$S/pasef_lambda_qpx.log" 2>&1
grep Overall "$S/pasef_lambda_qpx.log"
echo '== counts (previous build, then exact lambdas)'
/home/singjc/Documents/github/pyprophet/.venv/bin/python "$S/summarize.py" "$S/out_pasef_owned" "$S/out_pasef_lambda" "$S/out_astral_inference" "$S/out_astral_lambda"
