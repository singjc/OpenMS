#!/usr/bin/env bash
set -euo pipefail
S=/home/singjc/.jcode/scratch/qpx-output-20261127
OUT=$S/out_astral_qpx
/home/singjc/.jcode/scratch/qpx-conformance-env/bin/python "$S/validate_opendia_qpx.py" "$OUT/qpx" --json-output "$S/astral-conformance.json" > "$S/astral-conformance.log" 2>&1
cat "$S/astral-conformance.log"
echo 'JCODE_PROGRESS {"message":"Actual upstream QPX and every producer ID passed"}'
"$S/check-env/bin/python" "$S/check_outputs.py" "$OUT" --refdir /home/singjc/.jcode/scratch/perf-20261127/out_base > "$S/astral-independent.log" 2>&1
cat "$S/astral-independent.log"
echo 'JCODE_PROGRESS {"message":"Independent reports, matrices, precision, legacy and XIPM links passed"}'
"$S/read_native_xipm" "$OUT/LFQ_Astral_DIA_15min_50ng_Human_01.xipm" > "$S/native-xipm.log" 2>&1
cat "$S/native-xipm.log"
