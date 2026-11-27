#!/usr/bin/env bash
set -euo pipefail
S=/home/singjc/.jcode/scratch/entity-inference-20261127
P=/home/singjc/.jcode/scratch/perf-20261127
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA
OUT=$S/out_astral_inference
if test -e "$OUT"; then echo "Refusing to reuse existing output $OUT" >&2; exit 1; fi
mkdir -p "$OUT"
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8
sha256sum "$B" "$S/PLAN.md" > "$S/acceptance.sha256"
echo 'JCODE_PROGRESS {"message":"Astral QPX+both+legacy+XIPM started"}'
/usr/bin/time -v -o "$S/astral.time.txt" "$B" \
  -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw \
  -lib "$P/lib/library.pqp" -scoring_lib "$P/lib/library12_full.pqp" \
  -out_dir "$OUT" -threads 16 -seed 0 -output_format both -legacy_outputs -peakmaps_out \
  > "$S/astral.log" 2>&1
for file in "$P"/out_base/*; do
  name=$(basename "$file")
  case "$name" in
    *.precursors.tsv|precursors.tsv) actual=$OUT/legacy/$name ;;
    *.calibration.json) actual=$OUT/$name ;;
    *) continue ;;
  esac
  cmp "$file" "$actual"
  echo "BYTE_IDENTICAL $name"
done
if find "$OUT" -type d -name 'opendia_cache_*' | grep -q .; then echo 'FAIL cache directory remains'; exit 1; fi
if find "$OUT" -type f -name '*.partial*' | grep -q .; then echo 'FAIL partial file remains'; exit 1; fi
grep -E 'Elapsed|Maximum resident' "$S/astral.time.txt"
echo 'JCODE_PROGRESS {"message":"Astral identity gate complete"}'
