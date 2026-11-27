#!/usr/bin/env bash
# Perf round 2026-11-27: one OpenDIA arm on Astral raw01 seed0 (laptop profile, 16 threads), with GNU time and perf stat.
# Usage: run_arm.sh <arm> <binary>
set -euo pipefail
P=/home/singjc/.jcode/scratch/perf-20261127
arm=$1; bin=$2
out=$P/out_$arm
test ! -e "$out"
mkdir -p "$out"
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8
raw=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
sha256sum "$bin" > "$P/$arm.binary.sha256"
echo "JCODE_PROGRESS {\"message\":\"$arm started\"}"
perf stat -x, -o "$P/$arm.perfstat.csv" -e task-clock,cycles,instructions,cache-references,cache-misses,LLC-loads,LLC-load-misses,page-faults \
  /usr/bin/time -v -o "$P/$arm.time.txt" "$bin" -in "$raw" -lib $P/lib/library.pqp -scoring_lib $P/lib/library12_full.pqp \
  -out_dir "$out" -threads 16 -seed 0 > "$P/$arm.log" 2>&1
echo "JCODE_PROGRESS {\"message\":\"$arm done\"}"
grep -E "Elapsed|Maximum resident" "$P/$arm.time.txt"
