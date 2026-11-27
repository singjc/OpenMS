#!/usr/bin/env bash
set -euo pipefail
scratch=/home/singjc/.jcode/scratch/decoyprep-20261127
source=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment
build=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build
mkdir -p "$scratch/cli_ref" "$scratch/cli_owned"
cp "$build/src/tests/class_tests/openms/topp_OpenDIA_tests/preparation_prediction/preparation_targets.pqp" "$scratch/cli_targets.pqp"
set +e
"$scratch/OpenDIA.ref" -in "$source/src/tests/topp/OpenSwathWorkflow_17_input.mzML" -lib "$scratch/cli_targets.pqp" -out_dir "$scratch/cli_ref" -decoys predict -threads 4 > "$scratch/cli_ref.log" 2>&1
old_status=$?
"$build/bin/OpenDIA" -in "$source/src/tests/topp/OpenSwathWorkflow_17_input.mzML" -lib "$scratch/cli_targets.pqp" -out_dir "$scratch/cli_owned" -decoys predict -threads 4 > "$scratch/cli_owned.log" 2>&1
new_status=$?
set -e
echo "CLI exit codes: reference=$old_status owned=$new_status"
test "$old_status" = "$new_status"
cmp "$scratch/cli_ref/cli_targets.predicted_decoys.pqp" "$scratch/cli_owned/cli_targets.predicted_decoys.pqp"
echo 'PASS: fresh CLI predicted-decoy PQP bytes identical'
sha256sum "$scratch/cli_ref/cli_targets.predicted_decoys.pqp" "$scratch/cli_owned/cli_targets.predicted_decoys.pqp"
if find "$scratch/cli_ref" "$scratch/cli_owned" -name '*.tmp.pqp' | grep .; then exit 1; fi
echo 'PASS: prediction temporary files removed'
