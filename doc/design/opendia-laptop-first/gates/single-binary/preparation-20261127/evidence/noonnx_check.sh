#!/usr/bin/env bash
set -euo pipefail
source=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment
build=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build
scratch=/home/singjc/.jcode/scratch/decoyprep-20261127
# noonnx/OpenMS/config.h includes this build's generated OpenMS/config.h, then undefines WITH_ONNX.
# The overlay does not modify either source or generated configuration files.
g++ -std=c++23 -O1 -fopenmp -I"$scratch/noonnx" -I"$source/src/topp/OpenDIA/include" \
  -I"$source/src/openms/include" -I"$build/src/openms/include" \
  -isystem /home/singjc/Documents/github/openms-build/vcpkg_installed/x64-linux/include \
  "$scratch/noonnx_test.cpp" "$source/src/topp/OpenDIA/source/library/preparation/Preparation.cpp" \
  "$source/src/topp/OpenDIA/source/library/preparation/AssayPrediction.cpp" \
  "$source/src/topp/OpenDIA/source/library/preparation/DecoySequences.cpp" \
  -L"$build/lib" -Wl,-rpath,"$build/lib" -lOpenMS -o "$scratch/noonnx_test"
"$scratch/noonnx_test"
