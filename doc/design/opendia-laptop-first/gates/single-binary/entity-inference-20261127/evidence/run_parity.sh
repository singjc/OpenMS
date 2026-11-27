#!/usr/bin/env bash
set -euo pipefail
BUILD=${BUILD:-/home/singjc/Documents/github/OpenMS-r1-full-ninja-build}
PYTHON=${PYTHON:-/home/singjc/Documents/github/pyprophet/.venv/bin/python}
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
OUT=${1:-"$HERE/output"}
if ninja -C "$BUILD" -t targets all | grep -q '^OpenDIAEntityInferenceParity_test:'; then
  cmake --build "$BUILD" --target OpenDIAEntityInferenceParity_test -j 3
  EXE=$(find "$BUILD" -type f -name OpenDIAEntityInferenceParity_test -perm -u+x -print -quit)
else
  cmake --build "$BUILD" --target OpenDIA_impl -j 3
  "$PYTHON" "$HERE/build_standalone.py" --build "$BUILD"
  EXE="$HERE/OpenDIAEntityInferenceParity_test"
fi
test -n "$EXE"
"$EXE" "$OUT"
"$PYTHON" "$HERE/parity.py" "$OUT" | tee "$OUT/parity.log"
