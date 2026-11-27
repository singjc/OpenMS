#!/usr/bin/env bash
set -euo pipefail
root=/home/singjc/.jcode/scratch/decoyprep-20261127/guard
script=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/src/topp/OpenDIA/test/no_openswath.cmake
mkdir -p "$root/include/library/preparation" "$root/source/library/preparation"
printf '#pragma once\n' > "$root/include/library/preparation/Assay.h"
printf '#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>\n' > "$root/OpenDIA.cpp"
cmake -DSOURCE_DIR="$root" -P "$script"
printf '#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryPreparation.h>\n' >> "$root/OpenDIA.cpp"
if cmake -DSOURCE_DIR="$root" -P "$script"; then echo 'FAIL: forbidden CLI include accepted'; exit 1; fi
printf '#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>\n' > "$root/OpenDIA.cpp"
printf '#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>\n' > "$root/source/library/preparation/Bad.cpp"
if cmake -DSOURCE_DIR="$root" -P "$script"; then echo 'FAIL: nested forbidden core include accepted'; exit 1; fi
echo 'PASS: CLI boundary and nested core negative fixtures rejected'
