# OpenDIA uses OpenMS readers and low-level classes only, no OpenSWATH code (design rule). Fails if the OpenDIA core
# library sources (include/, source/) reference OpenSWATH headers or the OpenSwath:: namespace. The command-line tool
# may use only TransitionPQPFile at the PQP conversion edge for -decoys predict.
file(GLOB_RECURSE core "${SOURCE_DIR}/include/*.h" "${SOURCE_DIR}/source/*.cpp")
if(NOT core)
  message(FATAL_ERROR "No OpenDIA core files found under ${SOURCE_DIR}")
endif()
set(bad "")
foreach(f IN LISTS core)
  file(STRINGS "${f}" hits REGEX "OPENSWATH|OpenSwath::|PeakMapExtractor|TransitionPQPFile|DIAHelper")
  if(hits)
    list(APPEND bad "${f}: ${hits}")
  endif()
endforeach()
# Keep the CLI format edge narrow as library-building functionality grows.
file(STRINGS "${SOURCE_DIR}/OpenDIA.cpp" edge_includes
  REGEX "^#include .*OpenMS/(ANALYSIS/OPENSWATH|OPENSWATHALGO)/")
foreach(include IN LISTS edge_includes)
  if(NOT include STREQUAL "#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>")
    list(APPEND bad "OpenDIA.cpp: ${include}")
  endif()
endforeach()
if(bad)
  string(REPLACE ";" "\n" bad "${bad}")
  message(FATAL_ERROR "OpenSWATH code in the OpenDIA core:\n${bad}")
endif()
message(STATUS "OpenDIA core: no OpenSWATH references in ${SOURCE_DIR}/include and /source")
