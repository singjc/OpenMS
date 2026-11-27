# OpenDIA does not use OpenSWATH extraction, workflow or OSW I/O code. The one core exception is the
# file-format-independent LevelContextInference statistics API, shared with the original OpenDIA and PyProphet.
# The CLI also permits TransitionPQPFile at the PQP conversion edge for -decoys predict.
get_filename_component(SOURCE_DIR "${SOURCE_DIR}" REALPATH)
file(GLOB_RECURSE core "${SOURCE_DIR}/include/*.h" "${SOURCE_DIR}/source/*.cpp")
if(NOT core)
  message(FATAL_ERROR "No OpenDIA core files found under ${SOURCE_DIR}")
endif()
set(bad "")
foreach(f IN LISTS core)
  file(STRINGS "${f}" hits REGEX "OPENSWATH|OpenSwath::|PeakMapExtractor|TransitionPQPFile|DIAHelper")
  foreach(hit IN LISTS hits)
    if(f STREQUAL "${SOURCE_DIR}/source/results/Quantification.cpp" AND
       hit STREQUAL "#include <OpenMS/ANALYSIS/OPENSWATH/LevelContextInference.h>")
      continue()
    endif()
    list(APPEND bad "${f}: ${hit}")
  endforeach()
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
message(STATUS "OpenDIA: no OpenSWATH workflow/extraction dependencies; only shared LevelContextInference statistics and CLI PQP conversion are permitted")
