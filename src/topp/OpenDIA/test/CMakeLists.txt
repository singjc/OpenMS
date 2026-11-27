# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------

# Output contracts are tested through both native readers and direct Parquet inspection.
foreach(_output_test IN ITEMS Quantification EntityInferenceParity QpxExport ReportExport PeakMapExport Annotations)
  add_executable(OpenDIA${_output_test}_test ${_output_test}_test.cpp)
  target_link_libraries(OpenDIA${_output_test}_test PRIVATE OpenDIA_impl)
  add_test(NAME OpenDIA${_output_test}_test COMMAND OpenDIA${_output_test}_test ${CMAKE_CURRENT_BINARY_DIR}/OpenDIA${_output_test}_test_data)
  # The implementation uses headers privately. Reader-based output tests need the same Arrow headers.
  foreach(_output_arrow IN ITEMS ${OPENMS_ARROW_TARGET} ${OPENMS_PARQUET_TARGET} Arrow::arrow_shared Parquet::parquet_shared)
    if(TARGET ${_output_arrow})
      target_include_directories(OpenDIA${_output_test}_test SYSTEM PRIVATE $<TARGET_PROPERTY:${_output_arrow},INTERFACE_INCLUDE_DIRECTORIES>)
    endif()
  endforeach()
endforeach()

add_executable(OpenDIAPreparation_test Preparation_test.cpp)
target_link_libraries(OpenDIAPreparation_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAPreparation_test COMMAND OpenDIAPreparation_test ${CMAKE_CURRENT_BINARY_DIR}/preparation_sequences)
if(WITH_ONNX)
  add_test(NAME OpenDIAPreparation_prediction_test COMMAND OpenDIAPreparation_test ${CMAKE_CURRENT_BINARY_DIR}/preparation_prediction --prediction)
endif()

add_executable(OpenDIAWorkflow_test Workflow_test.cpp)
target_link_libraries(OpenDIAWorkflow_test PRIVATE OpenDIA_impl)
if(EXISTS "${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp")
  add_test(NAME OpenDIAWorkflow_test COMMAND OpenDIAWorkflow_test ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp)
else()
  add_test(NAME OpenDIAWorkflow_test COMMAND OpenDIAWorkflow_test ${CMAKE_CURRENT_BINARY_DIR})
endif()

add_executable(OpenDIASearch_test Search_test.cpp)
target_link_libraries(OpenDIASearch_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIASearch_test COMMAND OpenDIASearch_test)

add_executable(OpenDIATraces_test Traces_test.cpp)
target_link_libraries(OpenDIATraces_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIATraces_test COMMAND OpenDIATraces_test ${CMAKE_CURRENT_BINARY_DIR}/OpenDIATraces_test)

add_executable(OpenDIACalibration_test Calibration_test.cpp)
target_link_libraries(OpenDIACalibration_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIACalibration_test COMMAND OpenDIACalibration_test)

add_executable(OpenDIAPeakFeatures_test PeakFeatures_test.cpp)
target_link_libraries(OpenDIAPeakFeatures_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAPeakFeatures_test COMMAND OpenDIAPeakFeatures_test)

add_executable(OpenDIASpectrumCache_test SpectrumCache_test.cpp)
target_link_libraries(OpenDIASpectrumCache_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIASpectrumCache_test COMMAND OpenDIASpectrumCache_test)

add_executable(OpenDIAPeakPicker_test PeakPicker_test.cpp)
target_link_libraries(OpenDIAPeakPicker_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAPeakPicker_test COMMAND OpenDIAPeakPicker_test)

add_executable(OpenDIAInterference_test Interference_test.cpp)
target_link_libraries(OpenDIAInterference_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAInterference_test COMMAND OpenDIAInterference_test)

add_executable(OpenDIAScorer_test Scorer_test.cpp)
target_link_libraries(OpenDIAScorer_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAScorer_test COMMAND OpenDIAScorer_test)

add_executable(OpenDIAGBM_test GBM_test.cpp)
target_link_libraries(OpenDIAGBM_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAGBM_test COMMAND OpenDIAGBM_test)

add_executable(OpenDIAFragmentGrid_test FragmentGrid_test.cpp)
target_link_libraries(OpenDIAFragmentGrid_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAFragmentGrid_test COMMAND OpenDIAFragmentGrid_test)
add_executable(OpenDIAMzIndex_test MzIndex_test.cpp)
target_link_libraries(OpenDIAMzIndex_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAMzIndex_test COMMAND OpenDIAMzIndex_test)
add_executable(OpenDIARTRefine_test RTRefine_test.cpp)
target_link_libraries(OpenDIARTRefine_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIARTRefine_test COMMAND OpenDIARTRefine_test)
if(EXISTS "${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp")
  # Copies the PQP fixture into the build tree first, so no Parquet copy is ever written into src/tests.
  add_executable(OpenDIALibraryParquet_test LibraryParquet_test.cpp)
  target_link_libraries(OpenDIALibraryParquet_test PRIVATE OpenDIA_impl)
  add_test(NAME OpenDIALibraryParquet_test COMMAND OpenDIALibraryParquet_test ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp)
endif()
# Design rule: no OpenSWATH code in the OpenDIA core library.
add_test(NAME OpenDIA_no_openswath
  COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/.. -P ${CMAKE_CURRENT_SOURCE_DIR}/no_openswath.cmake)

if(EXISTS "${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML" AND
   EXISTS "${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp")
  # OpenDIA fails closed: a library too small to calibrate, and two inputs with one file name (outputs
  # would collide).
  add_test(NAME OpenDIA_fails_closed_on_tiny_library
    COMMAND $<TARGET_FILE:OpenDIA>
      -in ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML
      -lib ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
      -out_dir ${CMAKE_CURRENT_BINARY_DIR}/OpenDIA_tiny)
  set_tests_properties(OpenDIA_fails_closed_on_tiny_library PROPERTIES WILL_FAIL TRUE)
  add_test(NAME OpenDIA_rejects_duplicate_run_names
    COMMAND $<TARGET_FILE:OpenDIA>
      -in ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML
      -lib ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
      -out_dir ${CMAKE_CURRENT_BINARY_DIR}/OpenDIA_duplicate)
  set_tests_properties(OpenDIA_rejects_duplicate_run_names PROPERTIES WILL_FAIL TRUE)
endif()
