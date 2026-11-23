# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
endif()
project(OpenDIA)

if(NOT TARGET OpenMS::OpenMS_CLI)
  find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)
endif()

add_library(OpenDIA_impl STATIC
  source/OpenDIASearch.cpp
  source/OpenDIATraces.cpp
  source/OpenDIARTMap.cpp
  source/OpenDIACalibration.cpp
  source/OpenDIAStage2.cpp
  source/OpenDIAPeakPicker.cpp
  source/OpenDIAInterference.cpp
  source/OpenDIAScorer.cpp
  source/OpenDIAGBM.cpp
  source/OpenDIAFragmentGrid.cpp
  source/OpenDIAWorkflow.cpp
  source/OpenDIARTRefine.cpp
  source/OpenDIALibrary.cpp
  source/OpenDIALibraryParquet.cpp
  include/OpenDIALibrary.h
  include/OpenDIAMzIndex.h
  include/OpenDIASearch.h
  include/OpenDIATraces.h
  include/OpenDIARTMap.h
  include/OpenDIACalibration.h
  include/OpenDIAStage2.h
  include/OpenDIAPeakPicker.h
  include/OpenDIAInterference.h
  include/OpenDIAScorer.h
  include/OpenDIAGBM.h
  include/OpenDIAFragmentGrid.h
  include/OpenDIARTRefine.h
  include/OpenDIAWorkflow.h)
target_include_directories(OpenDIA_impl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(OpenDIA_impl PUBLIC OpenMS::OpenMS_CLI)
# OpenDIALibrary reads PQP files with the SQLite C API (the sqlite3 OpenMS itself is built with).
if(OPENMS_USE_VCPKG AND TARGET unofficial::sqlite3::sqlite3)
  target_link_libraries(OpenDIA_impl PRIVATE unofficial::sqlite3::sqlite3)
elseif(TARGET SQLite::SQLite3)
  target_link_libraries(OpenDIA_impl PRIVATE SQLite::SQLite3)
else()
  find_package(SQLite3 REQUIRED)
  target_link_libraries(OpenDIA_impl PRIVATE SQLite::SQLite3)
endif()
# OpenDIACalibration reads its JSON with nlohmann::json, a PRIVATE dependency of libOpenMS.
if(NOT TARGET nlohmann_json::nlohmann_json)
  find_package(nlohmann_json CONFIG REQUIRED)
endif()
target_link_libraries(OpenDIA_impl PRIVATE nlohmann_json::nlohmann_json)
# OpenDIALibraryParquet reads and writes Parquet with the Arrow/Parquet libOpenMS is built with. Inside the OpenMS
# build, Arrow is part of libOpenMS (statically linked and exported), so only its headers are added: linking the
# static archives again would put a second Arrow, with its own registries, into OpenDIA.
if(DEFINED OPENMS_ARROW_TARGET AND TARGET ${OPENMS_ARROW_TARGET} AND DEFINED OPENMS_PARQUET_TARGET AND TARGET ${OPENMS_PARQUET_TARGET})
  foreach(_opendia_arrow IN ITEMS ${OPENMS_ARROW_TARGET} ${OPENMS_PARQUET_TARGET})
    if(_opendia_arrow MATCHES "_static$")
      target_include_directories(OpenDIA_impl SYSTEM PRIVATE $<TARGET_PROPERTY:${_opendia_arrow},INTERFACE_INCLUDE_DIRECTORIES>)
      target_compile_definitions(OpenDIA_impl PRIVATE $<TARGET_PROPERTY:${_opendia_arrow},INTERFACE_COMPILE_DEFINITIONS>)
    else()
      target_link_libraries(OpenDIA_impl PRIVATE ${_opendia_arrow})
    endif()
  endforeach()
else()
  find_package(Arrow CONFIG REQUIRED)
  find_package(Parquet CONFIG REQUIRED)
  target_link_libraries(OpenDIA_impl PRIVATE Arrow::arrow_shared Parquet::parquet_shared)
endif()
if(COMMAND openms_add_library_compiler_flags)
  openms_add_library_compiler_flags(OpenDIA_impl)
endif()

add_executable(OpenDIA OpenDIA.cpp)
target_link_libraries(OpenDIA PRIVATE OpenDIA_impl)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(OpenDIA)
endif()

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(CTest)
  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()
endif()
