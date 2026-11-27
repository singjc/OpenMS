// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
#pragma once

#include <library/preparation/Assay.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS::OpenDIA::Preparation
{
  struct PredictionOptions
  {
    std::string model_directory;
    int threads{4};
    ProgressLogger::LogType log_type{ProgressLogger::NONE};
  };

  /**
    @brief Add pseudo-reversed, independently predicted decoys to a target library.

    Estimates instrument and collision energy from target intensities and maps predicted RT and mobility
    to the target library's scale. Retains the historical annotation exclusions, stable fragment ordering,
    short-assay policy and 80% compound/protein coverage requirement. No file I/O or analysis state is owned here.
    Requires WITH_ONNX. Targets are passed by value so callers can transfer their storage with std::move.
  */
  Assays predictDecoys(Assays targets, const PredictionOptions& options = {});
}
