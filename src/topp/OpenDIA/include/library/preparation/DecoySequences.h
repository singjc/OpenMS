// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
#pragma once

#include <library/preparation/Assay.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS::OpenDIA::Preparation
{
  /**
    @brief Generate pseudo-reversed decoy compound records, retaining legacy annotation eligibility.

    @param[in] targets Target compounds, proteins, and transitions.
    @param[in] log_type Progress reporting mode.
    @return Decoy compounds and referenced proteins, with no transitions (prediction supplies them).
  */
  Assays decoySequences(const Assays& targets, ProgressLogger::LogType log_type = ProgressLogger::NONE);
}
