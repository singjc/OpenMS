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
    @brief The OpenDIA pseudo-reverse decoy rule, in place.

    Reverses all residues but the last, moves residue modifications with their residues, then swaps a C-terminal K and
    R (any other last residue is replaced by one chosen by FNV-1a of the reversed sequence). @p sequence is the
    unmodified sequence; locations follow Modification (-1 N-terminus, size() C-terminus).
    @return false, with both arguments unchanged, if a modification is terminal or on the last residue (no decoy).
  */
  bool pseudoReverse(std::string& sequence, std::vector<Modification>& modifications);

  /// @p sequence with @p modifications in UniMod notation ("PEPC(UniMod:4)K"); modifications of one location in vector order.
  std::string uniModSequence(const std::string& sequence, const std::vector<Modification>& modifications);

  /**
    @brief Generate pseudo-reversed decoy compound records, retaining legacy annotation eligibility.

    @param[in] targets Target compounds, proteins, and transitions.
    @param[in] log_type Progress reporting mode.
    @return Decoy compounds and referenced proteins, with no transitions (prediction supplies them).
  */
  Assays decoySequences(const Assays& targets, ProgressLogger::LogType log_type = ProgressLogger::NONE);
}
