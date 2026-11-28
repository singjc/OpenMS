// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
#pragma once

#include <library/preparation/FastaLibrary.h>
#include <workflow/Workflow.h>

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /// Settings of a FASTA library (PLAN_FASTA).
  struct FastaLibraryOptions
  {
    std::string fasta;
    std::string cache_dir;                         ///< Where the cached libraries go.
    Preparation::CandidateOptions candidates;
    Preparation::PredictionSettings prediction;
    Preparation::FragmentOptions fragments;
    std::size_t prefilter_pairs = 0;               ///< K: pairs kept by the survey (0 = no prefilter).
    std::size_t panel_pairs = 40000;               ///< Calibration panel pairs of the survey library.
    bool survey_only = false;                      ///< Stop after the survey statistics (design stage).
  };

  struct FastaLibraryPaths
  {
    std::string candidates;   ///< Survey library: every pair, theoretical fragments (panel: predicted).
    std::string survey;       ///< Survey statistics TSV (empty without prefilter).
    std::string library;      ///< Final library (N fragments), empty with survey_only.
    std::string search;       ///< Final search library (top 6), empty with survey_only.
  };

  /**
    @brief Build (or reuse from FastaLibraryOptions::cache_dir) the fixed library of a FASTA for @p runs.

    1. <stem>.candidates.parquet: candidates with PeptDeep RT (and 1/K0) and survey fragments.
    2. With K > 0, <stem>.survey.tsv: per pair, the largest search co-elution of target and decoy over all runs
       (RunAnalysis::survey with @p run_options, calibration restricted to the panel pairs).
    3. <stem>.library.parquet and <stem>.library.search.parquet: predicted MS2 for the selected pairs.
    Each file has a .key next to it (inputs, settings, model and OpenMS versions); a file whose key matches is reused.
  */
  FastaLibraryPaths buildFastaLibrary(const FastaLibraryOptions& options, const std::vector<std::string>& runs,
                                      const RunOptions& run_options, const std::function<void(const std::string&)>& log = {});
}
