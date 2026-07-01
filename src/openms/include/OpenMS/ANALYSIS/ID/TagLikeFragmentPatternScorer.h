// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  /**
    @brief Score peptide fragment-pattern matches with DirectTag-like components.

    This scorer operates on FragmentIndex::SpectrumMatch records and combines
    matched-ion coverage, matched-intensity fraction, fragment m/z fidelity,
    complementary b/y cleavage support, and a simple random-match proxy into a
    single spectrum-level score.

    The implementation is intended to be reusable by staged DIA prefilters and
    other lightweight spectrum-to-peptide candidate filters without tying the
    logic to a specific workflow.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI TagLikeFragmentPatternScorer :
    public DefaultParamHandler
  {
public:
    /// Minimal spectrum-level context needed for scoring a candidate match.
    struct SpectrumContext
    {
      Size spectrum_peak_count{0};
      double total_spectrum_intensity{0.0};
      double spectrum_mz_span{1.0};
      double fragment_tolerance_da{0.0};
    };

    /// Minimal candidate-level context needed for scoring a candidate match.
    struct CandidateContext
    {
      Size theoretical_b_ions{0};
      Size theoretical_y_ions{0};
    };

    /// Full score breakdown for one candidate-spectrum match.
    struct Score
    {
      double total_score{0.0};
      double matched_intensity_fraction{0.0};
      double intensity_component{0.0};
      double mz_fidelity_component{0.0};
      double complementarity_component{0.0};
      double coverage_component{0.0};
      double poisson_component{0.0};
      double poisson_proxy{0.0};
      double mean_abs_mz_error_da{0.0};
      double longest_y_pct{0.0};
      Size matched_b_ions{0};
      Size matched_y_ions{0};
      Size longest_b_run{0};
      Size longest_y_run{0};
      Size complementary_pairs{0};
    };

    /// Default constructor
    TagLikeFragmentPatternScorer();

    /// Destructor
    ~TagLikeFragmentPatternScorer() override = default;

    /// Synchronize members with the parameter object.
    void updateMembers_() override;

    /**
      @brief Score one fragment-pattern match.

      @param[in] match FragmentIndex spectrum match
      @param[in] spectrum_context Spectrum-level context
      @param[in] candidate_context Candidate-level context
      @return Score breakdown and total score
    */
    Score score(const FragmentIndex::SpectrumMatch& match,
                const SpectrumContext& spectrum_context,
                const CandidateContext& candidate_context) const;

private:
    double intensity_scale_{1000.0};
    double coverage_weight_{1.0};
    double intensity_weight_{0.5};
    double mz_fidelity_weight_{1.0};
    double complementarity_weight_{0.75};
    double poisson_weight_{0.35};
    double coverage_run_bonus_weight_{0.35};
  };
} // namespace OpenMS
