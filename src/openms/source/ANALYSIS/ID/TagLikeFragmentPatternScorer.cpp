// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/TagLikeFragmentPatternScorer.h>

#include <algorithm>
#include <bit>
#include <cmath>

namespace OpenMS
{
  namespace
  {
    Size popcountOrdinalBits_(const std::array<std::uint64_t, 2>& words)
    {
      return static_cast<Size>(std::popcount(words[0]) + std::popcount(words[1]));
    }

    Size longestOrdinalRun_(const std::array<std::uint64_t, 2>& words)
    {
      Size longest_run = 0;
      Size current_run = 0;
      for (Size ordinal_index = 0; ordinal_index < 128; ++ordinal_index)
      {
        const std::uint64_t word = words[ordinal_index / 64];
        const bool is_set = (word & (std::uint64_t{1} << (ordinal_index % 64))) != 0;
        if (is_set)
        {
          ++current_run;
          longest_run = std::max(longest_run, current_run);
        }
        else
        {
          current_run = 0;
        }
      }
      return longest_run;
    }

    bool hasOrdinalBit_(const std::array<std::uint64_t, 2>& words, Size ordinal)
    {
      if (ordinal == 0 || ordinal > 128)
      {
        return false;
      }
      const Size zero_based = ordinal - 1;
      return (words[zero_based / 64] & (std::uint64_t{1} << (zero_based % 64))) != 0;
    }

    Size countComplementaryPairs_(const std::array<std::uint64_t, 2>& matched_b_ordinal_words,
                                  const std::array<std::uint64_t, 2>& matched_y_ordinal_words,
                                  Size theoretical_y_ions)
    {
      if (theoretical_y_ions == 0)
      {
        return 0;
      }

      Size complementary_pairs = 0;
      const Size max_b_ordinal = std::min<Size>(128, theoretical_y_ions);
      for (Size b_ordinal = 1; b_ordinal <= max_b_ordinal; ++b_ordinal)
      {
        if (!hasOrdinalBit_(matched_b_ordinal_words, b_ordinal))
        {
          continue;
        }

        const Size y_ordinal = theoretical_y_ions + 1 - b_ordinal;
        if (hasOrdinalBit_(matched_y_ordinal_words, y_ordinal))
        {
          ++complementary_pairs;
        }
      }
      return complementary_pairs;
    }

    double computePoissonProxy_(Size observed_matches,
                                Size theoretical_fragments,
                                Size spectrum_peak_count,
                                double fragment_tolerance_da,
                                double spectrum_mz_span)
    {
      if (observed_matches == 0 || theoretical_fragments == 0 || spectrum_peak_count == 0)
      {
        return 0.0;
      }

      const double bounded_mz_span = std::max(1.0, spectrum_mz_span);
      const double bounded_fragment_tolerance = std::max(1e-6, fragment_tolerance_da);
      const double expected_match_probability = std::clamp(
        (2.0 * bounded_fragment_tolerance * static_cast<double>(spectrum_peak_count)) / bounded_mz_span,
        1e-9, 0.95);
      const double lambda = std::max(
        1e-9, static_cast<double>(theoretical_fragments) * expected_match_probability);
      const double observed = static_cast<double>(observed_matches);
      if (observed <= lambda)
      {
        return 0.0;
      }
      return observed * std::log(observed / lambda) - (observed - lambda);
    }
  } // namespace

  TagLikeFragmentPatternScorer::TagLikeFragmentPatternScorer() :
    DefaultParamHandler("TagLikeFragmentPatternScorer")
  {
    defaults_.setValue("intensity_scale", intensity_scale_,
                       "Scale factor applied before log-transforming the matched-intensity fraction.");
    defaults_.setMinFloat("intensity_scale", 0.0);
    defaults_.setValue("coverage_weight", coverage_weight_,
                       "Weight applied to matched-ion coverage and ion-run support.");
    defaults_.setMinFloat("coverage_weight", 0.0);
    defaults_.setValue("intensity_weight", intensity_weight_,
                       "Weight applied to the matched-intensity fraction component.");
    defaults_.setMinFloat("intensity_weight", 0.0);
    defaults_.setValue("mz_fidelity_weight", mz_fidelity_weight_,
                       "Weight applied to the fragment m/z fidelity component.");
    defaults_.setMinFloat("mz_fidelity_weight", 0.0);
    defaults_.setValue("complementarity_weight", complementarity_weight_,
                       "Weight applied to complementary b/y cleavage support.");
    defaults_.setMinFloat("complementarity_weight", 0.0);
    defaults_.setValue("poisson_weight", poisson_weight_,
                       "Weight applied to the simple random-match proxy.");
    defaults_.setMinFloat("poisson_weight", 0.0);
    defaults_.setValue("coverage_run_bonus_weight", coverage_run_bonus_weight_,
                       "Bonus weight applied to longest matched-ion runs within the spectrum.");
    defaults_.setMinFloat("coverage_run_bonus_weight", 0.0);
    defaultsToParam_();
    updateMembers_();
  }

  void TagLikeFragmentPatternScorer::updateMembers_()
  {
    intensity_scale_ = static_cast<double>(param_.getValue("intensity_scale"));
    coverage_weight_ = static_cast<double>(param_.getValue("coverage_weight"));
    intensity_weight_ = static_cast<double>(param_.getValue("intensity_weight"));
    mz_fidelity_weight_ = static_cast<double>(param_.getValue("mz_fidelity_weight"));
    complementarity_weight_ = static_cast<double>(param_.getValue("complementarity_weight"));
    poisson_weight_ = static_cast<double>(param_.getValue("poisson_weight"));
    coverage_run_bonus_weight_ = static_cast<double>(param_.getValue("coverage_run_bonus_weight"));
  }

  TagLikeFragmentPatternScorer::Score TagLikeFragmentPatternScorer::score(
    const FragmentIndex::SpectrumMatch& match,
    const SpectrumContext& spectrum_context,
    const CandidateContext& candidate_context) const
  {
    Score score;
    score.matched_b_ions = popcountOrdinalBits_(match.matched_b_ordinal_words_);
    score.matched_y_ions = popcountOrdinalBits_(match.matched_y_ordinal_words_);
    score.longest_b_run = longestOrdinalRun_(match.matched_b_ordinal_words_);
    score.longest_y_run = longestOrdinalRun_(match.matched_y_ordinal_words_);
    score.longest_y_pct =
      candidate_context.theoretical_y_ions > 0 ?
      static_cast<double>(score.longest_y_run) / static_cast<double>(candidate_context.theoretical_y_ions) :
      0.0;

    score.matched_intensity_fraction =
      spectrum_context.total_spectrum_intensity > 0.0 ?
      static_cast<double>(match.matched_intensity_sum_) / spectrum_context.total_spectrum_intensity :
      0.0;

    if (match.num_matched_ > 0)
    {
      score.mean_abs_mz_error_da =
        static_cast<double>(match.matched_abs_mz_error_sum_) / static_cast<double>(match.num_matched_);
    }

    score.complementary_pairs = countComplementaryPairs_(
      match.matched_b_ordinal_words_, match.matched_y_ordinal_words_,
      candidate_context.theoretical_y_ions);

    const Size theoretical_fragments =
      candidate_context.theoretical_b_ions + candidate_context.theoretical_y_ions;
    score.poisson_proxy = computePoissonProxy_(
      score.matched_b_ions + score.matched_y_ions,
      theoretical_fragments,
      spectrum_context.spectrum_peak_count,
      spectrum_context.fragment_tolerance_da,
      spectrum_context.spectrum_mz_span);

    score.coverage_component =
      static_cast<double>(match.num_matched_) +
      coverage_run_bonus_weight_ *
        std::log1p(static_cast<double>(score.longest_b_run + score.longest_y_run));
    score.intensity_component =
      std::log1p(std::max(0.0, intensity_scale_) *
                 std::clamp(score.matched_intensity_fraction, 0.0, 1.0));
    score.mz_fidelity_component =
      spectrum_context.fragment_tolerance_da > 0.0 ?
      std::clamp(1.0 - (score.mean_abs_mz_error_da / spectrum_context.fragment_tolerance_da), 0.0, 1.0) :
      0.0;
    score.complementarity_component =
      std::min(candidate_context.theoretical_b_ions, candidate_context.theoretical_y_ions) > 0 ?
      static_cast<double>(score.complementary_pairs) /
        static_cast<double>(std::min(candidate_context.theoretical_b_ions, candidate_context.theoretical_y_ions)) :
      0.0;
    score.poisson_component = std::log1p(score.poisson_proxy);

    score.total_score =
      coverage_weight_ * score.coverage_component +
      intensity_weight_ * score.intensity_component +
      mz_fidelity_weight_ * score.mz_fidelity_component +
      complementarity_weight_ * score.complementarity_component +
      poisson_weight_ * score.poisson_component;
    return score;
  }
} // namespace OpenMS
