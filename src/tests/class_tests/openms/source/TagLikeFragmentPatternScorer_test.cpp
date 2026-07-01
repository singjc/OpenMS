// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/TagLikeFragmentPatternScorer.h>
#include <OpenMS/CONCEPT/ClassTest.h>

using namespace OpenMS;

START_TEST(TagLikeFragmentPatternScorer, "$Id$")

TagLikeFragmentPatternScorer* ptr = nullptr;

START_SECTION(TagLikeFragmentPatternScorer())
{
  ptr = new TagLikeFragmentPatternScorer();
  TEST_NOT_EQUAL(ptr, nullptr)
}
END_SECTION

START_SECTION((Score score(const FragmentIndex::SpectrumMatch&, const SpectrumContext&, const CandidateContext&) const))
{
  auto set_ordinal = [](std::array<std::uint64_t, 2>& words, Size ordinal)
  {
    const Size zero_based = ordinal - 1;
    words[zero_based / 64] |= (std::uint64_t{1} << (zero_based % 64));
  };

  FragmentIndex::SpectrumMatch complementary_match;
  complementary_match.num_matched_ = 4;
  complementary_match.matched_intensity_sum_ = 400.0f;
  complementary_match.matched_abs_mz_error_sum_ = 0.01f;
  set_ordinal(complementary_match.matched_b_ordinal_words_, 1);
  set_ordinal(complementary_match.matched_b_ordinal_words_, 2);
  set_ordinal(complementary_match.matched_y_ordinal_words_, 4);
  set_ordinal(complementary_match.matched_y_ordinal_words_, 3);

  FragmentIndex::SpectrumMatch non_complementary_match = complementary_match;
  non_complementary_match.matched_abs_mz_error_sum_ = 0.20f;
  non_complementary_match.matched_y_ordinal_words_ = {};
  set_ordinal(non_complementary_match.matched_y_ordinal_words_, 1);
  set_ordinal(non_complementary_match.matched_y_ordinal_words_, 2);

  const TagLikeFragmentPatternScorer::SpectrumContext spectrum_context{
    10, 1000.0, 100.0, 0.05
  };
  const TagLikeFragmentPatternScorer::CandidateContext candidate_context{
    4, 4
  };

  const auto complementary_score = ptr->score(
    complementary_match, spectrum_context, candidate_context);
  const auto non_complementary_score = ptr->score(
    non_complementary_match, spectrum_context, candidate_context);

  TEST_REAL_SIMILAR(complementary_score.matched_intensity_fraction, 0.4)
  TEST_REAL_SIMILAR(complementary_score.mean_abs_mz_error_da, 0.0025)
  TEST_EQUAL(complementary_score.matched_b_ions, 2)
  TEST_EQUAL(complementary_score.matched_y_ions, 2)
  TEST_EQUAL(complementary_score.longest_b_run, 2)
  TEST_EQUAL(complementary_score.longest_y_run, 2)
  TEST_EQUAL(complementary_score.complementary_pairs, 2)
  TEST_TRUE(complementary_score.complementarity_component > non_complementary_score.complementarity_component)
  TEST_TRUE(complementary_score.mz_fidelity_component > non_complementary_score.mz_fidelity_component)
  TEST_TRUE(complementary_score.total_score > non_complementary_score.total_score)
}
END_SECTION

delete ptr;
END_TEST
