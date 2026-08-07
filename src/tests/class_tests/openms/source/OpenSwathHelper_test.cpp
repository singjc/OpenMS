// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <boost/assign/std/vector.hpp>

///////////////////////////

START_TEST(OpenSwathHelper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace std;
using namespace OpenMS;
using namespace OpenSwath;

OpenSwathHelper* ptr = nullptr;
OpenSwathHelper* nullPointer = nullptr;

START_SECTION(OpenSwathHelper())
{
  ptr = new OpenSwathHelper();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~OpenSwathHelper())
{
  delete ptr;
}
END_SECTION

START_SECTION(static std::string computePrecursorId(const std::string& transition_group_id, int isotope))
{
  TEST_EQUAL(OpenSwathHelper::computePrecursorId("tr_gr2", 0), "tr_gr2_Precursor_i0")
  TEST_EQUAL(OpenSwathHelper::computePrecursorId("tr_gr2__test", 0), "tr_gr2__test_Precursor_i0")
}
END_SECTION

START_SECTION(static std::string computeTransitionGroupId(const std::string& precursor_id))
{
  TEST_EQUAL(OpenSwathHelper::computeTransitionGroupId("tr_gr2_Precursor_i0"), "tr_gr2")
  TEST_EQUAL(OpenSwathHelper::computeTransitionGroupId("tr_gr2__test_Precursor_i0"), "tr_gr2__test")
}
END_SECTION

START_SECTION(static void selectSwathTransitions(const OpenMS::TargetedExperiment &targeted_exp, OpenMS::TargetedExperiment &transition_exp_used, double min_upper_edge_dist, double lower, double upper))
{
  TargetedExperiment exp1;
  TargetedExperiment exp2;

  ReactionMonitoringTransition tr1;
  ReactionMonitoringTransition tr2;
  ReactionMonitoringTransition tr3;

  tr1.setPrecursorMZ(100.0);
  tr2.setPrecursorMZ(200.0);
  tr3.setPrecursorMZ(300.0);

  std::vector<ReactionMonitoringTransition> transitions;
  transitions.push_back(tr1);
  transitions.push_back(tr2);
  transitions.push_back(tr3);

  exp1.setTransitions(transitions);

  // select all transitions between 200 and 500
  OpenSwathHelper::selectSwathTransitions(exp1, exp2, 1.0, 199.9, 500);
  TEST_EQUAL(exp2.getTransitions().size(), 2)
}
END_SECTION

START_SECTION(static void selectSwathTransitions(const OpenSwath::LightTargetedExperiment &targeted_exp, OpenSwath::LightTargetedExperiment &transition_exp_used, double min_upper_edge_dist, double lower, double upper))
{
  LightTargetedExperiment exp1;
  LightTargetedExperiment exp2;

  LightTransition tr1;
  LightTransition tr2;
  LightTransition tr3;

  tr1.precursor_mz = 100.0;
  tr2.precursor_mz = 200.0;
  tr3.precursor_mz = 300.0;

  std::vector<LightTransition> transitions;
  transitions.push_back(tr1);
  transitions.push_back(tr2);
  transitions.push_back(tr3);

  exp1.transitions = transitions;

  // select all transitions between 200 and 500
  OpenSwathHelper::selectSwathTransitions(exp1, exp2, 1.0, 199.9, 500);
  TEST_EQUAL(exp2.getTransitions().size(), 2)
}
END_SECTION

START_SECTION((static bool pasefSwathMapMatchesPrecursor(const OpenSwath::SwathMap& swath_map, double precursor_mz, double precursor_im, double min_upper_edge_dist, bool include_upper_bound, double im_extraction_window)))
{
  SwathMap swath_map;
  swath_map.lower = 500.0;
  swath_map.upper = 525.0;
  swath_map.center = 512.5;
  swath_map.imLower = 0.60;
  swath_map.imUpper = 0.70;
  swath_map.ms1 = false;

  TEST_TRUE(OpenSwathHelper::pasefSwathMapMatchesPrecursor(swath_map, 510.0, 0.65, 1.0, false))
  TEST_FALSE(OpenSwathHelper::pasefSwathMapMatchesPrecursor(swath_map, 524.5, 0.65, 1.0, false))
  TEST_FALSE(OpenSwathHelper::pasefSwathMapMatchesPrecursor(swath_map, 525.0, 0.65, 0.0, false))
  TEST_TRUE(OpenSwathHelper::pasefSwathMapMatchesPrecursor(swath_map, 525.0, 0.65, 0.0, true))

  // Strict point matching rejects a target below the acquired IM map.
  TEST_FALSE(OpenSwathHelper::pasefSwathMapMatchesPrecursor(swath_map, 510.0, 0.58, 0.0, false))

  // A full 0.06 IM extraction interval centred at 0.58 overlaps [0.60, 0.70].
  TEST_TRUE(OpenSwathHelper::pasefSwathMapMatchesPrecursor(swath_map, 510.0, 0.58, 0.0, false, 0.06))

  // A zero-width boundary touch contains no extractable IM interval.
  TEST_FALSE(OpenSwathHelper::pasefSwathMapMatchesPrecursor(swath_map, 510.0, 0.57, 0.0, false, 0.06))

  // The same interval-overlap rule rescues targets above the map.
  TEST_TRUE(OpenSwathHelper::pasefSwathMapMatchesPrecursor(swath_map, 510.0, 0.72, 0.0, false, 0.06))

  // Non-positive widths retain strict point-matching behavior.
  TEST_FALSE(OpenSwathHelper::pasefSwathMapMatchesPrecursor(swath_map, 510.0, 0.58, 0.0, false, -1.0))

  // Preserve the historical inclusive upper IM bound used during calibration.
  TEST_FALSE(OpenSwathHelper::pasefSwathMapMatchesPrecursor(swath_map, 510.0, 0.70, 0.0, false))
  TEST_TRUE(OpenSwathHelper::pasefSwathMapMatchesPrecursor(swath_map, 510.0, 0.70, 0.0, true))
}
END_SECTION

START_SECTION((static int findBestPasefSwathMap(double precursor_mz, double precursor_im, double min_upper_edge_dist, const std::vector<OpenSwath::SwathMap>& swath_maps, bool include_upper_bound, double im_extraction_window)))
{
  vector<SwathMap> swath_maps(3);
  swath_maps[0].lower = 500.0;
  swath_maps[0].upper = 525.0;
  swath_maps[0].center = 512.5;
  swath_maps[0].imLower = 0.60;
  swath_maps[0].imUpper = 0.70;
  swath_maps[0].ms1 = false;

  swath_maps[1].lower = 500.0;
  swath_maps[1].upper = 525.0;
  swath_maps[1].center = 512.5;
  swath_maps[1].imLower = 0.65;
  swath_maps[1].imUpper = 0.75;
  swath_maps[1].ms1 = false;

  swath_maps[2].lower = 600.0;
  swath_maps[2].upper = 625.0;
  swath_maps[2].center = 612.5;
  swath_maps[2].imLower = 0.60;
  swath_maps[2].imUpper = 0.70;
  swath_maps[2].ms1 = false;

  // Preserve historical closest-IM-centre selection when multiple maps match.
  TEST_EQUAL(OpenSwathHelper::findBestPasefSwathMap(510.0, 0.69, 1.0, swath_maps, false), 1)

  // The configured full extraction interval rescues a target below the map.
  TEST_EQUAL(OpenSwathHelper::findBestPasefSwathMap(510.0, 0.58, 0.0, swath_maps, false, 0.06), 0)

  // A target interval that does not overlap any map remains unmatched.
  TEST_EQUAL(OpenSwathHelper::findBestPasefSwathMap(510.0, 0.56, 0.0, swath_maps, false, 0.06), -1)

  // Exact centre-distance ties retain the earliest map.
  // Use binary-exact fractions so the two calculated centres are exactly equal.
  vector<SwathMap> tied_maps(2);
  tied_maps[0] = swath_maps[0];
  tied_maps[0].imLower = 0.50;
  tied_maps[0].imUpper = 0.75;

  tied_maps[1] = swath_maps[0];
  tied_maps[1].imLower = 0.5625;
  tied_maps[1].imUpper = 0.6875;

  TEST_EQUAL(OpenSwathHelper::findBestPasefSwathMap(510.0, 0.625, 0.0, tied_maps, false), 0)
}
END_SECTION

START_SECTION((static void selectSwathTransitionsPasef(const OpenSwath::LightTargetedExperiment& transition_exp, std::vector<int>& tr_win_map, double min_upper_edge_dist, const std::vector<OpenSwath::SwathMap>& swath_maps, double im_extraction_window)))
{
  vector<SwathMap> swath_maps(1);
  swath_maps[0].lower = 500.0;
  swath_maps[0].upper = 525.0;
  swath_maps[0].center = 512.5;
  swath_maps[0].imLower = 0.60;
  swath_maps[0].imUpper = 0.70;
  swath_maps[0].ms1 = false;

  LightTargetedExperiment transition_exp;
  LightTransition rescued;
  rescued.precursor_mz = 510.0;
  rescued.precursor_im = 0.58;
  transition_exp.transitions.push_back(rescued);

  LightTransition unmatched;
  unmatched.precursor_mz = 510.0;
  unmatched.precursor_im = 0.56;
  transition_exp.transitions.push_back(unmatched);

  vector<int> tr_win_map;
  OpenSwathHelper::selectSwathTransitionsPasef(transition_exp, tr_win_map, 0.0, swath_maps, 0.06);
  TEST_EQUAL(tr_win_map.size(), 2)
  TEST_EQUAL(tr_win_map[0], 0)
  TEST_EQUAL(tr_win_map[1], -1)
}
END_SECTION

START_SECTION( (template < class TargetedExperimentT > static bool checkSwathMapAndSelectTransitions(const OpenMS::PeakMap &exp, const TargetedExperimentT &targeted_exp, TargetedExperimentT &transition_exp_used, double min_upper_edge_dist)))
{
  // tested above already
  NOT_TESTABLE
}
END_SECTION

START_SECTION(static void checkSwathMap(const OpenMS::PeakMap &swath_map, double &lower, double &upper))
{
  OpenMS::PeakMap swath_map;
  OpenMS::MSSpectrum spectrum;
  OpenMS::Precursor prec;
  std::vector<Precursor> precursors;
  prec.setMZ(250);
  prec.setIsolationWindowLowerOffset(50);
  prec.setIsolationWindowUpperOffset(50);
  precursors.push_back(prec);
  spectrum.setPrecursors(precursors);
  swath_map.addSpectrum(spectrum);

  double lower, upper, center;
  OpenSwathHelper::checkSwathMap(swath_map, lower, upper, center);

  TEST_REAL_SIMILAR(lower, 200);
  TEST_REAL_SIMILAR(upper, 300);
  TEST_REAL_SIMILAR(center, 250);
}
END_SECTION

START_SECTION((static std::pair<double,double> estimateRTRange(OpenSwath::LightTargetedExperiment & exp)))
{
  LightTargetedExperiment exp;

  LightCompound pep1;
  LightCompound pep2;
  LightCompound pep3;

  pep1.rt = -100.0;
  pep2.rt = 900.0;
  pep3.rt = 300.0;

  std::vector<LightCompound> peptides;
  peptides.push_back(pep1);
  peptides.push_back(pep2);
  peptides.push_back(pep3);

  exp.compounds = peptides;

  std::pair<double, double> range = OpenSwathHelper::estimateRTRange(exp);
  TEST_REAL_SIMILAR(range.first, -100)
  TEST_REAL_SIMILAR(range.second, 900)
}
END_SECTION

START_SECTION((static std::map<std::string, double> simple_find_best_feature(OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, 
        bool useQualCutoff = false, double qualCutoff = 0.0)))
{
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
