// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  /**
    @brief A helper class that is used by several OpenSWATH tools
  */
  class OPENMS_DLLAPI OpenSwathHelper
  {

public:
    /**
      @brief Strategy used to choose one map when multiple diaPASEF maps are eligible.

      CLOSEST_IM_CENTER preserves the historical behavior. MAXIMUM_IM_OVERLAP
      selects the map containing the largest fraction of the target-centered IM
      extraction interval and uses center distance as a tie-break.
    */
    enum class PasefMapSelectionStrategy
    {
      CLOSEST_IM_CENTER,
      MAXIMUM_IM_OVERLAP
    };

    /**
      @brief One candidate diaPASEF SWATH/IM window for a precursor.

      Used for both stage-internal assignment and debug/audit reporting.
    */
    struct PasefMapCandidate
    {
      int swath_map_index{-1};
      double mz_lower{0.0};
      double mz_upper{0.0};
      double mz_center{0.0};
      double im_lower{-1.0};
      double im_upper{-1.0};
      double im_center{-1.0};
      double upper_edge_distance{0.0};
      double mz_center_distance{0.0};
      double im_center_distance{0.0};
      /// Width of the overlap between the target-centered IM extraction interval and this map.
      double im_overlap_width{0.0};
      /// Fraction of the configured full IM extraction interval available in this map.
      double im_overlap_fraction{0.0};
      bool selected_best{false};
    };

    /**
      @brief All matching diaPASEF SWATH/IM windows for one precursor.
    */
    struct PasefMapMatch
    {
      int selected_swath_map_index{-1};
      std::vector<PasefMapCandidate> candidates;

      bool hasMatch() const
      {
        return selected_swath_map_index >= 0;
      }
    };

    /**
      @brief Compute unique precursor identifier

      Uses transition_group_id and isotope number to compute a unique precursor
      id of the form "groupID_Precursor_ix" where x is the isotope number, e.g.
      the monoisotopic precursor would become "groupID_Precursor_i0".

      @param[in] transition_group_id Unique id of the transition group (peptide/compound)
      @param[in] isotope Precursor isotope number

      @return Unique precursor identifier
    */
    static std::string computePrecursorId(const std::string& transition_group_id, int isotope)
    {
      return transition_group_id + "_Precursor_i" + StringUtils::toStr(isotope);
    }

    /**
      @brief Compute transition group id

      Uses the unique precursor identifier to compute the transition group id
      (peptide/compound identifier), reversing the operation performed by
      computePrecursorId().

      @param[in] precursor_id Precursor identifier as computed by computePrecursorId()

      @return Original transition group id
    */
    static std::string computeTransitionGroupId(const std::string& precursor_id)
    {
      std::vector<std::string> substrings;
      StringUtils::split(precursor_id, "_", substrings);

      if (substrings.size() == 3) return substrings[0];
      else if (substrings.size() > 3)
      {
        std::string r;
        for (Size k = 0; k < substrings.size() - 2; k++) r += substrings[k] + "_";
        return StringUtils::prefix(r, r.size() - 1);
      }
      return "";
    }

    /**
      @brief Select transitions between lower and upper and write them into the new TargetedExperiment

      Version for the OpenMS TargetedExperiment

      @param[in] targeted_exp Transition list for selection
      @param[out] selected_transitions Selected transitions for SWATH window
      @param[in] min_upper_edge_dist Distance in Th to the upper edge
      @param[in] lower Lower edge of SWATH window (in Th)
      @param[in] upper Upper edge of SWATH window (in Th)
    */
    static void selectSwathTransitions(const OpenMS::TargetedExperiment& targeted_exp,
                                       OpenMS::TargetedExperiment& selected_transitions,
                                       double min_upper_edge_dist,
                                       double lower, double upper);

    /**
      @brief Select transitions between lower and upper and write them into the new TargetedExperiment

      Version for the LightTargetedExperiment

      @param[in] targeted_exp Transition list for selection
      @param[out] selected_transitions Selected transitions for SWATH window
      @param[in] min_upper_edge_dist Distance in Th to the upper edge
      @param[in] lower Lower edge of SWATH window (in Th)
      @param[in] upper Upper edge of SWATH window (in Th)
    */
    static void selectSwathTransitions(const OpenSwath::LightTargetedExperiment& targeted_exp,
                                       OpenSwath::LightTargetedExperiment& selected_transitions,
                                       double min_upper_edge_dist,
                                       double lower, double upper);
    /**
     @brief Match transitions with their "best" window across m/z and ion mobility, save results in a vector.

     @param[in] transition_exp Transition list for selection
     @param[out] tr_win_map Mapping from transition (index) to the best matching entry in @p swath_maps
     @param[in] min_upper_edge_dist Distance in Th to the upper edge
     @param[in] swath_maps vector of SwathMap objects defining mz and im bounds
     @param[in] im_match_tolerance Symmetric half-width of the target-centered IM extraction interval
     @param[in] selection_strategy Strategy used to choose one eligible map
    */
    static void selectSwathTransitionsPasef(const OpenSwath::LightTargetedExperiment& transition_exp, std::vector<int>& tr_win_map,
                                       double min_upper_edge_dist, const std::vector< OpenSwath::SwathMap >& swath_maps,
                                       double im_match_tolerance = 0.0,
                                       PasefMapSelectionStrategy selection_strategy = PasefMapSelectionStrategy::CLOSEST_IM_CENTER);

    /**
      @brief Convert a full IM extraction window into a symmetric diaPASEF map-matching slack.

      OpenSWATH IM extraction windows are configured as full widths, while
      diaPASEF map matching compares a single precursor IM against window
      bounds. Using half the extraction window as a matching slack prevents
      hard drops when the library IM scale and the raw-data IM scale differ
      slightly, without changing the downstream extraction width.

      @param[in] im_extraction_window Full IM extraction window width
      @return Half-window matching slack, or 0 if the input window is non-positive/non-finite
    */
    static double computePasefMapMatchingImTolerance(double im_extraction_window);

    /**
      @brief Parse a user-facing diaPASEF map-selection strategy.

      Supported values are @c closest_im_center and @c maximum_im_overlap.
    */
    static PasefMapSelectionStrategy pasefMapSelectionStrategyFromString(const std::string& strategy);

    /// Convert a diaPASEF map-selection strategy to its user-facing string.
    static std::string pasefMapSelectionStrategyToString(PasefMapSelectionStrategy strategy);

    /**
      @brief Check whether a precursor falls into one diaPASEF SWATH/IM window.

      @param[in] swath_map SWATH/IM window to test
      @param[in] precursor_mz Precursor m/z
      @param[in] precursor_im Precursor ion mobility
      @param[in] min_upper_edge_dist Distance in Th to the upper edge
      @param[in] include_upper_bound If true, allow precursor values exactly on the upper m/z/IM boundary
      @param[in] im_match_tolerance Symmetric IM slack applied to the lower/upper window bounds during matching
      @return True if the precursor is inside the window
    */
    static bool pasefSwathMapContainsPrecursor(const OpenSwath::SwathMap& swath_map,
                                               double precursor_mz,
                                               double precursor_im,
                                               double min_upper_edge_dist,
                                               bool include_upper_bound = false,
                                               double im_match_tolerance = 0.0);

    /**
      @brief Collect all matching diaPASEF SWATH/IM windows for one precursor and mark the best match.

      Candidate eligibility is determined by m/z containment and by overlap of
      the target-centered IM extraction interval with the acquired map. The
      selected candidate is controlled by @p selection_strategy. Exact ties keep
      the earliest encountered SWATH map to preserve deterministic behavior.

      @param[in] precursor_mz Precursor m/z
      @param[in] precursor_im Precursor ion mobility
      @param[in] min_upper_edge_dist Distance in Th to the upper edge
      @param[in] swath_maps SWATH/IM windows
      @param[in] include_upper_bound If true, allow precursor values exactly on the upper m/z/IM boundary
      @param[in] im_match_tolerance Symmetric IM half-window used for candidate eligibility and overlap calculation
      @param[in] selection_strategy Strategy used to choose one candidate map
      @return Matching windows and the selected best match
    */
    static PasefMapMatch matchPasefSwathMaps(double precursor_mz,
                                             double precursor_im,
                                             double min_upper_edge_dist,
                                             const std::vector<OpenSwath::SwathMap>& swath_maps,
                                             bool include_upper_bound = false,
                                             double im_match_tolerance = 0.0,
                                             PasefMapSelectionStrategy selection_strategy = PasefMapSelectionStrategy::CLOSEST_IM_CENTER);

    /**
      @brief Get the lower / upper offset for this SWATH map and do some sanity checks

     
      Sanity check for the whole map:
       - all scans need to have exactly one precursor
       - all scans need to have the same MS levels (otherwise extracting an XIC
         from them makes no sense)
       - all scans need to have the same precursor isolation window (otherwise
         extracting an XIC from them makes no sense)

      @param[in] swath_map Input SWATH map to check
      @param[out] lower Lower edge of SWATH window (in Th)
      @param[out] upper Upper edge of SWATH window (in Th)
      @param[out] center Center of SWATH window (in Th)

      @throw IllegalArgument exception if the sanity checks fail.
    */
    static void checkSwathMap(const OpenMS::PeakMap& swath_map,
                              double& lower, double& upper, double& center);

    /**
      @brief Check the map and select transition in one function

      Computes lower and upper offset for the SWATH map and performs some
      sanity checks (see checkSwathMap()). Then selects transitions.

      @param[in] exp Input SWATH map to check
      @param[in] targeted_exp Transition list for selection
      @param[out] selected_transitions Selected transitions for SWATH window
      @param[in] min_upper_edge_dist Distance in Th to the upper edge
    */
    template <class TargetedExperimentT>
    static bool checkSwathMapAndSelectTransitions(const OpenMS::PeakMap& exp,
                                                  const TargetedExperimentT& targeted_exp,
                                                  TargetedExperimentT& selected_transitions,
                                                  double min_upper_edge_dist)
    {
      if (exp.empty() || exp[0].getPrecursors().empty())
      {
        std::cerr << "WARNING: File " << exp.getLoadedFilePath()
                  << " does not have any experiments or any precursors. Is it a SWATH map? "
                  << "I will move to the next map."
                  << std::endl;
        return false;
      }
      double upper, lower, center;
      OpenSwathHelper::checkSwathMap(exp, lower, upper, center);
      OpenSwathHelper::selectSwathTransitions(targeted_exp, selected_transitions, min_upper_edge_dist, lower, upper);
      if (selected_transitions.getTransitions().empty())
      {
        std::cerr << "WARNING: For File " << exp.getLoadedFilePath()
                  << " no transition were within the precursor window of " << lower << " to " << upper
                  << std::endl;
        return false;
      }
      return true;

    }

    /**
      @brief Computes the min and max retention time value
      
      Estimate the retention time span of a targeted experiment by returning
      the min/max values in retention time as a pair.

      @return A std::pair that contains (min,max)

    */
    static std::pair<double,double> estimateRTRange(const OpenSwath::LightTargetedExperiment & exp);

    /**
     * @brief Sample a subset of peptides uniformly across the RT range.
     *
     * Splits the RT span (min→max) into @p bins and randomly picks up to
     * @p peptides_per_bin compounds from each bin. Useful for on-the-fly
     * iRT calibration without external .irt files.
     *
     * @param[in] exp               Full LightTargetedExperiment (the input peptide query parameter assay list for targeted extraction)
     * @param[in] bins              Number of retention‐time bins (i.e. 10 bins across the RT range for linear iRT, 500-1000 bins across the RT for nonlinear iRT)
     * @param[in] peptides_per_bin  How many peptides to draw per bin (i.e. 5 peptides for linear iRT, 25 - 50 for non-linear iRT)
     * @param[in] seed              If non‐zero, used to seed the RNG (deterministic).
     *                              If zero, will use std::random_device for non-deterministic.
     * @param[in] sort_by_intensity     Whether to sort the assays by the highest cumulative intense transitions. This is useful for sampling the most intense peptides for iRTs.
     * @param[in] top_fraction     Only sample from the top N fraction of sorted assays to narrow down on only really intense peptides. This is useful for selecting a few "high quality" peptides to use for linear iRTs.
     * @param[in] priority_peptides Optional set of peptide sequences to prioritize during sampling. If provided, these peptides
     *                              will be sampled first if found in @p exp, before the remaining quota is filled with regular sampling.
     *                              Useful for ensuring common iRT peptides (e.g., from irtkit or cirtkit) are included when present.
     *
     * @return A new LightTargetedExperiment containing only the sampled
     *         compounds, their transitions, and associated proteins.
     */
    static OpenSwath::LightTargetedExperiment sampleExperiment(
      const OpenSwath::LightTargetedExperiment & exp,
      Size bins,
      Size peptides_per_bin,
      unsigned int seed = 0,
      bool sort_by_intensity = false,
      double top_fraction = 1.0,
      const std::unordered_set<std::string> & priority_peptides = std::unordered_set<std::string>());

    /**
      @brief Returns the feature with the highest score for each transition group.
      
      Simple method to extract the best feature for each transition group (e.g.
      for RT alignment). A quality cutoff can be used to skip some low-quality
      features altogether.

      @param[in] transition_group_map Input data containing the picked and scored map
      @param[in] useQualCutoff Whether to apply a quality cutoff to the data
      @param[in] qualCutoff What quality cutoff should be applied (all data above the cutoff will be kept)

      @return Result of the best scoring peaks (stored as map of peptide id and RT)

    */
    static std::map<std::string, double> simpleFindBestFeature(const OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, 
                                                               bool useQualCutoff = false,
                                                               double qualCutoff = 0.0);
  };

} // namespace OpenMS
