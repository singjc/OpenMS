// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <random>
#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace OpenMS
{
  namespace
  {
    constexpr double PASEF_MAP_SELECTION_EPSILON = 1e-12;

    bool isBetterPasefMapCandidate_(const OpenSwathHelper::PasefMapCandidate& candidate,
                                    const OpenSwathHelper::PasefMapCandidate& best,
                                    OpenSwathHelper::PasefMapSelectionStrategy strategy)
    {
      if (strategy == OpenSwathHelper::PasefMapSelectionStrategy::MAXIMUM_IM_OVERLAP)
      {
        const double overlap_delta = candidate.im_overlap_width - best.im_overlap_width;
        if (overlap_delta > PASEF_MAP_SELECTION_EPSILON)
        {
          return true;
        }
        if (overlap_delta < -PASEF_MAP_SELECTION_EPSILON)
        {
          return false;
        }
      }

      // CLOSEST_IM_CENTER is the primary criterion for the historical strategy
      // and the deterministic tie-break for equal maximum overlap. Exact ties
      // retain the earliest map because candidates are visited in map order.
      return candidate.im_center_distance + PASEF_MAP_SELECTION_EPSILON < best.im_center_distance;
    }
  }

  double OpenSwathHelper::computePasefMapMatchingImTolerance(double im_extraction_window)
  {
    return (std::isfinite(im_extraction_window) && im_extraction_window > 0.0) ?
      im_extraction_window / 2.0 :
      0.0;
  }

  OpenSwathHelper::PasefMapSelectionStrategy OpenSwathHelper::pasefMapSelectionStrategyFromString(const std::string& strategy)
  {
    if (strategy == "closest_im_center")
    {
      return PasefMapSelectionStrategy::CLOSEST_IM_CENTER;
    }
    if (strategy == "maximum_im_overlap")
    {
      return PasefMapSelectionStrategy::MAXIMUM_IM_OVERLAP;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Unsupported diaPASEF map-selection strategy", strategy);
  }

  std::string OpenSwathHelper::pasefMapSelectionStrategyToString(PasefMapSelectionStrategy strategy)
  {
    switch (strategy)
    {
      case PasefMapSelectionStrategy::CLOSEST_IM_CENTER:
        return "closest_im_center";
      case PasefMapSelectionStrategy::MAXIMUM_IM_OVERLAP:
        return "maximum_im_overlap";
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Unsupported diaPASEF map-selection strategy enum", std::to_string(static_cast<int>(strategy)));
  }

  bool OpenSwathHelper::pasefSwathMapContainsPrecursor(const OpenSwath::SwathMap& swath_map,
                                                       double precursor_mz,
                                                       double precursor_im,
                                                       double min_upper_edge_dist,
                                                       bool include_upper_bound,
                                                       double im_match_tolerance)
  {
    if (swath_map.ms1 || precursor_mz <= 0.0 || precursor_im < 0.0 ||
        swath_map.imLower < 0.0 || swath_map.imUpper < 0.0)
    {
      return false;
    }

    const bool mz_in_window = include_upper_bound ?
      (swath_map.lower < precursor_mz && precursor_mz <= swath_map.upper) :
      (swath_map.lower < precursor_mz && precursor_mz < swath_map.upper);
    if (!mz_in_window || std::fabs(swath_map.upper - precursor_mz) < min_upper_edge_dist)
    {
      return false;
    }

    const double effective_im_lower = swath_map.imLower - im_match_tolerance;
    const double effective_im_upper = swath_map.imUpper + im_match_tolerance;
    return include_upper_bound ?
      (effective_im_lower < precursor_im && precursor_im <= effective_im_upper) :
      (effective_im_lower < precursor_im && precursor_im < effective_im_upper);
  }

  OpenSwathHelper::PasefMapMatch OpenSwathHelper::matchPasefSwathMaps(double precursor_mz,
                                                                      double precursor_im,
                                                                      double min_upper_edge_dist,
                                                                      const std::vector<OpenSwath::SwathMap>& swath_maps,
                                                                      bool include_upper_bound,
                                                                      double im_match_tolerance,
                                                                      PasefMapSelectionStrategy selection_strategy)
  {
    PasefMapMatch match;
    SignedSize best_candidate_index = -1;

    for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
    {
      const auto& swath_map = swath_maps[static_cast<Size>(i)];
      if (!pasefSwathMapContainsPrecursor(swath_map, precursor_mz, precursor_im,
                                          min_upper_edge_dist, include_upper_bound,
                                          im_match_tolerance))
      {
        continue;
      }

      PasefMapCandidate candidate;
      candidate.swath_map_index = static_cast<int>(i);
      candidate.mz_lower = swath_map.lower;
      candidate.mz_upper = swath_map.upper;
      candidate.mz_center = swath_map.center;
      candidate.im_lower = swath_map.imLower;
      candidate.im_upper = swath_map.imUpper;
      candidate.im_center = (swath_map.imLower + swath_map.imUpper) / 2.0;
      candidate.upper_edge_distance = std::fabs(swath_map.upper - precursor_mz);
      candidate.mz_center_distance = std::fabs(swath_map.center - precursor_mz);
      candidate.im_center_distance = std::fabs(candidate.im_center - precursor_im);

      const double full_im_extraction_window = 2.0 * im_match_tolerance;
      if (full_im_extraction_window > 0.0 && std::isfinite(full_im_extraction_window))
      {
        const double extraction_im_lower = precursor_im - im_match_tolerance;
        const double extraction_im_upper = precursor_im + im_match_tolerance;
        const double overlap_lower = std::max(extraction_im_lower, swath_map.imLower);
        const double overlap_upper = std::min(extraction_im_upper, swath_map.imUpper);
        candidate.im_overlap_width = std::max(0.0, overlap_upper - overlap_lower);
        candidate.im_overlap_fraction = candidate.im_overlap_width / full_im_extraction_window;

        // A zero-width boundary touch contains no extractable IM interval.
        if (candidate.im_overlap_width <= 0.0)
        {
          continue;
        }
      }

      match.candidates.push_back(candidate);
      const SignedSize current_candidate_index = static_cast<SignedSize>(match.candidates.size()) - 1;
      if (best_candidate_index == -1 ||
          isBetterPasefMapCandidate_(match.candidates[static_cast<Size>(current_candidate_index)],
                                     match.candidates[static_cast<Size>(best_candidate_index)],
                                     selection_strategy))
      {
        best_candidate_index = current_candidate_index;
      }
    }

    if (best_candidate_index >= 0)
    {
      auto& best_candidate = match.candidates[static_cast<Size>(best_candidate_index)];
      best_candidate.selected_best = true;
      match.selected_swath_map_index = best_candidate.swath_map_index;
    }

    return match;
  }

  void OpenSwathHelper::selectSwathTransitions(const OpenMS::TargetedExperiment& targeted_exp,
                                               OpenMS::TargetedExperiment& transition_exp_used, double min_upper_edge_dist,
                                               double lower, double upper)
  {
    transition_exp_used.setPeptides(targeted_exp.getPeptides());
    transition_exp_used.setProteins(targeted_exp.getProteins());
    for (Size i = 0; i < targeted_exp.getTransitions().size(); i++)
    {
      ReactionMonitoringTransition tr = targeted_exp.getTransitions()[i];
      if (lower < tr.getPrecursorMZ() && tr.getPrecursorMZ() < upper &&
          std::fabs(upper - tr.getPrecursorMZ()) >= min_upper_edge_dist)
      {

         OPENMS_LOG_DEBUG << "Adding Precursor with m/z " << tr.getPrecursorMZ() <<  " to swath with mz lower of " << lower << " m/z upper of " << upper;
        transition_exp_used.addTransition(tr);
      }
    }
  }

  // For PASEF experiments it is possible to have DIA windows with the same m/z but different IM.
  // Select one eligible map using the configured strategy.
  void OpenSwathHelper::selectSwathTransitionsPasef(const OpenSwath::LightTargetedExperiment& transition_exp, std::vector<int>& tr_win_map,
                                               double min_upper_edge_dist, const std::vector< OpenSwath::SwathMap >& swath_maps,
                                               double im_match_tolerance,
                                               PasefMapSelectionStrategy selection_strategy)
  {
      OPENMS_PRECONDITION(std::any_of(transition_exp.transitions.begin(), transition_exp.transitions.end(), [](auto i){return i.getPrecursorIM()!=-1;}), "All transitions must have a valid IM value (not -1)");

      tr_win_map.resize(transition_exp.transitions.size(), -1);
      for (Size k = 0; k < transition_exp.transitions.size(); k++)
      {
        const OpenSwath::LightTransition& tr = transition_exp.transitions[k];
        const PasefMapMatch match = matchPasefSwathMaps(tr.getPrecursorMZ(), tr.getPrecursorIM(),
                                                        min_upper_edge_dist, swath_maps, false,
                                                        im_match_tolerance, selection_strategy);
        if (match.hasMatch())
        {
          tr_win_map[k] = match.selected_swath_map_index;
        }

        if (match.candidates.size() > 1)
        {
          const auto selected_it = std::find_if(match.candidates.begin(), match.candidates.end(),
                                                [](const PasefMapCandidate& candidate)
                                                {
                                                  return candidate.selected_best;
                                                });
          if (selected_it != match.candidates.end())
          {
            OPENMS_LOG_DEBUG << "For precursor IM " << tr.getPrecursorIM()
                             << " selecting SWATH map " << selected_it->swath_map_index
                             << " using " << pasefMapSelectionStrategyToString(selection_strategy)
                             << " (IM-center distance " << selected_it->im_center_distance
                             << ", usable IM overlap " << selected_it->im_overlap_width
                             << ", overlap fraction " << selected_it->im_overlap_fraction << ")"
                             << " among " << match.candidates.size() << " matching diaPASEF windows."
                             << std::endl;
          }
        }
      }
    }

  void OpenSwathHelper::checkSwathMap(const OpenMS::PeakMap& swath_map,
                                      double& lower, double& upper, double& center)
  {
    if (swath_map.empty() || swath_map[0].getPrecursors().empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Swath map has no Spectra");
    }
    const std::vector<Precursor>& first_prec = swath_map[0].getPrecursors();
    lower = first_prec[0].getMZ() - first_prec[0].getIsolationWindowLowerOffset();
    upper = first_prec[0].getMZ() + first_prec[0].getIsolationWindowUpperOffset();
    center = first_prec[0].getMZ();
    UInt expected_mslevel = swath_map[0].getMSLevel();

    for (Size index = 0; index < swath_map.size(); index++)
    {
      const std::vector<Precursor>& prec = swath_map[index].getPrecursors();
      if (prec.size() != 1)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Scan " + StringUtils::toStr(index) + " does not have exactly one precursor.");
      }
      if (swath_map[index].getMSLevel() != expected_mslevel)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Scan " + StringUtils::toStr(index) + " if of a different MS level than the first scan.");
      }
      if (
        fabs(prec[0].getMZ() - first_prec[0].getMZ()) > 0.1 ||
        fabs(prec[0].getIsolationWindowLowerOffset() - first_prec[0].getIsolationWindowLowerOffset()) > 0.1 ||
        fabs(prec[0].getIsolationWindowUpperOffset() - first_prec[0].getIsolationWindowUpperOffset()) > 0.1
        )
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Scan " + StringUtils::toStr(index) + " has a different precursor isolation window than the first scan.");
      }
    }
  }

  void OpenSwathHelper::selectSwathTransitions(const OpenSwath::LightTargetedExperiment& targeted_exp,
                                               OpenSwath::LightTargetedExperiment& transition_exp_used, double min_upper_edge_dist,
                                               double lower, double upper)
  {
    std::set<std::string> matching_compounds;
    for (Size i = 0; i < targeted_exp.transitions.size(); i++)
    {
      const OpenSwath::LightTransition& tr = targeted_exp.transitions[i];
      if (lower < tr.getPrecursorMZ() && tr.getPrecursorMZ() < upper &&
          std::fabs(upper - tr.getPrecursorMZ()) >= min_upper_edge_dist)
      {
        transition_exp_used.transitions.push_back(tr);
        matching_compounds.insert(tr.getPeptideRef());
      }
    }
    std::set<std::string> matching_proteins;
    for (Size i = 0; i < targeted_exp.compounds.size(); i++)
    {
      if (matching_compounds.contains(targeted_exp.compounds[i].id))
      {
        transition_exp_used.compounds.push_back( targeted_exp.compounds[i] );
        for (Size j = 0; j < targeted_exp.compounds[i].protein_refs.size(); j++)
        {
          matching_proteins.insert(targeted_exp.compounds[i].protein_refs[j]);
        }
      }
    }
    for (Size i = 0; i < targeted_exp.proteins.size(); i++)
    {
      if (matching_proteins.contains(targeted_exp.proteins[i].id))
      {
        transition_exp_used.proteins.push_back( targeted_exp.proteins[i] );
      }
    }
  }

  std::pair<double,double> OpenSwathHelper::estimateRTRange(const OpenSwath::LightTargetedExperiment & exp)
  {
    if (exp.getCompounds().empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Input list of targets is empty.");
    }
    double max = exp.getCompounds()[0].rt;
    double min = exp.getCompounds()[0].rt;
    for (Size i = 0; i < exp.getCompounds().size(); i++)
    {
      if (exp.getCompounds()[i].rt < min) min = exp.getCompounds()[i].rt;
      if (exp.getCompounds()[i].rt > max) max = exp.getCompounds()[i].rt;
    }
    return std::make_pair(min,max);
  }

  OpenSwath::LightTargetedExperiment
  OpenSwathHelper::sampleExperiment(
    const OpenSwath::LightTargetedExperiment & exp,
    Size bins,
    Size peptides_per_bin,
    unsigned int seed,
    bool sort_by_intensity,
    double top_fraction,
    const std::unordered_set<std::string> & priority_peptides)
  {
    OPENMS_PRECONDITION(bins >= 1, "bins must be >= 1");
    OPENMS_PRECONDITION(peptides_per_bin >= 1, "peptides_per_bin must be >= 1");
    OPENMS_PRECONDITION(!sort_by_intensity || (top_fraction > 0.0 && top_fraction <= 1.0),
                        "top_fraction must be in (0,1] when sort_by_intensity is true");

    // 0) initial candidate selection: exclude decoys 
    std::vector<OpenSwath::LightCompound> candidates;
    std::vector<OpenSwath::LightCompound> priority_candidates;

    std::unordered_set<std::string> good_ids;
    for (auto & tr : exp.getTransitions())
    {
      if (!tr.getDecoy())
        good_ids.insert(tr.getPeptideRef());
    }
    
    // Separate compounds into priority and regular candidates
    for (auto & cmp : exp.getCompounds())
    {
      if (good_ids.contains(cmp.id))
      {
        if (priority_peptides.contains(cmp.sequence))
        {
          priority_candidates.push_back(cmp);
        }
        else
        {
          candidates.push_back(cmp);
        }
      }
    }

    // 1) optionally sort by library intensities and trim to top fraction
    if (sort_by_intensity && top_fraction > 0.0 && top_fraction <= 1.0)
    {
      // sum intensities per peptide across all transitions
      std::unordered_map<std::string, double> intensity_sum;
      for (auto & tr : exp.getTransitions())
      {
        if (!tr.getDecoy())
        {
          intensity_sum[tr.getPeptideRef()] += tr.library_intensity;
        }
      }

      // sort regular candidates by descending sum
      std::sort(candidates.begin(), candidates.end(), [&](auto & a, auto & b) {
        return intensity_sum[a.id] > intensity_sum[b.id];
      });

      // trim to top N%
      Size max_keep = std::max<Size>(1, static_cast<Size>(candidates.size() * top_fraction));
      candidates.resize(std::min(max_keep, candidates.size()));
      
      // Also sort priority candidates by intensity (but don't trim them)
      std::sort(priority_candidates.begin(), priority_candidates.end(), [&](auto & a, auto & b) {
        return intensity_sum[a.id] > intensity_sum[b.id];
      });
    }

    // Combine all available candidates for sampling
    std::vector<OpenSwath::LightCompound> all_candidates;
    all_candidates.reserve(priority_candidates.size() + candidates.size());
    all_candidates.insert(all_candidates.end(), priority_candidates.begin(), priority_candidates.end());
    all_candidates.insert(all_candidates.end(), candidates.begin(), candidates.end());

    if (all_candidates.size() < 3)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Insufficient candidates for sampling: " + StringUtils::toStr(all_candidates.size()) +
                                         " found, minimum 3 required for meaningful iRT calibration.");
    }

    // 2) estimate RT range using all available candidates
    double rt_min = std::numeric_limits<double>::max();
    double rt_max = std::numeric_limits<double>::lowest();
    for (auto & cmp : all_candidates)
    {
      rt_min = std::min(rt_min, cmp.rt);
      rt_max = std::max(rt_max, cmp.rt);
    }
    double bin_width = (rt_max - rt_min) / static_cast<double>(bins);

    // 3) sample priority peptides first, then fill remaining quota with uniform sampling
    std::vector<OpenSwath::LightCompound> picked;
    std::unordered_set<std::string> picked_sequences; // Track which sequences we've already picked
    
    // Initialize OpenMS random shuffler
    Math::RandomShuffler rshuffler;
    rshuffler.seed(seed == 0 ? std::random_device{}() : seed);
    
    // First pass: sample priority peptides
    if (!priority_candidates.empty())
    {
      OPENMS_LOG_DEBUG << "Sampling " << priority_candidates.size() 
                      << " priority peptides from the input experiment" << std::endl;
      
      for (Size b = 0; b < bins; ++b)
      {
        double lo = rt_min + b * bin_width;
        double hi = (b + 1 == bins ? rt_max : lo + bin_width);
        std::vector<OpenSwath::LightCompound> bucket;
        
        for (auto & cmp : priority_candidates)
        {
          if (cmp.rt >= lo && cmp.rt < hi && !picked_sequences.contains(cmp.sequence))
            bucket.push_back(cmp);
        }
        
        if (!bucket.empty())
        {
          rshuffler.portable_random_shuffle(bucket.begin(), bucket.end());
          Size take = std::min(peptides_per_bin, bucket.size());
          for (Size i = 0; i < take; ++i)
          {
            picked.push_back(bucket[i]);
            picked_sequences.insert(bucket[i].sequence);
          }
        }
      }
      
      OPENMS_LOG_DEBUG << "Successfully sampled " << picked.size() 
                      << " priority peptides" << std::endl;
    }
    
    // Second pass: fill remaining quota with regular sampling
    Size total_quota = bins * peptides_per_bin;
    if (picked.size() < total_quota && !candidates.empty())
    {
      OPENMS_LOG_DEBUG << "Filling remaining quota (" << (total_quota - picked.size()) 
                      << " peptides) from regular candidates" << std::endl;
      
      for (Size b = 0; b < bins; ++b)
      {
        double lo = rt_min + b * bin_width;
        double hi = (b + 1 == bins ? rt_max : lo + bin_width);
        
        // Count how many priority peptides we already have in this bin
        Size priority_in_bin = 0;
        for (auto & p : picked)
        {
          if (p.rt >= lo && p.rt < hi)
            priority_in_bin++;
        }
        
        // Calculate how many more we need from this bin
        Size needed = (peptides_per_bin > priority_in_bin) ? (peptides_per_bin - priority_in_bin) : 0;
        
        if (needed > 0)
        {
          std::vector<OpenSwath::LightCompound> bucket;
          for (auto & cmp : candidates)
          {
            if (cmp.rt >= lo && cmp.rt < hi && !picked_sequences.contains(cmp.sequence))
              bucket.push_back(cmp);
          }
          
          if (!bucket.empty())
          {
            rshuffler.portable_random_shuffle(bucket.begin(), bucket.end());
            Size take = std::min(needed, bucket.size());
            for (Size i = 0; i < take; ++i)
            {
              picked.push_back(bucket[i]);
              picked_sequences.insert(bucket[i].sequence);
            }
          }
        }
      }
    }

    // 4) assemble output experiment
    OpenSwath::LightTargetedExperiment out_exp;
    out_exp.compounds = picked;

    // copy matching transitions, excluding decoys if requested
    std::unordered_set<std::string> pep_ids;
    for (auto & cmp : picked)
      pep_ids.insert(cmp.id);
    for (auto & tr : exp.getTransitions())
    {
      if (pep_ids.count(tr.getPeptideRef()) && (!tr.getDecoy()))
        out_exp.transitions.push_back(tr);
    }

    // copy associated proteins
    std::unordered_set<std::string> prot_ids;
    for (auto & cmp : picked)
      for (auto & pid : cmp.protein_refs)
        prot_ids.insert(pid);
    for (auto & prot : exp.getProteins())
    {
      if (prot_ids.contains(prot.id))
        out_exp.proteins.push_back(prot);
    }

    return out_exp;
  }

  std::map<std::string, double> OpenSwathHelper::simpleFindBestFeature(
      const OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map,
      bool useQualCutoff, double qualCutoff)
  {
    std::map<std::string, double> result;
    for (const auto & trgroup_it : transition_group_map)
    {
      if (trgroup_it.second.getFeatures().empty() ) {continue;}

      // Find the feature with the highest score
      auto bestf = trgroup_it.second.getBestFeature();

      // Skip if we did not find a feature or do not exceed a certain quality
      if (useQualCutoff && bestf.getOverallQuality() < qualCutoff )
      {
        continue;
      }

      // If we have a found a best feature, add it to the vector
      std::string pepref = trgroup_it.second.getTransitions()[0].getPeptideRef();
      result[ pepref ] = bestf.getRT();
    }
    return result;
  }
}
