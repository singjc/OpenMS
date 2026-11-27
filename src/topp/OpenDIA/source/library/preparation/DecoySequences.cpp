// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $

#include <library/preparation/DecoySequences.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string_view>
#include <unordered_set>

namespace OpenMS::OpenDIA::Preparation
{
  namespace
  {
    struct Ion
    {
      double mz;
      int charge;
      Size ordinal;
    };

    // Keep the cumulative summation and suffix subtraction order of the legacy ion series.
    std::map<std::string, Ion> ions(const AASequence& sequence, int precursor_charge)
    {
      std::map<std::string, Ion> result;
      const Size n = sequence.size();
      std::vector<double> cumulative(n + 1, 0.0);
      for (Size i = 0; i < n; ++i)
      {
        cumulative[i + 1] = cumulative[i] + sequence[i].getMonoWeight(Residue::Internal);
      }
      const double n_mod = sequence.hasNTerminalModification() ? sequence.getNTerminalModification()->getDiffMonoMass() : 0.0;
      const double c_mod = sequence.hasCTerminalModification() ? sequence.getCTerminalModification()->getDiffMonoMass() : 0.0;
      const double b_adj = Residue::getInternalToBIon().getMonoWeight();
      const double y_adj = Residue::getInternalToYIon().getMonoWeight();
      const int max_charge = std::min(4, precursor_charge > 0 ? precursor_charge : 1);
      for (const char type : {'b', 'y'})
      {
        for (int charge = 1; charge <= max_charge; ++charge)
        {
          const double protons = Constants::PROTON_MASS_U * static_cast<double>(charge);
          for (Size i = 1; i < n; ++i)
          {
            const double mass = type == 'b' ? cumulative[i] + n_mod : cumulative[n] - cumulative[n - i] + c_mod;
            const double mz = (mass + protons + (type == 'b' ? b_adj : y_adj)) / static_cast<double>(charge);
            result.emplace(std::string(1, type) + std::to_string(i) + "^" + std::to_string(charge),
                           Ion{Math::roundDecimal(mz, -4), charge, i});
          }
        }
      }
      return result;
    }

    bool eligible(const std::map<std::string, Ion>& target, const std::map<std::string, Ion>& decoy, double product_mz)
    {
      double closest = std::numeric_limits<double>::max();
      int best_charge = std::numeric_limits<int>::max();
      Size best_ordinal = std::numeric_limits<Size>::max();
      std::string best;
      constexpr double threshold = 0.025;
      constexpr double epsilon = threshold * 0.01;
      for (const auto& [name, ion] : target)
      {
        const double delta = std::fabs(ion.mz - product_mz);
        if (delta <= threshold &&
            (delta < closest - epsilon ||
             (delta <= closest + epsilon &&
              (ion.charge < best_charge || (ion.charge == best_charge && ion.ordinal < best_ordinal)))))
        {
          closest = delta;
          best_charge = ion.charge;
          best_ordinal = ion.ordinal;
          best = name;
        }
      }
      const auto found = decoy.find(best);
      return found != decoy.end() && found->second.mz > 0;
    }
  }

  Assays decoySequences(const Assays& targets, ProgressLogger::LogType log_type)
  {
    Assays result;
    std::unordered_set<std::string> sequences, excluded;
    std::map<std::string, const Compound*> target_map;
    for (const auto& compound : targets.compounds)
    {
      target_map[compound.id] = &compound;
      if (compound.isPeptide())
      {
        sequences.insert(compound.sequence + std::to_string(compound.charge));
      }
    }
    ProgressLogger logger;
    logger.setLogType(log_type);
    logger.startProgress(0, targets.compounds.size(), "Generating decoy peptides (Light)");
    Size progress = 0;
    for (const auto& target : targets.compounds)
    {
      logger.setProgress(++progress);
      Compound compound = target;
      compound.id = "DECOY_" + compound.id;
      for (auto& ref : compound.protein_refs) { ref = "DECOY_" + ref; }
      if (!compound.peptide_group_label.empty()) { compound.peptide_group_label = "DECOY_" + compound.peptide_group_label; }
      if (!compound.isPeptide())
      {
        result.compounds.push_back(std::move(compound));
        continue;
      }
      AASequence parsed;
      try
      {
        parsed = AASequence::fromString(compound.sequence);
      }
      catch (Exception::InvalidValue&)
      {
        excluded.insert(compound.id);
        continue;
      }
      std::string sequence = parsed.toUnmodifiedString();
      const int n = static_cast<int>(sequence.size());
      if (std::any_of(compound.modifications.begin(), compound.modifications.end(), [n](const Modification& mod)
          { return mod.location == -1 || mod.location == n || mod.location == n - 1; }))
      {
        excluded.insert(compound.id);
        continue;
      }
      if (!sequence.empty())
      {
        std::reverse(sequence.begin(), sequence.end() - 1);
        for (auto& mod : compound.modifications)
        {
          if (mod.location >= 0 && mod.location < n - 1) { mod.location = n - 2 - mod.location; }
        }
        if (sequence.back() == 'K') { sequence.back() = 'R'; }
        else if (sequence.back() == 'R') { sequence.back() = 'K'; }
        else
        {
          std::uint64_t hash = 14695981039346656037ULL;
          for (char c : sequence)
          {
            hash ^= static_cast<std::uint64_t>(c);
            hash *= 1099511628211ULL;
          }
          sequence.back() = std::string_view("ANDCEQGHILMFSTWYV")[hash % 17];
        }
      }
      compound.sequence.clear();
      for (int loc = -1; loc <= n; ++loc)
      {
        if (loc >= 0 && loc < n) { compound.sequence += sequence[loc]; }
        for (const auto& mod : compound.modifications)
        {
          if (mod.location == loc) { compound.sequence += "(UniMod:" + std::to_string(mod.unimod_id) + ")"; }
        }
      }
      // Reserve before annotation, even if a later annotation excludes this compound.
      if (!sequences.insert(compound.sequence + std::to_string(compound.charge)).second)
      {
        excluded.insert(compound.id);
        continue;
      }
      result.compounds.push_back(std::move(compound));
    }
    logger.endProgress();

    std::map<std::string, const Compound*> decoy_map;
    for (const auto& compound : result.compounds) { decoy_map[compound.id] = &compound; }
    std::map<std::string, std::vector<const Transition*>> groups;
    for (const auto& transition : targets.transitions) { groups[transition.peptide_ref].push_back(&transition); }
    logger.startProgress(0, groups.size(), "Generating decoy transitions (Light)");
    progress = 0;
    for (const auto& [ref, transitions] : groups)
    {
      logger.setProgress(++progress);
      const auto target = target_map.find(ref);
      const auto decoy = decoy_map.find("DECOY_" + ref);
      if (target == target_map.end() || decoy == decoy_map.end() || !target->second->isPeptide()) { continue; }
      AASequence target_sequence, decoy_sequence;
      try
      {
        target_sequence = AASequence::fromString(target->second->sequence);
        decoy_sequence = AASequence::fromString(decoy->second->sequence);
      }
      catch (Exception::InvalidValue&)
      {
        // Legacy parsing failure skips transitions, but does not exclude the compound.
        continue;
      }
      const auto target_ions = ions(target_sequence, target->second->charge);
      const auto decoy_ions = ions(decoy_sequence, decoy->second->charge);
      for (const auto* transition : transitions)
      {
        if (transition->detecting && !transition->decoy && !eligible(target_ions, decoy_ions, transition->product_mz))
        {
          excluded.insert(decoy->second->id);
        }
      }
    }
    logger.endProgress();
    std::erase_if(result.compounds, [&excluded](const Compound& compound) { return excluded.contains(compound.id); });
    std::unordered_set<std::string> protein_ids;
    for (const auto& compound : result.compounds)
    {
      for (const auto& ref : compound.protein_refs) { protein_ids.insert(ref); }
    }
    for (const auto& protein : targets.proteins)
    {
      Protein decoy = protein;
      decoy.id = "DECOY_" + protein.id;
      if (protein_ids.contains(decoy.id)) { result.proteins.push_back(std::move(decoy)); }
    }
    return result;
  }
}
