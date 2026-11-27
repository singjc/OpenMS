// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <results/EntityInference.h>
#include <results/Quantification.h>
#include <set>
#include <stdexcept>

namespace OpenMS::OpenDIA
{
namespace
{
  constexpr auto missing = UINT32_MAX;
  const double na = std::numeric_limits<double>::quiet_NaN();
  const double absent = -std::numeric_limits<double>::infinity();

  std::set<std::string> proteins(const std::string& text)
  {
    std::set<std::string> result;
    std::size_t begin = 0;
    while (begin < text.size())
    {
      const auto end = text.find(';', begin);
      auto name = text.substr(begin, end == std::string::npos ? end : end - begin);
      if (! name.empty() && name != "NA" && name != "N/A") result.insert(std::move(name));
      if (end == std::string::npos) break;
      begin = end + 1;
    }
    return result;
  }

  void evidence(QuantifiedEntity& entity, std::size_t run, double score)
  {
    if (entity.runs.empty() || entity.runs.back().run != run) entity.runs.push_back({run, score, na, 0, false});
    else
      entity.runs.back().score = std::max(entity.runs.back().score, score);
  }
} // namespace

Quantification quantify(const AnalysisLibrary& library, const std::vector<RunResult>& runs, const ReportingOptions& options)
{
  for (auto cut : {options.peakgroup_fdr, options.peptide_fdr, options.protein_fdr})
    if (! std::isfinite(cut) || cut < 0 || cut > 1) throw std::invalid_argument("Reporting FDR must be finite and in [0,1]");
  if (options.inference_context != "run-specific" && options.inference_context != "experiment-wide")
    throw std::invalid_argument("Inference context must be run-specific or experiment-wide");
  Quantification out;
  const auto n = library.precursors.size();
  out.precursor_peptide.assign(n, missing);
  out.precursor_protein_group.assign(n, missing);
  std::map<std::pair<bool, std::string>, std::uint32_t> families;
  std::vector<std::set<std::string>> family_proteins;
  for (std::uint32_t p = 0; p < n; ++p)
  {
    const auto& precursor = library.precursors[p];
    if (precursor.pair >= library.pairs) throw std::invalid_argument("Precursor pair outside library");
    std::string sequence;
    try
    {
      sequence = AASequence::fromString(precursor.sequence).toString();
    }
    catch (const Exception::BaseException&)
    {
      // Non-peptide library entries retain precursor-level reporting, but cannot define peptide families.
      continue;
    }
    if (sequence.empty()) continue;
    auto [it, added] = families.emplace(std::make_pair(precursor.decoy, sequence), static_cast<std::uint32_t>(out.peptides.size()));
    if (added)
    {
      out.peptides.push_back({sequence, {}, precursor.decoy, na, {}});
      family_proteins.emplace_back();
    }
    out.precursor_peptide[p] = it->second;
    if (precursor.proteins >= library.protein_sets.size()) throw std::invalid_argument("Protein set outside library");
    const auto names = proteins(library.protein_sets[precursor.proteins]);
    family_proteins[it->second].insert(names.begin(), names.end());
  }
  // Mapping is the union across every library entry of this modified peptide and class.
  std::map<std::pair<bool, std::string>, std::uint32_t> protein_ids;
  std::vector<std::uint32_t> family_group(family_proteins.size(), missing);
  for (std::size_t f = 0; f < family_proteins.size(); ++f)
  {
    out.peptides[f].proteins.assign(family_proteins[f].begin(), family_proteins[f].end());
    for (const auto& name : family_proteins[f])
    {
      const auto decoy = out.peptides[f].decoy;
      auto [it, added] = protein_ids.emplace(std::make_pair(decoy, name), static_cast<std::uint32_t>(out.protein_groups.size()));
      if (added) out.protein_groups.push_back({name, {name}, decoy, na, {}});
      if (family_proteins[f].size() == 1) family_group[f] = it->second;
    }
  }
  // Public IDs must remain unique even for a literal accession/sequence in both classes.
  // Reserve all actual names first so generated IDs cannot steal a real accession.
  auto unique_ids = [](auto& entities) {
    std::set<std::string> reserved, used;
    for (const auto& entity : entities)
      reserved.insert(entity.id);
    for (auto& entity : entities)
    {
      if (! entity.decoy) used.insert(entity.id);
    }
    for (auto& entity : entities)
    {
      if (! entity.decoy) continue;
      if (used.count(entity.id))
      {
        std::string id = "DECOY_" + entity.id;
        while (reserved.count(id) || used.count(id))
          id = "DECOY_" + id;
        entity.id = std::move(id);
      }
      used.insert(entity.id);
    }
  };
  unique_ids(out.peptides);
  unique_ids(out.protein_groups);
  for (std::size_t p = 0; p < n; ++p)
  {
    const auto e = out.precursor_peptide[p];
    if (e != missing) out.precursor_protein_group[p] = family_group[e];
  }

  std::vector<double> best(n, absent);
  std::vector<std::map<std::uint32_t, const RunWindow*>> winners(runs.size());
  for (std::size_t r = 0; r < runs.size(); ++r)
    for (const auto& w : runs[r].windows)
    {
      if (w.precursor >= n || ! std::isfinite(w.score) || ! std::isfinite(w.q_value) || w.q_value < 0 || w.q_value > 1 || ! std::isfinite(w.peak.area)
          || w.peak.area < 0)
        throw std::invalid_argument("Invalid quantification window index, score, q-value or area");
      if (w.score <= CrossFitScorer::NOT_KEPT) continue;
      best[w.precursor] = std::max(best[w.precursor], static_cast<double>(w.score));
      auto [it, added] = winners[r].emplace(w.precursor, &w);
      // A repeated candidate never contributes area twice. Equal scores retain the first input row.
      if (! added && w.score > it->second->score) it->second = &w;
    }

  // Preserve the existing CLI global_precursor_q algorithm and its ordering exactly.
  out.global_precursor_q.assign(n, 1.0);
  {
    std::vector<std::int64_t> pair_winner(library.pairs, -1);
    for (std::uint32_t p = 0; p < library.precursors.size(); ++p)
    {
      if (! std::isfinite(best[p])) continue;
      auto& w = pair_winner[library.precursors[p].pair];
      if (w < 0 || best[p] > best[w] || (best[p] == best[w] && library.precursors[p].decoy)) w = p;
    }
    std::vector<std::uint32_t> ranked;
    for (auto w : pair_winner)
      if (w >= 0) ranked.push_back(static_cast<std::uint32_t>(w));
    std::sort(ranked.begin(), ranked.end(), [&](std::uint32_t a, std::uint32_t b) {
      if (best[a] != best[b]) return best[a] > best[b];
      if (library.precursors[a].decoy != library.precursors[b].decoy) return library.precursors[a].decoy;
      return a < b;
    });
    std::vector<double> q(ranked.size());
    double t = 0.0, d = 0.0;
    for (std::size_t i = 0; i < ranked.size(); ++i)
    {
      (library.precursors[ranked[i]].decoy ? d : t) += 1.0;
      q[i] = std::min(1.0, (d + 1.0) / std::max(t, 1.0));
    }
    for (std::size_t i = q.size(); i-- > 1;)
      q[i - 1] = std::min(q[i - 1], q[i]);
    for (std::size_t i = 0; i < ranked.size(); ++i)
      out.global_precursor_q[ranked[i]] = q[i];
  }
  for (std::size_t r = 0; r < runs.size(); ++r)
    for (const auto& [p, w] : winners[r])
    {
      const auto e = out.precursor_peptide[p], g = out.precursor_protein_group[p];
      if (e != missing) evidence(out.peptides[e], r, w->score);
      if (g != missing) evidence(out.protein_groups[g], r, w->score);
    }
  EntityInference::estimate(out.peptides, runs.size(), "peptide", out.inference_warnings);
  EntityInference::estimate(out.protein_groups, runs.size(), "protein", out.inference_warnings);
  for (std::size_t r = 0; r < runs.size(); ++r)
    for (const auto& [p, w] : winners[r])
    {
      if (w->q_value > options.peakgroup_fdr || out.global_precursor_q[p] > options.peakgroup_fdr) continue;
      const auto e = out.precursor_peptide[p], g = out.precursor_protein_group[p];
      if (e == missing) continue;
      auto& peptide = out.peptides[e];
      auto& pr = *std::find_if(peptide.runs.begin(), peptide.runs.end(), [&](const auto& x) { return x.run == r; });
      pr.area += w->peak.area;
      if (! std::isfinite(pr.area)) throw std::invalid_argument("Peptide area overflow");
      pr.has_area = true;
      // Protein and peptide reporting cuts are independent. Writers apply each entity's own q cuts.
      if (g == missing) continue;
      auto& group = out.protein_groups[g];
      auto& gr = *std::find_if(group.runs.begin(), group.runs.end(), [&](const auto& x) { return x.run == r; });
      gr.area += w->peak.area;
      if (! std::isfinite(gr.area)) throw std::invalid_argument("Protein area overflow");
      gr.has_area = true;
    }
  return out;
}
} // namespace OpenMS::OpenDIA
