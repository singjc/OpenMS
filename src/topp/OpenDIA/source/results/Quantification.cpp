// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
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

  // Experimental reporting estimator, not an empirical calibration claim. A missing observed
  // counterpart has no estimable q (NaN), rather than being counted as a target victory.
  // Adjacent entities are the target and its library-derived decoy. Equal-score blocks share
  // their terminal (D+1)/T estimate before reverse cumulative minima are applied.
  std::vector<double> picked(const std::vector<double>& scores)
  {
    std::vector<double> q(scores.size(), na);
    std::vector<std::size_t> ranked;
    for (std::size_t i = 0; i + 1 < scores.size(); i += 2)
    {
      if (! std::isfinite(scores[i]) || ! std::isfinite(scores[i + 1])) continue;
      q[i] = q[i + 1] = 1.0;
      ranked.push_back(scores[i] > scores[i + 1] ? i : i + 1);
    }
    std::sort(ranked.begin(), ranked.end(), [&](auto a, auto b) { return scores[a] != scores[b] ? scores[a] > scores[b] : a < b; });
    double targets = 0, decoys = 0;
    std::vector<double> estimates(ranked.size());
    for (std::size_t begin = 0; begin < ranked.size();)
    {
      auto end = begin;
      while (end < ranked.size() && scores[ranked[end]] == scores[ranked[begin]])
      {
        (ranked[end] % 2 ? decoys : targets) += 1;
        ++end;
      }
      const auto estimate = std::min(1.0, (decoys + 1.0) / std::max(targets, 1.0));
      std::fill(estimates.begin() + begin, estimates.begin() + end, estimate);
      begin = end;
    }
    double minimum = 1.0;
    for (std::size_t i = ranked.size(); i-- > 0;)
    {
      minimum = std::min(minimum, estimates[i]);
      q[ranked[i]] = minimum;
    }
    return q;
  }

  void estimate(std::vector<QuantifiedEntity>& entities, std::size_t run_count)
  {
    std::vector<double> best(entities.size(), absent);
    std::vector<std::map<std::size_t, std::size_t>> observed(run_count);
    for (std::size_t e = 0; e < entities.size(); ++e)
    {
      for (std::size_t i = 0; i < entities[e].runs.size(); ++i)
      {
        const auto& r = entities[e].runs[i];
        best[e] = std::max(best[e], r.score);
        observed[r.run].emplace(e, i);
      }
    }
    const auto global = picked(best);
    for (std::size_t e = 0; e < entities.size(); ++e)
      entities[e].global_q_value = global[e];
    for (std::size_t r = 0; r < run_count; ++r)
    {
      std::vector<double> scores(entities.size(), absent);
      for (const auto& [e, i] : observed[r])
        scores[e] = entities[e].runs[i].score;
      const auto q = picked(scores);
      for (const auto& [e, i] : observed[r])
        entities[e].runs[i].q_value = q[e];
    }
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
  Quantification out;
  const auto n = library.precursors.size();
  out.precursor_peptide.assign(n, missing);
  out.precursor_protein_group.assign(n, missing);
  std::map<std::uint32_t, std::uint32_t> parents;
  std::map<std::string, std::uint32_t> families;
  std::vector<std::set<std::string>> family_proteins;
  for (std::uint32_t p = 0; p < n; ++p)
  {
    const auto& precursor = library.precursors[p];
    if (precursor.pair >= library.pairs) throw std::invalid_argument("Precursor pair outside library");
    if (precursor.decoy) continue;
    std::string sequence;
    try
    {
      sequence = AASequence::fromString(precursor.sequence).toUnmodifiedString();
    }
    catch (const Exception::BaseException&)
    {
      // Non-peptide library entries retain precursor-level reporting, but cannot define peptide families.
      continue;
    }
    if (sequence.empty()) continue;
    auto [it, added] = families.emplace(sequence, static_cast<std::uint32_t>(out.peptides.size()));
    if (added)
    {
      out.peptides.push_back({sequence, {}, false, na, {}});
      out.peptides.push_back({"DECOY_" + sequence, {}, true, na, {}});
      family_proteins.emplace_back();
    }
    out.precursor_peptide[p] = it->second;
    const auto [parent, inserted] = parents.emplace(precursor.pair, it->second);
    if (! inserted && parent->second != it->second) throw std::invalid_argument("Pair has conflicting target peptide families");
    if (precursor.proteins >= library.protein_sets.size()) throw std::invalid_argument("Protein set outside library");
    const auto names = proteins(library.protein_sets[precursor.proteins]);
    family_proteins[it->second / 2].insert(names.begin(), names.end());
  }
  for (std::size_t p = 0; p < n; ++p)
    if (library.precursors[p].decoy)
    {
      const auto parent = parents.find(library.precursors[p].pair);
      if (parent != parents.end()) out.precursor_peptide[p] = parent->second + 1;
    }

  // Infer indistinguishable proteins from ALL library target evidence, never from identified hits.
  std::map<std::string, std::set<std::size_t>> protein_evidence;
  for (std::size_t f = 0; f < family_proteins.size(); ++f)
  {
    for (const auto& name : family_proteins[f])
      protein_evidence[name].insert(f);
    for (auto e : {2 * f, 2 * f + 1})
      out.peptides[e].proteins.assign(family_proteins[f].begin(), family_proteins[f].end());
  }
  std::map<std::set<std::size_t>, std::vector<std::string>> groups;
  for (const auto& [name, evidence_set] : protein_evidence)
    groups[evidence_set].push_back(name);
  std::vector<std::set<std::uint32_t>> family_groups(family_proteins.size());
  std::set<std::string> reserved_group_ids;
  for (const auto& [evidence_set, names] : groups)
  {
    std::string id;
    for (const auto& name : names)
    {
      if (! id.empty()) id += ';';
      id += name;
    }
    reserved_group_ids.insert(std::move(id));
  }
  for (const auto& [evidence_set, names] : groups)
  {
    const auto g = static_cast<std::uint32_t>(out.protein_groups.size());
    std::string id;
    for (const auto& name : names)
    {
      if (! id.empty()) id += ';';
      id += name;
    }
    out.protein_groups.push_back({id, names, false, na, {}});
    std::string decoy_id = "DECOY_" + id;
    while (reserved_group_ids.count(decoy_id))
      decoy_id = "DECOY_" + decoy_id;
    reserved_group_ids.insert(decoy_id);
    out.protein_groups.push_back({decoy_id, names, true, na, {}});
    for (auto f : evidence_set)
      family_groups[f].insert(g);
  }
  for (std::size_t p = 0; p < n; ++p)
  {
    const auto e = out.precursor_peptide[p];
    if (e != missing && family_groups[e / 2].size() == 1) out.precursor_protein_group[p] = *family_groups[e / 2].begin() + e % 2;
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
  estimate(out.peptides, runs.size());
  estimate(out.protein_groups, runs.size());
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
