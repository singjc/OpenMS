// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <results/Quantification.h>
#include <set>
#include <stdexcept>

using namespace OpenMS::OpenDIA;

namespace
{
int failures = 0;
void check(bool condition, const std::string& message)
{
  if (! condition)
  {
    ++failures;
    std::cerr << "FAIL: " << message << '\n';
  }
}
void near(double actual, double expected, const std::string& message)
{
  check(std::isfinite(actual) && std::abs(actual - expected) < 1e-12, message);
}
template<typename F>
void invalid(F action, const std::string& message)
{
  try
  {
    action();
    check(false, message);
  }
  catch (const std::invalid_argument&)
  {
  }
}
void pair(AnalysisLibrary& lib, const std::string& sequence, const std::string& proteins)
{
  const auto p = static_cast<std::uint32_t>(lib.precursors.size());
  const auto group = static_cast<std::uint32_t>(lib.protein_sets.size());
  lib.protein_sets.push_back(proteins);
  for (int d = 0; d < 2; ++d)
  {
    AnalysisLibrary::Precursor precursor;
    precursor.sequence = d ? "EDITPEP" : sequence;
    precursor.proteins = group;
    precursor.decoy = d;
    precursor.charge = 2;
    precursor.pair = static_cast<std::uint32_t>(lib.pairs);
    precursor.decoy_partner = d ? UINT32_MAX : p + 1;
    lib.precursors.push_back(precursor);
  }
  ++lib.pairs;
}
void hit(RunResult& run, std::uint32_t p, float score, double area = 1, double q = 0)
{
  RunWindow w;
  w.precursor = p;
  w.score = score;
  w.q_value = q;
  w.peak.area = area;
  run.windows.push_back(w);
}
const QuantifiedEntity& peptide(const Quantification& result, std::size_t p)
{
  return result.peptides.at(result.precursor_peptide.at(p));
}
const QuantifiedEntity& protein(const Quantification& result, std::size_t p)
{
  return result.protein_groups.at(result.precursor_protein_group.at(p));
}
} // namespace

int main()
{
  const ReportingOptions permissive {1, 1, 1};
  AnalysisLibrary lib;
  pair(lib, "PEPTIDE", "P1");
  pair(lib, "PEPTIDK", "P2");
  pair(lib, "PEPTIDR", "P3");
  pair(lib, "PEPTIDA", "P4");
  pair(lib, "PEPTIDG", "P5");
  std::vector<RunResult> runs(2);
  // Independent hand-count oracle: picked ranks T10,T9,D8,T7,T6.
  // Prefix (D+1)/T: 1, 1/2, 1, 2/3, 1/2. Reverse minima: all 1/2.
  const float target[] = {10, 9, 1, 7, 6};
  const float decoy[] = {0, 0, 8, 0, 0};
  for (std::uint32_t i = 0; i < 5; ++i)
  {
    hit(runs[0], 2 * i, target[i], 10 + i);
    hit(runs[0], 2 * i + 1, decoy[i], 20 + i);
    hit(runs[1], 2 * i, 4, 30 + i);
    hit(runs[1], 2 * i + 1, 4, 40 + i);
  }
  auto result = quantify(lib, runs, permissive);
  for (std::size_t i = 0; i < 5; ++i)
  {
    const auto winner = 2 * i + (i == 2);
    near(peptide(result, winner).runs[0].q_value, .5, "hand-count peptide run q");
    near(protein(result, winner).runs[0].q_value, .5, "hand-count protein run q");
    near(peptide(result, winner).global_q_value, .5, "hand-count peptide global q");
    near(result.global_precursor_q[winner], .5, "legacy precursor hand-count q");
    near(peptide(result, 2 * i).runs[1].q_value, 1, "all-run target ties lose");
    near(peptide(result, 2 * i + 1).runs[1].q_value, 1, "all-run decoy ties win conservatively");
  }
  check(result.peptides.size() == 10 && result.protein_groups.size() == 10, "target and decoy entities separate");
  auto filtered = quantify(lib, runs, {.49, .49, .49});
  check(! peptide(filtered, 0).runs[0].has_area, "precursor global threshold suppresses area");
  near(peptide(filtered, 0).runs[0].q_value, .5, "q estimates precede reporting thresholds");

  // Equal-score block T10,T9 then T8,D8: terminal 2/3 must apply to BOTH score-8 winners,
  // not the anti-conservative intermediate T-only estimate 1/3.
  runs.resize(1);
  runs[0].windows.clear();
  for (std::uint32_t i = 0; i < 4; ++i)
  {
    hit(runs[0], 2 * i, i == 3 ? 0 : 10 - static_cast<float>(i));
    hit(runs[0], 2 * i + 1, i == 3 ? 8 : 0);
  }
  result = quantify(lib, runs, permissive);
  near(peptide(result, 4).global_q_value, 2.0 / 3, "equal-score block target conservatism");
  near(peptide(result, 7).global_q_value, 2.0 / 3, "equal-score block decoy same estimate");
  near(peptide(result, 0).global_q_value, .5, "monotonic reverse minimum");
  check(std::isnan(peptide(result, 8).global_q_value), "unobserved entity has NA global q");
  check(peptide(result, 8).runs.empty(), "unobserved entity has sparse empty runs");

  // Invalid/non-peptide entries must not invalidate unchanged precursor-level reporting.
  auto baseline = quantify(lib, runs, permissive);
  lib.precursors[0].sequence = "not-a-peptide!";
  lib.precursors[2].sequence.clear();
  auto unmappable = quantify(lib, runs, permissive);
  check(unmappable.precursor_peptide[0] == UINT32_MAX && unmappable.precursor_peptide[1] == UINT32_MAX
          && unmappable.precursor_peptide[2] == UINT32_MAX && unmappable.precursor_peptide[3] == UINT32_MAX,
        "invalid and empty target sequences and paired decoys have no family");
  check(unmappable.precursor_protein_group[0] == UINT32_MAX && unmappable.precursor_protein_group[2] == UINT32_MAX,
        "unmappable sequences have no protein group");
  check(std::memcmp(baseline.global_precursor_q.data(), unmappable.global_precursor_q.data(), baseline.global_precursor_q.size() * sizeof(double))
          == 0,
        "unmappable sequences preserve precursor q bytes exactly");
  lib.precursors[0].sequence = "PEPTIDE";
  lib.precursors[2].sequence = "PEPTIDK";
  AnalysisLibrary collisions;
  pair(collisions, "PEPTIDE", "A");
  pair(collisions, "PEPTIDK", "DECOY_A");
  auto collision_result = quantify(collisions, {}, permissive);
  std::set<std::string> ids;
  for (const auto& group : collision_result.protein_groups)
    ids.insert(group.id);
  check(ids.size() == collision_result.protein_groups.size(), "target and derived decoy group IDs never collide");
  check(protein(collision_result, 0).id == "A" && protein(collision_result, 2).id == "DECOY_A", "escaping decoy IDs preserves original target IDs");
  AnalysisLibrary mapping;
  pair(mapping, "MPEPTIDE", "A;B");
  pair(mapping, "M(Oxidation)PEPTIDE", "A;B");
  mapping.precursors[2].charge = 3;
  pair(mapping, "PEPTIDK", "A;B;C");
  pair(mapping, "PEPTIDR", "C");
  pair(mapping, "PEPTIDA", "NA");
  pair(mapping, "PEPTIDG", "");
  // Arbitrary decoy accession strings must have no influence on group inference.
  mapping.protein_sets.push_back("unrelated_decoy_accession");
  for (std::size_t p = 1; p < mapping.precursors.size(); p += 2)
    mapping.precursors[p].proteins = static_cast<std::uint32_t>(mapping.protein_sets.size() - 1);
  std::vector<RunResult> mapped_runs(1);
  for (std::uint32_t p = 0; p < mapping.precursors.size(); ++p)
    hit(mapped_runs[0], p, p % 2 ? 0 : 10, 2);
  hit(mapped_runs[0], 0, 9, 500);  // Lower candidate must not be counted.
  hit(mapped_runs[0], 0, 10, 700); // Equal repeated winner keeps the first row.
  result = quantify(mapping, mapped_runs, permissive);
  check(result.precursor_peptide[0] == result.precursor_peptide[2], "mods and charge collapse to parent AASequence");
  check(result.precursor_peptide[1] == result.precursor_peptide[3], "decoys use paired parent, not decoy sequence");
  check(protein(result, 0).proteins == std::vector<std::string>({"A", "B"}), "identical library evidence resolves indistinguishable proteins");
  check(result.precursor_protein_group[4] == UINT32_MAX, "shared peptide excluded from resolved groups");
  check(result.precursor_protein_group[5] == UINT32_MAX, "shared paired decoy excluded too");
  check(result.precursor_protein_group[8] == UINT32_MAX && result.precursor_protein_group[10] == UINT32_MAX,
        "NA and empty protein mappings excluded");
  near(peptide(result, 0).runs[0].area, 4, "sum raw areas across mods without repeated candidates");
  near(protein(result, 0).runs[0].area, 4, "group quantity excludes shared peptide");
  near(protein(result, 6).runs[0].area, 2, "second group excludes shared peptide");
  check(protein(result, 1).decoy && protein(result, 1).proteins == protein(result, 0).proteins, "decoy mirrors parent groups");
  auto strict_peptide = quantify(mapping, mapped_runs, {1, 0, 1});
  check(peptide(strict_peptide, 0).runs[0].has_area && protein(strict_peptide, 0).runs[0].has_area,
        "tightening peptide q does not suppress precursor-qualified protein area");
  near(protein(strict_peptide, 0).runs[0].area, protein(result, 0).runs[0].area, "protein area independent of peptide cutoff");
  near(protein(strict_peptide, 0).runs[0].q_value, protein(result, 0).runs[0].q_value, "protein run q independent of peptide cutoff");
  near(protein(strict_peptide, 0).global_q_value, protein(result, 0).global_q_value, "protein global q independent of peptide cutoff");
  auto run_filtered = mapped_runs;
  for (auto& w : run_filtered[0].windows)
    w.q_value = 0.9;
  auto run_cut = quantify(mapping, run_filtered, {0.5, 1, 1});
  check(! peptide(run_cut, 0).runs[0].has_area, "run precursor threshold suppresses area independently");
  near(peptide(run_cut, 0).global_q_value, peptide(result, 0).global_q_value, "run cutoff does not alter estimation");
  std::vector<RunResult> shared_only(1);
  hit(shared_only[0], 4, 100);
  hit(shared_only[0], 5, 0);
  auto shared_result = quantify(mapping, shared_only, permissive);
  check(shared_result.precursor_protein_group[4] == UINT32_MAX, "unobserved unique library evidence still resolves ambiguity");
  for (const auto& group : shared_result.protein_groups)
  {
    check(group.runs.empty() && std::isnan(group.global_q_value), "shared-only evidence cannot contribute group score or invent q");
  }

  // Counterpart exists in library but is absent in this run, and absent globally.
  mapped_runs[0].windows.clear();
  hit(mapped_runs[0], 0, 10, 0);
  result = quantify(mapping, mapped_runs, permissive);
  check(std::isnan(peptide(result, 0).runs[0].q_value) && std::isnan(peptide(result, 0).global_q_value), "missing observed counterpart gives NA");
  check(peptide(result, 0).runs[0].has_area && peptide(result, 0).runs[0].area == 0, "zero area is present, distinct from ineligible");
  check(protein(result, 0).runs[0].has_area, "missing peptide counterpart alone does not suppress precursor-qualified protein area");
  check(std::isnan(protein(result, 0).runs[0].q_value), "missing protein counterpart retains NA for writer's own matrix q gate");
  mapped_runs.resize(2);
  hit(mapped_runs[1], 1, 10);
  result = quantify(mapping, mapped_runs, permissive);
  near(peptide(result, 0).global_q_value, 1, "cross-run maxima tied target loses globally");
  check(std::isnan(peptide(result, 0).runs[0].q_value), "cross-run counterpart does not invent run estimate");
  // A missing library decoy also cannot invent q.
  mapping.precursors.pop_back();
  mapped_runs[0].windows.clear();
  mapped_runs[1].windows.clear();
  hit(mapped_runs[0], 10, 10);
  result = quantify(mapping, mapped_runs, permissive);
  check(std::isnan(peptide(result, 10).global_q_value), "missing library counterpart gives NA");

  for (double bad : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
  {
    invalid([&] { quantify(lib, {}, {bad, 1, 1}); }, "invalid precursor cutoff rejected");
    invalid([&] { quantify(lib, {}, {1, bad, 1}); }, "invalid peptide cutoff rejected");
    invalid([&] { quantify(lib, {}, {1, 1, bad}); }, "invalid protein cutoff rejected");
  }
  RunResult bad;
  hit(bad, 0, 10);
  for (double area : {-1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
  {
    bad.windows[0].peak.area = area;
    invalid([&] { quantify(lib, {bad}, permissive); }, "invalid numeric area rejected");
  }
  bad.windows[0].peak.area = 1;
  bad.windows[0].score = std::numeric_limits<float>::quiet_NaN();
  invalid([&] { quantify(lib, {bad}, permissive); }, "NaN score rejected");
  bad.windows[0].score = 10;
  bad.windows[0].q_value = -1;
  invalid([&] { quantify(lib, {bad}, permissive); }, "invalid run q rejected");
  bad.windows[0].q_value = 0;
  bad.windows[0].precursor = UINT32_MAX;
  invalid([&] { quantify(lib, {bad}, permissive); }, "invalid precursor index rejected");
  bad.windows[0].precursor = 0;
  bad.windows[0].score = CrossFitScorer::NOT_KEPT;
  result = quantify(lib, {bad}, permissive);
  check(peptide(result, 0).runs.empty(), "not-retained windows excluded");
  check(quantify(AnalysisLibrary {}, {}).peptides.empty(), "empty input supported");
  std::cout << "Quantification checks: " << failures << " failures\n";
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
