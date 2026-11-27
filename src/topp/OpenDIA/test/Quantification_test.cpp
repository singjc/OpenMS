// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $

#include <OpenMS/ANALYSIS/OPENSWATH/LevelContextInference.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <results/Quantification.h>
#include <set>
#include <stdexcept>

using namespace OpenMS;
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
  check(std::isfinite(actual) && std::isfinite(expected) && std::abs(actual - expected) < 1e-12, message);
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
void pair(AnalysisLibrary& lib,
          const std::string& sequence,
          const std::string& proteins,
          const std::string& decoy_sequence = "EDITPEP",
          const std::string& decoy_proteins = "rev_P")
{
  const auto p = static_cast<std::uint32_t>(lib.precursors.size());
  for (int d = 0; d < 2; ++d)
  {
    AnalysisLibrary::Precursor precursor;
    precursor.sequence = d ? decoy_sequence : sequence;
    precursor.proteins = static_cast<std::uint32_t>(lib.protein_sets.size());
    lib.protein_sets.push_back(d ? decoy_proteins : proteins);
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
std::string sequence(std::size_t i)
{
  std::string result = "PEPTIDE";
  for (int digit = 0; digit < 4; ++digit)
  {
    result += "ACDEFGHIKLMNPQRSTVWY"[i % 20];
    i /= 20;
  }
  return result;
}
// Cross-check only: OpenDIA's EntityInference must equal the shared OpenMS LevelContextInference exactly.
// The OpenDIA library itself does not depend on it.
void compareNative(const std::vector<QuantifiedEntity>& entities, InferenceLevel level, InferenceContext context, std::size_t run_id = 0)
{
  std::vector<LevelContextInputRow> rows;
  for (std::size_t e = 0; e < entities.size(); ++e)
  {
    double best = -std::numeric_limits<double>::infinity();
    for (const auto& run : entities[e].runs)
    {
      best = std::max(best, run.score);
      if (context == InferenceContext::Global || (context == InferenceContext::RunSpecific && run.run != run_id)) continue;
      rows.push_back({static_cast<Int64>(run.run), std::to_string(run.run) + ":" + std::to_string(e), static_cast<Int64>(e), entities[e].decoy,
                      run.score, context});
    }
    if (context == InferenceContext::Global && std::isfinite(best))
      rows.push_back({std::nullopt, std::to_string(e), static_cast<Int64>(e), entities[e].decoy, best, context});
  }
  LevelContextInferenceConfig config;
  config.level = level;
  config.context = context;
  // No paired competition and no changes to native ErrorEstimationConfig defaults.
  const auto expected = LevelContextInference::infer(rows, config);
  check(expected.size() == rows.size(), "native reference retains every entity including missing counterparts");
  for (const auto& row : expected)
  {
    const auto& entity = entities.at(static_cast<std::size_t>(row.entity_id));
    if (context == InferenceContext::Global)
    {
      near(entity.global_q_value, row.qvalue, "native global q");
      near(entity.global_p_value, row.pvalue, "native global p");
      near(entity.global_pep, row.pep, "native global PEP");
      continue;
    }
    const auto& run
      = *std::find_if(entity.runs.begin(), entity.runs.end(), [&](const auto& r) { return r.run == static_cast<std::size_t>(*row.run_id); });
    if (context == InferenceContext::RunSpecific)
    {
      near(run.q_value, row.qvalue, "native run q");
      near(run.p_value, row.pvalue, "native run p");
      near(run.pep, row.pep, "native run PEP");
    }
    else
    {
      near(run.experiment_q_value, row.qvalue, "native experiment q");
      near(run.experiment_p_value, row.pvalue, "native experiment p");
      near(run.experiment_pep, row.pep, "native experiment PEP");
    }
  }
}
} // namespace

int main()
{
  const ReportingOptions permissive {1, 1, 1};
  AnalysisLibrary stable;
  std::vector<RunResult> runs(3);
  for (std::uint32_t i = 0; i < 160; ++i)
  {
    pair(stable, sequence(i), "P" + std::to_string(i), sequence(i + 1000), "rev_P" + std::to_string(i));
    for (std::size_t r = 0; r < 2; ++r)
    {
      if (i != 1) hit(runs[r], 2 * i, static_cast<float>(i * .09 + r * ((i % 7) * .03)), 10 + i);
      if (i != 0) hit(runs[r], 2 * i + 1, static_cast<float>(i * .075 + r * ((i % 11) * .02)), 20 + i);
    }
  }
  // A target-only run must not discard either estimable run or pooled/global estimates.
  hit(runs[2], 0, 20, 0);
  auto result = quantify(stable, runs, permissive);
  check(result.peptides.size() == 320 && result.protein_groups.size() == 320, "actual target and decoy entities independent");
  for (auto level : {InferenceLevel::Peptide, InferenceLevel::Protein})
  {
    const auto& entities = level == InferenceLevel::Peptide ? result.peptides : result.protein_groups;
    compareNative(entities, level, InferenceContext::Global);
    compareNative(entities, level, InferenceContext::ExperimentWide);
    compareNative(entities, level, InferenceContext::RunSpecific, 0);
    compareNative(entities, level, InferenceContext::RunSpecific, 1);
  }
  check(std::isfinite(peptide(result, 0).runs[0].q_value), "missing paired decoy does not exclude target from estimable distribution");
  check(std::isfinite(peptide(result, 3).runs[0].q_value), "missing paired target does not exclude decoy from estimable distribution");
  check(std::isnan(peptide(result, 0).runs.back().q_value) && std::isnan(peptide(result, 0).runs.back().p_value)
          && std::isnan(peptide(result, 0).runs.back().pep),
        "unestimable run has NA statistics");
  check(result.inference_warnings.size() == 2, "one explicit warning per unestimable level/run");
  check(result.inference_warnings.front().find("run-specific run 2") != std::string::npos, "warning identifies failed context and run");
  check(peptide(result, 1).runs.empty() && std::isnan(peptide(result, 1).global_q_value), "unobserved entity stays sparse and NA");
  auto selected = quantify(stable, runs, {1, 1, 1, "experiment-wide"});
  near(peptide(selected, 0).runs[0].q_value, peptide(result, 0).runs[0].q_value, "selection does not redefine run q");
  auto filtered = quantify(stable, runs, {0, 0, 0});
  near(peptide(filtered, 0).global_p_value, peptide(result, 0).global_p_value, "all inference precedes reporting cuts");
  near(peptide(filtered, 0).runs[0].score, peptide(result, 0).runs[0].score, "entity scorer winner remains before q cuts");
  check(! peptide(filtered, 0).runs[0].has_area, "precursor cut suppresses area only");

  AnalysisLibrary mapping;
  pair(mapping, "MPEPTIDE", "A", "EDITPEPM", "actual_decoy_A");
  pair(mapping, "M(Oxidation)PEPTIDE", "B", "EDITPEPK", "actual_decoy_B");
  pair(mapping, "MPEPTIDE", "A", "EDITPEPM", "actual_decoy_A");
  mapping.precursors[4].charge = 3;
  pair(mapping, "PEPTIDK", "A;C", "EDITPEPR", "actual_decoy_C");
  pair(mapping, "PEPTIDK", "A", "EDITPEPA", "actual_decoy_D");
  pair(mapping, "PEPTIDR", "NA", "EDITPEPG", "");
  std::vector<RunResult> mapped(1);
  for (std::uint32_t p = 0; p < mapping.precursors.size(); ++p)
    hit(mapped[0], p, 10 + p, 2);
  hit(mapped[0], 0, 9, 500);
  hit(mapped[0], 0, 10, 700);
  result = quantify(mapping, mapped, permissive);
  check(result.precursor_peptide[0] != result.precursor_peptide[2], "modified peptides remain separate");
  check(result.precursor_peptide[0] == result.precursor_peptide[4], "charge states collapse within modified peptide");
  check(result.precursor_peptide[1] == result.precursor_peptide[5], "decoys aggregate by actual modified decoy sequence");
  check(result.precursor_peptide[1] != result.precursor_peptide[3], "distinct actual decoy sequences do not mirror target families");
  check(peptide(result, 0).proteins == std::vector<std::string> {"A"}, "actual target mapping");
  check(peptide(result, 1).proteins == std::vector<std::string> {"actual_decoy_A"}, "actual decoy mapping, not target mirror");
  check(result.precursor_protein_group[6] == UINT32_MAX && result.precursor_protein_group[8] == UINT32_MAX,
        "mapping union across all library entries makes shared peptide ineligible");
  check(result.precursor_protein_group[7] != UINT32_MAX, "unique actual decoy mapping independently eligible");
  check(result.precursor_protein_group[10] == UINT32_MAX && result.precursor_protein_group[11] == UINT32_MAX, "absent mappings ineligible");
  check(protein(result, 0).proteins == std::vector<std::string> {"A"} && protein(result, 1).proteins == std::vector<std::string> {"actual_decoy_A"},
        "singleton protein groups retain actual class-specific accessions");
  near(peptide(result, 0).runs[0].score, 14, "entity winner across charge states before q cuts");
  near(peptide(result, 0).runs[0].area, 4, "sum raw eligible areas with repeated candidates counted once");
  near(protein(result, 0).runs[0].area, 4, "protein area excludes shared peptide");
  near(protein(result, 2).runs[0].area, 2, "modified sequence independently maps to protein");
  auto strict = quantify(mapping, mapped, {1, 0, 0});
  near(protein(strict, 0).runs[0].area, 4, "entity q cuts do not alter raw area eligibility");
  auto run_filtered = mapped;
  for (auto& w : run_filtered[0].windows)
    w.q_value = .9;
  auto run_cut = quantify(mapping, run_filtered, {.5, 1, 1});
  check(! peptide(run_cut, 0).runs[0].has_area && ! protein(run_cut, 0).runs[0].has_area,
        "run precursor cutoff suppresses raw area independently of entity inference");
  near(peptide(run_cut, 0).runs[0].score, peptide(result, 0).runs[0].score, "run cut never changes entity scorer winner");
  // A non-observed library entry still changes peptide uniqueness.
  mapping.protein_sets[mapping.precursors[4].proteins] = "Z";
  mapped[0].windows.erase(std::remove_if(mapped[0].windows.begin(), mapped[0].windows.end(), [](const auto& w) { return w.precursor == 4; }),
                          mapped[0].windows.end());
  auto union_mapping = quantify(mapping, mapped, permissive);
  check(union_mapping.precursor_protein_group[0] == UINT32_MAX, "all same-class library mappings determine uniqueness before observed hits");

  AnalysisLibrary canonical;
  pair(canonical, "M(Oxidation)PEPTIDE", "P1", "EDITPEPM", "rev_P1");
  pair(canonical, "M(UniMod:35)PEPTIDE", "P1", "EDITPEPK", "rev_P2");
  const auto canonical_result = quantify(canonical, {}, permissive);
  check(canonical_result.precursor_peptide[0] == canonical_result.precursor_peptide[2],
        "equivalent modification spellings share canonical modified sequence");

  AnalysisLibrary collisions;
  pair(collisions, "PEPTIDE", "A", "PEPTIDE", "A");
  pair(collisions, "PEPTIDK", "DECOY_A", "EDITPEPK", "rev_B");
  result = quantify(collisions, {}, permissive);
  for (const auto* entities : {&result.peptides, &result.protein_groups})
  {
    std::set<std::string> ids;
    for (const auto& entity : *entities)
      ids.insert(entity.id);
    check(ids.size() == entities->size(), "class IDs remain collision-safe even when accession/sequence is shared");
  }
  check(result.precursor_peptide[0] != result.precursor_peptide[1], "same sequence in opposite classes remains independent");
  check(protein(result, 0).id == "A" && protein(result, 2).id == "DECOY_A", "actual target accession IDs reserved before collision escaping");
  check(protein(result, 1).proteins == std::vector<std::string> {"A"}, "escaped public ID never changes actual decoy accession");

  // Exact unchanged precursor oracle, including legacy rank order and tie behavior.
  AnalysisLibrary legacy;
  std::vector<RunResult> legacy_runs(1);
  const float targets[] = {10, 9, 1, 7, 6};
  const float decoys[] = {0, 0, 8, 0, 0};
  for (std::uint32_t i = 0; i < 5; ++i)
  {
    pair(legacy, sequence(i), "P", sequence(i + 50), "rev_P");
    hit(legacy_runs[0], 2 * i, targets[i]);
    hit(legacy_runs[0], 2 * i + 1, decoys[i]);
  }
  auto baseline = quantify(legacy, legacy_runs, permissive);
  for (std::size_t i = 0; i < 5; ++i)
  {
    near(baseline.global_precursor_q[2 * i + (i == 2)], .5, "unchanged legacy precursor winner q");
    near(baseline.global_precursor_q[2 * i + (i != 2)], 1, "unchanged legacy precursor loser q");
  }
  legacy.precursors[0].sequence = "not-a-peptide!";
  legacy.precursors[2].sequence.clear();
  result = quantify(legacy, legacy_runs, permissive);
  check(result.precursor_peptide[0] == UINT32_MAX && result.precursor_peptide[2] == UINT32_MAX, "invalid/empty sequences excluded");
  check(result.precursor_peptide[1] != UINT32_MAX && result.precursor_peptide[3] != UINT32_MAX,
        "valid decoys never depend on paired target validity");
  check(std::memcmp(result.global_precursor_q.data(), baseline.global_precursor_q.data(), baseline.global_precursor_q.size() * sizeof(double)) == 0,
        "mapping changes preserve precursor q bytes exactly");

  // Missing decoy distribution is explicitly unavailable, never a tuned fallback.
  std::vector<RunResult> tiny(1);
  hit(tiny[0], 0, 10, 0);
  result = quantify(collisions, tiny, permissive);
  check(std::isnan(peptide(result, 0).global_q_value) && std::isnan(peptide(result, 0).global_p_value) && std::isnan(peptide(result, 0).global_pep),
        "unestimable global fields all NA");
  check(std::isnan(peptide(result, 0).runs[0].experiment_q_value), "unestimable experiment fields NA");
  check(result.inference_warnings.size() == 6, "tiny target-only input warns for every failed level/context");
  check(peptide(result, 0).runs[0].has_area && peptide(result, 0).runs[0].area == 0, "zero raw area present despite inference unavailable");

  hit(tiny[0], 1, 0);
  const auto singleton = quantify(collisions, tiny, permissive);
  check(singleton.inference_warnings.size() == 6 && std::isnan(peptide(singleton, 0).runs[0].q_value),
        "one target and one decoy are unestimable for local-FDR KDE, not a crash or fallback");
  check(singleton.inference_warnings.front().find("fewer than two target scores") != std::string::npos,
        "singleton warning explains KDE prerequisite");

  AnalysisLibrary ties;
  std::vector<RunResult> tie_runs(1);
  for (std::uint32_t i = 0; i < 20; ++i)
  {
    pair(ties, sequence(i), "T" + std::to_string(i), sequence(i + 500), "rev_T" + std::to_string(i));
    hit(tie_runs[0], 2 * i, 5);
    hit(tie_runs[0], 2 * i + 1, 5);
  }
  const auto tied = quantify(ties, tie_runs, permissive);
  check(tied.inference_warnings.size() == 6, "all-tie scores warn in every level/context like upstream pi0 <= 0 error");
  check(! tied.inference_warnings.empty() && tied.inference_warnings.front().find("pi0 <= 0") != std::string::npos,
        "all-tie warning names pi0 failure");
  check(std::isnan(peptide(tied, 0).global_q_value) && std::isnan(peptide(tied, 0).global_pep) && std::isnan(peptide(tied, 0).runs[0].q_value)
          && std::isnan(peptide(tied, 0).runs[0].experiment_pep) && std::isnan(protein(tied, 1).runs[0].p_value),
        "all-tie statistics remain NA, not zero");

  for (double bad : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
  {
    invalid([&] { quantify(collisions, {}, {bad, 1, 1}); }, "invalid precursor cutoff rejected");
    invalid([&] { quantify(collisions, {}, {1, bad, 1}); }, "invalid peptide cutoff rejected");
    invalid([&] { quantify(collisions, {}, {1, 1, bad}); }, "invalid protein cutoff rejected");
  }
  invalid([&] { quantify(collisions, {}, {1, 1, 1, "global"}); }, "unsupported matrix context rejected");
  RunResult bad;
  hit(bad, 0, 10);
  for (double area : {-1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
  {
    bad.windows[0].peak.area = area;
    invalid([&] { quantify(collisions, {bad}, permissive); }, "invalid numeric area rejected");
  }
  bad.windows[0].peak.area = 1;
  bad.windows[0].score = std::numeric_limits<float>::quiet_NaN();
  invalid([&] { quantify(collisions, {bad}, permissive); }, "NaN score rejected");
  bad.windows[0].score = 10;
  bad.windows[0].q_value = -1;
  invalid([&] { quantify(collisions, {bad}, permissive); }, "invalid run q rejected");
  bad.windows[0].q_value = 0;
  bad.windows[0].precursor = UINT32_MAX;
  invalid([&] { quantify(collisions, {bad}, permissive); }, "invalid precursor index rejected");
  bad.windows[0].precursor = 0;
  bad.windows[0].score = CrossFitScorer::NOT_KEPT;
  result = quantify(collisions, {bad}, permissive);
  check(peptide(result, 0).runs.empty(), "not-retained windows excluded");
  check(quantify(AnalysisLibrary {}, {}).peptides.empty(), "empty input supported");
  std::cout << "Quantification checks: " << failures << " failures\n";
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
