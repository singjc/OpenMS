// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $

#include <OpenMS/MATH/STATISTICS/MultipleTesting.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <results/EntityInference.h>
#include <stdexcept>

namespace OpenMS::OpenDIA
{
namespace
{
  using Math::MultipleTesting;

  // PyProphet error_statistics defaults. pi0_lambda is np.arange(0.1, 0.5, 0.05), whose
  // binary64 values are not the decimal literals (e.g. 0.15000000000000002). Empirical p-values
  // are often exactly k/n, so p >= lambda at the boundary needs these exact values.
  const std::vector<double> PI0_LAMBDA {0.1, 0.15000000000000002, 0.20000000000000004, 0.25000000000000006,
                                        0.30000000000000004, 0.3500000000000001, 0.40000000000000013, 0.45000000000000007};
  constexpr double LFDR_ADJ = 1.5;
  constexpr double LFDR_EPS = 1e-8;

  /// One score of a distribution and where its statistics go.
  struct Slot
  {
    double score;
    bool decoy;
    double* p;
    double* q;
    double* pep;
  };

  /// Statistics per sorted target score, the lookup basis.
  struct Table
  {
    std::vector<double> cutoffs, p, q, pep;
  };

  /// Nearest cutoff, ties in distance to the higher one, equal cutoffs to the first (PyProphet lookup).
  std::size_t nearest(const std::vector<double>& cutoffs, double score)
  {
    const auto i = static_cast<std::size_t>(std::lower_bound(cutoffs.begin(), cutoffs.end(), score) - cutoffs.begin());
    if (i < cutoffs.size() && cutoffs[i] == score) return i;
    std::size_t best;
    if (i == 0) best = 0;
    else if (i == cutoffs.size()) best = cutoffs.size() - 1;
    else best = score - cutoffs[i - 1] < cutoffs[i] - score ? i - 1 : i;
    while (best > 0 && cutoffs[best - 1] == cutoffs[best]) --best;
    return best;
  }

  /// Empty string on success, otherwise why the distribution cannot be estimated.
  std::string build(const std::vector<Slot>& slots, Table& table)
  {
    std::vector<double> decoys;
    for (const auto& slot : slots)
    {
      if (!std::isfinite(slot.score)) continue;
      (slot.decoy ? decoys : table.cutoffs).push_back(slot.score);
    }
    if (table.cutoffs.empty()) return "no target scores";
    if (decoys.empty()) return "no decoy scores";
    // The probit local-FDR density needs a bandwidth, so at least two targets.
    if (table.cutoffs.size() < 2) return "fewer than two target scores";
    std::sort(table.cutoffs.begin(), table.cutoffs.end());
    std::sort(decoys.begin(), decoys.end());
    try
    {
      table.p = MultipleTesting::pEmp(table.cutoffs, decoys);
      const double pi0 = MultipleTesting::pi0Est(table.p, PI0_LAMBDA, MultipleTesting::Pi0Method::Bootstrap).pi0;
      // PyProphet rejects this (e.g. all-tied scores) rather than reporting q and PEP of 0.
      if (!(pi0 > 0)) return "estimated pi0 <= 0";
      table.q = MultipleTesting::qValue(table.p, pi0, false);
      table.pep = MultipleTesting::lfdr(table.p, pi0, true, true, MultipleTesting::LfdrTransform::Probit, LFDR_ADJ, LFDR_EPS);
    }
    catch (const std::invalid_argument& error)
    {
      // MultipleTesting reports unestimable inputs this way, e.g. a zero-width density grid.
      return error.what();
    }
    for (std::size_t i = 0; i < table.cutoffs.size(); ++i)
      if (!std::isfinite(table.p[i]) || !std::isfinite(table.q[i]) || !std::isfinite(table.pep[i])) return "non-finite error statistics";
    return {};
  }

  void infer(const std::vector<Slot>& slots, const std::string& label, std::vector<std::string>& warnings)
  {
    if (slots.empty()) return;
    Table table;
    if (const auto reason = build(slots, table); !reason.empty())
    {
      warnings.push_back(label + ": inference unavailable: " + reason);
      return;
    }
    for (const auto& slot : slots)
    {
      if (!std::isfinite(slot.score)) continue;
      const auto i = nearest(table.cutoffs, slot.score);
      *slot.p = table.p[i];
      *slot.q = table.q[i];
      *slot.pep = table.pep[i];
    }
  }
} // namespace

void EntityInference::estimate(std::vector<QuantifiedEntity>& entities, std::size_t run_count, const std::string& level,
                               std::vector<std::string>& warnings)
{
  std::vector<Slot> global, experiment;
  std::vector<std::vector<Slot>> by_run(run_count);
  for (auto& entity : entities)
  {
    double best = -std::numeric_limits<double>::infinity();
    for (auto& run : entity.runs)
    {
      if (run.run >= run_count) throw std::out_of_range("Entity inference run index");
      by_run[run.run].push_back({run.score, entity.decoy, &run.p_value, &run.q_value, &run.pep});
      experiment.push_back({run.score, entity.decoy, &run.experiment_p_value, &run.experiment_q_value, &run.experiment_pep});
      best = std::max(best, run.score);
    }
    if (std::isfinite(best)) global.push_back({best, entity.decoy, &entity.global_p_value, &entity.global_q_value, &entity.global_pep});
  }
  infer(global, level + " global", warnings);
  infer(experiment, level + " experiment-wide", warnings);
  for (std::size_t r = 0; r < run_count; ++r)
    infer(by_run[r], level + " run-specific run " + std::to_string(r), warnings);
}
} // namespace OpenMS::OpenDIA
