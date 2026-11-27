// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#pragma once

#include <results/Quantification.h>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief OpenDIA-owned peptide/protein error statistics in three contexts.

    Method of PyProphet and Rosenberger et al. 2017 (10.1038/nmeth.4398), defaults of PyProphet's
    error_statistics: empirical p-values against all decoy scores, bootstrap pi0 over
    lambda = 0.10..0.45 step 0.05, q-values, and probit local FDR (adj 1.5, eps 1e-8,
    truncated, monotone) as PEP. Numerical steps use OpenMS Math::MultipleTesting.

    Each context is one target/decoy distribution:
    - run-specific: the best score of each entity in one run, separately for every run;
    - experiment-wide: those entity/run scores, pooled over runs;
    - global: the best score of each entity over all runs.
    Every score, target or decoy, gets the statistics of the nearest target score.

    A distribution that cannot be estimated (no targets or decoys, fewer than two targets,
    pi0 <= 0 as rejected by PyProphet, or non-finite results) keeps its NaN fields and adds one
    warning. Other runs and contexts are unaffected. There is no fallback estimator.
  */
  class EntityInference
  {
  public:
    /// Fills p/q/PEP of all contexts in place from entity.runs[].score. level names warnings.
    static void estimate(std::vector<QuantifiedEntity>& entities, std::size_t run_count, const std::string& level,
                         std::vector<std::string>& warnings);
  };
}
