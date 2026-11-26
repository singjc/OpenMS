// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenDIAScorer.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief Deterministic histogram gradient boosting for binary logistic loss.

    Features are quantile binned per feature on the training rows (up to max_bins bins, uint8 codes). Trees grow
    leaf-wise up to max_leaf_nodes leaves; a split needs min_samples_leaf rows on each side. Split gain and leaf values
    use gradient and hessian sums with l2 regularisation (sums guarded by a small epsilon). The initial score is the
    weighted log-odds. Boosting stops early when the weighted validation logloss has not improved for `patience`
    rounds (all fitted trees are kept, as in sklearn). Histograms are built in parallel over features; every reduction
    has a fixed order, so the model does not depend on the thread count.
  */
  class HistGBM
  {
  public:
    struct Options
    {
      int max_iter = 300;
      double learning_rate = 0.05;
      int max_leaf_nodes = 31;
      int max_bins = 255;
      std::size_t min_samples_leaf = 20;
      double l2 = 0.0;
      int patience = 10;
      double tolerance = 1e-7;
      int threads = 1;
    };

    struct Node
    {
      int feature = -1;     ///< -1 for a leaf.
      float threshold = 0;  ///< Go left if x < threshold (NaN goes right).
      int left = -1, right = -1;
      double value = 0.0;   ///< Leaf output (shrinkage applied).
    };

    double init = 0.0;
    std::vector<std::vector<Node>> trees;

    /**
      @brief Optional restrictions (FDR-20261125). Empty vectors mean none and leave the fit unchanged.

      use_feature[j] == 0: the model never splits on feature j.
      monotone[j] = +1 / -1: the model's output is non-decreasing / non-increasing in feature j. Each leaf carries bounds
      on its weight; a split on a constrained feature must order the children's weights, and its children's bounds meet at
      the midpoint of their weights. Splits on other features pass the bounds down (the basic method of LightGBM). With a
      monotone constraint the split gain is computed from the bounded weights.
    */
    std::vector<std::uint8_t> use_feature;
    std::vector<std::int8_t> monotone;

    /**
      @brief Fit on rows `train` of the row-major matrix `x` (n x d); rows `valid` drive early stopping (may be empty).
      @param label 1 positive, 0 negative, indexed like `x` rows.
      @param weight Sample weight, indexed like `x` rows.
    */
    void fit(const float* x, std::size_t d, const std::vector<std::uint32_t>& train, const std::vector<std::uint32_t>& valid,
             const std::vector<std::uint8_t>& label, const std::vector<double>& weight, const Options& options);

    /// Raw margin (log-odds) of one row.
    double predict(const float* row) const;
  };

  /// GBM feature restrictions for gbmRescore (FDR-20261125; see HistGBM::use_feature and HistGBM::monotone). Empty = none.
  struct GBMVariant
  {
    std::vector<std::uint8_t> use_feature;  ///< Per feature column: 0 excludes it from the GBM.
    std::vector<std::int8_t> monotone;      ///< Per feature column: +1 non-decreasing, -1 non-increasing, 0 free.
    /**
      Fold-local training (FDR-20261126, opt-in; false = unchanged). Each outer fold f gets its own training chain on the
      rows outside f: positive selection, hard-negative mining and every round's training scores use only those rows, and
      the training rows are rescored between rounds by inner cross-fitting (each other fold scored by a model trained on
      the remaining one). Fold f is only ever predicted, never selected from or trained on. The initial network score is
      still the pooled out-of-fold network score. Needs at least 3 folds.
    */
    bool fold_local = false;
  };

  /**
    @brief Hard-negative GBM rescoring of CrossFitScorer windows.

    Up to `rounds` rounds (default 3): positives are the target winners at pooled paired TDC q <= `positive_q` (default 0.01) under the current score, negatives
    all rows of the top `hard_fraction` decoy precursors ranked by their winner score (0: all decoy winners). Per fold
    a HistGBM is trained on the other folds (negatives subsampled to at most 1.5 M, positive weight n_neg / n_pos,
    validation rows by pair_hash % 100 < 10) and scores its fold by the raw margin. Finally each precursor keeps its
    best window, scores are calibrated with CrossFitScorer::calibrateTail and q-values come from pooledPairedTDC.
    Only score and q_value of the result are filled. @p variant restricts the GBM's use of feature columns (default: none).
  */
  CrossFitScorer::Result gbmRescore(const CrossFitScorer::Input& input, const std::vector<float>& score, double hard_fraction,
                                    std::uint64_t seed, int threads, const std::function<void(const std::string&)>& log = {},
                                    int rounds = 3, double positive_q = 0.01, const GBMVariant& variant = {});

  /**
    @brief Blend the network and GBM calibrated scores per precursor, then recalibrate and recompete (FDR-20261125).

    For each precursor, b = w * network + (1 - w) * gbm. Here network is the calibrated score of the precursor's network-kept
    window, and gbm is the calibrated score of its GBM-kept window. b is placed on the GBM-kept row, and the other rows
    stay NOT_KEPT. Scores are then recalibrated with CrossFitScorer::calibrateTail, and q-values come from pooledPairedTDC.
    A precursor without a network-kept window keeps its GBM score as b. w must lie in [0, 1].
    Only score and q_value of the result are filled.
  */
  CrossFitScorer::Result blendScores(const CrossFitScorer::Input& input, const std::vector<float>& network, const std::vector<float>& gbm,
                                     double w);

  /**
    @brief Label-free bulk-null weighted target-decoy q-values (FDR-20261126, opt-in; CALIB section 9 ported to OpenDIA).

    The ranking is the pooled paired TDC of @p score (kept rows, as pooledPairedTDC) and is not changed; only the decoy count
    in q is reweighted. The bulk is every precursor winner whose score lies below its fold's median decoy-winner score, so it
    is almost all false targets and decoys. Per fold f a HistGBM (default options) learns target versus decoy on the bulk
    winners of the other folds (feature columns only, unit weights, validation rows by pair_hash % 100 < 10). Every decoy
    winner of f gets w = exp(margin) * n_decoy / n_target of that model (the false-target / decoy density ratio), clipped to
    [0.2, 5] and divided by the mean w of f's bulk decoys. Along the ranked list q = min(1, (sum of w over decoys + 1) /
    targets) with a reverse cumulative minimum. Returns q per row (1 for rows outside the ranked list).
    @throws std::invalid_argument With fewer than 2 folds or a fold without bulk targets or bulk decoys.
  */
  std::vector<double> bulkNullQ(const CrossFitScorer::Input& input, const std::vector<float>& score, int threads,
                                const std::function<void(const std::string&)>& log = {});
}
