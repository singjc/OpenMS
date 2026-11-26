// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIAGBM.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>

using namespace OpenMS::OpenDIA;

namespace
{
  int failures = 0;
  void check(bool ok, const std::string& what)
  {
    if (!ok)
    {
      ++failures;
      std::cerr << "FAIL: " << what << '\n';
    }
  }

  double auc(const std::vector<double>& s, const std::vector<std::uint8_t>& y)
  {
    std::vector<std::size_t> o(s.size());
    std::iota(o.begin(), o.end(), 0);
    std::sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return s[a] < s[b]; });
    double rank_sum = 0.0, np = 0.0;
    for (std::size_t i = 0; i < o.size(); ++i)
      if (y[o[i]])
      {
        rank_sum += static_cast<double>(i + 1);
        np += 1.0;
      }
    const double nn = static_cast<double>(s.size()) - np;
    return (rank_sum - np * (np + 1.0) / 2.0) / (np * nn);
  }

  struct Problem
  {
    std::size_t d = 0;
    std::vector<float> x;
    std::vector<std::uint8_t> y;
    std::vector<double> w;
  };

  // Nonlinear separable problem: y = 1 iff x0 * x1 > 0 or x2 > 1.5, plus noise features.
  Problem makeProblem(std::size_t n, std::size_t d, unsigned seed)
  {
    Problem p;
    p.d = d;
    p.x.resize(n * d);
    p.y.resize(n);
    p.w.assign(n, 1.0);
    std::mt19937_64 rng(seed);
    std::normal_distribution<float> nd(0.0f, 1.0f);
    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j < d; ++j) p.x[i * d + j] = nd(rng);
      p.y[i] = (p.x[i * d] * p.x[i * d + 1] > 0.0f) || p.x[i * d + 2] > 1.5f;
    }
    return p;
  }
}

int main()
{
  // (a) Separable synthetic problem.
  {
    const std::size_t n = 20000;
    auto p = makeProblem(n, 8, 1);
    std::vector<std::uint32_t> tr, va, te;
    for (std::uint32_t i = 0; i < n; ++i) (i % 10 < 7 ? tr : i % 10 < 8 ? va : te).push_back(i);
    HistGBM m;
    HistGBM::Options o;
    o.threads = 2;
    m.fit(p.x.data(), p.d, tr, va, p.y, p.w, o);
    std::vector<double> s;
    std::vector<std::uint8_t> y;
    for (auto i : te)
    {
      s.push_back(m.predict(p.x.data() + i * p.d));
      y.push_back(p.y[i]);
    }
    const double a = auc(s, y);
    std::cout << "separable AUC " << a << " trees " << m.trees.size() << '\n';
    check(a > 0.99, "separable AUC > 0.99");
  }

  // (b) Thread-count independence.
  {
    const std::size_t n = 30000;
    auto p = makeProblem(n, 12, 2);
    for (std::size_t i = 0; i < n; ++i) p.w[i] = p.y[i] ? 2.5 : 1.0;
    std::vector<std::uint32_t> tr, va;
    for (std::uint32_t i = 0; i < n; ++i) (i % 10 ? tr : va).push_back(i);
    HistGBM a, b;
    HistGBM::Options o;
    o.max_iter = 60;
    o.threads = 1;
    a.fit(p.x.data(), p.d, tr, va, p.y, p.w, o);
    o.threads = 4;
    b.fit(p.x.data(), p.d, tr, va, p.y, p.w, o);
    bool same = a.trees.size() == b.trees.size() && a.init == b.init;
    for (std::size_t i = 0; same && i < n; ++i) same = a.predict(p.x.data() + i * p.d) == b.predict(p.x.data() + i * p.d);
    check(same, "identical model for threads 1 and 4");
  }

  // (c) End-to-end gbmRescore on synthetic windows.
  {
    CrossFitScorer::Input in;
    const std::size_t pairs = 3000, wpp = 4, d = 6;
    in.dimension = d;
    std::mt19937_64 rng(7);
    std::normal_distribution<float> nd(0.0f, 1.0f);
    std::vector<float> weak;
    for (std::uint32_t pr = 0; pr < pairs; ++pr)
      for (int dec = 0; dec < 2; ++dec)
        for (std::size_t w = 0; w < wpp; ++w)
        {
          const bool signal = !dec && w == 0 && pr % 2 == 0;  // half the targets are present
          for (std::size_t j = 0; j < d; ++j) in.features.push_back(nd(rng));
          if (signal)
          {
            in.features[in.features.size() - d] += 3.0f;
            in.features[in.features.size() - d + 1] += 2.0f;
          }
          in.precursor.push_back(2 * pr + dec);
          in.pair.push_back(pr);
          in.decoy.push_back(static_cast<std::uint8_t>(dec));
          in.fold.push_back(static_cast<std::uint8_t>(pr % 3));
          in.pair_hash.push_back(CrossFitScorer::fnv1a("p" + std::to_string(pr)));
          // Weakly informative incoming score: a noisy view of feature 0.
          weak.push_back(in.features[in.features.size() - d] + 0.7f * nd(rng));
        }
    const auto count = [&](const std::vector<double>& q)
    {
      std::size_t c = 0;
      for (std::size_t r = 0; r < q.size(); ++r) c += !in.decoy[r] && q[r] <= 0.01;
      return c;
    };
    std::vector<std::uint32_t> all(in.rows());
    std::iota(all.begin(), all.end(), 0u);
    const auto comp = CrossFitScorer::pooledPairedTDC(in, weak, all);
    std::vector<double> q0(in.rows(), 1.0);
    for (std::size_t i = 0; i < comp.ranked.size(); ++i) q0[comp.ranked[i]] = comp.q[i];
    std::vector<float> flat(in.rows());
    for (auto& v : flat) v = nd(rng);
    const auto compf = CrossFitScorer::pooledPairedTDC(in, flat, all);
    std::vector<double> qf(in.rows(), 1.0);
    for (std::size_t i = 0; i < compf.ranked.size(); ++i) qf[compf.ranked[i]] = compf.q[i];
    for (double frac : {0.02, 0.0})
    {
      const auto res = gbmRescore(in, weak, frac, 1, 2, [](const std::string& l) { std::cout << l << '\n'; });
      const std::size_t c = count(res.q_value);
      std::cout << "gbmRescore frac " << frac << ": targets at q<=0.01 " << c << " (weak input " << count(q0) << ", uninformative "
                << count(qf) << ")\n";
      check(c > count(q0) && c > count(qf), "gbmRescore gains targets");
      std::size_t kept = 0;
      for (float v : res.score) kept += v > -1e29f;
      check(kept == 2 * pairs, "one kept window per precursor");
    }

    // (c2) Fold-local training (FDR-20261126): still gains targets, and the held-out fold's scores do not depend on the
    // incoming scores of that fold (a perturbation confined to fold 0 leaves fold 0's raw GBM margins unchanged).
    {
      GBMVariant fl;
      fl.fold_local = true;
      // A stronger incoming score (feature 0 itself) so every outer fold's training rows yield positives.
      std::vector<float> strong(in.rows());
      for (std::size_t r = 0; r < in.rows(); ++r) strong[r] = in.features[r * d];
      const auto res = gbmRescore(in, strong, 0.02, 1, 2, [](const std::string& l) { std::cout << l << '\n'; }, 3, 0.01, fl);
      const std::size_t c = count(res.q_value);
      std::cout << "gbmRescore fold-local: targets at q<=0.01 " << c << " (weak input " << count(q0) << ")\n";
      check(c > count(q0), "fold-local gbmRescore gains targets");
      std::size_t kept = 0;
      for (float v : res.score) kept += v > -1e29f;
      check(kept == 2 * pairs, "fold-local: one kept window per precursor");
      auto shifted = strong;
      for (std::size_t r = 0; r < in.rows(); ++r)
        if (in.fold[r] == 0) shifted[r] += in.decoy[r] ? 0.3f : -0.3f;
      // Raw margins are compared through the kept-window choice and the calibrated score order within fold 0. The
      // calibration of fold 0 uses fold 0's own decoys only, so identical margins give identical calibrated scores.
      const auto res2 = gbmRescore(in, shifted, 0.02, 1, 2, {}, 3, 0.01, fl);
      bool same0 = true;
      for (std::size_t r = 0; r < in.rows(); ++r)
        if (in.fold[r] == 0) same0 = same0 && res.score[r] == res2.score[r];
      check(same0, "fold-local: fold 0 scores are independent of fold 0's incoming scores");
      // Default (pooled) training does depend on them; this documents the leak the option removes.
      const auto p1 = gbmRescore(in, strong, 0.02, 1, 2);
      const auto p2 = gbmRescore(in, shifted, 0.02, 1, 2);
      bool pooled_same0 = true;
      for (std::size_t r = 0; r < in.rows(); ++r)
        if (in.fold[r] == 0) pooled_same0 = pooled_same0 && p1.score[r] == p2.score[r];
      check(!pooled_same0, "pooled training: fold 0 scores depend on fold 0's incoming scores");
      bool threw = false;
      try
      {
        auto two = in;
        for (auto& f : two.fold) f = static_cast<std::uint8_t>(f % 2);
        gbmRescore(two, weak, 0.02, 1, 2, {}, 3, 0.01, fl);
      }
      catch (const std::invalid_argument&) { threw = true; }
      check(threw, "fold-local needs at least 3 folds");
    }

    // (d) blendScores: kept rows follow the GBM, the order at w = 1 follows the network, bad w throws.
    {
      const auto net = CrossFitScorer::calibrateTail(in, [&]
      {
        // Network-like input: the best window per precursor keeps the weak score.
        std::vector<float> k(in.rows(), CrossFitScorer::NOT_KEPT);
        for (std::size_t r = 0; r < in.rows(); r += wpp) k[r] = weak[r];
        return k;
      }());
      const auto g = gbmRescore(in, weak, 0.02, 1, 2);
      const auto b0 = blendScores(in, net, g.score, 0.0);
      const auto b5 = blendScores(in, net, g.score, 0.5);
      const auto b1 = blendScores(in, net, g.score, 1.0);
      bool rows_ok = true;
      for (std::size_t r = 0; r < in.rows(); ++r)
        rows_ok = rows_ok && ((b5.score[r] > -1e29f) == (g.score[r] > -1e29f)) && ((b1.score[r] > -1e29f) == (g.score[r] > -1e29f));
      check(rows_ok, "blend keeps exactly the GBM-kept rows");
      // w = 0: same order of precursors as the GBM score (calibrateTail is monotone per fold).
      bool order0 = true;
      for (std::size_t a = 0; a < in.rows() && order0; a += 97)
        for (std::size_t c = 0; c < in.rows() && order0; c += 89)
          if (g.score[a] > -1e29f && g.score[c] > -1e29f && in.fold[a] == in.fold[c] && g.score[a] < g.score[c]) order0 = b0.score[a] <= b0.score[c];
      check(order0, "blend w = 0 keeps the GBM order within a fold");
      // w = 1: within a fold, precursor order follows the network score.
      std::vector<float> net_p(2 * pairs, 0.0f);
      for (std::size_t r = 0; r < in.rows(); ++r)
        if (net[r] > -1e29f) net_p[in.precursor[r]] = net[r];
      bool order1 = true;
      for (std::size_t a = 0; a < in.rows() && order1; ++a)
        for (std::size_t c = a + 1; c < std::min(in.rows(), a + 400) && order1; ++c)
          if (b1.score[a] > -1e29f && b1.score[c] > -1e29f && in.fold[a] == in.fold[c] && net_p[in.precursor[a]] < net_p[in.precursor[c]])
            order1 = b1.score[a] <= b1.score[c];
      check(order1, "blend w = 1 follows the network order within a fold");
      const auto b5b = blendScores(in, net, g.score, 0.5);
      check(b5b.score == b5.score && b5b.q_value == b5.q_value, "blend is deterministic");
      std::cout << "blendScores w 0 / 0.5 / 1: targets at q<=0.01 " << count(b0.q_value) << " / " << count(b5.q_value) << " / " << count(b1.q_value) << '\n';
      bool threw = false;
      try { blendScores(in, net, g.score, 1.5); } catch (const std::invalid_argument&) { threw = true; }
      check(threw, "blend rejects w outside [0, 1]");
    }
  }

  // (f) Feature restrictions (FDR-20261125): all-ones use_feature and all-zero monotone reproduce the plain fit,
  // excluded features are never split on, and a monotone feature gives a monotone prediction.
  {
    const std::size_t n = 20000, d = 6;
    auto p = makeProblem(n, d, 11);
    std::vector<std::uint32_t> tr, va;
    for (std::uint32_t i = 0; i < n; ++i) (i % 10 ? tr : va).push_back(i);
    HistGBM::Options o;
    o.threads = 2;
    o.max_iter = 40;
    HistGBM plain, ones;
    plain.fit(p.x.data(), d, tr, va, p.y, p.w, o);
    ones.use_feature.assign(d, 1);
    ones.fit(p.x.data(), d, tr, va, p.y, p.w, o);
    bool same = plain.trees.size() == ones.trees.size();
    for (std::size_t i = 0; same && i < n; ++i) same = plain.predict(p.x.data() + i * d) == ones.predict(p.x.data() + i * d);
    check(same, "use_feature all ones reproduces the plain fit");

    HistGBM ex;
    ex.use_feature.assign(d, 1);
    ex.use_feature[2] = 0;
    ex.fit(p.x.data(), d, tr, va, p.y, p.w, o);
    bool none = !ex.trees.empty();
    for (const auto& t : ex.trees)
      for (const auto& nd : t) none = none && nd.feature != 2;
    check(none, "excluded feature is never split on");

    // x2 > 1.5 raises the label, so +1 is the natural direction; -1 forces the reverse. Check both on a sweep of x2.
    for (int dir : {+1, -1})
    {
      HistGBM mono;
      mono.monotone.assign(d, 0);
      mono.monotone[2] = static_cast<std::int8_t>(dir);
      mono.fit(p.x.data(), d, tr, va, p.y, p.w, o);
      bool ok = !mono.trees.empty();
      for (std::size_t i = 0; ok && i < 200; ++i)
      {
        std::vector<float> row(p.x.begin() + i * d, p.x.begin() + (i + 1) * d);
        double prev = 0.0;
        for (int k = 0; k <= 60 && ok; ++k)
        {
          row[2] = -3.0f + 0.1f * k;
          const double v = mono.predict(row.data());
          if (k) ok = dir > 0 ? v >= prev : v <= prev;
          prev = v;
        }
      }
      check(ok, dir > 0 ? "monotone +1 prediction is non-decreasing" : "monotone -1 prediction is non-increasing");
    }
    HistGBM zero;
    zero.monotone.assign(d, 0);
    zero.fit(p.x.data(), d, tr, va, p.y, p.w, o);
    std::vector<double> s;
    for (std::size_t i = 0; i < n; i += 10) s.push_back(zero.predict(p.x.data() + i * d));
    std::vector<std::uint8_t> y;
    for (std::size_t i = 0; i < n; i += 10) y.push_back(p.y[i]);
    check(auc(s, y) > 0.97, "all-free monotone vector still fits");
    bool threw = false;
    try
    {
      HistGBM bad;
      bad.use_feature.assign(d - 1, 1);
      bad.fit(p.x.data(), d, tr, va, p.y, p.w, o);
    }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "wrong-size use_feature throws");
  }

  // Optional timing run: OPENDIA_GBM_BENCH_ROWS=<rows> (49 features, OPENDIA_GBM_BENCH_THREADS threads, default 8).
  if (const char* env = std::getenv("OPENDIA_GBM_BENCH_ROWS"))
  {
    const std::size_t n = std::strtoull(env, nullptr, 10);
    const char* te = std::getenv("OPENDIA_GBM_BENCH_THREADS");
    auto p = makeProblem(n, 49, 3);
    std::mt19937_64 rng(9);
    std::normal_distribution<float> nd(0.0f, 1.0f);
    for (auto& v : p.x) v += 0.8f * nd(rng);  // label noise via feature noise, so boosting runs long
    std::vector<std::uint32_t> tr, va;
    for (std::uint32_t i = 0; i < n; ++i) (i % 10 ? tr : va).push_back(i);
    HistGBM m;
    HistGBM::Options o;
    o.threads = te ? std::atoi(te) : 8;
    const auto t0 = std::chrono::steady_clock::now();
    m.fit(p.x.data(), p.d, tr, va, p.y, p.w, o);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "bench " << tr.size() << " x 49 threads " << o.threads << ": " << m.trees.size() << " trees in " << sec << " s\n";
  }

  if (failures) std::cerr << failures << " failure(s)\n";
  else std::cout << "all OpenDIAGBM tests passed\n";
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
