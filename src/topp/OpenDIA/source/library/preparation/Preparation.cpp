// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $

#include <library/preparation/Preparation.h>
#include <library/preparation/DecoySequences.h>
#include <library/preparation/AssayPrediction.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/config.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <tuple>
#include <utility>

#ifdef WITH_ONNX
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#endif

namespace OpenMS::OpenDIA::Preparation
{
#ifdef WITH_ONNX
  namespace
  {
    // Keep the arithmetic, iteration order and sampling of the original predicted-decoy path.
    using FragmentKey = std::tuple<std::string, int, int>;
    using TargetFragments = std::map<std::string, std::vector<const Transition*>>;

    FragmentKey fragmentKey(const Transition& transition)
    {
      return {transition.fragment_type, transition.fragment_nr, transition.fragment_charge};
    }

    std::optional<std::pair<double, double>> fitLine(const std::vector<double>& x, const std::vector<double>& y)
    {
      const Size n = x.size();
      if (n < 2) return std::nullopt;
      double mx = 0.0, my = 0.0;
      for (Size i = 0; i < n; ++i)
      {
        mx += x[i];
        my += y[i];
      }
      mx /= n;
      my /= n;
      double sxx = 0.0, sxy = 0.0;
      for (Size i = 0; i < n; ++i)
      {
        sxx += (x[i] - mx) * (x[i] - mx);
        sxy += (x[i] - mx) * (y[i] - my);
      }
      if (sxx <= 0.0) return std::nullopt;
      const double b = sxy / sxx;
      return std::make_pair(my - b * mx, b);
    }

    AssayPredictor::Config predictorConfig(const PredictionOptions& options, bool predict_ccs)
    {
      AssayPredictor::Config config;
      if (!options.model_directory.empty())
      {
        const std::string dir = options.model_directory + "/";
        config.rt_model_path = dir + "peptdeep_rt_dynamic.onnx";
        config.ms2_model_path = dir + "peptdeep_ms2_dynamic.onnx";
        config.ccs_model_path = dir + "peptdeep_ccs_dynamic.onnx";
      }
      config.intra_op_threads = std::max(1, options.threads);
      config.predict_ccs = predict_ccs;
      return config;
    }

    std::optional<PredictionInput> predictionInput(const Compound& compound, double nce, Int64 instrument_index)
    {
      PredictionInput precursor;
      try
      {
        precursor.peptide = AASequence::fromString(compound.sequence);
      }
      catch (const Exception::BaseException&)
      {
        return std::nullopt;
      }
      if (precursor.peptide.size() < 2 || compound.charge <= 0) return std::nullopt;
      precursor.id = compound.id;
      precursor.charge = compound.charge;
      precursor.nce = static_cast<float>(nce);
      precursor.instrument_index = instrument_index;
      return precursor;
    }

    Assays predictTolerant(AssayPredictor& predictor, const std::vector<PredictionInput>& batch, Size& rejected)
    {
      try
      {
        return predictor.predict(batch);
      }
      catch (const Exception::BaseException&)
      {
        Assays merged;
        for (const auto& precursor : batch)
        {
          try
          {
            auto one = predictor.predict({precursor});
            merged.compounds.insert(merged.compounds.end(), one.compounds.begin(), one.compounds.end());
            merged.transitions.insert(merged.transitions.end(), one.transitions.begin(), one.transitions.end());
          }
          catch (const Exception::BaseException&)
          {
            ++rejected;
          }
        }
        return merged;
      }
    }

    std::pair<double, std::string> predictionSetting(const Assays& targets, const TargetFragments& by_compound,
                                                    const PredictionOptions& options)
    {
      std::vector<const Compound*> annotated;
      for (const auto& compound : targets.compounds)
      {
        const auto it = by_compound.find(compound.id);
        if (!compound.isPeptide() || it == by_compound.end() || it->second.size() < 3) continue;
        bool by_only = true;
        for (const auto* transition : it->second)
        {
          by_only = by_only && transition->fragment_nr > 0 &&
            (transition->fragment_charge == 1 || transition->fragment_charge == 2) &&
            (transition->fragment_type == "b" || transition->fragment_type == "y");
        }
        if (by_only) annotated.push_back(&compound);
      }
      const Size sample_size = std::min<Size>(annotated.size(), 400);
      if (sample_size < 20)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Decoy prediction needs at least 20 targets with annotated b/y transitions to estimate the "
          "collision energy and instrument; set them explicitly instead.");
      }
      std::vector<const Compound*> sample;
      for (Size i = 0; i < sample_size; ++i) sample.push_back(annotated[i * annotated.size() / sample_size]);

      AssayPredictor predictor(predictorConfig(options, false));
      double best_cosine = -1.0;
      std::pair<double, std::string> best{20.0, std::string(ML::ALPHAPEPTDEEP_INSTRUMENTS.front())};
      for (const auto instrument_name : ML::ALPHAPEPTDEEP_INSTRUMENTS)
      {
        const std::string instrument(instrument_name);
        for (int nce = 20; nce <= 40; ++nce)
        {
          std::vector<PredictionInput> batch;
          for (const auto* compound : sample)
          {
            if (auto precursor = predictionInput(*compound, nce, ML::instrumentIndex(instrument)))
              batch.push_back(std::move(*precursor));
          }
          Size rejected = 0;
          const auto predicted = predictTolerant(predictor, batch, rejected);
          std::map<std::string, std::map<FragmentKey, double>> predicted_by_compound;
          for (const auto& transition : predicted.transitions)
            predicted_by_compound[transition.peptide_ref][fragmentKey(transition)] = transition.library_intensity;
          std::vector<double> cosines;
          for (const auto& [id, intensities] : predicted_by_compound)
          {
            double dot = 0.0, nl = 0.0, np = 0.0;
            for (const auto* transition : by_compound.at(id))
            {
              const auto it = intensities.find(fragmentKey(*transition));
              const double p = it == intensities.end() ? 0.0 : it->second;
              dot += transition->library_intensity * p;
              nl += transition->library_intensity * transition->library_intensity;
              np += p * p;
            }
            if (nl > 0.0 && np > 0.0) cosines.push_back(dot / std::sqrt(nl * np));
          }
          if (cosines.empty()) continue;
          double mean = 0.0;
          for (const double cosine : cosines) mean += cosine;
          mean /= cosines.size();
          if (mean > best_cosine)
          {
            best_cosine = mean;
            best = {nce, instrument};
          }
        }
      }
      if (options.log_type != ProgressLogger::NONE)
      {
        OPENMS_LOG_INFO << "Decoy prediction: instrument " << best.second << ", NCE " << best.first << " (mean cosine "
                        << best_cosine << " to the library intensities of " << sample_size << " targets)." << std::endl;
      }
      return best;
    }

    void predictDecoyAssays(const Assays& targets, Assays& decoys, const PredictionOptions& options)
    {
      TargetFragments by_compound;
      double product_lo = std::numeric_limits<double>::max();
      double product_hi = std::numeric_limits<double>::lowest();
      for (const auto& transition : targets.transitions)
      {
        if (transition.decoy || !transition.detecting) continue;
        by_compound[transition.peptide_ref].push_back(&transition);
        product_lo = std::min(product_lo, transition.product_mz);
        product_hi = std::max(product_hi, transition.product_mz);
      }
      std::map<std::string, const Compound*> target_of;
      for (const auto& compound : targets.compounds) target_of[compound.id] = &compound;
      const auto [nce, instrument] = predictionSetting(targets, by_compound, options);
      const Int64 instrument_index = ML::instrumentIndex(instrument);

      Size targets_with_im = 0;
      std::vector<const Compound*> usable;
      for (const auto& compound : targets.compounds)
      {
        targets_with_im += compound.drift_time > 0.0;
        if (compound.isPeptide() && std::isfinite(compound.rt)) usable.push_back(&compound);
      }
      const bool map_im = targets_with_im * 2 > targets.compounds.size();
      AssayPredictor predictor(predictorConfig(options, map_im));
      std::optional<std::pair<double, double>> rt_line, im_line;
      const Size fit_size = std::min<Size>(usable.size(), 2000);
      if (fit_size >= 20)
      {
        std::vector<PredictionInput> batch;
        for (Size i = 0; i < fit_size; ++i)
        {
          if (auto precursor = predictionInput(*usable[i * usable.size() / fit_size], nce, instrument_index))
            batch.push_back(std::move(*precursor));
        }
        Size rejected = 0;
        const auto predicted = predictTolerant(predictor, batch, rejected);
        std::vector<double> rt_x, rt_y, im_x, im_y;
        for (const auto& compound : predicted.compounds)
        {
          const auto* target = target_of.at(compound.id);
          rt_x.push_back(compound.rt);
          rt_y.push_back(target->rt);
          if (map_im && target->drift_time > 0.0 && compound.drift_time > 0.0)
          {
            im_x.push_back(compound.drift_time);
            im_y.push_back(target->drift_time);
          }
        }
        rt_line = fitLine(rt_x, rt_y);
        if (im_x.size() >= 20) im_line = fitLine(im_x, im_y);
      }
      if (!rt_line)
        OPENMS_LOG_WARN << "Decoy prediction: too few targets with a library RT to map predicted RT; decoys keep their target's RT." << std::endl;
      else if (options.log_type != ProgressLogger::NONE)
        OPENMS_LOG_INFO << "Decoy prediction: library RT = " << rt_line->first << " + " << rt_line->second << " x predicted RT." << std::endl;
      if (map_im && !im_line)
        OPENMS_LOG_WARN << "Decoy prediction: too few targets with ion mobility to map predicted 1/K0; decoys keep their target's." << std::endl;

      std::map<std::string, Compound*> decoy_of;
      for (auto& compound : decoys.compounds) decoy_of[compound.id] = &compound;
      std::vector<Transition> transitions;
      std::set<std::string> kept;
      Size rejected = 0, too_few = 0, short_kept = 0;
      ProgressLogger progress;
      progress.setLogType(options.log_type);
      progress.startProgress(0, decoys.compounds.size(), "Predicting decoy assays");
      constexpr Size BATCH = 5000;
      for (Size start = 0; start < decoys.compounds.size(); start += BATCH)
      {
        progress.setProgress(start);
        std::vector<PredictionInput> batch;
        for (Size i = start; i < std::min(start + BATCH, decoys.compounds.size()); ++i)
        {
          if (!decoys.compounds[i].isPeptide()) continue;
          if (auto precursor = predictionInput(decoys.compounds[i], nce, instrument_index)) batch.push_back(std::move(*precursor));
          else ++rejected;
        }
        if (batch.empty()) continue;
        auto predicted = predictTolerant(predictor, batch, rejected);
        std::map<std::string, std::vector<Transition*>> fragments;
        for (auto& transition : predicted.transitions) fragments[transition.peptide_ref].push_back(&transition);
        for (const auto& compound : predicted.compounds)
        {
          Compound* decoy = decoy_of.at(compound.id);
          const auto target_it = by_compound.find(compound.id.substr(6)); // DECOY_
          if (target_it == by_compound.end()) continue;
          const Size k = target_it->second.size();
          std::vector<Transition*> candidates;
          for (auto* transition : fragments[compound.id])
          {
            if (transition->library_intensity > 0.0 &&
                (transition->fragment_type == "b" || transition->fragment_type == "y") &&
                transition->fragment_charge >= 1 && transition->fragment_charge <= 4 &&
                transition->product_mz >= product_lo && transition->product_mz <= product_hi)
              candidates.push_back(transition);
          }
          // OpenDIA's historical predict_keep_short=true policy, including rejection of empty assays.
          if (candidates.empty())
          {
            ++too_few;
            continue;
          }
          short_kept += candidates.size() < k;
          std::stable_sort(candidates.begin(), candidates.end(),
                          [](const auto* a, const auto* b) { return a->library_intensity > b->library_intensity; });
          candidates.resize(std::min(k, candidates.size()));
          const double top = candidates.front()->library_intensity;
          if (rt_line) decoy->rt = rt_line->first + rt_line->second * compound.rt;
          if (im_line) decoy->drift_time = im_line->first + im_line->second * compound.drift_time;
          const double precursor_mz = candidates.front()->precursor_mz;
          for (const auto* candidate : candidates)
          {
            Transition transition = *candidate;
            transition.library_intensity = candidate->library_intensity / top;
            transition.precursor_mz = precursor_mz;
            transition.precursor_im = decoy->drift_time;
            transition.decoy = true;
            transition.detecting = true;
            transition.quantifying = true;
            transition.identifying = false;
            transitions.push_back(std::move(transition));
          }
          kept.insert(compound.id);
        }
      }
      progress.endProgress();
      decoys.compounds.erase(std::remove_if(decoys.compounds.begin(), decoys.compounds.end(),
        [&kept](const Compound& c) { return !kept.contains(c.id); }), decoys.compounds.end());
      std::set<std::string> proteins;
      for (const auto& compound : decoys.compounds) proteins.insert(compound.protein_refs.begin(), compound.protein_refs.end());
      decoys.proteins.erase(std::remove_if(decoys.proteins.begin(), decoys.proteins.end(),
        [&proteins](const Protein& p) { return !proteins.contains(p.id); }), decoys.proteins.end());
      decoys.transitions = std::move(transitions);
      if (options.log_type != ProgressLogger::NONE)
      {
        OPENMS_LOG_INFO << "Decoy prediction: " << kept.size() << " decoys predicted; " << too_few << " dropped with fewer fragments than "
                        << "their target inside " << product_lo << "-" << product_hi << " m/z; " << short_kept << " kept with fewer fragments; "
                        << rejected << " rejected by the models." << std::endl;
      }
    }
  }
#endif

  Assays predictDecoys(Assays targets, const PredictionOptions& options)
  {
    auto decoys = decoySequences(targets, options.log_type);
#ifdef WITH_ONNX
    predictDecoyAssays(targets, decoys, options);
    if (targets.compounds.empty() || targets.proteins.empty())
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "The input experiment has no compounds or proteins.");
    if (static_cast<double>(decoys.compounds.size()) / targets.compounds.size() < 0.8 ||
        static_cast<double>(decoys.proteins.size()) / targets.proteins.size() < 0.8)
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "The number of decoys for compounds/peptides or proteins is below the threshold of 80% of the number of targets.");
    targets.transitions.insert(targets.transitions.end(), std::make_move_iterator(decoys.transitions.begin()),
                                std::make_move_iterator(decoys.transitions.end()));
    targets.compounds.insert(targets.compounds.end(), std::make_move_iterator(decoys.compounds.begin()),
                             std::make_move_iterator(decoys.compounds.end()));
    targets.proteins.insert(targets.proteins.end(), std::make_move_iterator(decoys.proteins.begin()),
                            std::make_move_iterator(decoys.proteins.end()));
    return targets;
#else
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Predicted decoys require an OpenMS build configured with WITH_ONNX=ON.");
#endif
  }
}
