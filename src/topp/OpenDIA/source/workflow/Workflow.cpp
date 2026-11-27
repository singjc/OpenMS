// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <workflow/Workflow.h>

#include <extraction/FragmentGrid.h>
#include <scoring/GBM.h>
#include <features/Interference.h>
#include <calibration/RTRefine.h>
#include <features/PeakFeatures.h>
#include <io/SpectrumCache.h>
#include <extraction/Traces.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/SysInfo.h>

#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cmath>
#include <fstream>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <thread>
#include <unordered_set>

namespace OpenMS::OpenDIA
{
  namespace
  {
    constexpr std::size_t K12 = FragmentGrid::FRAGMENTS;

    [[noreturn]] void invalid(const std::string& message, const std::string& value)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, value);
    }

    /// A PQP transition ID as an integer; only canonical decimals (to_string gives the ID back) are accepted, so
    /// comparing the formatted integers compares the IDs.
    std::int64_t transitionNumber(std::string_view id)
    {
      std::int64_t v = 0;
      const auto [end, error] = std::from_chars(id.data(), id.data() + id.size(), v);
      if (error != std::errc() || end != id.data() + id.size() || std::to_string(v) != id)
        invalid("Transition IDs must be decimal integers", std::string(id));
      return v;
    }

    /// The transition IDs @p a < @p b as strings.
    bool transitionIdLess(std::int64_t a, std::int64_t b)
    {
      char x[24], y[24];
      const auto ex = std::to_chars(x, x + sizeof(x), a).ptr;
      const auto ey = std::to_chars(y, y + sizeof(y), b).ptr;
      return std::string_view(x, static_cast<std::size_t>(ex - x)) < std::string_view(y, static_cast<std::size_t>(ey - y));
    }

    bool safeText(std::string_view text) { return !text.empty() && text.find_first_of("\t\r\n") == std::string_view::npos; }

    bool sameValue(double a, double b) { return a == b || (std::isnan(a) && std::isnan(b)); }

    /// " rss_mb=<current> peak_mb=<peak>" for stage log lines (diagnostic only).
    std::string memoryNote()
    {
      std::size_t now = 0, peak = 0;
      SysInfo::getProcessMemoryConsumption(now);
      SysInfo::getProcessPeakMemoryConsumption(peak);
      return " rss_mb=" + std::to_string(now / 1024) + " peak_mb=" + std::to_string(peak / 1024);
    }

    std::string unmodified(const std::string& sequence)
    {
      try
      {
        return AASequence::fromString(sequence).toUnmodifiedString();
      }
      catch (const Exception::BaseException&)
      {
        std::string out;
        int depth = 0;
        for (char c : sequence)
        {
          if (c == '(' || c == '[') ++depth;
          else if (c == ')' || c == ']') --depth;
          else if (depth == 0 && c >= 'A' && c <= 'Z') out += c;
        }
        return out;
      }
    }

    double seconds(std::chrono::steady_clock::time_point since)
    {
      return std::chrono::duration<double>(std::chrono::steady_clock::now() - since).count();
    }

    /// Query list for StreamSearch from the library calibrated by @p cal (RTs shifted by @p rt_correction when given), in
    /// the order of the former calibrated-library sort (precursor m/z, product m/z, RT, transition ID, precursor ID).
    std::vector<TransitionQuery> makeQueries(const AnalysisLibrary& lib, const CalibrationRT& cal, const std::vector<double>& rt_correction)
    {
      const std::size_t np = lib.precursors.size();
      struct Coordinates
      {
        double seconds = 0.0, im = 0.0;
        bool has_rt = false, has_im = false;
      };
      std::vector<Coordinates> coordinates(np);
      for (std::size_t pi = 0; pi < np; ++pi)
      {
        const auto& p = lib.precursors[pi];
        double rt = p.library_rt, im = p.library_im;
        if (std::isfinite(rt)) rt = cal.map(rt);
        if (!rt_correction.empty() && std::isfinite(rt)) rt += rt_correction[pi];
        if (cal.has_im && std::isfinite(im) && im > 0.0) im = cal.im_intercept + cal.im_slope * im;
        coordinates[pi] = {rt, im, std::isfinite(rt), std::isfinite(im) && im > 0.0};
      }
      struct Key
      {
        double precursor_mz, product_mz, seconds;
        std::uint32_t precursor;
        bool has_rt;
        std::uint64_t fragment;
      };
      const double factor = 1.0 + cal.ms2_offset_ppm * 1e-6;
      std::vector<Key> keys;
      keys.reserve(lib.fragments.size());
      for (std::uint32_t pi = 0; pi < np; ++pi)
      {
        const auto& c = coordinates[pi];
        for (std::uint64_t f = lib.fragment_offset[pi]; f < lib.fragment_offset[pi + 1]; ++f)
        {
          double product = lib.fragments[f].product_mz;
          product *= factor;
          keys.push_back({lib.precursors[pi].mz, product, c.has_rt ? c.seconds : 0.0, pi, c.has_rt, f});
        }
      }
      std::stable_sort(keys.begin(), keys.end(), [&lib](const Key& l, const Key& r)
      {
        if (l.precursor_mz != r.precursor_mz) return l.precursor_mz < r.precursor_mz;
        if (l.product_mz != r.product_mz) return l.product_mz < r.product_mz;
        if (l.has_rt != r.has_rt) return !l.has_rt;
        if (l.has_rt && l.seconds != r.seconds) return l.seconds < r.seconds;
        const std::int64_t a = lib.fragments[l.fragment].transition, b = lib.fragments[r.fragment].transition;
        if (a != b) return transitionIdLess(a, b);
        return lib.precursors[l.precursor].id < lib.precursors[r.precursor].id;
      });
      std::vector<TransitionQuery> queries;
      queries.reserve(keys.size());
      for (const Key& k : keys)
      {
        const auto& c = coordinates[k.precursor];
        TransitionQuery q;
        q.precursor_id = lib.precursors[k.precursor].id;
        q.precursor_mz = k.precursor_mz;
        q.fragment_mz = k.product_mz;
        q.library_rt = c.seconds;
        q.has_library_rt = c.has_rt;
        q.precursor_im = c.im;
        q.has_precursor_im = c.has_im;
        q.is_decoy = lib.precursors[k.precursor].decoy;
        queries.push_back(q);
      }
      return queries;
    }

    /// Up to @p N fragments of precursor @p pi, most intense first (ties: lower m/z, then library order), with m/z
    /// multiplied by @p factor first (1 = library m/z).
    template<std::size_t N>
    std::size_t topFragments(const AnalysisLibrary& lib, std::uint32_t pi, double factor, std::array<double, N>& mz, std::array<double, N>& intensity)
    {
      const std::uint64_t begin = lib.fragment_offset[pi], end = lib.fragment_offset[pi + 1];
      std::vector<std::pair<double, double>> f;  // (m/z, intensity)
      f.reserve(end - begin);
      for (std::uint64_t i = begin; i < end; ++i)
      {
        double product = lib.fragments[i].product_mz;
        if (factor != 1.0) product *= factor;
        f.emplace_back(product, lib.fragments[i].library_intensity);
      }
      std::stable_sort(f.begin(), f.end(), [](const auto& a, const auto& b)
      {
        return a.second != b.second ? a.second > b.second : a.first < b.first;
      });
      mz.fill(0.0);
      intensity.fill(0.0);
      for (std::size_t k = 0; k < std::min(N, f.size()); ++k)
      {
        mz[k] = f[k].first;
        intensity[k] = f[k].second;
      }
      return f.size();
    }

    /// Forwards spectra and accumulates the time spent inside the wrapped consumer (the rest of a pass is reading).
    struct TimedConsumer final : Interfaces::IMSDataConsumer
    {
      explicit TimedConsumer(Interfaces::IMSDataConsumer& c) : inner(c) {}
      void consumeSpectrum(MSSpectrum& s) override
      {
        const auto t = std::chrono::steady_clock::now();
        inner.consumeSpectrum(s);
        consume_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
        ++spectra;
        peaks += s.size();
      }
      void consumeChromatogram(MSChromatogram& c) override { inner.consumeChromatogram(c); }
      void setExpectedSize(Size s, Size c) override { inner.setExpectedSize(s, c); }
      void setExperimentalSettings(const ExperimentalSettings& e) override { inner.setExperimentalSettings(e); }
      Interfaces::IMSDataConsumer& inner;
      double consume_s = 0.0;
      std::size_t spectra = 0, peaks = 0;
    };

    /// Forward every spectrum to two consumers in order.
    struct TeeConsumer final : Interfaces::IMSDataConsumer
    {
      TeeConsumer(Interfaces::IMSDataConsumer& a, Interfaces::IMSDataConsumer& b) : first(a), second(b) {}
      void consumeSpectrum(MSSpectrum& s) override { first.consumeSpectrum(s); second.consumeSpectrum(s); }
      void consumeChromatogram(MSChromatogram& c) override { first.consumeChromatogram(c); second.consumeChromatogram(c); }
      void setExpectedSize(Size s, Size c) override { first.setExpectedSize(s, c); second.setExpectedSize(s, c); }
      void setExperimentalSettings(const ExperimentalSettings& e) override { first.setExperimentalSettings(e); second.setExperimentalSettings(e); }
      Interfaces::IMSDataConsumer& first;
      Interfaces::IMSDataConsumer& second;
    };

    /// Stream a run into @p consumer; returns " read_s=.. consume_s=.. spectra=.. peaks=.." for the stage log.
    /// With a spectrum cache, the first call reads MS1 and MS2 from the run once and writes them to the cache (spectra
    /// of the other levels are not passed on); later calls replay the cache. Both deliver the same spectra in the same
    /// order as a direct read of @p levels.
    std::string readSpectra(const std::string& path, Interfaces::IMSDataConsumer& consumer, std::vector<Int> levels,
                            SpectrumCache* cache = nullptr)
    {
      const auto t = std::chrono::steady_clock::now();
      TimedConsumer timed(consumer);
      std::string source = "";
      if (cache != nullptr && cache->complete())
      {
        cache->replay(timed, std::vector<int>(levels.begin(), levels.end()));
        source = " source=cache";
      }
      else
      {
        FileHandler reader;
        PeakFileOptions options;
        if (cache != nullptr && cache->usable())
        {
          // Write-through: every MS1/MS2 spectrum is cached as read; only the requested levels reach the consumer.
          struct WriteThrough final : Interfaces::IMSDataConsumer
          {
            WriteThrough(SpectrumCache& c, Interfaces::IMSDataConsumer& i, std::vector<Int> l) : cache(c), inner(i), levels(std::move(l)) {}
            void consumeSpectrum(MSSpectrum& s) override
            {
              cache.append(s);
              if (std::find(levels.begin(), levels.end(), static_cast<Int>(s.getMSLevel())) != levels.end()) inner.consumeSpectrum(s);
            }
            void consumeChromatogram(MSChromatogram& c) override { inner.consumeChromatogram(c); }
            void setExpectedSize(Size s, Size c) override { inner.setExpectedSize(s, c); }
            void setExperimentalSettings(const ExperimentalSettings& e) override { inner.setExperimentalSettings(e); }
            SpectrumCache& cache;
            Interfaces::IMSDataConsumer& inner;
            std::vector<Int> levels;
          } through(*cache, timed, levels);
          options.setMSLevels({1, 2});
          reader.setOptions(options);
          try
          {
            reader.transformExperiment(path, &through, {}, ProgressLogger::NONE);
            cache->finishWriting();
          }
          catch (...)
          {
            cache->discard();
            throw;
          }
          source = " source=raw+cache_write cache_mb=" + std::to_string(cache->bytes() >> 20);
        }
        else
        {
          options.setMSLevels(levels);
          reader.setOptions(options);
          reader.transformExperiment(path, &timed, {}, ProgressLogger::NONE);
        }
      }
      const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
      return " read_s=" + std::to_string(total - timed.consume_s) + " consume_s=" + std::to_string(timed.consume_s) +
             " spectra=" + std::to_string(timed.spectra) + " peaks=" + std::to_string(timed.peaks) + source;
    }
  }

  /// Mobility half width (1/K0) of the trace tensor; also tells the matched filter the tensor's mobility bin width.
  constexpr double TRACE_IM_HALF_WIDTH = 0.03;

  std::vector<std::string> runFeatureNames()
  {
    std::vector<std::string> n;
    for (const auto& s : searchFeatureNames()) n.push_back(s);
    for (const auto& s : boundaryFeatureNames()) n.push_back(s);
    for (const auto& s : byIonFeatureNames()) n.push_back(s);
    for (const char* s : FragmentGrid::evidenceNames()) n.emplace_back(s);
    for (const char* s : FragmentGrid::deficitNames()) n.emplace_back(s);
    for (const char* s : InterferenceResolver::featureNames()) n.emplace_back(s);
    for (const auto& s : evidenceFeatureNames()) n.push_back(s);
    for (const auto& s : matchedFilterFeatureNames()) n.push_back(s);
    return n;
  }

  // ------------------------------------------------------------------------------------------------ library

  AnalysisLibrary AnalysisLibrary::read(const std::string& search_pqp, const std::string& scoring_pqp,
                                        const std::function<void(const std::string&)>& log)
  {
    AnalysisLibrary lib;
    // The annotation vocabulary is small. Intern it while streaming, and discard it once the scoring fragments
    // have been matched. No transition or precursor identity string is repeated per fragment.
    std::unordered_map<std::string, std::uint32_t> annotation_ids;
    std::vector<std::uint32_t> search_annotations;
    struct SearchReader final : LibraryVisitor
    {
      AnalysisLibrary& lib;
      std::unordered_map<std::string, std::uint32_t>& annotations;
      std::vector<std::uint32_t>& annotation_of;
      std::unordered_map<std::string, std::uint32_t> proteins;
      std::vector<std::uint32_t> precursor_of;

      SearchReader(AnalysisLibrary& l, std::unordered_map<std::string, std::uint32_t>& a, std::vector<std::uint32_t>& f) :
        lib(l), annotations(a), annotation_of(f) {}
      void expect(std::size_t n) override
      {
        lib.fragments.reserve(n);
        precursor_of.reserve(n);
        annotation_of.reserve(n);
        lib.precursors.reserve(n / 6 + 1);
        lib.index.reserve(n / 6 + 1);
      }
      void reset() override
      {
        lib = AnalysisLibrary();
        annotations.clear();
        annotation_of.clear();
        precursor_of.clear();
        proteins.clear();
      }
      void fragment(const FragmentView& f, const PrecursorView& c) override
      {
        auto it = lib.index.find(c.id);
        if (it == lib.index.end())
        {
          if (lib.precursors.size() >= UINT32_MAX) invalid("Too many library precursors", std::string(c.id));
          it = lib.index.emplace(std::string(c.id), static_cast<std::uint32_t>(lib.precursors.size())).first;
          Precursor p;
          p.id = it->first;
          p.sequence = c.sequence;
          p.charge = c.charge;
          p.library_rt = c.rt;
          p.library_im = c.drift_time;
          p.has_library_rt = std::isfinite(c.rt);
          p.mz = f.precursor_mz;
          p.decoy = f.decoy;
          std::string joined;
          for (const auto s : *c.proteins) { if (!joined.empty()) joined += ';'; joined += s; }
          auto [protein, inserted] = proteins.try_emplace(joined, static_cast<std::uint32_t>(proteins.size()));
          if (inserted) lib.protein_sets.push_back(std::move(joined));
          p.proteins = protein->second;
          lib.precursors.push_back(std::move(p));
        }
        const auto& p = lib.precursors[it->second];
        if (!sameValue(p.mz, f.precursor_mz) || p.decoy != f.decoy)
          invalid("Invalid or inconsistent library fragment", std::string(p.id));
        lib.fragments.push_back({f.product_mz, f.library_intensity, transitionNumber(f.transition_id)});
        precursor_of.push_back(it->second);
        const std::string a = fragmentAnnotation(f.fragment_type, f.fragment_nr, f.fragment_charge);
        annotation_of.push_back(annotations.try_emplace(a, static_cast<std::uint32_t>(annotations.size())).first->second);
      }
      void finish()
      {
        lib.fragment_offset.assign(lib.precursors.size() + 1, 0);
        for (const auto pi : precursor_of) ++lib.fragment_offset[pi + 1];
        std::partial_sum(lib.fragment_offset.begin(), lib.fragment_offset.end(), lib.fragment_offset.begin());
        // Preserve the relative file order within each precursor, even for interleaved fragments.
        if (!std::is_sorted(precursor_of.begin(), precursor_of.end()))
        {
          std::vector<Fragment> sorted(lib.fragments.size());
          std::vector<std::uint32_t> sorted_annotations(annotation_of.size());
          auto next = lib.fragment_offset;
          for (std::size_t i = 0; i < lib.fragments.size(); ++i)
          {
            const auto at = next[precursor_of[i]]++;
            sorted[at] = lib.fragments[i];
            sorted_annotations[at] = annotation_of[i];
          }
          lib.fragments.swap(sorted);
          annotation_of.swap(sorted_annotations);
        }
        std::vector<std::uint32_t>().swap(precursor_of);
        proteins.clear();
      }
    };
    auto traml_ids = [&]()
    {
      SearchReader reader(lib, annotation_ids, search_annotations);
      auto ids = visitAnyLibrary(search_pqp, reader, log);
      reader.finish();
      return ids;
    }();
    {
      // Keep the original map insertion and iteration order: duplicate TRAML names are resolved just as they were
      // when readPrecursorTraMLIDs fed the workflow and CalibrationCollector.
      const auto traml = traMLIDMap(traml_ids);
      std::vector<std::pair<std::string, std::string>>().swap(traml_ids);
      std::unordered_map<std::string, std::uint32_t> by_name;
      for (const auto& [id, name] : traml) by_name.emplace(name, lib.find(id));
      std::unordered_map<std::string_view, std::uint32_t> seed_partner;
      std::unordered_map<std::string, std::uint32_t> pair_index, sibling_index;
      for (std::uint32_t pi = 0; pi < lib.precursors.size(); ++pi)
      {
        const auto t = traml.find(std::string(lib.precursors[pi].id));
        if (t == traml.end()) invalid("Precursor without TRAML_ID", std::string(lib.precursors[pi].id));
        seed_partner.emplace(t->second, pi);
      }
      for (std::uint32_t pi = 0; pi < lib.precursors.size(); ++pi)
      {
        auto& p = lib.precursors[pi];
        const auto& name = traml.at(std::string(p.id));
        std::string key = name;
        if (p.decoy)
        {
          if (name.rfind("DECOY_", 0) != 0) invalid("Decoy TRAML_ID must be DECOY_<target TRAML_ID>", name);
          key = name.substr(6);
        }
        const std::string* sequence = &p.sequence;
        if (p.decoy)
        {
          const auto target = by_name.find(key);
          if (target != by_name.end() && target->second != UINT32_MAX) sequence = &lib.precursors[target->second].sequence;
        }
        const std::string family = unmodified(*sequence);
        const auto hash = CrossFitScorer::fnv1a(family);
        p.fold = static_cast<std::uint8_t>(hash % 3);
        p.rt_fold = static_cast<std::uint8_t>(hash % 5);
        p.family_seed = CrossFitScorer::fnv1a("v" + family);
        p.pair_seed = CrossFitScorer::fnv1a("v" + key);
        p.trace_safe = safeText(p.id) && safeText(family);
        p.pair = pair_index.try_emplace(key, static_cast<std::uint32_t>(pair_index.size())).first->second;
        const auto partner = seed_partner.find("DECOY_" + name);
        if (partner != seed_partner.end()) p.decoy_partner = partner->second;
        if (name.rfind("DECOY_", 0) != 0)
        {
          const auto decoy = by_name.find("DECOY_" + name);
          if (decoy != by_name.end()) lib.calibration_pairs.emplace_back(pi, decoy->second);
        }
        p.sibling_group = sibling_index.try_emplace((p.decoy ? "D" : "T") + p.sequence,
                                                   static_cast<std::uint32_t>(sibling_index.size())).first->second;
        if (p.sibling_group == lib.sibling_group_sizes.size()) lib.sibling_group_sizes.push_back(0);
        ++lib.sibling_group_sizes[p.sibling_group];
      }
      lib.pairs = pair_index.size();
      std::vector<int> pair_fold(lib.pairs, -1);
      for (const auto& p : lib.precursors) if (!p.decoy) pair_fold[p.pair] = p.fold;
      for (auto& p : lib.precursors) if (pair_fold[p.pair] >= 0) p.fold = static_cast<std::uint8_t>(pair_fold[p.pair]);
    }
    if (scoring_pqp.empty())
    {
      for (std::uint32_t pi = 0; pi < lib.precursors.size(); ++pi)
      {
        std::vector<std::uint64_t> list;
        for (auto f = lib.fragment_offset[pi]; f < lib.fragment_offset[pi + 1]; ++f) list.push_back(f);
        std::sort(list.begin(), list.end(), [&](auto a, auto b) { return lib.fragments[a].transition < lib.fragments[b].transition; });
        auto& p = lib.precursors[pi];
        p.scoring_count = static_cast<std::uint8_t>(std::min(K12, list.size()));
        for (std::size_t s = 0; s < p.scoring_count; ++s)
        {
          p.scoring_mz[s] = lib.fragments[list[s]].product_mz;
          p.scoring_library[s] = lib.fragments[list[s]].library_intensity;
          p.scoring_in_search[s] = true;
        }
      }
    }
    else
    {
      struct ScoringReader final : LibraryVisitor
      {
        AnalysisLibrary& lib;
        const std::unordered_map<std::string, std::uint32_t>& annotation_ids;
        const std::vector<std::uint32_t>& search_annotations;
        // Only the smallest twelve transition IDs are retained per precursor, never the scoring Library.
        std::vector<std::array<std::int64_t, K12>> ids;
        ScoringReader(AnalysisLibrary& l, const std::unordered_map<std::string, std::uint32_t>& a, const std::vector<std::uint32_t>& f) :
          lib(l), annotation_ids(a), search_annotations(f), ids(l.precursors.size()) {}
        void reset() override
        {
          for (auto& p : lib.precursors)
          {
            p.scoring_count = 0;
            p.scoring_mz.fill(0.0);
            p.scoring_library.fill(0.0);
            p.scoring_in_search.fill(false);
          }
        }
        void fragment(const FragmentView& f, const PrecursorView&) override
        {
          const auto pi = lib.find(f.precursor_id);
          if (pi == UINT32_MAX) return;
          auto& p = lib.precursors[pi];
          const auto id = transitionNumber(f.transition_id);
          std::size_t at = 0;
          while (at < p.scoring_count && ids[pi][at] <= id) ++at;
          if (at == K12) return;
          const auto count = std::min(K12, static_cast<std::size_t>(p.scoring_count) + 1);
          for (std::size_t s = count - 1; s > at; --s)
          {
            ids[pi][s] = ids[pi][s - 1];
            p.scoring_mz[s] = p.scoring_mz[s - 1];
            p.scoring_library[s] = p.scoring_library[s - 1];
            p.scoring_in_search[s] = p.scoring_in_search[s - 1];
          }
          ids[pi][at] = id;
          p.scoring_count = static_cast<std::uint8_t>(count);
          p.scoring_mz[at] = f.product_mz;
          p.scoring_library[at] = f.library_intensity;
          const auto a = annotation_ids.find(fragmentAnnotation(f.fragment_type, f.fragment_nr, f.fragment_charge));
          p.scoring_in_search[at] = a != annotation_ids.end() &&
            std::find(search_annotations.begin() + lib.fragment_offset[pi], search_annotations.begin() + lib.fragment_offset[pi + 1], a->second)
             != search_annotations.begin() + lib.fragment_offset[pi + 1];
        }
      };
      ScoringReader reader(lib, annotation_ids, search_annotations);
      visitAnyLibrary(scoring_pqp, reader, log);
    }
    for (auto& p : lib.precursors)
    {
      double maximum = 0.0;
      for (std::size_t s = 0; s < p.scoring_count; ++s) maximum = std::max(maximum, p.scoring_library[s]);
      for (std::size_t s = 0; s < p.scoring_count; ++s) p.scoring_library[s] /= std::max(maximum, 1e-12);
    }
    return lib;
  }

  // ------------------------------------------------------------------------------------------------ run

  RunAnalysis::RunAnalysis(AnalysisLibrary& library, RunOptions options, Log log) :
    library_(library), options_(std::move(options)), log_(std::move(log))
  {
  }

  RunResult RunAnalysis::run(const std::string& path)
  {
    RunResult result;
    result.run = path;
    const auto say = [&](const std::string& line)
    {
      result.log.push_back(line + memoryNote());
      if (log_) log_(result.log.back());
    };
    auto& lib = library_;
    const auto t_run = std::chrono::steady_clock::now();
    // Run-local spectrum cache (RunOptions::spectrum_cache_dir): written during the first read, replayed by the later
    // passes, removed when the run ends (also on an exception).
    struct CacheScope
    {
      std::unique_ptr<SpectrumCache>& cache;
      ~CacheScope() { cache.reset(); }
    } cache_scope{cache_};
    if (!options_.spectrum_cache_dir.empty())
    {
      const std::string dir = options_.spectrum_cache_dir + "/opendia_cache_" + File::getUniqueName(false);
      cache_ = std::make_unique<SpectrumCache>(dir);
      say("Spectrum cache: " + dir);
    }

    // 1. Calibration (uncalibrated library).
    const auto describe = [&lib](std::uint32_t pi, CalibrationPrecursor& out)
    {
      if (pi == UINT32_MAX || pi >= lib.precursors.size()) return false;
      const auto& p = lib.precursors[pi];
      out.id = p.id;
      out.precursor_mz = p.mz;
      out.library_rt = p.library_rt;
      out.library_im = p.library_im;
      out.is_decoy = p.decoy;
      out.fragments = topFragments(lib, pi, 1.0, out.fragment_mz, out.library_intensity);
      return out.fragments != 0;
    };
    std::vector<std::pair<std::string, double>> anchors;  // used anchors: precursor id, run RT
    {
      CalibrationOptions copt = options_.calibration;
      std::optional<CalibrationResult> calibrated;
      std::optional<CalibrationWindow> window;
      std::size_t passes = 0;
      while (!calibrated)
      {
        CalibrationCollector collector(lib.calibration_pairs, describe, copt, window);
        collector.setThreads(options_.threads);
        collector.setIonMobilityArrayAsVSSC(options_.im_array_as_vssc);
        const std::string io = readSpectra(path, collector, {2}, cache_.get());
        ++passes;
        say("Calibration: pass " + std::to_string(passes) + " collected panel=" + std::to_string(collector.panelSize()) + io);
        std::optional<CalibrationWindow> retry_window;
        try
        {
          calibrated = collector.calibrate(options_.window_panel_targets ? &retry_window : nullptr);
        }
        catch (const std::runtime_error& error)
        {
          if (!window && options_.window_panel_targets > copt.panel_targets && retry_window)
          {
            say(std::string("Calibration: whole-run pass failed (") + error.what() + "); windowed pass with up to " +
                std::to_string(options_.window_panel_targets) + " targets.");
            window = retry_window;
            copt.panel_targets = options_.window_panel_targets;
            continue;
          }
          throw;
        }
      }
      std::ostringstream json;
      writeCalibrationJSON(*calibrated, copt, json);
      result.calibration_json = json.str();
      std::istringstream in(result.calibration_json);
      result.calibration = readCalibrationRT(in, path);
      for (const auto& a : calibrated->anchors)
        if (a.used) anchors.emplace_back(a.precursor_id, a.run_rt);
      say("Calibration: passes=" + std::to_string(passes) + " used anchors=" + std::to_string(anchors.size()));
    }
    if (options_.im_tolerance_scale != 1.0 && result.calibration.has_im)
      result.calibration.im_tolerance *= options_.im_tolerance_scale;
    const CalibrationRT& cal = result.calibration;
    if (!cal.has_ms2) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Calibration has no fragment mass model", path);
    say("Calibration: rt_tolerance_s=" + std::to_string(cal.tolerance_s) + " ms2_offset_ppm=" + std::to_string(cal.ms2_offset_ppm) +
        " ms2_tolerance_ppm=" + std::to_string(cal.ms2_tolerance_ppm) + " boundary_half_widths=" + std::to_string(cal.boundary_half_widths) +
        " seconds=" + std::to_string(seconds(t_run)));

    // Steps 2-6.
    searchAndScore_(path, result, anchors, {}, cal.tolerance_s, say, t_run);
    if (!options_.rt_refine) return result;

    // Per-run ridge RT correction from the first pass's own q <= 0.01 target IDs, then a
    // second search, extraction and scoring pass with corrected RTs and a narrower RT window.
    std::vector<std::uint32_t> best(lib.precursors.size(), UINT32_MAX);
    for (std::uint32_t i = 0; i < result.windows.size(); ++i)
    {
      const auto& w = result.windows[i];
      auto& b = best[w.precursor];
      if (b == UINT32_MAX || w.score > result.windows[b].score) b = i;
    }
    std::vector<std::size_t> train;
    std::vector<double> y;
    std::vector<std::uint32_t> folds;
    for (std::uint32_t pi = 0; pi < lib.precursors.size(); ++pi)
    {
      const auto& p = lib.precursors[pi];
      if (p.decoy || !p.has_library_rt || best[pi] == UINT32_MAX || !(result.windows[best[pi]].q_value <= 0.01)) continue;
      train.push_back(pi);
      y.push_back(result.windows[best[pi]].peak.rt_apex - cal.map(p.library_rt));
      folds.push_back(p.rt_fold);
    }
    const auto refined = refineRT_(train, y, folds, cal.tolerance_s, "RT refine", say);
    if (refined.first.empty()) return result;
    // Release the first pass before the second.
    std::vector<RunWindow>().swap(result.windows);
    std::vector<float>().swap(result.features);
    result.feature_names.clear();
    searchAndScore_(path, result, anchors, refined.first, refined.second, say, t_run);
    return result;
  }

  std::pair<std::vector<double>, double> RunAnalysis::refineRT_(const std::vector<std::size_t>& train, const std::vector<double>& y,
                                                                const std::vector<std::uint32_t>& folds, double tolerance_s,
                                                                const std::string& label, const Log& say)
  {
    auto& lib = library_;
    if (train.size() < options_.rt_refine_min_ids)
    {
      say(label + ": skipped (" + std::to_string(train.size()) + " training IDs)");
      return {};
    }
    const RTRefineModel model(lib.precursors.size(), [&lib](std::size_t i) -> const std::string& { return lib.precursors[i].sequence; });
    const std::size_t dim = model.dimension();
    std::vector<double> x;
    x.reserve(train.size() * dim);
    std::vector<std::size_t> dense_train;
    dense_train.reserve(train.size());
    // Keep training order unchanged so every floating-point accumulation is identical.
    for (const auto i : train)
    {
      const auto& p = lib.precursors[i];
      const auto f = model.features(p.sequence, p.charge, p.mz);
      x.insert(x.end(), f.begin(), f.end());
      dense_train.push_back(dense_train.size());
    }
    auto refine = fitRTRefine(x, dim, dense_train, y, folds);
    std::vector<double>().swap(x);
    // Predict the full library with only one feature row resident at a time.
    refine.correction.resize(lib.precursors.size());
    for (std::size_t i = 0; i < lib.precursors.size(); ++i)
    {
      const auto& p = lib.precursors[i];
      const auto f = model.features(p.sequence, p.charge, p.mz);
      refine.correction[i] = refine.fit.predict(f.data());
    }
    const auto p95 = [](std::vector<double> v)
    {
      for (auto& e : v) e = std::abs(e);
      std::sort(v.begin(), v.end());
      // Linear interpolation between order statistics (numpy's default quantile).
      const double pos = 0.95 * static_cast<double>(v.size() - 1);
      const std::size_t lo = static_cast<std::size_t>(pos);
      const std::size_t hi = std::min(lo + 1, v.size() - 1);
      return v[lo] + (pos - static_cast<double>(lo)) * (v[hi] - v[lo]);
    };
    std::vector<double> residual(y.size());
    for (std::size_t r = 0; r < y.size(); ++r) residual[r] = y[r] - refine.out_of_fold[r];
    const double base = p95(y), refined = p95(residual);
    const double ratio = base > 0.0 ? std::clamp(refined / base, 0.3, 1.0) : 1.0;
    const double new_tol = tolerance_s * ratio;
    say(label + ": training IDs=" + std::to_string(train.size()) + ", alpha=" + std::to_string(refine.alpha) + ", p95 |residual| base=" +
        std::to_string(base) + " s refined=" + std::to_string(refined) + " s, rt_tolerance_s " + std::to_string(tolerance_s) + " -> " +
        std::to_string(new_tol));
    return {std::move(refine.correction), new_tol};
  }

  void RunAnalysis::searchAndScore_(const std::string& path, RunResult& result,
                                    const std::vector<std::pair<std::string, double>>& anchors, const std::vector<double>& rt_correction,
                                    double rt_tolerance_s, const Log& say, std::chrono::steady_clock::time_point t_run)
  {
    auto& lib = library_;
    const CalibrationRT& cal = result.calibration;

    // 2. Search with the calibrated library (RTs shifted by the per-precursor correction when given).
    const auto t_search = std::chrono::steady_clock::now();
    std::vector<std::uint32_t> search_to_library;
    std::vector<SearchCandidate> candidates;
    std::string search_io;
    {
      SearchOptions so;
      so.rt_window_seconds = rt_tolerance_s;
      so.fragment_tolerance_ppm = cal.ms2_tolerance_ppm;
      if (cal.has_im) so.ion_mobility_tolerance_vssc = cal.im_tolerance;
      so.ion_mobility_array_as_vssc = options_.im_array_as_vssc;
      so.picker.max_peaks = options_.max_peaks;
      so.picker.valley_fraction = options_.search_valley_fraction;
      so.picker.boundary_half_widths = cal.boundary_half_widths;
      so.min_transitions = options_.min_transitions;
      // Candidates are kept in memory.
      StreamSearch search(makeQueries(lib, cal, rt_correction), so);
      search_io = readSpectra(path, search, {2}, cache_.get());
      const auto t_finish = std::chrono::steady_clock::now();
      search.finish();
      search_io += " finish_s=" + std::to_string(seconds(t_finish));
      candidates = search.candidates();
      search_to_library.resize(search.precursorCount());
      for (std::size_t p = 0; p < search.precursorCount(); ++p) search_to_library[p] = lib.find(search.precursorId(p));
    }
    say("Search: candidates=" + std::to_string(candidates.size()) + " seconds=" + std::to_string(seconds(t_search)) + search_io);

    // 3. Population: every candidate of pairs whose best candidate has enough supporting fragments.
    std::vector<std::uint32_t> pair_support(lib.pairs, 0);
    for (const auto& c : candidates)
    {
      auto& s = pair_support[lib.precursors[search_to_library[c.precursor]].pair];
      s = std::max(s, c.transition_support);
    }
    for (std::uint32_t row = 0; row < candidates.size(); ++row)
    {
      const auto& c = candidates[row];
      const std::uint32_t pi = search_to_library[c.precursor];
      if (pair_support[lib.precursors[pi].pair] < options_.population_min_support) continue;
      RunWindow w;
      w.precursor = pi;
      w.candidate = row;
      w.peak = c;
      result.windows.push_back(w);
    }
    std::vector<SearchCandidate>().swap(candidates);
    auto& windows = result.windows;
    const std::size_t n = windows.size();

    // 4. One pass: 6-fragment traces + b/y ions (TraceCollector) and the 12-fragment grid (FragmentGridCollector).
    const auto t_extract = std::chrono::steady_clock::now();
    std::vector<std::uint32_t> local_of(lib.precursors.size(), UINT32_MAX);
    std::vector<std::uint32_t> local_precursors;
    for (const auto& w : windows)
    {
      if (local_of[w.precursor] == UINT32_MAX)
      {
        local_of[w.precursor] = static_cast<std::uint32_t>(local_precursors.size());
        local_precursors.push_back(w.precursor);
      }
    }
    std::vector<FragmentGridCollector::Precursor> grid_precursors(local_precursors.size());
    std::vector<InterferenceResolver::Precursor> resolver_precursors(local_precursors.size());
    // Calibrated precursor 1/K0 (0 without a library value or a mobility calibration), for the mobility-resolved grid.
    std::vector<double> calibrated_im(lib.precursors.size(), 0.0);
    if (cal.has_im)
    {
      for (std::size_t pi = 0; pi < lib.precursors.size(); ++pi)
      {
        const auto& p = lib.precursors[pi];
        if (std::isfinite(p.library_im) && p.library_im > 0.0)
          calibrated_im[pi] = cal.im_intercept + cal.im_slope * p.library_im;
      }
    }
    for (std::size_t l = 0; l < local_precursors.size(); ++l)
    {
      const auto& p = lib.precursors[local_precursors[l]];
      auto& g = grid_precursors[l];
      g.mz = p.mz;
      g.im = calibrated_im[local_precursors[l]];
      resolver_precursors[l].im = g.im;
      for (std::size_t s = 0; s < K12; ++s)
      {
        g.extraction_mz[s] = p.scoring_mz[s] > 0.0 ? p.scoring_mz[s] * (1.0 + cal.ms2_offset_ppm * 1e-6) : 0.0;
        g.library[s] = p.scoring_library[s];
        g.search[s] = p.scoring_in_search[s];
        resolver_precursors[l].fragment_mz[s] = p.scoring_mz[s];
        resolver_precursors[l].library[s] = p.scoring_library[s];
      }
    }
    // Per window, only the 12-fragment central areas outlive its extraction group (fragment deficit and noise floor);
    // the grid profile is used by the interference features, computed inside the group.
    std::vector<FragmentGrid::Values> core_areas(n);
    std::vector<float> interference(n * InterferenceResolver::FEATURES, 0.0f);
    std::vector<TracePrecursor> trace_precursors(local_precursors.size());
    for (std::size_t l = 0; l < local_precursors.size(); ++l)
    {
      const auto pi = local_precursors[l];
      const auto& p = lib.precursors[pi];
      if (!p.trace_safe) invalid("Invalid candidate identity or group", std::string(p.id));
      auto& t = trace_precursors[l];
      t.precursor_mz = p.mz;
      t.precursor_im = p.library_im;
      if (cal.has_im && std::isfinite(t.precursor_im) && t.precursor_im > 0.0)
        t.precursor_im = cal.im_intercept + cal.im_slope * t.precursor_im;
      t.charge = p.charge;
      t.is_decoy = p.decoy;
      const double factor = 1.0 + cal.ms2_offset_ppm * 1e-6;
      for (auto f = lib.fragment_offset[pi]; f < lib.fragment_offset[pi + 1]; ++f)
      {
        const auto& fragment = lib.fragments[f];
        const double mz = fragment.product_mz * factor;
        if (!std::isfinite(mz) || mz <= 0.0 || !std::isfinite(p.mz) || p.mz <= 0.0 ||
            !std::isfinite(fragment.library_intensity) || fragment.library_intensity < 0.0)
          invalid("Invalid or inconsistent library fragment", std::string(p.id));
      }
      // The 6-fragment tensor uses calibrated m/z and precursor 1/K0, resolved here without a Library copy.
      topFragments(lib, pi, factor, t.fragment_mz, t.library_intensity);
    }
    // Streaming tensors: a window's 6-fragment tensor is held only while its RT range is being read. The tensor-derived
    // features are computed when the window closes and the tensor is released, so tensor memory follows the windows
    // co-eluting at one time.
    PeakFeatureOptions feature_options;
    feature_options.rt_half_width = options_.grid_half_width_s;
    feature_options.core_half_width = options_.core_half_width_s;
    feature_options.im_half_width = TRACE_IM_HALF_WIDTH;
    std::vector<float> boundary(n * BOUNDARY_FEATURES, 0.0f);
    std::vector<float> evidence_block(n * EVIDENCE_FEATURES, 0.0f);
    std::vector<float> mf_block(n * MATCHED_FILTER_FEATURES, 0.0f);
    std::vector<std::array<std::uint16_t, 2>> byions(n);
    std::vector<std::array<float, FragmentGrid::EVIDENCE>> grid_evidence(n);
    const int threads = options_.threads;

    // Extraction groups: precursors contiguous in m/z with about equal window counts. Every per-window value depends only
    // on the spectra whose isolation windows contain the window's precursor m/z (traces, MS1 evidence, b/y ions and the
    // 12-fragment grid), so a group that receives all MS1 spectra and every MS2 spectrum with an isolation window
    // overlapping its m/z range computes exactly the values of a single pass, and several groups run concurrently.
    // Collector memory follows the groups in flight. Groups need the cache (one replay each); without it, one group.
    const bool partition = cache_ && cache_->complete();
    // Isolation windows of the run, as the 12-fragment grid keys them: every valid isolation window of a non-empty MS2
    // spectrum with a finite RT. With the cache they are known before extraction (from its index); otherwise the single
    // group's collector reports them.
    std::vector<std::pair<double, double>> isolation_keys;
    const auto set_isolation = [&]()
    {
      std::vector<double> mz(local_precursors.size());
      for (std::size_t l = 0; l < local_precursors.size(); ++l) mz[l] = lib.precursors[local_precursors[l]].mz;
      const auto isolation = FragmentGridCollector::isolationIndices(isolation_keys, mz);
      for (std::size_t l = 0; l < resolver_precursors.size(); ++l) resolver_precursors[l].isolation = isolation[l];
    };
    if (partition)
    {
      std::set<std::pair<double, double>> keys;
      const auto& isolation = cache_->precursors();
      for (const auto& e : cache_->index())
      {
        if (e.level != 2 || e.peaks == 0 || !std::isfinite(e.rt)) continue;
        for (std::size_t k = e.precursor_begin; k < std::size_t{e.precursor_begin} + e.precursor_count; ++k)
        {
          std::pair<double, double> key;
          if (FragmentGridCollector::isolationKey(isolation[k][0], isolation[k][1], isolation[k][2], key)) keys.insert(key);
        }
      }
      isolation_keys.assign(keys.begin(), keys.end());
      set_isolation();
    }
    std::size_t group_count = partition ? (options_.extraction_groups == 0 ? 16 : options_.extraction_groups) : 1;
    group_count = std::clamp<std::size_t>(group_count, 1, std::max<std::size_t>(1, local_precursors.size()));
    std::vector<std::vector<std::uint32_t>> group_members;      // global window indices, ascending
    std::vector<std::pair<double, double>> group_mz;            // precursor m/z range of each group
    {
      std::vector<std::uint32_t> by_mz(local_precursors.size());
      std::iota(by_mz.begin(), by_mz.end(), 0u);
      std::sort(by_mz.begin(), by_mz.end(), [&](std::uint32_t x, std::uint32_t y)
      {
        const double mx = lib.precursors[local_precursors[x]].mz, my = lib.precursors[local_precursors[y]].mz;
        return mx != my ? mx < my : x < y;
      });
      std::vector<std::uint32_t> windows_of(local_precursors.size(), 0);
      for (const auto& w : windows) ++windows_of[local_of[w.precursor]];
      std::vector<std::uint32_t> group_of(local_precursors.size(), 0);
      std::size_t acc = 0, g = 0;
      // Groups hold whole isolation windows (the precursor isolation index is nondecreasing in m/z), so the interference
      // features, which link windows of one isolation window only, are computed inside the group.
      std::int32_t last_isolation = -1;
      for (const std::uint32_t l : by_mz)
      {
        const std::int32_t iso = resolver_precursors[l].isolation;
        const bool may_cut = iso >= 0 && iso != last_isolation;
        if (may_cut && g + 1 < group_count && acc * group_count >= (g + 1) * n) ++g;
        if (iso >= 0) last_isolation = iso;
        group_of[l] = static_cast<std::uint32_t>(g);
        acc += windows_of[l];
      }
      group_count = local_precursors.empty() ? 1 : g + 1;
      group_members.assign(group_count, {});
      group_mz.assign(group_count, {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()});
      for (std::uint32_t i = 0; i < n; ++i)
      {
        const auto l = local_of[windows[i].precursor];
        group_members[group_of[l]].push_back(i);
        auto& r = group_mz[group_of[l]];
        r.first = std::min(r.first, lib.precursors[local_precursors[l]].mz);
        r.second = std::max(r.second, lib.precursors[local_precursors[l]].mz);
      }
    }
    std::size_t parallel = 1;
    if (partition)
    {
      parallel = options_.extraction_parallel == 0 ? static_cast<std::size_t>(std::max(1, threads / 4)) : options_.extraction_parallel;
      parallel = std::clamp<std::size_t>(parallel, 1, group_count);
    }
    const int group_threads = std::max(1, threads / static_cast<int>(parallel));

    struct GroupStats
    {
      std::size_t peak_open = 0, peak_held = 0, spectra = 0, peaks = 0;
      double consume_s = 0.0, interference_s = 0.0;
      std::string io;
      std::vector<std::pair<double, double>> isolation;  ///< Isolation windows seen (single group only).
    };
    std::vector<GroupStats> group_stats(group_count);
    std::mutex byion_setup;
    const auto extract_group = [&](std::size_t gi, GroupStats& st)
    {
      const auto& members = group_members[gi];
      // Group-local precursors, first seen in window order.
      std::vector<std::uint32_t> group_local(local_precursors.size(), UINT32_MAX), group_precursors;
      for (const auto i : members)
      {
        const auto l = local_of[windows[i].precursor];
        if (group_local[l] == UINT32_MAX)
        {
          group_local[l] = static_cast<std::uint32_t>(group_precursors.size());
          group_precursors.push_back(l);
        }
      }
      std::vector<FragmentGridCollector::Precursor> gp(group_precursors.size());
      std::vector<TracePrecursor> tp(group_precursors.size());
      for (std::size_t k = 0; k < group_precursors.size(); ++k)
      {
        gp[k] = grid_precursors[group_precursors[k]];
        tp[k] = trace_precursors[group_precursors[k]];
      }
      std::vector<InterferenceResolver::Window> gw(members.size());
      std::vector<double> apex(members.size());
      std::vector<std::uint32_t> precursor_of(members.size());
      for (std::size_t j = 0; j < members.size(); ++j)
      {
        const auto& src = windows[members[j]].peak;
        gw[j].precursor = group_local[local_of[windows[members[j]].precursor]];
        gw[j].rt_start = src.rt_start;
        gw[j].rt_apex = src.rt_apex;
        gw[j].rt_end = src.rt_end;
        apex[j] = src.rt_apex;
        precursor_of[j] = gw[j].precursor;
      }
      std::vector<std::uint32_t>().swap(group_local);
      std::vector<InterferenceResolver::Precursor> rp(group_precursors.size());
      for (std::size_t k = 0; k < group_precursors.size(); ++k) rp[k] = resolver_precursors[group_precursors[k]];
      FragmentGridCollector grids(gp, gw, options_.grid_half_width_s, options_.grid_extraction_ppm);
      // With a mobility calibration the 12-fragment grid (and with it the interference inputs) is extracted within the
      // calibrated 1/K0 tolerance of each precursor, as the search is; without one it pools all points.
      if (cal.has_im && cal.im_tolerance > 0.0) grids.enableIonMobility(cal.im_tolerance, options_.im_array_as_vssc);
      grids.setThreads(group_threads);
      TraceCollector::CloseCallback on_close = [&](const std::vector<std::uint32_t>& closed, const TraceCollector& tc)
      {
#pragma omp parallel num_threads(group_threads)
        {
          std::array<float, TraceCollector::SAMPLE_SIZE> view;
#pragma omp for schedule(dynamic, 256)
          for (std::ptrdiff_t c = 0; c < static_cast<std::ptrdiff_t>(closed.size()); ++c)
          {
            const std::size_t j = closed[static_cast<std::size_t>(c)];
            const std::size_t i = members[j];
            const float* tensor = tc.windowTensor(j, view);
            const auto b = boundaryFeatures(tensor, tc.libraryIntensities(j), tc.fragmentMZ(j), tc.evidence(j), feature_options);
            for (std::size_t k = 0; k < BOUNDARY_FEATURES; ++k) boundary[i * BOUNDARY_FEATURES + k] = static_cast<float>(b[k]);
            const auto x = evidenceFeatures(tensor, tc.libraryIntensities(j), tc.fragmentMZ(j), tc.evidence(j), feature_options);
            for (std::size_t k = 0; k < EVIDENCE_FEATURES; ++k) evidence_block[i * EVIDENCE_FEATURES + k] = static_cast<float>(x[k]);
            const auto m = matchedFilterFeatures(tensor, tc.libraryIntensities(j), tc.fragmentMZ(j), tc.evidence(j), feature_options);
            for (std::size_t k = 0; k < MATCHED_FILTER_FEATURES; ++k) mf_block[i * MATCHED_FILTER_FEATURES + k] = static_cast<float>(m[k]);
          }
        }
      };
      TraceCollector traces(apex, precursor_of, tp, 0.02, options_.grid_half_width_s, TRACE_IM_HALF_WIDTH, members.size(),
                            options_.trace_tensor_gb << 30, on_close);
      std::vector<double>().swap(apex);
      std::vector<std::uint32_t>().swap(precursor_of);
      std::vector<TracePrecursor>().swap(tp);
      traces.setIonMobilityArrayAsVSSC(options_.im_array_as_vssc);
      traces.setFragmentTolerancePpm(options_.trace_tolerance_ppm);
      traces.setThreads(group_threads);
      traces.enableEvidence(20.0, options_.core_half_width_s);
      {
        // Sequence parsing reaches the shared residue and modification databases; set up one group at a time.
        std::lock_guard<std::mutex> lock(byion_setup);
        traces.enableByIonCounts([&](std::size_t j) -> const std::string& { return lib.precursors[windows[members[j]].precursor].sequence; },
                                 [&](std::size_t j) -> std::uint64_t { return windows[members[j]].precursor; },
                                 options_.byion_ppm, options_.byion_intensity_min);
      }
      {
        TeeConsumer both(grids, traces);
        if (partition)
        {
          // MS1 always; MS2 when one of its isolation windows (as the trace and grid collectors test them) overlaps
          // the group's precursor m/z range.
          const double lo = group_mz[gi].first, hi = group_mz[gi].second;
          const auto& isolation = cache_->precursors();
          TimedConsumer timed(both);
          cache_->replay(timed, {1, 2}, 64, [&](const SpectrumCache::Entry& e)
          {
            if (e.level != 2) return true;
            for (std::size_t k = e.precursor_begin; k < std::size_t{e.precursor_begin} + e.precursor_count; ++k)
            {
              const auto& v = isolation[k];
              if (!(v[1] >= 0.0 && v[2] >= 0.0 && v[1] + v[2] > 0.0)) continue;
              if (v[0] - v[1] <= hi && v[0] + v[2] >= lo) return true;
            }
            return false;
          });
          st.consume_s = timed.consume_s;
          st.spectra = timed.spectra;
          st.peaks = timed.peaks;
        }
        else
        {
          st.io = readSpectra(path, both, {1, 2}, cache_.get());
        }
      }
      traces.finish();
      grids.finish();
      for (std::size_t j = 0; j < members.size(); ++j)
      {
        const std::size_t i = members[j];
        core_areas[i] = gw[j].core_area;
        grid_evidence[i] = grids.evidence()[j];
        byions[i] = traces.byIonCounts(j);
      }
      st.peak_open = grids.peakOpenGrids();
      st.peak_held = traces.peakHeldTensors();
      if (!partition)
      {
        // Single group: it saw every spectrum, so its collector knows the run's isolation windows.
        st.isolation = grids.isolationKeys();
        const auto isolation = grids.isolationIndices();
        for (std::size_t k = 0; k < rp.size(); ++k) rp[k].isolation = isolation[k];
      }
      // Interference within the group: every window's partners share its isolation window, which lies in this group.
      // Window and precursor order are those of the run restricted to the group, so the features are unchanged.
      const auto t_i = std::chrono::steady_clock::now();
      const auto group_interference = InterferenceResolver::compute(rp, gw, cal.ms2_tolerance_ppm,
                                                                    2.0 * options_.grid_half_width_s / InterferenceResolver::BINS, group_threads,
                                                                    nullptr, cal.has_im ? 2.0 * cal.im_tolerance : 0.0);
      st.interference_s = seconds(t_i);
      for (std::size_t j = 0; j < members.size(); ++j)
        for (std::size_t k = 0; k < InterferenceResolver::FEATURES; ++k)
          interference[members[j] * InterferenceResolver::FEATURES + k] = group_interference[j * InterferenceResolver::FEATURES + k];
    };
    {
      std::vector<std::exception_ptr> errors(group_count);
      std::atomic<std::size_t> next{0};
      std::atomic<bool> failed{false};
      const auto worker = [&]()
      {
        for (std::size_t g; !failed && (g = next++) < group_count;)
        {
          try
          {
            extract_group(g, group_stats[g]);
          }
          catch (...)
          {
            errors[g] = std::current_exception();
            failed = true;
          }
        }
      };
      if (parallel == 1)
      {
        worker();
      }
      else
      {
        std::vector<std::thread> pool;
        for (std::size_t t = 0; t < parallel; ++t) pool.emplace_back(worker);
        for (auto& t : pool) t.join();
      }
      for (const auto& e : errors)
        if (e) std::rethrow_exception(e);
    }
    std::vector<TracePrecursor>().swap(trace_precursors);
    if (!partition)
    {
      isolation_keys = std::move(group_stats[0].isolation);
      set_isolation();
    }
    std::string extract_io;
    std::size_t peak_open = 0, peak_held = 0;
    for (const auto& st : group_stats)
    {
      peak_open = std::max(peak_open, st.peak_open);
      peak_held = std::max(peak_held, st.peak_held);
    }
    if (partition)
    {
      double consume_s = 0.0;
      std::size_t spectra = 0, peaks = 0;
      for (const auto& st : group_stats)
      {
        consume_s += st.consume_s;
        spectra += st.spectra;
        peaks += st.peaks;
      }
      extract_io = " groups=" + std::to_string(group_count) + " parallel=" + std::to_string(parallel) + " group_threads=" +
                   std::to_string(group_threads) + " consume_s_sum=" + std::to_string(consume_s) + " spectra_sum=" + std::to_string(spectra) +
                   " peaks_sum=" + std::to_string(peaks) + " source=cache";
    }
    else
    {
      extract_io = group_stats[0].io;
    }
    say("Extraction: windows=" + std::to_string(n) + " precursors=" + std::to_string(local_precursors.size()) +
        " isolation_windows=" + std::to_string(isolation_keys.size()) + " peak_open_grids=" + std::to_string(peak_open) +
        " peak_held_tensors=" + std::to_string(peak_held) + " (max per group)" +
        " seconds=" + std::to_string(seconds(t_extract)) + extract_io);

    // 5. Features.
    const auto t_features = std::chrono::steady_clock::now();
    const auto names = runFeatureNames();
    const std::size_t d = names.size();
    result.feature_names = names;
    // Capacity for the sibling columns added after the first scorer pass, so widening never reallocates (untouched
    // capacity is not resident).
    result.features.reserve(n * (d + SIBLING_FEATURES));
    result.features.assign(n * d, 0.0f);
    double noise = 1.0;
    {
      // Median of the nonzero central areas of library fragments over all windows (FragmentGrid::noiseLevel()).
      std::vector<double> v;
      for (std::size_t i = 0; i < n; ++i)
      {
        const auto& lib12 = resolver_precursors[local_of[windows[i].precursor]].library;
        for (std::size_t f = 0; f < K12; ++f)
          if (core_areas[i][f] > 0.0 && lib12[f] > 0.0) v.push_back(core_areas[i][f]);
      }
      noise = FragmentGrid::noiseLevel(v);
    }
    double interference_s = 0.0;
    for (const auto& st : group_stats) interference_s += st.interference_s;
#pragma omp parallel num_threads(options_.threads)
    {
#pragma omp for schedule(dynamic, 1024)
      for (std::ptrdiff_t ii = 0; ii < static_cast<std::ptrdiff_t>(n); ++ii)
      {
        const std::size_t i = static_cast<std::size_t>(ii);
        const auto& w = windows[i];
        float* f = result.features.data() + i * d;
        std::size_t o = 0;
        SearchPeakFeatures peak;
        peak.transition_support = w.peak.transition_support;
        peak.matched_scans = w.peak.matched_scans;
        peak.peak_width_seconds = w.peak.rt_end - w.peak.rt_start;
        peak.abs_rt_error_seconds = lib.precursors[w.precursor].has_library_rt ? std::abs(w.peak.apex_rt_error_seconds) : 0.0;
        for (double v : searchFeatures(peak)) f[o++] = static_cast<float>(v);
        for (std::size_t k = 0; k < BOUNDARY_FEATURES; ++k) f[o++] = boundary[i * BOUNDARY_FEATURES + k];
        for (const auto v : byions[i]) f[o++] = static_cast<float>(static_cast<double>(v));
        for (float v : grid_evidence[i]) f[o++] = v;
        const auto& gp = grid_precursors[local_of[w.precursor]];
        for (float v : FragmentGrid::deficit(core_areas[i], gp.library, gp.search, noise)) f[o++] = v;
        for (std::size_t k = 0; k < InterferenceResolver::FEATURES; ++k) f[o++] = interference[i * InterferenceResolver::FEATURES + k];
        for (std::size_t k = 0; k < EVIDENCE_FEATURES; ++k) f[o++] = evidence_block[i * EVIDENCE_FEATURES + k];
        for (std::size_t k = 0; k < MATCHED_FILTER_FEATURES; ++k) f[o++] = mf_block[i * MATCHED_FILTER_FEATURES + k];
      }
    }
    say("Features: " + std::to_string(d) + " per window, noise=" + std::to_string(noise) + " seconds=" + std::to_string(seconds(t_features)) +
        " interference_seconds_sum=" + std::to_string(interference_s));
    // Everything below reads only result.features, windows and the library. Release the extraction state now (about
    // about 180 bytes per window) instead of holding it through scoring.
    std::vector<std::array<std::uint16_t, 2>>().swap(byions);
    std::vector<std::array<float, FragmentGrid::EVIDENCE>>().swap(grid_evidence);
    std::vector<FragmentGrid::Values>().swap(core_areas);
    std::vector<float>().swap(interference);
    std::vector<FragmentGridCollector::Precursor>().swap(grid_precursors);
    std::vector<InterferenceResolver::Precursor>().swap(resolver_precursors);
    std::vector<float>().swap(boundary);
    std::vector<float>().swap(evidence_block);
    std::vector<float>().swap(mf_block);

    // 6. Scorer.
    const auto t_score = std::chrono::steady_clock::now();
    CrossFitScorer::Input input;
    input.dimension = d;
    // The scorer input owns the feature table during scoring; result.features is restored when scoring ends.
    input.features = std::move(result.features);
    std::vector<std::uint32_t> dense_precursor(lib.precursors.size(), UINT32_MAX), dense_pair(lib.pairs, UINT32_MAX);
    std::uint32_t np = 0, npair = 0;
    for (const auto& w : windows)
    {
      const auto& p = lib.precursors[w.precursor];
      if (dense_precursor[w.precursor] == UINT32_MAX) dense_precursor[w.precursor] = np++;
      if (dense_pair[p.pair] == UINT32_MAX) dense_pair[p.pair] = npair++;
      input.precursor.push_back(dense_precursor[w.precursor]);
      input.pair.push_back(dense_pair[p.pair]);
      input.decoy.push_back(p.decoy);
      input.fold.push_back(p.fold);
      input.pair_hash.push_back(p.pair_seed);
    }
    // Seeds: the window nearest each used anchor (|apex - anchor run RT| <= anchor_apex_s) is positive; the decoy partner's
    // best search window (most supporting fragments, then most scans, then lowest row) and a sample of other decoy windows
    // are negative.
    std::unordered_map<std::uint32_t, double> anchor_rt;
    for (const auto& [id, rt] : anchors)
    {
      const auto it = lib.index.find(id);
      if (it != lib.index.end()) anchor_rt.emplace(it->second, rt);
    }
    std::vector<std::uint32_t> best_of(lib.precursors.size(), UINT32_MAX), anchor_window(lib.precursors.size(), UINT32_MAX);
    for (std::uint32_t i = 0; i < n; ++i)
    {
      const auto& w = windows[i];
      auto& b = best_of[w.precursor];
      if (b == UINT32_MAX || w.peak.transition_support > windows[b].peak.transition_support ||
          (w.peak.transition_support == windows[b].peak.transition_support && w.peak.matched_scans > windows[b].peak.matched_scans))
      {
        b = i;
      }
      const auto a = anchor_rt.find(w.precursor);
      if (a == anchor_rt.end()) continue;
      const double dt = std::abs(w.peak.rt_apex - a->second);
      auto& aw = anchor_window[w.precursor];
      if (dt <= options_.anchor_apex_s && (aw == UINT32_MAX || dt < std::abs(windows[aw].peak.rt_apex - a->second))) aw = i;
    }
    std::unordered_set<std::uint32_t> negative_rows;
    std::size_t positives = 0;
    for (std::uint32_t pi = 0; pi < lib.precursors.size(); ++pi)
    {
      if (anchor_window[pi] == UINT32_MAX) continue;
      const auto& p = lib.precursors[pi];
      input.seeds.push_back({anchor_window[pi], true, p.family_seed});
      ++positives;
      const auto partner = p.decoy_partner;
      if (partner != UINT32_MAX && best_of[partner] != UINT32_MAX && negative_rows.insert(best_of[partner]).second)
      {
        input.seeds.push_back({best_of[partner], false, lib.precursors[partner].family_seed});
      }
    }
    {
      std::vector<std::uint32_t> decoys;
      for (std::uint32_t i = 0; i < n; ++i)
        if (lib.precursors[windows[i].precursor].decoy && !negative_rows.contains(i)) decoys.push_back(i);
      // Deterministic sample (splitmix64 partial shuffle).
      std::uint64_t state = options_.seed + 13;
      const auto next = [&state]()
      {
        std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
      };
      const std::size_t take = std::min(options_.decoy_seeds, decoys.size());
      for (std::size_t k = 0; k < take; ++k)
      {
        std::swap(decoys[k], decoys[k + next() % (decoys.size() - k)]);
        const auto& p = lib.precursors[windows[decoys[k]].precursor];
        input.seeds.push_back({decoys[k], false, p.family_seed});
      }
    }
    say("Seeds: positives=" + std::to_string(positives) + " negatives=" + std::to_string(input.seeds.size() - positives));
    CrossFitScorer::Options so;
    so.seed = options_.seed;
    so.threads = options_.threads;
    auto scored = CrossFitScorer::run(input, so, [&](const std::string& l) { say("Scorer " + l); });
    {
      // Charge-state sibling rescoring (PLAN_SIB): siblings share the family fold, so first-pass scores are out-of-fold for each other.
      std::vector<std::uint32_t> group(n), prec(n);
      std::vector<double> rt(n);
      for (std::size_t i = 0; i < n; ++i)
      {
        prec[i] = windows[i].precursor;
        group[i] = lib.precursors[prec[i]].sibling_group;
        rt[i] = windows[i].peak.rt_apex;
      }
      const auto sib = siblingFeatures(group, prec, rt, scored.score, lib.sibling_group_sizes, CrossFitScorer::NOT_KEPT);
      std::vector<std::uint32_t>().swap(group);
      std::vector<double>().swap(rt);
      constexpr std::size_t S = SIBLING_FEATURES;
      const std::size_t d2 = d + S;
      // Widen in place, last row first: row i moves from i * d to i * d2 >= i * d, so no unread value is overwritten.
      auto& f = input.features;
      f.resize(n * d2);
      std::size_t co_t = 0, co_d = 0;
      for (std::size_t i = n; i-- > 0;)
      {
        std::copy_backward(f.data() + i * d, f.data() + i * d + d, f.data() + i * d2 + d);
        std::copy_n(sib.data() + i * S, S, f.data() + i * d2 + d);
        if (sib[i * S + 1] > 0.0f) (lib.precursors[prec[i]].decoy ? co_d : co_t) += 1;
      }
      for (const char* nm : {"s_n", "s_co", "s_best", "s_best_d"}) result.feature_names.emplace_back(nm);
      input.dimension = d2;
      scored = CrossFitScorer::run(input, so, [&](const std::string& l) { say("Scorer " + l); });
      say("Sibling features: windows with a co-eluting sibling target=" + std::to_string(co_t) + " decoy=" + std::to_string(co_d) + "; rescored");
    }
    {
      // Hard-negative GBM rescoring (PLAN_GBT) on the hardest 2 % of decoy precursors.
      for (std::size_t i = 0; i < n; ++i)
      {
        windows[i].network_score = scored.score[i];
        windows[i].network_q_value = scored.q_value[i];
      }
      auto gbm = gbmRescore(input, scored.score, 0.02, options_.seed, options_.threads, say, options_.gbm_rounds, options_.gbm_positive_q);
      scored.score = std::move(gbm.score);
      scored.q_value = std::move(gbm.q_value);
    }
    result.features = std::move(input.features);
    for (std::size_t i = 0; i < n; ++i)
    {
      windows[i].score = scored.score[i];
      windows[i].q_value = scored.q_value[i];
    }
    std::size_t at1 = 0;
    for (std::size_t i = 0; i < n; ++i) at1 += !lib.precursors[windows[i].precursor].decoy && windows[i].q_value <= 0.01;
    say("Scoring: target precursors at q <= 0.01: " + std::to_string(at1) + " seconds=" + std::to_string(seconds(t_score)) +
        "; run total seconds=" + std::to_string(seconds(t_run)));
  }
}
