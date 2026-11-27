// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#include <io/PeakMapExport.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/QPXIdentity.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace OpenMS::OpenDIA
{
  namespace
  {
    void checked(const arrow::Status& status)
    {
      if (!status.ok()) throw std::runtime_error("XIPM export: " + status.ToString());
    }
    template <typename T> T take(arrow::Result<T> result)
    {
      checked(result.status());
      return std::move(result).ValueOrDie();
    }
    struct Map
    {
      const RunWindow* window{};
      std::int64_t transition{-1};
      double target{};
      std::vector<double> mz, rt, im, intensity;
    };
    class Collector final : public Interfaces::IMSDataConsumer
    {
    public:
      Collector(const AnalysisLibrary& library, const RunResult& run, const PeakMapExportOptions& options,
                std::vector<const RunWindow*> selected, std::function<void(const Map&)> emit) :
        library_(library), run_(run), options_(options), pending_(std::move(selected)), emit_(std::move(emit))
      {
        std::stable_sort(pending_.begin(), pending_.end(), [](const auto* a, const auto* b)
        {
          return a->peak.rt_start < b->peak.rt_start;
        });
      }
      void finish()
      {
        for (const auto& map : maps_) emit_(map);
        maps_.clear();
        while (next_ < pending_.size())
        {
          activate_(pending_[next_++]);
          for (const auto& map : maps_) emit_(map);
          maps_.clear();
        }
      }
      void setExpectedSize(std::size_t, std::size_t) override {}
      void setExperimentalSettings(const ExperimentalSettings&) override {}
      void consumeChromatogram(MSChromatogram&) override {}
      void consumeSpectrum(MSSpectrum& spectrum) override
      {
        const auto level = spectrum.getMSLevel();
        if (level != 1 && level != 2) return;
        const double rt = spectrum.getRT();
        if (!std::isfinite(rt) || rt < last_rt_)
          throw std::runtime_error("XIPM export requires finite, nondecreasing spectrum RT in replay");
        last_rt_ = rt;
        // Retire complete intervals before activating new ones. Only overlapping windows retain points.
        const auto remaining = std::remove_if(maps_.begin(), maps_.end(), [&](const Map& map)
        {
          if (map.window->peak.rt_end + options_.rt_padding_s >= rt) return false;
          emit_(map);
          points_ -= map.mz.size();
          return true;
        });
        maps_.erase(remaining, maps_.end());
        while (next_ < pending_.size() && pending_[next_]->peak.rt_start - options_.rt_padding_s <= rt)
        {
          const auto* window = pending_[next_++];
          const auto begin = maps_.size();
          activate_(window);
          if (window->peak.rt_end + options_.rt_padding_s < rt)
          {
            for (std::size_t i = begin; i < maps_.size(); ++i) emit_(maps_[i]);
            maps_.resize(begin);
          }
        }
        if (!spectrum.isSorted()) spectrum.sortByPosition();
        const DataArrays::FloatDataArray* mobility = nullptr;
        if (spectrum.containsIMData())
        {
          const auto [index, unit] = spectrum.getIMData();
          if (spectrum.getFloatDataArrays()[index].size() == spectrum.size() &&
              (unit == DriftTimeUnit::VSSC || (options_.im_array_as_vssc && unit != DriftTimeUnit::CCS)))
          {
            mobility = &spectrum.getFloatDataArrays()[index];
          }
        }
        const double scalar = spectrum.getDriftTime();
        const bool scalar_valid = spectrum.getDriftTimeUnit() == DriftTimeUnit::VSSC && std::isfinite(scalar) && scalar > 0;
        for (auto& map : maps_)
        {
          if ((map.transition < 0 ? 1u : 2u) != level) continue;
          const auto& window = *map.window;
          if (spectrum.getRT() < window.peak.rt_start - options_.rt_padding_s ||
              spectrum.getRT() > window.peak.rt_end + options_.rt_padding_s) continue;
          const auto& precursor = library_.precursors.at(window.precursor);
          const double target_im = run_.calibration.has_im && precursor.library_im > 0 ?
            run_.calibration.im_intercept + run_.calibration.im_slope * precursor.library_im : -1;
          // Each precursor's isolation rectangle is tested separately, never their bounding box.
          const auto isolates = [&](double im)
          {
            for (const auto& p : spectrum.getPrecursors())
            {
              const double lo = p.getIsolationWindowLowerOffset(), hi = p.getIsolationWindowUpperOffset();
              if (!(lo >= 0 && hi >= 0 && lo + hi > 0) ||
                  !(precursor.mz >= p.getMZ() - lo && precursor.mz <= p.getMZ() + hi)) continue;
              if (p.getDriftTimeUnit() == DriftTimeUnit::VSSC && p.getDriftTime() > 0 &&
                  p.getDriftTimeWindowLowerOffset() + p.getDriftTimeWindowUpperOffset() > 0)
              {
                const double position = std::isfinite(im) ? im : target_im;
                if (position > 0 && !(position >= p.getDriftTime() - p.getDriftTimeWindowLowerOffset() &&
                                      position <= p.getDriftTime() + p.getDriftTimeWindowUpperOffset())) continue;
              }
              return true;
            }
            return false;
          };
          if (level == 2 && std::none_of(spectrum.getPrecursors().begin(), spectrum.getPrecursors().end(), [&](const auto& p)
          {
            const double lo = p.getIsolationWindowLowerOffset(), hi = p.getIsolationWindowUpperOffset();
            return lo >= 0 && hi >= 0 && lo + hi > 0 &&
              precursor.mz >= p.getMZ() - lo && precursor.mz <= p.getMZ() + hi;
          })) continue;
          const double ppm = level == 2 && run_.calibration.has_ms2 && run_.calibration.ms2_tolerance_ppm > 0 ?
            run_.calibration.ms2_tolerance_ppm : 10.0;
          const double tolerance = map.target * ppm * 1e-6;
          const auto first = std::lower_bound(spectrum.begin(), spectrum.end(), map.target - tolerance,
            [](const Peak1D& peak, double mz) { return peak.getMZ() < mz; });
          for (auto it = first; it != spectrum.end() && it->getMZ() <= map.target + tolerance; ++it)
          {
            const auto& peak = *it;
            const auto i = static_cast<std::size_t>(it - spectrum.begin());
            if (!std::isfinite(peak.getMZ()) || std::abs(peak.getMZ() - map.target) > tolerance) continue;
            const double im = mobility ? (*mobility)[i] : scalar_valid ? scalar : std::numeric_limits<double>::quiet_NaN();
            if (level == 2 && !isolates(im)) continue;
            if (target_im > 0 && std::isfinite(im) && std::abs(im - target_im) > run_.calibration.im_tolerance) continue;
            if (++points_ > 4000000) throw std::runtime_error("XIPM export resident point limit (4000000) exceeded, reduce RT padding");
            map.mz.push_back(peak.getMZ());
            map.rt.push_back(spectrum.getRT());
            map.im.push_back(im);
            map.intensity.push_back(peak.getIntensity());
          }
        }
      }
    private:
      void activate_(const RunWindow* window)
      {
        const auto index = window->precursor;
        if (!std::isfinite(library_.precursors[index].mz) || library_.precursors[index].mz <= 0)
          throw std::invalid_argument("XIPM export invalid precursor m/z");
        maps_.push_back({window, -1, library_.precursors[index].mz, {}, {}, {}, {}});
        for (auto f = library_.fragment_offset.at(index); f < library_.fragment_offset.at(index + 1); ++f)
        {
          const auto& fragment = library_.fragments.at(f);
          if (fragment.transition < 0) throw std::invalid_argument("XIPM export negative transition ID");
          const double factor = run_.calibration.has_ms2 ? 1 + run_.calibration.ms2_offset_ppm * 1e-6 : 1;
          if (!std::isfinite(fragment.product_mz * factor) || fragment.product_mz * factor <= 0)
            throw std::invalid_argument("XIPM export invalid calibrated fragment m/z");
          maps_.push_back({window, fragment.transition, fragment.product_mz * factor, {}, {}, {}, {}});
        }
      }
      const AnalysisLibrary& library_;
      const RunResult& run_;
      const PeakMapExportOptions& options_;
      std::vector<const RunWindow*> pending_;
      std::function<void(const Map&)> emit_;
      std::vector<Map> maps_;
      std::size_t next_{};
      std::size_t points_{};
      double last_rt_{-std::numeric_limits<double>::infinity()};
    };

    std::vector<std::int64_t> nativePrecursorIds(const AnalysisLibrary& library)
    {
      std::vector<std::int64_t> ids(library.precursors.size());
      std::vector<bool> numeric(ids.size(), false);
      std::unordered_set<std::int64_t> used;
      // Reserve all original numeric IDs first, independent of score/q selection.
      for (std::size_t i = 0; i < ids.size(); ++i)
      {
        const auto id = library.precursors[i].id;
        if (id.empty()) continue;
        const auto parsed = std::from_chars(id.data(), id.data() + id.size(), ids[i]);
        numeric[i] = parsed.ec == std::errc{} && parsed.ptr == id.data() + id.size();
        if (numeric[i] && !used.insert(ids[i]).second)
          throw std::invalid_argument("XIPM export duplicate numeric precursor ID: " + std::string(id));
      }
      std::int64_t next_free = 0;
      for (std::size_t i = 0; i < ids.size(); ++i)
      {
        if (numeric[i]) continue;
        ids[i] = static_cast<std::int64_t>(i);
        if (used.contains(ids[i]))
        {
          while (used.contains(next_free)) ++next_free;
          ids[i] = next_free++;
        }
        used.insert(ids[i]);
      }
      return ids;
    }

    std::shared_ptr<arrow::Table> tableFor(const std::shared_ptr<arrow::Schema>& schema, const Map& map,
                                         const AnalysisLibrary& library, const RunResult& run,
                                         const PeakMapExportOptions& options, const std::vector<std::int64_t>& precursor_ids)
    {
      const auto& w = *map.window;
      const auto& p = library.precursors.at(w.precursor);
      const auto id = precursor_ids.at(w.precursor);
      std::vector<std::shared_ptr<arrow::Array>> arrays;
      for (const auto& field : schema->fields())
      {
        const auto& name = field->name();
        std::unique_ptr<arrow::ArrayBuilder> builder = take(arrow::MakeBuilder(field->type()));
        bool assigned = false;
        const auto integer = [&](std::int64_t value) { checked(static_cast<arrow::Int64Builder&>(*builder).Append(value)); assigned = true; };
        const auto real = [&](double value) { checked(static_cast<arrow::DoubleBuilder&>(*builder).Append(value)); assigned = true; };
        const auto text = [&](const std::string& value) { checked(static_cast<arrow::StringBuilder&>(*builder).Append(value)); assigned = true; };
        const auto binary = [&](const std::vector<double>& values)
        {
          checked(static_cast<arrow::BinaryBuilder&>(*builder).Append(reinterpret_cast<const std::uint8_t*>(values.data()), static_cast<std::int32_t>(values.size() * sizeof(double))));
          assigned = true;
        };
        if (name == "RUN_ID") integer(options.run_id);
        else if (name == "SOURCE_FILE") text(run.run);
        else if (name == "MS_LEVEL") integer(map.transition < 0 ? 1 : 2);
        else if (name == "PEAKMAP_TYPE") text(map.transition < 0 ? "precursor" : "transition");
        else if (name == "PRECURSOR_ID") integer(id);
        else if (name == "TRANSITION_ID" && map.transition >= 0) integer(map.transition);
        else if (name == "MODIFIED_SEQUENCE") text(p.sequence);
        else if (name == "PRECURSOR_CHARGE") integer(p.charge);
        else if (name == "PRECURSOR_DECOY" || (name == "PRODUCT_DECOY" && map.transition >= 0)) integer(p.decoy);
        else if (name == "DETECTING_TRANSITION" && map.transition >= 0) integer(1);
        else if (name == "TARGET_MZ") real(map.target);
        else if (name == "TARGET_RT") real(w.peak.rt_apex);
        else if (name == "TARGET_ION_MOBILITY" && run.calibration.has_im && p.library_im > 0) real(run.calibration.im_intercept + run.calibration.im_slope * p.library_im);
        else if (name == "RT_START") real(w.peak.rt_start - options.rt_padding_s);
        else if (name == "RT_END") real(w.peak.rt_end + options.rt_padding_s);
        else if (name == "MZ_DATA") binary(map.mz);
        else if (name == "RT_DATA") binary(map.rt);
        else if (name == "MOBILITY_DATA") binary(map.im);
        else if (name == "INTENSITY_DATA") binary(map.intensity);
        else if (name.ends_with("_COMPRESSION")) integer(0);
        else if (name == "opendia_precursor_id") text(std::string(p.id));
        else if (name == "opendia_library_index") integer(w.precursor);
        else if (name == "candidate_id" || name == "opendia_candidate_row") integer(w.candidate);
        else if (name == "run_file_name") text(File::stemName(run.run));
        else if (name == "feature_id") integer(QPXIdentity::deriveId({QPXIdentity::Value{File::stemName(run.run)},
          QPXIdentity::Value{std::string(p.id)}, QPXIdentity::Value{static_cast<Int64>(w.candidate)}}));
        if (!assigned) checked(builder->AppendNull());
        arrays.push_back(take(builder->Finish()));
      }
      return arrow::Table::Make(schema, arrays);
    }
  }

  void exportPeakMaps(const AnalysisLibrary& library, const RunResult& result,
                      const PeakMapExportOptions& options, const SpectrumReplay& replay)
  {
    if (options.path.empty()) return;
    if (!std::isfinite(options.q_value) || options.q_value < 0 || options.q_value > 1 ||
        !std::isfinite(options.rt_padding_s) || options.rt_padding_s < 0 ||
        options.batch_precursors == 0 || options.batch_precursors > 32 ||
        options.run_id > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) || !replay)
      throw std::invalid_argument("XIPM export requires q in [0,1], finite nonnegative RT padding, batch size in [1,32], signed run ID and spectrum replay");
    const auto& calibration = result.calibration;
    if ((calibration.has_ms2 && (!std::isfinite(calibration.ms2_offset_ppm) ||
          !std::isfinite(calibration.ms2_tolerance_ppm) || calibration.ms2_tolerance_ppm <= 0)) ||
        (calibration.has_im && (!std::isfinite(calibration.im_intercept) || !std::isfinite(calibration.im_slope) ||
          !std::isfinite(calibration.im_tolerance) || calibration.im_tolerance < 0)))
      throw std::invalid_argument("XIPM export invalid calibrated m/z or mobility model");
    const auto precursor_ids = nativePrecursorIds(library);
    std::vector<const RunWindow*> winners(library.precursors.size(), nullptr);
    for (const auto& window : result.windows)
    {
      if (window.precursor >= winners.size()) throw std::invalid_argument("XIPM export invalid precursor index");
      if (!(window.score > CrossFitScorer::NOT_KEPT) || !std::isfinite(window.score)) continue;
      if (!std::isfinite(window.peak.rt_start) || !std::isfinite(window.peak.rt_end) || window.peak.rt_start > window.peak.rt_end)
        throw std::invalid_argument("XIPM export invalid candidate RT bounds");
      auto& best = winners[window.precursor];
      if (!best || window.score > best->score) best = &window;
    }
    std::vector<const RunWindow*> selected;
    for (const auto* winner : winners)
    {
      if (winner && winner->q_value >= 0 && winner->q_value <= options.q_value &&
          (options.include_decoys || !library.precursors[winner->precursor].decoy)) selected.push_back(winner);
    }
    auto fields = XIPMSchema::schema()->fields();
    fields.push_back(arrow::field("opendia_precursor_id", arrow::utf8()));
    fields.push_back(arrow::field("opendia_library_index", arrow::int64()));
    fields.push_back(arrow::field("candidate_id", arrow::int64()));
    fields.push_back(arrow::field("opendia_candidate_row", arrow::int64()));
    fields.push_back(arrow::field("run_file_name", arrow::utf8()));
    fields.push_back(arrow::field("feature_id", arrow::int64()));
    auto schema = arrow::schema(fields);
    static std::atomic<std::uint64_t> serial{};
    const auto temporary = options.path + ".partial." + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + "." + std::to_string(serial++);
    try
    {
      auto sink = take(arrow::io::FileOutputStream::Open(temporary));
      parquet::WriterProperties::Builder properties;
      properties.compression(parquet::Compression::ZSTD);
      auto writer = take(parquet::arrow::FileWriter::Open(*schema, arrow::default_memory_pool(), sink, properties.build()));
      std::vector<std::shared_ptr<arrow::Table>> completed;
      std::size_t completed_points = 0;
      const auto flush = [&]()
      {
        if (completed.empty()) return;
        const auto table = take(arrow::ConcatenateTables(completed));
        checked(writer->WriteTable(*table, table->num_rows()));
        completed.clear();
        completed_points = 0;
      };
      Collector collector(library, result, options, std::move(selected), [&](const Map& map)
      {
        // Bound both encoded binary buffers and row-group metadata. A single large map is written immediately.
        if (completed_points + map.mz.size() > 524288) flush();
        completed.push_back(tableFor(schema, map, library, result, options, precursor_ids));
        completed_points += map.mz.size();
        if (completed.size() >= options.batch_precursors * 7 || completed_points >= 524288) flush();
      });
      replay(collector, {1, 2});
      collector.finish();
      flush();
      checked(writer->Close());
      checked(sink->Close());
      std::filesystem::rename(temporary, options.path);
    }
    catch (...)
    {
      std::error_code ignored;
      std::filesystem::remove(temporary, ignored);
      throw;
    }
  }
}
