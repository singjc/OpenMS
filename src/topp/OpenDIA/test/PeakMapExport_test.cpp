// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#include <io/PeakMapExport.h>
#include <OpenMS/FORMAT/XIPMParquetFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>

using namespace OpenMS;
using namespace OpenMS::OpenDIA;
namespace
{
  void check(bool value, const char* message)
  {
    if (!value) throw std::runtime_error(message);
  }
  MSSpectrum spectrum(unsigned level, double rt, double mz, bool im)
  {
    MSSpectrum s;
    s.setMSLevel(level);
    s.setRT(rt);
    Peak1D peak;
    peak.setMZ(mz);
    peak.setIntensity(123.25f);
    s.push_back(peak);
    if (im)
    {
      s.setDriftTime(1.125);
      s.setDriftTimeUnit(DriftTimeUnit::VSSC);
    }
    OpenMS::Precursor off, on;
    off.setMZ(400);
    off.setIsolationWindowLowerOffset(5);
    off.setIsolationWindowUpperOffset(5);
    on.setMZ(505);
    on.setIsolationWindowLowerOffset(5);
    on.setIsolationWindowUpperOffset(5);
    s.setPrecursors({off, on});
    return s;
  }
}
int main(int argc, char** argv)
{
  const auto directory = std::filesystem::path(argc > 1 ? argv[1] : ".") / "PeakMapExport_fixture";
  std::filesystem::create_directories(directory);
  try
  {
    AnalysisLibrary library;
    library.precursors.resize(3);
    library.precursors[0].id = "42";
    library.precursors[1].id = "DECOY_42";
    library.precursors[2].id = "rejected";
    for (auto& p : library.precursors)
    {
      p.mz = 500;
      p.charge = 2;
      p.sequence = "PEPTIDE";
      p.library_im = 1.0;
    }
    library.precursors[1].decoy = true;
    library.fragment_offset = {0, 1, 2, 3};
    library.fragments = {{300, 1, 701}, {300, 1, 702}, {300, 1, 703}};
    RunResult run;
    run.run = "synthetic.mzML";
    run.calibration.has_ms2 = true;
    run.calibration.ms2_offset_ppm = 20;
    run.calibration.ms2_tolerance_ppm = 5;
    run.calibration.has_im = true;
    run.calibration.im_intercept = 0.125;
    run.calibration.im_tolerance = 0.01;
    for (unsigned i = 0; i < 3; ++i)
    {
      RunWindow w;
      w.precursor = i;
      w.candidate = 123 + i;
      w.score = 2;
      w.q_value = i == 2 ? 0.1 : 0.005;
      w.peak.rt_start = 10;
      w.peak.rt_end = 12;
      w.peak.rt_apex = 11;
      run.windows.push_back(w);
    }
    auto discarded = run.windows[0];
    discarded.precursor = 2;
    discarded.score = CrossFitScorer::NOT_KEPT;
    discarded.q_value = 0;
    run.windows.push_back(discarded);
    const double corrected = 300 * (1 + 20e-6);
    std::vector<MSSpectrum> raw = {
      spectrum(1, 9.0, 500.000123456789, false),
      spectrum(2, 10.1234567890123, corrected, true),
      spectrum(2, 12, corrected + 0.00001, false),
      spectrum(2, 8.9, corrected, true),
      spectrum(2, 11, 300, true)
    };
    auto wrong_im = spectrum(2, 11, corrected, true);
    wrong_im.setDriftTime(1.5);
    raw.push_back(wrong_im);
    auto off = spectrum(2, 11, corrected, true);
    off.getPrecursors().resize(1);
    raw.push_back(off);
    std::stable_sort(raw.begin(), raw.end(), [](const auto& a, const auto& b) { return a.getRT() < b.getRT(); });
    int replays = 0;
    SpectrumReplay replay = [&](Interfaces::IMSDataConsumer& consumer, const std::vector<int>& levels)
    {
      check(levels == std::vector<int>({1, 2}), "replay must request MS1 and MS2");
      ++replays;
      for (auto s : raw) consumer.consumeSpectrum(s);
    };
    PeakMapExportOptions options;
    options.path = (directory / "native.xipm").string();
    options.rt_padding_s = 1;
    options.run_id = 91;
    options.batch_precursors = 1;
    exportPeakMaps(library, run, options, replay);
    std::vector<XIPMParquetFile::XIPMPeakMap> maps;
    XIPMParquetFile(options.path).load(maps);
    std::vector<std::string> columns;
    XIPMParquetFile(options.path).getColumns(columns);
    for (const auto* name : {"feature_id", "opendia_precursor_id", "opendia_library_index", "opendia_candidate_row", "run_file_name"})
      check(std::find(columns.begin(), columns.end(), name) != columns.end(), "QPX linking and original-ID mapping columns");
    check(maps.size() == 2 && replays == 1, "only retained target maps");
    check(maps[0].ms_level == 1 && maps[0].precursor_id == 42 && maps[0].mz.size() == 1, "MS1 precursor metadata");
    check(maps[0].mz[0] == 500.000123456789 && maps[0].rt[0] == 9 && std::isnan(maps[0].ion_mobility[0]), "MS1 exact mz and absent IM");
    check(maps[1].transition_id == 701 && maps[1].mz.size() == 2, "multiple isolation windows, calibrated mz and IM filtering");
    check(maps[1].mz[0] == corrected && maps[1].rt[0] == 10.1234567890123 && maps[1].intensity[0] == 123.25 && maps[1].ion_mobility[0] == 1.125, "lossless native double roundtrip");
    check(std::isnan(maps[1].ion_mobility[1]) && maps[1].source_file == run.run && maps[1].run_id == 91, "missing IM and provenance");
    options.include_decoys = true;
    exportPeakMaps(library, run, options, replay);
    XIPMParquetFile(options.path).load(maps);
    check(maps.size() == 4 && replays == 2 && maps[2].precursor_decoy == 1 && maps[2].precursor_id == 1, "decoy flag and stable nonnumeric ID mapping");

    // The same consumer contract must accept cache replay without changing raw values.
    SpectrumCache cache((directory / "cache").string());
    for (const auto& s : raw) cache.append(s);
    cache.finishWriting();
    options.path = (directory / "cache.xipm").string();
    exportPeakMaps(library, run, options, [&](Interfaces::IMSDataConsumer& c, const std::vector<int>& levels) { cache.replay(c, levels, 0); });
    std::vector<XIPMParquetFile::XIPMPeakMap> cached;
    XIPMParquetFile(options.path).load(cached);
    check(cached.size() == maps.size() && cached[1].mz == maps[1].mz && cached[1].rt == maps[1].rt, "raw/cache replay equivalence");

    const auto bytes = std::filesystem::file_size(options.path);
    bool failed = false;
    try { exportPeakMaps(library, run, options, [](Interfaces::IMSDataConsumer&, const std::vector<int>&) { throw std::runtime_error("synthetic replay failure"); }); }
    catch (const std::runtime_error&) { failed = true; }
    check(failed && std::filesystem::file_size(options.path) == bytes, "failure preserves existing final");
    for (const auto& entry : std::filesystem::directory_iterator(directory)) check(entry.path().filename().string().find(".partial.") == std::string::npos, "partial cleanup");
    auto dense = spectrum(2, 11, corrected, false);
    dense.resize(100001);
    for (auto& peak : dense)
    {
      peak.setMZ(corrected);
      peak.setIntensity(1);
    }
    failed = false;
    try
    {
      exportPeakMaps(library, run, options, [&](Interfaces::IMSDataConsumer& consumer, const std::vector<int>&)
      {
        for (unsigned i = 0; i < 41; ++i) consumer.consumeSpectrum(dense);
      });
    }
    catch (const std::runtime_error& error)
    {
      failed = std::string(error.what()).find("resident point limit") != std::string::npos;
    }
    check(failed && std::filesystem::file_size(options.path) == bytes, "resident point guard preserves existing output");
    for (const auto& entry : std::filesystem::directory_iterator(directory)) check(entry.path().filename().string().find(".partial.") == std::string::npos, "point guard partial cleanup");
    for (int bad = 0; bad < 4; ++bad)
    {
      auto invalid = options;
      if (bad == 0) invalid.q_value = -1;
      if (bad == 1) invalid.rt_padding_s = -1;
      if (bad == 2) invalid.batch_precursors = 0;
      if (bad == 3) invalid.batch_precursors = 33;
      failed = false;
      try { exportPeakMaps(library, run, invalid, replay); } catch (const std::invalid_argument&) { failed = true; }
      check(failed, "invalid export options rejected");
    }
    // Mobility arrays use their declared unit, with the same explicit non-CCS override as Traces.
    auto array_scan = spectrum(2, 11, corrected, false);
    DataArrays::FloatDataArray im_array;
    IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);
    im_array.push_back(1.125f);
    array_scan.getFloatDataArrays().push_back(im_array);
    check(array_scan.containsIMData(), "recognized mobility array name");
    auto export_scan = [&](MSSpectrum scan)
    {
      exportPeakMaps(library, run, options, [&](Interfaces::IMSDataConsumer& c, const std::vector<int>&) { c.consumeSpectrum(scan); });
      XIPMParquetFile(options.path).load(maps);
    };
    options.include_decoys = false;
    options.path = (directory / "im_edges.xipm").string();
    export_scan(array_scan);
    check(maps[1].ion_mobility == std::vector<double>({1.125}), "VSSC mobility array extraction");
    array_scan.getFloatDataArrays()[0][0] = 1.5f;
    export_scan(array_scan);
    check(maps[1].mz.empty(), "off IM peak excluded");
    array_scan.getFloatDataArrays()[0].setName("Ion Mobility");
    export_scan(array_scan);
    check(maps[1].mz.size() == 1 && std::isnan(maps[1].ion_mobility[0]), "unitless mobility unavailable by default");
    options.im_array_as_vssc = true;
    export_scan(array_scan);
    check(maps[1].mz.empty(), "explicit unitless VSSC override filters off IM");
    options.im_array_as_vssc = false;
    auto rectangle = spectrum(2, 11, corrected, true);
    rectangle.getPrecursors()[1].setDriftTimeUnit(DriftTimeUnit::VSSC);
    rectangle.getPrecursors()[1].setDriftTime(1.25);
    rectangle.getPrecursors()[1].setDriftTimeWindowLowerOffset(0.125);
    rectangle.getPrecursors()[1].setDriftTimeWindowUpperOffset(0.0625);
    export_scan(rectangle);
    check(maps[1].mz.size() == 1, "inclusive diaPASEF lower mz and IM edge");
    rectangle.getPrecursors()[1].setDriftTime(1.5);
    export_scan(rectangle);
    check(maps[1].mz.empty(), "off diaPASEF mobility rectangle excluded");
    auto inside_rectangle = rectangle;
    inside_rectangle.getPrecursors()[1].setDriftTime(1.25);
    options.path = (directory / "rectangle_raw.xipm").string();
    exportPeakMaps(library, run, options, [&](Interfaces::IMSDataConsumer& c, const std::vector<int>&)
    {
      c.consumeSpectrum(inside_rectangle);
      c.consumeSpectrum(rectangle);
    });
    std::vector<XIPMParquetFile::XIPMPeakMap> rectangle_raw;
    XIPMParquetFile(options.path).load(rectangle_raw);
    SpectrumCache rectangle_cache((directory / "rectangle_cache").string());
    rectangle_cache.append(inside_rectangle);
    rectangle_cache.append(rectangle);
    rectangle_cache.finishWriting();
    options.path = (directory / "rectangle_cached.xipm").string();
    exportPeakMaps(library, run, options, [&](Interfaces::IMSDataConsumer& c, const std::vector<int>& levels)
    {
      rectangle_cache.replay(c, levels, 2);
    });
    XIPMParquetFile(options.path).load(maps);
    check(maps[1].mz.size() == 1 && maps[1].mz == rectangle_raw[1].mz && maps[1].rt == rectangle_raw[1].rt &&
          maps[1].ion_mobility == rectangle_raw[1].ion_mobility && maps[1].intensity == rectangle_raw[1].intensity,
          "cached scalar IM and precursor rectangles preserve exact raw XIPM selection");
    options.include_decoys = true;
    options.path = (directory / "empty_points.xipm").string();
    exportPeakMaps(library, run, options, [](Interfaces::IMSDataConsumer&, const std::vector<int>&) {});
    XIPMParquetFile(options.path).load(maps);
    check(maps.size() == 4 && maps[0].mz.empty() && maps[1].ion_mobility.empty(), "empty maps reader compatible");
    options.path = (directory / "empty_selection.xipm").string();
    options.q_value = 0;
    exportPeakMaps(library, run, options, replay);
    XIPMParquetFile(options.path).load(maps);
    check(maps.empty(), "empty selection reader compatible");
    AnalysisLibrary collision;
    collision.precursors.resize(2);
    collision.precursors[0] = library.precursors[0];
    collision.precursors[1] = library.precursors[0];
    collision.precursors[0].id = "1";
    collision.precursors[1].id = "nonnumeric";
    collision.fragment_offset = {0, 0, 0};
    RunResult collision_run;
    collision_run.run = run.run;
    for (unsigned i = 0; i < 2; ++i)
    {
      auto w = run.windows[0];
      w.precursor = i;
      collision_run.windows.push_back(w);
    }
    options.path = (directory / "id_collision.xipm").string();
    options.q_value = 0.01;
    exportPeakMaps(collision, collision_run, options, [](Interfaces::IMSDataConsumer&, const std::vector<int>&) {});
    XIPMParquetFile(options.path).load(maps);
    check(maps.size() == 2 && maps[0].precursor_id == 1 && maps[1].precursor_id != 1,
          "fallback library index cannot collide with original numeric precursor ID");
    collision.precursors[1].id = "01";
    failed = false;
    try { exportPeakMaps(collision, collision_run, options, replay); }
    catch (const std::invalid_argument& error) { failed = std::string(error.what()).find("duplicate numeric") != std::string::npos; }
    check(failed, "ambiguous original numeric IDs fail rather than silently merge maps");
    AnalysisLibrary many;
    RunResult many_run;
    many_run.run = run.run;
    many_run.calibration = run.calibration;
    std::vector<std::string> ids;
    ids.reserve(40);
    many.fragment_offset.push_back(0);
    for (unsigned i = 0; i < 40; ++i)
    {
      ids.push_back(std::to_string(1000 + i));
      auto p = library.precursors[0];
      p.id = ids.back();
      many.precursors.push_back(p);
      many.fragments.push_back({300, 1, 10000 + i});
      many.fragment_offset.push_back(i + 1);
      auto w = run.windows[0];
      w.precursor = i;
      w.peak.rt_start += i * 3;
      w.peak.rt_end += i * 3;
      w.peak.rt_apex += i * 3;
      many_run.windows.push_back(w);
    }
    options.path = (directory / "single_replay.xipm").string();
    options.q_value = 0.01;
    unsigned large_replays = 0;
    exportPeakMaps(many, many_run, options, [&](Interfaces::IMSDataConsumer& consumer, const std::vector<int>&)
    {
      ++large_replays;
      for (unsigned i = 0; i < 40; ++i)
      {
        auto scan = spectrum(2, 11 + i * 3, corrected, true);
        consumer.consumeSpectrum(scan);
      }
    });
    XIPMParquetFile(options.path).load(maps);
    check(large_replays == 1 && maps.size() == 80, "more than 32 precursors require exactly one replay");
    for (std::size_t i = 1; i < maps.size(); i += 2) check(maps[i].mz.size() == 1, "RT lifecycle expires and emits each interval");
    failed = false;
    try
    {
      exportPeakMaps(many, many_run, options, [&](Interfaces::IMSDataConsumer& consumer, const std::vector<int>&)
      {
        auto later = spectrum(2, 12, corrected, true);
        auto earlier = spectrum(2, 11, corrected, true);
        consumer.consumeSpectrum(later);
        consumer.consumeSpectrum(earlier);
      });
    }
    catch (const std::runtime_error& error) { failed = std::string(error.what()).find("nondecreasing") != std::string::npos; }
    check(failed, "nonmonotonic replay fails explicitly");
    std::filesystem::remove_all(directory);
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
