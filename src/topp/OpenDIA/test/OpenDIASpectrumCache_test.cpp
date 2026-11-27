// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIASpectrumCache.h>

#include <OpenMS/IONMOBILITY/IMDataArrayUtils.h>
#include <OpenMS/METADATA/Precursor.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  void require(bool ok, const std::string& message)
  {
    if (!ok)
    {
      std::cerr << "FAIL: " << message << '\n';
      std::exit(1);
    }
  }

  struct Collect final : Interfaces::IMSDataConsumer
  {
    std::vector<MSSpectrum> spectra;
    void consumeSpectrum(MSSpectrum& s) override { spectra.push_back(s); }
    void consumeChromatogram(MSChromatogram&) override {}
    void setExpectedSize(Size, Size) override {}
    void setExperimentalSettings(const ExperimentalSettings&) override {}
  };

  struct Thrower final : Interfaces::IMSDataConsumer
  {
    std::size_t seen = 0;
    void consumeSpectrum(MSSpectrum&) override
    {
      if (++seen == 3) throw std::runtime_error("consumer stop");
    }
    void consumeChromatogram(MSChromatogram&) override {}
    void setExpectedSize(Size, Size) override {}
    void setExperimentalSettings(const ExperimentalSettings&) override {}
  };

  MSSpectrum make(int level, double rt, std::size_t peaks, bool im_array, bool scalar_im, bool sorted)
  {
    MSSpectrum s;
    s.setMSLevel(level);
    s.setRT(rt);
    for (std::size_t i = 0; i < peaks; ++i)
    {
      Peak1D p;
      // Unsorted spectra keep their order in the cache (the collectors sort them themselves).
      p.setMZ(sorted ? 100.0 + 0.1234567891234 * i : 900.0 - 0.37 * i);
      p.setIntensity(static_cast<float>(1.0 + i * 3.25));
      s.push_back(p);
    }
    if (level == 2)
    {
      Precursor a;
      a.setMZ(500.25 + rt * 1e-3);
      a.setIsolationWindowLowerOffset(12.5);
      a.setIsolationWindowUpperOffset(12.75);
      Precursor b = a;
      b.setMZ(700.125);
      s.setPrecursors({a, b});
    }
    if (im_array)
    {
      DataArrays::FloatDataArray unrelated;
      unrelated.setName("charge-like");
      unrelated.assign(peaks, 7.0f);
      s.getFloatDataArrays().push_back(unrelated);
      DataArrays::FloatDataArray im;
      IMDataArrayUtils::setIMUnit(im, DriftTimeUnit::VSSC);
      for (std::size_t i = 0; i < peaks; ++i) im.push_back(0.8f + 0.001f * static_cast<float>(i));
      s.getFloatDataArrays().push_back(im);
    }
    if (scalar_im)
    {
      s.setDriftTime(1.0375);
      s.setDriftTimeUnit(DriftTimeUnit::VSSC);
    }
    return s;
  }

  void same(const MSSpectrum& a, const MSSpectrum& b, const std::string& what)
  {
    require(a.getRT() == b.getRT() && a.getMSLevel() == b.getMSLevel(), what + ": RT/level");
    require(a.size() == b.size(), what + ": size");
    for (std::size_t i = 0; i < a.size(); ++i)
      require(a[i].getMZ() == b[i].getMZ() && a[i].getIntensity() == b[i].getIntensity(), what + ": peak " + std::to_string(i));
    require(a.isSorted() == b.isSorted(), what + ": sortedness");
    require(a.getPrecursors().size() == b.getPrecursors().size(), what + ": precursor count");
    for (std::size_t i = 0; i < a.getPrecursors().size(); ++i)
    {
      const auto& p = a.getPrecursors()[i];
      const auto& q = b.getPrecursors()[i];
      require(p.getMZ() == q.getMZ() && p.getIsolationWindowLowerOffset() == q.getIsolationWindowLowerOffset() &&
                  p.getIsolationWindowUpperOffset() == q.getIsolationWindowUpperOffset(),
              what + ": precursor");
    }
    require(a.getDriftTime() == b.getDriftTime() && a.getDriftTimeUnit() == b.getDriftTimeUnit(), what + ": drift time");
    require(a.containsIMData() == b.containsIMData(), what + ": has IM");
    if (a.containsIMData())
    {
      const auto [ia, ua] = a.getIMData();
      const auto [ib, ub] = b.getIMData();
      require(ua == ub, what + ": IM unit");
      require(a.getFloatDataArrays()[ia].getName() == b.getFloatDataArrays()[ib].getName(), what + ": IM name");
      require(std::vector<float>(a.getFloatDataArrays()[ia].begin(), a.getFloatDataArrays()[ia].end()) ==
                  std::vector<float>(b.getFloatDataArrays()[ib].begin(), b.getFloatDataArrays()[ib].end()),
              what + ": IM values");
    }
  }
}

int main()
{
  const auto base = std::filesystem::temp_directory_path() / ("opendia_cache_test_" + std::to_string(::getpid()));
  std::vector<MSSpectrum> input;
  input.push_back(make(1, 10.0, 50, false, false, true));
  input.push_back(make(2, 10.5, 40, true, false, true));
  input.push_back(make(2, 10.5, 0, true, false, true));    // empty, same RT
  input.push_back(make(3, 11.0, 5, false, false, true));   // MS3: not cached
  input.push_back(make(2, 11.25, 30, false, true, false)); // scalar IM, unsorted
  input.push_back(make(1, 12.0, 20, true, false, true));
  {
    OpenDIA::SpectrumCache cache((base / "c1").string());
    require(std::filesystem::exists(base / "c1" / "spectra.bin"), "cache file created");
    for (const auto& s : input) cache.append(s);
    cache.finishWriting();
    require(cache.complete() && cache.index().size() == 5, "five MS1/MS2 records");
    require(cache.index()[1].isolation_low == 500.25 + 10.5e-3 - 12.5 && cache.index()[1].isolation_high == 700.125 + 12.75, "isolation span");
    for (std::size_t prefetch : {std::size_t{0}, std::size_t{1}, std::size_t{64}})
    {
      Collect both, ms2;
      cache.replay(both, {1, 2}, prefetch);
      cache.replay(ms2, {2}, prefetch);
      require(both.spectra.size() == 5 && ms2.spectra.size() == 3, "replay counts");
      const std::vector<std::size_t> kept{0, 1, 2, 4, 5};
      for (std::size_t i = 0; i < kept.size(); ++i) same(input[kept[i]], both.spectra[i], "replay " + std::to_string(i));
      same(input[1], ms2.spectra[0], "ms2 0");
      same(input[4], ms2.spectra[2], "ms2 2");
      // Filtered replay keeps the order of the records it passes.
      Collect odd;
      cache.replay(odd, {1, 2}, prefetch, [](const OpenDIA::SpectrumCache::Entry& e) { return e.level == 1 || e.peaks == 0; });
      require(odd.spectra.size() == 3, "filtered replay count");
      same(input[0], odd.spectra[0], "filtered 0");
      same(input[2], odd.spectra[1], "filtered 1");
      same(input[5], odd.spectra[2], "filtered 2");
    }
    // The index lists every record's precursors as stored.
    for (std::size_t r = 0; r < cache.index().size(); ++r)
    {
      const auto& e = cache.index()[r];
      const std::vector<std::size_t> kept{0, 1, 2, 4, 5};
      const auto& s = input[kept[r]];
      require(e.peaks == s.size() && e.precursor_count == s.getPrecursors().size(), "index peaks and precursors " + std::to_string(r));
      for (std::size_t k = 0; k < e.precursor_count; ++k)
      {
        const auto& v = cache.precursors()[e.precursor_begin + k];
        const auto& p = s.getPrecursors()[k];
        require(v[0] == p.getMZ() && v[1] == p.getIsolationWindowLowerOffset() && v[2] == p.getIsolationWindowUpperOffset(),
                "index precursor " + std::to_string(r));
      }
    }
    Thrower t;
    bool threw = false;
    try
    {
      cache.replay(t, {1, 2});
    }
    catch (const std::runtime_error&)
    {
      threw = true;
    }
    require(threw && t.seen == 3, "consumer exception is rethrown");
  }
  require(!std::filesystem::exists(base / "c1"), "cache directory removed on destruction");
  {
    OpenDIA::SpectrumCache cache((base / "c2").string());
    cache.append(input[0]);
    cache.discard();
    require(!cache.usable() && !cache.complete(), "discarded cache is unusable");
    bool threw = false;
    try
    {
      Collect c;
      cache.replay(c, {1});
    }
    catch (const std::logic_error&)
    {
      threw = true;
    }
    require(threw, "replay of an incomplete cache throws");
  }
  require(!std::filesystem::exists(base / "c2"), "discarded cache directory removed");
  std::filesystem::remove_all(base);
  std::cout << "OpenDIASpectrumCache_test passed\n";
  return 0;
}
