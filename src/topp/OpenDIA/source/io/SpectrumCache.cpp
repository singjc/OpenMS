// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <io/SpectrumCache.h>

#include <OpenMS/IONMOBILITY/IMDataArrayUtils.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <condition_variable>
#include <cstring>
#include <deque>
#include <exception>
#include <filesystem>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace OpenMS::OpenDIA
{
  namespace
  {
    constexpr std::uint32_t NO_IM = std::numeric_limits<std::uint32_t>::max();

    /// Fixed record header; each precursor is six doubles (m/z, isolation offsets, drift time, drift offsets)
    /// and one unit byte, followed by peak m/z (doubles), intensities (floats), and mobility values (floats).
    /// This cache is run-local, so no persistent-format compatibility is required.
    struct Header
    {
      double rt;
      double drift_time;
      std::uint32_t peaks;
      std::uint32_t im_name;   ///< Index into the name table, NO_IM without an ion-mobility array.
      std::uint32_t im_values; ///< Length of the ion-mobility array (may differ from peaks).
      std::uint16_t precursors;
      std::uint8_t level;
      std::uint8_t drift_unit;
    };
    static_assert(sizeof(Header) == 32, "packed cache record header");

    template <class T>
    void put(std::vector<char>& out, const T* data, std::size_t count)
    {
      const auto* bytes = reinterpret_cast<const char*>(data);
      out.insert(out.end(), bytes, bytes + count * sizeof(T));
    }
  }

  SpectrumCache::SpectrumCache(std::string directory) : directory_(std::move(directory))
  {
    std::error_code ec;
    if (std::filesystem::exists(directory_, ec)) throw std::runtime_error("SpectrumCache: directory already exists: " + directory_);
    std::filesystem::create_directories(directory_, ec);
    if (ec) throw std::runtime_error("SpectrumCache: cannot create " + directory_ + ": " + ec.message());
    file_ = (std::filesystem::path(directory_) / "spectra.bin").string();
    out_ = std::make_unique<std::ofstream>(file_, std::ios::binary | std::ios::trunc);
    if (!*out_)
    {
      std::filesystem::remove_all(directory_, ec);
      throw std::runtime_error("SpectrumCache: cannot open " + file_);
    }
  }

  SpectrumCache::~SpectrumCache()
  {
    out_.reset();
    std::error_code ec;
    std::filesystem::remove_all(directory_, ec);
  }

  void SpectrumCache::append(const MSSpectrum& spectrum)
  {
    if (complete_ || discarded_) throw std::logic_error("SpectrumCache: append after finishWriting() or discard().");
    const int level = static_cast<int>(spectrum.getMSLevel());
    if (level != 1 && level != 2) return;
    if (spectrum.size() > std::numeric_limits<std::uint32_t>::max() || spectrum.getPrecursors().size() > 65535)
      throw std::runtime_error("SpectrumCache: spectrum too large for the cache record format.");
    Header h{};
    h.rt = spectrum.getRT();
    h.drift_time = spectrum.getDriftTime();
    h.peaks = static_cast<std::uint32_t>(spectrum.size());
    h.level = static_cast<std::uint8_t>(level);
    h.drift_unit = static_cast<std::uint8_t>(spectrum.getDriftTimeUnit());
    h.precursors = static_cast<std::uint16_t>(spectrum.getPrecursors().size());
    h.im_name = NO_IM;
    h.im_values = 0;
    // The first float array that names an ion-mobility unit, exactly as MSSpectrum::getIMData() selects it.
    const DataArrays::FloatDataArray* im = nullptr;
    for (const auto& fda : spectrum.getFloatDataArrays())
    {
      DriftTimeUnit unit;
      if (IMDataArrayUtils::getIMUnit(fda, unit))
      {
        im = &fda;
        break;
      }
    }
    if (im != nullptr)
    {
      const std::string& name = im->getName();
      auto it = std::find(im_names_.begin(), im_names_.end(), name);
      if (it == im_names_.end()) it = im_names_.insert(im_names_.end(), name);
      h.im_name = static_cast<std::uint32_t>(it - im_names_.begin());
      h.im_values = static_cast<std::uint32_t>(im->size());
    }
    Entry e;
    e.offset = bytes_;
    e.rt = h.rt;
    e.level = h.level;
    e.peaks = h.peaks;
    e.precursor_begin = static_cast<std::uint32_t>(precursors_.size());
    e.precursor_count = h.precursors;
    double low = std::numeric_limits<double>::infinity(), high = -low;
    out_buffer_.clear();
    put(out_buffer_, &h, 1);
    for (const auto& p : spectrum.getPrecursors())
    {
      const double v[6] = {p.getMZ(), p.getIsolationWindowLowerOffset(), p.getIsolationWindowUpperOffset(),
                           p.getDriftTime(), p.getDriftTimeWindowLowerOffset(), p.getDriftTimeWindowUpperOffset()};
      const auto unit = static_cast<std::uint8_t>(p.getDriftTimeUnit());
      put(out_buffer_, v, 6);
      put(out_buffer_, &unit, 1);
      precursors_.push_back({v[0], v[1], v[2]});
      low = std::min(low, v[0] - v[1]);
      high = std::max(high, v[0] + v[2]);
    }
    if (level == 2 && h.precursors > 0)
    {
      e.isolation_low = low;
      e.isolation_high = high;
    }
    scratch_mz_.resize(h.peaks);
    scratch_intensity_.resize(h.peaks);
    for (std::size_t i = 0; i < h.peaks; ++i)
    {
      scratch_mz_[i] = spectrum[i].getMZ();
      scratch_intensity_[i] = spectrum[i].getIntensity();
    }
    put(out_buffer_, scratch_mz_.data(), h.peaks);
    put(out_buffer_, scratch_intensity_.data(), h.peaks);
    if (im != nullptr) put(out_buffer_, im->data(), im->size());
    out_->write(out_buffer_.data(), static_cast<std::streamsize>(out_buffer_.size()));
    if (!*out_) throw std::runtime_error("SpectrumCache: write failed for " + file_);
    bytes_ += out_buffer_.size();
    peaks_ += h.peaks;
    index_.push_back(e);
  }

  void SpectrumCache::finishWriting()
  {
    if (discarded_) throw std::logic_error("SpectrumCache: finishWriting() after discard().");
    out_->flush();
    if (!*out_) throw std::runtime_error("SpectrumCache: write failed for " + file_);
    out_.reset();
    std::vector<double>().swap(scratch_mz_);
    std::vector<float>().swap(scratch_intensity_);
    std::vector<char>().swap(out_buffer_);
    complete_ = true;
  }

  void SpectrumCache::discard()
  {
    discarded_ = true;
    complete_ = false;
    out_.reset();
    index_.clear();
    precursors_.clear();
    std::error_code ec;
    std::filesystem::remove(file_, ec);
  }

  void SpectrumCache::read_(std::ifstream& in, const Entry& entry, MSSpectrum& out) const
  {
    in.seekg(static_cast<std::streamoff>(entry.offset));
    Header h;
    in.read(reinterpret_cast<char*>(&h), sizeof(h));
    if (!in) throw std::runtime_error("SpectrumCache: cannot read record header for " + file_);
    std::vector<Precursor> precursors(h.precursors);
    for (auto& precursor : precursors)
    {
      double values[6];
      std::uint8_t unit;
      in.read(reinterpret_cast<char*>(values), sizeof(values));
      in.read(reinterpret_cast<char*>(&unit), sizeof(unit));
      if (!in) throw std::runtime_error("SpectrumCache: cannot read precursor isolation rectangle for " + file_);
      precursor.setMZ(values[0]);
      precursor.setIsolationWindowLowerOffset(values[1]);
      precursor.setIsolationWindowUpperOffset(values[2]);
      precursor.setDriftTime(values[3]);
      precursor.setDriftTimeWindowLowerOffset(values[4]);
      precursor.setDriftTimeWindowUpperOffset(values[5]);
      precursor.setDriftTimeUnit(static_cast<DriftTimeUnit>(unit));
    }
    std::vector<double> mz(h.peaks);
    std::vector<float> intensity(h.peaks);
    in.read(reinterpret_cast<char*>(mz.data()), static_cast<std::streamsize>(mz.size() * sizeof(double)));
    in.read(reinterpret_cast<char*>(intensity.data()), static_cast<std::streamsize>(intensity.size() * sizeof(float)));
    out = MSSpectrum();
    if (h.im_name != NO_IM)
    {
      DataArrays::FloatDataArray im;
      im.setName(im_names_.at(h.im_name));
      im.resize(h.im_values);
      in.read(reinterpret_cast<char*>(im.data()), static_cast<std::streamsize>(im.size() * sizeof(float)));
      out.getFloatDataArrays().push_back(std::move(im));
    }
    if (!in) throw std::runtime_error("SpectrumCache: read failed for " + file_);
    out.setRT(h.rt);
    out.setMSLevel(h.level);
    out.setDriftTime(h.drift_time);
    out.setDriftTimeUnit(static_cast<DriftTimeUnit>(h.drift_unit));
    out.setPrecursors(std::move(precursors));
    out.resize(h.peaks);
    for (std::size_t i = 0; i < h.peaks; ++i)
    {
      out[i].setMZ(mz[i]);
      out[i].setIntensity(intensity[i]);
    }
  }

  void SpectrumCache::replay(Interfaces::IMSDataConsumer& consumer, const std::vector<int>& levels, std::size_t prefetch,
                             const std::function<bool(const Entry&)>& keep) const
  {
    if (!complete_) throw std::logic_error("SpectrumCache: replay() before finishWriting().");
    std::vector<const Entry*> selected;
    for (const auto& e : index_)
      if (std::find(levels.begin(), levels.end(), static_cast<int>(e.level)) != levels.end() && (!keep || keep(e))) selected.push_back(&e);
    consumer.setExpectedSize(selected.size(), 0);
    std::ifstream in(file_, std::ios::binary);
    if (!in) throw std::runtime_error("SpectrumCache: cannot open " + file_);
    if (prefetch == 0)
    {
      MSSpectrum s;
      for (const Entry* e : selected)
      {
        read_(in, *e, s);
        consumer.consumeSpectrum(s);
      }
      return;
    }
    // Producer: decodes records in order into a bounded queue. Consumer: this thread, in the same order.
    std::mutex m;
    std::condition_variable cv;
    std::deque<MSSpectrum> queue;
    bool done = false, stop = false;
    std::exception_ptr producer_error;
    std::thread producer([&]()
    {
      try
      {
        for (const Entry* e : selected)
        {
          MSSpectrum s;
          read_(in, *e, s);
          std::unique_lock lock(m);
          cv.wait(lock, [&] { return stop || queue.size() < prefetch; });
          if (stop) return;
          queue.push_back(std::move(s));
          cv.notify_all();
        }
      }
      catch (...)
      {
        std::lock_guard lock(m);
        producer_error = std::current_exception();
      }
      std::lock_guard lock(m);
      done = true;
      cv.notify_all();
    });
    std::exception_ptr consumer_error;
    try
    {
      for (;;)
      {
        MSSpectrum s;
        {
          std::unique_lock lock(m);
          cv.wait(lock, [&] { return !queue.empty() || done; });
          if (queue.empty()) break;
          s = std::move(queue.front());
          queue.pop_front();
          cv.notify_all();
        }
        consumer.consumeSpectrum(s);
      }
    }
    catch (...)
    {
      consumer_error = std::current_exception();
      std::lock_guard lock(m);
      stop = true;
      cv.notify_all();
    }
    producer.join();
    if (consumer_error) std::rethrow_exception(consumer_error);
    if (producer_error) std::rethrow_exception(producer_error);
  }
}
