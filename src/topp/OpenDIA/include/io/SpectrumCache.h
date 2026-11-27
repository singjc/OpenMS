// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief Run-local spectrum cache: the first raw read writes every MS1 and MS2 spectrum to one flat binary file, later
    passes replay it instead of decoding the raw file again.

    Only the fields the OpenDIA collectors read are kept, at their in-memory precision: RT, MS level, precursor m/z,
    isolation offsets, precursor mobility rectangles (drift time, unit, and lower/upper offsets), peaks
    (m/z double, intensity float), spectrum drift time and its unit, and the first ion-mobility float
    array with its exact name (the name carries the unit). Records are written in arrival order, so a replay delivers
    the same spectra in the same order. An in-memory index (level, RT, isolation window, file offset) allows reading
    the records of selected isolation windows only.

    Replay decodes records on a prefetch thread into a bounded queue while the consumer runs.
    The cache directory is created on construction and removed (with its files) on destruction.
  */
  class SpectrumCache
  {
  public:
    struct Entry
    {
      std::uint64_t offset = 0;      ///< Byte offset of the record in the cache file.
      double rt = 0.0;
      double isolation_low = 0.0;    ///< Lowest isolation window bound over the precursors (MS2; 0 for MS1).
      double isolation_high = 0.0;   ///< Highest isolation window bound over the precursors (MS2; 0 for MS1).
      std::uint32_t peaks = 0;
      std::uint32_t precursor_begin = 0; ///< First entry in precursors() of this record.
      std::uint16_t precursor_count = 0;
      std::uint8_t level = 0;
    };
    /// Precursor m/z and isolation window lower/upper offsets, as stored for one record.
    using Isolation = std::array<double, 3>;

    /// Creates @p directory (and parents). @throws std::runtime_error if it cannot be created or the file not opened.
    explicit SpectrumCache(std::string directory);
    ~SpectrumCache();
    SpectrumCache(const SpectrumCache&) = delete;
    SpectrumCache& operator=(const SpectrumCache&) = delete;

    /// Append one spectrum while writing. @throws std::runtime_error on a write error (the cache is then unusable).
    void append(const MSSpectrum& spectrum);
    /// End writing; the cache becomes replayable. @throws std::runtime_error on a write error.
    void finishWriting();
    /// Mark the cache unusable (e.g. the writing pass was interrupted); its file is removed.
    void discard();

    bool complete() const { return complete_; }
    bool usable() const { return !discarded_; }
    const std::vector<Entry>& index() const { return index_; }
    /// Precursors of all records (Entry::precursor_begin / precursor_count).
    const std::vector<Isolation>& precursors() const { return precursors_; }
    std::uint64_t bytes() const { return bytes_; }
    std::uint64_t peaks() const { return peaks_; }
    const std::string& directory() const { return directory_; }

    /**
      @brief Feed the cached spectra whose MS level is in @p levels to @p consumer, in their original order.
      @param[in] prefetch Records decoded ahead on a separate thread (0 decodes on the calling thread).
      @param[in] keep If set, only records for which it returns true are fed (the order is kept).
      @throws std::logic_error if the cache is not complete; consumer exceptions are rethrown.
    */
    void replay(Interfaces::IMSDataConsumer& consumer, const std::vector<int>& levels, std::size_t prefetch = 64,
                const std::function<bool(const Entry&)>& keep = {}) const;

  private:
    void read_(std::ifstream& in, const Entry& entry, MSSpectrum& out) const;

    std::string directory_, file_;
    std::unique_ptr<std::ofstream> out_;
    std::vector<char> out_buffer_;
    std::vector<Entry> index_;
    std::vector<Isolation> precursors_;
    std::vector<std::string> im_names_;
    std::vector<double> scratch_mz_;
    std::vector<float> scratch_intensity_;
    std::uint64_t bytes_ = 0, peaks_ = 0;
    bool complete_ = false, discarded_ = false;
  };
}
