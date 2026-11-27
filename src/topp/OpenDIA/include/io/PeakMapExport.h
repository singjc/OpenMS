// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <workflow/Workflow.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <functional>

namespace OpenMS::OpenDIA
{
  struct PeakMapExportOptions
  {
    std::string path;
    double q_value{0.01};
    double rt_padding_s{15.0};
    bool include_decoys{false};
    bool im_array_as_vssc{false};
    std::uint64_t run_id{0};
    std::size_t batch_precursors{32};
  };
  using SpectrumReplay = std::function<void(Interfaces::IMSDataConsumer&, const std::vector<int>&)>;

  /// Export native OpenMS XIPMSchema maps, without an analysis dependency on external assay records.
  /// Called once after the final scoring pass, before the run-local spectrum cache is destroyed.
  void exportPeakMaps(const AnalysisLibrary& library, const RunResult& result,
                      const PeakMapExportOptions& options, const SpectrumReplay& replay);
}
