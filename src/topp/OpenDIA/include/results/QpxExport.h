// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <results/Quantification.h>
#include <map>
#include <string>

namespace OpenMS::OpenDIA
{
  using GeneAnnotations = std::map<std::string, std::string>;

  /// Stable producer identity for a candidate peak group, independent of float32 RT rounding.
  std::int64_t outputFeatureId(const std::string& run, const std::string& precursor, std::uint32_t candidate);

  /// Stream full-precision per-candidate scorer diagnostics before RunResult::features is released.
  void writeCandidateEvidence(const std::string& directory, const AnalysisLibrary& library, const RunResult& run);

  /// Canonical feature/pg, run/sample and provenance views with all candidate/decoy evidence.
  /// Producer-specific canonical columns use QPX x_opendia_* nullable strings. Exact typed diagnostics live separately.
  /// Feature IDs are supplied opaque IDs. Their actual producer composite is recorded independently of QPX's default.
  /// Missing observed precursor m/z and PEP remain null. No synthetic PSM or measured-fragment view is invented.
  void writeQpxDataset(const std::string& directory, const AnalysisLibrary& library,
                       const std::vector<RunResult>& runs, const Quantification& quantification,
                       const GeneAnnotations& genes, const ReportingOptions& options);
}
