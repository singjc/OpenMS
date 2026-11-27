// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <results/QpxExport.h>

namespace OpenMS::OpenDIA
{
  /// Full winner table (including decoys and pass flags) plus target-only, level-filtered raw MS2 area matrices.
  /// format is tsv, parquet, or both. Matrix missing values are null/NA, never imputed zeros.
  void writeReports(const std::string& directory, const AnalysisLibrary& library,
                    const std::vector<RunResult>& runs, const Quantification& quantification,
                    const GeneAnnotations& genes, const ReportingOptions& options,
                    const std::string& format = "tsv");
}
