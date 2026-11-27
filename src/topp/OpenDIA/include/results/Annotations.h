// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <results/QpxExport.h>

namespace OpenMS::OpenDIA
{
  /// Optional library gene annotations keyed by producer precursor ID. Missing input annotations stay absent.
  GeneAnnotations readGeneAnnotations(const std::string& library_path);
}
