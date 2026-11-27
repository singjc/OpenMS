// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>
#include <workflow/Workflow.h>

namespace OpenMS::OpenDIA
{
/// Reporting cuts only. They never change calibration, extraction or precursor scoring.
struct ReportingOptions
{
  double peakgroup_fdr {0.01};
  double peptide_fdr {0.01};
  double protein_fdr {0.01};
  std::string inference_context {"run-specific"};
};

/// Sparse, run-specific entity evidence. Area is the sum of eligible raw MS2 precursor areas, not LFQ.
struct EntityRun
{
  std::size_t run {};
  double score {};
  double q_value {std::numeric_limits<double>::quiet_NaN()}; ///< Run-specific q-value, NaN when inference is unavailable.
  double area {};
  bool has_area {false}; ///< Distinguishes an absent/ineligible quantity from a measured zero.
  double p_value {std::numeric_limits<double>::quiet_NaN()};
  double pep {std::numeric_limits<double>::quiet_NaN()};
  double experiment_q_value {std::numeric_limits<double>::quiet_NaN()};
  double experiment_p_value {std::numeric_limits<double>::quiet_NaN()};
  double experiment_pep {std::numeric_limits<double>::quiet_NaN()};
};

struct QuantifiedEntity
{
  std::string id;
  std::vector<std::string> proteins;
  bool decoy {false};
  double global_q_value {1.0};
  std::vector<EntityRun> runs;
  double global_p_value {std::numeric_limits<double>::quiet_NaN()};
  double global_pep {std::numeric_limits<double>::quiet_NaN()};
};

struct Quantification
{
  std::vector<double> global_precursor_q;
  std::vector<QuantifiedEntity> peptides;
  std::vector<QuantifiedEntity> protein_groups;
  std::vector<std::uint32_t> precursor_peptide; ///< UINT32_MAX for unmappable peptide sequences.
  /// UINT32_MAX when the modified peptide maps to multiple accessions or has no mapping.
  std::vector<std::uint32_t> precursor_protein_group;
  std::vector<std::string> inference_warnings;
};

/// Native default PyProphet level/context inference, before report filtering.
/// Modified sequences collapse charge states. Only uniquely mapped peptides score singleton proteins.
/// Entity areas require precursor cuts. Each matrix additionally applies only its corresponding level cuts.
Quantification quantify(const AnalysisLibrary& library, const std::vector<RunResult>& runs, const ReportingOptions& options = {});
} // namespace OpenMS::OpenDIA
