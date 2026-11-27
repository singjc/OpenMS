// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <workflow/Workflow.h>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /// Reporting cuts only. They never change calibration, extraction or precursor scoring.
  struct ReportingOptions
  {
    double peakgroup_fdr{0.01};
    double peptide_fdr{0.01};
    double protein_fdr{0.01};
  };

  /// Sparse, run-specific entity evidence. Area is the sum of eligible raw MS2 precursor areas, not LFQ.
  struct EntityRun
  {
    std::size_t run{};
    double score{};
    double q_value{1.0};  ///< NaN means that a target/decoy counterpart was not observed, so no estimate is available.
    double area{};
    bool has_area{false};  ///< Distinguishes an absent/ineligible quantity from a measured zero.
  };

  struct QuantifiedEntity
  {
    std::string id;
    std::vector<std::string> proteins;
    bool decoy{false};
    double global_q_value{1.0};
    std::vector<EntityRun> runs;
  };

  struct Quantification
  {
    std::vector<double> global_precursor_q;
    std::vector<QuantifiedEntity> peptides;
    std::vector<QuantifiedEntity> protein_groups;
    std::vector<std::uint32_t> precursor_peptide;  ///< UINT32_MAX for unmappable peptide sequences.
    /// UINT32_MAX when a precursor is shared across resolved protein groups or has no mapping.
    std::vector<std::uint32_t> precursor_protein_group;
  };

  /// No extraction/scorer changes. Experimental level-specific picked target-decoy estimates, before report filtering.
  /// Protein groups use all library target peptide evidence. Shared peptides are excluded from group scores/areas.
  /// Entity areas require precursor cuts. Each matrix additionally applies only its corresponding level cuts.
  Quantification quantify(const AnalysisLibrary& library, const std::vector<RunResult>& runs,
                          const ReportingOptions& options = {});
}
