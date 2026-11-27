// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
#pragma once

#include <limits>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA::Preparation
{
  /// Library-building records, independent of analysis-time library storage and PQP I/O types.
  struct Modification
  {
    int location{};
    int unimod_id{};
    bool operator==(const Modification&) const = default;
  };

  struct Compound
  {
    double drift_time{-1.0};
    double rt{std::numeric_limits<double>::quiet_NaN()};
    double rt_start{std::numeric_limits<double>::quiet_NaN()};
    double rt_end{std::numeric_limits<double>::quiet_NaN()};
    int charge{};
    std::string sequence;
    std::vector<std::string> protein_refs;
    std::string peptide_group_label;
    std::string gene_name;
    std::string id;
    std::string sum_formula;
    std::string compound_name;
    std::string label_type;
    std::string smiles;
    std::string adducts;
    std::vector<Modification> modifications;
    bool isPeptide() const { return compound_name.empty(); }
  };

  struct Transition
  {
    std::string transition_name;
    std::string peptide_ref;
    double library_intensity{};
    double product_mz{};
    double precursor_mz{};
    double precursor_im{-1.0};
    int fragment_charge{};
    int fragment_nr{-1};
    std::string fragment_type;
    std::vector<std::string> peptidoforms;
    bool decoy{false};
    bool detecting{true};
    bool quantifying{true};
    bool identifying{false};
  };

  struct Protein
  {
    std::string id;
    std::string sequence;
    std::string uniprot_id;
  };

  struct Assays
  {
    std::vector<Transition> transitions;
    std::vector<Compound> compounds;
    std::vector<Protein> proteins;
  };
}
