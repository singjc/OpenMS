// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
#pragma once

#include <library/preparation/Assay.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Types.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA::Preparation
{
  /** @brief A precursor with a caller-assigned ID, unique across the complete library. */
  struct PredictionInput
  {
    AASequence peptide;
    std::string id;
    int charge = 0;
    float nce = 30.0f;
    int64_t instrument_index = 0;
    std::vector<std::string> protein_refs;
  };

  /** @brief Predict owned assay records directly from the low-level PeptDeep models. */
  class AssayPredictor
  {
  public:
    struct Config
    {
      std::string rt_model_path;
      std::string ccs_model_path;
      std::string ms2_model_path;
      int intra_op_threads = 4;
      Size batch_size = 500;
      bool predict_ccs = true;
    };

    /// Construct using default models and inference settings.
    AssayPredictor();
    /** @brief Construct a predictor.
        @param[in] config Model paths and inference settings.
        @throws Exception::Precondition if ONNX support is unavailable.
    */
    explicit AssayPredictor(const Config& config);
    ~AssayPredictor();

    /** @brief Predict compounds, proteins and four b/y channels per cleavage.
        @param[in] precursors Nonempty precursor batch with unique IDs.
        @return Owned predicted assay records.
        @throws Exception::IllegalArgument for empty input, empty or duplicate IDs,
                non-positive charges, short peptides or empty protein IDs.
        @throws Exception::InvalidValue for unexpected inference output shapes.
        @throws Exception::Precondition if ONNX support is unavailable.
    */
    Assays predict(const std::vector<PredictionInput>& precursors);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
}
