// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
#include <library/preparation/AssayPrediction.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/config.h>

#ifdef WITH_ONNX
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepCCSInference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepRTInference.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h>

#include <unordered_set>
#include <utility>
#endif

namespace OpenMS::OpenDIA::Preparation
{
#ifdef WITH_ONNX
  namespace
  {
    constexpr Size PEPTDEEP_MS2_CHANNELS = 8;

    std::string resolveModelPath_(const std::string& configured_path, const char* model_name)
    {
      if (!configured_path.empty())
      {
        return configured_path;
      }
      // Prefer downloaded build-tree models, then the resolved OpenMS data path.
      return File::find(std::string("models/") + model_name,
                        {std::string(OPENMS_BINARY_PATH) + "/share/OpenMS"});
    }

    void addModifications_(const AASequence& peptide, Compound& compound)
    {
      if (peptide.hasNTerminalModification())
      {
        compound.modifications.push_back({-1, peptide.getNTerminalModification()->getUniModRecordId()});
      }
      for (Size i = 0; i < peptide.size(); ++i)
      {
        if (peptide[i].isModified())
        {
          compound.modifications.push_back(
            {static_cast<int>(i), peptide.getResidue(i).getModification()->getUniModRecordId()});
        }
      }
      if (peptide.hasCTerminalModification())
      {
        compound.modifications.push_back(
          {static_cast<int>(peptide.size()), peptide.getCTerminalModification()->getUniModRecordId()});
      }
    }
  }

  struct AssayPredictor::Impl
  {
    PeptDeepRTInference rt_predictor;
    std::unique_ptr<PeptDeepCCSInference> ccs_predictor;
    PeptDeepMS2Inference ms2_predictor;

    explicit Impl(const Config& config) :
      rt_predictor(resolveModelPath_(config.rt_model_path, "peptdeep_rt_dynamic.onnx"), config.intra_op_threads, config.batch_size),
      ms2_predictor(resolveModelPath_(config.ms2_model_path, "peptdeep_ms2_dynamic.onnx"), config.intra_op_threads, config.batch_size)
    {
      if (config.predict_ccs)
      {
        ccs_predictor = std::make_unique<PeptDeepCCSInference>(
          resolveModelPath_(config.ccs_model_path, "peptdeep_ccs_dynamic.onnx"), config.intra_op_threads, config.batch_size);
      }
    }
  };
#else
  struct AssayPredictor::Impl {};
#endif

  AssayPredictor::AssayPredictor() : AssayPredictor(Config{})
  {
  }

  AssayPredictor::AssayPredictor(const Config& config)
#ifdef WITH_ONNX
    : impl_(std::make_unique<Impl>(config))
#endif
  {
#ifndef WITH_ONNX
    (void)config;
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Predicted library preparation requires an OpenMS build configured with WITH_ONNX=ON.");
#endif
  }

  AssayPredictor::~AssayPredictor() = default;

  Assays AssayPredictor::predict(const std::vector<PredictionInput>& precursors)
  {
#ifndef WITH_ONNX
    (void)precursors;
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Predicted library preparation requires an OpenMS build configured with WITH_ONNX=ON.");
#else
    if (precursors.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor batch cannot be empty.");
    }

    std::vector<std::string> peptide_strings;
    std::vector<float> charges;
    std::vector<float> nces;
    std::vector<int64_t> instrument_indices;
    peptide_strings.reserve(precursors.size());
    charges.reserve(precursors.size());
    nces.reserve(precursors.size());
    instrument_indices.reserve(precursors.size());

    Size transition_count = 0;
    std::unordered_set<std::string> precursor_ids;
    precursor_ids.reserve(precursors.size());
    for (const auto& precursor : precursors)
    {
      if (precursor.id.empty())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor ID cannot be empty.");
      }
      if (!precursor_ids.insert(precursor.id).second)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Precursor IDs must be unique within a prediction batch: " + precursor.id);
      }
      if (precursor.charge <= 0)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor charge must be positive.");
      }
      if (precursor.peptide.size() < 2)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Peptide sequences must contain at least two residues.");
      }
      for (const auto& protein_ref : precursor.protein_refs)
      {
        if (protein_ref.empty())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Protein IDs cannot be empty.");
        }
      }
      peptide_strings.push_back(precursor.peptide.toString());
      charges.push_back(static_cast<float>(precursor.charge));
      nces.push_back(precursor.nce);
      instrument_indices.push_back(precursor.instrument_index);
      transition_count += (precursor.peptide.size() - 1) * 4;
    }

    const std::vector<float> predicted_rt = impl_->rt_predictor.predictRT(peptide_strings);
    const std::vector<std::vector<float>> predicted_ms2 =
      impl_->ms2_predictor.predictMS2(peptide_strings, charges, nces, instrument_indices);
    std::vector<float> predicted_ccs;
    if (impl_->ccs_predictor)
    {
      predicted_ccs = impl_->ccs_predictor->predictCCS(peptide_strings, charges);
    }
    if (predicted_rt.size() != precursors.size() || predicted_ms2.size() != precursors.size() ||
        (impl_->ccs_predictor && predicted_ccs.size() != precursors.size()))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "prediction count mismatch", "PeptDeep returned an unexpected number of predictions.");
    }

    Assays experiment;
    experiment.compounds.reserve(precursors.size());
    experiment.transitions.reserve(transition_count);
    std::unordered_set<std::string> protein_ids;
    for (Size precursor_index = 0; precursor_index < precursors.size(); ++precursor_index)
    {
      const auto& precursor = precursors[precursor_index];
      const Size cleavage_count = precursor.peptide.size() - 1;
      const Size expected_ms2_values = cleavage_count * PEPTDEEP_MS2_CHANNELS;
      if (predicted_ms2[precursor_index].size() != expected_ms2_values)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      std::to_string(predicted_ms2[precursor_index].size()),
                                      "Expected PeptDeep MS2 output size " + std::to_string(expected_ms2_values) +
                                      " for precursor " + precursor.id + ".");
      }
      const double precursor_mz = precursor.peptide.getMZ(precursor.charge);
      double precursor_im = -1.0;
      if (impl_->ccs_predictor)
      {
        precursor_im = IMTypes::ccsToOneOverK0(predicted_ccs[precursor_index], precursor_mz, precursor.charge);
      }
      Compound compound;
      compound.id = precursor.id;
      compound.sequence = precursor.peptide.toUniModString();
      compound.charge = precursor.charge;
      compound.rt = predicted_rt[precursor_index];
      compound.drift_time = precursor_im;
      compound.protein_refs = precursor.protein_refs;
      addModifications_(precursor.peptide, compound);
      experiment.compounds.push_back(std::move(compound));
      for (const auto& protein_ref : precursor.protein_refs)
      {
        if (protein_ids.insert(protein_ref).second)
        {
          Protein protein;
          protein.id = protein_ref;
          experiment.proteins.push_back(std::move(protein));
        }
      }

      const auto& ms2 = predicted_ms2[precursor_index];
      for (Size row = 0; row < cleavage_count; ++row)
      {
        const Size cleavage = row + 1;
        const Size y_ordinal = precursor.peptide.size() - cleavage;
        const AASequence b_fragment = precursor.peptide.getPrefix(cleavage);
        const AASequence y_fragment = precursor.peptide.getSuffix(y_ordinal);
        const Size offset = row * PEPTDEEP_MS2_CHANNELS;
        auto add_transition = [&](const std::string& fragment_type, Size ordinal, int fragment_charge,
                                  double product_mz, float intensity)
        {
          Transition transition;
          transition.transition_name = precursor.id + "_" + fragment_type + std::to_string(ordinal) +
                                       "^" + std::to_string(fragment_charge);
          transition.peptide_ref = precursor.id;
          transition.library_intensity = intensity;
          transition.product_mz = product_mz;
          transition.precursor_mz = precursor_mz;
          transition.precursor_im = precursor_im;
          transition.fragment_charge = fragment_charge;
          transition.fragment_nr = static_cast<int>(ordinal);
          transition.fragment_type = fragment_type;
          experiment.transitions.push_back(std::move(transition));
        };
        // The remaining four channels are modification-loss ions and are not emitted.
        add_transition("b", cleavage, 1, b_fragment.getMZ(1, Residue::BIon), ms2[offset]);
        add_transition("b", cleavage, 2, b_fragment.getMZ(2, Residue::BIon), ms2[offset + 1]);
        add_transition("y", y_ordinal, 1, y_fragment.getMZ(1, Residue::YIon), ms2[offset + 2]);
        add_transition("y", y_ordinal, 2, y_fragment.getMZ(2, Residue::YIon), ms2[offset + 3]);
      }
    }
    return experiment;
#endif
  }
}
