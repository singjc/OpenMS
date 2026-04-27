// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/FastaEvidenceFilter.h>

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionListEvidenceFilter.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{
  namespace
  {
    struct IndexPeptideInfo
    {
      std::string modified_peptide_sequence;
      bool decoy{false};
    };

    bool withinSwathWindow_(double precursor_mz, const OpenSwath::SwathMap& map)
    {
      return precursor_mz >= map.lower && precursor_mz <= map.upper;
    }

    bool parseIonAnnotation_(const String& annotation, std::string& product_type, int& ordinal)
    {
      static const std::regex ion_pattern("^([A-Za-z\\.']+)([0-9]+)");
      std::smatch match;
      const std::string annotation_std = annotation.c_str();
      if (!std::regex_search(annotation_std, match, ion_pattern))
      {
        return false;
      }
      product_type = match[1].str();
      ordinal = std::stoi(match[2].str());
      return true;
    }

    double halfTolerance_(double full_window, bool /*ppm*/, double fallback)
    {
      if (full_window > 0.0)
      {
        return std::max(full_window * 0.5, fallback);
      }
      return fallback;
    }
  }

  FastaEvidenceFilter::FastaEvidenceFilter() :
    DefaultParamHandler("FastaEvidenceFilter"),
    ProgressLogger()
  {
    defaults_.setValue("aggregation_method", "any",
                       "How to combine stage-1 evidence across multiple DIA runs.");
    defaults_.setValidStrings("aggregation_method", {"any", "all"});

    TransitionListEvidenceFilter stage1_filter;
    Param stage1_defaults = stage1_filter.getParameters();
    stage1_defaults.remove("enabled");
    defaults_.insert("Stage1:", stage1_defaults);
    defaults_.setValue("Stage1:precursor_batch_size", 50000,
                       "Maximum number of precursors materialized in one stage-1 transition batch.");
    defaults_.setMinInt("Stage1:precursor_batch_size", 1);

    defaults_.setValue("Stage2:mode", "qvalue",
                       "How to accept stage-2 confirmed peptides.");
    defaults_.setValidStrings("Stage2:mode", {"qvalue", "raw_score"});
    defaults_.setValue("Stage2:max_qvalue", 0.01,
                       "Maximum peptide-level q-value in Stage2:mode=qvalue.");
    defaults_.setMinFloat("Stage2:max_qvalue", 0.0);
    defaults_.setMaxFloat("Stage2:max_qvalue", 1.0);
    defaults_.setValue("Stage2:min_matched_ions", 5,
                       "Minimum matched fragment ions required in stage 2.");
    defaults_.setMinInt("Stage2:min_matched_ions", 0);
    defaults_.setValue("Stage2:decoys", "true",
                       "Generate target-decoy peptides for stage-2 q-value estimation.");
    defaults_.setValidStrings("Stage2:decoys", {"true", "false"});
    defaults_.setValue("Stage2:decoy_prefix", "DECOY_",
                       "Prefix added to generated stage-2 decoy protein accessions.");

    defaults_.setValue("Protein:min_confirmed_peptides", 1,
                       "Minimum number of confirmed peptides needed to keep a protein.");
    defaults_.setMinInt("Protein:min_confirmed_peptides", 1);
    defaults_.setValue("Protein:unique_peptides_only", "false",
                       "If true, only peptides mapping to one protein count toward protein support.");
    defaults_.setValidStrings("Protein:unique_peptides_only", {"true", "false"});

    defaults_.setValue("Export:export_fragments", "false",
                       "If true, emit one TSV row per precursor-fragment pair instead of one row per precursor.");
    defaults_.setValidStrings("Export:export_fragments", {"true", "false"});

    std::vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    std::vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);

    defaults_.setValue("SearchSpace:enzyme", "Trypsin", "The enzyme used for in-silico digestion.");
    defaults_.setValidStrings("SearchSpace:enzyme", ListUtils::create<std::string>(all_enzymes));
    defaults_.setValue("SearchSpace:specificity", "full",
                       "Required enzyme specificity at the peptide termini.");
    defaults_.setValidStrings("SearchSpace:specificity", {"full", "semi", "none"});
    defaults_.setValue("SearchSpace:missed_cleavages", 1,
                       "Number of missed cleavages allowed during digestion.");
    defaults_.setMinInt("SearchSpace:missed_cleavages", 0);
    defaults_.setValue("SearchSpace:min_size", 7,
                       "Minimum peptide length after digestion.");
    defaults_.setMinInt("SearchSpace:min_size", 1);
    defaults_.setValue("SearchSpace:max_size", 40,
                       "Maximum peptide length after digestion.");
    defaults_.setMinInt("SearchSpace:max_size", 1);
    defaults_.setValue("SearchSpace:min_mass", 100,
                       "Minimum peptide mass kept in the search space.");
    defaults_.setMinInt("SearchSpace:min_mass", 0);
    defaults_.setValue("SearchSpace:max_mass", 9000,
                       "Maximum peptide mass kept in the search space.");
    defaults_.setMinInt("SearchSpace:max_mass", 0);
    defaults_.setValue("SearchSpace:modifications:fixed",
                       std::vector<std::string>{"Carbamidomethyl (C)"},
                       "Fixed modifications specified with search-engine names.");
    defaults_.setValidStrings("SearchSpace:modifications:fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("SearchSpace:modifications:variable",
                       std::vector<std::string>{"Oxidation (M)"},
                       "Variable modifications specified with search-engine names.");
    defaults_.setValidStrings("SearchSpace:modifications:variable", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("SearchSpace:modifications:variable_max_per_peptide", 2,
                       "Maximum number of variable modifications per peptide candidate.");
    defaults_.setMinInt("SearchSpace:modifications:variable_max_per_peptide", 0);
    defaults_.setValue("SearchSpace:precursor:min_charge", 2,
                       "Minimum precursor charge exported to the peptide search space.");
    defaults_.setMinInt("SearchSpace:precursor:min_charge", 1);
    defaults_.setValue("SearchSpace:precursor:max_charge", 5,
                       "Maximum precursor charge exported to the peptide search space.");
    defaults_.setMinInt("SearchSpace:precursor:max_charge", 1);
    defaults_.setValue("SearchSpace:fragment:min_charge", 1,
                       "Minimum fragment charge generated for each precursor.");
    defaults_.setMinInt("SearchSpace:fragment:min_charge", 1);
    defaults_.setValue("SearchSpace:fragment:max_charge", 2,
                       "Maximum fragment charge generated for each precursor.");
    defaults_.setMinInt("SearchSpace:fragment:max_charge", 1);
    defaults_.setValue("SearchSpace:fragment:max_per_precursor", 6,
                       "Maximum number of theoretical fragments kept per precursor.");
    defaults_.setMinInt("SearchSpace:fragment:max_per_precursor", 1);
    defaults_.setValue("SearchSpace:fragment:min_ion_index", 2,
                       "Skip theoretical fragment ions with ordinal less than or equal to this value.");
    defaults_.setMinInt("SearchSpace:fragment:min_ion_index", 0);

    defaultsToParam_();
    updateMembers_();
  }

  void FastaEvidenceFilter::updateMembers_()
  {
    aggregation_method_ = param_.getValue("aggregation_method").toString();

    stage1_min_supported_precursors_ = static_cast<Size>(param_.getValue("Stage1:min_supported_precursors"));
    stage1_precursor_batch_size_ = static_cast<Size>(param_.getValue("Stage1:precursor_batch_size"));

    stage2_mode_ = param_.getValue("Stage2:mode").toString();
    stage2_max_qvalue_ = static_cast<double>(param_.getValue("Stage2:max_qvalue"));
    stage2_min_matched_ions_ = static_cast<Int>(param_.getValue("Stage2:min_matched_ions"));
    stage2_decoys_ = param_.getValue("Stage2:decoys").toString() == "true";
    stage2_decoy_prefix_ = param_.getValue("Stage2:decoy_prefix").toString();

    protein_min_confirmed_peptides_ = static_cast<Size>(param_.getValue("Protein:min_confirmed_peptides"));
    protein_unique_peptides_only_ = param_.getValue("Protein:unique_peptides_only").toString() == "true";

    export_fragments_ = param_.getValue("Export:export_fragments").toString() == "true";

    enzyme_ = param_.getValue("SearchSpace:enzyme").toString();
    enzyme_specificity_ = param_.getValue("SearchSpace:specificity").toString();
    peptide_missed_cleavages_ = static_cast<Int>(param_.getValue("SearchSpace:missed_cleavages"));
    peptide_min_size_ = static_cast<Int>(param_.getValue("SearchSpace:min_size"));
    peptide_max_size_ = static_cast<Int>(param_.getValue("SearchSpace:max_size"));
    peptide_min_mass_ = static_cast<Int>(param_.getValue("SearchSpace:min_mass"));
    peptide_max_mass_ = static_cast<Int>(param_.getValue("SearchSpace:max_mass"));
    modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("SearchSpace:modifications:fixed"));
    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("SearchSpace:modifications:variable"));
    max_variable_mods_per_peptide_ = static_cast<Int>(param_.getValue("SearchSpace:modifications:variable_max_per_peptide"));
    precursor_min_charge_ = static_cast<Int>(param_.getValue("SearchSpace:precursor:min_charge"));
    precursor_max_charge_ = static_cast<Int>(param_.getValue("SearchSpace:precursor:max_charge"));
    fragment_min_charge_ = static_cast<Int>(param_.getValue("SearchSpace:fragment:min_charge"));
    fragment_max_charge_ = static_cast<Int>(param_.getValue("SearchSpace:fragment:max_charge"));
    max_fragments_per_precursor_ = static_cast<Size>(param_.getValue("SearchSpace:fragment:max_per_precursor"));
    fragment_min_ion_index_ = static_cast<Int>(param_.getValue("SearchSpace:fragment:min_ion_index"));
  }

  std::string FastaEvidenceFilter::makeCanonicalPeptideKey(const std::string& modified_peptide_sequence, int precursor_charge)
  {
    return modified_peptide_sequence + "/" + std::to_string(precursor_charge);
  }

  bool FastaEvidenceFilter::hasDecoyPrefix_(const std::string& value, const std::string& decoy_prefix)
  {
    if (!decoy_prefix.empty() && value.rfind(decoy_prefix, 0) == 0)
    {
      return true;
    }
    return value.rfind("DECOY_", 0) == 0 || value.rfind("Decoy_", 0) == 0 || value.rfind("decoy_", 0) == 0;
  }

  std::string FastaEvidenceFilter::buildInternalKey_(const std::string& modified_peptide_sequence,
                                                     int precursor_charge,
                                                     bool decoy) const
  {
    const std::string canonical_key = makeCanonicalPeptideKey(modified_peptide_sequence, precursor_charge);
    if (!decoy)
    {
      return canonical_key;
    }
    return std::string(stage2_decoy_prefix_.c_str()) + canonical_key;
  }

  std::vector<FastaEvidenceFilter::FragmentRecord> FastaEvidenceFilter::buildTheoreticalFragments_(const AASequence& modified_sequence,
                                                                                                   int precursor_charge) const
  {
    std::vector<FragmentRecord> fragments;

    const int max_fragment_charge = std::min(fragment_max_charge_, precursor_charge);
    if (max_fragment_charge < fragment_min_charge_)
    {
      return fragments;
    }

    TheoreticalSpectrumGenerator spectrum_generator;
    Param tsg_params = spectrum_generator.getParameters();
    tsg_params.setValue("add_metainfo", "true");
    tsg_params.setValue("add_losses", "false");
    tsg_params.setValue("add_precursor_peaks", "false");
    tsg_params.setValue("isotope_model", "none");
    tsg_params.setValue("add_abundant_immonium_ions", "false");
    tsg_params.setValue("add_a_ions", "false");
    tsg_params.setValue("add_b_ions", "true");
    tsg_params.setValue("add_c_ions", "false");
    tsg_params.setValue("add_x_ions", "false");
    tsg_params.setValue("add_y_ions", "true");
    tsg_params.setValue("add_z_ions", "false");
    spectrum_generator.setParameters(tsg_params);

    MSSpectrum theoretical_spectrum;
    spectrum_generator.getSpectrum(theoretical_spectrum, modified_sequence, fragment_min_charge_, max_fragment_charge, precursor_charge);

    if (theoretical_spectrum.getStringDataArrays().empty() || theoretical_spectrum.getIntegerDataArrays().empty())
    {
      return fragments;
    }

    const auto& annotations = theoretical_spectrum.getStringDataArrays()[0];
    const auto& charges = theoretical_spectrum.getIntegerDataArrays()[0];
    const Size n = std::min(theoretical_spectrum.size(), std::min(annotations.size(), charges.size()));
    for (Size i = 0; i < n; ++i)
    {
      std::string product_type;
      int ordinal = 0;
      if (!parseIonAnnotation_(annotations[i], product_type, ordinal))
      {
        continue;
      }
      if (ordinal <= fragment_min_ion_index_)
      {
        continue;
      }

      FragmentRecord fragment;
      fragment.product_mz = theoretical_spectrum[i].getMZ();
      fragment.product_charge = charges[i];
      fragment.product_type = product_type;
      fragment.product_ordinal = ordinal;
      fragments.push_back(std::move(fragment));
    }

    std::sort(fragments.begin(), fragments.end(),
              [](const FragmentRecord& lhs, const FragmentRecord& rhs)
              {
                if (lhs.product_charge != rhs.product_charge) return lhs.product_charge < rhs.product_charge;
                if (lhs.product_type != rhs.product_type) return lhs.product_type < rhs.product_type;
                if (lhs.product_ordinal != rhs.product_ordinal) return lhs.product_ordinal < rhs.product_ordinal;
                return lhs.product_mz < rhs.product_mz;
              });

    if (fragments.size() > max_fragments_per_precursor_)
    {
      fragments.resize(max_fragments_per_precursor_);
    }

    return fragments;
  }

  std::vector<FastaEvidenceFilter::PeptideEntry> FastaEvidenceFilter::generatePeptideEntries(const std::vector<FASTAFile::FASTAEntry>& fasta_entries) const
  {
    return generatePeptideEntries_(fasta_entries, true);
  }

  std::vector<FastaEvidenceFilter::PeptideEntry> FastaEvidenceFilter::generatePeptideEntries_(const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                                                                               bool include_fragments) const
  {
    FragmentIndex fragment_index;
    Param fragment_index_params = fragment_index.getParameters();
    fragment_index_params.setValue("enzyme", enzyme_);
    fragment_index_params.setValue("peptide:enzyme_specificity", enzyme_specificity_);
    fragment_index_params.setValue("peptide:missed_cleavages", peptide_missed_cleavages_);
    fragment_index_params.setValue("peptide:min_size", peptide_min_size_);
    fragment_index_params.setValue("peptide:max_size", peptide_max_size_);
    fragment_index_params.setValue("peptide:min_mass", peptide_min_mass_);
    fragment_index_params.setValue("peptide:max_mass", peptide_max_mass_);
    fragment_index_params.setValue("modifications:fixed", ListUtils::create<std::string>(modifications_fixed_));
    fragment_index_params.setValue("modifications:variable", ListUtils::create<std::string>(modifications_variable_));
    fragment_index_params.setValue("modifications:variable_max_per_peptide", max_variable_mods_per_peptide_);
    fragment_index_params.setValue("precursor:min_charge", precursor_min_charge_);
    fragment_index_params.setValue("precursor:max_charge", precursor_max_charge_);
    fragment_index_params.setValue("fragment:max_charge", fragment_max_charge_);
    fragment_index_params.setValue("fragment:min_ion_index", fragment_min_ion_index_);
    fragment_index_params.setValue("fragment:min_mz", 0);
    fragment_index_params.setValue("fragment:max_mz", 2000);
    fragment_index.setParameters(fragment_index_params);
    fragment_index.buildPeptidesOnly(fasta_entries);

    std::map<std::string, PeptideEntry> peptide_map;
    for (const auto& peptide : fragment_index.getPeptides())
    {
      const AASequence modified_sequence = fragment_index.reconstructModifiedSequence(peptide, fasta_entries);
      const std::string modified_string = modified_sequence.toString().c_str();
      const std::string unmodified_string = modified_sequence.toUnmodifiedString().c_str();
      const std::string protein_ref = fasta_entries[peptide.protein_idx].identifier.c_str();
      const bool decoy = hasDecoyPrefix_(protein_ref, stage2_decoy_prefix_.c_str());

      for (int precursor_charge = precursor_min_charge_; precursor_charge <= precursor_max_charge_; ++precursor_charge)
      {
        const std::string internal_key = buildInternalKey_(modified_string, precursor_charge, decoy);
        auto [it, inserted] = peptide_map.emplace(internal_key, PeptideEntry{});
        PeptideEntry& entry = it->second;
        if (inserted)
        {
          entry.internal_key = internal_key;
          entry.canonical_key = makeCanonicalPeptideKey(modified_string, precursor_charge);
          entry.peptide_sequence = unmodified_string;
          entry.modified_peptide_sequence = modified_string;
          entry.precursor_charge = precursor_charge;
          entry.precursor_mz = modified_sequence.getMZ(precursor_charge);
          if (include_fragments)
          {
            entry.fragments = buildTheoreticalFragments_(modified_sequence, precursor_charge);
          }
          entry.decoy = decoy;
        }
        entry.protein_refs.push_back(protein_ref);
      }
    }

    std::vector<PeptideEntry> peptides;
    peptides.reserve(peptide_map.size());
    for (auto& item : peptide_map)
    {
      auto& refs = item.second.protein_refs;
      std::sort(refs.begin(), refs.end());
      refs.erase(std::unique(refs.begin(), refs.end()), refs.end());
      peptides.push_back(std::move(item.second));
    }

    std::sort(peptides.begin(), peptides.end(),
              [](const PeptideEntry& lhs, const PeptideEntry& rhs)
              {
                if (lhs.precursor_mz != rhs.precursor_mz) return lhs.precursor_mz < rhs.precursor_mz;
                if (lhs.precursor_charge != rhs.precursor_charge) return lhs.precursor_charge < rhs.precursor_charge;
                return lhs.modified_peptide_sequence < rhs.modified_peptide_sequence;
              });

    return peptides;
  }

  OpenSwath::LightTargetedExperiment FastaEvidenceFilter::buildStage1Experiment_(const std::vector<PeptideEntry>& peptides,
                                                                                 Size begin_idx,
                                                                                 Size end_idx) const
  {
    OpenSwath::LightTargetedExperiment experiment;
    end_idx = std::min(end_idx, peptides.size());
    for (Size peptide_idx = begin_idx; peptide_idx < end_idx; ++peptide_idx)
    {
      const auto& peptide = peptides[peptide_idx];
      OpenSwath::LightCompound compound;
      compound.id = peptide.canonical_key;
      compound.sequence = peptide.modified_peptide_sequence;
      compound.charge = peptide.precursor_charge;
      experiment.compounds.push_back(compound);

      const std::vector<FragmentRecord> fragments = peptide.fragments.empty()
        ? buildTheoreticalFragments_(AASequence::fromString(peptide.modified_peptide_sequence), peptide.precursor_charge)
        : peptide.fragments;
      Size intensity_rank = fragments.size();
      for (const auto& fragment : fragments)
      {
        OpenSwath::LightTransition transition;
        transition.transition_name = peptide.canonical_key + "_" + fragment.product_type + std::to_string(fragment.product_ordinal) + "_" + std::to_string(fragment.product_charge);
        transition.peptide_ref = peptide.canonical_key;
        transition.precursor_mz = peptide.precursor_mz;
        transition.product_mz = fragment.product_mz;
        transition.fragment_charge = static_cast<int8_t>(fragment.product_charge);
        transition.fragment_nr = static_cast<int16_t>(fragment.product_ordinal);
        transition.setFragmentType(fragment.product_type);
        transition.setLibraryIntensity(static_cast<double>(intensity_rank--));
        transition.setDetectingTransition(true);
        transition.setQuantifyingTransition(true);
        transition.setDecoy(false);
        experiment.transitions.push_back(std::move(transition));
      }
    }

    return experiment;
  }

  void FastaEvidenceFilter::populateFragments_(std::vector<PeptideEntry>& peptides) const
  {
    for (auto& peptide : peptides)
    {
      if (!peptide.fragments.empty())
      {
        continue;
      }
      peptide.fragments = buildTheoreticalFragments_(AASequence::fromString(peptide.modified_peptide_sequence),
                                                     peptide.precursor_charge);
    }
  }

  std::vector<FASTAFile::FASTAEntry> FastaEvidenceFilter::buildReducedTargetFasta_(const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                                                                   const std::vector<PeptideEntry>& supported_peptides) const
  {
    std::unordered_set<std::string> accessions;
    for (const auto& peptide : supported_peptides)
    {
      accessions.insert(peptide.protein_refs.begin(), peptide.protein_refs.end());
    }

    std::vector<FASTAFile::FASTAEntry> reduced_fasta;
    reduced_fasta.reserve(accessions.size());
    for (const auto& entry : fasta_entries)
    {
      const std::string accession = entry.identifier.c_str();
      if (accessions.find(accession) != accessions.end())
      {
        reduced_fasta.push_back(entry);
      }
    }
    return reduced_fasta;
  }

  std::vector<FASTAFile::FASTAEntry> FastaEvidenceFilter::buildDecoyDatabase_(const std::vector<FASTAFile::FASTAEntry>& target_fasta) const
  {
    std::vector<FASTAFile::FASTAEntry> decoy_fasta;
    decoy_fasta.reserve(target_fasta.size());

    DecoyGenerator decoy_generator;
    for (const auto& entry : target_fasta)
    {
      FASTAFile::FASTAEntry decoy_entry = entry;
      if (enzyme_specificity_ == "none")
      {
        decoy_entry.sequence = decoy_generator.reverseProtein(AASequence::fromString(entry.sequence)).toString();
      }
      else
      {
        decoy_entry.sequence = decoy_generator.reversePeptides(AASequence::fromString(entry.sequence), enzyme_).toString();
      }
      decoy_entry.identifier = stage2_decoy_prefix_ + entry.identifier;
      decoy_fasta.push_back(std::move(decoy_entry));
    }

    return decoy_fasta;
  }

  Param FastaEvidenceFilter::buildFragmentIndexParams_(const ChromExtractParams& ms1_params,
                                                       const ChromExtractParams& ms2_params) const
  {
    FragmentIndex fragment_index;
    Param params = fragment_index.getParameters();

    params.setValue("enzyme", enzyme_);
    params.setValue("peptide:enzyme_specificity", enzyme_specificity_);
    params.setValue("peptide:missed_cleavages", peptide_missed_cleavages_);
    params.setValue("peptide:min_size", peptide_min_size_);
    params.setValue("peptide:max_size", peptide_max_size_);
    params.setValue("peptide:min_mass", peptide_min_mass_);
    params.setValue("peptide:max_mass", peptide_max_mass_);
    params.setValue("modifications:fixed", ListUtils::create<std::string>(modifications_fixed_));
    params.setValue("modifications:variable", ListUtils::create<std::string>(modifications_variable_));
    params.setValue("modifications:variable_max_per_peptide", max_variable_mods_per_peptide_);
    params.setValue("precursor:min_charge", precursor_min_charge_);
    params.setValue("precursor:max_charge", precursor_max_charge_);
    params.setValue("fragment:max_charge", fragment_max_charge_);
    params.setValue("fragment:min_ion_index", fragment_min_ion_index_);
    params.setValue("fragment:min_mz", 0);
    params.setValue("fragment:max_mz", 2000);
    params.setValue("fragment:min_matched_ions", stage2_min_matched_ions_);
    params.setValue("precursor:isotope_error_min", 0);
    params.setValue("precursor:isotope_error_max", 0);
    params.setValue("scoring:max_candidates_per_spectrum", 500);
    params.setValue("decoys", "false");
    params.setValue("snes_enabled", "false");

    params.setValue("precursor:mass_tolerance_unit", ms1_params.ppm ? "ppm" : "Da");
    params.setValue("precursor:mass_tolerance_lower", halfTolerance_(ms1_params.mz_extraction_window, ms1_params.ppm, ms1_params.ppm ? 10.0 : 0.05));
    params.setValue("precursor:mass_tolerance_upper", halfTolerance_(ms1_params.mz_extraction_window, ms1_params.ppm, ms1_params.ppm ? 10.0 : 0.05));
    params.setValue("fragment:mass_tolerance_unit", ms2_params.ppm ? "ppm" : "Da");
    params.setValue("fragment:mass_tolerance", halfTolerance_(ms2_params.mz_extraction_window, ms2_params.ppm, ms2_params.ppm ? 10.0 : 0.05));

    return params;
  }

  FastaEvidenceFilter::Stage2ScoreBundle FastaEvidenceFilter::scoreStage2_(const std::vector<RunData>& runs,
                                                                           const std::vector<FASTAFile::FASTAEntry>& full_fasta,
                                                                           const std::vector<PeptideEntry>& candidates,
                                                                           const ChromExtractParams& ms1_params,
                                                                           const ChromExtractParams& ms2_params) const
  {
    Stage2ScoreBundle bundle;
    if (candidates.empty())
    {
      return bundle;
    }

    FragmentIndex fragment_index;
    Param fragment_index_params = buildFragmentIndexParams_(ms1_params, ms2_params);
    fragment_index.setParameters(fragment_index_params);
    fragment_index.build(full_fasta);

    std::unordered_set<std::string> allowed_keys;
    std::unordered_map<std::string, bool> key_is_decoy;
    allowed_keys.reserve(candidates.size());
    key_is_decoy.reserve(candidates.size());
    for (const auto& candidate : candidates)
    {
      allowed_keys.insert(candidate.internal_key);
      key_is_decoy[candidate.internal_key] = candidate.decoy;
    }

    const auto& index_peptides = fragment_index.getPeptides();
    std::vector<IndexPeptideInfo> peptide_info(index_peptides.size());
    for (Size peptide_index = 0; peptide_index < index_peptides.size(); ++peptide_index)
    {
      const AASequence modified_sequence = fragment_index.reconstructModifiedSequence(index_peptides[peptide_index], full_fasta);
      peptide_info[peptide_index].modified_peptide_sequence = modified_sequence.toString().c_str();
      const std::string accession = full_fasta[index_peptides[peptide_index].protein_idx].identifier.c_str();
      peptide_info[peptide_index].decoy = hasDecoyPrefix_(accession, stage2_decoy_prefix_.c_str());
    }

    for (const auto& run : runs)
    {
      for (const auto& swath_map : run.swath_maps)
      {
        if (swath_map.ms1 || !swath_map.sptr)
        {
          continue;
        }

        std::vector<const PeptideEntry*> map_candidates;
        map_candidates.reserve(candidates.size());
        for (const auto& candidate : candidates)
        {
          if (withinSwathWindow_(candidate.precursor_mz, swath_map))
          {
            map_candidates.push_back(&candidate);
          }
        }

        if (map_candidates.empty())
        {
          continue;
        }

        const size_t n_spectra = swath_map.sptr->getNrSpectra();
        for (size_t spectrum_index = 0; spectrum_index < n_spectra; ++spectrum_index)
        {
          OpenSwath::SpectrumPtr spectrum_ptr;
          if (run.pasef && swath_map.imLower >= 0.0 && swath_map.imUpper >= 0.0)
          {
            spectrum_ptr = swath_map.sptr->getSpectrumById(static_cast<int>(spectrum_index), swath_map.imLower, swath_map.imUpper);
          }
          else
          {
            spectrum_ptr = swath_map.sptr->getSpectrumById(static_cast<int>(spectrum_index));
          }

          if (!spectrum_ptr)
          {
            continue;
          }

          MSSpectrum openms_spectrum;
          OpenSwathDataAccessHelper::convertToOpenMSSpectrum(spectrum_ptr, openms_spectrum);
          openms_spectrum.setMSLevel(2);
          openms_spectrum.sortByPosition();
          if (openms_spectrum.empty())
          {
            continue;
          }

          for (const auto* candidate : map_candidates)
          {
            MSSpectrum pseudo_spectrum = openms_spectrum;
            Precursor precursor;
            precursor.setMZ(candidate->precursor_mz);
            precursor.setCharge(candidate->precursor_charge);
            pseudo_spectrum.setPrecursors({precursor});

            FragmentIndex::SpectrumMatchesTopN matches;
            fragment_index.querySpectrum(pseudo_spectrum, full_fasta, matches);
            for (const auto& match : matches.hits_)
            {
              if (match.peptide_idx_ >= peptide_info.size())
              {
                continue;
              }
              const auto& info = peptide_info[match.peptide_idx_];
              const std::string matched_key = buildInternalKey_(info.modified_peptide_sequence, match.precursor_charge_, info.decoy);
              if (allowed_keys.find(matched_key) == allowed_keys.end())
              {
                continue;
              }
              auto best_it = bundle.best_scores.find(matched_key);
              const double score = static_cast<double>(match.num_matched_);
              if (best_it == bundle.best_scores.end() || score > best_it->second)
              {
                bundle.best_scores[matched_key] = score;
              }
            }
          }
        }
      }
    }

    bundle.score_records.reserve(bundle.best_scores.size());
    for (const auto& score_item : bundle.best_scores)
    {
      bundle.score_records.push_back({score_item.first, score_item.second, key_is_decoy[score_item.first]});
    }
    return bundle;
  }

  std::unordered_map<std::string, double> FastaEvidenceFilter::computePeptideQValues(const std::vector<PeptideScoreRecord>& score_records)
  {
    std::vector<PeptideScoreRecord> sorted_records = score_records;
    std::sort(sorted_records.begin(), sorted_records.end(),
              [](const PeptideScoreRecord& lhs, const PeptideScoreRecord& rhs)
              {
                if (lhs.score != rhs.score) return lhs.score > rhs.score;
                if (lhs.decoy != rhs.decoy) return lhs.decoy && !rhs.decoy;
                return lhs.peptide_key < rhs.peptide_key;
              });

    std::vector<double> threshold_fdr(sorted_records.size(), 1.0);
    Size target_count = 0;
    Size decoy_count = 0;
    for (Size i = 0; i < sorted_records.size(); ++i)
    {
      if (sorted_records[i].decoy) ++decoy_count;
      else ++target_count;
      threshold_fdr[i] = target_count == 0 ? 1.0 : static_cast<double>(decoy_count) / static_cast<double>(target_count);
    }

    for (SignedSize i = static_cast<SignedSize>(threshold_fdr.size()) - 2; i >= 0; --i)
    {
      threshold_fdr[static_cast<Size>(i)] = std::min(threshold_fdr[static_cast<Size>(i)],
                                                     threshold_fdr[static_cast<Size>(i + 1)]);
    }

    std::unordered_map<std::string, double> qvalues;
    qvalues.reserve(sorted_records.size());
    for (Size i = 0; i < sorted_records.size(); ++i)
    {
      qvalues[sorted_records[i].peptide_key] = threshold_fdr[i];
    }
    return qvalues;
  }

  std::vector<FastaEvidenceFilter::PeptideEntry> FastaEvidenceFilter::selectConfirmedPeptides(
    const std::vector<PeptideEntry>& target_peptides,
    const std::unordered_map<std::string, double>& best_scores,
    const std::vector<PeptideScoreRecord>& score_records) const
  {
    std::unordered_map<std::string, double> qvalues;
    if (stage2_mode_ == "qvalue")
    {
      qvalues = computePeptideQValues(score_records);
    }

    std::vector<PeptideEntry> confirmed;
    for (const auto& peptide : target_peptides)
    {
      const auto score_it = best_scores.find(peptide.internal_key);
      if (score_it == best_scores.end() || score_it->second < static_cast<double>(stage2_min_matched_ions_))
      {
        continue;
      }

      bool keep = false;
      if (stage2_mode_ == "raw_score")
      {
        keep = true;
      }
      else
      {
        const auto qvalue_it = qvalues.find(peptide.internal_key);
        keep = qvalue_it != qvalues.end() && qvalue_it->second <= stage2_max_qvalue_;
      }

      if (keep)
      {
        confirmed.push_back(peptide);
      }
    }

    return confirmed;
  }

  std::unordered_set<std::string> FastaEvidenceFilter::selectSupportedProteins(const std::vector<PeptideEntry>& confirmed_peptides,
                                                                               Size min_confirmed_peptides,
                                                                               bool unique_only)
  {
    std::unordered_map<std::string, Size> protein_counts;
    for (const auto& peptide : confirmed_peptides)
    {
      if (unique_only && peptide.protein_refs.size() != 1)
      {
        continue;
      }
      for (const auto& protein_ref : peptide.protein_refs)
      {
        ++protein_counts[protein_ref];
      }
    }

    std::unordered_set<std::string> supported_proteins;
    for (const auto& item : protein_counts)
    {
      if (item.second >= min_confirmed_peptides)
      {
        supported_proteins.insert(item.first);
      }
    }
    return supported_proteins;
  }

  FastaEvidenceFilter::Result FastaEvidenceFilter::filter(const std::vector<RunData>& runs,
                                                          const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                                          const ChromExtractParams& ms1_params,
                                                          const ChromExtractParams& ms2_params,
                                                          int threads) const
  {
    if (runs.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "FastaEvidenceFilter requires at least one DIA run.");
    }
    if (fasta_entries.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "FastaEvidenceFilter requires a non-empty FASTA database.");
    }
    if (stage2_mode_ == "qvalue" && !stage2_decoys_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Stage2:mode=qvalue requires Stage2:decoys=true.");
    }

    const std::vector<PeptideEntry> all_target_peptides = generatePeptideEntries_(fasta_entries, false);
    if (all_target_peptides.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "No peptide precursors were generated from the FASTA database.");
    }

    Param stage1_params = param_.copy("Stage1:", true);
    stage1_params.remove("precursor_batch_size");
    stage1_params.setValue("enabled", "false");

    std::unordered_map<std::string, Size> supported_run_counts;
    const Size batch_size = std::max<Size>(1, stage1_precursor_batch_size_);
    const Size num_batches = (all_target_peptides.size() + batch_size - 1) / batch_size;
    for (Size batch_idx = 0; batch_idx < num_batches; ++batch_idx)
    {
      const Size begin_idx = batch_idx * batch_size;
      const Size end_idx = std::min(begin_idx + batch_size, all_target_peptides.size());
      OPENMS_LOG_INFO << "Stage 1: processing precursor batch " << (batch_idx + 1)
                      << "/" << num_batches << " (" << (end_idx - begin_idx)
                      << " precursors)" << std::endl;
      const OpenSwath::LightTargetedExperiment stage1_experiment =
        buildStage1Experiment_(all_target_peptides, begin_idx, end_idx);

      for (const auto& run : runs)
      {
        TransitionListEvidenceFilter stage1_filter;
        stage1_filter.setParameters(stage1_params);
        stage1_filter.setLogType(getLogType());
        const auto stage1_result = stage1_filter.filter(run.swath_maps, stage1_experiment, ms1_params, ms2_params, run.pasef, threads);

        std::unordered_set<std::string> run_supported;
        for (const auto& compound : stage1_result.filtered_targets.getCompounds())
        {
          run_supported.insert(compound.id);
        }
        for (const auto& compound_id : run_supported)
        {
          ++supported_run_counts[compound_id];
        }
      }
    }

    std::unordered_set<std::string> selected_target_keys;
    for (const auto& item : supported_run_counts)
    {
      if ((aggregation_method_ == "any" && item.second > 0) ||
          (aggregation_method_ == "all" && item.second == runs.size()))
      {
        selected_target_keys.insert(item.first);
      }
    }

    if (selected_target_keys.size() < stage1_min_supported_precursors_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "FastaEvidenceFilter retained only " + String(selected_target_keys.size()) +
                                       " stage-1 precursors, fewer than Stage1:min_supported_precursors=" +
                                       String(stage1_min_supported_precursors_) + ".");
    }

    std::vector<PeptideEntry> selected_target_peptides;
    selected_target_peptides.reserve(selected_target_keys.size());
    for (const auto& peptide : all_target_peptides)
    {
      if (selected_target_keys.find(peptide.canonical_key) != selected_target_keys.end())
      {
        selected_target_peptides.push_back(peptide);
      }
    }

    const std::vector<FASTAFile::FASTAEntry> reduced_target_fasta = buildReducedTargetFasta_(fasta_entries, selected_target_peptides);
    if (reduced_target_fasta.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Stage 1 did not retain any proteins in the input FASTA order.");
    }

    std::vector<PeptideEntry> stage2_candidates = selected_target_peptides;
    std::vector<FASTAFile::FASTAEntry> full_fasta = reduced_target_fasta;
    if (stage2_decoys_)
    {
      std::vector<FASTAFile::FASTAEntry> decoy_fasta = buildDecoyDatabase_(reduced_target_fasta);
      std::vector<PeptideEntry> decoy_peptides = generatePeptideEntries_(decoy_fasta, false);
      stage2_candidates.insert(stage2_candidates.end(), decoy_peptides.begin(), decoy_peptides.end());
      full_fasta.insert(full_fasta.end(), decoy_fasta.begin(), decoy_fasta.end());
    }

    const Stage2ScoreBundle stage2_scores = scoreStage2_(runs, full_fasta, stage2_candidates, ms1_params, ms2_params);
    std::vector<PeptideEntry> confirmed_target_peptides = selectConfirmedPeptides(selected_target_peptides,
                                                                                  stage2_scores.best_scores,
                                                                                  stage2_scores.score_records);

    const std::unordered_set<std::string> supported_proteins = selectSupportedProteins(confirmed_target_peptides,
                                                                                        protein_min_confirmed_peptides_,
                                                                                        protein_unique_peptides_only_);

    std::vector<FASTAFile::FASTAEntry> filtered_fasta;
    filtered_fasta.reserve(supported_proteins.size());
    for (const auto& entry : fasta_entries)
    {
      const std::string accession = entry.identifier.c_str();
      if (supported_proteins.find(accession) != supported_proteins.end())
      {
        filtered_fasta.push_back(entry);
      }
    }

    if (export_fragments_)
    {
      populateFragments_(confirmed_target_peptides);
    }

    Result result;
    result.filtered_fasta = std::move(filtered_fasta);
    result.confirmed_peptides = confirmed_target_peptides;
    result.stage1_supported_precursors = selected_target_peptides.size();
    result.stage2_confirmed_precursors = confirmed_target_peptides.size();
    result.retained_proteins = result.filtered_fasta.size();
    result.summary = "FastaEvidenceFilter retained " + String(result.stage2_confirmed_precursors) +
                     " peptide precursors and " + String(result.retained_proteins) +
                     " proteins after stage 2.";
    return result;
  }
}
