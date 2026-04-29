// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/SYSTEM/File.h>

#include <memory>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{
  class AASequence;

  /**
    @brief Prefilter a FASTA database by DIA raw-data evidence.

    The algorithm digests a FASTA database into a theoretical peptide/fragment
    space, uses TransitionListEvidenceFilter as a fast stage-1 DIA prefilter,
    and confirms the stage-1 survivors with a second pass based on
    FragmentIndex::querySpectrum over the same DIA MS2 spectra.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI FastaEvidenceFilter :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// One theoretical product ion exported for a precursor.
    struct FragmentRecord
    {
      double product_mz{0.0};
      int product_charge{0};
      std::string product_type;
      int product_ordinal{0};
    };

    /// One peptide precursor candidate in the generated search space.
    struct PeptideEntry
    {
      std::string internal_key;
      std::string canonical_key;
      std::string peptide_sequence;
      std::string modified_peptide_sequence;
      double precursor_mz{0.0};
      int precursor_charge{0};
      std::vector<std::string> protein_refs;
      std::map<std::string, std::string> protein_gene_names_by_accession;
      std::vector<FragmentRecord> fragments;
      bool decoy{false};
    };

    /// One stage-2 score record used for target-decoy ranking and q-value estimation.
    struct PeptideScoreRecord
    {
      std::string peptide_key;
      double score{0.0};
      bool decoy{false};
    };

    /// One scored stage-2 candidate exported for score-distribution inspection.
    struct Stage2CandidateScore
    {
      std::string peptide_key;
      std::string peptide_sequence;
      std::string modified_peptide_sequence;
      double precursor_mz{0.0};
      int precursor_charge{0};
      std::vector<std::string> protein_refs;
      std::map<std::string, std::string> protein_gene_names_by_accession;
      bool decoy{false};
      std::string source_file;
      std::string native_spectrum_id;
      Size best_matched_ions{0};
      Size supporting_spectra{0};
      Size supporting_runs{0};
      double best_spectrum_matched_intensity_fraction{0.0};
      Size best_spectrum_matched_b_ions{0};
      Size best_spectrum_matched_y_ions{0};
      Size best_spectrum_longest_b_run{0};
      Size best_spectrum_longest_y_run{0};
      double best_spectrum_longest_y_pct{0.0};
      double best_spectrum_poisson_proxy{0.0};
      double best_spectrum_score{0.0};
      double top_run_score_1{0.0};
      double top_run_score_2{0.0};
      double top_run_score_3{0.0};
      double composite_score{0.0};
      double qvalue{-1.0};
      bool accepted{false};
    };

    /// One DIA run loaded into OpenSWATH map containers.
    struct RunData
    {
      std::vector<OpenSwath::SwathMap> swath_maps;
      std::vector<String> swath_map_sources;
      bool pasef{false};
      /// Keeps per-run cached mzML temp files alive for lightClone()-based readers.
      std::shared_ptr<File::TempDir> cache_dir_guard;
    };

    /// Full filtering result for one invocation.
    struct Result
    {
      std::vector<FASTAFile::FASTAEntry> filtered_fasta;
      std::vector<PeptideEntry> confirmed_peptides;
      std::vector<Stage2CandidateScore> stage2_candidate_scores;
      Size stage1_supported_precursors{0};
      Size stage2_confirmed_precursors{0};
      Size retained_proteins{0};
      String summary;
    };

    /** @name Constructors and Destructors
    */
    //@{
    /// Default constructor
    FastaEvidenceFilter();

    /// Destructor
    ~FastaEvidenceFilter() override = default;
    //@}

    /// Synchronize members with the parameter object.
    void updateMembers_() override;

    /**
      @brief Run the full two-stage FASTA evidence filter.

      @param[in] runs Loaded DIA runs
      @param[in] fasta_entries Input FASTA entries
      @param[in] ms1_params MS1 extraction parameters used for stage 1 and stage-2 precursor tolerances
      @param[in] ms2_params MS2 extraction parameters used for stage 1 and stage-2 fragment tolerances
      @param[in] threads Number of threads for stage-1 raw evidence filtering
      @return Filtered proteins and confirmed peptide precursors
    */
    Result filter(const std::vector<RunData>& runs,
                  const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                  const ChromExtractParams& ms1_params,
                  const ChromExtractParams& ms2_params,
                  int threads = 1) const;

    /**
      @brief Digest a FASTA database into unique peptide precursors.

      The peptide entries are collapsed by modified peptide sequence plus
      precursor charge, while preserving many-to-many protein mappings.

      @param[in] fasta_entries FASTA entries to digest
      @return Unique peptide precursor entries sorted by precursor m/z
    */
    std::vector<PeptideEntry> generatePeptideEntries(const std::vector<FASTAFile::FASTAEntry>& fasta_entries) const;

    /**
      @brief Build the canonical target key used throughout the workflow.

      @param[in] modified_peptide_sequence Modified peptide sequence string
      @param[in] precursor_charge Precursor charge state
      @return Canonical key in the form "<modified_sequence>/<charge>"
    */
    static std::string makeCanonicalPeptideKey(const std::string& modified_peptide_sequence, int precursor_charge);

    /**
      @brief Estimate monotone peptide-level q-values from best target/decoy scores.

      @param[in] score_records One best-score record per peptide key
      @param[in] total_target_candidates Total number of target peptide candidates in the stage-2 library
      @param[in] total_decoy_candidates Total number of decoy peptide candidates in the stage-2 library
      @return Mapping from peptide key to empirical q-value
    */
    static std::unordered_map<std::string, double> computePeptideQValues(const std::vector<PeptideScoreRecord>& score_records,
                                                                         Size total_target_candidates,
                                                                         Size total_decoy_candidates);

    /**
      @brief Select confirmed target peptides under the configured stage-2 rule.

      @param[in] target_peptides Target peptide entries from stage 1
      @param[in] best_matched_ions Best stage-2 matched-ion counts keyed by peptide key
      @param[in] score_records Target and decoy best-score records
      @param[in] total_decoy_candidates Total number of stage-2 decoy candidates
      @return Confirmed target peptides
    */
    std::vector<PeptideEntry> selectConfirmedPeptides(const std::vector<PeptideEntry>& target_peptides,
                                                      const std::unordered_map<std::string, Size>& best_matched_ions,
                                                      const std::vector<PeptideScoreRecord>& score_records,
                                                      Size total_decoy_candidates) const;

    /**
      @brief Collapse confirmed peptides to supported proteins.

      @param[in] confirmed_peptides Confirmed target peptides
      @param[in] min_confirmed_peptides Minimal number of peptides needed per protein
      @param[in] unique_only Whether shared peptides should be ignored
      @return Retained protein accessions
    */
    static std::unordered_set<std::string> selectSupportedProteins(const std::vector<PeptideEntry>& confirmed_peptides,
                                                                   Size min_confirmed_peptides,
                                                                   bool unique_only);

private:
    struct Stage2ScoreBundle
    {
      std::unordered_map<std::string, Size> best_matched_ions;
      std::unordered_map<std::string, Stage2CandidateScore> candidate_scores;
      std::vector<PeptideScoreRecord> score_records;
    };

    std::vector<PeptideEntry> generatePeptideEntries_(const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                                      bool include_fragments) const;

    std::vector<FragmentRecord> buildTheoreticalFragments_(const AASequence& modified_sequence,
                                                           int precursor_charge) const;

    OpenSwath::LightTargetedExperiment buildStage1Experiment_(const std::vector<PeptideEntry>& peptides,
                                                              Size begin_idx,
                                                              Size end_idx,
                                                              int threads) const;

    void populateFragments_(std::vector<PeptideEntry>& peptides) const;

    std::vector<FASTAFile::FASTAEntry> buildReducedTargetFasta_(const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                                                const std::vector<PeptideEntry>& supported_peptides) const;

    std::vector<PeptideEntry> buildStage2DecoyPeptides_(const std::vector<PeptideEntry>& target_peptides) const;

    Stage2ScoreBundle scoreStage2_(const std::vector<RunData>& runs,
                                   const std::vector<PeptideEntry>& candidates,
                                   const ChromExtractParams& ms1_params,
                                   const ChromExtractParams& ms2_params,
                                   int threads) const;

    Param buildFragmentIndexParams_(const ChromExtractParams& ms1_params,
                                    const ChromExtractParams& ms2_params) const;

    std::string buildInternalKey_(const std::string& modified_peptide_sequence,
                                  int precursor_charge,
                                  bool decoy) const;

    static bool hasDecoyPrefix_(const std::string& value, const std::string& decoy_prefix);

    String aggregation_method_{"any"};
    Size stage1_min_supported_precursors_{1};
    String stage2_mode_{"qvalue"};
    double stage2_max_qvalue_{0.01};
    Int stage2_min_matched_ions_{5};
    bool stage2_decoys_{true};
    String stage2_decoy_prefix_{"DECOY_"};
    Size protein_min_confirmed_peptides_{1};
    bool protein_unique_peptides_only_{false};
    bool export_fragments_{false};
    bool export_stage2_scores_{false};
    String enzyme_{"Trypsin"};
    String enzyme_specificity_{"full"};
    Int peptide_missed_cleavages_{1};
    Int peptide_min_size_{7};
    Int peptide_max_size_{40};
    Int peptide_min_mass_{100};
    Int peptide_max_mass_{9000};
    StringList modifications_fixed_;
    StringList modifications_variable_;
    Int max_variable_mods_per_peptide_{2};
    Int precursor_min_charge_{2};
    Int precursor_max_charge_{5};
    Int fragment_min_charge_{1};
    Int fragment_max_charge_{2};
    Size max_fragments_per_precursor_{6};
    Int fragment_min_ion_index_{2};
    Size stage1_precursor_batch_size_{50000};
  };
}
