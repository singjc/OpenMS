// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/**
  FASTA in silico library (PLAN_FASTA): candidates from a FASTA with one pseudo-reversed decoy per target, cheap
  coordinates (PeptDeep RT, CCS) for all of them, a survey library with theoretical fragments for the data-driven
  prefilter, and the final library with predicted fragment intensities for the kept pairs only.
*/
namespace OpenMS::OpenDIA::Preparation
{
  /// FASTA digestion and precursor settings (PLAN_FASTA defaults).
  struct CandidateOptions
  {
    std::string enzyme = "Trypsin/P";
    int missed_cleavages = 1;
    int min_length = 7;
    int max_length = 30;
    int min_charge = 2;
    int max_charge = 4;
    double min_mz = 400.0;
    double max_mz = 1200.0;
    std::vector<std::string> fixed_mods{"Carbamidomethyl (C)"};
    std::vector<std::string> variable_mods{"Oxidation (M)"};
    int max_variable_mods = 1;
  };

  /// One peptidoform: a target or its pseudo-reversed decoy.
  struct Peptidoform
  {
    std::string sequence;            ///< UniMod notation ("PEPC(UniMod:4)K").
    std::uint32_t proteins = 0;      ///< Index into CandidateSet::protein_sets (decoys: DECOY_ prefixed accessions).
    std::uint8_t missed_cleavages = 0;
    std::uint8_t variable_mods = 0;
    std::uint8_t length = 0;
    bool decoy = false;
    float rt = 0.0f;                 ///< Library RT (100 x PeptDeep RT), set by predictCoordinates.
  };

  /// One precursor: a peptidoform at one charge.
  struct CandidatePrecursor
  {
    std::uint32_t peptidoform = 0;
    double mz = 0.0;
    float im = -1.0f;                ///< 1/K0 from predicted CCS, -1 without.
    std::int8_t charge = 0;
  };

  /**
    @brief The candidates of a FASTA library: complete target-decoy pairs only.

    Precursor 2k is the target and 2k + 1 the decoy of pair k (same charge). The precursor ID is the decimal index; the
    TRAML_ID is "<target sequence>/<charge>", with DECOY_ for the decoy, so the analysis pairs them by name.
  */
  struct CandidateSet
  {
    std::vector<Peptidoform> peptidoforms;
    std::vector<CandidatePrecursor> precursors;
    std::vector<std::vector<std::string>> protein_sets;
    std::size_t proteins_read = 0;
    std::size_t target_peptides = 0;     ///< Distinct target peptides (unmodified) with at least one pair.
    std::size_t dropped_no_decoy = 0;    ///< Target peptidoforms without an admissible decoy (terminal modification).
    std::size_t dropped_collision = 0;   ///< Target peptidoforms whose decoy equals a target or an earlier decoy.
    std::size_t skipped_residues = 0;    ///< Peptides with residues outside the 20 standard amino acids.

    std::size_t pairs() const { return precursors.size() / 2; }
    std::string tramlId(std::size_t precursor) const;
  };

  /**
    @brief Digest @p fasta, apply modifications and add one pseudo-reversed decoy (pseudoReverse) per target peptidoform.

    A target whose decoy is impossible or equals a target or an earlier decoy peptidoform is dropped with all its
    charges. A pair is kept at every charge at which the target m/z is in range. Peptides are visited in the order of
    first appearance in the FASTA and list every protein they occur in.
    @throws Exception::InvalidValue for an unknown modification, a terminal one, or one without UniMod ID or formula.
  */
  CandidateSet enumerateCandidates(const std::string& fasta, const CandidateOptions& options,
                                   const std::function<void(const std::string&)>& log = {});

  /// Internal cleavage sites of @p unmodified for @p enzyme.
  int missedCleavages(const std::string& enzyme, const std::string& unmodified);

  /// A theoretical or predicted fragment.
  struct FragmentIon
  {
    double mz = 0.0;
    float intensity = 0.0f;
    std::int16_t ordinal = 0;
    std::int8_t charge = 1;
    char type = 'y';
  };

  /// b and y ion and precursor masses of one peptidoform (cumulative residue masses).
  class IonLadder
  {
  public:
    explicit IonLadder(const AASequence& peptide);
    std::size_t size() const { return cumulative_.size() - 1; }
    double ion(char type, std::size_t ordinal, int charge) const;
    double precursor(int charge) const;

  private:
    std::vector<double> cumulative_;
  };

  /**
    @brief The 6 theoretical survey fragments (PLAN_FASTA step 3): y1+ of ordinal 3 .. n-1 in [@p min_mz, @p max_mz],
    the 6 closest to ordinal 7 (lower first on ties), topped up with b1+ (ordinal 3 .. n-1, closest to 5, lower first).
    Intensity 1, in selection order.
  */
  std::vector<FragmentIon> surveyFragments(const IonLadder& ladder, double min_mz, double max_mz);

  /**
    @brief Final fragments from a PeptDeep MS2 prediction (8 channels per cleavage, b/y at z = 1, 2 first): ordinal >= 3,
    fragment charge <= max(1, min(2, z - 1)), m/z in range, intensity > 0; the top @p count by intensity (ties: lower
    m/z), in that order.
  */
  std::vector<FragmentIon> predictedFragments(const IonLadder& ladder, int precursor_charge, const std::vector<float>& ms2,
                                              double min_mz, double max_mz, std::size_t count);

  /**
    @brief The calibration panel: pairs whose target has no variable modification, charge 2, no missed cleavage and
    length 7-20, at an even stride over target RT (ties: pair index), at most @p count. @return Pair indices.
  */
  std::vector<std::uint32_t> calibrationPanel(const CandidateSet& candidates, std::size_t count);

  /**
    @brief The pairs kept by the prefilter: the @p keep largest @p statistic > 0, ties by @p tie then index, returned
    sorted by index. @p keep = 0 keeps every pair.
  */
  std::vector<std::uint32_t> selectPairs(const std::vector<double>& statistic, const std::vector<std::uint64_t>& tie, std::size_t keep);

  /**
    @brief PLAN_SURVEY +pool: per pair, the largest @p statistic over the pairs with the same target modified sequence
    (charge siblings) plus @p protein_weight x the protein evidence, the largest over the pair's proteins of the mean
    of that protein's top 3 pair statistics. @p sequence and @p proteins are the target's, per pair (label-free
    attributes: the result is a symmetric function of each pair's max(target, decoy)). A pair with an empty @p sequence
    (absent from the survey library) gets 0 and contributes nothing; a present pair with statistic 0 may be lifted by its
    siblings and proteins.
  */
  std::vector<double> pooledStatistic(const std::vector<double>& statistic, const std::vector<std::string>& sequence,
                                      const std::vector<std::vector<std::string>>& proteins, double protein_weight = 0.5);

  /// PeptDeep settings of a FASTA library.
  struct PredictionSettings
  {
    std::string model_directory;   ///< peptdeep_{rt,ccs,ms2}_dynamic.onnx; empty = installed models.
    std::string instrument = "QE";
    float nce = 30.0f;
    bool predict_im = false;
    int threads = 4;
  };

  /// The resolved path of PeptDeep model @p name ("peptdeep_rt_dynamic.onnx") for @p settings.
  std::string modelPath(const PredictionSettings& settings, const std::string& name);

  /// Fragment settings of the survey and final libraries.
  struct FragmentOptions
  {
    double min_mz = 200.0;
    double max_mz = 1800.0;
    std::size_t count = 12;        ///< Final library fragments per precursor (at most 99).
    std::size_t min_count = 3;     ///< A pair is dropped if a member has fewer.
    /// PLAN_SURVEY V1/V3: survey fragments are the top 6 PeptDeep fragments of every pair instead of theoretical
    /// y ions (slower: MS2 for every candidate). Default off.
    bool predicted_survey = false;
  };

  /// Library RT (100 x PeptDeep RT) for every peptidoform and, with @p settings.predict_im, 1/K0 for every precursor.
  void predictCoordinates(CandidateSet& candidates, const PredictionSettings& settings,
                          const std::function<void(const std::string&)>& log = {});

  /**
    @brief Write the survey library: every pair with survey fragments, except the @p panel pairs, which come first and
    have their top 6 predicted fragments (so the identification-free calibration sees intensities). The footer key
    opendia.fasta.panel_precursors holds the number of leading panel precursors.
  */
  void writeSurveyLibrary(const CandidateSet& candidates, const std::vector<std::uint32_t>& panel, const PredictionSettings& settings,
                          const FragmentOptions& fragments, const std::string& path,
                          const std::vector<std::pair<std::string, std::string>>& metadata,
                          const std::function<void(const std::string&)>& log = {});

  /// Leading panel precursors of a survey or final FASTA library (0 if the footer key is missing).
  std::size_t surveyPanelPrecursors(const std::string& path);

  struct FinalLibraryStats
  {
    std::size_t pairs_in = 0, pairs_out = 0, fragments = 0, panel_pairs = 0;
  };

  /**
    @brief Write the final library from a survey library: PeptDeep MS2 for both members of every pair k with
    @p keep(k) (precursor IDs 2k and 2k + 1), FragmentOptions::count fragments into @p path and the top 6 into @p search_path (same
    IDs). Pairs with a member below FragmentOptions::min_count fragments are dropped.
    The kept panel pairs of the survey come first; the footer key opendia.fasta.panel_precursors holds their precursor
    count, so the final search calibrates on them only (as the survey does).
  */
  FinalLibraryStats writeFinalLibrary(const std::string& survey, const std::function<bool(std::size_t)>& keep,
                                      const PredictionSettings& settings, const FragmentOptions& fragments, const std::string& path,
                                      const std::string& search_path, const std::vector<std::pair<std::string, std::string>>& metadata,
                                      const std::function<void(const std::string&)>& log = {});
}
