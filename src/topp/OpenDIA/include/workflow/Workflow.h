// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <calibration/Calibration.h>
#include <scoring/Scorer.h>
#include <extraction/Search.h>

#include <library/Library.h>
#include <io/SpectrumCache.h>

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief The library of one OpenDIA analysis: precursors with their target-decoy pairing, peptide family, fold, the
    search fragments (up to 6) and the scoring fragments (up to 12).

    Decoys are named DECOY_<target TRAML_ID>. A precursor's family is its parent target's unmodified sequence; the fold
    is FNV-1a 64 of the family mod 3, so a target, its decoy and all charge states and modified forms of one peptide
    share a fold.
  */
  struct AnalysisLibrary
  {
    struct Precursor
    {
      std::string_view id;     ///< PQP precursor ID (views the key in index).
      std::string sequence;    ///< Modified sequence.
      std::uint32_t proteins = 0;  ///< Index into protein_sets.
      int charge = 0;
      double mz = 0.0;
      double library_rt = 0.0;
      double library_im = -1.0;  ///< Library 1/K0 (uncalibrated), -1 if absent.
      bool has_library_rt = false;
      bool decoy = false;
      bool trace_safe = false;   ///< ID and family are non-empty and free of tabs and line breaks.
      std::uint8_t fold = 0;     ///< FNV-1a of the family mod 3 (the target's fold for a pair).
      std::uint8_t rt_fold = 0;  ///< FNV-1a of the family mod 5 (RT refinement folds).
      std::uint32_t pair = 0;    ///< Dense pair index.
      std::uint32_t decoy_partner = UINT32_MAX;  ///< First precursor named DECOY_<TRAML_ID>, UINT32_MAX if none.
      std::uint32_t sibling_group = 0;           ///< Same decoy flag and modified sequence.
      std::uint64_t family_seed = 0;             ///< FNV-1a of "v" + family.
      std::uint64_t pair_seed = 0;               ///< FNV-1a of "v" + pair key (the target's TRAML_ID).
      std::array<double, 12> scoring_mz{};        ///< Library m/z of the scoring fragments (0 = absent).
      std::array<double, 12> scoring_library{};   ///< Relative to the precursor's most intense scoring fragment.
      std::array<bool, 12> scoring_in_search{};   ///< Also one of the search fragments (same annotation).
      std::uint8_t scoring_count = 0;
    };

    /// One search fragment (library m/z and intensity, PQP transition ID).
    struct Fragment
    {
      double product_mz = 0.0;
      double library_intensity = 0.0;
      std::int64_t transition = 0;
    };

    /// Hash of IDs for heterogeneous (std::string_view) lookup.
    struct IdHash
    {
      using is_transparent = void;
      std::size_t operator()(std::string_view s) const noexcept { return std::hash<std::string_view>{}(s); }
    };

    std::vector<Precursor> precursors;
    /// Search fragments of precursor i: fragments[fragment_offset[i], fragment_offset[i + 1]), in library order.
    std::vector<Fragment> fragments;
    std::vector<std::uint64_t> fragment_offset;
    /// Distinct protein columns of precursors.tsv (references joined by ';').
    std::vector<std::string> protein_sets;
    /// PQP precursor ID -> index. Owns the IDs (nodes are stable, so Precursor::id stays valid when moved).
    std::unordered_map<std::string, std::uint32_t, IdHash, std::equal_to<>> index;
    /// Calibration candidates: (target, DECOY_ partner or UINT32_MAX) in precursor order, as CalibrationCollector pairs them.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> calibration_pairs;
    std::vector<std::uint32_t> sibling_group_sizes;
    std::size_t pairs = 0;

    AnalysisLibrary() = default;
    AnalysisLibrary(AnalysisLibrary&&) = default;
    AnalysisLibrary& operator=(AnalysisLibrary&&) = default;
    AnalysisLibrary(const AnalysisLibrary&) = delete;  ///< Precursor::id views index.
    AnalysisLibrary& operator=(const AnalysisLibrary&) = delete;

    /// Index of precursor @p id, UINT32_MAX if unknown.
    std::uint32_t find(std::string_view id) const
    {
      const auto it = index.find(id);
      return it == index.end() ? UINT32_MAX : it->second;
    }

    /**
      @brief Read the search library and, if given, the scoring library (same precursor IDs; up to 12 fragments each,
      in transition ID order). Without a scoring library the search fragments are the scoring fragments.

      Either may be a PQP or an OpenDIA Parquet library; a PQP is converted to Parquet once (visitAnyLibrary), and
      @p log receives what was read or written. Parquet inputs and valid caches are streamed into this table, with no
      resident Library. The first PQP read still materializes a Library while writing its automatic cache.
      @throws Exception::InvalidValue On a decoy TRAML_ID without the DECOY_ prefix, a transition ID that is not a
      canonical decimal integer, or fragments of one precursor with different precursor m/z or decoy flags.
    */
    static AnalysisLibrary read(const std::string& search_pqp, const std::string& scoring_pqp = "",
                                const std::function<void(const std::string&)>& log = {});
  };

  /// Settings of one run.
  struct RunOptions
  {
    CalibrationOptions calibration;
    std::size_t window_panel_targets = 300000;  ///< Windowed second calibration pass (0 disables).
    std::size_t max_peaks = 5;
    double search_valley_fraction = 0.5;        ///< XicPickerOptions::valley_fraction of the candidate search.
    /// Multiplies the calibrated 1/K0 half width used by the search, the grid, the trace tensor and the interference
    /// resolver (PLAN_IMNARROW D pass: 0.6; 1 = as calibrated). The written calibration JSON keeps the calibrated value.
    double im_tolerance_scale = 0.6;
    double trace_tolerance_ppm = 10.0;  ///< Trace-tensor fragment half width in ppm (PLAN_TRPPM D pass: 10). Must be positive.
    int gbm_rounds = 3;                 ///< GBM rescore rounds (PRECHECK_GBMFIX; default 3 = PLAN_GBT).
    double gbm_positive_q = 0.01;       ///< GBM positives: target winners at this pooled q (default 0.01 = PLAN_GBT).
    std::size_t min_transitions = 1;
    std::size_t population_min_support = 3;     ///< Pairs whose best candidate has fewer supporting fragments are dropped.
    double grid_half_width_s = 15.0;
    double grid_extraction_ppm = 17.5;
    double byion_ppm = 10.0;
    double byion_intensity_min = 300.0;
    double core_half_width_s = 4.0;
    double anchor_apex_s = 5.0;                 ///< An anchor labels the window of its precursor with apex within this.
    std::size_t decoy_seeds = 6000;
    std::uint64_t seed = 0;
    int threads = 1;
    bool im_array_as_vssc = false;
    std::size_t trace_tensor_gb = 64;  ///< Hard cap on the 6-fragment trace tensor (raw floats), in GiB.
    /// Per-run ridge RT correction and a second search/extract/score pass (PLAN_RTR D pass; -no_rt_refine turns it off).
    bool rt_refine = true;
    std::size_t rt_refine_min_ids = 500;  ///< Fewer training IDs skip the refinement.
    /// Parent directory of a run-local spectrum cache: the first raw read writes MS1/MS2 there and the later passes
    /// replay it (same spectra, same order; outputs unchanged). Removed after the run. Empty = read the raw file each pass.
    std::string spectrum_cache_dir;
    /// Extraction with a spectrum cache runs per group of precursors (contiguous in m/z, about equal windows each); each
    /// group replays only the cached spectra that can isolate its precursors. Per-window results are unchanged, and
    /// extraction memory follows the groups in flight. 0 = auto (16). Without a cache extraction is one group.
    std::size_t extraction_groups = 0;
    /// Groups extracted concurrently, each with threads / this many threads. 0 = auto (threads / 4, at least 1).
    std::size_t extraction_parallel = 0;
  };

  /// One scored candidate window of the run (the population), with its features and result.
  struct RunWindow
  {
    std::uint32_t precursor = 0;  ///< AnalysisLibrary::precursors index.
    std::uint32_t candidate = 0;  ///< Search candidate row.
    SearchCandidate peak;
    float score = CrossFitScorer::NOT_KEPT;
    double q_value = 1.0;
    float network_score = CrossFitScorer::NOT_KEPT;  ///< Calibrated network score before the GBM rescore (diagnostic).
    double network_q_value = 1.0;                     ///< Pooled paired TDC q of network_score (diagnostic).
  };

  struct RunResult
  {
    std::string run;
    CalibrationRT calibration;
    std::string calibration_json;  ///< Written calibration.
    std::vector<RunWindow> windows;
    std::vector<std::string> feature_names;
    std::vector<float> features;  ///< windows x feature_names, row-major.
    std::vector<std::string> log;
  };

  /**
    @brief The OpenDIA analysis of one DIA run, from raw spectra to scored precursors.

    1. Calibration (OpenDIACalibration): library RT, fragment m/z and 1/K0 to run coordinates and tolerances,
       identification-free; with the windowed second pass when the whole-run pass finds too few anchors.
    2. Candidate search (StreamSearch) with the calibrated library: up to max_peaks peaks per precursor.
    3. Population: every candidate of every pair whose best candidate has >= population_min_support fragments.
    4. One pass over the MS2 spectra collects, per population window, the 6-fragment trace tensor and the apex b/y ion
       counts (TraceCollector) and the 12-fragment grid (FragmentGridCollector).
    5. Features (24): search (4), boundary (5), b/y ions (2), 12-fragment evidence (6), fragment deficit (3), shared-signal
       interference (4).
    6. CrossFitScorer: three folds, seeds = calibration anchor windows (positives), their decoys' best windows and a
       random sample of decoy windows (negatives); q-values by pooled paired target-decoy competition.

    Steps 1, 2 and 4 each read the spectra once (calibration may take up to three passes).
  */
  class RunAnalysis
  {
  public:
    using Log = std::function<void(const std::string&)>;
    RunAnalysis(AnalysisLibrary& library, RunOptions options, Log log = {});

    /// Analyse one run (raw, mzML, ...). The library is calibrated for each run on the fly (it is not changed).
    RunResult run(const std::string& path);

  private:
    /// Steps 2-6 into @p result (calibration set). @p rt_correction (seconds per precursor, empty = none) is added to
    /// the calibrated library RT.
    void searchAndScore_(const std::string& path, RunResult& result,
                         const std::vector<std::pair<std::string, double>>& anchors, const std::vector<double>& rt_correction,
                         double rt_tolerance_s, const Log& say, std::chrono::steady_clock::time_point t_run);
    /// Fit the ridge RT correction on @p train (library precursor indices, targets @p y in seconds, CV folds @p folds).
    /// Returns the per-precursor correction and the narrowed RT tolerance, or an empty correction when there are fewer
    /// than rt_refine_min_ids training rows.
    std::pair<std::vector<double>, double> refineRT_(const std::vector<std::size_t>& train, const std::vector<double>& y,
                                                     const std::vector<std::uint32_t>& folds, double tolerance_s,
                                                     const std::string& label, const Log& say);

    AnalysisLibrary& library_;
    RunOptions options_;
    Log log_;
    std::unique_ptr<SpectrumCache> cache_; ///< Only during run().
  };

  /// Feature names of RunAnalysis before the sibling block, in column order.
  std::vector<std::string> runFeatureNames();
}
