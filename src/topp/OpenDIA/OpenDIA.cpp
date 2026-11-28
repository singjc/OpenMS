// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <workflow/Workflow.h>
#include <workflow/FastaWorkflow.h>
#include <library/preparation/Preparation.h>
#include <results/Annotations.h>
#include <results/QpxExport.h>
#include <results/ReportExport.h>

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/SysInfo.h>

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <unistd.h>

using namespace OpenMS;

namespace
{
  namespace Preparation = OpenMS::OpenDIA::Preparation;

  // The PQP reader/writer is the only preparation boundary using external assay records.
  template <typename To, typename From>
  void moveCompoundMetadata(To& to, From& from)
  {
    to.id = std::move(from.id);
    to.sequence = std::move(from.sequence);
    to.charge = from.charge;
    to.rt = from.rt;
    to.rt_start = from.rt_start;
    to.rt_end = from.rt_end;
    to.drift_time = from.drift_time;
    to.protein_refs = std::move(from.protein_refs);
    to.peptide_group_label = std::move(from.peptide_group_label);
    to.gene_name = std::move(from.gene_name);
    to.sum_formula = std::move(from.sum_formula);
    to.compound_name = std::move(from.compound_name);
    to.label_type = std::move(from.label_type);
    to.smiles = std::move(from.smiles);
    to.adducts = std::move(from.adducts);
    to.modifications.reserve(from.modifications.size());
    for (const auto& mod : from.modifications) to.modifications.push_back({mod.location, mod.unimod_id});
  }

  template <typename To, typename From>
  void moveTransitionMetadata(To& to, From& from)
  {
    to.transition_name = std::move(from.transition_name);
    to.peptide_ref = std::move(from.peptide_ref);
    to.library_intensity = from.library_intensity;
    to.product_mz = from.product_mz;
    to.precursor_mz = from.precursor_mz;
    to.precursor_im = from.precursor_im;
    to.fragment_charge = from.fragment_charge;
    to.fragment_nr = from.fragment_nr;
    to.peptidoforms = std::move(from.peptidoforms);
  }

  Preparation::Assays fromPQP(OpenSwath::LightTargetedExperiment&& source)
  {
    Preparation::Assays result;
    result.compounds.reserve(source.compounds.size());
    result.transitions.reserve(source.transitions.size());
    result.proteins.reserve(source.proteins.size());
    for (auto& c : source.compounds)
    {
      Preparation::Compound compound;
      moveCompoundMetadata(compound, c);
      result.compounds.push_back(std::move(compound));
    }
    for (auto& t : source.transitions)
    {
      Preparation::Transition transition;
      moveTransitionMetadata(transition, t);
      transition.fragment_type = t.getFragmentType();
      transition.decoy = t.getDecoy();
      transition.detecting = t.isDetectingTransition();
      transition.quantifying = t.isQuantifyingTransition();
      transition.identifying = t.isIdentifyingTransition();
      result.transitions.push_back(std::move(transition));
    }
    for (auto& p : source.proteins)
      result.proteins.push_back({std::move(p.id), std::move(p.sequence), std::move(p.uniprot_id)});
    return result;
  }

  OpenSwath::LightTargetedExperiment toPQP(Preparation::Assays&& source)
  {
    OpenSwath::LightTargetedExperiment result;
    result.compounds.reserve(source.compounds.size());
    result.transitions.reserve(source.transitions.size());
    result.proteins.reserve(source.proteins.size());
    for (auto& c : source.compounds)
    {
      OpenSwath::LightCompound compound;
      moveCompoundMetadata(compound, c);
      result.compounds.push_back(std::move(compound));
    }
    for (auto& t : source.transitions)
    {
      OpenSwath::LightTransition transition;
      moveTransitionMetadata(transition, t);
      transition.setFragmentType(t.fragment_type);
      transition.setDecoy(t.decoy);
      transition.setDetectingTransition(t.detecting);
      transition.setQuantifyingTransition(t.quantifying);
      transition.setIdentifyingTransition(t.identifying);
      result.transitions.push_back(std::move(transition));
    }
    for (auto& p : source.proteins)
    {
      OpenSwath::LightProtein protein;
      protein.id = std::move(p.id);
      protein.sequence = std::move(p.sequence);
      protein.uniprot_id = std::move(p.uniprot_id);
      result.proteins.push_back(std::move(protein));
    }
    return result;
  }
}

/**
  @page TOPP_OpenDIA OpenDIA
  @brief Complete OpenDIA analysis of DIA runs: calibration, candidate search, features, scoring and q-values.

  For every input run: identification-free calibration of library RT, fragment m/z and 1/K0; candidate search with the
  calibrated library; one more pass over the spectra for the scoring evidence (6-fragment traces and apex b/y ions,
  12-fragment grids, shared-signal interference); a cross-fitted, self-trained scorer per run; and q-values by pooled
  paired target-decoy competition. Across runs, precursors get a global q-value from their best run score.

  Output directory: qpx/ holds the canonical QPX evidence dataset. precursors.tsv (or Parquet) contains every
  retained precursor winner across runs. Precursor, peptide and protein-group matrices apply the configured
  run/global q-value cuts. Calibration JSON remains a sidecar. Optional XIPM files contain measured peak maps.
  Legacy TSVs are available under legacy/ for exact comparisons.
*/
class TOPPOpenDIA final : public TOPPBase
{
public:
  TOPPOpenDIA() : TOPPBase("OpenDIA", "Complete OpenDIA analysis of DIA runs: raw spectra to precursor q-values.", {}, false) {}

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", {}, "DIA runs (Thermo raw, mzML, ...)");
    registerInputFile_("lib", "<file>", "", "Library (PQP or OpenDIA .parquet) with decoys named DECOY_<target TRAML_ID> (search fragments, "
                       "up to 6 per precursor). A PQP is converted once to <stem>.opendia.parquet next to it, which later runs read instead. "
                       "Either -lib or -fasta", false);
    setValidFormats_("lib", {"pqp", "parquet"});
    registerInputFile_("fasta", "<file>", "", "Instead of -lib: build one fixed library from this FASTA for all runs (PLAN_FASTA). "
                       "In silico digestion with one pseudo-reversed decoy per target peptidoform, PeptDeep RT (and CCS), a "
                       "survey of every run with theoretical fragments that keeps the -prefilter_pairs best target-decoy pairs (by "
                       "the larger co-elution of target and decoy, so both are treated alike), then PeptDeep fragment intensities "
                       "for the kept pairs only. Cached in -library_cache. Needs a WITH_ONNX build", false);
    setValidFormats_("fasta", {"fasta"});
    registerStringOption_("enzyme", "<name>", "Trypsin/P", "FASTA library: digestion enzyme (OpenMS enzyme name)", false);
    registerIntOption_("missed_cleavages", "<n>", 1, "FASTA library: maximum missed cleavages", false);
    setMinInt_("missed_cleavages", 0);
    registerIntOption_("min_length", "<n>", 7, "FASTA library: minimum peptide length", false);
    setMinInt_("min_length", 4);
    registerIntOption_("max_length", "<n>", 30, "FASTA library: maximum peptide length", false);
    setMinInt_("max_length", 4);
    registerIntOption_("min_charge", "<z>", 2, "FASTA library: minimum precursor charge", false);
    setMinInt_("min_charge", 1);
    registerIntOption_("max_charge", "<z>", 4, "FASTA library: maximum precursor charge", false);
    setMinInt_("max_charge", 1);
    registerDoubleOption_("min_mz", "<m/z>", 400.0, "FASTA library: minimum precursor m/z", false);
    registerDoubleOption_("max_mz", "<m/z>", 1200.0, "FASTA library: maximum precursor m/z", false);
    registerStringList_("fixed_mods", "<mods>", {"Carbamidomethyl (C)"}, "FASTA library: fixed modifications (UniMod names, residue modifications only)", false);
    registerStringList_("variable_mods", "<mods>", {"Oxidation (M)"}, "FASTA library: variable modifications (UniMod names, residue "
                        "modifications only: a terminal modification leaves no pseudo-reversed decoy)", false);
    registerIntOption_("max_variable_mods", "<n>", 1, "FASTA library: maximum variable modifications per peptide", false);
    setMinInt_("max_variable_mods", 0);
    registerIntOption_("fragments", "<n>", 12, "FASTA library: predicted fragments per precursor in the library (b/y, ordinal >= 3); "
                       "the search uses the top 6, scoring the top 12", false);
    setMinInt_("fragments", 3);
    setMaxInt_("fragments", 99);
    registerDoubleOption_("fragment_min_mz", "<m/z>", 200.0, "FASTA library: minimum fragment m/z", false);
    registerDoubleOption_("fragment_max_mz", "<m/z>", 1800.0, "FASTA library: maximum fragment m/z", false);
    registerIntOption_("prefilter_pairs", "<n>", 800000, "FASTA library: target-decoy pairs kept by the survey prefilter and predicted "
                       "(0 = no prefilter, every candidate pair is predicted; default K* of DESIGN_K)", false);
    setMinInt_("prefilter_pairs", 0);
    registerStringOption_("instrument", "<name>", "auto", "FASTA library: PeptDeep instrument (auto: timsTOF for .d runs, SciexTOF for "
                          ".wiff, QE otherwise)", false);
    setValidStrings_("instrument", {"auto", "QE", "Lumos", "timsTOF", "SciexTOF", "ThermoTOF"});
    registerDoubleOption_("nce", "<nce>", 30.0, "FASTA library: PeptDeep normalized collision energy", false);
    registerStringOption_("library_im", "<mode>", "auto", "FASTA library: predict precursor 1/K0 from PeptDeep CCS (auto: for .d runs)", false);
    setValidStrings_("library_im", {"auto", "on", "off"});
    registerStringOption_("library_cache", "<dir>", "", "FASTA library: directory of the cached libraries (default: out_dir). A library "
                          "built earlier with the same FASTA, settings, runs, models and OpenMS revision is reused", false);
    registerStringOption_("model_dir", "<dir>", "", "FASTA library: directory with peptdeep_{rt,ms2,ccs}_dynamic.onnx; empty uses the "
                          "installed models", false, true);
    registerFlag_("survey_only", "Advanced: FASTA library: stop after the survey and write <fasta stem>.survey.tsv (design stage)", true);
    registerFlag_("survey_predicted", "Advanced: FASTA library: survey with the top 6 predicted fragments of every pair instead of theoretical "
                  "y ions (MS2 prediction for every candidate; PLAN_SURVEY V1)", true);
    registerFlag_("survey_normalized", "Advanced: FASTA library: survey co-elution on fragments scaled to their own maximum (PLAN_SURVEY V2)", true);
    registerFlag_("survey_features", "Advanced: FASTA library: also write <fasta stem>.survey.features.tsv (best survey candidate per "
                  "precursor and run)", true);
    registerInputFile_("scoring_lib", "<file>", "", "Optional library (PQP or OpenDIA .parquet) with up to 12 scoring fragments per precursor "
                       "(same precursor IDs); default: the search fragments", false);
    setValidFormats_("scoring_lib", {"pqp", "parquet"});
    registerOutputDir_("out_dir", "<dir>", "", "Output directory");
    registerDoubleOption_("q", "<q>", 0.01, "Run and global peak-group q-value cut for quantification matrices", false);
    setMinFloat_("q", 0.0);
    setMaxFloat_("q", 1.0);
    registerDoubleOption_("peptide_q", "<q>", 0.01, "Selected-context and global peptide q-value cut for the peptide matrix", false);
    setMinFloat_("peptide_q", 0.0);
    setMaxFloat_("peptide_q", 1.0);
    registerDoubleOption_("protein_q", "<q>", 0.01, "Selected-context and global protein q-value cut for the protein matrix", false);
    setMinFloat_("protein_q", 0.0);
    setMaxFloat_("protein_q", 1.0);
    registerStringOption_("inference_context", "<context>", "run-specific", "Peptide/protein reporting context, combined with their global q-value cuts; all contexts are exported", false);
    setValidStrings_("inference_context", {"run-specific", "experiment-wide"});
    registerStringOption_("output_format", "<format>", "tsv", "High-level tables and matrix format (canonical QPX is always Parquet)", false);
    setValidStrings_("output_format", {"tsv", "parquet", "both"});
    registerFlag_("legacy_outputs", "Also write the original byte-compatible precursor TSVs in legacy/");
    registerFlag_("peakmaps_out", "Also export measured MS1/MS2 peak maps as <run>.xipm (additional spectrum replay)");
    registerDoubleOption_("peakmaps_q", "<q>", 0.01, "Run q-value cut for optional XIPM peak maps", false);
    setMinFloat_("peakmaps_q", 0.0);
    setMaxFloat_("peakmaps_q", 1.0);
    registerDoubleOption_("peakmaps_rt_padding", "<seconds>", 15.0, "RT padding on each side of the winning peak for optional XIPM export", false);
    setMinFloat_("peakmaps_rt_padding", 0.0);
    registerFlag_("peakmaps_decoys", "Include decoys passing peakmaps_q in optional XIPM export");
    registerFlag_("features_out", "Also write <run>.features.tsv with every scored window and its features");
    registerIntOption_("seed", "<n>", 0, "Scorer seed", false);
    setMinInt_("seed", 0);
    registerFlag_("im_array_as_vssc", "Interpret a spectrum ion-mobility array whose unit is not CCS as 1/K0");
    registerIntOption_("max_peaks", "<n>", 5, "Advanced: candidate peaks kept per precursor by the search, best fragment co-elution first", false, true);
    setMinInt_("max_peaks", 1);
    registerDoubleOption_("valley_fraction", "<fraction>", 0.5, "Advanced: in the search, two maxima stay separate candidate peaks only if the valley "
                          "between them is at most this fraction of the lower one (1 keeps every smoothed maximum)", false, true);
    setMinFloat_("valley_fraction", 0.0);
    setMaxFloat_("valley_fraction", 1.0);
    registerDoubleOption_("im_tolerance_scale", "<factor>", 0.6, "Advanced: multiply the calibrated 1/K0 half width used for search and extraction "
                          "(PLAN_IMNARROW: 0.6; 1 = as calibrated)", false, true);
    setMinFloat_("im_tolerance_scale", 0.05);
    setMaxFloat_("im_tolerance_scale", 4.0);
    registerDoubleOption_("trace_tolerance_ppm", "<ppm>", 10.0, "Advanced: fragment m/z half width of the 6-fragment trace tensor in ppm "
                          "(PLAN_TRPPM: 10, +18 % IDs at 1 % true FDP on the r ~ 7 HeLa entrapment)", false, true);
    setMinFloat_("trace_tolerance_ppm", 0.0);
    registerIntOption_("gbm_rounds", "<n>", 3, "Advanced, experimental (PRECHECK_GBMFIX): hard-negative GBM rescore rounds", false, true);
    setMinInt_("gbm_rounds", 1);
    registerDoubleOption_("gbm_positive_q", "<q>", 0.01, "Advanced, experimental (PRECHECK_GBMFIX): GBM positives are target winners at this pooled q", false, true);
    setMinFloat_("gbm_positive_q", 0.0001);
    setMaxFloat_("gbm_positive_q", 1.0);
    registerIntOption_("trace_tensor_gb", "<GiB>", 64, "Advanced: hard cap on the 6-fragment trace tensor in GiB (about 6.1 KiB per candidate window with ion mobility); "
                       "raise it for very large libraries, memory permitting", false, true);
    setMinInt_("trace_tensor_gb", 1);
    registerFlag_("no_rt_refine", "Advanced: skip the per-run RT refinement (PLAN_RTR): by default a ridge RT correction on sequence composition is fitted from the run's own q<=0.01 IDs after a first pass, then search, extraction and scoring run again with corrected RTs and a narrower RT window (about 1.8-1.9x wall time)", true);
    registerIntOption_("rt_refine_min_ids", "<n>", 500, "Advanced: skip the RT refinement with fewer training IDs", false, true);
    setMinInt_("rt_refine_min_ids", 1);
    registerFlag_("no_spectrum_cache", "Advanced: read the run file on every pass instead of caching its MS1/MS2 spectra once in a temporary "
                  "directory (outputs are identical either way; the cache needs about 12 bytes per peak of free disk space)", true);
    registerStringOption_("spectrum_cache_dir", "<dir>", "", "Advanced: parent directory of the temporary spectrum cache (default: out_dir); "
                          "a fast local disk is best. The cache is removed after each run", false, true);
    registerIntOption_("extraction_groups", "<n>", 0, "Advanced: with the spectrum cache, extract features in this many precursor m/z groups, "
                       "each replaying only the spectra that can isolate its precursors (outputs are identical; fewer windows are held at "
                       "once). 0 = auto (16)", false, true);
    setMinInt_("extraction_groups", 0);
    registerIntOption_("extraction_parallel", "<n>", 0, "Advanced: extraction groups run concurrently, sharing -threads (0 = auto, threads / 4)", false, true);
    setMinInt_("extraction_parallel", 0);
    registerStringOption_("decoys", "<mode>", "library",
                          "library: use the library's decoys. predict: replace them by pseudo-reversed decoys predicted like the "
                          "targets (own RT, ion mobility and fragment intensities, PeptDeep; DESIGN_OPENDIA_DECOYS), written to "
                          "<out_dir>/<lib>.predicted_decoys.pqp. Needs a WITH_ONNX build; for predicted libraries", false);
    setValidStrings_("decoys", {"library", "predict"});
    registerStringOption_("decoy_model_dir", "<dir>", "",
                          "decoys predict: directory with peptdeep_{rt,ms2,ccs}_dynamic.onnx; empty uses the installed models", false, true);
    registerStringOption_("decoy_cache", "<dir>", "",
                          "decoys predict: directory of the predicted library (default: out_dir). A library predicted earlier from the "
                          "same input library, models and OpenMS revision is reused instead of predicted again, so a study pays the "
                          "prediction once even when every run is a separate call", false, true);
  }

  static std::string stem_(const std::string& path)
  {
    return File::stemName(path);
  }

  static void number_(std::string& row, double v)
  {
    char buffer[32];
    const auto r = std::to_chars(buffer, buffer + sizeof(buffer), v);
    row.append(1, '\t').append(buffer, r.ptr);
  }

  /// -decoys predict: drop the library's decoys and give every target a pseudo-reversed decoy predicted like a target
  /// (library/preparation, settings estimated from the library). Returns the path of the target/decoy PQP it writes to
  /// @p cache_dir. The library is predicted once per
  /// input: <stem>.predicted_decoys.pqp.key holds the SHA-1 of the input library, the model directory and the OpenMS
  /// revision, and a library with a matching key is reused. Prediction writes to a temporary name and renames, so a
  /// concurrent or interrupted call never leaves a partial library under the final name.
  /// TransitionPQPFile remains the format boundary; sequence generation and inference use OpenDIA-owned records.
  std::string predictDecoys_(const std::string& lib, const std::string& cache_dir)
  {
    const std::string stem = stem_(lib);
    const std::string out = cache_dir + "/" + stem + ".predicted_decoys.pqp";
    const std::string key_file = out + ".key";
    const std::string key = "lib_sha1=" + FileHandler::computeFileHash(lib) + "\nmodel_dir=" + getStringOption_("decoy_model_dir") +
                            "\nopenms=" + VersionInfo::getVersion() + " " + VersionInfo::getRevision() + "\n";
    {
      std::ifstream in(key_file);
      const std::string cached((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      if (cached == key && std::filesystem::exists(out))
      {
        writeLogInfo_("Decoy prediction: reusing " + out + " (same input library, models and OpenMS revision)");
        return out;
      }
    }
    std::filesystem::create_directories(cache_dir);
    OpenSwath::LightTargetedExperiment all;
    // The written library has canonical precursor IDs (the output's precursor_id refers to it); sequence, charge and
    // proteins are those of the input.
    TransitionPQPFile().convertPQPToTargetedExperiment(lib.c_str(), all);
    OpenSwath::LightTargetedExperiment targets;
    std::set<std::string> decoy_ids;
    for (const auto& t : all.transitions)
    {
      if (t.getDecoy()) decoy_ids.insert(t.peptide_ref);
    }
    for (auto& c : all.compounds)
    {
      if (!decoy_ids.contains(c.id)) targets.compounds.push_back(std::move(c));
    }
    for (auto& t : all.transitions)
    {
      if (!t.getDecoy()) targets.transitions.push_back(std::move(t));
    }
    std::set<std::string> used;
    for (const auto& c : targets.compounds) used.insert(c.protein_refs.begin(), c.protein_refs.end());
    for (auto& p : all.proteins)
    {
      if (used.contains(p.id)) targets.proteins.push_back(std::move(p));
    }
    all = OpenSwath::LightTargetedExperiment();
    writeLogInfo_("Predicting decoys for " + std::to_string(targets.compounds.size()) + " target precursors (" +
                  std::to_string(decoy_ids.size()) + " library decoys dropped)");

    const std::string tag = "." + std::to_string(::getpid()) + ".tmp.pqp";
    const std::string targets_pqp = cache_dir + "/" + stem + ".targets" + tag;
    const std::string partial = cache_dir + "/" + stem + ".predicted_decoys" + tag;
    std::filesystem::remove(targets_pqp);
    std::filesystem::remove(partial);
    std::filesystem::remove(key_file);
    // Remove the temporaries also when prediction throws (e.g. too few annotated targets).
    struct Cleanup
    {
      std::vector<std::string> paths;
      ~Cleanup()
      {
        std::error_code ignored;
        for (const auto& p : paths) std::filesystem::remove(p, ignored);
      }
    } cleanup{{targets_pqp, partial}};
    TransitionPQPFile().convertLightTargetedExperimentToPQP(targets_pqp.c_str(), targets);
    targets = OpenSwath::LightTargetedExperiment();

    // Preserve the canonical-ID and numeric round trip of the previous preparation path.
    TransitionPQPFile().convertPQPToTargetedExperiment(targets_pqp.c_str(), targets);
    auto assays = fromPQP(std::move(targets));
    targets = OpenSwath::LightTargetedExperiment();
    Preparation::PredictionOptions options;
    options.model_directory = getStringOption_("decoy_model_dir");
    options.threads = std::max(1, getIntOption_("threads"));
    options.log_type = log_type_;
    auto predicted = toPQP(Preparation::predictDecoys(std::move(assays), options));
    TransitionPQPFile().convertLightTargetedExperimentToPQP(partial.c_str(), predicted);
    std::filesystem::rename(partial, out);
    std::ofstream(key_file) << key;
    return out;
  }

  ExitCodes main_(int, const char**) override
  {
    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<std::string> inputs = getStringList_("in");
    const std::string out_dir = getOutputDirOption("out_dir");
    const double q_cut = getDoubleOption_("q");
    if (inputs.empty())
    {
      writeLogError_("No input runs.");
      return ILLEGAL_PARAMETERS;
    }
    for (const auto value : {q_cut, getDoubleOption_("peptide_q"), getDoubleOption_("protein_q"), getDoubleOption_("peakmaps_q")})
    {
      if (!std::isfinite(value) || value < 0.0 || value > 1.0)
      {
        writeLogError_("Output q-value cuts must be finite and in [0, 1].");
        return ILLEGAL_PARAMETERS;
      }
    }
    if (!std::isfinite(getDoubleOption_("peakmaps_rt_padding")))
    {
      writeLogError_("Peak-map RT padding must be finite.");
      return ILLEGAL_PARAMETERS;
    }
    std::filesystem::create_directories(out_dir);
    std::map<std::string, int> stems;
    for (const auto& in : inputs)
    {
      if (++stems[stem_(in)] > 1)
      {
        writeLogError_("Two input runs share the file name " + stem_(in) + "; outputs would collide.");
        return ILLEGAL_PARAMETERS;
      }
    }

    OpenDIA::RunOptions options;
    options.threads = std::max(1, getIntOption_("threads"));
    options.seed = static_cast<std::uint64_t>(getIntOption_("seed"));
    options.im_array_as_vssc = getFlag_("im_array_as_vssc");
    options.max_peaks = static_cast<std::size_t>(getIntOption_("max_peaks"));
    options.search_valley_fraction = getDoubleOption_("valley_fraction");
    options.im_tolerance_scale = getDoubleOption_("im_tolerance_scale");
    options.trace_tolerance_ppm = getDoubleOption_("trace_tolerance_ppm");
    options.gbm_rounds = getIntOption_("gbm_rounds");
    options.gbm_positive_q = getDoubleOption_("gbm_positive_q");
    if (!(options.search_valley_fraction > 0.0))
    {
      writeLogError_("-valley_fraction must be in (0, 1].");
      return ILLEGAL_PARAMETERS;
    }
    if (!(options.trace_tolerance_ppm > 0.0))
    {
      writeLogError_("-trace_tolerance_ppm must be positive.");
      return ILLEGAL_PARAMETERS;
    }
    options.trace_tensor_gb = static_cast<std::size_t>(getIntOption_("trace_tensor_gb"));
    options.rt_refine = !getFlag_("no_rt_refine");
    options.rt_refine_min_ids = static_cast<std::size_t>(getIntOption_("rt_refine_min_ids"));
    options.calibration.panel_targets = 40000;
    if (!getFlag_("no_spectrum_cache"))
    {
      const std::string cache_dir = getStringOption_("spectrum_cache_dir");
      options.spectrum_cache_dir = cache_dir.empty() ? out_dir : cache_dir;
    }
    options.extraction_groups = static_cast<std::size_t>(getIntOption_("extraction_groups"));
    options.extraction_parallel = static_cast<std::size_t>(getIntOption_("extraction_parallel"));

    std::string lib_path = getStringOption_("lib");
    std::string scoring_path = getStringOption_("scoring_lib");
    const std::string fasta = getStringOption_("fasta");
    if (lib_path.empty() == fasta.empty())
    {
      writeLogError_("Give either -lib or -fasta.");
      return ILLEGAL_PARAMETERS;
    }
    if (!fasta.empty())
    {
      if (!scoring_path.empty() || getStringOption_("decoys") != "library")
      {
        writeLogError_("-fasta builds its own decoys and scoring fragments: -scoring_lib and -decoys predict do not apply.");
        return ILLEGAL_PARAMETERS;
      }
      OpenDIA::FastaLibraryOptions fo;
      fo.fasta = fasta;
      const std::string cache = getStringOption_("library_cache");
      fo.cache_dir = cache.empty() ? out_dir : cache;
      auto& c = fo.candidates;
      c.enzyme = getStringOption_("enzyme");
      c.missed_cleavages = getIntOption_("missed_cleavages");
      c.min_length = getIntOption_("min_length");
      c.max_length = getIntOption_("max_length");
      c.min_charge = getIntOption_("min_charge");
      c.max_charge = getIntOption_("max_charge");
      c.min_mz = getDoubleOption_("min_mz");
      c.max_mz = getDoubleOption_("max_mz");
      const StringList fixed = getStringList_("fixed_mods"), variable = getStringList_("variable_mods");
      c.fixed_mods.assign(fixed.begin(), fixed.end());
      c.variable_mods.assign(variable.begin(), variable.end());
      c.max_variable_mods = getIntOption_("max_variable_mods");
      // Instrument and mobility from the run type: all runs must agree.
      const auto kind = [](const std::string& run)
      {
        std::string e = std::filesystem::path(run).extension().string();
        std::transform(e.begin(), e.end(), e.begin(), [](unsigned char x) { return static_cast<char>(std::tolower(x)); });
        return e == ".d" ? std::string("timsTOF") : e == ".wiff" || e == ".wiff2" ? std::string("SciexTOF") : std::string("QE");
      };
      std::set<std::string> kinds;
      for (const auto& in : inputs) kinds.insert(kind(in));
      auto& p = fo.prediction;
      p.model_directory = getStringOption_("model_dir");
      p.instrument = getStringOption_("instrument");
      if (p.instrument == "auto")
      {
        if (kinds.size() != 1)
        {
          writeLogError_("-instrument auto needs runs of one kind (.d, .wiff or other); set -instrument.");
          return ILLEGAL_PARAMETERS;
        }
        p.instrument = *kinds.begin();
      }
      p.nce = static_cast<float>(getDoubleOption_("nce"));
      const std::string im = getStringOption_("library_im");
      p.predict_im = im == "on" || (im == "auto" && kinds.size() == 1 && *kinds.begin() == "timsTOF");
      p.threads = options.threads;
      fo.fragments.count = static_cast<std::size_t>(getIntOption_("fragments"));
      fo.fragments.min_mz = getDoubleOption_("fragment_min_mz");
      fo.fragments.max_mz = getDoubleOption_("fragment_max_mz");
      fo.prefilter_pairs = static_cast<std::size_t>(getIntOption_("prefilter_pairs"));
      fo.panel_pairs = options.calibration.panel_targets;
      fo.survey_only = getFlag_("survey_only");
      fo.fragments.predicted_survey = getFlag_("survey_predicted");
      fo.survey_features = getFlag_("survey_features");
      options.survey_normalized_coelution = getFlag_("survey_normalized");
      writeLogInfo_("FASTA library: " + fasta + " instrument=" + p.instrument + " nce=" + std::to_string(p.nce) + " im=" +
                    (p.predict_im ? "on" : "off") + " prefilter_pairs=" + std::to_string(fo.prefilter_pairs));
      const auto paths = OpenDIA::buildFastaLibrary(fo, inputs, options, [this](const std::string& line) { writeLogInfo_(line); });
      if (fo.survey_only)
      {
        writeLogInfo_("Survey only: " + paths.survey);
        return EXECUTION_OK;
      }
      lib_path = paths.search;
      scoring_path = paths.library;
    }
    if (getStringOption_("decoys") == "predict")
    {
      if (!scoring_path.empty())
      {
        writeLogError_("-decoys predict does not support -scoring_lib (its decoys would not match the predicted ones).");
        return ILLEGAL_PARAMETERS;
      }
      if (OpenDIA::isParquetLibrary(lib_path))
      {
        writeLogError_("-decoys predict needs a PQP -lib (the prediction format boundary uses the PQP reader/writer).");
        return ILLEGAL_PARAMETERS;
      }
      const std::string cache = getStringOption_("decoy_cache");
      lib_path = predictDecoys_(lib_path, cache.empty() ? out_dir : cache);
    }
    writeLogInfo_("Reading library " + lib_path);
    OpenDIA::AnalysisLibrary library =
      OpenDIA::AnalysisLibrary::read(lib_path, scoring_path, [this](const std::string& line) { writeLogInfo_(line); });
    if (!fasta.empty())
    {
      // As in the survey: only the panel pairs (the leading precursors of the FASTA library) calibrate.
      const std::size_t panel = OpenDIA::Preparation::surveyPanelPrecursors(lib_path);
      std::erase_if(library.calibration_pairs, [&](const auto& pair) { return pair.first >= panel; });
      writeLogInfo_("FASTA library: calibration pairs=" + std::to_string(library.calibration_pairs.size()));
    }
    std::size_t library_rss = 0;
    SysInfo::getProcessMemoryConsumption(library_rss);
    writeLogInfo_("Library: precursors=" + std::to_string(library.precursors.size()) + " pairs=" + std::to_string(library.pairs) +
                  " rss_mb=" + std::to_string(library_rss / 1024));

    const auto genes = OpenDIA::readGeneAnnotations(lib_path);
    const OpenDIA::ReportingOptions reporting{q_cut, getDoubleOption_("peptide_q"), getDoubleOption_("protein_q"), getStringOption_("inference_context")};
    const bool legacy_outputs = getFlag_("legacy_outputs");
    const std::string legacy_dir = out_dir + "/legacy";
    const std::string qpx_dir = out_dir + "/qpx";
    std::filesystem::create_directories(qpx_dir);
    ParamXMLFile().store(out_dir + "/parameters.ini", getParam_());
    if (legacy_outputs) std::filesystem::create_directories(legacy_dir);

    // Best (calibrated) score of each precursor per run, for the global q-value.
    struct Hit
    {
      std::size_t run;
      std::uint32_t window;
      float score;
      double q;
    };
    std::vector<std::vector<Hit>> hits(library.precursors.size());
    std::vector<OpenDIA::RunResult> results;
    for (std::size_t r = 0; r < inputs.size(); ++r)
    {
      const std::string name = stem_(inputs[r]);
      writeLogInfo_("Run " + std::to_string(r + 1) + "/" + std::to_string(inputs.size()) + ": " + inputs[r]);
      options.peakmaps_path = getFlag_("peakmaps_out") ? out_dir + "/" + name + ".xipm" : "";
      options.peakmaps_q = getDoubleOption_("peakmaps_q");
      options.peakmaps_rt_padding_s = getDoubleOption_("peakmaps_rt_padding");
      options.peakmaps_include_decoys = getFlag_("peakmaps_decoys");
      options.peakmaps_run_id = r;
      OpenDIA::RunAnalysis analysis(library, options, [&](const std::string& line) { writeLogInfo_("  " + line); });
      OpenDIA::RunResult result = analysis.run(inputs[r]);
      {
        std::ofstream json(out_dir + "/" + name + ".calibration.json");
        json << result.calibration_json;
      }
      if (legacy_outputs && !writeRun_(legacy_dir + "/" + name + ".precursors.tsv", name, library, result))
      {
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      if (getFlag_("features_out") && !writeFeatures_(out_dir + "/" + name + ".features.tsv", library, result))
      {
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      for (std::uint32_t i = 0; i < result.windows.size(); ++i)
      {
        const auto& w = result.windows[i];
        if (w.score > OpenDIA::CrossFitScorer::NOT_KEPT) hits[w.precursor].push_back({r, i, w.score, w.q_value});
      }
      const auto evidence_start = std::chrono::steady_clock::now();
      OpenDIA::writeCandidateEvidence(qpx_dir, library, result);
      writeLogInfo_("Candidate evidence: windows=" + std::to_string(result.windows.size()) + " seconds=" +
                    std::to_string(std::chrono::duration<double>(std::chrono::steady_clock::now() - evidence_start).count()));
      // Candidate diagnostics have been persisted. Retain windows for cross-run reporting.
      std::vector<float>().swap(result.features);
      results.push_back(std::move(result));
    }

    const auto reporting_start = std::chrono::steady_clock::now();
    const auto quantification = OpenDIA::quantify(library, results, reporting);
    for (const auto& warning : quantification.inference_warnings)
      writeLogInfo_("WARNING: entity inference unavailable: " + warning);
    const auto& global_q = quantification.global_precursor_q;
    OpenDIA::writeQpxDataset(qpx_dir, library, results, quantification, genes, reporting);
    OpenDIA::writeReports(out_dir, library, results, quantification, genes, reporting, getStringOption_("output_format"));
    writeLogInfo_("QPX and reports: seconds=" +
                  std::to_string(std::chrono::duration<double>(std::chrono::steady_clock::now() - reporting_start).count()));

    std::ofstream all;
    if (legacy_outputs)
    {
      all.open(legacy_dir + "/precursors.tsv");
      all << "run\tprecursor_id\tsequence\tcharge\tproteins\trt_apex\trt_start\trt_end\tarea\tscore\tq_value\tglobal_q_value\n";
    }
    std::size_t rows = 0, global_targets = 0;
    for (std::uint32_t p = 0; p < library.precursors.size(); ++p)
    {
      const auto& pre = library.precursors[p];
      if (pre.decoy) continue;
      if (global_q[p] <= q_cut) ++global_targets;
      for (const auto& h : hits[p])
      {
        if (h.q > q_cut || global_q[p] > q_cut) continue;
        const auto& w = results[h.run].windows[h.window];
        const auto& proteins = library.protein_sets[pre.proteins];
        std::string row = stem_(inputs[h.run]) + '\t' + std::string(pre.id) + '\t' + pre.sequence + '\t' + std::to_string(pre.charge) + '\t' + proteins;
        for (double v : {w.peak.rt_apex, w.peak.rt_start, w.peak.rt_end, w.peak.area, static_cast<double>(w.score), w.q_value, global_q[p]}) number_(row, v);
        if (legacy_outputs) all << row << '\n';
        ++rows;
      }
    }
    if (legacy_outputs)
    {
      all.close();
      if (!all) return CANNOT_WRITE_OUTPUT_FILE;
    }
    writeLogInfo_("Combined: " + std::to_string(global_targets) + " target precursors at global q <= " + std::to_string(q_cut) + ", " +
                  std::to_string(rows) + " run-level rows; seconds=" +
                  std::to_string(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()));
    return EXECUTION_OK;
  }

  /// One row per precursor winner (the kept window), with its run q-value.
  static bool writeRun_(const std::string& path, const std::string& run, const OpenDIA::AnalysisLibrary& library, const OpenDIA::RunResult& result)
  {
    std::ofstream out(path);
    out << "run\tcandidate_row\tprecursor_id\tis_decoy\tsequence\tcharge\trt_apex\trt_start\trt_end\tarea\ttransition_support\tscore\tq_value\n";
    for (const auto& w : result.windows)
    {
      if (!(w.score > OpenDIA::CrossFitScorer::NOT_KEPT)) continue;
      const auto& p = library.precursors[w.precursor];
      std::string row = run + '\t' + std::to_string(w.candidate) + '\t' + std::string(p.id) + '\t' + (p.decoy ? "1" : "0") + '\t' + p.sequence + '\t' + std::to_string(p.charge);
      for (double v : {w.peak.rt_apex, w.peak.rt_start, w.peak.rt_end, w.peak.area}) number_(row, v);
      row += '\t' + std::to_string(w.peak.transition_support);
      for (double v : {static_cast<double>(w.score), w.q_value}) number_(row, v);
      out << row << '\n';
    }
    out.close();
    return static_cast<bool>(out);
  }

  /// Every scored window with its features (candidate_row = search feature table row).
  static bool writeFeatures_(const std::string& path, const OpenDIA::AnalysisLibrary& library, const OpenDIA::RunResult& result)
  {
    std::ofstream out(path);
    out << "candidate_row\tprecursor_id\tis_decoy\tfold\trt_apex";
    for (const auto& n : result.feature_names) out << '\t' << n;
    out << "\tscore\tq_value\tnetwork_score\tnetwork_q_value\n";
    const std::size_t d = result.feature_names.size();
    for (std::size_t i = 0; i < result.windows.size(); ++i)
    {
      const auto& w = result.windows[i];
      const auto& p = library.precursors[w.precursor];
      std::string row = std::to_string(w.candidate) + '\t' + std::string(p.id) + '\t' + (p.decoy ? "1" : "0") + '\t' + std::to_string(p.fold);
      number_(row, w.peak.rt_apex);
      for (std::size_t k = 0; k < d; ++k) number_(row, result.features[i * d + k]);
      number_(row, w.score);
      number_(row, w.q_value);
      number_(row, w.network_score);
      number_(row, w.network_q_value);
      out << row << '\n';
    }
    out.close();
    return static_cast<bool>(out);
  }
};

int main(int argc, const char** argv)
{
  TOPPOpenDIA tool;
  return tool.main(argc, argv);
}
