// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
#include <workflow/FastaWorkflow.h>

#include <scoring/Scorer.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>

namespace OpenMS::OpenDIA
{
  namespace
  {
    using Log = std::function<void(const std::string&)>;
    namespace fs = std::filesystem;

    std::string readFile(const std::string& path)
    {
      std::ifstream in(path, std::ios::binary);
      return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
    }

    /// The file is reused if it exists and its .key holds @p key.
    bool cached(const std::string& path, const std::string& key)
    {
      std::error_code error;
      return fs::is_regular_file(path, error) && readFile(path + ".key") == key;
    }

    /// Write @p key next to @p path (temporary name, then rename).
    void writeKey(const std::string& path, const std::string& key)
    {
      const std::string partial = path + ".key." + std::to_string(::getpid()) + ".tmp";
      {
        std::ofstream out(partial, std::ios::binary);
        out << key;
        if (!out) throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, partial);
      }
      fs::rename(partial, path + ".key");
    }

    std::string join(const std::vector<std::string>& v)
    {
      std::string out;
      for (const auto& s : v) out += (out.empty() ? "" : ",") + s;
      return out;
    }

    std::string number(double v)
    {
      std::ostringstream s;
      s.precision(17);
      s << v;
      return s.str();
    }

    std::string modelKey(const Preparation::PredictionSettings& settings, const std::string& name)
    {
      return name + "=" + FileHandler::computeFileHash(Preparation::modelPath(settings, name)) + "\n";
    }

    /// The pair statistic and FNV-1a of the target TRAML_ID per pair, from a survey TSV.
    void readSurvey(const std::string& path, std::vector<double>& statistic, std::vector<std::uint64_t>& tie)
    {
      statistic.clear();
      tie.clear();
      std::ifstream in(path);
      std::string line;
      std::getline(in, line);
      while (std::getline(in, line))
      {
        std::istringstream row(line);
        std::string field[5];
        for (auto& x : field) std::getline(row, x, '\t');
        statistic.push_back(std::stod(field[4]));
        tie.push_back(CrossFitScorer::fnv1a(field[1]));
      }
    }
  }

  FastaLibraryPaths buildFastaLibrary(const FastaLibraryOptions& options, const std::vector<std::string>& runs, const RunOptions& run_options,
                                      const Log& log)
  {
    const auto say = [&](const std::string& line) { if (log) log(line); };
    const auto& c = options.candidates;
    const auto& p = options.prediction;
    const auto& f = options.fragments;
    fs::create_directories(options.cache_dir);
    const std::string stem = (fs::path(options.cache_dir) / fs::path(options.fasta).stem()).string();
    FastaLibraryPaths paths;
    paths.candidates = stem + ".candidates.parquet";

    // 1. Candidates with coordinates and survey fragments.
    std::string key = "opendia_fasta_library=1\nfasta_sha1=" + FileHandler::computeFileHash(options.fasta) + "\nenzyme=" + c.enzyme +
                      "\nmissed_cleavages=" + std::to_string(c.missed_cleavages) + "\nlength=" + std::to_string(c.min_length) + "-" +
                      std::to_string(c.max_length) + "\ncharge=" + std::to_string(c.min_charge) + "-" + std::to_string(c.max_charge) +
                      "\nmz=" + number(c.min_mz) + "-" + number(c.max_mz) + "\nfixed=" + join(c.fixed_mods) + "\nvariable=" +
                      join(c.variable_mods) + "\nmax_variable=" + std::to_string(c.max_variable_mods) + "\ninstrument=" + p.instrument +
                      "\nnce=" + number(p.nce) + "\nim=" + (p.predict_im ? "1" : "0") + "\nfragment_mz=" + number(f.min_mz) + "-" +
                      number(f.max_mz) + "\npanel_pairs=" + std::to_string(options.panel_pairs) + "\n" +
                      (f.predicted_survey ? "predicted_survey=1\n" : "") +
                      modelKey(p, "peptdeep_rt_dynamic.onnx") + modelKey(p, "peptdeep_ms2_dynamic.onnx") +
                      (p.predict_im ? modelKey(p, "peptdeep_ccs_dynamic.onnx") : std::string()) + "openms=" + VersionInfo::getVersion() +
                      " " + VersionInfo::getRevision() + "\n";
    const std::string candidates_key = key;
    if (cached(paths.candidates, key))
    {
      say("FASTA library: reusing " + paths.candidates);
    }
    else
    {
      auto set = Preparation::enumerateCandidates(options.fasta, c, log);
      if (set.pairs() == 0) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The FASTA gives no target-decoy pairs", options.fasta);
      Preparation::predictCoordinates(set, p, log);
      const auto panel = Preparation::calibrationPanel(set, options.panel_pairs);
      Preparation::writeSurveyLibrary(set, panel, p, f, paths.candidates, {{"opendia.fasta.source", fs::absolute(options.fasta).string()}}, log);
      writeKey(paths.candidates, key);
    }

    // 2. Survey of every run (K > 0).
    std::vector<double> statistic;
    std::vector<std::uint64_t> tie;
    std::string survey_key;
    if (options.prefilter_pairs > 0 || options.survey_only)
    {
      paths.survey = stem + ".survey.tsv";
      survey_key = candidates_key + "survey=1\nmax_peaks=" + std::to_string(run_options.max_peaks) + "\nvalley_fraction=" +
                   number(run_options.search_valley_fraction) + "\nim_tolerance_scale=" + number(run_options.im_tolerance_scale) +
                   "\nmin_transitions=" + std::to_string(run_options.min_transitions) + "\nim_array_as_vssc=" +
                   (run_options.im_array_as_vssc ? "1" : "0") + "\npanel_targets=" + std::to_string(run_options.calibration.panel_targets) +
                   "\nwindow_panel_targets=" + std::to_string(run_options.window_panel_targets) + "\n" +
                   (run_options.survey_normalized_coelution ? "normalized_coelution=1\n" : "");
      for (const auto& run : runs)
      {
        std::error_code error;
        const auto canonical = fs::weakly_canonical(run, error).string();
        std::uintmax_t size = 0;
        if (fs::is_directory(run, error))
        {
          for (const auto& e : fs::recursive_directory_iterator(run, error))
            if (e.is_regular_file(error)) size += e.file_size(error);
        }
        else size = fs::file_size(run, error);
        survey_key += "run=" + canonical + " size=" + std::to_string(size) + "\n";
      }
      if (cached(paths.survey, survey_key))
      {
        say("FASTA survey: reusing " + paths.survey);
      }
      else
      {
        AnalysisLibrary lib = AnalysisLibrary::read(paths.candidates, "", log);
        // Only the panel pairs (the leading precursors) calibrate: theoretical fragments carry no intensities.
        const std::size_t panel = Preparation::surveyPanelPrecursors(paths.candidates);
        std::erase_if(lib.calibration_pairs, [&](const auto& pair) { return pair.first >= panel; });
        say("FASTA survey: library precursors=" + std::to_string(lib.precursors.size()) + " calibration pairs=" +
            std::to_string(lib.calibration_pairs.size()));
        // Pair k is precursors 2k and 2k + 1; pairs without survey fragments are missing from the library (S = 0).
        std::vector<std::size_t> pair_of(lib.precursors.size());
        std::size_t pairs = 0;
        for (std::size_t i = 0; i < lib.precursors.size(); ++i)
        {
          pair_of[i] = std::stoull(std::string(lib.precursors[i].id)) / 2;
          pairs = std::max(pairs, pair_of[i] + 1);
        }
        std::vector<double> best_target(pairs, 0.0), best_decoy(pairs, 0.0);
        std::ofstream features;
        if (options.survey_features)
        {
          features.open(stem + ".survey.features.tsv");
          features.precision(9);
          features << "run\tprecursor_id\tdecoy\tcharge\tcoelution\tsupport\tmatched_scans\tapex_rt_error_s\tfwhm_s\tapex_intensity\tarea\tpeaks_in_window\n";
        }
        for (std::size_t r = 0; r < runs.size(); ++r)
        {
          say("FASTA survey: run " + std::to_string(r + 1) + "/" + std::to_string(runs.size()) + ": " + runs[r]);
          RunAnalysis analysis(lib, run_options, [&](const std::string& line) { say("  " + line); });
          const auto s = analysis.survey(runs[r]);
          for (std::size_t i = 0; i < lib.precursors.size(); ++i)
          {
            auto& best = lib.precursors[i].decoy ? best_decoy[pair_of[i]] : best_target[pair_of[i]];
            best = std::max(best, s.coelution[i]);
            const auto& c = s.best[i];
            if (features.is_open() && c.precursor != UINT32_MAX)
              features << r << '\t' << lib.precursors[i].id << '\t' << lib.precursors[i].decoy << '\t' << lib.precursors[i].charge << '\t'
                       << c.coelution << '\t' << c.transition_support << '\t' << c.matched_scans << '\t' << c.apex_rt_error_seconds << '\t'
                       << c.fwhm_seconds << '\t' << c.apex_intensity << '\t' << c.area << '\t' << c.peaks_in_window << '\n';
          }
        }
        std::vector<std::string> traml(pairs);
        for (std::size_t i = 0; i < lib.precursors.size(); ++i)
          if (!lib.precursors[i].decoy) traml[pair_of[i]] = lib.precursors[i].sequence + "/" + std::to_string(lib.precursors[i].charge);
        const std::string partial = paths.survey + "." + std::to_string(::getpid()) + ".tmp";
        {
          std::ofstream out(partial);
          out << "pair\ttraml_id\ttarget_coelution\tdecoy_coelution\tstatistic\n";
          out.precision(9);
          for (std::size_t k = 0; k < pairs; ++k)
          {
            out << k << '\t' << traml[k] << '\t' << best_target[k] << '\t' << best_decoy[k] << '\t' << std::max(best_target[k], best_decoy[k]) << '\n';
          }
          if (!out) throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, partial);
        }
        fs::rename(partial, paths.survey);
        writeKey(paths.survey, survey_key);
      }
      // Read back from the TSV (9 significant digits) also when just written, so a reused survey selects the same pairs.
      readSurvey(paths.survey, statistic, tie);
      if (options.survey_only) return paths;
    }

    // 3. Final library for the selected pairs.
    paths.library = stem + ".library.parquet";
    paths.search = stem + ".library.search.parquet";
    const std::string final_key = (survey_key.empty() ? candidates_key : survey_key) + "final=2\nprefilter_pairs=" +
                                  std::to_string(options.prefilter_pairs) + "\nfragments=" + std::to_string(f.count) + "\nmin_fragments=" +
                                  std::to_string(f.min_count) + "\n" + (options.survey_pool ? "survey_pool=1\n" : "");
    if (cached(paths.library, final_key) && cached(paths.search, final_key))
    {
      say("FASTA library: reusing " + paths.library);
      return paths;
    }
    std::vector<bool> keep;
    if (options.prefilter_pairs > 0)
    {
      if (options.survey_pool)
      {
        // PLAN_SURVEY +pool: the target's modified sequence and proteins per pair, from the survey library.
        struct Targets final : LibraryVisitor
        {
          std::vector<std::string> sequence;
          std::vector<std::vector<std::string>> proteins;
          std::string_view last;
          void reset() override { sequence.clear(); proteins.clear(); last = {}; }
          void fragment(const FragmentView& f, const PrecursorView& c) override
          {
            if (c.id == last) return;
            last = c.id;
            if (f.decoy) return;
            const std::size_t k = std::stoull(std::string(c.id)) / 2;
            if (k >= sequence.size()) { sequence.resize(k + 1); proteins.resize(k + 1); }
            sequence[k] = c.sequence;
            proteins[k].assign(c.proteins->begin(), c.proteins->end());
          }
        } targets;
        visitParquetLibrary(paths.candidates, targets);
        targets.sequence.resize(statistic.size());
        targets.proteins.resize(statistic.size());
        statistic = Preparation::pooledStatistic(statistic, targets.sequence, targets.proteins);
        say("FASTA prefilter: pooled statistic (charge siblings + 0.5 x protein top-3)");
      }
      const auto selected = Preparation::selectPairs(statistic, tie, options.prefilter_pairs);
      keep.assign(statistic.size(), false);
      for (const auto k : selected) keep[k] = true;
      std::size_t positive = static_cast<std::size_t>(std::count_if(statistic.begin(), statistic.end(), [](double s) { return s > 0.0; }));
      say("FASTA prefilter: kept " + std::to_string(selected.size()) + " of " + std::to_string(statistic.size()) + " pairs (" +
          std::to_string(positive) + " with a survey candidate)");
    }
    const auto stats = Preparation::writeFinalLibrary(
      paths.candidates, [&](std::size_t k) { return keep.empty() || (k < keep.size() && keep[k]); }, p, f, paths.library, paths.search,
      {{"opendia.fasta.source", fs::absolute(options.fasta).string()}, {"opendia.fasta.prefilter_pairs", std::to_string(options.prefilter_pairs)}},
      log);
    if (stats.pairs_out == 0) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The FASTA library has no pairs", paths.library);
    writeKey(paths.library, final_key);
    writeKey(paths.search, final_key);
    return paths;
  }
}
