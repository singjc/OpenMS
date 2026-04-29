// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/OpenSwathBase.h>
#include <OpenMS/ANALYSIS/OPENSWATH/FastaEvidenceFilter.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <iomanip>
#include <memory>
#include <vector>

using namespace OpenMS;

/**
@page TOPP_FastaEvidenceFilter FastaEvidenceFilter

@brief Filters a FASTA database by staged DIA raw-data evidence.

The tool performs a two-stage DIA-only FASTA prefilter:

1. digest the FASTA into a theoretical precursor/fragment search space and keep
   precursors with quick raw-data evidence via TransitionListEvidenceFilter
2. confirm the stage-1 survivors with FragmentIndex-based matched-ion counts

Outputs are a filtered FASTA database and a peptide TSV intended for downstream
transition-library or peptide-property predictors.

@experimental This tool is experimental and intended as a prefilter for large
predicted search spaces, not as a full identification search engine.

@ingroup TOPP
*/

/// @cond TOPPCLASSES
class TOPPFastaEvidenceFilter :
  public TOPPOpenSwathBase
{
public:
  TOPPFastaEvidenceFilter() :
    TOPPOpenSwathBase("FastaEvidenceFilter",
                      "Filter a FASTA database by staged DIA evidence.",
                      true)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input mzML/sqMass files. By default each file is treated as one DIA run.");
    StringList in_formats = {"mzML", "mzXML", "sqMass"};
#ifdef WITH_OPENTIMS
    in_formats.push_back("d");
#endif
    setValidFormats_("in", in_formats);

    registerInputFile_("database", "<file>", "", "Input protein database in FASTA format.");

    registerOutputFile_("out_fasta", "<file>", "", "Filtered output FASTA file.");
    setValidFormats_("out_fasta", {"fasta"});
    registerOutputFile_("out_peptides", "<file>", "", "Filtered peptide precursor table in TSV format.");
    setValidFormats_("out_peptides", {"tsv"});
    registerOutputFile_("out_stage2_scores", "<file>", "",
                        "Optional Stage-2 target/decoy score table in TSV format for score-distribution inspection.",
                        false);
    setValidFormats_("out_stage2_scores", {"tsv"});

    registerStringOption_("aggregation_method", "<any|all>", "any",
                          "How to combine stage-1 evidence across multiple DIA runs.",
                          false, true);
    setValidStrings_("aggregation_method", {"any", "all"});

    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0,
                          "Minimal distance to the upper edge of a SWATH window to still consider a precursor, in Thomson.",
                          false, true);
    setMinFloat_("min_upper_edge_dist", 0.0);
    registerFlag_("pasef", "Data is PASEF data. If omitted, the tool auto-detects ion mobility SWATH windows.");
    registerFlag_("prm", "Data is targeted DIA / PRM-like data with potentially overlapping DIA windows.", true);
    registerFlag_("force", "Override SWATH window gap/overlap sanity checks.", true);
    registerInputFile_("swath_windows_file", "<file>", "",
                       "Optional tab-separated file containing SWATH windows for extraction: lower_offset upper_offset. The first line is a header and will be skipped.",
                       false);
    registerFlag_("sort_swath_maps", "Sort input SWATH files when matching to SWATH windows from swath_windows_file.", true);

    registerFlag_("split_file_input", "Treat all input files as one split SWATH run, with each file containing one SWATH window.", true);
    registerStringOption_("readOptions", "<normal|cache>", "normal",
                          "Whether to run directly on input data or cache data to disk first. If 'cache', set tempDirectory as needed.",
                          false, true);
    setValidStrings_("readOptions", {"normal", "cache"});
    registerStringOption_("tempDirectory", "<tmp>", File::getTempDirectory(), "Temporary directory for cached data.", false, true);
    registerFlag_("keep_cached_files", "If set, do not remove cached files created in tempDirectory.", false);

    registerDoubleOption_("mz_extraction_window_ms1", "<double>", 50.0,
                          "MS1 precursor m/z extraction window full width, in Thomson or ppm.",
                          false, true);
    setMinFloat_("mz_extraction_window_ms1", 0.0);
    registerStringOption_("mz_extraction_window_ms1_unit", "<Th|ppm>", "ppm", "Unit of mz_extraction_window_ms1.", false, true);
    setValidStrings_("mz_extraction_window_ms1_unit", {"Th", "ppm"});
    registerDoubleOption_("im_extraction_window_ms1", "<double>", -1.0,
                          "MS1 ion mobility extraction window full width. -1 disables MS1 IM filtering.",
                          false, true);

    registerDoubleOption_("mz_extraction_window_ms2", "<double>", 50.0,
                          "MS2 fragment m/z extraction window full width, in Thomson or ppm.",
                          false, true);
    setMinFloat_("mz_extraction_window_ms2", 0.0);
    registerStringOption_("mz_extraction_window_ms2_unit", "<Th|ppm>", "ppm", "Unit of mz_extraction_window_ms2.", false, true);
    setValidStrings_("mz_extraction_window_ms2_unit", {"Th", "ppm"});
    registerDoubleOption_("im_extraction_window_ms2", "<double>", -1.0,
                          "MS2 ion mobility extraction window full width. -1 disables MS2 IM filtering.",
                          false, true);

    registerSubsection_("Stage1", "Quick DIA evidence prefilter parameters.");
    registerSubsection_("Stage2", "FragmentIndex-based confirmation parameters.");
    registerSubsection_("Protein", "Protein-level collapsing parameters.");
    registerSubsection_("SearchSpace", "In-silico digestion and precursor/fragment generation parameters.");
    registerSubsection_("Export", "Peptide TSV export parameters.");
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    Param defaults = FastaEvidenceFilter().getParameters();
    if (name == "Stage1") return defaults.copy("Stage1:", true);
    if (name == "Stage2") return defaults.copy("Stage2:", true);
    if (name == "Protein") return defaults.copy("Protein:", true);
    if (name == "SearchSpace") return defaults.copy("SearchSpace:", true);
    if (name == "Export") return defaults.copy("Export:", true);
    return Param();
  }

  ExitCodes main_(int, const char**) override
  {
    const StringList in_files = getStringList_("in");
    if (in_files.empty())
    {
      writeLogError_("Error: No input raw data files provided.");
      return ILLEGAL_PARAMETERS;
    }

    const String database_file = getStringOption_("database");
    const String out_fasta_file = getStringOption_("out_fasta");
    const String out_peptides_file = getStringOption_("out_peptides");
    const String out_stage2_scores_file = getStringOption_("out_stage2_scores");

    if (File::isDirectory(database_file))
    {
      writeLogError_("Error: Parameter '-database' must point to a FASTA file, not a directory: '" + database_file + "'.");
      return ILLEGAL_PARAMETERS;
    }
    if (File::exists(out_fasta_file) && File::isDirectory(out_fasta_file))
    {
      writeLogError_("Error: Parameter '-out_fasta' must point to a FASTA file, not a directory: '" + out_fasta_file + "'.");
      return ILLEGAL_PARAMETERS;
    }
    if (File::exists(out_peptides_file) && File::isDirectory(out_peptides_file))
    {
      writeLogError_("Error: Parameter '-out_peptides' must point to a TSV file, not a directory: '" + out_peptides_file + "'.");
      return ILLEGAL_PARAMETERS;
    }
    if (!out_stage2_scores_file.empty() && File::exists(out_stage2_scores_file) && File::isDirectory(out_stage2_scores_file))
    {
      writeLogError_("Error: Parameter '-out_stage2_scores' must point to a TSV file, not a directory: '" + out_stage2_scores_file + "'.");
      return ILLEGAL_PARAMETERS;
    }

    std::vector<FASTAFile::FASTAEntry> fasta_db;
    FASTAFile().load(database_file, fasta_db);
    if (fasta_db.empty())
    {
      writeLogError_("Error: FASTA database is empty.");
      return INPUT_FILE_EMPTY;
    }

    ChromExtractParams ms1_params = makeChromExtractParams_(
      getDoubleOption_("mz_extraction_window_ms1"),
      getStringOption_("mz_extraction_window_ms1_unit") == "ppm",
      getDoubleOption_("im_extraction_window_ms1"));
    ms1_params.min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");

    ChromExtractParams ms2_params = makeChromExtractParams_(
      getDoubleOption_("mz_extraction_window_ms2"),
      getStringOption_("mz_extraction_window_ms2_unit") == "ppm",
      getDoubleOption_("im_extraction_window_ms2"));
    ms2_params.min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");

    const bool split_file_input = getFlag_("split_file_input");
    const String readoptions = getStringOption_("readOptions");
    const String tmp_dir = File::absolutePath(getStringOption_("tempDirectory")).ensureLastChar('/');
    const bool keep_cached_files = getFlag_("keep_cached_files");
    const bool force = getFlag_("force");
    const bool sort_swath_maps = getFlag_("sort_swath_maps");
    const bool prm = getFlag_("prm");
    const String swath_windows_file = getStringOption_("swath_windows_file");

    std::vector<StringList> run_groups;
    if (split_file_input)
    {
      run_groups.push_back(in_files);
    }
    else
    {
      for (const auto& file : in_files)
      {
        run_groups.push_back(StringList{file});
      }
    }

    std::vector<FastaEvidenceFilter::RunData> runs;
    runs.reserve(run_groups.size());
    for (const auto& run_files : run_groups)
    {
      String per_run_tmp = tmp_dir;
      std::shared_ptr<File::TempDir> per_run_temp_dir;
      if (readoptions == "cache")
      {
        per_run_temp_dir = std::make_shared<File::TempDir>(tmp_dir, keep_cached_files);
        per_run_tmp = per_run_temp_dir->getPath();
      }

      std::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
      std::vector<OpenSwath::SwathMap> swath_maps;
      std::vector<String> swath_map_sources;
      if (!loadSwathFiles(run_files, exp_meta, swath_maps, swath_map_sources, split_file_input,
                          per_run_tmp, readoptions, swath_windows_file,
                          ms2_params.min_upper_edge_dist, force, sort_swath_maps, prm))
      {
        writeLogError_("Error: Failed to load DIA input files.");
        return PARSE_ERROR;
      }

      bool run_pasef = getFlag_("pasef");
      if (!run_pasef)
      {
        run_pasef = std::any_of(swath_maps.begin(), swath_maps.end(),
                                [](const OpenSwath::SwathMap& map)
                                {
                                  return !map.ms1 && map.imLower >= 0.0 && map.imUpper >= 0.0;
                                });
      }

      FastaEvidenceFilter::RunData run;
      run.swath_maps = std::move(swath_maps);
      run.swath_map_sources = std::move(swath_map_sources);
      run.pasef = run_pasef;
      run.cache_dir_guard = per_run_temp_dir;
      runs.push_back(std::move(run));
    }

    FastaEvidenceFilter algorithm;
    Param algorithm_params = algorithm.getParameters();
    algorithm_params.setValue("aggregation_method", getStringOption_("aggregation_method"));
    algorithm_params.update(getParam_().copy("Stage1:"));
    algorithm_params.update(getParam_().copy("Stage2:"));
    algorithm_params.update(getParam_().copy("Protein:"));
    algorithm_params.update(getParam_().copy("SearchSpace:"));
    algorithm_params.update(getParam_().copy("Export:"));
    algorithm_params.setValue("Export:export_stage2_scores", out_stage2_scores_file.empty() ? "false" : "true");
    algorithm.setParameters(algorithm_params);
    algorithm.setLogType(log_type_);

    const auto result = algorithm.filter(runs, fasta_db, ms1_params, ms2_params, static_cast<int>(getIntOption_("threads")));
    FASTAFile().store(out_fasta_file, result.filtered_fasta);
    writePeptideTable_(out_peptides_file, result.confirmed_peptides,
                       getParam_().getValue("Export:export_fragments").toString() == "true");
    if (!out_stage2_scores_file.empty())
    {
      writeStage2ScoreTable_(out_stage2_scores_file, result.stage2_candidate_scores);
    }

    OPENMS_LOG_INFO << result.summary << "\n";
    return EXECUTION_OK;
  }

private:
  static ChromExtractParams makeChromExtractParams_(double mz_window, bool ppm, double im_window)
  {
    ChromExtractParams params;
    params.min_upper_edge_dist = 0.0;
    params.mz_extraction_window = mz_window;
    params.im_extraction_window = im_window;
    params.ppm = ppm;
    params.extraction_function = "tophat";
    params.rt_extraction_window = -1.0;
    params.extra_rt_extract = 0.0;
    return params;
  }

  static void writePeptideTable_(const String& filename,
                                 const std::vector<FastaEvidenceFilter::PeptideEntry>& peptides,
                                 bool export_fragments)
  {
    std::ofstream out(filename.c_str());
    if (!out)
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    out << "peptide_sequence\tmodified_peptide_sequence\tprecursor_mz\tprecursor_charge\tprotein_accession\tgene_name\tproduct_mz\tproduct_charge\tproduct_type\tproduct_ordinal\n";
    out << std::fixed << std::setprecision(6);
    for (const auto& peptide : peptides)
    {
      std::vector<std::string> gene_names;
      gene_names.reserve(peptide.protein_refs.size());
      for (const auto& protein_ref : peptide.protein_refs)
      {
        const auto gene_name_it = peptide.protein_gene_names_by_accession.find(protein_ref);
        gene_names.push_back(gene_name_it != peptide.protein_gene_names_by_accession.end() ? gene_name_it->second : std::string{});
      }
      const String protein_accessions = ListUtils::concatenate(peptide.protein_refs, ";");
      const String joined_gene_names = ListUtils::concatenate(gene_names, ";");

      if (!export_fragments || peptide.fragments.empty())
      {
        out << peptide.peptide_sequence << '\t'
            << peptide.modified_peptide_sequence << '\t'
            << peptide.precursor_mz << '\t'
            << peptide.precursor_charge << '\t'
            << protein_accessions << '\t'
            << joined_gene_names << '\t'
            << '\t' << '\t' << '\t' << '\n';
        continue;
      }

      for (const auto& fragment : peptide.fragments)
      {
        out << peptide.peptide_sequence << '\t'
            << peptide.modified_peptide_sequence << '\t'
            << peptide.precursor_mz << '\t'
            << peptide.precursor_charge << '\t'
            << protein_accessions << '\t'
            << joined_gene_names << '\t'
            << fragment.product_mz << '\t'
            << fragment.product_charge << '\t'
            << fragment.product_type << '\t'
            << fragment.product_ordinal << '\n';
      }
    }
  }

  static void writeStage2ScoreTable_(const String& filename,
                                     const std::vector<FastaEvidenceFilter::Stage2CandidateScore>& scores)
  {
    std::ofstream out(filename.c_str());
    if (!out)
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    out << "protein_accession\tgene_name\tpeptide_sequence\tmodified_peptide_sequence\tprecursor_mz\tprecursor_charge\tdecoy\tsource_file\tnative_spectrum_id\tbest_matched_ions\tsupporting_spectra\tsupporting_runs\tstrong_supporting_spectra\tstrong_supporting_runs\truns_with_streak_ge_2\truns_with_streak_ge_3\tbest_spectrum_matched_intensity_fraction\tbest_spectrum_matched_b_ions\tbest_spectrum_matched_y_ions\tbest_spectrum_longest_b_run\tbest_spectrum_longest_y_run\tbest_spectrum_longest_y_pct\tbest_spectrum_poisson_proxy\tbest_spectrum_score\tbest_run_streak_length\tbest_run_streak_score\ttop_run_score_1\ttop_run_score_2\ttop_run_score_3\ttop_run_streak_length_1\ttop_run_streak_length_2\ttop_run_streak_length_3\tbest_local_rank\tbest_local_pvalue\tcombined_pvalue\tcomposite_score\tqvalue\taccepted\n";
    out << std::fixed << std::setprecision(6);
    for (const auto& score : scores)
    {
      std::vector<std::string> gene_names;
      gene_names.reserve(score.protein_refs.size());
      for (const auto& protein_ref : score.protein_refs)
      {
        const auto gene_name_it = score.protein_gene_names_by_accession.find(protein_ref);
        gene_names.push_back(gene_name_it != score.protein_gene_names_by_accession.end() ? gene_name_it->second : std::string{});
      }
      const String protein_accessions = ListUtils::concatenate(score.protein_refs, ";");
      const String joined_gene_names = ListUtils::concatenate(gene_names, ";");

      out << protein_accessions << '\t'
          << joined_gene_names << '\t'
          << score.peptide_sequence << '\t'
          << score.modified_peptide_sequence << '\t'
          << score.precursor_mz << '\t'
          << score.precursor_charge << '\t'
          << (score.decoy ? 1 : 0) << '\t'
          << score.source_file << '\t'
          << score.native_spectrum_id << '\t'
          << score.best_matched_ions << '\t'
          << score.supporting_spectra << '\t'
          << score.supporting_runs << '\t'
          << score.strong_supporting_spectra << '\t'
          << score.strong_supporting_runs << '\t'
          << score.runs_with_streak_ge_2 << '\t'
          << score.runs_with_streak_ge_3 << '\t'
          << score.best_spectrum_matched_intensity_fraction << '\t'
          << score.best_spectrum_matched_b_ions << '\t'
          << score.best_spectrum_matched_y_ions << '\t'
          << score.best_spectrum_longest_b_run << '\t'
          << score.best_spectrum_longest_y_run << '\t'
          << score.best_spectrum_longest_y_pct << '\t'
          << score.best_spectrum_poisson_proxy << '\t'
          << score.best_spectrum_score << '\t'
          << score.best_run_streak_length << '\t'
          << score.best_run_streak_score << '\t'
          << score.top_run_score_1 << '\t'
          << score.top_run_score_2 << '\t'
          << score.top_run_score_3 << '\t'
          << score.top_run_streak_length_1 << '\t'
          << score.top_run_streak_length_2 << '\t'
          << score.top_run_streak_length_3 << '\t'
          << score.best_local_rank << '\t';

      if (score.best_local_pvalue >= 0.0)
      {
        out << score.best_local_pvalue;
      }
      out << '\t';
      if (score.combined_pvalue >= 0.0)
      {
        out << score.combined_pvalue;
      }
      out << '\t'
          << score.composite_score << '\t';

      if (score.qvalue >= 0.0)
      {
        out << score.qvalue;
      }
      out << '\t'
          << (score.accepted ? 1 : 0) << '\n';
    }
  }
};
/// @endcond

int main(int argc, const char** argv)
{
  TOPPFastaEvidenceFilter tool;
  return tool.main(argc, argv);
}
