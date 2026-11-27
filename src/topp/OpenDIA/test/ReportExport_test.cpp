// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $

#include <results/ReportExport.h>
#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/file_reader.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace OpenMS::OpenDIA;
namespace fs = std::filesystem;

namespace
{
  void require(bool condition, const std::string& message)
  {
    if (!condition) throw std::runtime_error(message);
  }

  template <typename T> T unwrap(arrow::Result<T> result)
  {
    if (!result.ok()) throw std::runtime_error(result.status().ToString());
    return std::move(result).ValueOrDie();
  }

  void check(const arrow::Status& status)
  {
    require(status.ok(), status.ToString());
  }

  // Parse quoted TSV, including embedded line breaks and doubled quotes. Keeping an independent
  // reader here checks data recovery, rather than just matching a writer-specific escaped string.
  std::vector<std::vector<std::string>> tsv(const fs::path& path)
  {
    std::ifstream stream(path, std::ios::binary);
    require(stream.good(), "read TSV " + path.string());
    std::vector<std::vector<std::string>> rows;
    std::vector<std::string> row;
    std::string cell;
    bool quoted = false;
    char c;
    while (stream.get(c))
    {
      if (c == '"')
      {
        if (quoted && stream.peek() == '"') { stream.get(c); cell += '"'; }
        else quoted = !quoted;
      }
      else if (!quoted && (c == '\t' || c == '\n'))
      {
        row.push_back(cell);
        cell.clear();
        if (c == '\n') { rows.push_back(row); row.clear(); }
      }
      else cell += c;
    }
    require(!quoted && row.empty() && cell.empty(), "complete quoted TSV");
    return rows;
  }

  std::size_t column(const std::vector<std::string>& header, const std::string& name)
  {
    for (std::size_t i = 0; i < header.size(); ++i) if (header[i] == name) return i;
    throw std::runtime_error("Missing TSV column " + name);
  }

  std::shared_ptr<arrow::Table> readParquet(const fs::path& path)
  {
    auto file = unwrap(arrow::io::ReadableFile::Open(path.string()));
    auto reader = unwrap(parquet::arrow::OpenFile(file, arrow::default_memory_pool()));
    std::shared_ptr<arrow::Table> table;
    check(reader->ReadTable(&table));
    return unwrap(table->CombineChunks());
  }

  std::shared_ptr<arrow::DoubleArray> doubles(const std::shared_ptr<arrow::Table>& table, const std::string& name)
  {
    return std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName(name)->chunk(0));
  }

  std::shared_ptr<arrow::StringArray> strings(const std::shared_ptr<arrow::Table>& table, const std::string& name)
  {
    return std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName(name)->chunk(0));
  }

  void noPartials(const fs::path& directory)
  {
    for (const auto& entry : fs::directory_iterator(directory))
      require(entry.path().filename().string().find(".partial.") == std::string::npos, "no stale partial files");
  }

  void addPrecursor(AnalysisLibrary& library, const std::string& id, bool decoy)
  {
    const auto p = static_cast<std::uint32_t>(library.precursors.size());
    const auto it = library.index.emplace(id, p).first;
    AnalysisLibrary::Precursor precursor;
    precursor.id = it->first;
    precursor.sequence = "PEPM(Oxidation)IDE";
    precursor.charge = 2;
    precursor.decoy = decoy;
    precursor.proteins = 0;
    library.precursors.push_back(precursor);
  }

  RunWindow window(std::uint32_t p, float score, double q, double area)
  {
    RunWindow result;
    result.precursor = p;
    result.candidate = p + 10;
    result.score = score;
    result.q_value = q;
    result.peak.area = area;
    result.peak.rt_start = 10.123456789012345;
    result.peak.rt_apex = 11.123456789012345;
    result.peak.rt_end = 12.123456789012345;
    return result;
  }

  void exercise(const fs::path& root)
  {
    const double na = std::numeric_limits<double>::quiet_NaN();
    AnalysisLibrary library;
    library.protein_sets = {"P1;P2"};
    addPrecursor(library, "p0", false);
    addPrecursor(library, "p1", false);
    addPrecursor(library, "DECOY_p0", true);
    addPrecursor(library, "p3", false);
    std::vector<RunResult> runs(2);
    runs[0].run = "/input/runA.mzML.gz";
    runs[1].run = "/input/runB.d";
    const double precise = 123456.78901234567;
    runs[0].windows = {window(3, 6, na, 8), window(2, 5, 0, 0), window(1, 7, 0.001, precise),
                       window(0, 10, 0.01, 0), window(0, 1, 0.001, 999),
                       window(1, CrossFitScorer::NOT_KEPT, 0, 111)};
    runs[1].windows = {window(0, 9, 0.02, precise), window(1, 8, 0.001, 22)};
    Quantification quantification;
    quantification.global_precursor_q = {0.01, 0.02, 0, na};
    quantification.precursor_peptide = {0, 1, 2, 3};
    quantification.precursor_protein_group = {0, UINT32_MAX, 1, UINT32_MAX};
    quantification.peptides = {
      {"PEPMIDE", {"P1", "P2"}, false, 0.01, {{0, 10, 0.02, 0, true}, {1, 9, 0.01, 0, false}}},
      {"OTHER", {"P1"}, false, 0.02, {{0, 7, 0.01, precise, true}, {1, 8, 0.02, 22, true}}},
      {"DECOY_PEPMIDE", {"P1"}, true, 0, {{0, 5, 0, 0, true}}},
      {"UNKNOWN", {}, false, na, {{0, 6, na, 0, false}}},
      {"HIGH_GLOBAL", {}, false, 0.021, {{0, 6, 0.01, 9, true}}}
    };
    quantification.protein_groups = {
      {"GROUP1", {"P1", "P2"}, false, 0.03, {{0, 10, 0.03, 0, true}, {1, 9, 0.031, 19, true}}},
      {"DECOY_P1;P2", {"P1", "P2"}, true, 0, {{0, 5, 0, 0, true}}},
      {"MISSING", {}, false, 0, {{1, 3, 0, 0, false}}}
    };
    // Distinct contexts make accidental relabeling and cross-level gating observable.
    quantification.peptides[0].id = "PEPM(Oxidation)IDE";
    auto& peptide_evidence = quantification.peptides[0].runs[0];
    peptide_evidence.p_value = 0;
    peptide_evidence.pep = 0.12345678901234567;
    peptide_evidence.experiment_q_value = 0.5;
    peptide_evidence.experiment_p_value = 0.25;
    peptide_evidence.experiment_pep = 0;
    quantification.peptides[0].global_p_value = 0;
    quantification.peptides[0].global_pep = 0.125;
    quantification.peptides[1].runs[0].q_value = 0.5;
    quantification.peptides[1].runs[1].q_value = 0.5;
    quantification.peptides[1].runs[0].experiment_q_value = 0.01;
    quantification.peptides[4].runs[0].experiment_q_value = 0;
    quantification.protein_groups[0].runs[0].experiment_q_value = 0.02;
    const GeneAnnotations genes{{"p0", "GENE\t\"quoted\"\nsecond line"}, {"p1", "NA"}};
    const ReportingOptions options{0.01, 0.02, 0.03};
    const auto both = root / "both";
    writeReports(both.string(), library, runs, quantification, genes, options, "both");
    auto precursor = tsv(both / "precursors.tsv");
    require(precursor.size() == 7, "all six retained winners including decoy, no q filtering or duplicate candidates");
    const auto& header = precursor.front();
    auto value = [&](std::size_t row, const std::string& name) -> const std::string& { return precursor.at(row).at(column(header, name)); };
    require(value(1, "precursor_id") == "p0" && value(2, "run") == "runB", "library then run deterministic order");
    require(value(1, "sequence") == "PEPM(Oxidation)IDE" && value(1, "unmodified_sequence") == "PEPMIDE", "source and unmodified sequence");
    require(value(1, "genes") == genes.at("p0") && value(3, "genes") == "NA", "quoted gene round trip");
    require(value(1, "feature_id") == std::to_string(outputFeatureId(runs[0].run, "p0", 10)), "stable producer feature ID");
    require(std::stod(value(3, "area")) == precise, "full precision TSV area");
    require(std::stod(value(1, "rt_apex")) == runs[0].windows[3].peak.rt_apex, "full precision TSV RT");
    require(value(1, "pass_peakgroup") == "true" && value(1, "pass_peptide") == "true" && value(1, "pass_protein_group") == "true", "inclusive independent thresholds");
    require(value(2, "pass_peakgroup") == "false" && value(2, "pass_peptide") == "false" && value(2, "pass_protein_group") == "false", "cumulative feature gates do not alter independent matrix gates");
    require(value(3, "protein_group_run_q_value") == "NA" && value(3, "protein_group_global_q_value") == "NA", "ambiguous protein is unknown, not zero q");
    require(value(5, "decoy") == "true" && value(5, "pass_peakgroup") == "false", "decoy retained only in long table");
    require(value(6, "run_q_value") == "NA" && value(6, "genes") == "NA", "missing q and genes stay missing");
    require(value(1, "quantity_type") == "raw_ms2_area", "quantity semantics");

    auto matrix = tsv(both / "precursor_matrix.tsv");
    require(matrix.size() == 2 && matrix[1][0] == "p0", "matrix excludes decoys and rows without eligible cells");
    require(matrix[1][column(matrix[0], "runA")] == "0" && matrix[1][column(matrix[0], "runB")] == "NA", "real zero versus absent TSV");
    auto peptide = tsv(both / "peptide_matrix.tsv");
    require(peptide.size() == 2 && peptide[1][0] == "PEPM(Oxidation)IDE", "peptide input entity order and level-specific gate");
    require(peptide[1][2] == "0" && peptide[1][3] == "NA", "use provided entity area and has_area");
    auto protein = tsv(both / "protein_group_matrix.tsv");
    require(protein.size() == 2 && protein[1][0] == "P1;P2" && protein[1][1] == "P1;P2", "complete protein group membership");
    require(protein[1][2] == "0" && protein[1][3] == "NA", "protein independent threshold");
    auto long_peptide = tsv(both / "peptides.tsv");
    require(long_peptide.size() == 8, "entity long table retains failed and decoy observations");
    require(long_peptide[2][column(long_peptide[0], "area")] == "NA", "missing entity area not fabricated");

    auto context_options = options;
    context_options.inference_context = "experiment-wide";
    const auto context_output = root / "experiment_context";
    writeReports(context_output.string(), library, runs, quantification, genes, context_options, "both");
    const auto context_precursors = tsv(context_output / "precursors.tsv");
    require(context_precursors[1][column(context_precursors[0], "pass_peptide")] == "false", "experiment q gates peptide flag");
    require(context_precursors[1][column(context_precursors[0], "pass_protein_group")] == "true", "protein is not gated by peptide experiment q");
    require(context_precursors[1][column(context_precursors[0], "inference_context")] == "experiment-wide", "selected context explicit");
    require(context_precursors[1][column(context_precursors[0], "peptide_run_q_value")] == precursor[1][column(precursor[0], "peptide_run_q_value")], "run q never relabeled");
    const auto context_entities = tsv(context_output / "peptides.tsv");
    require(context_entities.back()[column(context_entities[0], "pass")] == "false", "global entity cut still applies in experiment context");
    const auto context_matrix = tsv(context_output / "peptide_matrix.tsv");
    require(context_matrix.size() == 2 && context_matrix[1][0] == "OTHER", "experiment context selects divergent peptide evidence");
    require(tsv(context_output / "precursor_matrix.tsv") == matrix, "precursor matrix context invariant");
    require(tsv(context_output / "protein_group_matrix.tsv") == protein, "protein matrix independent of peptide context failure");
    const auto entity_table = readParquet(both / "peptides.parquet");
    const auto precursor_table = readParquet(both / "precursors.parquet");
    for (const auto& entry : std::vector<std::pair<std::string, double>>{{"run_p_value", 0}, {"run_pep", 0.12345678901234567}, {"experiment_q_value", 0.5}, {"experiment_p_value", 0.25}, {"experiment_pep", 0}, {"global_p_value", 0}, {"global_pep", 0.125}})
    {
      require(doubles(entity_table, entry.first)->Value(0) == entry.second, "exact entity Parquet statistic " + entry.first);
      require(doubles(precursor_table, "peptide_" + entry.first)->Value(0) == entry.second, "exact precursor Parquet statistic " + entry.first);
      require(std::stod(long_peptide[1][column(long_peptide[0], entry.first)]) == entry.second, "exact TSV statistic " + entry.first);
    }
    require(doubles(entity_table, "run_p_value")->IsNull(1), "unavailable p stays Parquet null");
    require(long_peptide[2][column(long_peptide[0], "run_p_value")] == "NA", "unavailable p stays TSV NA");

    auto invalid_context = options;
    invalid_context.inference_context = "global";
    bool context_failed = false;
    try { writeReports((root / "bad_context").string(), library, runs, quantification, {}, invalid_context, "both"); }
    catch (const std::invalid_argument&) { context_failed = true; }
    require(context_failed && !fs::exists(root / "bad_context"), "invalid context rejected before output");

    auto pm = readParquet(both / "precursor_matrix.parquet");
    require(pm->num_rows() == 1, "Parquet matrix row count");
    require(!doubles(pm, "runA")->IsNull(0) && doubles(pm, "runA")->Value(0) == 0 && doubles(pm, "runB")->IsNull(0), "Parquet zero versus null");
    auto pp = readParquet(both / "precursors.parquet");
    require(pp->num_rows() == 6 && doubles(pp, "area")->Value(2) == precise, "Parquet full precision and all winner rows");
    require(strings(pp, "genes")->GetString(0) == genes.at("p0") && strings(pp, "genes")->IsNull(4), "Parquet string and gene null");
    require(doubles(pp, "protein_group_global_q_value")->IsNull(2), "Parquet ambiguous protein q null");
    require(doubles(pp, "run_q_value")->IsNull(5), "Parquet unknown q null");
    auto ep = readParquet(both / "peptides.parquet");
    require(doubles(ep, "area")->IsNull(1), "Parquet entity has_area false null");
    noPartials(both);

    for (const std::string format : {"tsv", "parquet"})
    {
      const auto output = root / format;
      writeReports(output.string(), library, runs, quantification, {}, options, format);
      require(fs::exists(output / ("precursors." + format)), "requested format exists");
      require(!fs::exists(output / (format == "tsv" ? "precursors.parquet" : "precursors.tsv")), "unrequested format absent");
      noPartials(output);
    }
    AnalysisLibrary empty_library;
    const auto empty = root / "empty";
    writeReports(empty.string(), empty_library, {}, {}, {}, {}, "both");
    for (const std::string name : {"precursors", "precursor_matrix", "peptide_matrix", "protein_group_matrix", "peptides", "protein_groups"})
    {
      require(tsv(empty / (name + ".tsv")).size() == 1, "empty TSV has schema only");
      require(readParquet(empty / (name + ".parquet"))->num_rows() == 0, "empty Parquet has schema only");
    }

    bool failed = false;
    try { writeReports((root / "invalid").string(), library, runs, quantification, {}, options, "json"); }
    catch (const std::invalid_argument&) { failed = true; }
    require(failed && !fs::exists(root / "invalid"), "format rejection before output");

    // Fail late, after all precursor files have been completed. Prior published reports remain
    // untouched and no partial files survive the exception.
    const auto before = precursor;
    auto bad = quantification;
    bad.peptides[0].runs.push_back({999, 1, 0, 0, true});
    failed = false;
    try { writeReports(both.string(), library, runs, bad, genes, options, "both"); }
    catch (const std::out_of_range&) { failed = true; }
    require(failed && tsv(both / "precursors.tsv") == before, "failed generation preserves published output");
    noPartials(both);

    const auto rename_failure = root / "rename_failure";
    fs::create_directories(rename_failure / "precursors.tsv");
    failed = false;
    try { writeReports(rename_failure.string(), library, runs, quantification, {}, options, "both"); }
    catch (const fs::filesystem_error&) { failed = true; }
    require(failed && fs::is_directory(rename_failure / "precursors.tsv"), "publication errors propagate without deleting destination");
    noPartials(rename_failure);

    failed = false;
    auto invalid_options = options;
    invalid_options.peptide_fdr = na;
    try { writeReports((root / "bad_threshold").string(), library, runs, quantification, {}, invalid_options, "tsv"); }
    catch (const std::invalid_argument&) { failed = true; }
    require(failed && !fs::exists(root / "bad_threshold"), "invalid threshold rejected before output");

    auto duplicate_runs = runs;
    duplicate_runs[1].run = "/elsewhere/runA.raw";
    failed = false;
    try { writeReports((root / "duplicate").string(), library, duplicate_runs, quantification, {}, options, "tsv"); }
    catch (const std::invalid_argument&) { failed = true; }
    require(failed && !fs::exists(root / "duplicate"), "ambiguous run stems rejected before output");

    auto independent_levels = quantification;
    independent_levels.peptides[0].global_q_value = 0.021;
    const auto independent = root / "independent_levels";
    writeReports(independent.string(), library, runs, independent_levels, {}, options, "both");
    const auto independent_table = tsv(independent / "precursors.tsv");
    require(independent_table[1][column(independent_table[0], "pass_peakgroup")] == "true"
         && independent_table[1][column(independent_table[0], "pass_peptide")] == "false"
         && independent_table[1][column(independent_table[0], "pass_protein_group")] == "true",
            "protein pass requires peak and protein cuts, never the independent peptide cut");
    require(tsv(independent / "protein_group_matrix.tsv").size() == 2, "peptide q failure does not hide eligible protein row");
    noPartials(independent);

    const auto malformed = root / "malformed_sequence";
    library.precursors[0].sequence = "PEP(UnknownModificationThatDoesNotExist)TIDE";
    writeReports(malformed.string(), library, runs, quantification, {}, options, "both");
    const auto malformed_table = readParquet(malformed / "precursors.parquet");
    require(strings(malformed_table, "sequence")->GetString(0) == library.precursors[0].sequence
         && strings(malformed_table, "unmodified_sequence")->IsNull(0), "unparseable source sequence retained with nullable normalized metadata");
    library.precursors[0].sequence.clear();
    const auto empty_sequence = root / "empty_sequence";
    writeReports(empty_sequence.string(), library, runs, quantification, {}, options, "parquet");
    require(strings(readParquet(empty_sequence / "precursors.parquet"), "unmodified_sequence")->IsNull(0), "empty normalized sequence is null");
    library.precursors[0].sequence = "PEPM(Oxidation)IDE";
    noPartials(malformed);

    auto reserved_runs = runs;
    reserved_runs[0].run = "/input/proteins.mzML.gz";
    failed = false;
    try { writeReports((root / "reserved").string(), library, reserved_runs, quantification, {}, options, "both"); }
    catch (const std::invalid_argument&) { failed = true; }
    require(failed && !fs::exists(root / "reserved"), "run stem cannot collide with matrix metadata");

    AnalysisLibrary colliding_library;
    colliding_library.protein_sets = {"P1"};
    addPrecursor(colliding_library, "same", false);
    addPrecursor(colliding_library, "same", false);
    RunResult colliding_run;
    colliding_run.run = "/input/collision.mzML.gz";
    colliding_run.windows = {window(0, 1, 0, 1), window(1, 1, 0, 1)};
    colliding_run.windows[1].candidate = colliding_run.windows[0].candidate;
    failed = false;
    const auto collision_output = root / "collision";
    try { writeReports(collision_output.string(), colliding_library, {colliding_run}, {}, {}, options, "both"); }
    catch (const std::runtime_error&) { failed = true; }
    require(failed && !fs::exists(collision_output / "precursors.tsv"), "duplicate/colliding feature IDs fail before publication");
    noPartials(collision_output);

    // Force multiple bounded row groups without a dense full matrix. Entity order stays intact.
    Quantification many;
    for (int i = 0; i < 8200; ++i) many.peptides.push_back({"PEPTIDE_" + std::to_string(i), {}, false, 0, {{0, 1, 0, double(i), true}}});
    const auto bounded = root / "bounded";
    RunResult sparse_run;
    sparse_run.run = runs[0].run;
    writeReports(bounded.string(), empty_library, {sparse_run}, many, {}, options, "parquet");
    auto reader = parquet::ParquetFileReader::OpenFile((bounded / "peptide_matrix.parquet").string());
    const auto meta = reader->metadata();
    require(meta->num_rows() == 8200 && meta->num_row_groups() >= 3, "multiple bounded row groups");
    for (int i = 0; i < meta->num_row_groups(); ++i) require(meta->RowGroup(i)->num_rows() <= 4096, "row group cap");
    const auto wide = readParquet(bounded / "peptide_matrix.parquet");
    require(strings(wide, "modified_sequence")->GetString(8199) == "PEPTIDE_8199" && doubles(wide, "runA")->Value(8199) == 8199, "multi-batch row order and numeric round trip");
    noPartials(bounded);
  }
}

int main()
{
  const auto root = fs::temp_directory_path() / ("OpenDIA_ReportExport_test_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  try
  {
    fs::create_directories(root);
    exercise(root);
    fs::remove_all(root);
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << "\nTest data: " << root << '\n';
    fs::remove_all(root);
    return 1;
  }
}
