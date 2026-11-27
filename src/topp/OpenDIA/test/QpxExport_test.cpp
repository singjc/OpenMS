// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/QPXIdentity.h>
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <parquet/arrow/reader.h>
#include <parquet/file_reader.h>
#include <results/QpxExport.h>
#include <stdexcept>

using namespace OpenMS;
using namespace OpenMS::OpenDIA;
namespace fs = std::filesystem;
namespace
{
void require(bool condition, const std::string& message)
{
  if (! condition) throw std::runtime_error(message);
}
template<typename T>
T take(arrow::Result<T> result)
{
  if (! result.ok()) throw std::runtime_error(result.status().ToString());
  return std::move(result).ValueOrDie();
}
std::shared_ptr<arrow::Table> read(const fs::path& path)
{
  auto input = take(arrow::io::ReadableFile::Open(path.string()));
  auto reader = take(parquet::arrow::OpenFile(input, arrow::default_memory_pool()));
  std::shared_ptr<arrow::Table> table;
  auto status = reader->ReadTable(&table);
  require(status.ok(), status.ToString());
  require(table->ValidateFull().ok(), "invalid linked-reader table");
  return take(table->CombineChunks());
}
template<typename Array>
std::shared_ptr<Array> column(const std::shared_ptr<arrow::Table>& table, const std::string& name)
{
  auto col = table->GetColumnByName(name);
  require(col && col->num_chunks() == 1, "missing column " + name);
  return std::static_pointer_cast<Array>(col->chunk(0));
}
double extensionNumber(const std::shared_ptr<arrow::Table>& table, const std::string& name, std::int64_t row)
{ return std::stod(column<arrow::StringArray>(table, "x_opendia_" + name)->GetString(row)); }
std::shared_ptr<arrow::FloatArray> intensities(const std::shared_ptr<arrow::Table>& table)
{
  auto values = std::static_pointer_cast<arrow::StructArray>(column<arrow::ListArray>(table, "intensities")->values());
  return std::static_pointer_cast<arrow::FloatArray>(values->field(1));
}
bool exact(double a, double b)
{ return std::memcmp(&a, &b, sizeof(double)) == 0; }
AnalysisLibrary library()
{
  AnalysisLibrary lib;
  lib.index.emplace("original_42", 0);
  lib.index.emplace("DECOY_original_42", 1);
  for (const auto& id : {"original_42", "DECOY_original_42"})
  {
    AnalysisLibrary::Precursor p;
    p.id = lib.index.find(id)->first;
    p.sequence = "PEPM(Oxidation)IDEK";
    p.charge = 2;
    p.mz = 512.1234567890123;
    p.has_library_rt = true;
    p.library_rt = 17.123456789;
    p.decoy = lib.precursors.size() == 1;
    lib.precursors.push_back(p);
  }
  return lib;
}
RunResult run()
{
  RunResult result;
  result.run = "/raw/acquisition.mzML";
  result.calibration_json = "{\"rt\":1.23456789012345}";
  result.feature_names = {"scorer_a", "scorer_b"};
  for (std::uint32_t i = 0; i < 3; ++i)
  {
    RunWindow w;
    w.precursor = i == 2 ? 1 : 0;
    w.candidate = i;
    w.peak.rt_apex = 1000.123456789 + static_cast<double>(i) * 1e-8;
    w.peak.rt_start = 998.123456789;
    w.peak.rt_end = 1002.987654321;
    w.peak.area = i == 0 ? 0.0 : 123456789.12345679;
    w.score = i == 1 ? 10.0f : 2.0f;
    w.q_value = 0.003 + i * 0.001;
    w.network_score = 1.2345678f + i;
    w.network_q_value = 0.023 + i * 0.001;
    result.windows.push_back(w);
    result.features.push_back(std::nextafter(0.1234567f, 1.0f));
    result.features.push_back(-2.345678f - i);
  }
  return result;
}
} // namespace

int main(int argc, char** argv)
{
  try
  {
    require(argc > 1, "test needs an output directory");
    const fs::path output = argv[1];
    fs::create_directories(output);
    auto lib = library();
    auto sample = run();
    require(outputFeatureId(sample.run, "original_42", 0)
              == QPXIdentity::deriveId(
                {QPXIdentity::Value {std::string("acquisition")}, QPXIdentity::Value {std::string("original_42")}, QPXIdentity::Value {Int64 {0}}}),
            "identity composite");
    require(outputFeatureId(sample.run, "original_42", 0) != outputFeatureId(sample.run, "original_42", 1), "candidate identity collision");
    writeCandidateEvidence(output.string(), lib, sample);
    const auto candidate_path = output / "opendia/acquisition.candidates.parquet";
    auto candidates = read(candidate_path);
    require(candidates->num_rows() == 3, "all candidate rows retained");
    require(candidates->schema()->GetFieldByName("scorer_a")->type()->id() == arrow::Type::FLOAT, "scorer float32 schema");
    for (int i = 0; i < 3; ++i)
    {
      require(column<arrow::FloatArray>(candidates, "scorer_a")->Value(i) == sample.features[i * 2], "scorer float bits");
      require(exact(column<arrow::DoubleArray>(candidates, "opendia_area")->Value(i), sample.windows[i].peak.area), "candidate area precision");
      require(exact(column<arrow::DoubleArray>(candidates, "opendia_rt_apex")->Value(i), sample.windows[i].peak.rt_apex), "candidate RT precision");
    }
    require(! column<arrow::BooleanArray>(candidates, "opendia_winner")->Value(0)
              && column<arrow::BooleanArray>(candidates, "opendia_winner")->Value(1),
            "winner flag");
    auto footer = candidates->schema()->metadata();
    require(take(footer->Get("identity_composite")) == "run_file_name,opendia_precursor_id,opendia_candidate_row", "identity metadata");
    require(take(footer->Get("opendia_calibration_json")) == sample.calibration_json, "calibration metadata");

    Quantification quant;
    quant.global_precursor_q = {0.123, 0.234};
    quant.precursor_peptide = {0, UINT32_MAX};
    quant.precursor_protein_group = {0, UINT32_MAX};
    QuantifiedEntity peptide;
    peptide.id = "PEPMIDEK";
    peptide.global_q_value = 0.345;
    peptide.runs.push_back({0, 10, 0.456, 0, false});
    quant.peptides.push_back(peptide);
    QuantifiedEntity group;
    group.id = "P1;P2";
    group.proteins = {"P1", "P2"};
    group.global_q_value = 0.567;
    group.runs.push_back({0, 10, 0.678, 0, false});
    quant.protein_groups.push_back(group);
    GeneAnnotations genes {{"original_42", "GENE1; GENE2;GENE1;"}};
    writeQpxDataset(output.string(), lib, {sample}, quant, genes, {});
    auto features = read(output / "opendia.feature.parquet");
    auto feature_footer = features->schema()->metadata();
    require(take(feature_footer->Get("identity_source")) == "supplied", "canonical producer supplied identity");
    require(take(feature_footer->Get("opendia_identity_composite")) == "run_file_name,opendia_precursor_id,opendia_candidate_row",
            "producer identity retained");
    require(features->schema()->GetFieldByName("x_opendia_area")->type()->id() == arrow::Type::STRING,
            "canonical extensions are standard nullable strings");
    require(exact(extensionNumber(features, "area", 1), sample.windows[1].peak.area), "canonical string area exact double");
    require(exact(extensionNumber(features, "rt_apex", 1), sample.windows[1].peak.rt_apex), "canonical string apex exact double");
    require(features->num_rows() == 3, "canonical all candidates including decoys");
    require(column<arrow::BooleanArray>(features, "is_decoy")->Value(2), "decoy retained");
    require(column<arrow::Int64Array>(features, "feature_id")->Value(0) != column<arrow::Int64Array>(features, "feature_id")->Value(1),
            "rounded RT no collision");
    require(column<arrow::FloatArray>(features, "rt")->Value(0) == column<arrow::FloatArray>(features, "rt")->Value(1), "fixture RT rounds together");
    require(intensities(features)->IsNull(0) && ! intensities(features)->IsNull(1), "nonwinner intensity null");
    require(column<arrow::FloatArray>(features, "observed_mz")->null_count() == 3, "no fabricated observed m/z");
    require(column<arrow::DoubleArray>(features, "posterior_error_probability")->null_count() == 3, "no fabricated PEP");
    require(column<arrow::StringArray>(features, "sequence")->GetString(0) == "PEPMIDEK", "unmodified sequence");
    require(column<arrow::StringArray>(features, "peptidoform")->GetString(0).find("UNIMOD:35") != std::string::npos, "canonical ProForma");
    require(column<arrow::ListArray>(features, "pg_accessions")->IsNull(2), "unresolved group does not use ambiguous union");
    require(column<arrow::ListArray>(features, "gg_names")->IsNull(2), "gene map keyed by precursor ID");
    require(column<arrow::ListArray>(features, "gg_names")->value_length(0) == 2, "gene annotations split and deduplicated");
    require(column<arrow::FloatArray>(features, "predicted_rt")->null_count() == 3, "source iRT is not predicted seconds");
    require(exact(extensionNumber(features, "library_rt", 0), lib.precursors[0].library_rt), "raw library RT preserved");
    require(column<arrow::DoubleArray>(features, "global_qvalue")->Value(0) == 0.123, "global precursor q independent");
    require(column<arrow::DoubleArray>(features, "precursor_qvalue")->Value(0) == sample.windows[0].q_value, "run precursor q independent");
    require(column<arrow::DoubleArray>(features, "peptide_qvalue")->Value(0) == 0.456, "peptide q independent");
    require(extensionNumber(features, "pg_qvalue", 0) == 0.678, "protein q independent");
    require(exact(extensionNumber(features, "mz", 0), lib.precursors[0].mz), "mz exact double");
    require(exact(extensionNumber(features, "rt_start", 0), sample.windows[0].peak.rt_start), "boundary exact double");
    require(exact(extensionNumber(features, "rt_end", 0), sample.windows[0].peak.rt_end), "end exact double");
    auto pg = read(output / "opendia.pg.parquet");
    require(pg->num_rows() == 1 && column<arrow::FloatArray>(pg, "intensity")->IsNull(0), "PG missing area is null");
    require(column<arrow::DoubleArray>(pg, "global_qvalue")->Value(0) == 0.567 && column<arrow::DoubleArray>(pg, "pg_qvalue")->Value(0) == 0.678,
            "PG global/run q");
    auto run_view = read(output / "opendia.run.parquet");
    require(column<arrow::StringArray>(run_view, "opendia_calibration_json")->GetString(0) == sample.calibration_json, "run calibration");
    auto sample_view = read(output / "opendia.sample.parquet");
    require(sample_view->num_rows() == 1 && column<arrow::StringArray>(sample_view, "sample_accession")->GetString(0) == "opendia_run:acquisition",
            "technical sample mapping");
    require(column<arrow::ListArray>(run_view, "samples")->value_length(0) == 1, "run maps one LFQ sample");
    require(column<arrow::StringArray>(run_view, "instrument")->IsNull(0), "unavailable instrument remains null");
    require(fs::exists(output / "opendia.provenance.parquet") && ! fs::exists(output / "opendia.psm.parquet")
              && ! fs::exists(output / "opendia.mz.parquet"),
            "honest dataset views");

    // Target and decoy entities legitimately share resolved protein accessions.
    auto paired = quant;
    auto decoy_group = group;
    decoy_group.id = "DECOY_" + group.id;
    decoy_group.decoy = true;
    decoy_group.global_q_value = std::numeric_limits<double>::quiet_NaN();
    decoy_group.runs[0].q_value = std::numeric_limits<double>::quiet_NaN();
    paired.protein_groups.push_back(decoy_group);
    paired.precursor_protein_group[1] = 1;
    QuantifiedEntity prefix_collision;
    prefix_collision.id = "DECOY_P1";
    prefix_collision.proteins = {"DECOY_P1"};
    paired.protein_groups.push_back(prefix_collision);
    writeQpxDataset((output / "paired").string(), lib, {sample}, paired, {}, {});
    auto paired_pg = read(output / "paired/opendia.pg.parquet");
    require(paired_pg->num_rows() == 2, "target and decoy PG retained");
    require(column<arrow::Int64Array>(paired_pg, "pg_id")->Value(0) != column<arrow::Int64Array>(paired_pg, "pg_id")->Value(1),
            "decoy-aware protein identity");
    require(column<arrow::StringArray>(paired_pg, "anchor_protein")->GetString(1) == "DECOY_DECOY_P1",
            "derived decoy accession avoids real target collision");
    auto paired_features = read(output / "paired/opendia.feature.parquet");
    require(column<arrow::StringArray>(paired_features, "anchor_protein")->GetString(2) == "DECOY_DECOY_P1", "feature decoy membership matches PG");
    require(column<arrow::DoubleArray>(paired_pg, "global_qvalue")->IsNull(1) && column<arrow::DoubleArray>(paired_pg, "pg_qvalue")->IsNull(1),
            "unestimable q values are null, not NaN");

    // A positive-quantity, fully resolved target/decoy fixture for the unmodified Python conformance gate.
    auto conformant = paired;
    conformant.protein_groups.pop_back();
    for (auto& entity : conformant.protein_groups)
    {
      entity.global_q_value = 0.01;
      entity.runs[0].q_value = 0.01;
      entity.runs[0].area = 123456789.12345679;
      entity.runs[0].has_area = true;
    }
    writeQpxDataset((output / "conformance").string(), lib, {sample}, conformant, {}, {});

    // Zero is a measured area, not a missing value.
    quant.protein_groups[0].runs[0].has_area = true;
    writeQpxDataset((output / "zero").string(), lib, {sample}, quant, {}, {});
    auto zero_pg = read(output / "zero/opendia.pg.parquet");
    require(! column<arrow::FloatArray>(zero_pg, "intensity")->IsNull(0) && column<arrow::FloatArray>(zero_pg, "intensity")->Value(0) == 0,
            "zero PG area retained");

    RunResult empty;
    empty.run = "empty.mzML";
    writeCandidateEvidence((output / "empty").string(), lib, empty);
    writeQpxDataset((output / "empty").string(), lib, {empty}, {}, {}, {});
    require(read(output / "empty/opendia.feature.parquet")->num_rows() == 0, "empty canonical schema");
    require(read(output / "empty/opendia.pg.parquet")->num_rows() == 0, "empty PG schema");
    require(read(output / "empty/opendia/empty.candidates.parquet")->num_rows() == 0, "empty scorer schema");

    writeQpxDataset((output / "no_runs").string(), lib, {}, {}, {}, {});
    require(read(output / "no_runs/opendia.feature.parquet")->num_rows() == 0, "no-run feature schema");
    require(read(output / "no_runs/opendia.run.parquet")->num_rows() == 0, "no-run metadata schema");

    // A score excluded from competition cannot receive a canonical intensity.
    auto excluded = sample;
    for (auto& window : excluded.windows)
      window.score = CrossFitScorer::NOT_KEPT;
    writeQpxDataset((output / "excluded").string(), lib, {excluded}, quant, {}, {});
    auto excluded_features = read(output / "excluded/opendia.feature.parquet");
    require(column<arrow::StringArray>(excluded_features, "x_opendia_score_kept")->GetString(0) == "false", "diagnostic score mask explicit");
    auto score_lists = column<arrow::ListArray>(excluded_features, "additional_scores");
    auto score_structs = std::static_pointer_cast<arrow::StructArray>(score_lists->values());
    auto score_names = std::static_pointer_cast<arrow::StringArray>(score_structs->field(0));
    for (std::int64_t i = 0; i < score_names->length(); ++i)
      require(score_names->GetString(i) != "opendia_score" && score_names->GetString(i) != "opendia_qvalue",
              "masked sentinel is not a canonical margin or q-value");
    require(intensities(excluded_features)->null_count() == 3, "unscored candidates have no canonical intensity");

    // A row group boundary must not change values or drop rows.
    auto large = sample;
    large.windows.clear();
    large.features.clear();
    for (std::uint32_t i = 0; i < 8201; ++i)
    {
      auto w = sample.windows[0];
      w.candidate = i;
      large.windows.push_back(w);
      large.features.insert(large.features.end(), {0.25f, -0.5f});
    }
    writeCandidateEvidence((output / "large").string(), lib, large);
    auto parquet_reader = parquet::ParquetFileReader::OpenFile((output / "large/opendia/acquisition.candidates.parquet").string());
    require(parquet_reader->metadata()->num_rows() == 8201 && parquet_reader->metadata()->num_row_groups() == 3, "bounded row groups");
    const auto original_charge = lib.precursors[0].charge;
    for (const int unsupported_charge : {0, 32768})
    {
      lib.precursors[0].charge = unsupported_charge;
      bool charge_rejected = false;
      try { writeQpxDataset(output.string(), lib, {sample}, quant, {}, {}); }
      catch (const std::runtime_error&) { charge_rejected = true; }
      require(charge_rejected, "unsupported charge rejected before canonical int16 narrowing");
      require(read(output / "opendia.feature.parquet")->num_rows() == 3, "invalid charge preserves previous complete view");
    }
    lib.precursors[0].charge = original_charge;

    auto duplicate = sample;
    duplicate.windows[1].candidate = duplicate.windows[0].candidate;
    const auto duplicate_directory = output / "duplicate";
    writeCandidateEvidence(duplicate_directory.string(), lib, sample);
    bool rejected = false;
    try
    {
      writeCandidateEvidence(duplicate_directory.string(), lib, duplicate);
    }
    catch (const std::runtime_error&)
    {
      rejected = true;
    }
    require(rejected, "duplicate identity rejected rather than dropping row");
    require(read(duplicate_directory / "opendia/acquisition.candidates.parquet")->num_rows() == 3,
            "failed replacement preserves previous complete output");
    for (const auto& file : fs::recursive_directory_iterator(output))
      require(file.path().filename().string().find(".partial.") == std::string::npos, "failed export leaves no partial file");
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
