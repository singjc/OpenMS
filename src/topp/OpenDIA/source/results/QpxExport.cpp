// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/QPXDatasetExport.h>
#include <OpenMS/FORMAT/QPXIdentity.h>
#include <algorithm>
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
#include <parquet/arrow/writer.h>
#include <results/QpxExport.h>
#include <set>
#include <stdexcept>
#include <unordered_set>

namespace OpenMS::OpenDIA
{
namespace
{
  constexpr std::size_t BATCH_SIZE = 4096;
  constexpr const char* IDENTITY = "run_file_name,opendia_precursor_id,opendia_candidate_row";
  using Scalar = std::shared_ptr<arrow::Scalar>;
  using Row = std::map<std::string, Scalar>;

  void check(const arrow::Status& status)
  {
    if (! status.ok()) throw std::runtime_error("OpenDIA QPX: " + status.ToString());
  }
  template<typename T>
  T unwrap(arrow::Result<T> value)
  {
    if (! value.ok()) throw std::runtime_error("OpenDIA QPX: " + value.status().ToString());
    return std::move(value).ValueOrDie();
  }
  Scalar text(const std::string& value)
  { return std::make_shared<arrow::StringScalar>(value); }
  Scalar real(double value)
  { return std::isfinite(value) ? std::make_shared<arrow::DoubleScalar>(value) : arrow::MakeNullScalar(arrow::float64()); }
  Scalar single(double value)
  {
    const float narrowed = static_cast<float>(value);
    if (std::isfinite(value) && !std::isfinite(narrowed))
      throw std::runtime_error("OpenDIA QPX: finite numeric value exceeds the float32 range");
    return std::make_shared<arrow::FloatScalar>(narrowed);
  }
  Scalar boolean(bool value)
  { return std::make_shared<arrow::BooleanScalar>(value); }
  Scalar integer(std::int64_t value)
  { return std::make_shared<arrow::Int64Scalar>(value); }

  Scalar canonicalExtension(const Scalar& value)
  {
    if (! value->is_valid) return arrow::MakeNullScalar(arrow::utf8());
    if (value->type->id() == arrow::Type::DOUBLE) return text(QPXDatasetExport::formatNumber(static_cast<const arrow::DoubleScalar&>(*value).value));
    if (value->type->id() == arrow::Type::FLOAT) return text(QPXDatasetExport::formatNumber(static_cast<const arrow::FloatScalar&>(*value).value));
    return text(value->ToString());
  }
  Row canonicalRow(const Row& source, bool feature)
  {
    Row row;
    for (const auto& [name, value] : source)
    {
      if (name.starts_with("opendia_")) row.emplace("x_" + name, canonicalExtension(value));
      else if (feature && (name == "peptide_global_qvalue" || name == "pg_qvalue"))
        row.emplace("x_opendia_" + name, canonicalExtension(value));
      else
        row.emplace(name, value);
    }
    return row;
  }
  std::shared_ptr<arrow::DataType> featureIntensitiesType()
  { return arrow::list(arrow::struct_({arrow::field("label", arrow::utf8(), false), arrow::field("intensity", arrow::float32(), true)})); }
  Scalar list(const std::shared_ptr<arrow::DataType>& type, const std::vector<Scalar>& values)
  {
    auto builder = unwrap(arrow::MakeBuilder(type->field(0)->type(), arrow::default_memory_pool()));
    for (const auto& value : values)
      check(builder->AppendScalar(*value));
    return std::make_shared<arrow::ListScalar>(unwrap(builder->Finish()), type);
  }
  Scalar strings(const std::vector<std::string>& values)
  {
    std::vector<Scalar> items;
    for (const auto& value : values)
      items.push_back(text(value));
    return list(arrow::list(arrow::utf8()), items);
  }
  Scalar structure(const std::shared_ptr<arrow::DataType>& type, const std::vector<Scalar>& values)
  { return std::make_shared<arrow::StructScalar>(values, type); }
  Scalar scores(const RunWindow& window, bool winner)
  {
    auto type = QPXFeatureSchema::additionalScoresType();
    auto item = type->field(0)->type();
    std::vector<Scalar> values;
    const auto add = [&](const std::string& name, double value, bool higher) {
      if (std::isfinite(value)) values.push_back(structure(item, {text(name), real(value), boolean(higher)}));
    };
    if (window.score > CrossFitScorer::NOT_KEPT) add("opendia_score", window.score, true);
    if (window.score > CrossFitScorer::NOT_KEPT) add("opendia_qvalue", window.q_value, false);
    if (window.network_score > CrossFitScorer::NOT_KEPT) add("opendia_network_score", window.network_score, true);
    if (window.network_score > CrossFitScorer::NOT_KEPT) add("opendia_network_qvalue", window.network_q_value, false);
    add("opendia_winner", winner ? 1.0 : 0.0, true);
    add("opendia_score_kept", std::isfinite(window.score) && window.score > CrossFitScorer::NOT_KEPT ? 1.0 : 0.0, true);
    return list(type, values);
  }

  std::string runName(const std::string& run)
  { return ArrowIOHelpers::qpxRunFileName(run); }
  std::vector<std::size_t> winners(const AnalysisLibrary& library, const RunResult& run)
  {
    std::vector<std::size_t> result(library.precursors.size(), std::numeric_limits<std::size_t>::max());
    for (std::size_t i = 0; i < run.windows.size(); ++i)
    {
      const auto& w = run.windows[i];
      if (w.precursor >= result.size()) throw std::runtime_error("OpenDIA QPX: precursor index out of range");
      if (! std::isfinite(w.score) || w.score <= CrossFitScorer::NOT_KEPT) continue;
      auto& best = result[w.precursor];
      if (best == std::numeric_limits<std::size_t>::max() || w.score > run.windows[best].score) best = i;
    }
    return result;
  }

  std::map<std::string, std::string> metadata(const ReportingOptions& options, const std::string& view)
  {
    return {
      {"qpx_view", view},
      {"opendia_schema_version", "1"},
      {"opendia_algorithm", "OpenDIA precursor scoring; pooled paired precursor competition; PyProphet-adapted level-context inference"},
      {"opendia_peakgroup_fdr", QPXDatasetExport::formatNumber(options.peakgroup_fdr)},
      {"opendia_peptide_fdr", QPXDatasetExport::formatNumber(options.peptide_fdr)},
      {"opendia_protein_fdr", QPXDatasetExport::formatNumber(options.protein_fdr)},
      {"opendia_rt_unit", "second"},
      {"opendia_library_rt_unit", "source_library_scale"},
      {"opendia_masked_score_sentinel", QPXDatasetExport::formatNumber(CrossFitScorer::NOT_KEPT)},
      {"opendia_openms_version", VersionInfo::getVersion()},
      {"opendia_openms_revision", VersionInfo::getRevision()},
      {"opendia_canonical_extension_encoding",
       "nullable x_opendia_* strings; shortest round-trip decimals preserve source binary64; typed candidate companion retains scorer binary32"},
      {"opendia_entity_estimator_status", "shared LevelContextInference; empirical calibration remains dataset-dependent"},
      {"opendia_entity_inference_reference", "10.1038/nmeth.4398"},
      {"opendia_entity_inference_context", options.inference_context},
      {"opendia_entity_inference_contexts", "run-specific,experiment-wide,global"},
      {"opendia_peptide_identity", "canonical modified sequence and target/decoy class; charge collapsed"},
      {"opendia_protein_identity", "supplied accession and target/decoy class; only uniquely mapped modified peptides"},
      {"opendia_entity_error_model", "empirical p-values; bootstrap pi0; lambda=0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45; FDR; probit local FDR adj=1.5 eps=1e-8"},
      {"opendia_area_semantics", "raw MS2 precursor area; not LFQ"}};
  }

  // Only one row group is resident. We deliberately do not use the older native value validator:
  // current QPX permits null observed m/z and a producer-declared candidate identity composite.
  class Writer
  {
  public:
    Writer(const std::filesystem::path& path,
           std::shared_ptr<arrow::Schema> schema,
           const std::string& view,
           const std::map<std::string, std::string>& extra):
        path_(path),
        schema_(std::move(schema))
    {
      static std::atomic<std::uint64_t> serial {0};
      temporary_ = path_;
      temporary_
        += ".partial." + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + "." + std::to_string(serial.fetch_add(1));
      std::filesystem::create_directories(path.parent_path());
      auto footer = ArrowIOHelpers::qpxFileMetadata(view + "_file", {}, extra);
      if (! footer) throw std::runtime_error("OpenDIA QPX: cannot create metadata");
      auto mutable_footer = footer->Copy();
      if (view == "feature")
      {
        check(mutable_footer->Set("identity_source", "supplied"));
        check(mutable_footer->Set("identity_composite", "run_file_name,cv_params"));
        check(mutable_footer->Set("opendia_identity_composite", IDENTITY));
      }
      else if (view == "candidates")
      {
        check(mutable_footer->Set("identity_composite", IDENTITY));
        check(mutable_footer->Set("primary_key", "feature_id"));
      }
      else if (view == "pg")
        check(mutable_footer->Set("identity_composite", "pg_accessions,grouped_runs,label"));
      schema_ = schema_->WithMetadata(mutable_footer);
      try
      {
        output_ = unwrap(arrow::io::FileOutputStream::Open(temporary_.string()));
        parquet::WriterProperties::Builder properties;
        properties.compression(parquet::Compression::ZSTD);
        parquet::ArrowWriterProperties::Builder arrow_properties;
        arrow_properties.store_schema();
        writer_
          = unwrap(parquet::arrow::FileWriter::Open(*schema_, arrow::default_memory_pool(), output_, properties.build(), arrow_properties.build()));
        reset();
      }
      catch (...)
      {
        writer_.reset();
        if (output_) (void)output_->Close();
        std::error_code error;
        std::filesystem::remove(temporary_, error);
        throw;
      }
    }
    ~Writer()
    {
      if (! closed_)
      {
        writer_.reset();
        if (output_) (void)output_->Close();
        std::error_code error;
        std::filesystem::remove(temporary_, error);
      }
    }
    void append(const Row& row)
    {
      for (int i = 0; i < schema_->num_fields(); ++i)
      {
        const auto it = row.find(schema_->field(i)->name());
        if (it == row.end()) check(builders_[i]->AppendNull());
        else
          check(builders_[i]->AppendScalar(*it->second));
      }
      if (++rows_ == BATCH_SIZE) flush();
    }
    void close()
    {
      // An empty dataset must still contain a typed, readable Parquet file.
      if (rows_ || ! written_) flush();
      check(writer_->Close());
      check(output_->Close());
      std::filesystem::rename(temporary_, path_);
      closed_ = true;
    }
    void appendTable(const std::shared_ptr<arrow::Table>& table)
    {
      if (! table) throw std::runtime_error("OpenDIA QPX: null native view");
      for (std::int64_t i = 0; i < table->num_rows(); ++i)
      {
        Row row;
        for (int j = 0; j < table->num_columns(); ++j)
          row.emplace(table->field(j)->name(), unwrap(table->column(j)->GetScalar(i)));
        append(row);
      }
    }

  private:
    void reset()
    {
      builders_.clear();
      for (const auto& field : schema_->fields())
        builders_.push_back(unwrap(arrow::MakeBuilder(field->type(), arrow::default_memory_pool())));
    }
    void flush()
    {
      std::vector<std::shared_ptr<arrow::Array>> arrays;
      for (auto& builder : builders_)
        arrays.push_back(unwrap(builder->Finish()));
      auto table = arrow::Table::Make(schema_, arrays, static_cast<std::int64_t>(rows_));
      check(table->ValidateFull());
      check(writer_->WriteTable(*table, BATCH_SIZE));
      written_ = true;
      rows_ = 0;
      reset();
    }
    std::filesystem::path path_;
    std::filesystem::path temporary_;
    std::shared_ptr<arrow::Schema> schema_;
    std::shared_ptr<arrow::io::FileOutputStream> output_;
    std::unique_ptr<parquet::arrow::FileWriter> writer_;
    std::vector<std::unique_ptr<arrow::ArrayBuilder>> builders_;
    std::size_t rows_ {};
    bool written_ {false};
    bool closed_ {false};
  };

  std::vector<std::shared_ptr<arrow::Field>> extensions()
  {
    return {arrow::field("opendia_precursor_id", arrow::utf8(), false),
            arrow::field("opendia_candidate_row", arrow::int64(), false),
            arrow::field("opendia_rt_apex", arrow::float64(), false),
            arrow::field("opendia_rt_start", arrow::float64(), false),
            arrow::field("opendia_rt_end", arrow::float64(), false),
            arrow::field("opendia_area", arrow::float64(), false),
            arrow::field("opendia_mz", arrow::float64(), false),
            arrow::field("opendia_library_rt", arrow::float64()),
            arrow::field("opendia_score_kept", arrow::boolean(), false),
            arrow::field("opendia_network_score_kept", arrow::boolean(), false),
            arrow::field("opendia_score", arrow::float32(), false),
            arrow::field("opendia_qvalue", arrow::float64(), false),
            arrow::field("opendia_network_score", arrow::float32(), false),
            arrow::field("opendia_network_qvalue", arrow::float64(), false),
            arrow::field("opendia_winner", arrow::boolean(), false)};
  }
  Row candidateRow(const AnalysisLibrary::Precursor& p, const RunResult& run, const RunWindow& w, bool winner, std::unordered_set<std::int64_t>& ids)
  {
    const auto id = outputFeatureId(run.run, std::string(p.id), w.candidate);
    if (! ids.insert(id).second) throw std::runtime_error("OpenDIA QPX: duplicate/colliding candidate feature_id " + std::to_string(id));
    return {{"feature_id", integer(id)},
            {"run_file_name", text(runName(run.run))},
            {"opendia_precursor_id", text(std::string(p.id))},
            {"opendia_candidate_row", integer(w.candidate)},
            {"opendia_rt_apex", real(w.peak.rt_apex)},
            {"opendia_rt_start", real(w.peak.rt_start)},
            {"opendia_rt_end", real(w.peak.rt_end)},
            {"opendia_area", real(w.peak.area)},
            {"opendia_mz", real(p.mz)},
            {"opendia_library_rt", p.has_library_rt ? real(p.library_rt) : arrow::MakeNullScalar(arrow::float64())},
            {"opendia_score_kept", boolean(std::isfinite(w.score) && w.score > CrossFitScorer::NOT_KEPT)},
            {"opendia_network_score_kept", boolean(std::isfinite(w.network_score) && w.network_score > CrossFitScorer::NOT_KEPT)},
            {"opendia_score", single(w.score)},
            {"opendia_qvalue", real(w.q_value)},
            {"opendia_network_score", single(w.network_score)},
            {"opendia_network_qvalue", real(w.network_q_value)},
            {"opendia_winner", boolean(winner)},
            {"is_decoy", boolean(p.decoy)}};
  }
  std::string decoyPrefix(const Quantification& quant)
  {
    std::set<std::string> targets, reserved;
    for (const auto& group : quant.protein_groups)
    {
      reserved.insert(group.proteins.begin(), group.proteins.end());
      if (! group.decoy) targets.insert(group.proteins.begin(), group.proteins.end());
    }
    std::string prefix = "DECOY_";
    for (;;)
    {
      bool collision = false;
      for (const auto& group : quant.protein_groups)
        if (group.decoy)
          for (const auto& protein : group.proteins)
            if (targets.contains(protein) && reserved.contains(prefix + protein)) collision = true;
      if (! collision) return prefix;
      prefix += "DECOY_";
    }
  }
  std::vector<std::string> groupAccessions(const QuantifiedEntity& group, const std::string& prefix,
                                         const std::set<std::string>& targets)
  {
    auto proteins = group.proteins;
    if (group.decoy)
      for (auto& protein : proteins)
        if (targets.contains(protein)) protein = prefix + protein;
    return proteins;
  }
  std::vector<std::string> geneNames(const std::string& raw)
  {
    std::vector<std::string> names;
    std::size_t begin = 0;
    while (begin <= raw.size())
    {
      const auto end = raw.find(';', begin);
      const auto token = raw.substr(begin, end == std::string::npos ? end : end - begin);
      const auto first = token.find_first_not_of(" \t\r\n");
      if (first != std::string::npos)
      {
        const auto name = token.substr(first, token.find_last_not_of(" \t\r\n") - first + 1);
        if (std::find(names.begin(), names.end(), name) == names.end()) names.push_back(name);
      }
      if (end == std::string::npos) break;
      begin = end + 1;
    }
    return names;
  }
  const EntityRun* entityRun(const QuantifiedEntity& entity, std::size_t run)
  {
    for (const auto& evidence : entity.runs)
      if (evidence.run == run) return &evidence;
    return nullptr;
  }
} // namespace

std::int64_t outputFeatureId(const std::string& run, const std::string& precursor, std::uint32_t candidate)
{
  return QPXIdentity::deriveId(
    {QPXIdentity::Value {runName(run)}, QPXIdentity::Value {precursor}, QPXIdentity::Value {static_cast<Int64>(candidate)}});
}

void writeCandidateEvidence(const std::string& directory, const AnalysisLibrary& library, const RunResult& run)
{
  if (run.features.size() != run.windows.size() * run.feature_names.size())
    throw std::runtime_error("OpenDIA QPX: scorer feature matrix shape mismatch");
  auto fields = extensions();
  fields.insert(fields.begin(), {arrow::field("feature_id", arrow::int64(), false), arrow::field("run_file_name", arrow::utf8(), false),
                                 arrow::field("is_decoy", arrow::boolean(), false)});
  std::set<std::string> names;
  for (const auto& field : fields)
    names.insert(field->name());
  for (const auto& name : run.feature_names)
  {
    if (! names.insert(name).second) throw std::runtime_error("OpenDIA QPX: duplicate/reserved scorer feature name " + name);
    fields.push_back(arrow::field(name, arrow::float32(), false));
  }
  auto meta = metadata({}, "opendia_candidates");
  meta["opendia_calibration_json"] = run.calibration_json.empty() ? "null" : run.calibration_json;
  Writer writer(std::filesystem::path(directory) / "opendia" / (runName(run.run) + ".candidates.parquet"), arrow::schema(fields), "candidates", meta);
  const auto best = winners(library, run);
  std::unordered_set<std::int64_t> ids;
  for (std::size_t i = 0; i < run.windows.size(); ++i)
  {
    const auto& w = run.windows[i];
    auto row = candidateRow(library.precursors.at(w.precursor), run, w, best[w.precursor] == i, ids);
    for (std::size_t j = 0; j < run.feature_names.size(); ++j)
      row.emplace(run.feature_names[j], single(run.features[i * run.feature_names.size() + j]));
    writer.append(row);
  }
  writer.close();
}

void writeQpxDataset(const std::string& directory,
                     const AnalysisLibrary& library,
                     const std::vector<RunResult>& runs,
                     const Quantification& quant,
                     const GeneAnnotations& genes,
                     const ReportingOptions& options)
{
  std::set<std::string> target_accessions;
  for (const auto& group : quant.protein_groups)
    if (!group.decoy) target_accessions.insert(group.proteins.begin(), group.proteins.end());
  const auto datasetMetadata = [&](const std::string& view) {
    auto result = metadata(options, view);
    std::string warnings;
    for (const auto& warning : quant.inference_warnings)
    {
      if (!warnings.empty()) warnings += "\n";
      warnings += warning;
    }
    result["opendia_entity_inference_warnings"] = warnings;
    result["opendia_decoy_accession_policy"] = "source accessions preserved; target/decoy collisions only are prefixed";
    return result;
  };
  auto fields = QPXFeatureSchema::schema()->fields();
  for (auto& field : fields)
  {
    if (field->name() == "observed_mz" || field->name() == "calculated_mz") field = field->WithNullable(true);
    else if (field->name() == "intensities")
      field = arrow::field("intensities", featureIntensitiesType());
  }
  for (const auto& field : extensions())
    fields.push_back(arrow::field("x_" + field->name(), arrow::utf8()));
  for (const auto& name : {"global_qvalue", "precursor_qvalue", "peptide_qvalue"})
    fields.push_back(arrow::field(name, arrow::float64()));
  fields.push_back(arrow::field("x_opendia_peptide_global_qvalue", arrow::utf8()));
  fields.push_back(arrow::field("x_opendia_pg_qvalue", arrow::utf8()));
  fields.push_back(arrow::field("x_opendia_peptide_experiment_qvalue", arrow::utf8()));
  fields.push_back(arrow::field("x_opendia_pg_experiment_qvalue", arrow::utf8()));
  fields.push_back(arrow::field("pg_ids", arrow::list(arrow::int64())));
  const auto decoy_prefix = decoyPrefix(quant);
  auto feature_metadata = datasetMetadata("feature");
  feature_metadata["opendia_decoy_accession_prefix"] = decoy_prefix;
  Writer feature(std::filesystem::path(directory) / "opendia.feature.parquet", arrow::schema(fields), "feature", feature_metadata);
  std::unordered_set<std::int64_t> ids;
  std::set<std::string> run_names;
  for (std::size_t r = 0; r < runs.size(); ++r)
  {
    const auto& run = runs[r];
    if (! run_names.insert(runName(run.run)).second) throw std::runtime_error("OpenDIA QPX: duplicate run_file_name " + runName(run.run));
    const auto best = winners(library, run);
    for (std::size_t i = 0; i < run.windows.size(); ++i)
    {
      const auto& w = run.windows[i];
      const auto& p = library.precursors.at(w.precursor);
      const bool winner = best[w.precursor] == i;
      auto row = candidateRow(p, run, w, winner, ids);
      const auto sequence = AASequence::fromString(p.sequence);
      row["sequence"] = text(sequence.toUnmodifiedString());
      row["peptidoform"] = text(ProForma::toString(ProForma::fromAASequence(sequence), ProForma::WriteMode::CANONICAL));
      if (p.charge <= 0 || p.charge > std::numeric_limits<std::int16_t>::max())
        throw std::runtime_error("OpenDIA QPX: precursor charge is outside the positive int16 range");
      row["charge"] = std::make_shared<arrow::Int16Scalar>(static_cast<std::int16_t>(p.charge));
      if (p.charge > 0) row["calculated_mz"] = single(sequence.getMZ(p.charge));
      row["scan"] = list(arrow::list(arrow::int32()), {});
      auto cv_type = QPXFeatureSchema::cvParamsType();
      row["cv_params"] = list(cv_type, {structure(cv_type->field(0)->type(), {text("opendia_precursor_id"), text(std::string(p.id))}),
                                        structure(cv_type->field(0)->type(), {text("opendia_candidate_row"), text(std::to_string(w.candidate))})});
      row["rt"] = single(w.peak.rt_apex);
      row["rt_start"] = single(w.peak.rt_start);
      row["rt_stop"] = single(w.peak.rt_end);
      // The source library coordinate may be iRT. No calibrated predicted RT is available here.
      row["additional_scores"] = scores(w, winner);
      // Keep the known LFQ channel even without a quantity. Current QPX permits null nested intensity.
      auto intensity_type = featureIntensitiesType();
      row["intensities"]
        = list(intensity_type,
               {structure(intensity_type->field(0)->type(), {text("LFQ"), winner ? single(w.peak.area) : arrow::MakeNullScalar(arrow::float32())})});
      if (std::isfinite(w.score) && w.score > CrossFitScorer::NOT_KEPT) row["precursor_qvalue"] = real(w.q_value);
      if (w.precursor < quant.global_precursor_q.size()) row["global_qvalue"] = real(quant.global_precursor_q[w.precursor]);
      if (w.precursor < quant.precursor_peptide.size() && quant.precursor_peptide[w.precursor] < quant.peptides.size())
      {
        const auto& peptide = quant.peptides[quant.precursor_peptide[w.precursor]];
        row["peptide_global_qvalue"] = real(peptide.global_q_value);
        if (const auto* evidence = entityRun(peptide, r))
        {
          row["peptide_qvalue"] = real(evidence->q_value);
          row["opendia_peptide_experiment_qvalue"] = real(evidence->experiment_q_value);
        }
      }
      if (w.precursor < quant.precursor_protein_group.size() && quant.precursor_protein_group[w.precursor] < quant.protein_groups.size())
      {
        const auto& group = quant.protein_groups[quant.precursor_protein_group[w.precursor]];
        std::vector<Scalar> accessions;
        auto type = QPXFeatureSchema::pgAccessionsType();
        auto item = type->field(0)->type();
        const auto proteins = groupAccessions(group, decoy_prefix, target_accessions);
        for (const auto& protein : proteins)
          accessions.push_back(structure(item, {text(protein), arrow::MakeNullScalar(arrow::int32()), arrow::MakeNullScalar(arrow::int32()),
                                                arrow::MakeNullScalar(arrow::utf8()), arrow::MakeNullScalar(arrow::utf8())}));
        row["pg_accessions"] = list(type, accessions);
        if (! proteins.empty()) row["anchor_protein"] = text(proteins.front());
        row["pg_global_qvalue"] = real(group.global_q_value);
        if (const auto* evidence = entityRun(group, r))
        {
          row["pg_qvalue"] = real(evidence->q_value);
          row["opendia_pg_experiment_qvalue"] = real(evidence->experiment_q_value);
          const auto pg_id = QPXIdentity::pgId(proteins, {runName(run.run)}, std::string("LFQ"));
          row["pg_ids"] = list(arrow::list(arrow::int64()), {integer(pg_id)});
        }
        row["unique"] = boolean(group.proteins.size() == 1);
      }
      const auto gene = genes.find(std::string(p.id));
      if (gene != genes.end() && ! gene->second.empty()) row["gg_names"] = strings(geneNames(gene->second));
      feature.append(canonicalRow(row, true));
    }
  }
  feature.close();

  auto pg_fields = QPXPgSchema::schema()->fields();
  pg_fields.push_back(arrow::field("x_opendia_decoy_class", arrow::utf8()));
  pg_fields.push_back(arrow::field("x_opendia_area", arrow::utf8()));
  pg_fields.push_back(arrow::field("x_opendia_score", arrow::utf8()));
  pg_fields.push_back(arrow::field("x_opendia_experiment_qvalue", arrow::utf8()));
  auto pg_metadata = datasetMetadata("pg");
  pg_metadata["opendia_decoy_accession_prefix"] = decoy_prefix;
  Writer pg(std::filesystem::path(directory) / "opendia.pg.parquet", arrow::schema(pg_fields), "pg", pg_metadata);
  ids.clear();
  for (const auto& group : quant.protein_groups)
  {
    if (group.proteins.empty()) throw std::runtime_error("OpenDIA QPX: resolved protein group has no accession");
    for (const auto& evidence : group.runs)
    {
      if (evidence.run >= runs.size()) throw std::runtime_error("OpenDIA QPX: protein group run index out of range");
      const std::vector<std::string> grouped_runs {runName(runs[evidence.run].run)};
      const std::string decoy_class = group.decoy ? "decoy" : "target";
      const auto proteins = groupAccessions(group, decoy_prefix, target_accessions);
      const auto id = QPXIdentity::pgId(proteins, grouped_runs, std::string("LFQ"));
      if (! ids.insert(id).second) throw std::runtime_error("OpenDIA QPX: duplicate/colliding pg_id");
      Row row {{"pg_id", integer(id)},
               {"opendia_decoy_class", text(decoy_class)},
               {"pg_accessions", strings(proteins)},
               {"anchor_protein", text(proteins.front())},
               {"grouped_runs", strings(grouped_runs)},
               {"global_qvalue", real(group.global_q_value)},
               {"pg_qvalue", real(evidence.q_value)},
               {"label", text("LFQ")},
               {"is_decoy", boolean(group.decoy)},
               {"opendia_score", real(evidence.score)},
               {"opendia_experiment_qvalue", real(evidence.experiment_q_value)},
               {"peptides", list(QPXPgSchema::peptidesType(), {})}};
      if (evidence.has_area)
      {
        row["intensity"] = single(evidence.area);
        row["opendia_area"] = real(evidence.area);
      }
      pg.append(canonicalRow(row, false));
    }
  }
  pg.close();

  // Producer companion keeps every context's full-precision statistics without duplicating
  // them across millions of candidate features. A null run denotes the global context.
  Writer inference(std::filesystem::path(directory) / "opendia" / "inference.parquet",
                   arrow::schema({arrow::field("level", arrow::utf8(), false),
                                  arrow::field("entity_id", arrow::utf8(), false),
                                  arrow::field("is_decoy", arrow::boolean(), false),
                                  arrow::field("run_file_name", arrow::utf8()),
                                  arrow::field("context", arrow::utf8(), false),
                                  arrow::field("score", arrow::float64()),
                                  arrow::field("p_value", arrow::float64()),
                                  arrow::field("q_value", arrow::float64()),
                                  arrow::field("pep", arrow::float64())}),
                   "inference", datasetMetadata("opendia_inference"));
  const auto appendInference = [&](const std::vector<QuantifiedEntity>& entities, const std::string& level) {
    for (const auto& entity : entities)
    {
      double best = -std::numeric_limits<double>::infinity();
      for (const auto& evidence : entity.runs)
      {
        if (evidence.run >= runs.size()) throw std::runtime_error("OpenDIA QPX: entity run index out of range");
        best = std::max(best, evidence.score);
        const Row base{{"level", text(level)}, {"entity_id", text(entity.id)}, {"is_decoy", boolean(entity.decoy)},
                       {"run_file_name", text(runName(runs[evidence.run].run))}, {"score", real(evidence.score)}};
        auto row = base;
        row["context"] = text("run-specific");
        row["p_value"] = real(evidence.p_value);
        row["q_value"] = real(evidence.q_value);
        row["pep"] = real(evidence.pep);
        inference.append(row);
        row = base;
        row["context"] = text("experiment-wide");
        row["p_value"] = real(evidence.experiment_p_value);
        row["q_value"] = real(evidence.experiment_q_value);
        row["pep"] = real(evidence.experiment_pep);
        inference.append(row);
      }
      if (!entity.runs.empty())
        inference.append({{"level", text(level)}, {"entity_id", text(entity.id)}, {"is_decoy", boolean(entity.decoy)},
                          {"context", text("global")}, {"score", real(best)}, {"p_value", real(entity.global_p_value)},
                          {"q_value", real(entity.global_q_value)}, {"pep", real(entity.global_pep)}});
    }
  };
  appendInference(quant.peptides, "peptide");
  appendInference(quant.protein_groups, "protein");
  inference.close();

  std::vector<QPXDatasetExport::Run> run_rows;
  std::vector<QPXDatasetExport::Sample> sample_rows;
  std::vector<Scalar> calibration;
  for (const auto& run : runs)
  {
    QPXDatasetExport::Run row;
    row.accession = runName(run.run);
    row.run_file_name = runName(run.run);
    row.file_name = std::filesystem::path(run.run).filename().string();
    row.acquisition_method = "DIA";
    const std::string sample_accession = "opendia_run:" + runName(run.run);
    QPXDatasetExport::RunSample sample;
    sample.accession = sample_accession;
    sample.label = "LFQ";
    row.samples.push_back(std::move(sample));
    QPXDatasetExport::Sample technical_sample;
    technical_sample.accession = sample_accession;
    sample_rows.push_back(std::move(technical_sample));
    run_rows.push_back(std::move(row));
    calibration.push_back(run.calibration_json.empty() ? arrow::MakeNullScalar(arrow::utf8()) : text(run.calibration_json));
  }
  auto run_table = QPXDatasetExport::runView(run_rows);
  if (! run_table) throw std::runtime_error("OpenDIA QPX: cannot create run view");
  auto calibration_builder = unwrap(arrow::MakeBuilder(arrow::utf8(), arrow::default_memory_pool()));
  for (const auto& value : calibration)
    check(calibration_builder->AppendScalar(*value));
  run_table = unwrap(run_table->AddColumn(run_table->num_columns(), arrow::field("opendia_calibration_json", arrow::utf8()),
                                          std::make_shared<arrow::ChunkedArray>(unwrap(calibration_builder->Finish()))));
  QPXDatasetExport::ProvenanceStep step;
  step.category = "quantification";
  step.name = "opendia";
  step.tool = "OpenDIA";
  step.output_views = {"feature", "pg", "run", "sample"};
  for (const auto& [key, value] : datasetMetadata("provenance"))
    step.parameters.emplace_back(key, value);
  auto run_metadata = datasetMetadata("run");
  run_metadata["opendia_sample_accession_source"] = "producer-generated technical run identifier; biological sample unknown";
  Writer run_writer(std::filesystem::path(directory) / "opendia.run.parquet", run_table->schema(), "run", run_metadata);
  run_writer.appendTable(run_table);
  run_writer.close();
  auto sample_table = QPXDatasetExport::sampleView(sample_rows);
  if (! sample_table) throw std::runtime_error("OpenDIA QPX: cannot create sample view");
  auto sample_metadata = datasetMetadata("sample");
  sample_metadata["opendia_sample_accession_source"] = "producer-generated technical run identifier; biological sample unknown";
  Writer sample_writer(std::filesystem::path(directory) / "opendia.sample.parquet", sample_table->schema(), "sample", sample_metadata);
  sample_writer.appendTable(sample_table);
  sample_writer.close();
  auto provenance_table = QPXDatasetExport::provenanceView({step});
  if (! provenance_table) throw std::runtime_error("OpenDIA QPX: cannot create provenance view");
  Writer provenance_writer(std::filesystem::path(directory) / "opendia.provenance.parquet", provenance_table->schema(), "provenance",
                           datasetMetadata("provenance"));
  provenance_writer.appendTable(provenance_table);
  provenance_writer.close();
}
} // namespace OpenMS::OpenDIA
