// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $

#include <results/ReportExport.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <stdexcept>
#include <unordered_set>
#include <variant>

namespace OpenMS::OpenDIA
{
  namespace
  {
    namespace fs = std::filesystem;
    using Cell = std::variant<std::monostate, std::string, double, std::int64_t, bool>;
    using Row = std::vector<Cell>;
    const double na = std::numeric_limits<double>::quiet_NaN();

    void check(const arrow::Status& status)
    {
      if (!status.ok()) throw std::runtime_error("OpenDIA report: " + status.ToString());
    }

    template <typename T> T unwrap(arrow::Result<T> result)
    {
      check(result.status());
      return std::move(result).ValueOrDie();
    }

    Cell number(double value)
    {
      return std::isfinite(value) ? Cell(value) : Cell{};
    }

    std::string escaped(const std::string& text)
    {
      if (text != "NA" && text.find_first_of("\t\r\n\"") == std::string::npos) return text;
      std::string result = "\"";
      for (const char c : text)
      {
        result += c;
        if (c == '"') result += '"';
      }
      return result + '"';
    }

    // Files become visible only after all writers have closed successfully. Every failure removes
    // its unique partial files, including failures during Arrow initialization or final close.
    struct Publication
    {
      std::vector<std::pair<fs::path, fs::path>> paths;
      ~Publication()
      {
        for (const auto& path : paths)
        {
          std::error_code error;
          fs::remove(path.first, error);
        }
      }
      fs::path add(const fs::path& destination)
      {
        static std::atomic<unsigned long long> serial{0};
        auto partial = destination;
        partial += ".partial." + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())
                 + "." + std::to_string(serial++);
        paths.emplace_back(partial, destination);
        return partial;
      }
      void publish()
      {
        for (const auto& path : paths) fs::rename(path.first, path.second);
      }
    };

    class Table
    {
    public:
      Table(Publication& publication, const fs::path& base, std::vector<std::shared_ptr<arrow::Field>> fields,
            bool tsv, bool parquet) : schema_(arrow::schema(std::move(fields)))
      {
        // A wide matrix still has bounded buffering, not a run-count times 4096 allocation floor.
        batch_limit_ = std::max<std::size_t>(1, std::min<std::size_t>(4096, 1048576 / schema_->num_fields()));
        if (tsv)
        {
          text_.exceptions(std::ios::badbit | std::ios::failbit);
          text_.open(publication.add(base.string() + ".tsv"), std::ios::binary);
          text_.imbue(std::locale::classic());
          text_ << std::setprecision(std::numeric_limits<double>::max_digits10);
          for (int i = 0; i < schema_->num_fields(); ++i)
          {
            if (i) text_ << '\t';
            text_ << escaped(schema_->field(i)->name());
          }
          text_ << '\n';
        }
        if (parquet)
        {
          output_ = unwrap(arrow::io::FileOutputStream::Open(publication.add(base.string() + ".parquet").string()));
          parquet::WriterProperties::Builder properties;
          properties.compression(parquet::Compression::ZSTD);
          writer_ = unwrap(parquet::arrow::FileWriter::Open(*schema_, arrow::default_memory_pool(), output_, properties.build()));
          for (const auto& field : schema_->fields()) builders_.push_back(unwrap(arrow::MakeBuilder(field->type())));
        }
      }

      void append(const Row& row)
      {
        if (row.size() != static_cast<std::size_t>(schema_->num_fields())) throw std::logic_error("Report column count mismatch");
        for (std::size_t i = 0; i < row.size(); ++i)
        {
          const auto& cell = row[i];
          if (text_.is_open())
          {
            if (i) text_ << '\t';
            if (std::holds_alternative<std::monostate>(cell)) text_ << "NA";
            else if (const auto* s = std::get_if<std::string>(&cell)) text_ << escaped(*s);
            else if (const auto* d = std::get_if<double>(&cell)) text_ << *d;
            else if (const auto* n = std::get_if<std::int64_t>(&cell)) text_ << *n;
            else text_ << (std::get<bool>(cell) ? "true" : "false");
          }
          if (writer_)
          {
            auto& builder = *builders_[i];
            if (std::holds_alternative<std::monostate>(cell)) check(builder.AppendNull());
            else if (const auto* s = std::get_if<std::string>(&cell)) check(static_cast<arrow::StringBuilder&>(builder).Append(*s));
            else if (const auto* d = std::get_if<double>(&cell)) check(static_cast<arrow::DoubleBuilder&>(builder).Append(*d));
            else if (const auto* n = std::get_if<std::int64_t>(&cell)) check(static_cast<arrow::Int64Builder&>(builder).Append(*n));
            else check(static_cast<arrow::BooleanBuilder&>(builder).Append(std::get<bool>(cell)));
          }
        }
        if (text_.is_open()) text_ << '\n';
        if (++buffered_ >= batch_limit_) flush();
      }

      void close()
      {
        flush();
        if (writer_) check(writer_->Close());
        if (output_) check(output_->Close());
        if (text_.is_open()) text_.close();
      }

    private:
      void flush()
      {
        if (writer_ && buffered_)
        {
          std::vector<std::shared_ptr<arrow::Array>> arrays;
          for (auto& builder : builders_) arrays.push_back(unwrap(builder->Finish()));
          check(writer_->WriteTable(*arrow::Table::Make(schema_, arrays), buffered_));
        }
        buffered_ = 0;
      }
      std::shared_ptr<arrow::Schema> schema_;
      std::ofstream text_;
      std::shared_ptr<arrow::io::FileOutputStream> output_;
      std::unique_ptr<parquet::arrow::FileWriter> writer_;
      std::vector<std::unique_ptr<arrow::ArrayBuilder>> builders_;
      std::size_t buffered_ = 0, batch_limit_ = 4096;
    };

    auto textField(const std::string& name) { return arrow::field(name, arrow::utf8()); }
    auto realField(const std::string& name) { return arrow::field(name, arrow::float64()); }
    auto intField(const std::string& name) { return arrow::field(name, arrow::int64()); }
    auto boolField(const std::string& name) { return arrow::field(name, arrow::boolean()); }

    bool passes(double run_q, double global_q, double threshold)
    {
      return std::isfinite(run_q) && std::isfinite(global_q) && run_q >= 0 && global_q >= 0
          && run_q <= threshold && global_q <= threshold;
    }

    const QuantifiedEntity* entity(const std::vector<QuantifiedEntity>& entities,
                                  const std::vector<std::uint32_t>& mapping, std::size_t precursor)
    {
      return precursor < mapping.size() && mapping[precursor] < entities.size() ? &entities[mapping[precursor]] : nullptr;
    }

    const EntityRun* evidence(const QuantifiedEntity* e, std::size_t run)
    {
      if (e) for (const auto& r : e->runs) if (r.run == run) return &r;
      return nullptr;
    }

    double selectedQ(const EntityRun& run, const ReportingOptions& options)
    {
      return options.inference_context == "experiment-wide" ? run.experiment_q_value : run.q_value;
    }

    void appendStatistics(Row& row, const QuantifiedEntity* entity, const EntityRun* run)
    {
      for (double value : {run ? run->q_value : na, run ? run->p_value : na, run ? run->pep : na,
                           run ? run->experiment_q_value : na, run ? run->experiment_p_value : na,
                           run ? run->experiment_pep : na, entity ? entity->global_q_value : na,
                           entity ? entity->global_p_value : na, entity ? entity->global_pep : na})
      {
        row.push_back(number(value));
      }
    }

    void appendStatisticFields(std::vector<std::shared_ptr<arrow::Field>>& fields, const std::string& prefix)
    {
      for (const std::string context : {"run", "experiment", "global"})
      {
        for (const std::string statistic : {"q_value", "p_value", "pep"})
        {
          fields.push_back(arrow::field(prefix + context + "_" + statistic, arrow::float64()));
        }
      }
    }

    std::string membership(const std::vector<std::string>& proteins)
    {
      std::string result;
      for (const auto& protein : proteins)
      {
        if (!result.empty()) result += ';';
        result += protein;
      }
      return result;
    }

    Row precursorMetadata(const AnalysisLibrary& library, std::size_t p, const GeneAnnotations& genes)
    {
      const auto& precursor = library.precursors[p];
      Cell sequence;
      // Preserve source text even when its modification syntax cannot be interpreted. Unknown
      // normalized metadata is null, never guessed from an unrelated entity or stripped by regex.
      try
      {
        const auto parsed = OpenMS::AASequence::fromString(precursor.sequence).toUnmodifiedString();
        if (!parsed.empty()) sequence = parsed;
      }
      catch (const OpenMS::Exception::BaseException&)
      {
      }
      const auto gene = genes.find(std::string(precursor.id));
      return {std::string(precursor.id), precursor.sequence, sequence, std::int64_t(precursor.charge),
              precursor.proteins < library.protein_sets.size() ? Cell(library.protein_sets[precursor.proteins]) : Cell{},
              gene != genes.end() && !gene->second.empty() ? Cell(gene->second) : Cell{}};
    }
  }

  void writeReports(const std::string& directory, const AnalysisLibrary& library,
                    const std::vector<RunResult>& runs, const Quantification& quantification,
                    const GeneAnnotations& genes, const ReportingOptions& options, const std::string& format)
  {
    if (format != "tsv" && format != "parquet" && format != "both") throw std::invalid_argument("Report format must be tsv, parquet, or both");
    for (const double threshold : {options.peakgroup_fdr, options.peptide_fdr, options.protein_fdr})
      if (!std::isfinite(threshold) || threshold < 0 || threshold > 1) throw std::invalid_argument("Report FDR thresholds must be in [0,1]");
    if (options.inference_context != "run-specific" && options.inference_context != "experiment-wide")
      throw std::invalid_argument("Report inference context must be run-specific or experiment-wide");
    std::vector<std::string> run_names;
    for (const auto& run : runs)
    {
      const auto name = ArrowIOHelpers::qpxRunFileName(run.run);
      if (name.empty() || std::find(run_names.begin(), run_names.end(), name) != run_names.end())
        throw std::invalid_argument("Report run stems must be nonempty and unique");
      for (const char* reserved : {"precursor_id", "sequence", "unmodified_sequence", "charge", "proteins", "genes", "protein_group", "modified_sequence"})
      {
        if (name == reserved) throw std::invalid_argument("Report run stem collides with matrix metadata column: " + name);
      }
      run_names.push_back(name);
    }

    // Index only observed retained winners, never the library by run Cartesian product.
    std::map<std::uint32_t, std::map<std::size_t, const RunWindow*>> winners;
    for (std::size_t r = 0; r < runs.size(); ++r)
      for (const auto& window : runs[r].windows)
      {
        if (window.precursor >= library.precursors.size()) throw std::out_of_range("Report precursor index");
        if (!std::isfinite(window.score) || window.score <= CrossFitScorer::NOT_KEPT) continue;
        auto& winner = winners[window.precursor][r];
        if (!winner || window.score > winner->score) winner = &window;
      }
    fs::create_directories(directory);
    Publication publication;
    const bool tsv = format != "parquet", parquet = format != "tsv";
    std::vector<std::shared_ptr<arrow::Field>> metadata{
      textField("precursor_id"), textField("sequence"), textField("unmodified_sequence"),
      intField("charge"), textField("proteins"), textField("genes")};
    auto long_fields = metadata;
    long_fields.insert(long_fields.begin(), {intField("feature_id"), textField("run")});
    const std::vector<std::shared_ptr<arrow::Field>> tail{
      realField("rt_apex"), realField("rt_start"), realField("rt_end"), realField("area"), textField("quantity_type"),
      realField("discriminant_score"), realField("run_q_value"), realField("global_q_value"), boolField("decoy"),
      textField("inference_context")};
    long_fields.insert(long_fields.end(), tail.begin(), tail.end());
    appendStatisticFields(long_fields, "peptide_");
    appendStatisticFields(long_fields, "protein_group_");
    long_fields.insert(long_fields.end(), {boolField("pass_peakgroup"), boolField("pass_peptide"), boolField("pass_protein_group")});
    Table precursors(publication, fs::path(directory) / "precursors", long_fields, tsv, parquet);
    auto matrix_fields = metadata;
    for (const auto& name : run_names) matrix_fields.push_back(realField(name));
    Table precursor_matrix(publication, fs::path(directory) / "precursor_matrix", matrix_fields, tsv, parquet);
    std::unordered_set<std::int64_t> feature_ids;
    for (const auto& [p, observed] : winners)
    {
      const auto& precursor = library.precursors[p];
      const auto* peptide = entity(quantification.peptides, quantification.precursor_peptide, p);
      const auto* group = entity(quantification.protein_groups, quantification.precursor_protein_group, p);
      const double global = p < quantification.global_precursor_q.size() ? quantification.global_precursor_q[p] : na;
      const auto values = precursorMetadata(library, p, genes);
      auto matrix = values;
      matrix.resize(metadata.size() + runs.size());
      bool eligible = false;
      for (const auto& [r, window] : observed)
      {
        const auto* pe = evidence(peptide, r);
        const auto* pg = evidence(group, r);
        const bool pass_peak = !precursor.decoy && passes(window->q_value, global, options.peakgroup_fdr);
        const bool pass_peptide = pass_peak && peptide && !peptide->decoy && pe
                               && passes(selectedQ(*pe, options), peptide->global_q_value, options.peptide_fdr);
        const bool pass_group = pass_peak && group && !group->decoy && pg
                             && passes(selectedQ(*pg, options), group->global_q_value, options.protein_fdr);
        const auto feature_id = outputFeatureId(runs[r].run, std::string(precursor.id), window->candidate);
        if (!feature_ids.insert(feature_id).second)
        {
          throw std::runtime_error("OpenDIA report: duplicate/colliding feature_id " + std::to_string(feature_id));
        }
        Row row{feature_id, run_names[r]};
        row.insert(row.end(), values.begin(), values.end());
        const Row rest{number(window->peak.rt_apex), number(window->peak.rt_start), number(window->peak.rt_end),
                       number(window->peak.area), std::string("raw_ms2_area"), number(window->score), number(window->q_value),
                       number(global), precursor.decoy, options.inference_context};
        row.insert(row.end(), rest.begin(), rest.end());
        // Preserve native context statistics even when another context supplies the gate.
        appendStatistics(row, peptide, pe);
        appendStatistics(row, group, pg);
        row.insert(row.end(), {pass_peak, pass_peptide, pass_group});
        precursors.append(row);
        if (pass_peak && std::isfinite(window->peak.area))
        {
          matrix[metadata.size() + r] = window->peak.area;
          eligible = true;
        }
      }
      if (eligible) precursor_matrix.append(matrix);
    }
    precursors.close();
    precursor_matrix.close();

    // Quantification supplies areas after precursor run/global cuts. Apply only the selected
    // entity context and global entity cuts here, without cross-gating peptide and protein levels.
    const auto write_entities = [&](const std::vector<QuantifiedEntity>& entities, const std::string& stem,
                                    const std::string& id_column, double threshold)
    {
      std::vector<std::shared_ptr<arrow::Field>> fields{textField(id_column), textField("proteins")};
      auto wide_fields = fields;
      for (const auto& name : run_names) wide_fields.push_back(realField(name));
      fields.insert(fields.end(), {textField("run"), realField("discriminant_score")});
      appendStatisticFields(fields, "");
      fields.insert(fields.end(), {boolField("decoy"), realField("area"),
                                  textField("quantity_type"), textField("inference_context"), boolField("pass")});
      Table long_table(publication, fs::path(directory) / (stem == "peptide" ? "peptides" : "protein_groups"), fields, tsv, parquet);
      Table matrix(publication, fs::path(directory) / (stem + "_matrix"), wide_fields, tsv, parquet);
      for (const auto& e : entities)
      {
        const std::string proteins = membership(e.proteins);
        Row wide{stem == "protein_group" && !proteins.empty() ? proteins : e.id, proteins};
        wide.resize(2 + runs.size());
        bool eligible = false;
        for (const auto& r : e.runs)
        {
          if (r.run >= runs.size()) throw std::out_of_range("Report entity run index");
          const bool pass = !e.decoy && passes(selectedQ(r, options), e.global_q_value, threshold);
          const Cell area = r.has_area ? number(r.area) : Cell{};
          Row row{e.id, proteins, run_names[r.run], number(r.score)};
          appendStatistics(row, &e, &r);
          row.insert(row.end(), {e.decoy, area, std::string("raw_ms2_area"), options.inference_context, pass});
          long_table.append(row);
          if (pass && !std::holds_alternative<std::monostate>(area))
          {
            wide[2 + r.run] = area;
            eligible = true;
          }
        }
        if (eligible) matrix.append(wide);
      }
      long_table.close();
      matrix.close();
    };
    write_entities(quantification.peptides, "peptide", "modified_sequence", options.peptide_fdr);
    write_entities(quantification.protein_groups, "protein_group", "protein_group", options.protein_fdr);
    publication.publish();
  }
}
