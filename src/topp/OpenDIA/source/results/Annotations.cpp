// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#include <arrow/api.h>
#include <arrow/io/file.h>
#include <library/Library.h>
#include <memory>
#include <numeric>
#include <parquet/arrow/reader.h>
#include <results/Annotations.h>
#include <set>
#include <sqlite3.h>
#include <stdexcept>

namespace OpenMS::OpenDIA
{
namespace
{
  void add(std::map<std::string, std::set<std::string>>& values, const std::string& id, const std::string& gene)
  {
    if (id.empty() || gene.empty() || gene == "NA") return;
    std::size_t begin = 0;
    while (begin < gene.size())
    {
      const auto end = gene.find(';', begin);
      const auto name = gene.substr(begin, end == std::string::npos ? end : end - begin);
      if (! name.empty() && name != "NA") values[id].insert(name);
      if (end == std::string::npos) break;
      begin = end + 1;
    }
  }
} // namespace

GeneAnnotations readGeneAnnotations(const std::string& path)
{
  std::map<std::string, std::set<std::string>> values;
  if (isParquetLibrary(path))
  {
    const auto opened = arrow::io::ReadableFile::Open(path);
    if (! opened.ok()) throw std::runtime_error("Cannot read output annotations: " + opened.status().ToString());
    auto reader_result = parquet::arrow::OpenFile(*opened, arrow::default_memory_pool());
    if (! reader_result.ok()) throw std::runtime_error(reader_result.status().ToString());
    auto reader = std::move(*reader_result);
    std::shared_ptr<arrow::Schema> schema;
    auto status = reader->GetSchema(&schema);
    if (! status.ok()) throw std::runtime_error(status.ToString());
    const int id_column = schema->GetFieldIndex("precursor_id");
    int gene_column = schema->GetFieldIndex("gene_name");
    if (gene_column < 0) gene_column = schema->GetFieldIndex("genes");
    // Existing version-1 library caches did not retain genes. Do not invent an annotation.
    if (id_column < 0 || gene_column < 0) return {};
    std::vector<int> row_groups(reader->num_row_groups());
    std::iota(row_groups.begin(), row_groups.end(), 0);
    auto batches_result = reader->GetRecordBatchReader(row_groups, {id_column, gene_column});
    if (! batches_result.ok()) throw std::runtime_error(batches_result.status().ToString());
    auto batches = std::move(*batches_result);
    std::shared_ptr<arrow::RecordBatch> batch;
    while (true)
    {
      status = batches->ReadNext(&batch);
      if (! status.ok()) throw std::runtime_error(status.ToString());
      if (! batch) break;
      if (batch->column(0)->type_id() != arrow::Type::STRING || batch->column(1)->type_id() != arrow::Type::STRING)
        throw std::runtime_error("Library output annotation columns must be strings");
      const auto& ids = static_cast<const arrow::StringArray&>(*batch->column(0));
      const auto& genes = static_cast<const arrow::StringArray&>(*batch->column(1));
      for (std::int64_t i = 0; i < batch->num_rows(); ++i)
        if (! ids.IsNull(i) && ! genes.IsNull(i)) add(values, ids.GetString(i), genes.GetString(i));
    }
  }
  else
  {
    sqlite3* raw = nullptr;
    if (sqlite3_open_v2(path.c_str(), &raw, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK)
    {
      const std::string error = raw ? sqlite3_errmsg(raw) : "allocation failure";
      if (raw) sqlite3_close(raw);
      throw std::runtime_error("Cannot read library output annotations: " + error);
    }
    const std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
    const auto prepare = [&](const char* query) {
      sqlite3_stmt* statement = nullptr;
      if (sqlite3_prepare_v2(db.get(), query, -1, &statement, nullptr) != SQLITE_OK)
        throw std::runtime_error("Cannot read library output annotations: " + std::string(sqlite3_errmsg(db.get())));
      return std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)>(statement, sqlite3_finalize);
    };
    auto tables = prepare("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name IN ('GENE','PEPTIDE_GENE_MAPPING')");
    if (sqlite3_step(tables.get()) != SQLITE_ROW || sqlite3_column_int(tables.get(), 0) != 2) return {};
    // Keep PQP PRECURSOR_ID, the producer key used by AnalysisLibrary. TRAML_ID is a
    // pairing alias and must not replace the output precursor identifier.
    auto rows = prepare("SELECT CAST(M.PRECURSOR_ID AS TEXT), G.GENE_NAME FROM PRECURSOR_PEPTIDE_MAPPING M "
                        "JOIN PEPTIDE_GENE_MAPPING P ON M.PEPTIDE_ID=P.PEPTIDE_ID JOIN GENE G ON P.GENE_ID=G.ID "
                        "ORDER BY M.PRECURSOR_ID, G.GENE_NAME");
    int status = SQLITE_OK;
    while ((status = sqlite3_step(rows.get())) == SQLITE_ROW)
    {
      if (sqlite3_column_type(rows.get(), 0) == SQLITE_NULL || sqlite3_column_type(rows.get(), 1) == SQLITE_NULL) continue;
      add(values, reinterpret_cast<const char*>(sqlite3_column_text(rows.get(), 0)),
          reinterpret_cast<const char*>(sqlite3_column_text(rows.get(), 1)));
    }
    if (status != SQLITE_DONE) throw std::runtime_error(sqlite3_errmsg(db.get()));
  }
  GeneAnnotations result;
  for (const auto& [id, names] : values)
  {
    auto& text = result[id];
    for (const auto& name : names)
    {
      if (! text.empty()) text += ';';
      text += name;
    }
  }
  return result;
}
} // namespace OpenMS::OpenDIA
