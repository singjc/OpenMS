// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#include <arrow/api.h>
#include <arrow/io/file.h>
#include <filesystem>
#include <iostream>
#include <parquet/arrow/writer.h>
#include <results/Annotations.h>
#include <sqlite3.h>
#include <stdexcept>

using namespace OpenMS::OpenDIA;
namespace fs = std::filesystem;
namespace
{
void require(bool ok, const std::string& message)
{
  if (! ok) throw std::runtime_error(message);
}
void sql(sqlite3* db, const char* query)
{
  char* error = nullptr;
  if (sqlite3_exec(db, query, nullptr, nullptr, &error) != SQLITE_OK)
  {
    const std::string message = error ? error : "SQLite failure";
    sqlite3_free(error);
    throw std::runtime_error(message);
  }
}
void parquetFile(const fs::path& path, const std::string& gene_column, bool wrong_type = false)
{
  arrow::StringBuilder ids, genes;
  for (const auto& id : {"0042", "0042", "decoy-original", "null-gene", "empty-gene"})
    require(ids.Append(id).ok(), "append id");
  require(genes.Append("Z;A").ok(), "append genes");
  require(genes.Append("A;B;NA").ok(), "append duplicate genes");
  require(genes.Append("D").ok(), "append decoy gene");
  require(genes.AppendNull().ok(), "append null gene");
  require(genes.Append("").ok(), "append empty gene");
  std::shared_ptr<arrow::Array> id_array, gene_array;
  require(ids.Finish(&id_array).ok() && genes.Finish(&gene_array).ok(), "finish arrays");
  // An unrelated leading column distinguishes projection from row-group selection.
  std::vector<std::shared_ptr<arrow::Field>> fields {arrow::field("unrelated", arrow::utf8()), arrow::field("precursor_id", arrow::utf8())};
  std::vector<std::shared_ptr<arrow::Array>> arrays {gene_array, id_array};
  if (wrong_type)
  {
    fields[1] = arrow::field("precursor_id", arrow::int64());
    const auto null_ids = arrow::MakeArrayOfNull(arrow::int64(), 5);
    require(null_ids.ok(), "make malformed ID array");
    arrays[1] = *null_ids;
  }
  if (! gene_column.empty())
  {
    fields.push_back(arrow::field(gene_column, arrow::utf8()));
    arrays.push_back(gene_array);
  }
  const auto opened = arrow::io::FileOutputStream::Open(path.string());
  require(opened.ok(), "open parquet fixture");
  require(parquet::arrow::WriteTable(*arrow::Table::Make(arrow::schema(fields), arrays), arrow::default_memory_pool(), *opened, 2).ok(),
          "write fixture");
  require((*opened)->Close().ok(), "close fixture");
}
} // namespace
int main(int argc, char** argv)
{
  try
  {
    require(argc > 1, "Expected fixture directory");
    fs::create_directories(argv[1]);
    const fs::path dir(argv[1]);
    const auto pqp = dir / "annotations.pqp";
    fs::remove(pqp);
    sqlite3* raw = nullptr;
    require(sqlite3_open(pqp.c_str(), &raw) == SQLITE_OK, "create PQP");
    const std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
    sql(db.get(), "CREATE TABLE PRECURSOR(ID INTEGER, TRAML_ID TEXT);"
                  "CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID INTEGER, PEPTIDE_ID INTEGER);"
                  "INSERT INTO PRECURSOR VALUES(42,'not-the-producer-id'),(43,'DECOY_not-the-producer-id');"
                  "INSERT INTO PRECURSOR_PEPTIDE_MAPPING VALUES(42,1),(43,2);");
    require(readGeneAnnotations(pqp.string()).empty(), "PQP without optional genes stays absent");
    sql(db.get(), "CREATE TABLE GENE(ID INTEGER, GENE_NAME TEXT);"
                  "CREATE TABLE PEPTIDE_GENE_MAPPING(PEPTIDE_ID INTEGER, GENE_ID INTEGER);"
                  "INSERT INTO GENE VALUES(1,'Z'),(2,'A'),(3,'A;B'),(4,NULL),(5,'NA'),(6,''),(7,'DECOY_G');"
                  "INSERT INTO PEPTIDE_GENE_MAPPING VALUES(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(2,7);");
    const auto genes = readGeneAnnotations(pqp.string());
    require(genes.size() == 2 && genes.at("42") == "A;B;Z", "PQP genes deduplicate and sort under original precursor ID");
    require(genes.at("43") == "DECOY_G" && ! genes.count("not-the-producer-id"), "preserve decoy and producer IDs without inferring aliases");
    for (const auto& column : {std::string("gene_name"), std::string("genes"), std::string()})
    {
      const auto file = dir / ("annotations_" + column + ".parquet");
      parquetFile(file, column);
      const auto read = readGeneAnnotations(file.string());
      if (column.empty()) require(read.empty(), "old parquet without gene metadata stays absent");
      else
      {
        require(read.size() == 2 && read.at("0042") == "A;B;Z", "Parquet preserves exact string precursor ID and merges row groups");
        require(read.at("decoy-original") == "D", "Parquet keeps explicit decoy annotation");
      }
    }
    const auto malformed = dir / "bad_type.parquet";
    parquetFile(malformed, "genes", true);
    bool rejected = false;
    try
    {
      readGeneAnnotations(malformed.string());
    }
    catch (const std::runtime_error&)
    {
      rejected = true;
    }
    require(rejected, "malformed annotation types rejected");
    rejected = false;
    try
    {
      readGeneAnnotations((dir / "does_not_exist.pqp").string());
    }
    catch (const std::runtime_error&)
    {
      rejected = true;
    }
    require(rejected, "missing input rejected without creating a database");
    require(! fs::exists(dir / "does_not_exist.pqp"), "read-only loading never creates a missing PQP");
    std::cout << "Annotation checks passed\n";
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
