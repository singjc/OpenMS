// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

// The OpenDIA Parquet library: one row per fragment, in the order OpenDIA uses it. A PQP is converted on first use.

#include <library/Library.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <unistd.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <filesystem>
#include <numeric>
#include <unordered_set>

namespace OpenMS::OpenDIA
{
  namespace
  {
    constexpr const char* FORMAT_KEY = "opendia.library.format";
    constexpr const char* FORMAT_VALUE = "opendia_library";
    constexpr const char* SCHEMA_KEY = "opendia.library.schema_version";
    constexpr const char* SCHEMA_VALUE = "1";
    constexpr const char* FRAGMENTS_KEY = "opendia.library.fragments";
    // The PQP a cached copy was made from: its size, modification time and the version of the PQP reading it went
    // through (bump READER_VALUE whenever readLibrary changes what it returns, so old copies are rebuilt).
    constexpr const char* SOURCE_SIZE_KEY = "opendia.source.size";
    constexpr const char* SOURCE_MTIME_KEY = "opendia.source.mtime";
    constexpr const char* READER_KEY = "opendia.source.reader_version";
    constexpr const char* READER_VALUE = "1";
    constexpr std::int64_t ROW_GROUP = 1 << 20;

    // Column order of schema version 1.
    enum Column
    {
      PRECURSOR_INDEX,   ///< Position of the precursor in Library::traml_ids (the PQP PRECURSOR order).
      PRECURSOR_ID,
      TRAML_ID,
      DECOY,
      PRECURSOR_MZ,
      PRECURSOR_CHARGE,
      LIBRARY_RT,
      PRECURSOR_IM,
      SEQUENCE,
      PROTEINS,
      TRANSITION_ID,     ///< Null on a row without fragment.
      PRODUCT_MZ,
      LIBRARY_INTENSITY,
      FRAGMENT_TYPE,
      FRAGMENT_ORDINAL,
      FRAGMENT_CHARGE,
      FRAGMENT_RANK,     ///< 1 = most intense fragment of the precursor (ties in library order).
      COLUMNS
    };

    std::shared_ptr<arrow::Schema> schemaV1()
    {
      return arrow::schema({
        arrow::field("precursor_index", arrow::int32(), false),
        arrow::field("precursor_id", arrow::utf8(), false),
        arrow::field("traml_id", arrow::utf8(), false),
        arrow::field("decoy", arrow::boolean()),
        arrow::field("precursor_mz", arrow::float64()),
        arrow::field("precursor_charge", arrow::int32()),
        arrow::field("library_rt", arrow::float64()),
        arrow::field("precursor_im", arrow::float64()),
        arrow::field("sequence", arrow::utf8()),
        arrow::field("proteins", arrow::list(arrow::utf8())),
        arrow::field("transition_id", arrow::utf8()),
        arrow::field("product_mz", arrow::float64()),
        arrow::field("library_intensity", arrow::float64()),
        arrow::field("fragment_type", arrow::utf8()),
        arrow::field("fragment_ordinal", arrow::int16()),
        arrow::field("fragment_charge", arrow::int8()),
        arrow::field("fragment_rank", arrow::int32()),
      });
    }

    [[noreturn]] void notWritable(const std::string& path, const arrow::Status& status)
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path + ": " + status.ToString());
    }

    [[noreturn]] void notReadable(const std::string& path, const std::string& message)
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path + ": " + message);
    }

    void check(const arrow::Status& status, const std::string& path)
    {
      if (!status.ok()) notWritable(path, status);
    }

    /// The builders of one record batch.
    struct Batch
    {
      arrow::Int32Builder precursor_index;
      arrow::StringBuilder precursor_id, traml_id, sequence, transition_id, fragment_type;
      arrow::BooleanBuilder decoy;
      arrow::DoubleBuilder precursor_mz, library_rt, precursor_im, product_mz, library_intensity;
      arrow::Int32Builder precursor_charge, fragment_rank;
      arrow::Int16Builder fragment_ordinal;
      arrow::Int8Builder fragment_charge;
      std::shared_ptr<arrow::StringBuilder> protein = std::make_shared<arrow::StringBuilder>();
      arrow::ListBuilder proteins{arrow::default_memory_pool(), protein};
      std::int64_t rows = 0;

      void precursor(std::int32_t index, const std::string& id, const std::string& traml, const std::string& path)
      {
        check(precursor_index.Append(index), path);
        check(precursor_id.Append(id), path);
        check(traml_id.Append(traml), path);
        ++rows;
      }

      void fragment(const LibraryTransition& t, const LibraryCompound& c, std::int32_t rank, const std::string& path)
      {
        check(decoy.Append(t.decoy), path);
        check(precursor_mz.Append(t.precursor_mz), path);
        check(precursor_charge.Append(c.charge), path);
        check(library_rt.Append(c.rt), path);
        check(precursor_im.Append(t.precursor_im), path);
        check(sequence.Append(c.sequence), path);
        check(proteins.Append(), path);
        for (const auto& p : c.protein_refs) check(protein->Append(p), path);
        check(transition_id.Append(t.transition_name), path);
        check(product_mz.Append(t.product_mz), path);
        check(library_intensity.Append(t.library_intensity), path);
        check(fragment_type.Append(t.fragment_type), path);
        check(fragment_ordinal.Append(t.fragment_nr), path);
        check(fragment_charge.Append(t.fragment_charge), path);
        check(fragment_rank.Append(rank), path);
      }

      void noFragment(const std::string& path)
      {
        check(decoy.AppendNull(), path);
        check(precursor_mz.AppendNull(), path);
        check(precursor_charge.AppendNull(), path);
        check(library_rt.AppendNull(), path);
        check(precursor_im.AppendNull(), path);
        check(sequence.AppendNull(), path);
        check(proteins.AppendNull(), path);
        check(transition_id.AppendNull(), path);
        check(product_mz.AppendNull(), path);
        check(library_intensity.AppendNull(), path);
        check(fragment_type.AppendNull(), path);
        check(fragment_ordinal.AppendNull(), path);
        check(fragment_charge.AppendNull(), path);
        check(fragment_rank.AppendNull(), path);
      }

      std::shared_ptr<arrow::RecordBatch> finish(const std::shared_ptr<arrow::Schema>& schema, const std::string& path)
      {
        std::vector<std::shared_ptr<arrow::Array>> arrays(COLUMNS);
        check(precursor_index.Finish(&arrays[PRECURSOR_INDEX]), path);
        check(precursor_id.Finish(&arrays[PRECURSOR_ID]), path);
        check(traml_id.Finish(&arrays[TRAML_ID]), path);
        check(decoy.Finish(&arrays[DECOY]), path);
        check(precursor_mz.Finish(&arrays[PRECURSOR_MZ]), path);
        check(precursor_charge.Finish(&arrays[PRECURSOR_CHARGE]), path);
        check(library_rt.Finish(&arrays[LIBRARY_RT]), path);
        check(precursor_im.Finish(&arrays[PRECURSOR_IM]), path);
        check(sequence.Finish(&arrays[SEQUENCE]), path);
        check(proteins.Finish(&arrays[PROTEINS]), path);
        check(transition_id.Finish(&arrays[TRANSITION_ID]), path);
        check(product_mz.Finish(&arrays[PRODUCT_MZ]), path);
        check(library_intensity.Finish(&arrays[LIBRARY_INTENSITY]), path);
        check(fragment_type.Finish(&arrays[FRAGMENT_TYPE]), path);
        check(fragment_ordinal.Finish(&arrays[FRAGMENT_ORDINAL]), path);
        check(fragment_charge.Finish(&arrays[FRAGMENT_CHARGE]), path);
        check(fragment_rank.Finish(&arrays[FRAGMENT_RANK]), path);
        const auto batch = arrow::RecordBatch::Make(schema, rows, std::move(arrays));
        rows = 0;
        return batch;
      }
    };

    std::unordered_map<std::string, std::string> metadataOf(const arrow::Schema& schema)
    {
      std::unordered_map<std::string, std::string> out;
      if (const auto& m = schema.metadata())
      {
        for (std::int64_t i = 0; i < m->size(); ++i) out.emplace(m->key(i), m->value(i));
      }
      return out;
    }

    std::unique_ptr<parquet::arrow::FileReader> openReader(const std::string& path)
    {
      auto file = arrow::io::ReadableFile::Open(path);
      if (!file.ok()) notReadable(path, file.status().ToString());
      auto reader = parquet::arrow::OpenFile(*file, arrow::default_memory_pool());
      if (!reader.ok()) notReadable(path, reader.status().ToString());
      return std::move(*reader);
    }

    std::string view(const arrow::StringArray& a, std::int64_t i) { return a.IsNull(i) ? std::string() : std::string(a.GetView(i)); }

    template <typename Array> auto get(const Array& a, std::int64_t i, decltype(a.Value(0)) missing)
    {
      return a.IsNull(i) ? missing : a.Value(i);
    }

    /// PQP size and modification time (ns), the key of its Parquet copy.
    std::pair<std::string, std::string> sourceKey(const std::string& pqp)
    {
      return {std::to_string(std::filesystem::file_size(pqp)),
              std::to_string(std::filesystem::last_write_time(pqp).time_since_epoch().count())};
    }
  }

  bool isParquetLibrary(const std::string& path)
  {
    std::string extension = std::filesystem::path(path).extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return extension == ".parquet";
  }

  std::string parquetCachePath(const std::string& pqp)
  {
    const std::filesystem::path p(pqp);
    return (p.parent_path() / (p.stem().string() + ".opendia.parquet")).string();
  }

  void writeParquetLibrary(const Library& library, const std::string& path, const std::vector<std::pair<std::string, std::string>>& metadata)
  {
    // Precursor index (first position in traml_ids) and compound of every precursor ID.
    std::unordered_map<std::string, std::int32_t> precursor_index;
    precursor_index.reserve(library.traml_ids.size());
    for (std::size_t i = 0; i < library.traml_ids.size(); ++i) precursor_index.try_emplace(library.traml_ids[i].first, static_cast<std::int32_t>(i));
    std::unordered_map<std::string, std::uint32_t> compound_index;
    compound_index.reserve(library.compounds.size());
    for (std::size_t i = 0; i < library.compounds.size(); ++i) compound_index.try_emplace(library.compounds[i].id, static_cast<std::uint32_t>(i));

    // Intensity rank within each precursor (1 = most intense; ties in library order).
    const auto& transitions = library.transitions;
    std::vector<std::int32_t> rank(transitions.size(), 0);
    std::vector<std::uint32_t> compound_of(transitions.size());
    std::vector<bool> covered(library.traml_ids.size(), false);
    {
      std::vector<std::vector<std::uint32_t>> members(library.compounds.size());
      for (std::size_t i = 0; i < transitions.size(); ++i)
      {
        const auto c = compound_index.find(transitions[i].peptide_ref);
        const auto p = precursor_index.find(transitions[i].peptide_ref);
        if (c == compound_index.end() || p == precursor_index.end())
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Fragment of a precursor without compound or TRAML_ID",
                                        transitions[i].peptide_ref);
        compound_of[i] = c->second;
        covered[p->second] = true;
        members[c->second].push_back(static_cast<std::uint32_t>(i));
      }
      for (auto& m : members)
      {
        std::stable_sort(m.begin(), m.end(), [&](std::uint32_t a, std::uint32_t b) { return transitions[a].library_intensity > transitions[b].library_intensity; });
        for (std::size_t r = 0; r < m.size(); ++r) rank[m[r]] = static_cast<std::int32_t>(r + 1);
      }
    }

    // Host and process make the temporary name unique also when several nodes convert the same PQP on a shared disk.
    char host[256] = {};
    ::gethostname(host, sizeof(host) - 1);
    const std::string partial = path + "." + host + "." + std::to_string(::getpid()) + ".tmp";
    std::error_code ignored;
    std::filesystem::remove(partial, ignored);
    try
    {
      std::vector<std::string> keys{FORMAT_KEY, SCHEMA_KEY, FRAGMENTS_KEY};
      std::vector<std::string> values{FORMAT_VALUE, SCHEMA_VALUE, std::to_string(transitions.size())};
      for (const auto& [k, v] : metadata)
      {
        keys.push_back(k);
        values.push_back(v);
      }
      const auto schema = schemaV1()->WithMetadata(arrow::key_value_metadata(keys, values));
      auto out = arrow::io::FileOutputStream::Open(partial);
      if (!out.ok()) notWritable(partial, out.status());
      const auto properties = parquet::WriterProperties::Builder()
                                .compression(parquet::Compression::ZSTD)
                                ->max_row_group_length(ROW_GROUP)
                                ->build();
      const auto arrow_properties = parquet::ArrowWriterProperties::Builder().store_schema()->build();
      auto writer = parquet::arrow::FileWriter::Open(*schema, arrow::default_memory_pool(), *out, properties, arrow_properties);
      if (!writer.ok()) notWritable(partial, writer.status());

      Batch batch;
      const auto flush = [&](bool force)
      {
        if (batch.rows == 0 || (!force && batch.rows < ROW_GROUP)) return;
        check((*writer)->WriteRecordBatch(*batch.finish(schema, partial)), partial);
      };
      for (std::size_t i = 0; i < transitions.size(); ++i)
      {
        const auto& t = transitions[i];
        const std::int32_t p = precursor_index.at(t.peptide_ref);
        batch.precursor(p, t.peptide_ref, library.traml_ids[p].second, partial);
        batch.fragment(t, library.compounds[compound_of[i]], rank[i], partial);
        flush(false);
      }
      // Precursors without fragments (and repeated PRECURSOR IDs) keep their TRAML_ID, after all fragments.
      for (std::size_t p = 0; p < library.traml_ids.size(); ++p)
      {
        if (covered[p] && precursor_index.at(library.traml_ids[p].first) == static_cast<std::int32_t>(p)) continue;
        batch.precursor(static_cast<std::int32_t>(p), library.traml_ids[p].first, library.traml_ids[p].second, partial);
        batch.noFragment(partial);
        flush(false);
      }
      flush(true);
      check((*writer)->Close(), partial);
      check((*out)->Close(), partial);
      std::filesystem::rename(partial, path);
    }
    catch (...)
    {
      std::filesystem::remove(partial, ignored);
      throw;
    }
  }

  struct ParquetLibraryWriter::Impl
  {
    std::string path, partial;
    std::shared_ptr<arrow::Schema> schema;
    std::shared_ptr<arrow::io::FileOutputStream> out;
    std::unique_ptr<parquet::arrow::FileWriter> writer;
    Batch batch;
    bool closed = false;

    void flush(bool force)
    {
      if (batch.rows == 0 || (!force && batch.rows < ROW_GROUP)) return;
      check(writer->WriteRecordBatch(*batch.finish(schema, partial)), partial);
    }
  };

  ParquetLibraryWriter::ParquetLibraryWriter(const std::string& path, const std::vector<std::pair<std::string, std::string>>& metadata) :
    impl_(std::make_unique<Impl>())
  {
    char host[256] = {};
    ::gethostname(host, sizeof(host) - 1);
    impl_->path = path;
    impl_->partial = path + "." + host + "." + std::to_string(::getpid()) + ".tmp";
    std::error_code ignored;
    std::filesystem::remove(impl_->partial, ignored);
    // The fragment count is not known up front; it is only a capacity hint for the readers, so it is left out.
    std::vector<std::string> keys{FORMAT_KEY, SCHEMA_KEY};
    std::vector<std::string> values{FORMAT_VALUE, SCHEMA_VALUE};
    for (const auto& [k, v] : metadata)
    {
      keys.push_back(k);
      values.push_back(v);
    }
    impl_->schema = schemaV1()->WithMetadata(arrow::key_value_metadata(keys, values));
    auto out = arrow::io::FileOutputStream::Open(impl_->partial);
    if (!out.ok()) notWritable(impl_->partial, out.status());
    impl_->out = *out;
    const auto properties = parquet::WriterProperties::Builder().compression(parquet::Compression::ZSTD)->max_row_group_length(ROW_GROUP)->build();
    const auto arrow_properties = parquet::ArrowWriterProperties::Builder().store_schema()->build();
    auto writer = parquet::arrow::FileWriter::Open(*impl_->schema, arrow::default_memory_pool(), impl_->out, properties, arrow_properties);
    if (!writer.ok()) notWritable(impl_->partial, writer.status());
    impl_->writer = std::move(*writer);
  }

  ParquetLibraryWriter::~ParquetLibraryWriter()
  {
    if (impl_ && !impl_->closed)
    {
      if (impl_->writer) (void)impl_->writer->Close();
      if (impl_->out) (void)impl_->out->Close();
      std::error_code ignored;
      std::filesystem::remove(impl_->partial, ignored);
    }
  }

  void ParquetLibraryWriter::add(const Precursor& p, const std::vector<Fragment>& fragments)
  {
    if (fragments.empty())
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor without fragments", std::string(p.id));
    if (precursors_ >= static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Too many library precursors", std::string(p.id));
    std::vector<std::uint32_t> order(fragments.size());
    std::iota(order.begin(), order.end(), 0u);
    std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) { return fragments[a].intensity > fragments[b].intensity; });
    std::vector<std::int32_t> rank(fragments.size());
    for (std::size_t r = 0; r < order.size(); ++r) rank[order[r]] = static_cast<std::int32_t>(r + 1);
    auto& b = impl_->batch;
    const auto& path = impl_->partial;
    const auto index = static_cast<std::int32_t>(precursors_);
    for (std::size_t i = 0; i < fragments.size(); ++i)
    {
      const auto& f = fragments[i];
      check(b.precursor_index.Append(index), path);
      check(b.precursor_id.Append(p.id.data(), static_cast<std::int32_t>(p.id.size())), path);
      check(b.traml_id.Append(p.traml_id.data(), static_cast<std::int32_t>(p.traml_id.size())), path);
      ++b.rows;
      check(b.decoy.Append(p.decoy), path);
      check(b.precursor_mz.Append(p.mz), path);
      check(b.precursor_charge.Append(p.charge), path);
      check(b.library_rt.Append(p.rt), path);
      check(b.precursor_im.Append(p.im), path);
      check(b.sequence.Append(p.sequence.data(), static_cast<std::int32_t>(p.sequence.size())), path);
      check(b.proteins.Append(), path);
      for (const auto& protein : *p.proteins) check(b.protein->Append(protein), path);
      check(b.transition_id.Append(f.transition_id), path);
      check(b.product_mz.Append(f.product_mz), path);
      check(b.library_intensity.Append(f.intensity), path);
      check(b.fragment_type.Append(f.type.data(), static_cast<std::int32_t>(f.type.size())), path);
      check(b.fragment_ordinal.Append(f.ordinal), path);
      check(b.fragment_charge.Append(f.charge), path);
      check(b.fragment_rank.Append(rank[i]), path);
    }
    ++precursors_;
    fragments_ += fragments.size();
    impl_->flush(false);
  }

  void ParquetLibraryWriter::close()
  {
    impl_->flush(true);
    check(impl_->writer->Close(), impl_->partial);
    check(impl_->out->Close(), impl_->partial);
    std::filesystem::rename(impl_->partial, impl_->path);
    impl_->closed = true;
  }

  std::unordered_map<std::string, std::string> readParquetLibraryMetadata(const std::string& path)
  {
    std::error_code error;
    if (!std::filesystem::is_regular_file(path, error)) return {};
    try
    {
      const auto reader = openReader(path);
      std::shared_ptr<arrow::Schema> schema;
      if (!reader->GetSchema(&schema).ok() || !schema) return {};
      return metadataOf(*schema);
    }
    catch (const Exception::BaseException&)
    {
      return {};
    }
  }

  std::vector<std::pair<std::string, std::string>> visitParquetLibrary(const std::string& path, LibraryVisitor& visitor)
  {
    const auto reader = openReader(path);
    // Single-threaded decode: on the 600k-precursor Astral libraries threaded decoding saved 1.7 s but left 0.62 GB
    // more peak resident memory (worker-thread malloc arenas).
    reader->set_use_threads(false);
    std::shared_ptr<arrow::Schema> schema;
    const auto status = reader->GetSchema(&schema);
    if (!status.ok()) notReadable(path, status.ToString());
    const auto metadata = metadataOf(*schema);
    const auto format = metadata.find(FORMAT_KEY);
    const auto version = metadata.find(SCHEMA_KEY);
    if (format == metadata.end() || format->second != FORMAT_VALUE)
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path + ": not an OpenDIA Parquet library", FORMAT_KEY);
    if (version == metadata.end() || version->second != SCHEMA_VALUE)
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path + ": unsupported OpenDIA library schema version",
                                    version == metadata.end() ? "" : version->second);
    const auto expected = schemaV1();
    if (schema->num_fields() != expected->num_fields())
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path + ": unexpected number of columns", std::to_string(schema->num_fields()));
    for (int c = 0; c < expected->num_fields(); ++c)
    {
      if (schema->field(c)->name() != expected->field(c)->name() || !schema->field(c)->type()->Equals(*expected->field(c)->type()))
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path + ": unexpected column", schema->field(c)->ToString());
    }

    if (const auto n = metadata.find(FRAGMENTS_KEY); n != metadata.end()) visitor.expect(std::stoull(n->second));
    std::vector<std::pair<std::string, std::string>> traml_ids;
    std::vector<bool> has_traml;
    std::vector<std::string_view> protein_views;
    for (int g = 0; g < reader->num_row_groups(); ++g)
    {
      auto read = reader->ReadRowGroup(g);
      if (!read.ok()) notReadable(path, read.status().ToString());
      const std::shared_ptr<arrow::Table> table = std::move(*read);
      arrow::TableBatchReader batches(*table);
      std::shared_ptr<arrow::RecordBatch> batch;
      while (true)
      {
        const auto next = batches.ReadNext(&batch);
        if (!next.ok()) notReadable(path, next.ToString());
        if (!batch) break;
        const auto column = [&](int c) { return batch->column(c).get(); };
        const auto& index = static_cast<const arrow::Int32Array&>(*column(PRECURSOR_INDEX));
        const auto& id = static_cast<const arrow::StringArray&>(*column(PRECURSOR_ID));
        const auto& traml = static_cast<const arrow::StringArray&>(*column(TRAML_ID));
        const auto& decoy = static_cast<const arrow::BooleanArray&>(*column(DECOY));
        const auto& precursor_mz = static_cast<const arrow::DoubleArray&>(*column(PRECURSOR_MZ));
        const auto& charge = static_cast<const arrow::Int32Array&>(*column(PRECURSOR_CHARGE));
        const auto& rt = static_cast<const arrow::DoubleArray&>(*column(LIBRARY_RT));
        const auto& im = static_cast<const arrow::DoubleArray&>(*column(PRECURSOR_IM));
        const auto& sequence = static_cast<const arrow::StringArray&>(*column(SEQUENCE));
        const auto& proteins = static_cast<const arrow::ListArray&>(*column(PROTEINS));
        const auto& protein = static_cast<const arrow::StringArray&>(*proteins.values());
        const auto& transition = static_cast<const arrow::StringArray&>(*column(TRANSITION_ID));
        const auto& product_mz = static_cast<const arrow::DoubleArray&>(*column(PRODUCT_MZ));
        const auto& intensity = static_cast<const arrow::DoubleArray&>(*column(LIBRARY_INTENSITY));
        const auto& type = static_cast<const arrow::StringArray&>(*column(FRAGMENT_TYPE));
        const auto& ordinal = static_cast<const arrow::Int16Array&>(*column(FRAGMENT_ORDINAL));
        const auto& fragment_charge = static_cast<const arrow::Int8Array&>(*column(FRAGMENT_CHARGE));
        // Every row of a precursor carries the same precursor fields, so its protein list is decoded once per run of rows.
        std::int32_t proteins_of = -1;
        for (std::int64_t r = 0; r < batch->num_rows(); ++r)
        {
          const std::int32_t p = index.Value(r);
          if (p < 0) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path + ": negative precursor_index", std::to_string(p));
          if (static_cast<std::size_t>(p) >= traml_ids.size())
          {
            traml_ids.resize(p + 1);
            has_traml.resize(p + 1, false);
          }
          if (!has_traml[p])
          {
            traml_ids[p] = {view(id, r), view(traml, r)};
            has_traml[p] = true;
          }
          if (transition.IsNull(r)) continue;
          FragmentView f;
          f.precursor_id = id.IsNull(r) ? std::string_view() : id.GetView(r);
          f.transition_id = transition.GetView(r);
          f.fragment_type = type.IsNull(r) ? std::string_view() : type.GetView(r);
          f.library_intensity = get(intensity, r, 0.0);
          f.product_mz = get(product_mz, r, 0.0);
          f.precursor_mz = get(precursor_mz, r, 0.0);
          f.precursor_im = get(im, r, -1.0);
          f.fragment_charge = get(fragment_charge, r, std::int8_t{0});
          f.fragment_nr = get(ordinal, r, std::int16_t{-1});
          f.decoy = get(decoy, r, false);
          if (p != proteins_of)
          {
            protein_views.clear();
            if (!proteins.IsNull(r))
            {
              for (std::int64_t k = proteins.value_offset(r); k < proteins.value_offset(r + 1); ++k)
                protein_views.push_back(protein.IsNull(k) ? std::string_view() : protein.GetView(k));
            }
            proteins_of = p;
          }
          PrecursorView c;
          c.id = f.precursor_id;
          c.sequence = sequence.IsNull(r) ? std::string_view() : sequence.GetView(r);
          c.rt = get(rt, r, 0.0);
          c.drift_time = f.precursor_im;
          c.charge = get(charge, r, 0);
          c.proteins = &protein_views;
          visitor.fragment(f, c);
        }
      }
    }
    if (std::find(has_traml.begin(), has_traml.end(), false) != has_traml.end())
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path + ": precursor_index has gaps", "");
    // Hand the decoded row-group buffers back to the system; the pool otherwise keeps them resident.
    arrow::default_memory_pool()->ReleaseUnused();
    return traml_ids;
  }

  Library readParquetLibrary(const std::string& path)
  {
    LibraryBuilder builder;
    auto traml_ids = visitParquetLibrary(path, builder);
    builder.library.traml_ids = std::move(traml_ids);
    return std::move(builder.library);
  }

  std::vector<std::pair<std::string, std::string>> visitAnyLibrary(const std::string& path, LibraryVisitor& visitor,
                                                                   const std::function<void(const std::string&)>& log,
                                                                   std::size_t cache_min_transitions)
  {
    const auto say = [&](const std::string& line) { if (log) log(line); };
    if (isParquetLibrary(path)) return visitParquetLibrary(path, visitor);

    const std::string cache = parquetCachePath(path);
    const auto [size, mtime] = sourceKey(path);
    const auto metadata = readParquetLibraryMetadata(cache);
    const auto has = [&](const char* key, const std::string& value)
    {
      const auto it = metadata.find(key);
      return it != metadata.end() && it->second == value;
    };
    if (has(FORMAT_KEY, FORMAT_VALUE) && has(SCHEMA_KEY, SCHEMA_VALUE) && has(READER_KEY, READER_VALUE) && has(SOURCE_SIZE_KEY, size) &&
        has(SOURCE_MTIME_KEY, mtime))
    {
      try
      {
        auto traml_ids = visitParquetLibrary(cache, visitor);
        say("Library: reading " + cache + " (Parquet copy of " + path + ")");
        return traml_ids;
      }
      catch (const std::exception& e)
      {
        visitor.reset();
        say("Library: cannot read " + cache + " (" + e.what() + "); reading the PQP");
      }
    }

    Library lib = readLibrary(path);
    if (lib.transitions.size() >= cache_min_transitions)
    {
      const auto start = std::chrono::steady_clock::now();
      try
      {
        writeParquetLibrary(lib, cache, {{SOURCE_SIZE_KEY, size}, {SOURCE_MTIME_KEY, mtime}, {READER_KEY, READER_VALUE},
                                         {"opendia.source.path", std::filesystem::absolute(path).string()}});
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        say("Library: wrote Parquet copy " + cache + " (" + std::to_string(std::filesystem::file_size(cache) >> 20) + " MiB in " +
            std::to_string(static_cast<long>(seconds + 0.5)) + " s); later runs read it instead of the PQP");
      }
      catch (const std::exception& e)
      {
        say(std::string("Library: no Parquet copy written (") + e.what() + "); reading the PQP each run");
      }
    }
    visitLibrary(lib, visitor);
    return std::move(lib.traml_ids);
  }

  Library readAnyLibrary(const std::string& path, const std::function<void(const std::string&)>& log, std::size_t cache_min_transitions)
  {
    if (!isParquetLibrary(path))
    {
      // Without a usable Parquet copy the PQP is read into a Library anyway: return it instead of rebuilding it.
      const std::string cache = parquetCachePath(path);
      const auto [size, mtime] = sourceKey(path);
      const auto metadata = readParquetLibraryMetadata(cache);
      const auto has = [&](const char* key, const std::string& value)
      {
        const auto it = metadata.find(key);
        return it != metadata.end() && it->second == value;
      };
      if (!(has(FORMAT_KEY, FORMAT_VALUE) && has(SCHEMA_KEY, SCHEMA_VALUE) && has(READER_KEY, READER_VALUE) && has(SOURCE_SIZE_KEY, size) &&
            has(SOURCE_MTIME_KEY, mtime)))
      {
        const auto say = [&](const std::string& line) { if (log) log(line); };
        Library lib = readLibrary(path);
        if (lib.transitions.size() < cache_min_transitions) return lib;
        const auto start = std::chrono::steady_clock::now();
        try
        {
          writeParquetLibrary(lib, cache, {{SOURCE_SIZE_KEY, size}, {SOURCE_MTIME_KEY, mtime}, {READER_KEY, READER_VALUE},
                                           {"opendia.source.path", std::filesystem::absolute(path).string()}});
          const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
          say("Library: wrote Parquet copy " + cache + " (" + std::to_string(std::filesystem::file_size(cache) >> 20) + " MiB in " +
              std::to_string(static_cast<long>(seconds + 0.5)) + " s); later runs read it instead of the PQP");
        }
        catch (const std::exception& e)
        {
          say(std::string("Library: no Parquet copy written (") + e.what() + "); reading the PQP each run");
        }
        return lib;
      }
    }
    LibraryBuilder builder;
    auto traml_ids = visitAnyLibrary(path, builder, log, cache_min_transitions);
    builder.library.traml_ids = std::move(traml_ids);
    return std::move(builder.library);
  }
}
