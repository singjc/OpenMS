// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
// OpenDIA Parquet library: PQP round trip, the automatic Parquet copy of a PQP (reuse, staleness, size threshold),
// precursors without fragments, repeated precursor IDs and rejected inputs. Writes only below argv[1].
#include <library/Library.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace OpenMS;
using namespace OpenMS::OpenDIA;
namespace fs = std::filesystem;

namespace
{
  void require(bool value, const std::string& message)
  {
    if (!value) throw std::runtime_error(message);
  }

  bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

  void requireEqual(const Library& a, const Library& b, const std::string& what)
  {
    require(a.transitions.size() == b.transitions.size(), what + ": transition count");
    require(a.compounds.size() == b.compounds.size(), what + ": compound count");
    require(a.traml_ids == b.traml_ids, what + ": TRAML_IDs");
    for (std::size_t i = 0; i < a.transitions.size(); ++i)
    {
      const auto& x = a.transitions[i];
      const auto& y = b.transitions[i];
      require(x.transition_name == y.transition_name && x.peptide_ref == y.peptide_ref && same(x.library_intensity, y.library_intensity) &&
                same(x.product_mz, y.product_mz) && same(x.precursor_mz, y.precursor_mz) && same(x.precursor_im, y.precursor_im) &&
                x.fragment_charge == y.fragment_charge && x.fragment_nr == y.fragment_nr && x.fragment_type == y.fragment_type &&
                x.decoy == y.decoy,
              what + ": transition " + std::to_string(i));
    }
    for (std::size_t i = 0; i < a.compounds.size(); ++i)
    {
      const auto& x = a.compounds[i];
      const auto& y = b.compounds[i];
      require(x.id == y.id && same(x.rt, y.rt) && same(x.drift_time, y.drift_time) && x.charge == y.charge && x.sequence == y.sequence &&
                x.protein_refs == y.protein_refs,
              what + ": compound " + std::to_string(i));
    }
  }

  template <typename Action> void rejects(Action action, const std::string& what)
  {
    bool rejected = false;
    try
    {
      action();
    }
    catch (const Exception::BaseException&)
    {
      rejected = true;
    }
    require(rejected, what + " was accepted");
  }

  LibraryTransition fragment(const std::string& id, const std::string& precursor, double mz, double intensity, bool decoy = false)
  {
    LibraryTransition t;
    t.transition_name = id;
    t.peptide_ref = precursor;
    t.product_mz = mz;
    t.library_intensity = intensity;
    t.precursor_mz = 500.25;
    t.precursor_im = -1.0;
    t.fragment_type = "y";
    t.fragment_nr = 3;
    t.fragment_charge = 1;
    t.decoy = decoy;
    return t;
  }

  /// Two precursors with fragments (one decoy, two proteins, a NaN RT), one without fragments, a repeated PRECURSOR ID.
  Library synthetic()
  {
    Library lib;
    LibraryCompound a;
    a.id = "1";
    a.rt = 12.5;
    a.charge = 2;
    a.sequence = "PEPT(UniMod:21)IDE";
    a.protein_refs = {"P1", "P2"};
    LibraryCompound d;
    d.id = "2";
    d.rt = std::nan("");
    d.drift_time = 0.9;
    d.charge = 3;
    d.sequence = "EDITPEP";
    lib.compounds = {a, d};
    lib.transitions = {fragment("10", "1", 300.1, 5.0), fragment("11", "1", 400.2, 9.0), fragment("20", "2", 350.0, 1.0, true),
                       fragment("12", "1", 410.0, 9.0)};
    lib.transitions[2].precursor_im = 0.9;
    lib.transitions[2].fragment_type = "";
    lib.transitions[2].fragment_nr = -1;
    lib.traml_ids = {{"1", "pep_1"}, {"3", "pep_without_fragments"}, {"2", "DECOY_pep_1"}, {"1", "pep_1_again"}};
    return lib;
  }
}

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cerr << "usage: OpenDIALibraryParquet_test <scratch dir> <PQP fixture>\n";
    return 2;
  }
  const fs::path dir = fs::path(argv[1]) / "OpenDIALibraryParquet_test";
  fs::remove_all(dir);
  fs::create_directories(dir);

  // Synthetic library: everything survives, including precursors without fragments and a repeated precursor ID.
  {
    const Library lib = synthetic();
    const std::string path = (dir / "synthetic.parquet").string();
    writeParquetLibrary(lib, path, {{"note", "synthetic"}});
    require(fs::exists(path), "synthetic library not written");
    require(isParquetLibrary(path) && isParquetLibrary("X.PARQUET") && !isParquetLibrary("x.pqp"), "isParquetLibrary");
    for (const auto& entry : fs::directory_iterator(dir)) require(entry.path().extension() != ".tmp", "temporary file left behind");
    const Library back = readParquetLibrary(path);
    requireEqual(lib, back, "synthetic");
    // The TRAML map keeps the first entry of a repeated precursor ID, like readPrecursorTraMLIDs.
    const auto map = traMLIDMap(back.traml_ids);
    require(map.size() == 3 && map.at("1") == "pep_1" && map.at("3") == "pep_without_fragments", "TRAML map");
    const auto metadata = readParquetLibraryMetadata(path);
    require(metadata.at("note") == "synthetic" && metadata.at("opendia.library.schema_version") == "1" &&
              metadata.at("opendia.library.fragments") == "4",
            "metadata");
    require(readParquetLibraryMetadata((dir / "missing.parquet").string()).empty(), "metadata of a missing file");
    // readAnyLibrary reads a .parquet directly.
    requireEqual(lib, readAnyLibrary(path), "readAnyLibrary(.parquet)");
  }

  // Rejected inputs: a fragment of an unknown precursor, a non-OpenDIA Parquet file, a non-Parquet file.
  {
    Library bad = synthetic();
    bad.transitions[0].peptide_ref = "unknown";
    const std::string path = (dir / "bad.parquet").string();
    rejects([&] { writeParquetLibrary(bad, path); }, "fragment of an unknown precursor");
    require(!fs::exists(path), "rejected library was written");
    const std::string text = (dir / "text.parquet").string();
    std::ofstream(text) << "not parquet";
    rejects([&] { readParquetLibrary(text); }, "a text file");
  }

  // PQP fixture: the Parquet copy holds exactly what the PQP reader returns.
  const fs::path pqp = dir / "library.pqp";
  fs::copy_file(argv[2], pqp);
  const Library reference = readLibrary(pqp.string());
  require(!reference.transitions.empty(), "empty fixture");
  {
    const auto map = readPrecursorTraMLIDs(pqp.string());
    const auto rebuilt = traMLIDMap(reference.traml_ids);
    require(map == rebuilt, "TRAML map differs from readPrecursorTraMLIDs");
    auto a = map.begin();
    auto b = rebuilt.begin();
    for (; a != map.end(); ++a, ++b) require(a->first == b->first, "TRAML map iteration order differs");
  }
  const std::string cache = parquetCachePath(pqp.string());
  require(cache == (dir / "library.opendia.parquet").string(), "cache path");

  // Below the size threshold nothing is written.
  std::vector<std::string> log;
  const auto record = [&](const std::string& line) { log.push_back(line); };
  requireEqual(reference, readAnyLibrary(pqp.string(), record), "small PQP");
  require(!fs::exists(cache), "Parquet copy of a small library written");

  // First read writes the copy, second read uses it.
  log.clear();
  requireEqual(reference, readAnyLibrary(pqp.string(), record, 0), "PQP, writing the copy");
  require(fs::exists(cache) && log.size() == 1 && log[0].find("wrote Parquet copy") != std::string::npos, "copy not written");
  requireEqual(reference, readParquetLibrary(cache), "Parquet copy");
  const auto written = fs::last_write_time(cache);
  log.clear();
  requireEqual(reference, readAnyLibrary(pqp.string(), record, 0), "PQP, reusing the copy");
  require(log.size() == 1 && log[0].find("reading " + cache) != std::string::npos && fs::last_write_time(cache) == written, "copy not reused");

  // A changed PQP (new modification time) makes the copy stale: it is rewritten.
  fs::last_write_time(pqp, fs::last_write_time(pqp) + std::chrono::seconds(5));
  log.clear();
  requireEqual(reference, readAnyLibrary(pqp.string(), record, 0), "touched PQP");
  require(log.size() == 1 && log[0].find("wrote Parquet copy") != std::string::npos, "stale copy reused");

  // A damaged copy is ignored and replaced.
  std::ofstream(cache, std::ios::trunc) << "damaged";
  log.clear();
  requireEqual(reference, readAnyLibrary(pqp.string(), record, 0), "damaged copy");
  require(log.size() == 1 && log[0].find("wrote Parquet copy") != std::string::npos, "damaged copy not replaced");
  requireEqual(reference, readParquetLibrary(cache), "replaced copy");

  // An unwritable directory: the PQP is read, no copy is made.
  {
    const fs::path locked = dir / "locked";
    fs::create_directories(locked);
    fs::copy_file(pqp, locked / "library.pqp");
    fs::permissions(locked, fs::perms::owner_read | fs::perms::owner_exec);
    const bool writable = std::ofstream((locked / "probe").string()).good();
    if (!writable)
    {
      log.clear();
      requireEqual(reference, readAnyLibrary((locked / "library.pqp").string(), record, 0), "read-only directory");
      require(log.size() == 1 && log[0].find("no Parquet copy written") != std::string::npos, "read-only directory not reported");
      require(!fs::exists(locked / "library.opendia.parquet"), "copy in a read-only directory");
    }
    fs::permissions(locked, fs::perms::owner_all);
  }

  fs::remove_all(dir);
  std::cout << "OpenDIALibraryParquet_test passed\n";
  return 0;
}
