// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIAWorkflow.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace OpenMS::OpenDIA;
namespace fs = std::filesystem;

namespace
{
void require(bool condition, const std::string& message)
{
  if (! condition) throw std::runtime_error(message);
}

bool same(double a, double b)
{
  return a == b || (std::isnan(a) && std::isnan(b));
}

// Independent, materialized oracle: group by precursor in first-fragment order, then sort the complete
// scoring list by numeric transition ID before truncation. This is deliberately not a bounded insertion list.
void compareLegacy(const AnalysisLibrary& actual, const Library& search, const Library& scoring)
{
  std::vector<std::string> order;
  std::unordered_map<std::string, std::vector<const LibraryTransition*>> search_groups, scoring_groups;
  for (const auto& t : search.transitions)
  {
    if (! search_groups.count(t.peptide_ref)) order.push_back(t.peptide_ref);
    search_groups[t.peptide_ref].push_back(&t);
  }
  for (const auto& t : scoring.transitions)
    scoring_groups[t.peptide_ref].push_back(&t);
  require(actual.precursors.size() == order.size(), "legacy precursor count");
  require(actual.fragment_offset.size() == order.size() + 1 && actual.fragment_offset.front() == 0
            && actual.fragment_offset.back() == search.transitions.size(),
          "legacy CSR extent");
  for (std::size_t i = 0; i < order.size(); ++i)
  {
    const auto& id = order[i];
    const auto& p = actual.precursors[i];
    require(p.id == id && actual.find(id) == i, "first-fragment precursor order and index");
    const auto c = std::find_if(search.compounds.begin(), search.compounds.end(), [&](const auto& x) { return x.id == id; });
    require(c != search.compounds.end(), "oracle missing compound");
    std::string proteins;
    for (const auto& protein : c->protein_refs)
    {
      if (! proteins.empty()) proteins += ';';
      proteins += protein;
    }
    require(p.sequence == c->sequence && p.charge == c->charge && same(p.library_rt, c->rt) && same(p.library_im, c->drift_time)
              && p.has_library_rt == std::isfinite(c->rt) && actual.protein_sets.at(p.proteins) == proteins,
            "legacy precursor metadata");
    const auto& fragments = search_groups.at(id);
    require(p.mz == fragments.front()->precursor_mz && p.decoy == fragments.front()->decoy, "legacy precursor identity");
    require(actual.fragment_offset[i + 1] - actual.fragment_offset[i] == fragments.size(), "legacy CSR group size");
    for (std::size_t f = 0; f < fragments.size(); ++f)
    {
      const auto& got = actual.fragments[actual.fragment_offset[i] + f];
      const auto& expected = *fragments[f];
      require(got.transition == std::stoll(expected.transition_name) && got.product_mz == expected.product_mz
                && got.library_intensity == expected.library_intensity,
              "stable interleaved fragment grouping");
    }
    auto list = scoring_groups[id];
    std::stable_sort(list.begin(), list.end(),
                     [](const auto* a, const auto* b) { return std::stoll(a->transition_name) < std::stoll(b->transition_name); });
    list.resize(std::min<std::size_t>(12, list.size()));
    require(p.scoring_count == list.size(), "legacy top12 count");
    double maximum = 1e-12;
    for (const auto* t : list)
      maximum = std::max(maximum, t->library_intensity);
    for (std::size_t s = 0; s < list.size(); ++s)
    {
      const bool in_search
        = std::any_of(fragments.begin(), fragments.end(), [&](const auto* t) { return t->getAnnotation() == list[s]->getAnnotation(); });
      require(p.scoring_mz[s] == list[s]->product_mz && p.scoring_library[s] == list[s]->library_intensity / maximum
                && p.scoring_in_search[s] == in_search,
              "legacy top12 numeric order, normalization and annotation membership");
    }
    for (std::size_t s = list.size(); s < 12; ++s)
    {
      require(p.scoring_mz[s] == 0.0 && p.scoring_library[s] == 0.0 && ! p.scoring_in_search[s], "unused scoring slots cleared");
    }
  }
  require(actual.find("not-a-precursor") == UINT32_MAX, "unknown precursor lookup");
}

std::uint64_t fnv(const std::string& text)
{
  std::uint64_t hash = 14695981039346656037ULL;
  for (const unsigned char c : text)
    hash = (hash ^ c) * 1099511628211ULL;
  return hash;
}

Library synthetic(bool scoring, unsigned seed)
{
  Library library;
  const std::vector<std::string> ids = {"101", "201", "102", "202", "103", "203", "104"};
  const std::vector<std::string> sequences = {"PEPTIDE", "EDITPEP", "PEPTIDE", "EDITPEP", "PEPT(UniMod:21)IDE", "EDITPEP", "OTHER"};
  const std::vector<std::string> names = {"pep_z2", "DECOY_pep_z2", "pep_z3", "DECOY_pep_z3", "pep_mod", "DECOY_pep_mod", "alone"};
  for (std::size_t i = 0; i < ids.size(); ++i)
  {
    LibraryCompound c;
    c.id = ids[i];
    c.sequence = sequences[i];
    c.charge = i == 2 || i == 3 ? 3 : 2;
    c.rt = i == 6 ? std::nan("") : 10.0 + i;
    c.drift_time = i == 6 ? -1.0 : 0.8 + i * 0.01;
    if (i < 4) c.protein_refs = {"P1", "P2"};
    else if (i < 6)
      c.protein_refs = {"P2", "P1"};
    library.compounds.push_back(c);
    library.traml_ids.emplace_back(ids[i], names[i]);
    const int count = scoring ? (i == 6 ? 0 : 19) : (i == 6 ? 2 : 6);
    for (int f = 0; f < count; ++f)
    {
      LibraryTransition t;
      t.transition_name = std::to_string(1000 * f + i + 1);
      t.peptide_ref = ids[i];
      t.product_mz = 200.0 + 20 * i + f;
      t.precursor_mz = 500.0 + i;
      t.precursor_im = c.drift_time;
      // Maximum intensity is outside the top12. It must not enter the retained-fragment normalization.
      t.library_intensity = f == 18 ? 100000.0 : (i == 5 ? 0.0 : f + 1.0);
      t.fragment_type = f % 2 ? "b" : "y";
      t.fragment_nr = static_cast<std::int16_t>(f + 1);
      t.fragment_charge = f % 3 == 0 ? 2 : 1;
      t.decoy = i % 2 == 1;
      library.transitions.push_back(t);
    }
  }
  // Compounds without fragments must not become compact precursors.
  LibraryCompound empty;
  empty.id = "999";
  library.compounds.push_back(empty);
  library.traml_ids.emplace_back("999", "empty");
  std::mt19937 random(seed);
  std::shuffle(library.transitions.begin(), library.transitions.end(), random);
  // Make traversal order explicitly unlike compound and numeric transition order.
  return library;
}

void checkIdentity(const AnalysisLibrary& library)
{
  require(library.pairs == 4, "dense target-decoy pair count");
  require(library.protein_sets.size() == 3, "protein columns interned including empty and reordered references");
  require(library.calibration_pairs.size() == 3, "calibration target-decoy pairs");
  for (const auto& [target, decoy] : library.calibration_pairs)
  {
    require(! library.precursors[target].decoy && library.precursors[decoy].decoy && library.precursors[target].decoy_partner == decoy,
            "calibration partner indices");
  }
  for (const auto& [target_id, decoy_id, key] :
       std::vector<std::tuple<std::string, std::string, std::string>> {{"101", "201", "pep_z2"}, {"102", "202", "pep_z3"}, {"103", "203", "pep_mod"}})
  {
    const auto target = library.find(target_id), decoy = library.find(decoy_id);
    const auto& t = library.precursors[target];
    const auto& d = library.precursors[decoy];
    require(t.pair == d.pair && t.decoy_partner == decoy && d.decoy_partner == UINT32_MAX, "pair and partner identity");
    for (const auto index : {target, decoy})
    {
      const auto& p = library.precursors[index];
      require(p.fold == fnv("PEPTIDE") % 3 && p.rt_fold == fnv("PEPTIDE") % 5 && p.family_seed == fnv("vPEPTIDE") && p.pair_seed == fnv("v" + key)
                && p.trace_safe,
              "unmodified target family hashes, folds and seeds");
    }
  }
  const auto& a = library.precursors[library.find("101")];
  const auto& b = library.precursors[library.find("102")];
  const auto& mod = library.precursors[library.find("103")];
  require(a.sibling_group == b.sibling_group && a.sibling_group != mod.sibling_group && library.sibling_group_sizes.at(a.sibling_group) == 2
            && library.sibling_group_sizes.at(mod.sibling_group) == 1,
          "siblings share modified sequence, not unmodified family");
  require(library.sibling_group_sizes.at(library.precursors[library.find("201")].sibling_group) == 3, "decoy sibling group size");
  require(a.proteins == b.proteins && a.proteins != mod.proteins, "protein reference order and interning");
  require(library.find("999") == UINT32_MAX, "fragmentless precursor omitted");
}
} // namespace

int main(int argc, char** argv)
{
  if (argc < 2) return 2;
  const fs::path dir = fs::path(argv[1]) / "OpenDIAWorkflow_test";
  fs::remove_all(dir);
  fs::create_directories(dir);
  for (unsigned seed : {7U, 29U, 101U})
  {
    const auto search = synthetic(false, seed), scoring = synthetic(true, seed + 1);
    const auto search_path = (dir / "search.parquet").string(), scoring_path = (dir / "scoring.parquet").string();
    writeParquetLibrary(search, search_path);
    writeParquetLibrary(scoring, scoring_path);
    auto compact = AnalysisLibrary::read(search_path, scoring_path);
    compareLegacy(compact, search, scoring);
    checkIdentity(compact);
    auto search_only = AnalysisLibrary::read(search_path);
    compareLegacy(search_only, search, search);
    checkIdentity(search_only);
    // Move construction and assignment must keep every string_view attached to its owning index node.
    const char* id_address = compact.precursors.front().id.data();
    AnalysisLibrary moved(std::move(compact));
    require(moved.precursors.front().id.data() == id_address, "move construction retains ID allocation");
    search_only = std::move(moved);
    compact = AnalysisLibrary();
    moved = AnalysisLibrary();
    require(search_only.precursors.front().id.data() == id_address, "move assignment retains ID allocation");
    compareLegacy(search_only, search, scoring);
    checkIdentity(search_only);
    for (const auto& p : search_only.precursors)
    {
      require(p.id.data() == search_only.index.find(p.id)->first.data(), "ID views own index after source destruction/reset");
    }
  }
  if (argc > 2)
  {
    const auto pqp = (dir / "fixture.pqp").string();
    fs::copy_file(argv[2], pqp);
    const auto legacy = readLibrary(pqp);
    const auto parquet = (dir / "fixture.parquet").string();
    writeParquetLibrary(legacy, parquet);
    compareLegacy(AnalysisLibrary::read(pqp), legacy, legacy);
    compareLegacy(AnalysisLibrary::read(parquet), legacy, legacy);
    compareLegacy(AnalysisLibrary::read(pqp, parquet), legacy, legacy);
    // Obtain genuine cache metadata, then make the last transition fail after earlier fragments were accepted.
    // AnalysisLibrary must reset both search and scoring visitor state before retrying from the original PQP.
    readAnyLibrary(pqp, {}, 0);
    const auto cache = parquetCachePath(pqp);
    const auto metadata = readParquetLibraryMetadata(cache);
    // A Codon experiment can link immutable baseline PQPs into its own directory. Automatic conversion must write
    // beside the link, not beside the resolved baseline source. Symlink creation needs extra privileges on Windows.
    const auto linked = (dir / "linked.pqp").string();
    std::error_code symlink_error;
    fs::create_symlink(pqp, linked, symlink_error);
    if (!symlink_error)
    {
      const auto original_cache_time = fs::last_write_time(cache);
      const auto linked_library = readAnyLibrary(linked, {}, 0);
      require(fs::exists(parquetCachePath(linked)), "automatic cache must be adjacent to the PQP symlink");
      require(fs::last_write_time(cache) == original_cache_time, "symlink conversion must not rewrite the baseline cache");
      require(fs::equivalent(pqp, linked) && linked_library.transitions.size() == legacy.transitions.size(), "linked PQP parity");
      compareLegacy(AnalysisLibrary::read(linked), legacy, legacy);
    }
    std::vector<std::pair<std::string, std::string>> entries(metadata.begin(), metadata.end());
    auto damaged = legacy;
    require(damaged.transitions.size() > 1, "PQP fixture needs multiple transitions for partial fallback");
    damaged.transitions.back().transition_name = "not-an-integer";
    for (const bool scoring_failure : {false, true})
    {
      writeParquetLibrary(damaged, cache, entries);
      std::vector<std::string> messages;
      const auto log = [&](const std::string& line) { messages.push_back(line); };
      auto compact = scoring_failure ? AnalysisLibrary::read(parquet, pqp, log) : AnalysisLibrary::read(pqp, "", log);
      compareLegacy(compact, legacy, legacy);
      require(std::any_of(messages.begin(), messages.end(), [](const auto& line) { return line.find("cannot read") != std::string::npos; }),
              "partial visitor failure did not take PQP fallback");
    }
  }
  fs::remove_all(dir);
  std::cout << "OpenDIAWorkflow_test passed\n";
  return 0;
}
