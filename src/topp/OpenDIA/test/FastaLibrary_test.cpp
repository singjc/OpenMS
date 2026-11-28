// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// FASTA library (PLAN_FASTA): candidate enumeration, the shared decoy rule, fragment rules, the calibration panel,
// pair selection and, with WITH_ONNX, the survey and final library round trip.
#include <library/Library.h>
#include <library/preparation/DecoySequences.h>
#include <library/preparation/FastaLibrary.h>
#include <workflow/Workflow.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/config.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using namespace OpenMS;
namespace P = OpenMS::OpenDIA::Preparation;
namespace fs = std::filesystem;

namespace
{
  void require(bool condition, const std::string& message)
  {
    if (!condition) throw std::runtime_error(message);
  }

  bool near(double a, double b, double tolerance = 1e-6) { return std::fabs(a - b) <= tolerance; }

  std::string writeFasta(const fs::path& dir)
  {
    const auto path = dir / "test.fasta";
    std::ofstream out(path);
    // P1 and P2 share SAMPLEPEPTIDEK; P2 has a C (fixed mod) and an M (variable mod); P3 has a non-standard residue.
    out << ">sp|P00001|ONE_HUMAN first\nMKSAMPLEPEPTIDEKAGHVLLCDERSTWYK\n";
    out << ">sp|P00002|TWO_HUMAN second\nSAMPLEPEPTIDEKQQMVNPGTLLRAC\n";
    out << ">tr|Q00003|THREE_HUMAN third\nGGGXGGGGGKAAAAAAAAR\n";
    return path.string();
  }

  void decoyRuleTests()
  {
    // The shared rule: reverse all but the last residue, modifications move along, K <-> R.
    std::string s = "PEPCTIDEK";
    std::vector<P::Modification> mods{{3, 4}};
    require(P::pseudoReverse(s, mods), "decoy for a residue modification");
    require(s == "EDITCPEPR", "pseudo-reverse sequence: " + s);
    require(mods.size() == 1 && mods[0].location == 4, "modification follows its residue");
    require(P::uniModSequence(s, mods) == "EDITC(UniMod:4)PEPR", "UniMod notation");
    std::string t = "PEPTIDEK";
    std::vector<P::Modification> last{{7, 259}}, nterm{{-1, 1}}, cterm{{8, 2}};
    require(!P::pseudoReverse(t, last) && t == "PEPTIDEK", "last-residue modification has no decoy, input unchanged");
    require(!P::pseudoReverse(t, nterm), "N-terminal modification has no decoy");
    require(!P::pseudoReverse(t, cterm), "C-terminal modification has no decoy");
    std::string other = "ACDEFGHI";
    std::vector<P::Modification> none;
    require(P::pseudoReverse(other, none) && other.substr(0, 7) == "HGFEDCA", "non-tryptic C terminus is replaced");

    // decoySequences uses the same rule (no transitions: sequence stage only).
    P::Assays targets;
    P::Compound c;
    c.id = "t";
    c.sequence = "PEPC(UniMod:4)TIDEK";
    c.charge = 2;
    c.modifications = {{3, 4}};
    targets.compounds.push_back(c);
    const auto decoys = P::decoySequences(targets);
    require(decoys.compounds.size() == 1 && decoys.compounds[0].sequence == "EDITC(UniMod:4)PEPR", "decoySequences uses the shared rule");
  }

  void ladderTests()
  {
    const auto peptide = AASequence::fromString("PEPC(Carbamidomethyl)TIDEK");
    const P::IonLadder ladder(peptide);
    for (int z = 1; z <= 3; ++z) require(near(ladder.precursor(z), peptide.getMZ(z)), "precursor m/z");
    for (std::size_t o = 1; o < peptide.size(); ++o)
    {
      for (int z = 1; z <= 2; ++z)
      {
        require(near(ladder.ion('b', o, z), peptide.getPrefix(o).getMZ(z, Residue::BIon)), "b ion m/z");
        require(near(ladder.ion('y', o, z), peptide.getSuffix(o).getMZ(z, Residue::YIon)), "y ion m/z");
      }
    }
  }

  void fragmentRuleTests()
  {
    // Survey fragments: y1+ ordinal 3..n-1, 6 closest to 7 (lower first), then b1+ closest to 5.
    const P::IonLadder long_ladder(AASequence::fromString("AAAAAAAAAAAAAAAAAAAK"));  // n = 20
    const auto survey = P::surveyFragments(long_ladder, 0.0, 1e9);
    require(survey.size() == 6, "six survey fragments");
    const std::vector<int> expected{7, 6, 8, 5, 9, 4};
    for (std::size_t i = 0; i < 6; ++i) require(survey[i].type == 'y' && survey[i].ordinal == expected[i], "survey y order");
    const P::IonLadder short_ladder(AASequence::fromString("AAAAAK"));  // n = 6: y3..y5 only, then b5, b4, b3
    const auto top_up = P::surveyFragments(short_ladder, 0.0, 1e9);
    require(top_up.size() == 6, "survey top-up to six");
    require(top_up[0].type == 'y' && top_up[2].type == 'y' && top_up[3].type == 'b' && top_up[3].ordinal == 5 && top_up[5].ordinal == 3,
            "survey b top-up order");
    // m/z range applies.
    const auto ranged = P::surveyFragments(long_ladder, 0.0, long_ladder.ion('y', 6, 1) + 0.01);
    require(std::all_of(ranged.begin(), ranged.end(), [&](const auto& f) { return f.mz <= long_ladder.ion('y', 6, 1) + 0.01; }), "survey m/z range");

    // Predicted fragments: 8 channels per cleavage; z <= max(1, min(2, z_prec - 1)); ordinal >= 3; top N by intensity.
    const std::size_t n = 10;
    const P::IonLadder ladder(AASequence::fromString("AAAAAAAAAK"));
    std::vector<float> ms2((n - 1) * 8, 0.0f);
    for (std::size_t row = 0; row + 1 < n; ++row)
    {
      ms2[row * 8 + 0] = 0.1f * static_cast<float>(row + 1);  // b z1
      ms2[row * 8 + 1] = 0.9f;                                 // b z2
      ms2[row * 8 + 2] = 0.2f * static_cast<float>(row + 1);  // y z1
      ms2[row * 8 + 4] = 5.0f;                                 // modloss: never used
    }
    const auto z2 = P::predictedFragments(ladder, 2, ms2, 0.0, 1e9, 99);
    require(std::all_of(z2.begin(), z2.end(), [](const auto& f) { return f.charge == 1 && f.ordinal >= 3 && f.intensity > 0.0f; }),
            "charge-2 precursor: singly charged fragments of ordinal >= 3");
    require(z2.size() == 2 * 7, "b3..b9 and y3..y9 at z = 1");
    require(std::is_sorted(z2.begin(), z2.end(), [](const auto& a, const auto& b) { return a.intensity > b.intensity; }), "intensity order");
    const auto z3 = P::predictedFragments(ladder, 3, ms2, 0.0, 1e9, 99);
    require(std::any_of(z3.begin(), z3.end(), [](const auto& f) { return f.charge == 2; }), "charge-3 precursor admits z = 2 fragments");
    require(P::predictedFragments(ladder, 3, ms2, 0.0, 1e9, 4).size() == 4, "top-N cut");
    bool threw = false;
    try { P::predictedFragments(ladder, 2, std::vector<float>(3), 0.0, 1e9, 6); }
    catch (const Exception::InvalidValue&) { threw = true; }
    require(threw, "wrong MS2 size is rejected");
  }

  void candidateTests(const fs::path& dir)
  {
    const auto fasta = writeFasta(dir);
    P::CandidateOptions o;
    o.min_length = 6;
    o.min_charge = 2;
    o.max_charge = 3;
    o.min_mz = 300.0;
    o.max_mz = 1500.0;
    const auto set = P::enumerateCandidates(fasta, o);
    require(set.proteins_read == 3, "proteins read");
    require(set.skipped_residues >= 1, "non-standard residue skipped");
    require(set.precursors.size() % 2 == 0 && set.pairs() > 0, "complete pairs");
    std::set<std::string> targets, decoys;
    for (std::size_t k = 0; k < set.pairs(); ++k)
    {
      const auto& t = set.precursors[2 * k];
      const auto& d = set.precursors[2 * k + 1];
      const auto& tf = set.peptidoforms[t.peptidoform];
      const auto& df = set.peptidoforms[d.peptidoform];
      require(!tf.decoy && df.decoy && t.charge == d.charge, "pair layout: target, decoy, same charge");
      require(t.mz >= o.min_mz && t.mz <= o.max_mz, "target m/z in range");
      require(set.tramlId(2 * k + 1) == "DECOY_" + set.tramlId(2 * k), "TRAML_IDs pair by name");
      // The decoy is the shared rule applied to the target.
      const auto parsed = AASequence::fromString(tf.sequence);
      std::string reversed = parsed.toUnmodifiedString();
      std::vector<P::Modification> mods;
      for (Size i = 0; i < parsed.size(); ++i)
        if (parsed[i].isModified()) mods.push_back({static_cast<int>(i), parsed[i].getModification()->getUniModRecordId()});
      require(P::pseudoReverse(reversed, mods) && P::uniModSequence(reversed, mods) == df.sequence, "decoy is the shared rule");
      targets.insert(tf.sequence);
      decoys.insert(df.sequence);
      require(tf.sequence.find("C") == std::string::npos || tf.sequence.find("C(UniMod:4)") != std::string::npos, "fixed carbamidomethyl");
      for (const auto& p : set.protein_sets[df.proteins]) require(p.rfind("DECOY_", 0) == 0, "decoy proteins are DECOY_");
    }
    for (const auto& d : decoys) require(!targets.count(d), "no decoy equals a target");
    // Shared peptide lists both proteins; accessions are parsed from sp|...| headers.
    bool shared = false, oxidised = false;
    for (std::size_t k = 0; k < set.pairs(); ++k)
    {
      const auto& f = set.peptidoforms[set.precursors[2 * k].peptidoform];
      if (f.sequence == "SAMPLEPEPTIDEK")
      {
        const auto& proteins = set.protein_sets[f.proteins];
        shared = proteins == std::vector<std::string>{"P00001", "P00002"};
        require(f.missed_cleavages == 0, "no missed cleavage");
      }
      oxidised = oxidised || f.sequence.find("M(UniMod:35)") != std::string::npos;
    }
    require(shared, "shared peptide maps to both proteins");
    require(oxidised, "variable oxidation enumerated");
    require(P::missedCleavages("Trypsin/P", "AAKAAR") == 1 && P::missedCleavages("Trypsin/P", "AAAAR") == 0, "missed cleavages");

    // Terminal modifications are rejected up front.
    P::CandidateOptions bad = o;
    bad.variable_mods = {"Acetyl (N-term)"};
    bool threw = false;
    try { P::enumerateCandidates(fasta, bad); }
    catch (const Exception::InvalidValue&) { threw = true; }
    require(threw, "terminal modification rejected");
  }

  void panelAndSelectionTests()
  {
    P::CandidateSet set;
    // 10 pairs: even k eligible (charge 2, unmodified, no missed cleavage, length 10); RT descending in k.
    for (std::uint32_t k = 0; k < 10; ++k)
    {
      P::Peptidoform t, d;
      t.sequence = "AAAAAAAAAK";
      t.length = 10;
      t.rt = 100.0f - static_cast<float>(k);
      t.variable_mods = k % 2;
      d = t;
      d.decoy = true;
      set.peptidoforms.push_back(t);
      set.peptidoforms.push_back(d);
      set.precursors.push_back({2 * k, 500.0, -1.0f, 2});
      set.precursors.push_back({2 * k + 1, 500.0, -1.0f, 2});
    }
    const auto panel = P::calibrationPanel(set, 3);
    require(panel.size() == 3, "panel size");
    // Eligible sorted by RT: 8, 6, 4, 2, 0 -> stride picks i * 5 / 3 = 0, 1, 3 -> 8, 6, 2.
    require(panel == std::vector<std::uint32_t>({8, 6, 2}), "panel at an even RT stride");
    require(P::calibrationPanel(set, 100).size() == 5, "panel capped at the eligible pairs");

    const std::vector<double> s{0.0, 2.0, 1.0, 2.0, 3.0};
    const std::vector<std::uint64_t> tie{0, 7, 0, 5, 0};
    require(P::selectPairs(s, tie, 0) == std::vector<std::uint32_t>({0, 1, 2, 3, 4}), "K = 0 keeps every pair");
    require(P::selectPairs(s, tie, 2) == std::vector<std::uint32_t>({3, 4}), "top K, ties by hash");
    require(P::selectPairs(s, tie, 10) == std::vector<std::uint32_t>({1, 2, 3, 4}), "S = 0 pairs are never kept");
  }

  void writerTests(const fs::path& dir)
  {
    const auto path = (dir / "writer.parquet").string();
    {
      OpenDIA::ParquetLibraryWriter writer(path, {{"test.key", "1"}});
      const std::vector<std::string> proteins{"P1"}, decoy_proteins{"DECOY_P1"};
      OpenDIA::ParquetLibraryWriter::Precursor p;
      p.id = "0";
      p.traml_id = "PEPTIDEK/2";
      p.sequence = "PEPTIDEK";
      p.mz = 465.0;
      p.rt = 1234.5;
      p.charge = 2;
      p.proteins = &proteins;
      writer.add(p, {{"0", 300.0, 0.5, "y", 3, 1}, {"1", 400.0, 1.0, "y", 4, 1}, {"2", 500.0, 0.25, "b", 4, 1}});
      p.id = "1";
      p.traml_id = "DECOY_PEPTIDEK/2";
      p.sequence = "EDITPEPR";
      p.decoy = true;
      p.proteins = &decoy_proteins;
      writer.add(p, {{"100", 310.0, 1.0, "y", 3, 1}, {"101", 410.0, 0.5, "y", 4, 1}, {"102", 510.0, 0.25, "b", 4, 1}});
      require(!fs::exists(path), "file appears only on close");
      writer.close();
    }
    const auto lib = OpenDIA::readParquetLibrary(path);
    require(lib.compounds.size() == 2 && lib.transitions.size() == 6, "writer round trip sizes");
    require(lib.traml_ids[1].second == "DECOY_PEPTIDEK/2", "TRAML_ID round trip");
    require(lib.transitions[3].decoy && lib.transitions[3].peptide_ref == "1", "decoy flag round trip");
    require(OpenDIA::readParquetLibraryMetadata(path).at("test.key") == "1", "footer metadata");
    const auto analysis = OpenDIA::AnalysisLibrary::read(path);
    require(analysis.pairs == 1 && analysis.calibration_pairs.size() == 1, "analysis library pairs the written library");
    {
      OpenDIA::ParquetLibraryWriter abandoned((dir / "abandoned.parquet").string());
    }
    for (const auto& e : fs::directory_iterator(dir))
      require(e.path().string().find(".tmp") == std::string::npos, "abandoned writer leaves no temporary file");
  }

#ifdef WITH_ONNX
  void predictionTests(const fs::path& dir)
  {
    fs::create_directories(dir);
    const auto fasta = writeFasta(dir);
    P::CandidateOptions o;
    o.min_length = 6;
    o.min_mz = 300.0;
    o.max_mz = 1500.0;
    auto set = P::enumerateCandidates(fasta, o);
    P::PredictionSettings settings;
    settings.predict_im = true;
    settings.instrument = "timsTOF";
    settings.threads = 2;
    P::predictCoordinates(set, settings);
    for (const auto& p : set.precursors) require(p.im > 0.3f && p.im < 2.0f, "predicted 1/K0");
    for (const auto& f : set.peptidoforms) require(std::isfinite(f.rt), "predicted RT");
    P::FragmentOptions fragments;
    const auto panel = P::calibrationPanel(set, 2);
    const auto survey = (dir / "survey.parquet").string();
    P::writeSurveyLibrary(set, panel, settings, fragments, survey, {});
    const std::size_t lead = P::surveyPanelPrecursors(survey);
    require(lead % 2 == 0 && lead <= 2 * panel.size(), "panel precursors lead the survey library");
    const auto survey_lib = OpenDIA::readParquetLibrary(survey);
    std::map<std::string, std::size_t> per_precursor;
    for (const auto& t : survey_lib.transitions) ++per_precursor[t.peptide_ref];
    for (const auto& [id, count] : per_precursor) require(count <= 6, "at most 6 survey fragments");

    const auto library = (dir / "final.parquet").string(), search = (dir / "final.search.parquet").string();
    const auto stats = P::writeFinalLibrary(survey, [](std::size_t k) { return k % 2 == 0; }, settings, fragments, library, search, {});
    require(stats.pairs_out > 0 && stats.pairs_out <= stats.pairs_in, "final library pairs");
    const auto final_lib = OpenDIA::readParquetLibrary(library);
    const auto search_lib = OpenDIA::readParquetLibrary(search);
    require(final_lib.compounds.size() == 2 * stats.pairs_out && search_lib.compounds.size() == final_lib.compounds.size(), "same precursors");
    std::map<std::string, std::size_t> fragments_of;
    for (const auto& t : final_lib.transitions)
    {
      ++fragments_of[t.peptide_ref];
      require(std::stoull(t.peptide_ref) / 2 % 2 == 0, "only kept pairs");
      require(t.fragment_nr >= 3 && t.product_mz >= fragments.min_mz && t.product_mz <= fragments.max_mz, "fragment rules");
    }
    for (const auto& [id, count] : fragments_of) require(count >= 3 && count <= 12, "3 to 12 fragments");
    const auto analysis = OpenDIA::AnalysisLibrary::read(search, library);
    require(analysis.pairs == stats.pairs_out && analysis.calibration_pairs.size() == stats.pairs_out, "analysis pairs");
  }
#endif
}

int main(int argc, char** argv)
{
  try
  {
    require(argc >= 2, "usage: FastaLibrary_test SCRATCH");
    const fs::path dir(argv[1]);
    fs::remove_all(dir);
    fs::create_directories(dir);
    decoyRuleTests();
    ladderTests();
    fragmentRuleTests();
    candidateTests(dir);
    panelAndSelectionTests();
    writerTests(dir);
#ifdef WITH_ONNX
    predictionTests(dir / "prediction");
#endif
    std::cout << "FASTA library tests passed\n";
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
