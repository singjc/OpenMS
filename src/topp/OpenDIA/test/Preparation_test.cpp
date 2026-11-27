// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
#include <library/preparation/AssayPrediction.h>
#include <library/preparation/DecoySequences.h>
#include <library/preparation/Preparation.h>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryPreparation.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/config.h>
#ifdef WITH_ONNX
#include <OpenMS/ANALYSIS/OPENSWATH/PeptDeepLibraryPredictor.h>
#include <sqlite3.h>
#endif
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <memory>

using namespace OpenMS;
namespace P = OpenMS::OpenDIA::Preparation;

namespace
{
  void require(bool condition, const std::string& message)
  {
    if (!condition) throw std::runtime_error("Preparation parity: " + message);
  }

  // Explicit field lists deliberately do not use the light records' ID-only equality operators.
#define COMPOUND_FIELDS(X) X(drift_time) X(rt) X(rt_start) X(rt_end) X(charge) X(sequence) X(protein_refs) X(peptide_group_label) X(gene_name) X(id) X(sum_formula) X(compound_name) X(label_type) X(smiles) X(adducts)
#define TRANSITION_FIELDS(X) X(transition_name) X(peptide_ref) X(library_intensity) X(product_mz) X(precursor_mz) X(precursor_im) X(fragment_charge) X(fragment_nr) X(peptidoforms)
#define PROTEIN_FIELDS(X) X(id) X(sequence) X(uniprot_id)

  P::Assays owned(const OpenSwath::LightTargetedExperiment& input)
  {
    P::Assays result;
    for (const auto& s : input.compounds)
    {
      P::Compound d;
#define COPY(field) d.field = s.field;
      COMPOUND_FIELDS(COPY)
      for (const auto& m : s.modifications) d.modifications.push_back({m.location, m.unimod_id});
      result.compounds.push_back(std::move(d));
    }
    for (const auto& s : input.transitions)
    {
      P::Transition d;
      TRANSITION_FIELDS(COPY)
      d.fragment_type = s.getFragmentType();
      d.decoy = s.getDecoy();
      d.detecting = s.isDetectingTransition();
      d.quantifying = s.isQuantifyingTransition();
      d.identifying = s.isIdentifyingTransition();
      result.transitions.push_back(std::move(d));
    }
    for (const auto& s : input.proteins)
    {
      P::Protein d;
      PROTEIN_FIELDS(COPY)
      result.proteins.push_back(std::move(d));
    }
#undef COPY
    return result;
  }

  OpenSwath::LightTargetedExperiment light(const P::Assays& input)
  {
    OpenSwath::LightTargetedExperiment result;
    for (const auto& s : input.compounds)
    {
      OpenSwath::LightCompound d;
#define COPY(field) d.field = s.field;
      COMPOUND_FIELDS(COPY)
      for (const auto& m : s.modifications) d.modifications.push_back({m.location, m.unimod_id});
      result.compounds.push_back(std::move(d));
    }
    for (const auto& s : input.transitions)
    {
      OpenSwath::LightTransition d;
      TRANSITION_FIELDS(COPY)
      d.setFragmentType(s.fragment_type);
      d.setDecoy(s.decoy);
      d.setDetectingTransition(s.detecting);
      d.setQuantifyingTransition(s.quantifying);
      d.setIdentifyingTransition(s.identifying);
      result.transitions.push_back(std::move(d));
    }
    for (const auto& s : input.proteins)
    {
      OpenSwath::LightProtein d;
      PROTEIN_FIELDS(COPY)
      result.proteins.push_back(std::move(d));
    }
#undef COPY
    return result;
  }

  template <typename T> bool same(const T& a, const T& b) { return a == b; }
  bool same(double a, double b) { return a == b || (std::isnan(a) && std::isnan(b)); }

  void equal(const P::Assays& a, const P::Assays& b, const std::string& context)
  {
    require(a.compounds.size() == b.compounds.size(), context + " compound count");
    require(a.transitions.size() == b.transitions.size(), context + " transition count");
    require(a.proteins.size() == b.proteins.size(), context + " protein count");
#define CHECK(field) require(same(x.field, y.field), context + " " + std::to_string(i) + " " #field);
    for (size_t i = 0; i < a.compounds.size(); ++i)
    {
      const auto& x = a.compounds[i];
      const auto& y = b.compounds[i];
      COMPOUND_FIELDS(CHECK)
      CHECK(modifications)
    }
    for (size_t i = 0; i < a.transitions.size(); ++i)
    {
      const auto& x = a.transitions[i];
      const auto& y = b.transitions[i];
      TRANSITION_FIELDS(CHECK)
      CHECK(fragment_type) CHECK(decoy) CHECK(detecting) CHECK(quantifying) CHECK(identifying)
    }
    for (size_t i = 0; i < a.proteins.size(); ++i)
    {
      const auto& x = a.proteins[i];
      const auto& y = b.proteins[i];
      PROTEIN_FIELDS(CHECK)
    }
#undef CHECK
  }
#undef COMPOUND_FIELDS
#undef TRANSITION_FIELDS
#undef PROTEIN_FIELDS

  P::Assays oracle(const P::Assays& targets)
  {
    MRMDecoy generator;
    generator.setLogType(ProgressLogger::NONE);
    OpenSwath::LightTargetedExperiment decoys;
    generator.generateDecoysLight(light(targets), decoys, "pseudo-reverse", 1.0, true, "DECOY_",
                                  30, 0.5, 0.0, 20.0, 0.025, {"b", "y"}, {1, 2, 3, 4}, false, false, -4, 1);
    decoys.transitions.clear();
    return owned(decoys);
  }

  void add(P::Assays& assays, const std::string& id, const std::string& sequence, int charge = 2,
           std::vector<P::Modification> mods = {})
  {
    P::Compound c;
    c.id = id;
    c.sequence = sequence;
    c.charge = charge;
    c.modifications = std::move(mods);
    c.protein_refs = {"P2", "P1"};
    c.rt = 123.25;
    c.rt_start = 120.0;
    c.rt_end = 130.0;
    c.drift_time = 1.125;
    c.peptide_group_label = "group_" + id;
    c.gene_name = "GENE";
    c.label_type = "light";
    assays.compounds.push_back(c);
    const auto peptide = AASequence::fromString(sequence);
    P::Transition t;
    t.transition_name = id + "_y3";
    t.peptide_ref = id;
    t.precursor_mz = peptide.getMZ(charge);
    t.product_mz = peptide.getSuffix(3).getMonoWeight(Residue::YIon, 1);
    t.library_intensity = 42.5;
    t.precursor_im = c.drift_time;
    t.fragment_type = "y";
    t.fragment_nr = 3;
    t.fragment_charge = 1;
    t.peptidoforms = {sequence};
    t.identifying = true;
    t.quantifying = false;
    assays.transitions.push_back(t);
  }

  bool contains(const P::Assays& a, const std::string& id)
  {
    return std::any_of(a.compounds.begin(), a.compounds.end(), [&](const auto& c) { return c.id == "DECOY_" + id; });
  }

  void sequenceTests()
  {
    P::Assays targets;
    targets.proteins = {{"P2", "SECOND", "U2"}, {"unused", "UNUSED", "U0"}, {"P1", "FIRST", "U1"}};
    add(targets, "lysine", "PEPTIDEK");
    add(targets, "arginine", "GILGFVFTLR");
    add(targets, "nontryptic", "ACDEFGHI");
    add(targets, "internal", "AC(UniMod:4)DM(UniMod:35)EFGK", 2, {{1, 4}, {3, 35}});
    add(targets, "nterm", "(UniMod:1)PEPTIDEK", 2, {{-1, 1}});
    add(targets, "cterm", "PEPTIDEK.(UniMod:2)", 2, {{8, 2}});
    add(targets, "lastres", "PEPTIDEK(UniMod:259)", 2, {{7, 259}});
    add(targets, "charge3", "PEPTIDEK", 3);
    add(targets, "invalid", "AGHPEPTIDER");
    targets.transitions.back().product_mz = 1.0;
    add(targets, "filtered", "ASDFGHK");
    auto ignored = targets.transitions.back();
    ignored.product_mz = 1.0;
    ignored.detecting = false;
    ignored.transition_name += "_nondetecting";
    targets.transitions.push_back(ignored);
    ignored.detecting = true;
    ignored.decoy = true;
    ignored.transition_name += "_decoy";
    targets.transitions.push_back(ignored);
    // The first duplicate reserves its sequence before its invalid annotation rejects it.
    add(targets, "reserve_bad", "VAGHPEPTIDER");
    targets.transitions.back().product_mz = 1.0;
    add(targets, "reserve_later", "VAGHPEPTIDER");
    const auto actual = P::decoySequences(targets);
    equal(actual, oracle(targets), "sequence eligibility, metadata and order");
    equal(owned(light(targets)), targets, "test converter roundtrip");
    require(actual.transitions.empty(), "sequence stage must not retain transitions");
    for (const auto& id : {"lysine", "arginine", "nontryptic", "internal", "charge3", "filtered"})
      require(contains(actual, id), std::string("eligible fixture missing: ") + id);
    for (const auto& id : {"nterm", "cterm", "lastres", "invalid", "reserve_bad", "reserve_later"})
      require(!contains(actual, id), std::string("excluded fixture survived: ") + id);
    equal(P::decoySequences({}), oracle({}), "empty input");

    P::Assays collision;
    collision.proteins = targets.proteins;
    add(collision, "original", "PEPTIDEK");
    add(collision, "existing_target", "EDITPEPR");
    add(collision, "other_charge", "PEPTIDEK", 3);
    const auto collision_result = P::decoySequences(collision);
    equal(collision_result, oracle(collision), "target collision is charge-specific");
    require(!contains(collision_result, "original"), "decoy colliding with a target survived");
    require(contains(collision_result, "other_charge"), "collision leaked across precursor charges");
  }

#ifdef WITH_ONNX
  std::vector<char> bytes(const std::filesystem::path& path)
  {
    std::ifstream stream(path, std::ios::binary);
    require(stream.good(), "cannot read " + path.string());
    return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
  }

  // SQL diagnostics preserve column types and exact float bits, unlike a rounded text dump.
  std::vector<std::string> sqlRows(const std::string& path)
  {
    sqlite3* raw = nullptr;
    require(sqlite3_open_v2(path.c_str(), &raw, SQLITE_OPEN_READONLY, nullptr) == SQLITE_OK,
            "cannot open SQL diagnostic database " + path);
    const std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
    auto statement = [&](const std::string& sql)
    {
      sqlite3_stmt* stmt = nullptr;
      require(sqlite3_prepare_v2(db.get(), sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK,
              "cannot prepare SQL diagnostic: " + sql);
      return std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)>(stmt, sqlite3_finalize);
    };
    auto tables = statement("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name");
    std::vector<std::string> rows;
    while (sqlite3_step(tables.get()) == SQLITE_ROW)
    {
      const std::string table = reinterpret_cast<const char*>(sqlite3_column_text(tables.get(), 0));
      auto query = statement("SELECT * FROM \"" + table + "\"");
      std::vector<std::string> table_rows;
      while (sqlite3_step(query.get()) == SQLITE_ROW)
      {
        std::string row = table + ":";
        for (int c = 0; c < sqlite3_column_count(query.get()); ++c)
        {
          const int type = sqlite3_column_type(query.get(), c);
          row += std::to_string(type) + ":";
          if (type == SQLITE_FLOAT)
          {
            const double value = sqlite3_column_double(query.get(), c);
            row.append(reinterpret_cast<const char*>(&value), sizeof(value));
          }
          else if (type != SQLITE_NULL)
          {
            const auto* data = sqlite3_column_blob(query.get(), c);
            row.append(static_cast<const char*>(data), sqlite3_column_bytes(query.get(), c));
          }
          row += '\0';
        }
        table_rows.push_back(std::move(row));
      }
      std::sort(table_rows.begin(), table_rows.end());
      rows.insert(rows.end(), table_rows.begin(), table_rows.end());
    }
    return rows;
  }

  void predictionTests(const std::filesystem::path& scratch, const std::string& models, bool predict_ccs)
  {
    std::filesystem::create_directories(scratch);
    PeptDeepLibraryPredictor::Config old_config;
    P::AssayPredictor::Config config;
    if (!models.empty())
    {
      old_config.rt_model_path = config.rt_model_path = (std::filesystem::path(models) / "peptdeep_rt_dynamic.onnx").string();
      old_config.ccs_model_path = config.ccs_model_path = (std::filesystem::path(models) / "peptdeep_ccs_dynamic.onnx").string();
      old_config.ms2_model_path = config.ms2_model_path = (std::filesystem::path(models) / "peptdeep_ms2_dynamic.onnx").string();
    }
    old_config.batch_size = config.batch_size = 500;
    old_config.predict_ccs = config.predict_ccs = predict_ccs;
    old_config.intra_op_threads = config.intra_op_threads = 4;
    std::vector<PeptDeepLibraryPrecursor> old_inputs;
    std::vector<P::PredictionInput> inputs;
    const std::vector<std::string> sequences = {"PEPTIDEK", "GILGFVFTLR", "ACDEFGHIK", "VAGHPEPTIDER", "DASDFGHK", "LGGNEQVTR", "M(UniMod:35)PEPTIDER", "AC(UniMod:4)DEFGK", "VLGHPEPTIDER", "AGFAGDDAPR", "DLGEEHFK", "TVTAMDVVYALK"};
    for (const auto& sequence : sequences)
    {
      for (int charge : {2, 3})
      {
        const std::string id = "fixture_" + std::to_string(inputs.size());
        const auto peptide = AASequence::fromString(sequence);
        old_inputs.push_back({peptide, id, charge, 30.0f, 0, {"P1", "P2"}});
        inputs.push_back({peptide, id, charge, 30.0f, 0, {"P1", "P2"}});
      }
    }
    PeptDeepLibraryPredictor legacy(old_config);
    P::AssayPredictor predictor(config);
    auto targets = legacy.predict(old_inputs);
    equal(predictor.predict(inputs), owned(targets), "predictor exact fields");
    if (!predict_ccs)
    {
      for (const auto& c : targets.compounds) require(c.drift_time == -1.0, "CCS-disabled compound mobility");
      for (const auto& t : targets.transitions) require(t.precursor_im == -1.0, "CCS-disabled transition mobility");
    }
    // Terminal modifications exercise predictor metadata, but are intentionally absent
    // from the generator fixture because pseudo-reversal rejects terminal modifications.
    std::vector<PeptDeepLibraryPrecursor> terminal_old;
    std::vector<P::PredictionInput> terminal_inputs;
    for (const auto& sequence : {"(UniMod:1)PEPTIDEK", "PEPTIDEK.(UniMod:2)", "PEPTIDEK(UniMod:259)"})
    {
      const auto peptide = AASequence::fromString(sequence);
      const std::string id = "terminal_" + std::to_string(terminal_inputs.size());
      terminal_old.push_back({peptide, id, 2, 30.0f, 0, {"P1", "P2"}});
      terminal_inputs.push_back({peptide, id, 2, 30.0f, 0, {"P1", "P2"}});
    }
    const auto terminal_expected = owned(legacy.predict(terminal_old));
    equal(predictor.predict(terminal_inputs), terminal_expected, "terminal-mod predictor exact fields");
    require(terminal_expected.compounds[0].modifications == std::vector<P::Modification>{{-1, 1}},
            "N-terminal predictor modification location");
    require(terminal_expected.compounds[1].modifications == std::vector<P::Modification>{{8, 2}},
            "C-terminal predictor modification location");
    require(terminal_expected.compounds[2].modifications == std::vector<P::Modification>{{7, 259}},
            "last-residue predictor modification location");
    // A library assay is a selected subset, not the predictor's exhaustive channel grid.
    // In particular, b1 is not annotatable by the historical MRM ion series. The no-IM
    // fixture requests six fragments to also exercise retention of short decoy assays.
    const size_t fragment_count = predict_ccs ? 3 : 6;
    std::vector<OpenSwath::LightTransition> selected;
    for (const auto& compound : targets.compounds)
    {
      std::vector<OpenSwath::LightTransition> candidates;
      for (const auto& transition : targets.transitions)
      {
        if (transition.peptide_ref == compound.id && transition.fragment_charge == 1 &&
            transition.fragment_nr >= 2 && transition.product_mz >= 350.0 && transition.product_mz <= 2000.0)
        {
          candidates.push_back(transition);
        }
      }
      std::stable_sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b)
      {
        return a.library_intensity > b.library_intensity;
      });
      require(candidates.size() >= fragment_count, "prediction fixture has too few usable fragments");
      candidates.resize(fragment_count);
      selected.insert(selected.end(), candidates.begin(), candidates.end());
    }
    targets.transitions = std::move(selected);
    // Generate the input anew on every invocation and read it back to canonicalize PQP IDs.
    const auto input_path = (scratch / "preparation_targets.pqp").string();
    const auto old_path = (scratch / "preparation_legacy.pqp").string();
    const auto new_path = (scratch / "preparation_owned.pqp").string();
    for (const auto& path : {input_path, old_path, new_path}) std::filesystem::remove(path);
    TransitionPQPFile writer;
    writer.setLogType(ProgressLogger::NONE);
    writer.convertLightTargetedExperimentToPQP(input_path.c_str(), targets);
    OpenSwath::LightTargetedExperiment canonical;
    writer.convertPQPToTargetedExperiment(input_path.c_str(), canonical);
    const auto sequence_decoys = oracle(owned(canonical));
    require(sequence_decoys.compounds.size() * 5 >= canonical.compounds.size() * 4,
            "fixture annotation eligibility below 80%: " + std::to_string(sequence_decoys.compounds.size()) +
            "/" + std::to_string(canonical.compounds.size()));
    OpenSwathLibraryPreparation preparation;
    preparation.setLogType(ProgressLogger::NONE);
    OpenSwathLibraryPreparation::DecoyGeneratorParameters options;
    options.method = "pseudo-reverse";
    options.predict_decoys = true;
    options.predict_keep_short = true;
    options.predict_threads = 4;
    options.predict_model_directory = models;
    preparation.generateDecoys(input_path, FileTypes::PQP, old_path, FileTypes::PQP, options);
    P::PredictionOptions prediction_options;
    prediction_options.model_directory = models;
    prediction_options.threads = 4;
    prediction_options.log_type = ProgressLogger::NONE;
    const auto result = P::predictDecoys(owned(canonical), prediction_options);
    if (!predict_ccs)
    {
      size_t short_assays = 0;
      for (const auto& compound : result.compounds)
      {
        size_t fragments = 0;
        for (const auto& transition : result.transitions)
        {
          fragments += transition.decoy && transition.peptide_ref == compound.id;
        }
        short_assays += fragments > 0 && fragments < fragment_count;
      }
      require(short_assays > 0, "no-IM fixture did not exercise keep-short decoy policy");
    }
    writer.convertLightTargetedExperimentToPQP(new_path.c_str(), light(result));
    if (bytes(old_path) != bytes(new_path))
    {
      require(sqlRows(old_path) == sqlRows(new_path), "fresh generated PQP bytes and SQL rows differ");
      throw std::runtime_error("Preparation parity: fresh generated PQP bytes differ although SQL rows agree");
    }
    if (predict_ccs)
    {
      // Exercise both terms of the default coverage gate, without weakening its threshold.
      for (const bool protein_failure : {false, true})
      {
        auto failure_targets = owned(canonical);
        const std::string kind = protein_failure ? "protein" : "compound";
        if (protein_failure)
        {
          failure_targets.proteins.push_back({"excluded_compound_protein", "EXCLUDED", "EXCLUDED"});
          failure_targets.compounds.front().protein_refs.push_back("excluded_compound_protein");
          for (auto& transition : failure_targets.transitions)
          {
            const auto it = std::find_if(failure_targets.compounds.begin(), failure_targets.compounds.end(),
              [&](const auto& c) { return c.id == transition.peptide_ref; });
            // PQP protein mappings are peptide-level, shared by both precursor charges.
            if (it != failure_targets.compounds.end() && it->sequence == failure_targets.compounds.front().sequence)
            {
              transition.product_mz = 1.0;
            }
          }
        }
        else
        {
          for (auto& transition : failure_targets.transitions)
          {
            // Reject a third of the compounds at annotation, keeping enough other targets
            // for instrument/NCE estimation and making the failure specifically coverage.
            const auto it = std::find_if(failure_targets.compounds.begin(), failure_targets.compounds.end(),
              [&](const auto& c) { return c.id == transition.peptide_ref; });
            if (std::distance(failure_targets.compounds.begin(), it) < 8) transition.product_mz = 1.0;
          }
        }
        const auto failure_input = (scratch / ("coverage_" + kind + "_targets.pqp")).string();
        const auto failure_output = (scratch / ("coverage_" + kind + "_legacy.pqp")).string();
        std::filesystem::remove(failure_input);
        std::filesystem::remove(failure_output);
        writer.convertLightTargetedExperimentToPQP(failure_input.c_str(), light(failure_targets));
        OpenSwath::LightTargetedExperiment failure_canonical;
        writer.convertPQPToTargetedExperiment(failure_input.c_str(), failure_canonical);
        require(failure_canonical.proteins.size() == failure_targets.proteins.size(), "coverage fixture lost proteins");
        const auto eligible = oracle(owned(failure_canonical));
        require((eligible.compounds.size() * 5 >= failure_canonical.compounds.size() * 4) == protein_failure,
                "coverage fixture does not isolate the intended compound/protein gate");
        require((eligible.proteins.size() * 5 < failure_canonical.proteins.size() * 4) == protein_failure,
                "coverage fixture does not isolate the intended protein gate");
        auto coverageFailure = [&](auto&& action, const std::string& implementation)
        {
          bool threw = false;
          try
          {
            action();
          }
          catch (const Exception::IllegalArgument& error)
          {
            const std::string message = error.getMessage();
            require(message.find("below the threshold of 80") != std::string::npos,
                    implementation + " threw a different IllegalArgument: " + message);
            threw = true;
          }
          require(threw, implementation + " failed to reject low " + kind + " coverage");
        };
        coverageFailure([&]
        {
          preparation.generateDecoys(failure_input, FileTypes::PQP, failure_output, FileTypes::PQP, options);
        }, "legacy");
        coverageFailure([&]
        {
          P::predictDecoys(owned(failure_canonical), prediction_options);
        }, "owned");
      }
    }
  }
#endif
}

int main(int argc, char** argv)
{
  try
  {
    require(argc >= 2, "usage: Preparation_test SCRATCH [--prediction [MODEL_DIRECTORY]]");
    std::filesystem::create_directories(argv[1]);
    sequenceTests();
    if (argc > 2)
    {
      require(std::string(argv[2]) == "--prediction", "unknown option");
#ifdef WITH_ONNX
      predictionTests(argv[1], argc > 3 ? argv[3] : "", true);
      predictionTests(std::filesystem::path(argv[1]) / "no_im", argc > 3 ? argv[3] : "", false);
#else
      throw std::runtime_error("Preparation parity: --prediction requires WITH_ONNX");
#endif
    }
    std::cout << "Preparation parity passed\n";
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
