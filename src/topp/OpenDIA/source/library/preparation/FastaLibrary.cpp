// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
#include <library/preparation/FastaLibrary.h>

#include <library/Library.h>
#include <library/preparation/DecoySequences.h>

#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/config.h>
#include <OpenMS/openms_data_path.h>

#ifdef WITH_ONNX
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepCCSInference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepRTInference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#endif

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <map>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace OpenMS::OpenDIA::Preparation
{
  namespace
  {
    using Log = std::function<void(const std::string&)>;

    void say(const Log& log, const std::string& line)
    {
      if (log) log(line);
    }

    double since(std::chrono::steady_clock::time_point t)
    {
      return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
    }

    [[noreturn]] void invalid(const std::string& message, const std::string& value)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, value);
    }

    [[noreturn]] void needsOnnx()
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "A FASTA library needs an OpenMS build configured with WITH_ONNX=ON (PeptDeep).");
    }

    bool standardResidues(std::string_view s)
    {
      return std::all_of(s.begin(), s.end(), [](char c) { return std::string_view("ACDEFGHIKLMNPQRSTVWY").find(c) != std::string_view::npos; });
    }

    /// "sp|P12345|NAME_HUMAN" -> "P12345"; any other identifier as is.
    std::string accession(const std::string& identifier)
    {
      const auto a = identifier.find('|');
      if (a == std::string::npos) return identifier;
      const auto b = identifier.find('|', a + 1);
      const std::string db = identifier.substr(0, a);
      if (b == std::string::npos || (db != "sp" && db != "tr")) return identifier;
      return identifier.substr(a + 1, b - a - 1);
    }

    ModifiedPeptideGenerator::MapToResidueType modifications(const std::vector<std::string>& names)
    {
      StringList list(names.begin(), names.end());
      const auto map = ModifiedPeptideGenerator::getModifications(list);
      for (const auto& [mod, residue] : map.val)
      {
        if (mod->getTermSpecificity() != ResidueModification::ANYWHERE || residue == nullptr)
          invalid("Terminal modifications are not supported (their peptides get no pseudo-reversed decoy)", mod->getFullId());
        if (mod->getUniModRecordId() <= 0 || mod->getDiffFormula().isEmpty())
          invalid("Modifications need a UniMod ID and a formula (PeptDeep encodes them by composition)", mod->getFullId());
      }
      return map;
    }

    std::vector<Modification> residueModifications(const AASequence& peptide)
    {
      std::vector<Modification> mods;
      for (Size i = 0; i < peptide.size(); ++i)
      {
        if (peptide[i].isModified()) mods.push_back({static_cast<int>(i), peptide[i].getModification()->getUniModRecordId()});
      }
      return mods;
    }

    constexpr const char* PANEL_KEY = "opendia.fasta.panel_precursors";
    constexpr std::size_t PANEL_FRAGMENTS = 6;
    constexpr std::size_t PANEL_MIN_FRAGMENTS = 5;  ///< CalibrationOptions::strict_min_fragments.
    constexpr std::size_t SEARCH_FRAGMENTS = 6;

#ifdef WITH_ONNX
    /// PeptDeep MS2 for a batch of (sequence, charge), with one setting.
    std::vector<std::vector<float>> predictMS2(PeptDeepMS2Inference& model, const std::vector<std::string>& sequences,
                                               const std::vector<int>& charges, const PredictionSettings& settings)
    {
      std::vector<std::string> peptides;
      peptides.reserve(sequences.size());
      for (const auto& s : sequences) peptides.push_back(AASequence::fromString(s).toString());
      const std::vector<float> z(charges.begin(), charges.end());
      const std::vector<float> nce(sequences.size(), settings.nce);
      const std::vector<std::int64_t> instrument(sequences.size(), ML::instrumentIndex(settings.instrument));
      auto out = model.predictMS2(peptides, z, nce, instrument);
      if (out.size() != sequences.size()) invalid("PeptDeep returned an unexpected number of MS2 predictions", std::to_string(out.size()));
      return out;
    }
#endif
  }

  std::string CandidateSet::tramlId(std::size_t precursor) const
  {
    const auto& p = precursors[precursor];
    const auto& target = precursors[precursor & ~std::size_t{1}];
    std::string id = peptidoforms[target.peptidoform].sequence + "/" + std::to_string(p.charge);
    return (precursor & 1) ? "DECOY_" + id : id;
  }

  int missedCleavages(const std::string& enzyme, const std::string& unmodified)
  {
    ProteaseDigestion digestion;
    digestion.setEnzyme(enzyme);
    digestion.setMissedCleavages(0);
    std::vector<std::pair<Size, Size>> out;
    digestion.digestUnmodified(unmodified, out, 1, 0);
    return out.empty() ? 0 : static_cast<int>(out.size()) - 1;
  }

  CandidateSet enumerateCandidates(const std::string& fasta, const CandidateOptions& options, const Log& log)
  {
    const auto start = std::chrono::steady_clock::now();
    if (options.min_length < 2 || options.max_length < options.min_length || options.max_length > 60)
      invalid("Peptide lengths must satisfy 2 <= min <= max <= 60", std::to_string(options.min_length) + "-" + std::to_string(options.max_length));
    if (options.min_charge < 1 || options.max_charge < options.min_charge || options.max_charge > 10)
      invalid("Precursor charges must satisfy 1 <= min <= max <= 10", std::to_string(options.min_charge) + "-" + std::to_string(options.max_charge));
    const auto fixed = modifications(options.fixed_mods);
    const auto variable = modifications(options.variable_mods);

    CandidateSet set;
    // Distinct peptides in order of first appearance, with the proteins they occur in.
    std::vector<std::string> peptides;
    std::vector<std::vector<std::uint32_t>> peptide_proteins;
    std::vector<std::string> proteins;
    {
      ProteaseDigestion digestion;
      digestion.setEnzyme(options.enzyme);
      digestion.setMissedCleavages(static_cast<Size>(options.missed_cleavages));
      std::unordered_map<std::string, std::uint32_t> index;
      std::unordered_set<std::string> skipped;
      FASTAFile file;
      file.readStart(fasta);
      FASTAFile::FASTAEntry entry;
      std::vector<std::string_view> products;
      while (file.readNext(entry))
      {
        const auto protein = static_cast<std::uint32_t>(proteins.size());
        proteins.push_back(accession(entry.identifier));
        products.clear();
        digestion.digestUnmodified(entry.sequence, products, static_cast<Size>(options.min_length), static_cast<Size>(options.max_length));
        for (const auto product : products)
        {
          if (!standardResidues(product))
          {
            skipped.emplace(product);
            continue;
          }
          auto [it, inserted] = index.try_emplace(std::string(product), static_cast<std::uint32_t>(peptides.size()));
          if (inserted)
          {
            peptides.emplace_back(product);
            peptide_proteins.emplace_back();
          }
          auto& list = peptide_proteins[it->second];
          if (list.empty() || list.back() != protein) list.push_back(protein);
        }
      }
      set.proteins_read = proteins.size();
      set.skipped_residues = skipped.size();
    }
    say(log, "FASTA: proteins=" + std::to_string(proteins.size()) + " peptides=" + std::to_string(peptides.size()) +
             " skipped_nonstandard=" + std::to_string(set.skipped_residues) + " seconds=" + std::to_string(since(start)));

    // Target peptidoforms (all of them first, so a decoy colliding with any target is detected).
    struct Target
    {
      std::uint32_t peptide;
      std::string sequence;
      std::vector<Modification> mods;
      std::uint8_t variable;
    };
    std::vector<Target> targets;
    std::unordered_set<std::string> target_sequences;
    for (std::uint32_t i = 0; i < peptides.size(); ++i)
    {
      AASequence base = AASequence::fromString(peptides[i]);
      ModifiedPeptideGenerator::applyFixedModifications(fixed, base);
      const auto fixed_count = residueModifications(base).size();
      std::vector<AASequence> forms;
      ModifiedPeptideGenerator::applyVariableModifications(variable, base, static_cast<Size>(options.max_variable_mods), forms, true);
      for (const auto& form : forms)
      {
        Target t{i, {}, residueModifications(form), 0};
        t.variable = static_cast<std::uint8_t>(t.mods.size() - fixed_count);
        t.sequence = uniModSequence(peptides[i], t.mods);
        if (target_sequences.insert(t.sequence).second) targets.push_back(std::move(t));
      }
    }

    // Protein sets (targets and DECOY_ sets), deduplicated.
    std::map<std::vector<std::uint32_t>, std::uint32_t> target_set, decoy_set;
    const auto protein_set = [&](const std::vector<std::uint32_t>& members, bool decoy)
    {
      auto& sets = decoy ? decoy_set : target_set;
      auto [it, inserted] = sets.try_emplace(members, static_cast<std::uint32_t>(set.protein_sets.size()));
      if (inserted)
      {
        std::vector<std::string> names;
        for (const auto m : members) names.push_back(decoy ? "DECOY_" + proteins[m] : proteins[m]);
        set.protein_sets.push_back(std::move(names));
      }
      return it->second;
    };

    ProteaseDigestion cleavage;
    cleavage.setEnzyme(options.enzyme);
    cleavage.setMissedCleavages(0);
    std::unordered_set<std::string> decoy_sequences;
    std::vector<bool> peptide_used(peptides.size(), false);
    std::vector<std::pair<Size, Size>> sites;
    for (auto& t : targets)
    {
      std::string reversed = peptides[t.peptide];
      std::vector<Modification> mods = t.mods;
      if (!pseudoReverse(reversed, mods))
      {
        ++set.dropped_no_decoy;
        continue;
      }
      std::string decoy = uniModSequence(reversed, mods);
      if (target_sequences.count(decoy) || !decoy_sequences.insert(decoy).second)
      {
        ++set.dropped_collision;
        continue;
      }
      const AASequence target_peptide = AASequence::fromString(t.sequence);
      const IonLadder ladder(target_peptide);
      sites.clear();
      cleavage.digestUnmodified(peptides[t.peptide], sites, 1, 0);
      Peptidoform tp, dp;
      tp.sequence = std::move(t.sequence);
      tp.proteins = protein_set(peptide_proteins[t.peptide], false);
      tp.missed_cleavages = static_cast<std::uint8_t>(sites.empty() ? 0 : sites.size() - 1);
      tp.variable_mods = t.variable;
      tp.length = static_cast<std::uint8_t>(peptides[t.peptide].size());
      dp = tp;
      dp.sequence = std::move(decoy);
      dp.proteins = protein_set(peptide_proteins[t.peptide], true);
      dp.decoy = true;
      const IonLadder decoy_ladder(AASequence::fromString(dp.sequence));
      bool added = false;
      for (int z = options.min_charge; z <= options.max_charge; ++z)
      {
        const double mz = ladder.precursor(z);
        if (mz < options.min_mz || mz > options.max_mz) continue;
        if (!added)
        {
          set.peptidoforms.push_back(tp);
          set.peptidoforms.push_back(dp);
          added = true;
        }
        const auto tpi = static_cast<std::uint32_t>(set.peptidoforms.size() - 2);
        set.precursors.push_back({tpi, mz, -1.0f, static_cast<std::int8_t>(z)});
        set.precursors.push_back({tpi + 1, decoy_ladder.precursor(z), -1.0f, static_cast<std::int8_t>(z)});
      }
      if (added) peptide_used[t.peptide] = true;
    }
    set.target_peptides = static_cast<std::size_t>(std::count(peptide_used.begin(), peptide_used.end(), true));
    say(log, "FASTA candidates: target_peptidoforms=" + std::to_string(targets.size()) + " pairs=" + std::to_string(set.pairs()) +
             " target_peptides=" + std::to_string(set.target_peptides) + " dropped_no_decoy=" + std::to_string(set.dropped_no_decoy) +
             " dropped_collision=" + std::to_string(set.dropped_collision) + " seconds=" + std::to_string(since(start)));
    return set;
  }

  IonLadder::IonLadder(const AASequence& peptide) : cumulative_(peptide.size() + 1, 0.0)
  {
    if (peptide.hasNTerminalModification() || peptide.hasCTerminalModification())
      invalid("Terminal modifications are not supported", peptide.toString());
    for (Size i = 0; i < peptide.size(); ++i) cumulative_[i + 1] = cumulative_[i] + peptide[i].getMonoWeight(Residue::Internal);
  }

  double IonLadder::ion(char type, std::size_t ordinal, int charge) const
  {
    static const double b_adj = Residue::getInternalToBIon().getMonoWeight();
    static const double y_adj = Residue::getInternalToYIon().getMonoWeight();
    const std::size_t n = size();
    const double mass = type == 'b' ? cumulative_[ordinal] + b_adj : cumulative_[n] - cumulative_[n - ordinal] + y_adj;
    return (mass + Constants::PROTON_MASS_U * charge) / charge;
  }

  double IonLadder::precursor(int charge) const
  {
    static const double full = Residue::getInternalToFull().getMonoWeight();
    return (cumulative_.back() + full + Constants::PROTON_MASS_U * charge) / charge;
  }

  std::vector<FragmentIon> surveyFragments(const IonLadder& ladder, double min_mz, double max_mz)
  {
    const std::size_t n = ladder.size();
    std::vector<FragmentIon> out;
    const auto pick = [&](char type, int centre)
    {
      std::vector<FragmentIon> ions;
      for (std::size_t o = 3; o + 1 <= n; ++o)
      {
        const double mz = ladder.ion(type, o, 1);
        if (mz >= min_mz && mz <= max_mz) ions.push_back({mz, 1.0f, static_cast<std::int16_t>(o), 1, type});
      }
      std::stable_sort(ions.begin(), ions.end(), [centre](const FragmentIon& a, const FragmentIon& b)
      {
        const int da = std::abs(a.ordinal - centre), db = std::abs(b.ordinal - centre);
        return da != db ? da < db : a.ordinal < b.ordinal;
      });
      for (const auto& ion : ions)
      {
        if (out.size() == SEARCH_FRAGMENTS) break;
        out.push_back(ion);
      }
    };
    pick('y', 7);
    pick('b', 5);
    return out;
  }

  std::vector<FragmentIon> predictedFragments(const IonLadder& ladder, int precursor_charge, const std::vector<float>& ms2,
                                              double min_mz, double max_mz, std::size_t count)
  {
    const std::size_t n = ladder.size();
    constexpr std::size_t channels = 8;
    if (n < 2 || ms2.size() != (n - 1) * channels) invalid("Unexpected PeptDeep MS2 output size", std::to_string(ms2.size()));
    const int max_charge = std::max(1, std::min(2, precursor_charge - 1));
    std::vector<FragmentIon> ions;
    for (std::size_t row = 0; row + 1 < n; ++row)
    {
      const std::size_t b = row + 1, y = n - 1 - row;
      const std::pair<char, std::size_t> kinds[4] = {{'b', b}, {'b', b}, {'y', y}, {'y', y}};
      for (std::size_t c = 0; c < 4; ++c)
      {
        const int z = static_cast<int>(c % 2) + 1;
        const auto [type, ordinal] = kinds[c];
        const float intensity = ms2[row * channels + c];
        if (ordinal < 3 || z > max_charge || !(intensity > 0.0f)) continue;
        const double mz = ladder.ion(type, ordinal, z);
        if (mz < min_mz || mz > max_mz) continue;
        ions.push_back({mz, intensity, static_cast<std::int16_t>(ordinal), static_cast<std::int8_t>(z), type});
      }
    }
    std::sort(ions.begin(), ions.end(), [](const FragmentIon& a, const FragmentIon& b)
    {
      return a.intensity != b.intensity ? a.intensity > b.intensity : a.mz < b.mz;
    });
    if (ions.size() > count) ions.resize(count);
    return ions;
  }

  std::vector<std::uint32_t> calibrationPanel(const CandidateSet& candidates, std::size_t count)
  {
    std::vector<std::uint32_t> eligible;
    for (std::uint32_t k = 0; k < candidates.pairs(); ++k)
    {
      const auto& p = candidates.precursors[2 * k];
      const auto& f = candidates.peptidoforms[p.peptidoform];
      if (p.charge == 2 && f.variable_mods == 0 && f.missed_cleavages == 0 && f.length >= 7 && f.length <= 20) eligible.push_back(k);
    }
    const auto rt = [&](std::uint32_t k) { return candidates.peptidoforms[candidates.precursors[2 * k].peptidoform].rt; };
    std::stable_sort(eligible.begin(), eligible.end(), [&](auto a, auto b) { return rt(a) != rt(b) ? rt(a) < rt(b) : a < b; });
    const std::size_t take = std::min(count, eligible.size());
    std::vector<std::uint32_t> panel;
    panel.reserve(take);
    for (std::size_t i = 0; i < take; ++i) panel.push_back(eligible[i * eligible.size() / take]);
    return panel;
  }

  std::vector<std::uint32_t> selectPairs(const std::vector<double>& statistic, const std::vector<std::uint64_t>& tie, std::size_t keep)
  {
    std::vector<std::uint32_t> pairs;
    if (keep == 0)
    {
      pairs.resize(statistic.size());
      std::iota(pairs.begin(), pairs.end(), 0u);
      return pairs;
    }
    for (std::uint32_t k = 0; k < statistic.size(); ++k)
      if (statistic[k] > 0.0) pairs.push_back(k);
    const auto better = [&](std::uint32_t a, std::uint32_t b)
    {
      if (statistic[a] != statistic[b]) return statistic[a] > statistic[b];
      if (tie[a] != tie[b]) return tie[a] < tie[b];
      return a < b;
    };
    if (pairs.size() > keep)
    {
      std::nth_element(pairs.begin(), pairs.begin() + static_cast<std::ptrdiff_t>(keep), pairs.end(), better);
      pairs.resize(keep);
    }
    std::sort(pairs.begin(), pairs.end());
    return pairs;
  }

  std::string modelPath(const PredictionSettings& settings, const std::string& name)
  {
    if (!settings.model_directory.empty()) return settings.model_directory + "/" + name;
    return File::find("models/" + name, {std::string(OPENMS_BINARY_PATH) + "/share/OpenMS"});
  }

  void predictCoordinates(CandidateSet& candidates, const PredictionSettings& settings, const Log& log)
  {
#ifndef WITH_ONNX
    (void)candidates;
    (void)settings;
    (void)log;
    needsOnnx();
#else
    constexpr std::size_t CHUNK = 200000;
    auto start = std::chrono::steady_clock::now();
    {
      PeptDeepRTInference model(modelPath(settings, "peptdeep_rt_dynamic.onnx"), settings.threads);
      std::vector<std::string> batch;
      for (std::size_t begin = 0; begin < candidates.peptidoforms.size(); begin += CHUNK)
      {
        const std::size_t end = std::min(candidates.peptidoforms.size(), begin + CHUNK);
        batch.clear();
        for (std::size_t i = begin; i < end; ++i) batch.push_back(AASequence::fromString(candidates.peptidoforms[i].sequence).toString());
        const auto rt = model.predictRT(batch);
        if (rt.size() != batch.size()) invalid("PeptDeep returned an unexpected number of RT predictions", std::to_string(rt.size()));
        for (std::size_t i = begin; i < end; ++i) candidates.peptidoforms[i].rt = 100.0f * rt[i - begin];
      }
    }
    say(log, "FASTA RT: peptidoforms=" + std::to_string(candidates.peptidoforms.size()) + " seconds=" + std::to_string(since(start)));
    if (!settings.predict_im) return;
    start = std::chrono::steady_clock::now();
    PeptDeepCCSInference model(modelPath(settings, "peptdeep_ccs_dynamic.onnx"), settings.threads);
    std::vector<std::string> batch;
    std::vector<float> charges;
    for (std::size_t begin = 0; begin < candidates.precursors.size(); begin += CHUNK)
    {
      const std::size_t end = std::min(candidates.precursors.size(), begin + CHUNK);
      batch.clear();
      charges.clear();
      for (std::size_t i = begin; i < end; ++i)
      {
        const auto& p = candidates.precursors[i];
        batch.push_back(AASequence::fromString(candidates.peptidoforms[p.peptidoform].sequence).toString());
        charges.push_back(static_cast<float>(p.charge));
      }
      const auto ccs = model.predictCCS(batch, charges);
      if (ccs.size() != batch.size()) invalid("PeptDeep returned an unexpected number of CCS predictions", std::to_string(ccs.size()));
      for (std::size_t i = begin; i < end; ++i)
      {
        auto& p = candidates.precursors[i];
        p.im = static_cast<float>(IMTypes::ccsToOneOverK0(ccs[i - begin], p.mz, p.charge));
      }
    }
    say(log, "FASTA 1/K0: precursors=" + std::to_string(candidates.precursors.size()) + " seconds=" + std::to_string(since(start)));
#endif
  }

  namespace
  {
    ParquetLibraryWriter::Fragment row(std::size_t precursor, std::size_t slot, const FragmentIon& ion)
    {
      ParquetLibraryWriter::Fragment f;
      f.transition_id = std::to_string(precursor * 100 + slot);
      f.product_mz = ion.mz;
      f.intensity = ion.intensity;
      f.type = ion.type == 'b' ? "b" : "y";
      f.ordinal = ion.ordinal;
      f.charge = ion.charge;
      return f;
    }

    std::vector<ParquetLibraryWriter::Fragment> rows(std::size_t precursor, const std::vector<FragmentIon>& ions, std::size_t count)
    {
      std::vector<ParquetLibraryWriter::Fragment> out;
      for (std::size_t s = 0; s < std::min(count, ions.size()); ++s) out.push_back(row(precursor, s, ions[s]));
      return out;
    }
  }

  void writeSurveyLibrary(const CandidateSet& candidates, const std::vector<std::uint32_t>& panel, const PredictionSettings& settings,
                          const FragmentOptions& fragments, const std::string& path,
                          const std::vector<std::pair<std::string, std::string>>& metadata, const Log& log)
  {
#ifndef WITH_ONNX
    (void)candidates; (void)panel; (void)settings; (void)fragments; (void)path; (void)metadata; (void)log;
    needsOnnx();
#else
    const auto start = std::chrono::steady_clock::now();
    // Panel pairs with predicted fragments; a pair with a member below the calibration's 5 fragments is not a panel pair.
    std::vector<std::vector<FragmentIon>> panel_ions(2 * panel.size());
    std::vector<bool> is_panel(candidates.pairs(), false);
    std::vector<std::uint32_t> panel_pairs;
    {
      PeptDeepMS2Inference model(modelPath(settings, "peptdeep_ms2_dynamic.onnx"), settings.threads);
      std::vector<std::string> sequences;
      std::vector<int> charges;
      for (const auto k : panel)
      {
        for (std::size_t m = 0; m < 2; ++m)
        {
          const auto& p = candidates.precursors[2 * k + m];
          sequences.push_back(candidates.peptidoforms[p.peptidoform].sequence);
          charges.push_back(p.charge);
        }
      }
      const auto ms2 = sequences.empty() ? std::vector<std::vector<float>>() : predictMS2(model, sequences, charges, settings);
      for (std::size_t i = 0; i < panel.size(); ++i)
      {
        bool ok = true;
        for (std::size_t m = 0; m < 2; ++m)
        {
          const IonLadder ladder(AASequence::fromString(sequences[2 * i + m]));
          panel_ions[2 * i + m] = predictedFragments(ladder, charges[2 * i + m], ms2[2 * i + m], fragments.min_mz, fragments.max_mz, PANEL_FRAGMENTS);
          ok = ok && panel_ions[2 * i + m].size() >= PANEL_MIN_FRAGMENTS;
        }
        if (ok)
        {
          is_panel[panel[i]] = true;
          panel_pairs.push_back(static_cast<std::uint32_t>(i));
        }
      }
    }
    auto meta = metadata;
    meta.emplace_back(PANEL_KEY, std::to_string(2 * panel_pairs.size()));
    ParquetLibraryWriter writer(path, meta);
    const auto add = [&](std::size_t precursor, const std::vector<ParquetLibraryWriter::Fragment>& rows_)
    {
      const auto& p = candidates.precursors[precursor];
      const auto& f = candidates.peptidoforms[p.peptidoform];
      const std::string id = std::to_string(precursor), traml = candidates.tramlId(precursor);
      ParquetLibraryWriter::Precursor out;
      out.id = id;
      out.traml_id = traml;
      out.sequence = f.sequence;
      out.mz = p.mz;
      out.rt = f.rt;
      out.im = p.im;
      out.charge = p.charge;
      out.decoy = f.decoy;
      out.proteins = &candidates.protein_sets[f.proteins];
      writer.add(out, rows_);
    };
    for (const auto i : panel_pairs)
    {
      for (std::size_t m = 0; m < 2; ++m) add(2 * panel[i] + m, rows(2 * panel[i] + m, panel_ions[2 * i + m], PANEL_FRAGMENTS));
    }
    std::size_t dropped = 0;
    for (std::size_t k = 0; k < candidates.pairs(); ++k)
    {
      if (is_panel[k]) continue;
      std::vector<FragmentIon> ions[2];
      for (std::size_t m = 0; m < 2; ++m)
      {
        const auto& f = candidates.peptidoforms[candidates.precursors[2 * k + m].peptidoform];
        ions[m] = surveyFragments(IonLadder(AASequence::fromString(f.sequence)), fragments.min_mz, fragments.max_mz);
      }
      if (ions[0].empty() || ions[1].empty())
      {
        ++dropped;
        continue;
      }
      for (std::size_t m = 0; m < 2; ++m) add(2 * k + m, rows(2 * k + m, ions[m], SEARCH_FRAGMENTS));
    }
    writer.close();
    say(log, "FASTA survey library: " + path + " precursors=" + std::to_string(writer.precursors()) + " panel_pairs=" +
             std::to_string(panel_pairs.size()) + " dropped_without_fragments=" + std::to_string(dropped) + " seconds=" + std::to_string(since(start)));
#endif
  }

  std::size_t surveyPanelPrecursors(const std::string& path)
  {
    const auto meta = readParquetLibraryMetadata(path);
    const auto it = meta.find(PANEL_KEY);
    return it == meta.end() ? 0 : static_cast<std::size_t>(std::stoull(it->second));
  }

  FinalLibraryStats writeFinalLibrary(const std::string& survey, const std::function<bool(std::size_t)>& keep,
                                      const PredictionSettings& settings, const FragmentOptions& fragments, const std::string& path,
                                      const std::string& search_path, const std::vector<std::pair<std::string, std::string>>& metadata,
                                      const Log& log)
  {
#ifndef WITH_ONNX
    (void)survey; (void)keep; (void)settings; (void)fragments; (void)path; (void)search_path; (void)metadata; (void)log;
    needsOnnx();
#else
    if (fragments.count < 1 || fragments.count > 99) invalid("Library fragments per precursor must be in 1-99", std::to_string(fragments.count));
    const auto start = std::chrono::steady_clock::now();
    struct Record
    {
      std::size_t id = 0;
      std::string sequence;
      std::vector<std::string> proteins;
      double mz = 0.0, rt = 0.0, im = -1.0;
      int charge = 0;
    };
    // The kept precursors of the survey library (both members of a kept pair), one record each.
    struct Collector final : LibraryVisitor
    {
      const std::function<bool(std::size_t)>& keep;
      std::vector<Record> records;
      std::size_t last = SIZE_MAX;
      explicit Collector(const std::function<bool(std::size_t)>& k) : keep(k) {}
      void reset() override { records.clear(); last = SIZE_MAX; }
      void fragment(const FragmentView&, const PrecursorView& c) override
      {
        std::size_t id = 0;
        const auto r = std::from_chars(c.id.data(), c.id.data() + c.id.size(), id);
        if (r.ec != std::errc() || r.ptr != c.id.data() + c.id.size()) invalid("FASTA survey library precursor IDs must be decimal", std::string(c.id));
        if (id == last || !keep(id / 2)) return;
        last = id;
        Record rec;
        rec.id = id;
        rec.sequence = c.sequence;
        for (const auto p : *c.proteins) rec.proteins.emplace_back(p);
        rec.rt = c.rt;
        rec.im = c.drift_time;
        rec.charge = c.charge;
        records.push_back(std::move(rec));
      }
    } collector(keep);
    // Precursor m/z is recomputed from the sequence and the decoy flag is implied by the ID (odd = decoy).
    visitParquetLibrary(survey, collector);
    auto& records = collector.records;
    std::sort(records.begin(), records.end(), [](const Record& a, const Record& b) { return a.id < b.id; });
    // Complete pairs only.
    std::vector<Record> pairs;
    pairs.reserve(records.size());
    for (std::size_t i = 0; i + 1 < records.size();)
    {
      if (records[i].id % 2 == 0 && records[i + 1].id == records[i].id + 1)
      {
        pairs.push_back(std::move(records[i]));
        pairs.push_back(std::move(records[i + 1]));
        i += 2;
      }
      else ++i;
    }
    std::vector<Record>().swap(records);
    FinalLibraryStats stats;
    stats.pairs_in = pairs.size() / 2;
    say(log, "FASTA final library: predicting MS2 for " + std::to_string(stats.pairs_in) + " pairs");

    ParquetLibraryWriter scoring(path, metadata), search(search_path, metadata);
    PeptDeepMS2Inference model(modelPath(settings, "peptdeep_ms2_dynamic.onnx"), settings.threads);
    constexpr std::size_t CHUNK = 50000;  // precursors (even: pairs stay together)
    std::vector<std::string> sequences;
    std::vector<int> charges;
    for (std::size_t begin = 0; begin < pairs.size(); begin += CHUNK)
    {
      const std::size_t end = std::min(pairs.size(), begin + CHUNK);
      sequences.clear();
      charges.clear();
      for (std::size_t i = begin; i < end; ++i)
      {
        sequences.push_back(pairs[i].sequence);
        charges.push_back(pairs[i].charge);
      }
      const auto ms2 = predictMS2(model, sequences, charges, settings);
      for (std::size_t i = begin; i < end; i += 2)
      {
        std::vector<FragmentIon> ions[2];
        for (std::size_t m = 0; m < 2; ++m)
        {
          const IonLadder ladder(AASequence::fromString(pairs[i + m].sequence));
          pairs[i + m].mz = ladder.precursor(pairs[i + m].charge);
          ions[m] = predictedFragments(ladder, pairs[i + m].charge, ms2[i + m - begin], fragments.min_mz, fragments.max_mz, fragments.count);
        }
        if (ions[0].size() < fragments.min_count || ions[1].size() < fragments.min_count) continue;
        for (std::size_t m = 0; m < 2; ++m)
        {
          const auto& r = pairs[i + m];
          const std::string id = std::to_string(r.id);
          // TRAML_ID from the target's sequence, as in the survey library.
          const std::string target_traml = pairs[i].sequence + "/" + std::to_string(pairs[i].charge);
          const std::string traml = m == 0 ? target_traml : "DECOY_" + target_traml;
          ParquetLibraryWriter::Precursor out;
          out.id = id;
          out.traml_id = traml;
          out.sequence = r.sequence;
          out.mz = r.mz;
          out.rt = r.rt;
          out.im = r.im;
          out.charge = r.charge;
          out.decoy = m == 1;
          out.proteins = &r.proteins;
          scoring.add(out, rows(r.id, ions[m], fragments.count));
          search.add(out, rows(r.id, ions[m], SEARCH_FRAGMENTS));
        }
        ++stats.pairs_out;
      }
      say(log, "FASTA final library: " + std::to_string(end / 2) + "/" + std::to_string(stats.pairs_in) + " pairs predicted, seconds=" +
               std::to_string(since(start)));
    }
    stats.fragments = scoring.fragments();
    scoring.close();
    search.close();
    say(log, "FASTA final library: " + path + " pairs=" + std::to_string(stats.pairs_out) + " (dropped below " +
             std::to_string(fragments.min_count) + " fragments: " + std::to_string(stats.pairs_in - stats.pairs_out) + ") fragments=" +
             std::to_string(stats.fragments) + " seconds=" + std::to_string(since(start)));
    return stats;
#endif
  }
}
