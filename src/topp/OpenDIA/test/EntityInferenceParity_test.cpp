// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $

#include <results/Quantification.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace OpenMS::OpenDIA;

namespace
{
std::string sequence(unsigned n)
{
  const std::string alphabet = "ACDEFGHIKLNPQRSTVWY";
  std::string s = "M";
  for (unsigned j = 0; j < 4; ++j)
  {
    s += alphabet[n % alphabet.size()];
    n /= alphabet.size();
  }
  return s + "PEPTK";
}
void require(bool pass, const char* message)
{
  if (!pass) throw std::runtime_error(message);
}
}

int main(int argc, char** argv)
{
  try
  {
    require(argc == 2, "Supply output directory");
    const std::filesystem::path out(argv[1]);
    std::filesystem::create_directories(out);
    AnalysisLibrary library;
    std::vector<RunResult> runs(3);
    std::ofstream precursors(out / "precursors.tsv"), windows(out / "windows.tsv");
    precursors << "precursor\tsequence\tdecoy\tcharge\tproteins\n";
    windows << "run\tprecursor\tscore\tarea\n" << std::setprecision(17);
    for (unsigned i = 0; i < 480; ++i)
    {
      for (unsigned d = 0; d < 2; ++d)
      {
        for (unsigned charge = 2; charge <= (i % 9 == 0 ? 3U : 2U); ++charge)
        {
          AnalysisLibrary::Precursor p;
          // Identical strings in opposite classes deliberately test class-qualified identity.
          p.sequence = sequence(i / 2);
          if (i % 2 == 1) p.sequence.insert(1, "(Oxidation)");
          p.decoy = d;
          p.charge = charge;
          p.pair = static_cast<std::uint32_t>(library.pairs++);
          p.decoy_partner = UINT32_MAX; // Inference must not require paired counterparts.
          // Some accessions collide across classes too. Class must remain part of identity.
          const std::string prefix = d && i % 17 != 0 ? "DECOY_" : "";
          std::string protein = prefix + "P" + std::to_string(i / 2);
          if (i % 23 == 0) protein += ";" + (d ? std::string("DECOY_") : std::string()) + "SHARED";
          p.proteins = static_cast<std::uint32_t>(library.protein_sets.size());
          library.protein_sets.push_back(protein);
          const auto index = static_cast<std::uint32_t>(library.precursors.size());
          library.precursors.push_back(p);
          precursors << index << '\t' << p.sequence << '\t' << d << '\t' << charge << '\t' << protein << '\n';
          for (unsigned r = 0; r < runs.size(); ++r)
          {
            if ((i + r * 3) % 31 == 0 || (i == 479 && d == 1)) continue;
            // Reproducible heterogeneous distributions with explicit tied score blocks.
            double z = 0;
            for (unsigned k = 0; k < 6; ++k)
              z += ((i * (37 + 12 * k) + r * (53 + k * 8) + d * 71 + charge * 11) % 367) / 367.0;
            const float score = static_cast<float>(std::round((z - 3 + (d ? 0 : 0.65 + .10 * r) + .2 * r) * 16) / 16 + 5);
            for (unsigned candidate = 0; candidate < 2; ++candidate)
            {
              RunWindow w;
              w.precursor = index;
              w.score = score - candidate * .5f;
              w.q_value = 0;
              w.peak.area = 10 + i + charge + candidate * 10000;
              runs[r].windows.push_back(w);
              windows << r << '\t' << index << '\t' << w.score << '\t' << w.peak.area << '\n';
            }
          }
        }
      }
    }
    const auto result = quantify(library, runs, {1, 1, 1});
    std::ofstream aggregate(out / "native.tsv");
    aggregate << "level\tcontext\trun\tentity\tdecoy\tscore\tp\tq\tpep\n" << std::setprecision(17);
    for (unsigned level = 0; level < 2; ++level)
    {
      const auto& entities = level ? result.protein_groups : result.peptides;
      const auto& mapping = level ? result.precursor_protein_group : result.precursor_peptide;
      for (std::size_t e = 0; e < entities.size(); ++e)
      {
        const auto& entity = entities[e];
        std::string key;
        if (level)
        {
          require(entity.proteins.size() == 1, "No fabricated indistinguishable protein groups");
          key = entity.proteins.front();
        }
        else
        {
          auto p = std::find(mapping.begin(), mapping.end(), e);
          require(p != mapping.end(), "Entity missing precursor mapping");
          key = library.precursors[static_cast<std::size_t>(p - mapping.begin())].sequence;
        }
        double global_score = -std::numeric_limits<double>::infinity();
        for (const auto& run : entity.runs)
        {
          global_score = std::max(global_score, run.score);
          aggregate << (level ? "protein" : "peptide") << "\trun-specific\t" << run.run << '\t' << key << '\t' << entity.decoy << '\t' << run.score << '\t' << run.p_value << '\t' << run.q_value << '\t' << run.pep << '\n';
          aggregate << (level ? "protein" : "peptide") << "\texperiment-wide\t" << run.run << '\t' << key << '\t' << entity.decoy << '\t' << run.score << '\t' << run.experiment_p_value << '\t' << run.experiment_q_value << '\t' << run.experiment_pep << '\n';
        }
        if (!entity.runs.empty())
          aggregate << (level ? "protein" : "peptide") << "\tglobal\t-1\t" << key << '\t' << entity.decoy << '\t' << global_score << '\t' << entity.global_p_value << '\t' << entity.global_q_value << '\t' << entity.global_pep << '\n';
      }
    }
    // Degenerate input must fail closed with warnings, rather than inventing estimates.
    for (auto& run : runs) for (auto& w : run.windows) w.score = 1;
    const auto degenerate = quantify(library, runs, {1, 1, 1});
    require(!degenerate.inference_warnings.empty(), "Unestimable contexts must warn");
    for (const auto& entity : degenerate.peptides)
    {
      require(std::isnan(entity.global_q_value), "Unestimable global q must be NA");
      for (const auto& run : entity.runs)
        require(std::isnan(run.q_value) && std::isnan(run.p_value) && std::isnan(run.pep), "Unestimable run metrics must be NA");
    }
    std::cout << "Wrote independent PyProphet parity fixture to " << out << '\n';
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
