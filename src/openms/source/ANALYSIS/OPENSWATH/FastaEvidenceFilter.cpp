// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/FastaEvidenceFilter.h>

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionListEvidenceFilter.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <boost/math/distributions/chi_squared.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <future>
#include <fstream>
#include <limits>
#include <map>
#include <iomanip>
#include <mutex>
#include <numeric>
#include <sstream>
#include <regex>
#include <set>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
  namespace
  {
    bool withinSwathWindow_(double precursor_mz, const OpenSwath::SwathMap& map)
    {
      return precursor_mz >= map.lower && precursor_mz <= map.upper;
    }

    bool parseIonAnnotation_(const String& annotation, std::string& product_type, int& ordinal)
    {
      static const std::regex ion_pattern("^([A-Za-z\\.']+)([0-9]+)");
      std::smatch match;
      const std::string annotation_std = annotation.c_str();
      if (!std::regex_search(annotation_std, match, ion_pattern))
      {
        return false;
      }
      product_type = match[1].str();
      ordinal = std::stoi(match[2].str());
      return true;
    }

    double halfTolerance_(double full_window, bool /*ppm*/, double fallback)
    {
      if (full_window > 0.0)
      {
        return std::max(full_window * 0.5, fallback);
      }
      return fallback;
    }

    std::uint64_t makePeptideChargeKey_(Size peptide_index, std::uint16_t precursor_charge)
    {
      return (static_cast<std::uint64_t>(peptide_index) << 16) |
             static_cast<std::uint64_t>(precursor_charge);
    }

    int resolveThreadCount_(int threads)
    {
#ifdef _OPENMP
      const int max_threads = std::max(1, omp_get_max_threads());
      if (threads < 0)
      {
        return max_threads;
      }
      if (threads == 0)
      {
        return 1;
      }
      return std::max(1, std::min(threads, max_threads));
#else
      (void)threads;
      return 1;
#endif
    }

    struct MzCoverageInterval
    {
      double lower{0.0};
      double upper{0.0};
    };

    std::vector<MzCoverageInterval> buildSwathMzCoverage_(const std::vector<FastaEvidenceFilter::RunData>& runs,
                                                          double min_upper_edge_dist)
    {
      std::vector<MzCoverageInterval> intervals;
      for (const auto& run : runs)
      {
        for (const auto& swath_map : run.swath_maps)
        {
          if (swath_map.ms1)
          {
            continue;
          }

          const double adjusted_upper = swath_map.upper - min_upper_edge_dist;
          if (swath_map.lower >= adjusted_upper)
          {
            continue;
          }
          intervals.push_back({swath_map.lower, adjusted_upper});
        }
      }

      if (intervals.empty())
      {
        return intervals;
      }

      std::sort(intervals.begin(), intervals.end(),
                [](const MzCoverageInterval& lhs, const MzCoverageInterval& rhs)
                {
                  if (lhs.lower != rhs.lower) return lhs.lower < rhs.lower;
                  return lhs.upper < rhs.upper;
                });

      std::vector<MzCoverageInterval> merged;
      merged.reserve(intervals.size());
      for (const auto& interval : intervals)
      {
        if (merged.empty() || interval.lower > merged.back().upper)
        {
          merged.push_back(interval);
        }
        else
        {
          merged.back().upper = std::max(merged.back().upper, interval.upper);
        }
      }
      return merged;
    }

    std::vector<FastaEvidenceFilter::PeptideEntry> filterPeptidesByMzCoverage_(
      const std::vector<FastaEvidenceFilter::PeptideEntry>& peptides,
      const std::vector<MzCoverageInterval>& coverage)
    {
      if (coverage.empty())
      {
        return {};
      }

      std::vector<FastaEvidenceFilter::PeptideEntry> filtered;
      filtered.reserve(peptides.size());
      Size interval_index = 0;
      for (const auto& peptide : peptides)
      {
        while (interval_index < coverage.size() &&
               peptide.precursor_mz > coverage[interval_index].upper)
        {
          ++interval_index;
        }
        if (interval_index == coverage.size())
        {
          break;
        }
        if (peptide.precursor_mz > coverage[interval_index].lower &&
            peptide.precursor_mz <= coverage[interval_index].upper)
        {
          filtered.push_back(peptide);
        }
      }
      return filtered;
    }

    std::string serializeProteinGenePairs_(const FastaEvidenceFilter::PeptideEntry& peptide)
    {
      const auto sanitize_metadata_token = [](const std::string& value) -> std::string
      {
        std::string sanitized;
        sanitized.reserve(value.size());
        for (const unsigned char c : value)
        {
          if (c == '\r' || c == '\n' || c == '\t')
          {
            sanitized.push_back(' ');
          }
          else
          {
            sanitized.push_back(static_cast<char>(c));
          }
        }

        const Size begin = sanitized.find_first_not_of(' ');
        if (begin == std::string::npos)
        {
          return {};
        }
        const Size end = sanitized.find_last_not_of(' ');
        return sanitized.substr(begin, end - begin + 1);
      };

      std::vector<std::string> pairs;
      pairs.reserve(peptide.protein_refs.size());
      for (const auto& protein_ref : peptide.protein_refs)
      {
        auto gene_name_it = peptide.protein_gene_names_by_accession.find(protein_ref);
        const std::string gene_name =
          gene_name_it != peptide.protein_gene_names_by_accession.end() ?
          gene_name_it->second : std::string();
        pairs.push_back(protein_ref + "=" + sanitize_metadata_token(gene_name));
      }
      return ListUtils::concatenate(pairs, ";").c_str();
    }

    void parseProteinGenePairs_(const std::string& serialized_pairs,
                                std::vector<std::string>& protein_refs,
                                std::map<std::string, std::string>& gene_names_by_accession)
    {
      protein_refs.clear();
      gene_names_by_accession.clear();
      if (serialized_pairs.empty())
      {
        return;
      }

      Size begin = 0;
      while (begin <= serialized_pairs.size())
      {
        const Size end = serialized_pairs.find(';', begin);
        const std::string token = serialized_pairs.substr(
          begin, end == std::string::npos ? std::string::npos : end - begin);
        if (!token.empty())
        {
          const Size separator = token.find('=');
          const std::string accession = token.substr(0, separator);
          const std::string gene_name =
            separator == std::string::npos ? std::string() : token.substr(separator + 1);
          protein_refs.push_back(accession);
          if (!accession.empty())
          {
            gene_names_by_accession[accession] = gene_name;
          }
        }

        if (end == std::string::npos)
        {
          break;
        }
        begin = end + 1;
      }
    }

    void appendShardEntryBuffer_(std::string& buffer, const FastaEvidenceFilter::PeptideEntry& peptide)
    {
      std::ostringstream line_buffer;
      line_buffer << std::quoted(peptide.canonical_key) << ' '
                  << std::quoted(peptide.internal_key) << ' '
                  << std::quoted(peptide.peptide_sequence) << ' '
                  << std::quoted(peptide.modified_peptide_sequence) << ' '
                  << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << peptide.precursor_mz << ' '
                  << peptide.precursor_charge << ' '
                  << std::quoted(serializeProteinGenePairs_(peptide))
                  << '\n';
      buffer += line_buffer.str();
    }

    void mergeShardEntryLine_(std::map<std::string, FastaEvidenceFilter::PeptideEntry>& peptide_map,
                              const std::string& line,
                              const String& shard_path,
                              Size line_number)
    {
      if (line.empty())
      {
        return;
      }

      std::istringstream line_stream(line);
      std::string canonical_key;
      std::string internal_key;
      std::string peptide_sequence;
      std::string modified_peptide_sequence;
      std::string serialized_pairs;
      double precursor_mz = 0.0;
      int precursor_charge = 0;
      if (!(line_stream >> std::quoted(canonical_key) >>
            std::quoted(internal_key) >>
            std::quoted(peptide_sequence) >>
            std::quoted(modified_peptide_sequence) >>
            precursor_mz >>
            precursor_charge >>
            std::quoted(serialized_pairs)))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    line.c_str(),
                                    "Invalid sharded peptide entry at line " +
                                    String(line_number) + " in '" + shard_path + "'.");
      }

      auto [it, inserted] = peptide_map.emplace(internal_key, FastaEvidenceFilter::PeptideEntry{});
      auto& peptide = it->second;
      if (inserted)
      {
        peptide.canonical_key = canonical_key;
        peptide.internal_key = internal_key;
        peptide.peptide_sequence = peptide_sequence;
        peptide.modified_peptide_sequence = modified_peptide_sequence;
        peptide.precursor_mz = precursor_mz;
        peptide.precursor_charge = precursor_charge;
      }

      std::vector<std::string> protein_refs;
      std::map<std::string, std::string> gene_names_by_accession;
      parseProteinGenePairs_(serialized_pairs, protein_refs, gene_names_by_accession);
      peptide.protein_refs.insert(peptide.protein_refs.end(), protein_refs.begin(), protein_refs.end());
      for (const auto& gene_name_item : gene_names_by_accession)
      {
        std::string& stored_gene_name = peptide.protein_gene_names_by_accession[gene_name_item.first];
        if (stored_gene_name.empty())
        {
          stored_gene_name = gene_name_item.second;
        }
      }
    }

    std::vector<FastaEvidenceFilter::PeptideEntry> loadMergedShardEntries_(const String& shard_path)
    {
      std::ifstream input(shard_path.c_str());
      if (!input.good())
      {
        return {};
      }

      std::map<std::string, FastaEvidenceFilter::PeptideEntry> peptide_map;
      std::string line;
      Size line_number = 0;
      while (std::getline(input, line))
      {
        ++line_number;
        mergeShardEntryLine_(peptide_map, line, shard_path, line_number);
      }

      std::vector<FastaEvidenceFilter::PeptideEntry> peptides;
      peptides.reserve(peptide_map.size());
      for (auto& item : peptide_map)
      {
        auto& refs = item.second.protein_refs;
        std::sort(refs.begin(), refs.end());
        refs.erase(std::unique(refs.begin(), refs.end()), refs.end());
        peptides.push_back(std::move(item.second));
      }

      std::sort(peptides.begin(), peptides.end(),
                [](const FastaEvidenceFilter::PeptideEntry& lhs,
                   const FastaEvidenceFilter::PeptideEntry& rhs)
                {
                  if (lhs.precursor_mz != rhs.precursor_mz) return lhs.precursor_mz < rhs.precursor_mz;
                  if (lhs.precursor_charge != rhs.precursor_charge) return lhs.precursor_charge < rhs.precursor_charge;
                  return lhs.modified_peptide_sequence < rhs.modified_peptide_sequence;
                });
      return peptides;
    }

    std::string formatRetentionRatio_(Size retained, Size total)
    {
      std::ostringstream os;
      os << retained << " of " << total;
      if (total > 0)
      {
        const double retained_pct = 100.0 * static_cast<double>(retained) / static_cast<double>(total);
        const double reduced_pct = 100.0 - retained_pct;
        os << " (" << std::fixed << std::setprecision(1) << retained_pct
           << "%, reduced " << reduced_pct << "%)";
      }
      return os.str();
    }

    std::string extractGeneName_(const FASTAFile::FASTAEntry& entry)
    {
      const auto extract_from_text = [](const String& text) -> std::string
      {
        const Size gene_pos = text.find("GN=");
        if (gene_pos == String::npos)
        {
          return {};
        }

        const Size gene_begin = gene_pos + 3;
        const auto whitespace_it = std::find_if(text.begin() + static_cast<SignedSize>(gene_begin),
                                                text.end(),
                                                [](const char c)
                                                {
                                                  return std::isspace(static_cast<unsigned char>(c)) != 0;
                                                });
        const Size gene_end = whitespace_it == text.end() ?
                              text.size() :
                              static_cast<Size>(std::distance(text.begin(), whitespace_it));

        return text.substr(gene_begin, gene_end - gene_begin).c_str();
      };

      std::string gene_name = extract_from_text(entry.description);
      if (gene_name.empty())
      {
        gene_name = extract_from_text(entry.identifier);
      }
      return gene_name;
    }

    constexpr Size STAGE2_AGGREGATED_RUN_SCORES = 3;
    constexpr Size STAGE2_UNASSIGNED_RUN_ID = std::numeric_limits<Size>::max();
    constexpr std::array<double, STAGE2_AGGREGATED_RUN_SCORES> STAGE2_TOP_RUN_WEIGHTS{
      1.0, 0.75, 0.5};

    struct Stage2CandidateStats
    {
      Size best_matched_ions{0};
      Size supporting_spectra{0};
      std::uint64_t supporting_run_mask{0};
      Size strong_supporting_spectra{0};
      std::uint64_t strong_supporting_run_mask{0};
      std::uint64_t streak_ge_2_run_mask{0};
      std::uint64_t streak_ge_3_run_mask{0};
      double best_spectrum_matched_intensity_fraction{0.0};
      Size best_spectrum_matched_b_ions{0};
      Size best_spectrum_matched_y_ions{0};
      Size best_spectrum_longest_b_run{0};
      Size best_spectrum_longest_y_run{0};
      double best_spectrum_longest_y_pct{0.0};
      double best_spectrum_poisson_proxy{0.0};
      double best_spectrum_score{0.0};
      Size best_run_streak_length{0};
      double best_run_streak_score{0.0};
      std::array<double, STAGE2_AGGREGATED_RUN_SCORES> top_run_scores{};
      std::array<Size, STAGE2_AGGREGATED_RUN_SCORES> top_run_streak_lengths{};
      std::array<Size, STAGE2_AGGREGATED_RUN_SCORES> top_run_ids{
        STAGE2_UNASSIGNED_RUN_ID, STAGE2_UNASSIGNED_RUN_ID, STAGE2_UNASSIGNED_RUN_ID};
      std::string best_source_file;
      std::string best_native_spectrum_id;
    };

    struct Stage2RunStreakState
    {
      Size last_supported_spectrum_index{std::numeric_limits<Size>::max()};
      Size current_streak_length{0};
      Size best_streak_length{0};
      double best_spectrum_score{0.0};
      double best_streak_score{0.0};
    };

    struct Stage2ScoredSpectrumCandidate
    {
      Size candidate_id{0};
      FragmentIndex::SpectrumMatch match;
      double matched_intensity_fraction{0.0};
      double spectrum_score{0.0};
    };

    struct Stage2LowerOrderObservation
    {
      Size candidate_id{0};
      Size run_index{0};
      std::uint16_t precursor_charge{0};
      Size rank{0};
      Size matched_ions{0};
      double matched_intensity_fraction{0.0};
      double spectrum_score{0.0};
      std::string source_file;
      std::string native_spectrum_id;
      bool used_for_scoring{false};
      bool used_for_null{false};
      double local_pvalue{-1.0};
    };

    bool betterStage2SpectrumMatch_(const FragmentIndex::SpectrumMatch& lhs,
                                    const FragmentIndex::SpectrumMatch& rhs)
    {
      if (lhs.num_matched_ != rhs.num_matched_)
      {
        return lhs.num_matched_ > rhs.num_matched_;
      }
      if (lhs.matched_intensity_sum_ != rhs.matched_intensity_sum_)
      {
        return lhs.matched_intensity_sum_ > rhs.matched_intensity_sum_;
      }

      const auto abs_iso_lhs = lhs.isotope_error_ < 0 ? -lhs.isotope_error_ : lhs.isotope_error_;
      const auto abs_iso_rhs = rhs.isotope_error_ < 0 ? -rhs.isotope_error_ : rhs.isotope_error_;
      if (abs_iso_lhs != abs_iso_rhs)
      {
        return abs_iso_lhs < abs_iso_rhs;
      }
      if (lhs.isotope_error_ != rhs.isotope_error_)
      {
        return lhs.isotope_error_ < rhs.isotope_error_;
      }
      return lhs.precursor_charge_ < rhs.precursor_charge_;
    }

    inline void setStage2OrdinalBit_(std::array<std::uint64_t, 2>& words, Size ordinal)
    {
      if (ordinal == 0 || ordinal > 128)
      {
        return;
      }
      const Size zero_based = ordinal - 1;
      words[zero_based / 64] |= (std::uint64_t{1} << (zero_based % 64));
    }

    Size popcountStage2OrdinalBits_(const std::array<std::uint64_t, 2>& words)
    {
      return static_cast<Size>(std::popcount(words[0]) + std::popcount(words[1]));
    }

    Size longestStage2OrdinalRun_(const std::array<std::uint64_t, 2>& words)
    {
      Size longest_run = 0;
      Size current_run = 0;
      for (Size ordinal_index = 0; ordinal_index < 128; ++ordinal_index)
      {
        const std::uint64_t word = words[ordinal_index / 64];
        const bool is_set = (word & (std::uint64_t{1} << (ordinal_index % 64))) != 0;
        if (is_set)
        {
          ++current_run;
          longest_run = std::max(longest_run, current_run);
        }
        else
        {
          current_run = 0;
        }
      }
      return longest_run;
    }

    double computeStage2PoissonProxy_(Size observed_matches,
                                      Size theoretical_fragments,
                                      Size spectrum_peak_count,
                                      double fragment_tolerance_da,
                                      double spectrum_mz_span)
    {
      if (observed_matches == 0 || theoretical_fragments == 0 || spectrum_peak_count == 0)
      {
        return 0.0;
      }

      const double bounded_mz_span = std::max(1.0, spectrum_mz_span);
      const double bounded_fragment_tolerance = std::max(1e-6, fragment_tolerance_da);
      const double expected_match_probability = std::clamp(
        (2.0 * bounded_fragment_tolerance * static_cast<double>(spectrum_peak_count)) / bounded_mz_span,
        1e-9, 0.95);
      const double lambda = std::max(
        1e-9, static_cast<double>(theoretical_fragments) * expected_match_probability);
      const double observed = static_cast<double>(observed_matches);
      if (observed <= lambda)
      {
        return 0.0;
      }
      return observed * std::log(observed / lambda) - (observed - lambda);
    }

    double computeStage2SpectrumScore_(const FragmentIndex::SpectrumMatch& match,
                                       double matched_intensity_fraction)
    {
      const double bounded_fraction = std::clamp(matched_intensity_fraction, 0.0, 1.0);
      return static_cast<double>(match.num_matched_) +
             0.5 * std::log1p(1000.0 * bounded_fraction);
    }

    void computeStage2SpectrumDiagnostics_(const FragmentIndex::SpectrumMatch& match,
                                           Size theoretical_b_ions,
                                           Size theoretical_y_ions,
                                           Size spectrum_peak_count,
                                           double spectrum_mz_span,
                                           double fragment_tolerance_da,
                                           double& poisson_proxy,
                                           double& longest_y_pct,
                                           Size& matched_b_ions,
                                           Size& matched_y_ions,
                                           Size& longest_b_run,
                                           Size& longest_y_run)
    {
      matched_b_ions = popcountStage2OrdinalBits_(match.matched_b_ordinal_words_);
      matched_y_ions = popcountStage2OrdinalBits_(match.matched_y_ordinal_words_);
      longest_b_run = longestStage2OrdinalRun_(match.matched_b_ordinal_words_);
      longest_y_run = longestStage2OrdinalRun_(match.matched_y_ordinal_words_);
      longest_y_pct =
        theoretical_y_ions > 0 ?
        static_cast<double>(longest_y_run) / static_cast<double>(theoretical_y_ions) :
        0.0;

      poisson_proxy = computeStage2PoissonProxy_(
        matched_b_ions + matched_y_ions,
        theoretical_b_ions + theoretical_y_ions,
        spectrum_peak_count,
        fragment_tolerance_da,
        spectrum_mz_span);
    }

    double computeStage2RunStreakScore_(double best_spectrum_score, Size streak_length)
    {
      if (best_spectrum_score <= 0.0)
      {
        return 0.0;
      }
      return best_spectrum_score + 0.35 * std::log1p(static_cast<double>(streak_length));
    }

    void updateStage2RunStreak_(Stage2RunStreakState& state,
                                Size spectrum_index,
                                double spectrum_score)
    {
      if (state.current_streak_length > 0 &&
          state.last_supported_spectrum_index + 1 == spectrum_index)
      {
        ++state.current_streak_length;
      }
      else
      {
        state.current_streak_length = 1;
      }
      state.last_supported_spectrum_index = spectrum_index;
      state.best_streak_length = std::max(state.best_streak_length, state.current_streak_length);
      state.best_spectrum_score = std::max(state.best_spectrum_score, spectrum_score);

      const double streak_score =
        computeStage2RunStreakScore_(state.best_spectrum_score, state.best_streak_length);
      state.best_streak_score = std::max(state.best_streak_score, streak_score);
    }

    void sortStage2RunSummaries_(std::array<Size, STAGE2_AGGREGATED_RUN_SCORES>& run_ids,
                                 std::array<double, STAGE2_AGGREGATED_RUN_SCORES>& run_scores,
                                 std::array<Size, STAGE2_AGGREGATED_RUN_SCORES>& run_streak_lengths)
    {
      for (Size outer_idx = 0; outer_idx < run_scores.size(); ++outer_idx)
      {
        for (Size inner_idx = outer_idx + 1; inner_idx < run_scores.size(); ++inner_idx)
        {
          const bool swap_entries =
            run_scores[inner_idx] > run_scores[outer_idx] ||
            (run_scores[inner_idx] == run_scores[outer_idx] &&
             run_streak_lengths[inner_idx] > run_streak_lengths[outer_idx]) ||
            (run_scores[inner_idx] == run_scores[outer_idx] &&
             run_streak_lengths[inner_idx] == run_streak_lengths[outer_idx] &&
             run_ids[outer_idx] == STAGE2_UNASSIGNED_RUN_ID &&
             run_ids[inner_idx] != STAGE2_UNASSIGNED_RUN_ID);
          if (!swap_entries)
          {
            continue;
          }
          std::swap(run_scores[outer_idx], run_scores[inner_idx]);
          std::swap(run_streak_lengths[outer_idx], run_streak_lengths[inner_idx]);
          std::swap(run_ids[outer_idx], run_ids[inner_idx]);
        }
      }
    }

    void upsertStage2RunSummary_(Stage2CandidateStats& stats,
                                 Size run_index,
                                 double run_score,
                                 Size run_streak_length)
    {
      if (run_index < 64)
      {
        stats.supporting_run_mask |= (std::uint64_t{1} << run_index);
        if (run_streak_length >= 2)
        {
          stats.streak_ge_2_run_mask |= (std::uint64_t{1} << run_index);
        }
        if (run_streak_length >= 3)
        {
          stats.streak_ge_3_run_mask |= (std::uint64_t{1} << run_index);
        }
      }

      if (run_score > stats.best_run_streak_score ||
          (run_score == stats.best_run_streak_score &&
           run_streak_length > stats.best_run_streak_length))
      {
        stats.best_run_streak_score = run_score;
        stats.best_run_streak_length = run_streak_length;
      }

      for (Size run_slot = 0; run_slot < stats.top_run_ids.size(); ++run_slot)
      {
        if (stats.top_run_ids[run_slot] != run_index)
        {
          continue;
        }
        stats.top_run_scores[run_slot] = std::max(stats.top_run_scores[run_slot], run_score);
        stats.top_run_streak_lengths[run_slot] =
          std::max(stats.top_run_streak_lengths[run_slot], run_streak_length);
        sortStage2RunSummaries_(stats.top_run_ids, stats.top_run_scores, stats.top_run_streak_lengths);
        return;
      }

      for (Size run_slot = 0; run_slot < stats.top_run_ids.size(); ++run_slot)
      {
        if (stats.top_run_ids[run_slot] != STAGE2_UNASSIGNED_RUN_ID)
        {
          continue;
        }
        stats.top_run_ids[run_slot] = run_index;
        stats.top_run_scores[run_slot] = run_score;
        stats.top_run_streak_lengths[run_slot] = run_streak_length;
        sortStage2RunSummaries_(stats.top_run_ids, stats.top_run_scores, stats.top_run_streak_lengths);
        return;
      }

      if (run_score > stats.top_run_scores.back())
      {
        stats.top_run_ids.back() = run_index;
        stats.top_run_scores.back() = run_score;
        stats.top_run_streak_lengths.back() = run_streak_length;
        sortStage2RunSummaries_(stats.top_run_ids, stats.top_run_scores, stats.top_run_streak_lengths);
      }
    }

    Size countStage2SupportingRuns_(const Stage2CandidateStats& stats)
    {
      const Size masked_support = static_cast<Size>(std::popcount(stats.supporting_run_mask));
      if (masked_support > 0)
      {
        return masked_support;
      }

      Size tracked_support = 0;
      for (const Size run_id : stats.top_run_ids)
      {
        if (run_id != STAGE2_UNASSIGNED_RUN_ID)
        {
          ++tracked_support;
        }
      }
      return tracked_support;
    }

    Size countStage2StrongSupportingRuns_(const Stage2CandidateStats& stats)
    {
      return static_cast<Size>(std::popcount(stats.strong_supporting_run_mask));
    }

    Size countStage2RunsWithMinStreak_(const Stage2CandidateStats& stats, Size min_streak_length)
    {
      if (min_streak_length <= 1)
      {
        return countStage2SupportingRuns_(stats);
      }

      const std::uint64_t streak_mask =
        min_streak_length <= 2 ? stats.streak_ge_2_run_mask :
        stats.streak_ge_3_run_mask;
      const Size masked_support = static_cast<Size>(std::popcount(streak_mask));
      if (masked_support > 0)
      {
        return masked_support;
      }

      Size tracked_support = 0;
      for (Size run_slot = 0; run_slot < stats.top_run_ids.size(); ++run_slot)
      {
        if (stats.top_run_ids[run_slot] != STAGE2_UNASSIGNED_RUN_ID &&
            stats.top_run_streak_lengths[run_slot] >= min_streak_length)
        {
          ++tracked_support;
        }
      }
      return tracked_support;
    }

    double composeStage2PeptideScore_(const Stage2CandidateStats& stats)
    {
      double score = 0.0;
      for (Size i = 0; i < stats.top_run_scores.size(); ++i)
      {
        score += STAGE2_TOP_RUN_WEIGHTS[i] * stats.top_run_scores[i];
      }
      score += 1.5 * static_cast<double>(countStage2StrongSupportingRuns_(stats));
      score += 0.35 * std::log1p(static_cast<double>(stats.strong_supporting_spectra));
      score += 0.1 * std::log1p(static_cast<double>(countStage2SupportingRuns_(stats)));
      return score;
    }

    double empiricalStage2TailPValue_(const std::vector<double>& sorted_null_scores,
                                      double observed_score)
    {
      if (sorted_null_scores.empty())
      {
        return 1.0;
      }

      const auto first_ge = std::lower_bound(
        sorted_null_scores.begin(), sorted_null_scores.end(), observed_score);
      const Size exceed_count = static_cast<Size>(sorted_null_scores.end() - first_ge);
      return (static_cast<double>(exceed_count) + 1.0) /
             (static_cast<double>(sorted_null_scores.size()) + 1.0);
    }

    double combineStage2RunPValues_(const std::vector<double>& run_pvalues)
    {
      if (run_pvalues.empty())
      {
        return 1.0;
      }
      if (run_pvalues.size() == 1)
      {
        return std::clamp(run_pvalues.front(), 0.0, 1.0);
      }

      double fisher_statistic = 0.0;
      for (double pvalue : run_pvalues)
      {
        fisher_statistic += -2.0 * std::log(std::clamp(pvalue, 1e-300, 1.0));
      }

      const boost::math::chi_squared fisher_null(2.0 * static_cast<double>(run_pvalues.size()));
      const double combined_pvalue = boost::math::cdf(boost::math::complement(fisher_null, fisher_statistic));
      if (!std::isfinite(combined_pvalue))
      {
        return 1.0;
      }
      return std::clamp(combined_pvalue, 0.0, 1.0);
    }

    void mergeStage2CandidateStats_(Stage2CandidateStats& destination,
                                    const Stage2CandidateStats& source)
    {
      destination.best_matched_ions = std::max(destination.best_matched_ions, source.best_matched_ions);
      destination.supporting_spectra += source.supporting_spectra;
      destination.supporting_run_mask |= source.supporting_run_mask;
      destination.strong_supporting_spectra += source.strong_supporting_spectra;
      destination.strong_supporting_run_mask |= source.strong_supporting_run_mask;
      destination.streak_ge_2_run_mask |= source.streak_ge_2_run_mask;
      destination.streak_ge_3_run_mask |= source.streak_ge_3_run_mask;
      if (source.best_spectrum_score > destination.best_spectrum_score)
      {
        destination.best_spectrum_matched_intensity_fraction = source.best_spectrum_matched_intensity_fraction;
        destination.best_spectrum_matched_b_ions = source.best_spectrum_matched_b_ions;
        destination.best_spectrum_matched_y_ions = source.best_spectrum_matched_y_ions;
        destination.best_spectrum_longest_b_run = source.best_spectrum_longest_b_run;
        destination.best_spectrum_longest_y_run = source.best_spectrum_longest_y_run;
        destination.best_spectrum_longest_y_pct = source.best_spectrum_longest_y_pct;
        destination.best_spectrum_poisson_proxy = source.best_spectrum_poisson_proxy;
        destination.best_spectrum_score = source.best_spectrum_score;
        destination.best_source_file = source.best_source_file;
        destination.best_native_spectrum_id = source.best_native_spectrum_id;
      }
      if (source.best_run_streak_score > destination.best_run_streak_score ||
          (source.best_run_streak_score == destination.best_run_streak_score &&
           source.best_run_streak_length > destination.best_run_streak_length))
      {
        destination.best_run_streak_score = source.best_run_streak_score;
        destination.best_run_streak_length = source.best_run_streak_length;
      }
      for (Size run_slot = 0; run_slot < source.top_run_ids.size(); ++run_slot)
      {
        if (source.top_run_ids[run_slot] == STAGE2_UNASSIGNED_RUN_ID)
        {
          continue;
        }
        upsertStage2RunSummary_(
          destination,
          source.top_run_ids[run_slot],
          source.top_run_scores[run_slot],
          source.top_run_streak_lengths[run_slot]);
      }
    }
  }

  FastaEvidenceFilter::FastaEvidenceFilter() :
    DefaultParamHandler("FastaEvidenceFilter"),
    ProgressLogger()
  {
    defaults_.setValue("aggregation_method", "any",
                       "How to combine stage-1 evidence across multiple DIA runs.");
    defaults_.setValidStrings("aggregation_method", {"any", "all"});

    TransitionListEvidenceFilter stage1_filter;
    Param stage1_defaults = stage1_filter.getParameters();
    stage1_defaults.remove("enabled");
    defaults_.insert("Stage1:", stage1_defaults);
    defaults_.setValue("Stage1:precursor_batch_size", 50000,
                       "Maximum number of precursors materialized in one stage-1 transition batch.");
    defaults_.setMinInt("Stage1:precursor_batch_size", 1);
    defaults_.setValue("Stage1:max_concurrent_runs", 0,
                       "Maximum number of DIA runs filtered concurrently for one stage-1 batch. 0 uses the automatic limit derived from -threads.");
    defaults_.setMinInt("Stage1:max_concurrent_runs", 0);

    defaults_.setValue("Stage2:mode", "lower_order_null",
                       "How to accept stage-2 confirmed peptides.");
    defaults_.setValidStrings("Stage2:mode", {"raw_score", "lower_order_null"});
    defaults_.setValue("Stage2:max_qvalue", 0.01,
                       "Maximum peptide-level q-value in Stage2:mode=lower_order_null.");
    defaults_.setMinFloat("Stage2:max_qvalue", 0.0);
    defaults_.setMaxFloat("Stage2:max_qvalue", 1.0);
    defaults_.setValue("Stage2:min_matched_ions", 5,
                       "Minimum matched fragment ions required in stage 2.");
    defaults_.setMinInt("Stage2:min_matched_ions", 0);
    defaults_.setValue("Stage2:strong_min_matched_ions", 6,
                       "Minimum matched fragment ions needed for a spectrum to count as strong stage-2 support.");
    defaults_.setMinInt("Stage2:strong_min_matched_ions", 0);
    defaults_.setValue("Stage2:strong_min_intensity_fraction", 0.05,
                       "Minimum matched-intensity fraction needed for a spectrum to count as strong stage-2 support.");
    defaults_.setMinFloat("Stage2:strong_min_intensity_fraction", 0.0);
    defaults_.setMaxFloat("Stage2:strong_min_intensity_fraction", 1.0);
    defaults_.setValue("Stage2:lower_order_min_rank", 5,
                       "Lowest spectrum-rank position contributing to the decoy-free lower-order null model.");
    defaults_.setMinInt("Stage2:lower_order_min_rank", 2);
    defaults_.setValue("Stage2:lower_order_max_rank", 10,
                       "Highest spectrum-rank position contributing to the decoy-free lower-order null model.");
    defaults_.setMinInt("Stage2:lower_order_max_rank", 2);
    defaults_.setValue("Stage2:lower_order_scored_ranks", 3,
                       "Maximum number of top-ranked spectrum candidates per spectrum kept for decoy-free lower-order scoring.");
    defaults_.setMinInt("Stage2:lower_order_scored_ranks", 1);
    defaults_.setValue("Stage2:lower_order_min_null_scores", 256,
                       "Minimum number of lower-order null scores needed before a charge-specific null model is used.");
    defaults_.setMinInt("Stage2:lower_order_min_null_scores", 1);
    defaults_.setValue("Protein:min_confirmed_peptides", 1,
                       "Minimum number of confirmed peptides needed to keep a protein.");
    defaults_.setMinInt("Protein:min_confirmed_peptides", 1);
    defaults_.setValue("Protein:unique_peptides_only", "false",
                       "If true, only peptides mapping to one protein count toward protein support.");
    defaults_.setValidStrings("Protein:unique_peptides_only", {"true", "false"});

    defaults_.setValue("Export:export_fragments", "false",
                       "If true, emit one TSV row per precursor-fragment pair instead of one row per precursor.");
    defaults_.setValidStrings("Export:export_fragments", {"true", "false"});
    defaults_.setValue("Export:export_stage2_scores", "false",
                       "If true, keep scored Stage-2 export rows for optional TSV export.");
    defaults_.setValidStrings("Export:export_stage2_scores", {"true", "false"});
    defaults_.setValue("Export:modified_sequence_format", "unimod_accession",
                       "Format used when exporting modified peptide sequences to TSV outputs.");
    defaults_.setValidStrings("Export:modified_sequence_format", {"unimod_accession", "codename"});

    std::vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    std::vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);

    defaults_.setValue("SearchSpace:enzyme", "Trypsin", "The enzyme used for in-silico digestion.");
    defaults_.setValidStrings("SearchSpace:enzyme", ListUtils::create<std::string>(all_enzymes));
    defaults_.setValue("SearchSpace:specificity", "full",
                       "Required enzyme specificity at the peptide termini.");
    defaults_.setValidStrings("SearchSpace:specificity", {"full", "semi", "none"});
    defaults_.setValue("SearchSpace:missed_cleavages", 1,
                       "Number of missed cleavages allowed during digestion.");
    defaults_.setMinInt("SearchSpace:missed_cleavages", 0);
    defaults_.setValue("SearchSpace:min_size", 7,
                       "Minimum peptide length after digestion.");
    defaults_.setMinInt("SearchSpace:min_size", 1);
    defaults_.setValue("SearchSpace:max_size", 40,
                       "Maximum peptide length after digestion.");
    defaults_.setMinInt("SearchSpace:max_size", 1);
    defaults_.setValue("SearchSpace:min_mass", 100,
                       "Minimum peptide mass kept in the search space.");
    defaults_.setMinInt("SearchSpace:min_mass", 0);
    defaults_.setValue("SearchSpace:max_mass", 9000,
                       "Maximum peptide mass kept in the search space.");
    defaults_.setMinInt("SearchSpace:max_mass", 0);
    defaults_.setValue("SearchSpace:modifications:fixed",
                       std::vector<std::string>{"Carbamidomethyl (C)"},
                       "Fixed modifications specified with search-engine names.");
    defaults_.setValidStrings("SearchSpace:modifications:fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("SearchSpace:modifications:variable",
                       std::vector<std::string>{"Oxidation (M)"},
                       "Variable modifications specified with search-engine names.");
    defaults_.setValidStrings("SearchSpace:modifications:variable", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("SearchSpace:modifications:variable_max_per_peptide", 2,
                       "Maximum number of variable modifications per peptide candidate.");
    defaults_.setMinInt("SearchSpace:modifications:variable_max_per_peptide", 0);
    defaults_.setValue("SearchSpace:precursor:min_charge", 2,
                       "Minimum precursor charge exported to the peptide search space.");
    defaults_.setMinInt("SearchSpace:precursor:min_charge", 1);
    defaults_.setValue("SearchSpace:precursor:max_charge", 5,
                       "Maximum precursor charge exported to the peptide search space.");
    defaults_.setMinInt("SearchSpace:precursor:max_charge", 1);
    defaults_.setValue("SearchSpace:fragment:min_charge", 1,
                       "Minimum fragment charge generated for each precursor.");
    defaults_.setMinInt("SearchSpace:fragment:min_charge", 1);
    defaults_.setValue("SearchSpace:fragment:max_charge", 2,
                       "Maximum fragment charge generated for each precursor.");
    defaults_.setMinInt("SearchSpace:fragment:max_charge", 1);
    defaults_.setValue("SearchSpace:fragment:max_per_precursor", 6,
                       "Maximum number of theoretical fragments kept per precursor.");
    defaults_.setMinInt("SearchSpace:fragment:max_per_precursor", 1);
    defaults_.setValue("SearchSpace:fragment:min_ion_index", 2,
                       "Skip theoretical fragment ions with ordinal less than or equal to this value.");
    defaults_.setMinInt("SearchSpace:fragment:min_ion_index", 0);
    defaults_.setValue("SearchSpace:sharding:max_proteins_per_chunk", 1000,
                       "Maximum number of FASTA proteins materialized per in-memory peptide-generation chunk before sharding to disk. Set to 0 to disable sharding.");
    defaults_.setMinInt("SearchSpace:sharding:max_proteins_per_chunk", 0);
    defaults_.setValue("SearchSpace:sharding:num_shards", 128,
                       "Number of on-disk precursor shards used when sharded Stage-1 search-space generation is enabled.");
    defaults_.setMinInt("SearchSpace:sharding:num_shards", 1);
    defaults_.setValue("SearchSpace:sharding:temp_directory", File::getTempDirectory(),
                       "Base directory used for temporary shard files when sharded Stage-1 search-space generation is enabled.");
    defaults_.setValue("SearchSpace:sharding:keep_temporary_files", "false",
                       "If true, keep temporary sharded Stage-1 files instead of removing them automatically.");
    defaults_.setValidStrings("SearchSpace:sharding:keep_temporary_files", {"true", "false"});

    defaultsToParam_();
    updateMembers_();
  }

  void FastaEvidenceFilter::updateMembers_()
  {
    aggregation_method_ = param_.getValue("aggregation_method").toString();

    stage1_min_supported_precursors_ = static_cast<Size>(param_.getValue("Stage1:min_supported_precursors"));
    stage1_precursor_batch_size_ = static_cast<Size>(param_.getValue("Stage1:precursor_batch_size"));
    stage1_max_concurrent_runs_ = static_cast<Size>(param_.getValue("Stage1:max_concurrent_runs"));

    stage2_mode_ = param_.getValue("Stage2:mode").toString();
    stage2_max_qvalue_ = static_cast<double>(param_.getValue("Stage2:max_qvalue"));
    stage2_min_matched_ions_ = static_cast<Int>(param_.getValue("Stage2:min_matched_ions"));
    stage2_strong_min_matched_ions_ = static_cast<Int>(param_.getValue("Stage2:strong_min_matched_ions"));
    stage2_strong_min_intensity_fraction_ = static_cast<double>(param_.getValue("Stage2:strong_min_intensity_fraction"));
    stage2_lower_order_min_rank_ = static_cast<Int>(param_.getValue("Stage2:lower_order_min_rank"));
    stage2_lower_order_max_rank_ = static_cast<Int>(param_.getValue("Stage2:lower_order_max_rank"));
    stage2_lower_order_scored_ranks_ = static_cast<Int>(param_.getValue("Stage2:lower_order_scored_ranks"));
    stage2_lower_order_min_null_scores_ = static_cast<Int>(param_.getValue("Stage2:lower_order_min_null_scores"));

    protein_min_confirmed_peptides_ = static_cast<Size>(param_.getValue("Protein:min_confirmed_peptides"));
    protein_unique_peptides_only_ = param_.getValue("Protein:unique_peptides_only").toString() == "true";

    export_fragments_ = param_.getValue("Export:export_fragments").toString() == "true";
    export_stage2_scores_ = param_.getValue("Export:export_stage2_scores").toString() == "true";

    enzyme_ = param_.getValue("SearchSpace:enzyme").toString();
    enzyme_specificity_ = param_.getValue("SearchSpace:specificity").toString();
    peptide_missed_cleavages_ = static_cast<Int>(param_.getValue("SearchSpace:missed_cleavages"));
    peptide_min_size_ = static_cast<Int>(param_.getValue("SearchSpace:min_size"));
    peptide_max_size_ = static_cast<Int>(param_.getValue("SearchSpace:max_size"));
    peptide_min_mass_ = static_cast<Int>(param_.getValue("SearchSpace:min_mass"));
    peptide_max_mass_ = static_cast<Int>(param_.getValue("SearchSpace:max_mass"));
    modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("SearchSpace:modifications:fixed"));
    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("SearchSpace:modifications:variable"));
    max_variable_mods_per_peptide_ = static_cast<Int>(param_.getValue("SearchSpace:modifications:variable_max_per_peptide"));
    precursor_min_charge_ = static_cast<Int>(param_.getValue("SearchSpace:precursor:min_charge"));
    precursor_max_charge_ = static_cast<Int>(param_.getValue("SearchSpace:precursor:max_charge"));
    fragment_min_charge_ = static_cast<Int>(param_.getValue("SearchSpace:fragment:min_charge"));
    fragment_max_charge_ = static_cast<Int>(param_.getValue("SearchSpace:fragment:max_charge"));
    max_fragments_per_precursor_ = static_cast<Size>(param_.getValue("SearchSpace:fragment:max_per_precursor"));
    fragment_min_ion_index_ = static_cast<Int>(param_.getValue("SearchSpace:fragment:min_ion_index"));
    search_space_max_proteins_per_chunk_ = static_cast<Size>(param_.getValue("SearchSpace:sharding:max_proteins_per_chunk"));
    search_space_num_shards_ = static_cast<Size>(param_.getValue("SearchSpace:sharding:num_shards"));
    search_space_sharding_temp_directory_ =
      File::absolutePath(param_.getValue("SearchSpace:sharding:temp_directory").toString()).ensureLastChar('/');
    search_space_sharding_keep_temporary_files_ =
      param_.getValue("SearchSpace:sharding:keep_temporary_files").toString() == "true";
  }

  std::string FastaEvidenceFilter::makeCanonicalPeptideKey(const std::string& modified_peptide_sequence, int precursor_charge)
  {
    return modified_peptide_sequence + "/" + std::to_string(precursor_charge);
  }

  bool FastaEvidenceFilter::hasDecoyPrefix_(const std::string& value, const std::string& decoy_prefix)
  {
    if (!decoy_prefix.empty() && value.rfind(decoy_prefix, 0) == 0)
    {
      return true;
    }
    return value.rfind("DECOY_", 0) == 0 || value.rfind("Decoy_", 0) == 0 || value.rfind("decoy_", 0) == 0;
  }

  std::string FastaEvidenceFilter::buildInternalKey_(const std::string& modified_peptide_sequence,
                                                     int precursor_charge,
                                                     bool decoy) const
  {
    const std::string canonical_key = makeCanonicalPeptideKey(modified_peptide_sequence, precursor_charge);
    if (!decoy)
    {
      return canonical_key;
    }
    return std::string(stage2_decoy_prefix_.c_str()) + canonical_key;
  }

  std::vector<FastaEvidenceFilter::FragmentRecord> FastaEvidenceFilter::buildTheoreticalFragments_(const AASequence& modified_sequence,
                                                                                                   int precursor_charge) const
  {
    std::vector<FragmentRecord> fragments;

    const int max_fragment_charge = std::min(fragment_max_charge_, precursor_charge);
    if (max_fragment_charge < fragment_min_charge_)
    {
      return fragments;
    }

    TheoreticalSpectrumGenerator spectrum_generator;
    Param tsg_params = spectrum_generator.getParameters();
    tsg_params.setValue("add_metainfo", "true");
    tsg_params.setValue("add_losses", "false");
    tsg_params.setValue("add_precursor_peaks", "false");
    tsg_params.setValue("isotope_model", "none");
    tsg_params.setValue("add_abundant_immonium_ions", "false");
    tsg_params.setValue("add_a_ions", "false");
    tsg_params.setValue("add_b_ions", "true");
    tsg_params.setValue("add_c_ions", "false");
    tsg_params.setValue("add_x_ions", "false");
    tsg_params.setValue("add_y_ions", "true");
    tsg_params.setValue("add_z_ions", "false");
    spectrum_generator.setParameters(tsg_params);

    MSSpectrum theoretical_spectrum;
    spectrum_generator.getSpectrum(theoretical_spectrum, modified_sequence, fragment_min_charge_, max_fragment_charge, precursor_charge);

    if (theoretical_spectrum.getStringDataArrays().empty() || theoretical_spectrum.getIntegerDataArrays().empty())
    {
      return fragments;
    }

    const auto& annotations = theoretical_spectrum.getStringDataArrays()[0];
    const auto& charges = theoretical_spectrum.getIntegerDataArrays()[0];
    const Size n = std::min(theoretical_spectrum.size(), std::min(annotations.size(), charges.size()));
    for (Size i = 0; i < n; ++i)
    {
      std::string product_type;
      int ordinal = 0;
      if (!parseIonAnnotation_(annotations[i], product_type, ordinal))
      {
        continue;
      }
      if (ordinal <= fragment_min_ion_index_)
      {
        continue;
      }

      FragmentRecord fragment;
      fragment.product_mz = theoretical_spectrum[i].getMZ();
      fragment.product_charge = charges[i];
      fragment.product_type = product_type;
      fragment.product_ordinal = ordinal;
      fragments.push_back(std::move(fragment));
    }

    std::sort(fragments.begin(), fragments.end(),
              [](const FragmentRecord& lhs, const FragmentRecord& rhs)
              {
                if (lhs.product_charge != rhs.product_charge) return lhs.product_charge < rhs.product_charge;
                if (lhs.product_type != rhs.product_type) return lhs.product_type < rhs.product_type;
                if (lhs.product_ordinal != rhs.product_ordinal) return lhs.product_ordinal < rhs.product_ordinal;
                return lhs.product_mz < rhs.product_mz;
              });

    if (fragments.size() > max_fragments_per_precursor_)
    {
      fragments.resize(max_fragments_per_precursor_);
    }

    return fragments;
  }

  std::vector<FastaEvidenceFilter::PeptideEntry> FastaEvidenceFilter::generatePeptideEntries(const std::vector<FASTAFile::FASTAEntry>& fasta_entries) const
  {
    return generatePeptideEntries_(fasta_entries, true);
  }

  std::vector<FastaEvidenceFilter::PeptideEntry> FastaEvidenceFilter::generatePeptideEntries_(const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                                                                               bool include_fragments) const
  {
    FragmentIndex fragment_index;
    Param fragment_index_params = fragment_index.getParameters();
    fragment_index_params.setValue("enzyme", enzyme_);
    fragment_index_params.setValue("peptide:enzyme_specificity", enzyme_specificity_);
    fragment_index_params.setValue("peptide:missed_cleavages", peptide_missed_cleavages_);
    fragment_index_params.setValue("peptide:min_size", peptide_min_size_);
    fragment_index_params.setValue("peptide:max_size", peptide_max_size_);
    fragment_index_params.setValue("peptide:min_mass", peptide_min_mass_);
    fragment_index_params.setValue("peptide:max_mass", peptide_max_mass_);
    fragment_index_params.setValue("modifications:fixed", ListUtils::create<std::string>(modifications_fixed_));
    fragment_index_params.setValue("modifications:variable", ListUtils::create<std::string>(modifications_variable_));
    fragment_index_params.setValue("modifications:variable_max_per_peptide", max_variable_mods_per_peptide_);
    fragment_index_params.setValue("precursor:min_charge", precursor_min_charge_);
    fragment_index_params.setValue("precursor:max_charge", precursor_max_charge_);
    fragment_index_params.setValue("fragment:max_charge", fragment_max_charge_);
    fragment_index_params.setValue("fragment:min_ion_index", fragment_min_ion_index_);
    fragment_index_params.setValue("fragment:min_mz", 0);
    fragment_index_params.setValue("fragment:max_mz", 2000);
    fragment_index.setParameters(fragment_index_params);
    fragment_index.buildPeptidesOnly(fasta_entries);

    std::map<std::string, PeptideEntry> peptide_map;
    for (const auto& peptide : fragment_index.getPeptides())
    {
      const AASequence modified_sequence = fragment_index.reconstructModifiedSequence(peptide, fasta_entries);
      const std::string modified_string = modified_sequence.toString().c_str();
      const std::string unmodified_string = modified_sequence.toUnmodifiedString().c_str();
      const std::string protein_ref = fasta_entries[peptide.protein_idx].identifier.c_str();
      const std::string gene_name = extractGeneName_(fasta_entries[peptide.protein_idx]);
      const bool decoy = hasDecoyPrefix_(protein_ref, stage2_decoy_prefix_.c_str());

      for (int precursor_charge = precursor_min_charge_; precursor_charge <= precursor_max_charge_; ++precursor_charge)
      {
        const std::string internal_key = buildInternalKey_(modified_string, precursor_charge, decoy);
        auto [it, inserted] = peptide_map.emplace(internal_key, PeptideEntry{});
        PeptideEntry& entry = it->second;
        if (inserted)
        {
          entry.internal_key = internal_key;
          entry.canonical_key = makeCanonicalPeptideKey(modified_string, precursor_charge);
          entry.peptide_sequence = unmodified_string;
          entry.modified_peptide_sequence = modified_string;
          entry.precursor_charge = precursor_charge;
          entry.precursor_mz = modified_sequence.getMZ(precursor_charge);
          if (include_fragments)
          {
            entry.fragments = buildTheoreticalFragments_(modified_sequence, precursor_charge);
          }
          entry.decoy = decoy;
        }
        entry.protein_refs.push_back(protein_ref);
        std::string& stored_gene_name = entry.protein_gene_names_by_accession[protein_ref];
        if (stored_gene_name.empty())
        {
          stored_gene_name = gene_name;
        }
      }
    }

    std::vector<PeptideEntry> peptides;
    peptides.reserve(peptide_map.size());
    for (auto& item : peptide_map)
    {
      auto& refs = item.second.protein_refs;
      std::sort(refs.begin(), refs.end());
      refs.erase(std::unique(refs.begin(), refs.end()), refs.end());
      peptides.push_back(std::move(item.second));
    }

    std::sort(peptides.begin(), peptides.end(),
              [](const PeptideEntry& lhs, const PeptideEntry& rhs)
              {
                if (lhs.precursor_mz != rhs.precursor_mz) return lhs.precursor_mz < rhs.precursor_mz;
                if (lhs.precursor_charge != rhs.precursor_charge) return lhs.precursor_charge < rhs.precursor_charge;
                return lhs.modified_peptide_sequence < rhs.modified_peptide_sequence;
              });

    return peptides;
  }

  OpenSwath::LightTargetedExperiment FastaEvidenceFilter::buildStage1Experiment_(const std::vector<PeptideEntry>& peptides,
                                                                                 Size begin_idx,
                                                                                 Size end_idx,
                                                                                 int threads) const
  {
    OpenSwath::LightTargetedExperiment experiment;
    end_idx = std::min(end_idx, peptides.size());
    if (begin_idx >= end_idx)
    {
      return experiment;
    }

    const Size peptide_count = end_idx - begin_idx;
    experiment.compounds.reserve(peptide_count);
    experiment.transitions.reserve(peptide_count * max_fragments_per_precursor_);

    struct Stage1Chunk
    {
      std::vector<OpenSwath::LightCompound> compounds;
      std::vector<OpenSwath::LightTransition> transitions;
    };

    const int thread_count = std::min<int>(resolveThreadCount_(threads), static_cast<int>(peptide_count));
    if (thread_count <= 1 || peptide_count < 1024)
    {
      Stage1Chunk chunk;
      chunk.compounds.reserve(peptide_count);
      chunk.transitions.reserve(peptide_count * max_fragments_per_precursor_);
      for (Size peptide_idx = begin_idx; peptide_idx < end_idx; ++peptide_idx)
      {
        const auto& peptide = peptides[peptide_idx];
        OpenSwath::LightCompound compound;
        compound.id = peptide.canonical_key;
        compound.sequence = peptide.modified_peptide_sequence;
        compound.charge = peptide.precursor_charge;
        chunk.compounds.push_back(std::move(compound));

        std::vector<FragmentRecord> generated_fragments;
        const std::vector<FragmentRecord>* fragments = &peptide.fragments;
        if (fragments->empty())
        {
          generated_fragments = buildTheoreticalFragments_(AASequence::fromString(peptide.modified_peptide_sequence),
                                                           peptide.precursor_charge);
          fragments = &generated_fragments;
        }

        Size intensity_rank = fragments->size();
        for (const auto& fragment : *fragments)
        {
          OpenSwath::LightTransition transition;
          transition.transition_name = peptide.canonical_key + "_" + fragment.product_type + std::to_string(fragment.product_ordinal) + "_" + std::to_string(fragment.product_charge);
          transition.peptide_ref = peptide.canonical_key;
          transition.precursor_mz = peptide.precursor_mz;
          transition.product_mz = fragment.product_mz;
          transition.fragment_charge = static_cast<int8_t>(fragment.product_charge);
          transition.fragment_nr = static_cast<int16_t>(fragment.product_ordinal);
          transition.setFragmentType(fragment.product_type);
          transition.setLibraryIntensity(static_cast<double>(intensity_rank--));
          transition.setDetectingTransition(true);
          transition.setQuantifyingTransition(true);
          transition.setDecoy(false);
          chunk.transitions.push_back(std::move(transition));
        }
      }
      experiment.compounds = std::move(chunk.compounds);
      experiment.transitions = std::move(chunk.transitions);
      return experiment;
    }

    std::vector<Stage1Chunk> chunks(static_cast<Size>(thread_count));
#ifdef _OPENMP
    #pragma omp parallel num_threads(thread_count)
#endif
    {
#ifdef _OPENMP
      const int thread_id = omp_get_thread_num();
#else
      const int thread_id = 0;
#endif
      const Size local_begin = begin_idx + (peptide_count * static_cast<Size>(thread_id)) / static_cast<Size>(thread_count);
      const Size local_end = begin_idx + (peptide_count * static_cast<Size>(thread_id + 1)) / static_cast<Size>(thread_count);
      auto& chunk = chunks[static_cast<Size>(thread_id)];
      const Size local_count = local_end - local_begin;
      chunk.compounds.reserve(local_count);
      chunk.transitions.reserve(local_count * max_fragments_per_precursor_);

      for (Size peptide_idx = local_begin; peptide_idx < local_end; ++peptide_idx)
      {
        const auto& peptide = peptides[peptide_idx];
        OpenSwath::LightCompound compound;
        compound.id = peptide.canonical_key;
        compound.sequence = peptide.modified_peptide_sequence;
        compound.charge = peptide.precursor_charge;
        chunk.compounds.push_back(std::move(compound));

        std::vector<FragmentRecord> generated_fragments;
        const std::vector<FragmentRecord>* fragments = &peptide.fragments;
        if (fragments->empty())
        {
          generated_fragments = buildTheoreticalFragments_(AASequence::fromString(peptide.modified_peptide_sequence),
                                                           peptide.precursor_charge);
          fragments = &generated_fragments;
        }

        Size intensity_rank = fragments->size();
        for (const auto& fragment : *fragments)
        {
          OpenSwath::LightTransition transition;
          transition.transition_name = peptide.canonical_key + "_" + fragment.product_type + std::to_string(fragment.product_ordinal) + "_" + std::to_string(fragment.product_charge);
          transition.peptide_ref = peptide.canonical_key;
          transition.precursor_mz = peptide.precursor_mz;
          transition.product_mz = fragment.product_mz;
          transition.fragment_charge = static_cast<int8_t>(fragment.product_charge);
          transition.fragment_nr = static_cast<int16_t>(fragment.product_ordinal);
          transition.setFragmentType(fragment.product_type);
          transition.setLibraryIntensity(static_cast<double>(intensity_rank--));
          transition.setDetectingTransition(true);
          transition.setQuantifyingTransition(true);
          transition.setDecoy(false);
          chunk.transitions.push_back(std::move(transition));
        }
      }
    }

    for (auto& chunk : chunks)
    {
      experiment.compounds.insert(experiment.compounds.end(),
                                  std::make_move_iterator(chunk.compounds.begin()),
                                  std::make_move_iterator(chunk.compounds.end()));
      experiment.transitions.insert(experiment.transitions.end(),
                                    std::make_move_iterator(chunk.transitions.begin()),
                                    std::make_move_iterator(chunk.transitions.end()));
    }

    return experiment;
  }

  void FastaEvidenceFilter::populateFragments_(std::vector<PeptideEntry>& peptides) const
  {
    for (auto& peptide : peptides)
    {
      if (!peptide.fragments.empty())
      {
        continue;
      }
      peptide.fragments = buildTheoreticalFragments_(AASequence::fromString(peptide.modified_peptide_sequence),
                                                     peptide.precursor_charge);
    }
  }

  std::vector<FASTAFile::FASTAEntry> FastaEvidenceFilter::buildReducedTargetFasta_(const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                                                                   const std::vector<PeptideEntry>& supported_peptides) const
  {
    std::unordered_set<std::string> accessions;
    for (const auto& peptide : supported_peptides)
    {
      accessions.insert(peptide.protein_refs.begin(), peptide.protein_refs.end());
    }

    std::vector<FASTAFile::FASTAEntry> reduced_fasta;
    reduced_fasta.reserve(accessions.size());
    for (const auto& entry : fasta_entries)
    {
      const std::string accession = entry.identifier.c_str();
      if (accessions.find(accession) != accessions.end())
      {
        reduced_fasta.push_back(entry);
      }
    }
    return reduced_fasta;
  }

  Param FastaEvidenceFilter::buildFragmentIndexParams_(const ChromExtractParams& ms1_params,
                                                       const ChromExtractParams& ms2_params) const
  {
    FragmentIndex fragment_index;
    Param params = fragment_index.getParameters();

    params.setValue("enzyme", enzyme_);
    params.setValue("peptide:enzyme_specificity", enzyme_specificity_);
    params.setValue("peptide:missed_cleavages", peptide_missed_cleavages_);
    params.setValue("peptide:min_size", peptide_min_size_);
    params.setValue("peptide:max_size", peptide_max_size_);
    params.setValue("peptide:min_mass", peptide_min_mass_);
    params.setValue("peptide:max_mass", peptide_max_mass_);
    params.setValue("modifications:fixed", ListUtils::create<std::string>(modifications_fixed_));
    params.setValue("modifications:variable", ListUtils::create<std::string>(modifications_variable_));
    params.setValue("modifications:variable_max_per_peptide", max_variable_mods_per_peptide_);
    params.setValue("precursor:min_charge", precursor_min_charge_);
    params.setValue("precursor:max_charge", precursor_max_charge_);
    params.setValue("fragment:max_charge", fragment_max_charge_);
    params.setValue("fragment:min_ion_index", fragment_min_ion_index_);
    params.setValue("fragment:min_mz", 0);
    params.setValue("fragment:max_mz", 2000);
    params.setValue("fragment:min_matched_ions", stage2_min_matched_ions_);
    params.setValue("precursor:isotope_error_min", 0);
    params.setValue("precursor:isotope_error_max", 0);
    params.setValue("scoring:max_candidates_per_spectrum", 500);
    params.setValue("decoys", "false");
    params.setValue("snes_enabled", "false");

    params.setValue("precursor:mass_tolerance_unit", ms1_params.ppm ? "ppm" : "Da");
    params.setValue("precursor:mass_tolerance_lower", halfTolerance_(ms1_params.mz_extraction_window, ms1_params.ppm, ms1_params.ppm ? 10.0 : 0.05));
    params.setValue("precursor:mass_tolerance_upper", halfTolerance_(ms1_params.mz_extraction_window, ms1_params.ppm, ms1_params.ppm ? 10.0 : 0.05));
    params.setValue("fragment:mass_tolerance_unit", ms2_params.ppm ? "ppm" : "Da");
    params.setValue("fragment:mass_tolerance", halfTolerance_(ms2_params.mz_extraction_window, ms2_params.ppm, ms2_params.ppm ? 10.0 : 0.05));

    return params;
  }

  FastaEvidenceFilter::Stage2ScoreBundle FastaEvidenceFilter::scoreStage2_(const std::vector<RunData>& runs,
                                                                           const std::vector<PeptideEntry>& candidates,
                                                                           const ChromExtractParams& ms1_params,
                                                                           const ChromExtractParams& ms2_params,
                                                                           int threads) const
  {
    Stage2ScoreBundle bundle;
    if (candidates.empty())
    {
      return bundle;
    }

    FragmentIndex fragment_index;
    Param fragment_index_params = buildFragmentIndexParams_(ms1_params, ms2_params);
    fragment_index.setParameters(fragment_index_params);
    std::vector<FragmentIndex::ExplicitPeptide> explicit_candidates;
    explicit_candidates.reserve(candidates.size());
    for (Size candidate_id = 0; candidate_id < candidates.size(); ++candidate_id)
    {
      explicit_candidates.push_back({
        candidates[candidate_id].modified_peptide_sequence,
        static_cast<float>(
          candidates[candidate_id].precursor_mz * static_cast<double>(candidates[candidate_id].precursor_charge) -
          (static_cast<double>(candidates[candidate_id].precursor_charge - 1) * Constants::PROTON_MASS_U)),
        static_cast<UInt32>(candidate_id)});
    }
    fragment_index.buildFromPeptideSequences(explicit_candidates);
    const auto& index_peptides = fragment_index.getPeptides();

    struct Stage2MapJob
    {
      const RunData* run{nullptr};
      const OpenSwath::SwathMap* swath_map{nullptr};
      Size run_index{0};
      std::string source_file;
      std::vector<FragmentIndex::PrecursorRangeQuery> precursor_queries;
      std::unordered_set<Size> allowed_candidate_ids;
      Size spectrum_count{0};
      Size candidate_count{0};
      Size query_units_per_spectrum{0};
    };

    struct Stage2CandidateRef
    {
      const PeptideEntry* candidate{nullptr};
      Size candidate_id{0};
    };

    std::vector<Stage2CandidateRef> sorted_candidates;
    sorted_candidates.reserve(candidates.size());
    std::vector<Size> theoretical_b_ions_by_candidate(candidates.size(), 0);
    std::vector<Size> theoretical_y_ions_by_candidate(candidates.size(), 0);
    for (Size candidate_id = 0; candidate_id < candidates.size(); ++candidate_id)
    {
      const Size peptide_length = candidates[candidate_id].peptide_sequence.size();
      const Size total_terminal_ions =
        peptide_length > static_cast<Size>(fragment_min_ion_index_ + 1) ?
        peptide_length - static_cast<Size>(fragment_min_ion_index_ + 1) :
        0;
      theoretical_b_ions_by_candidate[candidate_id] = total_terminal_ions;
      theoretical_y_ions_by_candidate[candidate_id] = total_terminal_ions;
      sorted_candidates.push_back({&candidates[candidate_id], candidate_id});
    }
    std::sort(sorted_candidates.begin(), sorted_candidates.end(),
              [](const Stage2CandidateRef& lhs, const Stage2CandidateRef& rhs)
              {
                return lhs.candidate->precursor_mz < rhs.candidate->precursor_mz;
              });

    std::vector<Stage2MapJob> jobs;
    Size total_spectra = 0;
    SignedSize total_query_units = 0;
    for (Size run_index = 0; run_index < runs.size(); ++run_index)
    {
      const auto& run = runs[run_index];
      for (Size swath_map_index = 0; swath_map_index < run.swath_maps.size(); ++swath_map_index)
      {
        const auto& swath_map = run.swath_maps[swath_map_index];
        if (swath_map.ms1 || !swath_map.sptr)
        {
          continue;
        }

        const auto lower_it = std::lower_bound(
          sorted_candidates.begin(), sorted_candidates.end(), swath_map.lower,
          [](const Stage2CandidateRef& candidate_ref, double lower_bound)
          {
            return candidate_ref.candidate->precursor_mz < lower_bound;
          });
        const auto upper_it = std::upper_bound(
          sorted_candidates.begin(), sorted_candidates.end(), swath_map.upper,
          [](double upper_bound, const Stage2CandidateRef& candidate_ref)
          {
            return upper_bound < candidate_ref.candidate->precursor_mz;
          });
        if (lower_it == upper_it)
        {
          continue;
        }

        const Size n_spectra = swath_map.sptr->getNrSpectra();
        if (n_spectra == 0)
        {
          continue;
        }

        Stage2MapJob job;
        job.run = &run;
        job.swath_map = &swath_map;
        job.run_index = run_index;
        if (swath_map_index < run.swath_map_sources.size())
        {
          job.source_file = run.swath_map_sources[swath_map_index].c_str();
        }
        job.spectrum_count = n_spectra;
        job.allowed_candidate_ids.reserve(static_cast<Size>(std::distance(lower_it, upper_it)));
        std::map<std::uint16_t, std::pair<float, float>> precursor_mass_bounds_by_charge;
        for (auto candidate_it = lower_it; candidate_it != upper_it; ++candidate_it)
        {
          const auto* candidate = candidate_it->candidate;
          job.allowed_candidate_ids.insert(candidate_it->candidate_id);

          const std::uint16_t precursor_charge = static_cast<std::uint16_t>(candidate->precursor_charge);
          const float precursor_mass = static_cast<float>(
            candidate->precursor_mz * static_cast<double>(precursor_charge) -
            (static_cast<double>(precursor_charge - 1) * Constants::PROTON_MASS_U));
          auto range_it = precursor_mass_bounds_by_charge.find(precursor_charge);
          if (range_it == precursor_mass_bounds_by_charge.end())
          {
            precursor_mass_bounds_by_charge.emplace(precursor_charge,
                                                   std::make_pair(precursor_mass, precursor_mass));
          }
          else
          {
            range_it->second.first = std::min(range_it->second.first, precursor_mass);
            range_it->second.second = std::max(range_it->second.second, precursor_mass);
          }
        }

        job.candidate_count = job.allowed_candidate_ids.size();
        job.precursor_queries.reserve(precursor_mass_bounds_by_charge.size());
        for (const auto& item : precursor_mass_bounds_by_charge)
        {
          FragmentIndex::PrecursorRangeQuery query;
          query.precursor_mass_lower = item.second.first;
          query.precursor_mass_upper = item.second.second;
          query.precursor_charge = item.first;
          job.precursor_queries.push_back(query);
        }
        job.query_units_per_spectrum = job.precursor_queries.size();
        total_spectra += n_spectra;
        total_query_units += static_cast<SignedSize>(n_spectra * job.query_units_per_spectrum);
        jobs.push_back(std::move(job));
      }
    }

    if (jobs.empty())
    {
      return bundle;
    }

    const int thread_count = std::max(1, threads);
    const SignedSize total_progress = std::max<SignedSize>(1, total_query_units);
    const SignedSize progress_step = std::max<SignedSize>(1, total_progress / 100);
    const Size per_thread_reserve = std::min<Size>(
      std::max<Size>(4096, candidates.size() / (static_cast<Size>(thread_count) * 16) + 1),
      static_cast<Size>(65536));
    const bool use_lower_order_null = stage2_mode_ == "lower_order_null";
    const Size lower_order_min_rank = static_cast<Size>(std::max<Int>(2, stage2_lower_order_min_rank_));
    const Size lower_order_max_rank = static_cast<Size>(std::max<Int>(
      stage2_lower_order_max_rank_, stage2_lower_order_min_rank_));
    const Size lower_order_scored_ranks = static_cast<Size>(std::max<Int>(1, stage2_lower_order_scored_ranks_));
    const Size lower_order_min_null_scores = static_cast<Size>(std::max<Int>(1, stage2_lower_order_min_null_scores_));
    const Size per_thread_observation_reserve = std::min<Size>(
      std::max<Size>(
        1024,
        ((total_spectra / static_cast<Size>(thread_count)) + 1) * lower_order_scored_ranks),
      static_cast<Size>(262144));

    OPENMS_LOG_INFO << "Stage 2: scoring " << jobs.size() << " SWATH maps across "
                    << total_spectra << " spectra and " << total_progress
                    << " precursor-range queries with up to " << thread_count
                    << " threads." << std::endl;
    startProgress(0, total_progress, "Stage 2: scoring fragment-index candidates");

    std::atomic<SignedSize> processed_queries{0};
    std::atomic<SignedSize> next_progress_update{progress_step};
    std::mutex progress_mutex;
    std::vector<std::unordered_map<Size, Stage2CandidateStats>> local_candidate_stats(static_cast<Size>(thread_count));
    std::vector<std::unordered_map<int, std::vector<double>>> local_lower_order_null_scores_by_charge(
      static_cast<Size>(thread_count));
    std::vector<std::vector<Stage2LowerOrderObservation>> local_lower_order_observations(
      static_cast<Size>(thread_count));
    for (auto& local_scores : local_candidate_stats)
    {
      local_scores.reserve(per_thread_reserve);
    }
    if (use_lower_order_null)
    {
      for (auto& observations : local_lower_order_observations)
      {
        observations.reserve(per_thread_observation_reserve);
      }
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) num_threads(thread_count)
#endif
    for (SignedSize job_index = 0; job_index < static_cast<SignedSize>(jobs.size()); ++job_index)
    {
#ifdef _OPENMP
      const int thread_id = omp_get_thread_num();
#else
      const int thread_id = 0;
#endif
      auto& local_scores = local_candidate_stats[static_cast<Size>(thread_id)];
      auto& local_null_scores_by_charge = local_lower_order_null_scores_by_charge[static_cast<Size>(thread_id)];
      auto& local_observations = local_lower_order_observations[static_cast<Size>(thread_id)];
      std::unordered_map<Size, Stage2RunStreakState> run_streak_states;
      const auto& job = jobs[static_cast<Size>(job_index)];
      FragmentIndex::SpectrumMatchesTopN matches;
      if (export_stage2_scores_)
      {
        run_streak_states.reserve(std::min(job.candidate_count, per_thread_reserve));
      }

      for (Size spectrum_index = 0; spectrum_index < job.spectrum_count; ++spectrum_index)
      {
        OpenSwath::SpectrumPtr spectrum_ptr;
        if (job.run->pasef && job.swath_map->imLower >= 0.0 && job.swath_map->imUpper >= 0.0)
        {
          spectrum_ptr = job.swath_map->sptr->getSpectrumById(static_cast<int>(spectrum_index), job.swath_map->imLower, job.swath_map->imUpper);
        }
        else
        {
          spectrum_ptr = job.swath_map->sptr->getSpectrumById(static_cast<int>(spectrum_index));
        }

        if (spectrum_ptr)
        {
          std::string native_spectrum_id;
          if (export_stage2_scores_)
          {
            native_spectrum_id = job.swath_map->sptr->getSpectrumMetaById(
              static_cast<int>(spectrum_index)).id;
          }

          MSSpectrum openms_spectrum;
          OpenSwathDataAccessHelper::convertToOpenMSSpectrum(spectrum_ptr, openms_spectrum);
          openms_spectrum.setMSLevel(2);
          openms_spectrum.sortByPosition();
          if (!openms_spectrum.empty())
          {
            if (native_spectrum_id.empty())
            {
              native_spectrum_id = openms_spectrum.getNativeID().c_str();
            }
            const double total_spectrum_intensity = std::accumulate(
              openms_spectrum.begin(), openms_spectrum.end(), 0.0,
              [](double sum, const Peak1D& peak)
              {
                return sum + peak.getIntensity();
              });
            const double spectrum_mz_span =
              openms_spectrum.size() > 1 ?
              std::max(1e-6, openms_spectrum.back().getMZ() - openms_spectrum.front().getMZ()) :
              1.0;
            const double fragment_tolerance_reference_mz =
              openms_spectrum.size() > 1 ?
              0.5 * (openms_spectrum.front().getMZ() + openms_spectrum.back().getMZ()) :
              openms_spectrum.front().getMZ();
            const double fragment_tolerance_da =
              ms2_params.ppm ?
              Math::ppmToMass<double>(
                halfTolerance_(ms2_params.mz_extraction_window, true, 10.0),
                fragment_tolerance_reference_mz) :
              halfTolerance_(ms2_params.mz_extraction_window, false, 0.05);
            matches.clear();
            fragment_index.querySpectrum(openms_spectrum, job.precursor_queries, matches, false);
            std::unordered_map<Size, FragmentIndex::SpectrumMatch> spectrum_best_matches;
            spectrum_best_matches.reserve(matches.hits_.size());
            for (const auto& match : matches.hits_)
            {
              if (match.peptide_idx_ >= index_peptides.size())
              {
                continue;
              }

              const Size candidate_id = index_peptides[match.peptide_idx_].protein_idx;
              if (candidate_id >= candidates.size())
              {
                continue;
              }

              const auto& candidate = candidates[candidate_id];
              if (static_cast<std::uint16_t>(candidate.precursor_charge) != match.precursor_charge_ ||
                  job.allowed_candidate_ids.find(candidate_id) == job.allowed_candidate_ids.end())
              {
                continue;
              }

              auto best_match_it = spectrum_best_matches.find(candidate_id);
              if (best_match_it == spectrum_best_matches.end() ||
                  betterStage2SpectrumMatch_(match, best_match_it->second))
              {
                spectrum_best_matches[candidate_id] = match;
              }
            }

            std::vector<Stage2ScoredSpectrumCandidate> ranked_candidates;
            ranked_candidates.reserve(spectrum_best_matches.size());
            for (const auto& spectrum_match_item : spectrum_best_matches)
            {
              const auto& match = spectrum_match_item.second;
              const double matched_intensity_fraction =
                total_spectrum_intensity > 0.0 ?
                static_cast<double>(match.matched_intensity_sum_) / total_spectrum_intensity :
                0.0;
              ranked_candidates.push_back({
                spectrum_match_item.first,
                match,
                matched_intensity_fraction,
                computeStage2SpectrumScore_(match, matched_intensity_fraction)
              });
            }

            std::sort(ranked_candidates.begin(), ranked_candidates.end(),
                      [](const Stage2ScoredSpectrumCandidate& lhs,
                         const Stage2ScoredSpectrumCandidate& rhs)
                      {
                        if (lhs.spectrum_score != rhs.spectrum_score) return lhs.spectrum_score > rhs.spectrum_score;
                        return betterStage2SpectrumMatch_(lhs.match, rhs.match);
                      });

            for (const auto& scored_candidate : ranked_candidates)
            {
              const Size candidate_id = scored_candidate.candidate_id;
              const auto& match = scored_candidate.match;
              auto& stats = local_scores[candidate_id];
              stats.best_matched_ions = std::max(stats.best_matched_ions,
                                                 static_cast<Size>(match.num_matched_));
              ++stats.supporting_spectra;
              if (job.run_index < 64)
              {
                stats.supporting_run_mask |= (std::uint64_t{1} << job.run_index);
              }

              const double matched_intensity_fraction = scored_candidate.matched_intensity_fraction;
              const double spectrum_score = scored_candidate.spectrum_score;
              const bool strong_support =
                static_cast<Int>(match.num_matched_) >= stage2_strong_min_matched_ions_ &&
                matched_intensity_fraction >= stage2_strong_min_intensity_fraction_;
              if (strong_support)
              {
                ++stats.strong_supporting_spectra;
                if (job.run_index < 64)
                {
                  stats.strong_supporting_run_mask |= (std::uint64_t{1} << job.run_index);
                }
              }
              upsertStage2RunSummary_(stats, job.run_index, spectrum_score, 1);
              if (export_stage2_scores_)
              {
                updateStage2RunStreak_(
                  run_streak_states[candidate_id], spectrum_index, spectrum_score);
              }
              if (spectrum_score > stats.best_spectrum_score)
              {
                stats.best_spectrum_score = spectrum_score;
                stats.best_spectrum_matched_intensity_fraction = matched_intensity_fraction;
                if (export_stage2_scores_)
                {
                  double poisson_proxy = 0.0;
                  double longest_y_pct = 0.0;
                  Size matched_b_ions = 0;
                  Size matched_y_ions = 0;
                  Size longest_b_run = 0;
                  Size longest_y_run = 0;
                  computeStage2SpectrumDiagnostics_(
                    match,
                    theoretical_b_ions_by_candidate[candidate_id],
                    theoretical_y_ions_by_candidate[candidate_id],
                    openms_spectrum.size(),
                    spectrum_mz_span,
                    fragment_tolerance_da,
                    poisson_proxy,
                    longest_y_pct,
                    matched_b_ions,
                    matched_y_ions,
                    longest_b_run,
                    longest_y_run);
                  stats.best_spectrum_matched_b_ions = matched_b_ions;
                  stats.best_spectrum_matched_y_ions = matched_y_ions;
                  stats.best_spectrum_longest_b_run = longest_b_run;
                  stats.best_spectrum_longest_y_run = longest_y_run;
                  stats.best_spectrum_longest_y_pct = longest_y_pct;
                  stats.best_spectrum_poisson_proxy = poisson_proxy;
                  stats.best_source_file = job.source_file;
                  stats.best_native_spectrum_id = native_spectrum_id;
                }
              }
            }

            if (use_lower_order_null)
            {
              for (Size rank_idx = 0; rank_idx < ranked_candidates.size(); ++rank_idx)
              {
                const auto& scored_candidate = ranked_candidates[rank_idx];
                const auto& candidate = candidates[scored_candidate.candidate_id];
                const Size rank = rank_idx + 1;
                const bool used_for_null =
                  rank >= lower_order_min_rank && rank <= lower_order_max_rank;
                const bool used_for_scoring = rank <= lower_order_scored_ranks;

                if (used_for_null)
                {
                  local_null_scores_by_charge[candidate.precursor_charge].push_back(
                    scored_candidate.spectrum_score);
                }
                if (used_for_scoring || (export_stage2_scores_ && used_for_null))
                {
                  local_observations.push_back({
                    scored_candidate.candidate_id,
                    job.run_index,
                    static_cast<std::uint16_t>(candidate.precursor_charge),
                    rank,
                    static_cast<Size>(scored_candidate.match.num_matched_),
                    scored_candidate.matched_intensity_fraction,
                    scored_candidate.spectrum_score,
                    job.source_file,
                    native_spectrum_id,
                    used_for_scoring,
                    used_for_null
                  });
                }
              }
            }
          }
        }

        const SignedSize done = processed_queries.fetch_add(
          static_cast<SignedSize>(job.query_units_per_spectrum)) +
          static_cast<SignedSize>(job.query_units_per_spectrum);
        SignedSize next_update = next_progress_update.load(std::memory_order_relaxed);
        while (done >= next_update && next_update <= total_progress)
        {
          const SignedSize following_update =
            next_update == total_progress ? total_progress + 1 :
            std::min(total_progress, next_update + progress_step);
          if (next_progress_update.compare_exchange_weak(next_update, following_update))
          {
            std::lock_guard<std::mutex> lock(progress_mutex);
            setProgress(std::min(done, total_progress));
            break;
          }
        }
      }

      for (const auto& streak_item : run_streak_states)
      {
        auto& stats = local_scores[streak_item.first];
        if (streak_item.second.best_streak_score > stats.best_run_streak_score ||
            (streak_item.second.best_streak_score == stats.best_run_streak_score &&
             streak_item.second.best_streak_length > stats.best_run_streak_length))
        {
          stats.best_run_streak_score = streak_item.second.best_streak_score;
          stats.best_run_streak_length = streak_item.second.best_streak_length;
        }

        if (job.run_index < 64)
        {
          if (streak_item.second.best_streak_length >= 2)
          {
            stats.streak_ge_2_run_mask |= (std::uint64_t{1} << job.run_index);
          }
          if (streak_item.second.best_streak_length >= 3)
          {
            stats.streak_ge_3_run_mask |= (std::uint64_t{1} << job.run_index);
          }
        }

        for (Size run_slot = 0; run_slot < stats.top_run_ids.size(); ++run_slot)
        {
          if (stats.top_run_ids[run_slot] == job.run_index)
          {
            stats.top_run_streak_lengths[run_slot] = std::max(
              stats.top_run_streak_lengths[run_slot],
              streak_item.second.best_streak_length);
            break;
          }
        }
      }
    }

    endProgress();

    std::unordered_map<Size, Stage2CandidateStats> merged_candidate_stats;
    merged_candidate_stats.reserve(candidates.size());
    for (const auto& local_scores : local_candidate_stats)
    {
      for (const auto& item : local_scores)
      {
        mergeStage2CandidateStats_(merged_candidate_stats[item.first], item.second);
      }
    }

    std::unordered_map<Size, double> lower_order_combined_pvalues_by_candidate;
    std::unordered_map<Size, double> lower_order_best_local_pvalues_by_candidate;
    std::unordered_map<Size, Size> lower_order_best_local_ranks_by_candidate;
    if (use_lower_order_null)
    {
      std::unordered_map<int, std::vector<double>> null_scores_by_charge;
      std::vector<double> pooled_null_scores;
      std::vector<Stage2LowerOrderObservation> observations;

      for (Size thread_slot = 0; thread_slot < static_cast<Size>(thread_count); ++thread_slot)
      {
        for (auto& item : local_lower_order_null_scores_by_charge[thread_slot])
        {
          pooled_null_scores.insert(
            pooled_null_scores.end(), item.second.begin(), item.second.end());
          auto& destination = null_scores_by_charge[item.first];
          destination.insert(destination.end(), item.second.begin(), item.second.end());
        }
        auto& thread_observations = local_lower_order_observations[thread_slot];
        observations.insert(
          observations.end(), thread_observations.begin(), thread_observations.end());
      }

      std::sort(pooled_null_scores.begin(), pooled_null_scores.end());
      for (auto& item : null_scores_by_charge)
      {
        std::sort(item.second.begin(), item.second.end());
      }

      for (auto& observation : observations)
      {
        const auto charge_null_it = null_scores_by_charge.find(static_cast<int>(observation.precursor_charge));
        const std::vector<double>* active_null = &pooled_null_scores;
        if (charge_null_it != null_scores_by_charge.end() &&
            charge_null_it->second.size() >= lower_order_min_null_scores)
        {
          active_null = &charge_null_it->second;
        }

        observation.local_pvalue = empiricalStage2TailPValue_(*active_null, observation.spectrum_score);
      }

      std::vector<Size> scored_observation_indices;
      scored_observation_indices.reserve(observations.size());
      for (Size observation_idx = 0; observation_idx < observations.size(); ++observation_idx)
      {
        if (observations[observation_idx].used_for_scoring)
        {
          scored_observation_indices.push_back(observation_idx);
        }
      }

      std::sort(scored_observation_indices.begin(), scored_observation_indices.end(),
                [&observations](Size lhs_idx, Size rhs_idx)
                {
                  const auto& lhs = observations[lhs_idx];
                  const auto& rhs = observations[rhs_idx];
                  if (lhs.candidate_id != rhs.candidate_id) return lhs.candidate_id < rhs.candidate_id;
                  if (lhs.run_index != rhs.run_index) return lhs.run_index < rhs.run_index;
                  if (lhs.local_pvalue != rhs.local_pvalue) return lhs.local_pvalue < rhs.local_pvalue;
                  return lhs.rank < rhs.rank;
                });

      for (Size scored_idx = 0; scored_idx < scored_observation_indices.size();)
      {
        const Size candidate_id = observations[scored_observation_indices[scored_idx]].candidate_id;
        std::vector<double> run_pvalues;
        double best_local_pvalue = 1.0;
        Size best_local_rank = 0;

        while (scored_idx < scored_observation_indices.size() &&
               observations[scored_observation_indices[scored_idx]].candidate_id == candidate_id)
        {
          const auto& first_observation = observations[scored_observation_indices[scored_idx]];
          const Size run_index = first_observation.run_index;
          double min_run_pvalue = first_observation.local_pvalue;
          Size min_run_rank = first_observation.rank;
          ++scored_idx;
          while (scored_idx < scored_observation_indices.size())
          {
            const auto& current_observation = observations[scored_observation_indices[scored_idx]];
            if (current_observation.candidate_id != candidate_id ||
                current_observation.run_index != run_index)
            {
              break;
            }

            if (current_observation.local_pvalue < min_run_pvalue ||
                (current_observation.local_pvalue == min_run_pvalue &&
                 current_observation.rank < min_run_rank))
            {
              min_run_pvalue = current_observation.local_pvalue;
              min_run_rank = current_observation.rank;
            }
            ++scored_idx;
          }

          run_pvalues.push_back(min_run_pvalue);
          if (min_run_pvalue < best_local_pvalue ||
              (min_run_pvalue == best_local_pvalue &&
               (best_local_rank == 0 || min_run_rank < best_local_rank)))
          {
            best_local_pvalue = min_run_pvalue;
            best_local_rank = min_run_rank;
          }
        }

        lower_order_combined_pvalues_by_candidate[candidate_id] = combineStage2RunPValues_(run_pvalues);
        lower_order_best_local_pvalues_by_candidate[candidate_id] = best_local_pvalue;
        lower_order_best_local_ranks_by_candidate[candidate_id] = best_local_rank;
      }

      if (export_stage2_scores_)
      {
        bundle.observation_scores.reserve(observations.size());
        for (const auto& observation : observations)
        {
          bundle.observation_scores.push_back({
            candidates[observation.candidate_id].internal_key,
            observation.source_file,
            observation.native_spectrum_id,
            observation.run_index,
            observation.rank,
            observation.used_for_scoring,
            observation.used_for_null,
            observation.matched_ions,
            observation.matched_intensity_fraction,
            observation.spectrum_score,
            observation.local_pvalue
          });
        }
      }
    }

    bundle.best_matched_ions.reserve(merged_candidate_stats.size());
    bundle.peptide_pvalues.reserve(lower_order_combined_pvalues_by_candidate.size());
    if (export_stage2_scores_)
    {
      bundle.candidate_scores.reserve(merged_candidate_stats.size());
    }
    for (const auto& score_item : merged_candidate_stats)
    {
      const auto& candidate = candidates[score_item.first];
      const double composite_score = composeStage2PeptideScore_(score_item.second);
      bundle.best_matched_ions[candidate.internal_key] = score_item.second.best_matched_ions;
      if (use_lower_order_null)
      {
        const auto combined_pvalue_it = lower_order_combined_pvalues_by_candidate.find(score_item.first);
        if (combined_pvalue_it != lower_order_combined_pvalues_by_candidate.end())
        {
          bundle.peptide_pvalues[candidate.internal_key] = combined_pvalue_it->second;
        }
      }
      if (export_stage2_scores_)
      {
        Stage2CandidateScore score_row;
        score_row.peptide_key = candidate.internal_key;
        score_row.peptide_sequence = candidate.peptide_sequence;
        score_row.modified_peptide_sequence = candidate.modified_peptide_sequence;
        score_row.precursor_mz = candidate.precursor_mz;
        score_row.precursor_charge = candidate.precursor_charge;
        score_row.protein_refs = candidate.protein_refs;
        score_row.protein_gene_names_by_accession = candidate.protein_gene_names_by_accession;
        score_row.decoy = candidate.decoy;
        score_row.source_file = score_item.second.best_source_file;
        score_row.native_spectrum_id = score_item.second.best_native_spectrum_id;
        score_row.best_matched_ions = score_item.second.best_matched_ions;
        score_row.supporting_spectra = score_item.second.supporting_spectra;
        score_row.supporting_runs = countStage2SupportingRuns_(score_item.second);
        score_row.strong_supporting_spectra = score_item.second.strong_supporting_spectra;
        score_row.strong_supporting_runs = countStage2StrongSupportingRuns_(score_item.second);
        score_row.runs_with_streak_ge_2 = countStage2RunsWithMinStreak_(score_item.second, 2);
        score_row.runs_with_streak_ge_3 = countStage2RunsWithMinStreak_(score_item.second, 3);
        score_row.best_spectrum_matched_intensity_fraction = score_item.second.best_spectrum_matched_intensity_fraction;
        score_row.best_spectrum_matched_b_ions = score_item.second.best_spectrum_matched_b_ions;
        score_row.best_spectrum_matched_y_ions = score_item.second.best_spectrum_matched_y_ions;
        score_row.best_spectrum_longest_b_run = score_item.second.best_spectrum_longest_b_run;
        score_row.best_spectrum_longest_y_run = score_item.second.best_spectrum_longest_y_run;
        score_row.best_spectrum_longest_y_pct = score_item.second.best_spectrum_longest_y_pct;
        score_row.best_spectrum_poisson_proxy = score_item.second.best_spectrum_poisson_proxy;
        score_row.best_spectrum_score = score_item.second.best_spectrum_score;
        score_row.best_run_streak_length = score_item.second.best_run_streak_length;
        score_row.best_run_streak_score = score_item.second.best_run_streak_score;
        score_row.top_run_score_1 = score_item.second.top_run_scores[0];
        score_row.top_run_score_2 = score_item.second.top_run_scores[1];
        score_row.top_run_score_3 = score_item.second.top_run_scores[2];
        score_row.top_run_streak_length_1 = score_item.second.top_run_streak_lengths[0];
        score_row.top_run_streak_length_2 = score_item.second.top_run_streak_lengths[1];
        score_row.top_run_streak_length_3 = score_item.second.top_run_streak_lengths[2];
        const auto best_local_rank_it = lower_order_best_local_ranks_by_candidate.find(score_item.first);
        if (best_local_rank_it != lower_order_best_local_ranks_by_candidate.end())
        {
          score_row.best_local_rank = best_local_rank_it->second;
        }
        const auto best_local_pvalue_it = lower_order_best_local_pvalues_by_candidate.find(score_item.first);
        if (best_local_pvalue_it != lower_order_best_local_pvalues_by_candidate.end())
        {
          score_row.best_local_pvalue = best_local_pvalue_it->second;
        }
        const auto combined_pvalue_it = lower_order_combined_pvalues_by_candidate.find(score_item.first);
        if (combined_pvalue_it != lower_order_combined_pvalues_by_candidate.end())
        {
          score_row.combined_pvalue = combined_pvalue_it->second;
        }
        score_row.composite_score = composite_score;
        bundle.candidate_scores.emplace(score_row.peptide_key, std::move(score_row));
      }
    }
    return bundle;
  }

  std::unordered_map<std::string, double> FastaEvidenceFilter::computeBenjaminiHochbergQValues(
    const std::unordered_map<std::string, double>& peptide_pvalues)
  {
    std::vector<std::pair<std::string, double>> sorted_pvalues;
    sorted_pvalues.reserve(peptide_pvalues.size());
    for (const auto& item : peptide_pvalues)
    {
      sorted_pvalues.push_back(item);
    }

    std::sort(sorted_pvalues.begin(), sorted_pvalues.end(),
              [](const auto& lhs, const auto& rhs)
              {
                if (lhs.second != rhs.second) return lhs.second < rhs.second;
                return lhs.first < rhs.first;
              });

    std::vector<double> qvalue_buffer(sorted_pvalues.size(), 1.0);
    const double num_tests = static_cast<double>(sorted_pvalues.size());
    for (Size i = 0; i < sorted_pvalues.size(); ++i)
    {
      const double rank = static_cast<double>(i + 1);
      qvalue_buffer[i] = std::clamp(sorted_pvalues[i].second * num_tests / rank, 0.0, 1.0);
    }

    for (SignedSize i = static_cast<SignedSize>(qvalue_buffer.size()) - 2; i >= 0; --i)
    {
      qvalue_buffer[static_cast<Size>(i)] = std::min(
        qvalue_buffer[static_cast<Size>(i)],
        qvalue_buffer[static_cast<Size>(i + 1)]);
    }

    std::unordered_map<std::string, double> qvalues;
    qvalues.reserve(sorted_pvalues.size());
    for (Size i = 0; i < sorted_pvalues.size(); ++i)
    {
      qvalues[sorted_pvalues[i].first] = qvalue_buffer[i];
    }
    return qvalues;
  }

  std::vector<FastaEvidenceFilter::PeptideEntry> FastaEvidenceFilter::selectConfirmedPeptides(
    const std::vector<PeptideEntry>& target_peptides,
    const std::unordered_map<std::string, Size>& best_matched_ions,
    const std::unordered_map<std::string, double>* peptide_qvalues) const
  {
    std::vector<PeptideEntry> confirmed;
    for (const auto& peptide : target_peptides)
    {
      const auto matched_ions_it = best_matched_ions.find(peptide.internal_key);
      if (matched_ions_it == best_matched_ions.end() ||
          matched_ions_it->second < static_cast<Size>(stage2_min_matched_ions_))
      {
        continue;
      }

      bool keep = false;
      if (stage2_mode_ == "raw_score")
      {
        keep = true;
      }
      else
      {
        if (peptide_qvalues != nullptr)
        {
          const auto qvalue_it = peptide_qvalues->find(peptide.internal_key);
          keep = qvalue_it != peptide_qvalues->end() &&
                 qvalue_it->second <= stage2_max_qvalue_;
        }
      }

      if (keep)
      {
        confirmed.push_back(peptide);
      }
    }

    return confirmed;
  }

  std::unordered_set<std::string> FastaEvidenceFilter::selectSupportedProteins(const std::vector<PeptideEntry>& confirmed_peptides,
                                                                               Size min_confirmed_peptides,
                                                                               bool unique_only)
  {
    std::unordered_map<std::string, Size> protein_counts;
    for (const auto& peptide : confirmed_peptides)
    {
      if (unique_only && peptide.protein_refs.size() != 1)
      {
        continue;
      }
      for (const auto& protein_ref : peptide.protein_refs)
      {
        ++protein_counts[protein_ref];
      }
    }

    std::unordered_set<std::string> supported_proteins;
    for (const auto& item : protein_counts)
    {
      if (item.second >= min_confirmed_peptides)
      {
        supported_proteins.insert(item.first);
      }
    }
    return supported_proteins;
  }

  FastaEvidenceFilter::Result FastaEvidenceFilter::filter(const std::vector<RunData>& runs,
                                                          const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                                          const ChromExtractParams& ms1_params,
                                                          const ChromExtractParams& ms2_params,
                                                          int threads) const
  {
    if (runs.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "FastaEvidenceFilter requires at least one DIA run.");
    }
    if (fasta_entries.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "FastaEvidenceFilter requires a non-empty FASTA database.");
    }
    Param stage1_params = param_.copy("Stage1:", true);
    stage1_params.remove("max_concurrent_runs");
    stage1_params.remove("precursor_batch_size");
    stage1_params.setValue("enabled", "false");
    const String stage1_evidence_sources = stage1_params.getValue("evidence_sources").toString();
    const Size batch_size = std::max<Size>(1, stage1_precursor_batch_size_);
    const int thread_count = std::max(1, threads);
    const Size requested_max_concurrent_runs =
      stage1_max_concurrent_runs_ > 0 ?
      stage1_max_concurrent_runs_ :
      static_cast<Size>(thread_count);
    const Size max_concurrent_runs = std::max<Size>(
      1,
      std::min<Size>(
        runs.size(),
        std::min<Size>(requested_max_concurrent_runs, static_cast<Size>(thread_count))));
    const bool stage1_apply_ms2_prefilter = stage1_evidence_sources == "ms2";
    const std::vector<MzCoverageInterval> stage1_coverage =
      stage1_apply_ms2_prefilter ?
      buildSwathMzCoverage_(runs, ms2_params.min_upper_edge_dist) :
      std::vector<MzCoverageInterval>{};
    const bool use_search_space_sharding =
      search_space_max_proteins_per_chunk_ > 0 &&
      fasta_entries.size() > search_space_max_proteins_per_chunk_;

    Size total_target_precursors = 0;
    Size total_stage1_candidate_precursors = 0;
    Size estimated_stage1_target_precursors = 0;
    Size total_stage1_jobs = 0;
    std::vector<PeptideEntry> selected_target_peptides;
    std::vector<PeptideEntry> all_target_peptides;
    std::vector<PeptideEntry> stage1_candidate_peptides_storage;
    const std::vector<PeptideEntry>* stage1_candidate_peptides = nullptr;
    std::shared_ptr<File::TempDir> stage1_shard_dir;
    std::vector<String> shard_paths;

    if (use_search_space_sharding)
    {
      OPENMS_LOG_INFO << "Stage 1: enabling sharded search-space generation with up to "
                      << search_space_max_proteins_per_chunk_ << " proteins per chunk across "
                      << search_space_num_shards_ << " precursor hash shards." << std::endl;
      stage1_shard_dir = std::make_shared<File::TempDir>(search_space_sharding_temp_directory_,
                                                         search_space_sharding_keep_temporary_files_);
      shard_paths.resize(search_space_num_shards_);
      for (Size shard_idx = 0; shard_idx < search_space_num_shards_; ++shard_idx)
      {
        shard_paths[shard_idx] = stage1_shard_dir->getPath() + "/stage1_precursors_shard_" +
                                 String(shard_idx) + ".tsv";
      }

      const Size num_chunks =
        (fasta_entries.size() + search_space_max_proteins_per_chunk_ - 1) /
        search_space_max_proteins_per_chunk_;
      for (Size chunk_idx = 0; chunk_idx < num_chunks; ++chunk_idx)
      {
        const Size chunk_begin = chunk_idx * search_space_max_proteins_per_chunk_;
        const Size chunk_end = std::min(chunk_begin + search_space_max_proteins_per_chunk_,
                                        fasta_entries.size());
        std::vector<FASTAFile::FASTAEntry> fasta_chunk(
          fasta_entries.begin() + static_cast<SignedSize>(chunk_begin),
          fasta_entries.begin() + static_cast<SignedSize>(chunk_end));
        const auto chunk_peptides = generatePeptideEntries_(fasta_chunk, false);
        estimated_stage1_target_precursors += chunk_peptides.size();

        std::vector<std::string> shard_buffers(search_space_num_shards_);
        for (const auto& peptide : chunk_peptides)
        {
          const Size shard_idx =
            std::hash<std::string>{}(peptide.canonical_key) % search_space_num_shards_;
          appendShardEntryBuffer_(shard_buffers[shard_idx], peptide);
        }

        for (Size shard_idx = 0; shard_idx < search_space_num_shards_; ++shard_idx)
        {
          if (shard_buffers[shard_idx].empty())
          {
            continue;
          }

          std::ofstream output(shard_paths[shard_idx].c_str(), std::ios::app);
          output << shard_buffers[shard_idx];
        }
      }

      total_stage1_jobs =
        ((estimated_stage1_target_precursors + batch_size - 1) / batch_size) * runs.size();
      OPENMS_LOG_INFO << "Stage 1: filtering sharded precursor space across "
                      << runs.size() << " DIA runs (estimated up to " << total_stage1_jobs
                      << " TransitionListEvidenceFilter jobs, up to "
                      << max_concurrent_runs << " concurrent runs)." << std::endl;
    }
    else
    {
      all_target_peptides = generatePeptideEntries_(fasta_entries, false);
      if (all_target_peptides.empty())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "No peptide precursors were generated from the FASTA database.");
      }
      total_target_precursors = all_target_peptides.size();
      OPENMS_LOG_INFO << "Original target space: "
                      << total_target_precursors << " peptide precursors, "
                      << fasta_entries.size() << " proteins." << std::endl;

      stage1_candidate_peptides = &all_target_peptides;
      if (stage1_apply_ms2_prefilter)
      {
        stage1_candidate_peptides_storage = filterPeptidesByMzCoverage_(all_target_peptides, stage1_coverage);
        stage1_candidate_peptides = &stage1_candidate_peptides_storage;
        OPENMS_LOG_INFO << "Stage 1: prefiltered " << stage1_candidate_peptides->size()
                        << " of " << total_target_precursors
                        << " target precursors by DIA SWATH m/z coverage before batching."
                        << std::endl;
      }

      total_stage1_candidate_precursors = stage1_candidate_peptides->size();
      const Size num_batches = (stage1_candidate_peptides->size() + batch_size - 1) / batch_size;
      total_stage1_jobs = num_batches * runs.size();
      OPENMS_LOG_INFO << "Stage 1: filtering " << num_batches << " precursor batches across "
                      << runs.size() << " DIA runs (" << total_stage1_jobs
                      << " TransitionListEvidenceFilter jobs, up to "
                      << max_concurrent_runs << " concurrent runs)." << std::endl;
    }

    const auto stage1_begin = std::chrono::steady_clock::now();
    startProgress(0, static_cast<SignedSize>(std::max<Size>(1, total_stage1_jobs)),
                  "Stage 1: filtering precursor batches");
    Size processed_stage1_jobs = 0;
    bool stage1_progress_started = true;

    const auto run_stage1_batches =
      [&](const std::vector<PeptideEntry>& peptide_pool,
          const String& batch_context) -> std::unordered_map<std::string, Size>
      {
        std::unordered_map<std::string, Size> supported_run_counts;
        const Size num_batches = (peptide_pool.size() + batch_size - 1) / batch_size;
        for (Size batch_idx = 0; batch_idx < num_batches; ++batch_idx)
        {
          const Size begin_idx = batch_idx * batch_size;
          const Size end_idx = std::min(begin_idx + batch_size, peptide_pool.size());
          const auto batch_begin = std::chrono::steady_clock::now();
          const OpenSwath::LightTargetedExperiment stage1_experiment =
            buildStage1Experiment_(peptide_pool, begin_idx, end_idx, thread_count);

          for (Size run_offset = 0; run_offset < runs.size(); run_offset += max_concurrent_runs)
          {
            const Size active_runs = std::min(max_concurrent_runs, runs.size() - run_offset);
            const int base_threads_per_run = std::max(1, thread_count / static_cast<int>(active_runs));
            const int extra_threads = thread_count % static_cast<int>(active_runs);
            std::vector<std::future<std::unordered_set<std::string>>> futures;
            futures.reserve(active_runs);

            for (Size local_run_index = 0; local_run_index < active_runs; ++local_run_index)
            {
              const Size run_index = run_offset + local_run_index;
              const int run_threads = base_threads_per_run + (static_cast<int>(local_run_index) < extra_threads ? 1 : 0);
              futures.push_back(std::async(std::launch::async,
                [&, run_index, run_threads]()
                {
                  TransitionListEvidenceFilter stage1_filter;
                  stage1_filter.setParameters(stage1_params);
                  stage1_filter.setLogType(ProgressLogger::NONE);
                  const auto stage1_result = stage1_filter.filter(
                    runs[run_index].swath_maps, stage1_experiment, ms1_params, ms2_params,
                    runs[run_index].pasef, run_threads);

                  std::unordered_set<std::string> run_supported;
                  for (const auto& compound : stage1_result.filtered_targets.getCompounds())
                  {
                    run_supported.insert(compound.id);
                  }
                  return run_supported;
                }));
            }

            for (Size local_run_index = 0; local_run_index < active_runs; ++local_run_index)
            {
              const auto run_supported = futures[local_run_index].get();
              for (const auto& compound_id : run_supported)
              {
                ++supported_run_counts[compound_id];
              }

              ++processed_stage1_jobs;
              const SignedSize capped_progress = std::min<SignedSize>(
                static_cast<SignedSize>(processed_stage1_jobs),
                static_cast<SignedSize>(std::max<Size>(1, total_stage1_jobs)));
              setProgress(capped_progress);
            }
          }

          const double batch_elapsed_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - batch_begin).count();
          const double stage1_elapsed_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - stage1_begin).count();
          OPENMS_LOG_INFO << "Stage 1: finished precursor batch " << (batch_idx + 1)
                          << "/" << num_batches << " (" << (end_idx - begin_idx)
                          << " precursors)" << batch_context << " in "
                          << batch_elapsed_seconds << " s (overall "
                          << stage1_elapsed_seconds << " s)." << std::endl;
        }
        return supported_run_counts;
      };
    try
    {
      if (use_search_space_sharding)
      {
        for (Size shard_idx = 0; shard_idx < shard_paths.size(); ++shard_idx)
        {
          if (!File::exists(shard_paths[shard_idx]))
          {
            continue;
          }

          auto shard_peptides = loadMergedShardEntries_(shard_paths[shard_idx]);
          if (shard_peptides.empty())
          {
            continue;
          }
          total_target_precursors += shard_peptides.size();

          const std::vector<PeptideEntry>* shard_stage1_candidates = &shard_peptides;
          std::vector<PeptideEntry> filtered_shard_peptides;
          if (stage1_apply_ms2_prefilter)
          {
            filtered_shard_peptides = filterPeptidesByMzCoverage_(shard_peptides, stage1_coverage);
            shard_stage1_candidates = &filtered_shard_peptides;
          }
          total_stage1_candidate_precursors += shard_stage1_candidates->size();
          if (shard_stage1_candidates->empty())
          {
            continue;
          }

          const auto shard_supported_run_counts = run_stage1_batches(
            *shard_stage1_candidates,
            " in shard " + String(shard_idx + 1) + "/" + String(shard_paths.size()));
          for (const auto& peptide : *shard_stage1_candidates)
          {
            const auto support_it = shard_supported_run_counts.find(peptide.canonical_key);
            if (support_it == shard_supported_run_counts.end())
            {
              continue;
            }
            if ((aggregation_method_ == "any" && support_it->second > 0) ||
                (aggregation_method_ == "all" && support_it->second == runs.size()))
            {
              selected_target_peptides.push_back(peptide);
            }
          }
        }
      }
      else
      {
        const auto supported_run_counts = run_stage1_batches(*stage1_candidate_peptides, "");
        std::unordered_set<std::string> selected_target_keys;
        for (const auto& item : supported_run_counts)
        {
          if ((aggregation_method_ == "any" && item.second > 0) ||
              (aggregation_method_ == "all" && item.second == runs.size()))
          {
            selected_target_keys.insert(item.first);
          }
        }

        selected_target_peptides.reserve(selected_target_keys.size());
        for (const auto& peptide : all_target_peptides)
        {
          if (selected_target_keys.find(peptide.canonical_key) != selected_target_keys.end())
          {
            selected_target_peptides.push_back(peptide);
          }
        }
      }

      if (use_search_space_sharding)
      {
        OPENMS_LOG_INFO << "Original target space: "
                        << total_target_precursors << " peptide precursors, "
                        << fasta_entries.size() << " proteins." << std::endl;
        if (stage1_apply_ms2_prefilter)
        {
          OPENMS_LOG_INFO << "Stage 1: prefiltered " << total_stage1_candidate_precursors
                          << " of " << total_target_precursors
                          << " target precursors by DIA SWATH m/z coverage before batching."
                          << std::endl;
        }
      }

      if (total_stage1_jobs > 0)
      {
        setProgress(static_cast<SignedSize>(std::max<Size>(1, total_stage1_jobs)));
      }
      endProgress();
      stage1_progress_started = false;
    }
    catch (...)
    {
      if (stage1_progress_started)
      {
        endProgress();
      }
      throw;
    }
    OPENMS_LOG_INFO << "Stage 1: completed in "
                    << std::chrono::duration<double>(std::chrono::steady_clock::now() - stage1_begin).count()
                    << " s." << std::endl;

    if (selected_target_peptides.size() < stage1_min_supported_precursors_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "FastaEvidenceFilter retained only " + String(selected_target_peptides.size()) +
                                       " stage-1 precursors, fewer than Stage1:min_supported_precursors=" +
                                       String(stage1_min_supported_precursors_) + ".");
    }

    std::sort(selected_target_peptides.begin(), selected_target_peptides.end(),
              [](const PeptideEntry& lhs, const PeptideEntry& rhs)
              {
                if (lhs.precursor_mz != rhs.precursor_mz) return lhs.precursor_mz < rhs.precursor_mz;
                if (lhs.precursor_charge != rhs.precursor_charge) return lhs.precursor_charge < rhs.precursor_charge;
                return lhs.modified_peptide_sequence < rhs.modified_peptide_sequence;
              });

    const std::vector<FASTAFile::FASTAEntry> reduced_target_fasta = buildReducedTargetFasta_(fasta_entries, selected_target_peptides);
    if (reduced_target_fasta.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Stage 1 did not retain any proteins in the input FASTA order.");
    }
    OPENMS_LOG_INFO << "Stage 1 retained " << formatRetentionRatio_(selected_target_peptides.size(), total_target_precursors)
                    << " peptide precursors and "
                    << formatRetentionRatio_(reduced_target_fasta.size(), fasta_entries.size())
                    << " proteins." << std::endl;

    std::vector<PeptideEntry> stage2_candidates = selected_target_peptides;
    std::ostringstream stage2_space_message;
    stage2_space_message << "Stage 2 scoring space: "
                         << selected_target_peptides.size() << " target peptide precursors ("
                         << stage2_candidates.size() << " total), "
                         << reduced_target_fasta.size() << " target proteins";
    if (stage2_mode_ == "lower_order_null")
    {
      stage2_space_message << "; decoy-free lower-order null over spectrum ranks "
                           << stage2_lower_order_min_rank_ << "-" << stage2_lower_order_max_rank_
                           << " with top " << stage2_lower_order_scored_ranks_
                           << " scored ranks per spectrum";
    }
    stage2_space_message << ".";
    OPENMS_LOG_INFO << stage2_space_message.str() << std::endl;

    const Stage2ScoreBundle stage2_scores = scoreStage2_(runs, stage2_candidates, ms1_params, ms2_params, threads);
    std::unordered_map<std::string, double> stage2_qvalues;
    if (stage2_mode_ == "lower_order_null")
    {
      stage2_qvalues = computeBenjaminiHochbergQValues(stage2_scores.peptide_pvalues);
    }
    std::vector<PeptideEntry> confirmed_target_peptides = selectConfirmedPeptides(selected_target_peptides,
                                                                                  stage2_scores.best_matched_ions,
                                                                                  stage2_mode_ == "raw_score" ? nullptr : &stage2_qvalues);

    const std::unordered_set<std::string> supported_proteins = selectSupportedProteins(confirmed_target_peptides,
                                                                                        protein_min_confirmed_peptides_,
                                                                                        protein_unique_peptides_only_);

    std::vector<FASTAFile::FASTAEntry> filtered_fasta;
    filtered_fasta.reserve(supported_proteins.size());
    for (const auto& entry : fasta_entries)
    {
      const std::string accession = entry.identifier.c_str();
      if (supported_proteins.find(accession) != supported_proteins.end())
      {
        filtered_fasta.push_back(entry);
      }
    }

    OPENMS_LOG_INFO << "Stage 2 retained "
                    << formatRetentionRatio_(confirmed_target_peptides.size(), selected_target_peptides.size())
                    << " peptide precursors and "
                    << formatRetentionRatio_(filtered_fasta.size(), reduced_target_fasta.size())
                    << " proteins from the Stage 1 retained space." << std::endl;

    if (export_fragments_)
    {
      populateFragments_(confirmed_target_peptides);
    }

    Result result;
    result.filtered_fasta = std::move(filtered_fasta);
    result.confirmed_peptides = confirmed_target_peptides;
    if (export_stage2_scores_)
    {
      if (stage2_mode_ == "lower_order_null" && !stage2_scores.observation_scores.empty())
      {
        result.stage2_candidate_scores.reserve(stage2_scores.observation_scores.size());
        for (const auto& observation : stage2_scores.observation_scores)
        {
          const auto summary_it = stage2_scores.candidate_scores.find(observation.peptide_key);
          if (summary_it == stage2_scores.candidate_scores.end())
          {
            continue;
          }

          Stage2CandidateScore score_row = summary_it->second;
          score_row.source_file = observation.source_file;
          score_row.native_spectrum_id = observation.native_spectrum_id;
          score_row.run_index = observation.run_index;
          score_row.spectrum_rank = observation.spectrum_rank;
          score_row.used_for_scoring = observation.used_for_scoring;
          score_row.used_for_null = observation.used_for_null;
          score_row.observation_matched_ions = observation.observation_matched_ions;
          score_row.observation_matched_intensity_fraction = observation.observation_matched_intensity_fraction;
          score_row.observation_score = observation.observation_score;
          score_row.local_pvalue = observation.local_pvalue;

          const auto matched_ions_it = stage2_scores.best_matched_ions.find(observation.peptide_key);
          const bool passes_matched_ions =
            matched_ions_it != stage2_scores.best_matched_ions.end() &&
            matched_ions_it->second >= static_cast<Size>(stage2_min_matched_ions_);

          const auto qvalue_it = stage2_qvalues.find(observation.peptide_key);
          if (qvalue_it != stage2_qvalues.end())
          {
            score_row.qvalue = qvalue_it->second;
          }
          score_row.accepted = passes_matched_ions &&
                               score_row.qvalue >= 0.0 &&
                               score_row.qvalue <= stage2_max_qvalue_;
          result.stage2_candidate_scores.push_back(std::move(score_row));
        }
      }
      else
      {
        result.stage2_candidate_scores.reserve(stage2_scores.candidate_scores.size());
        for (const auto& item : stage2_scores.candidate_scores)
        {
          Stage2CandidateScore score_row = item.second;
          const auto matched_ions_it = stage2_scores.best_matched_ions.find(item.first);
          const bool passes_matched_ions =
            matched_ions_it != stage2_scores.best_matched_ions.end() &&
            matched_ions_it->second >= static_cast<Size>(stage2_min_matched_ions_);
          if (stage2_mode_ == "lower_order_null")
          {
            const auto qvalue_it = stage2_qvalues.find(item.first);
            if (qvalue_it != stage2_qvalues.end())
            {
              score_row.qvalue = qvalue_it->second;
            }
            score_row.accepted = passes_matched_ions &&
                                 score_row.qvalue >= 0.0 &&
                                 score_row.qvalue <= stage2_max_qvalue_;
          }
          else
          {
            score_row.accepted = passes_matched_ions;
          }
          result.stage2_candidate_scores.push_back(std::move(score_row));
        }
        std::sort(result.stage2_candidate_scores.begin(), result.stage2_candidate_scores.end(),
                  [](const Stage2CandidateScore& lhs, const Stage2CandidateScore& rhs)
                  {
                    if (lhs.composite_score != rhs.composite_score) return lhs.composite_score > rhs.composite_score;
                    return lhs.peptide_key < rhs.peptide_key;
                  });
      }
    }
    result.stage1_supported_precursors = selected_target_peptides.size();
    result.stage2_confirmed_precursors = confirmed_target_peptides.size();
    result.retained_proteins = result.filtered_fasta.size();
    result.summary = "FastaEvidenceFilter retained " +
                     String(formatRetentionRatio_(result.stage2_confirmed_precursors, total_target_precursors).c_str()) +
                     " peptide precursors and " +
                     String(formatRetentionRatio_(result.retained_proteins, fasta_entries.size()).c_str()) +
                     " proteins after stage 2.";
    OPENMS_LOG_INFO << "Final retained space: "
                    << formatRetentionRatio_(result.stage2_confirmed_precursors, total_target_precursors)
                    << " peptide precursors, "
                    << formatRetentionRatio_(result.retained_proteins, fasta_entries.size())
                    << " proteins." << std::endl;
    return result;
  }
}
