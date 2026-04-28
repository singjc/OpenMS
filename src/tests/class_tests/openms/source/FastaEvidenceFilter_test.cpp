// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/FastaEvidenceFilter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace OpenMS;
using namespace OpenSwath;
using namespace std;

namespace
{
  ChromExtractParams makeExtractParams(double mz_window, bool ppm = false, double im_window = -1.0)
  {
    ChromExtractParams params;
    params.min_upper_edge_dist = 0.0;
    params.mz_extraction_window = mz_window;
    params.im_extraction_window = im_window;
    params.ppm = ppm;
    params.extraction_function = "tophat";
    params.rt_extraction_window = -1.0;
    params.extra_rt_extract = 0.0;
    return params;
  }

  MSSpectrum makeSpectrum(double rt, const vector<pair<double, double>>& peaks)
  {
    MSSpectrum spectrum;
    spectrum.setRT(rt);
    for (const auto& peak_data : peaks)
    {
      Peak1D peak;
      peak.setMZ(peak_data.first);
      peak.setIntensity(peak_data.second);
      spectrum.push_back(peak);
    }
    spectrum.sortByPosition();
    return spectrum;
  }

  SwathMap makeSwathMap(bool ms1,
                        double lower,
                        double upper,
                        const vector<MSSpectrum>& spectra)
  {
    PeakMap* peak_map = new PeakMap;
    for (const auto& spectrum : spectra)
    {
      peak_map->addSpectrum(spectrum);
    }
    shared_ptr<PeakMap> exp(peak_map);

    SwathMap map(lower, upper, (lower + upper) * 0.5, ms1);
    map.sptr = SpectrumAccessPtr(new SpectrumAccessOpenMS(exp));
    return map;
  }

  FASTAFile::FASTAEntry makeFastaEntry(const String& identifier, const String& sequence, const String& description = "")
  {
    FASTAFile::FASTAEntry entry;
    entry.identifier = identifier;
    entry.description = description;
    entry.sequence = sequence;
    return entry;
  }

  FastaEvidenceFilter makeFilterForSingleChargePeptides()
  {
    FastaEvidenceFilter filter;
    Param params = filter.getParameters();
    params.setValue("SearchSpace:missed_cleavages", 0);
    params.setValue("SearchSpace:precursor:min_charge", 2);
    params.setValue("SearchSpace:precursor:max_charge", 2);
    params.setValue("SearchSpace:fragment:min_charge", 1);
    params.setValue("SearchSpace:fragment:max_charge", 1);
    params.setValue("SearchSpace:fragment:max_per_precursor", 4);
    filter.setParameters(params);
    return filter;
  }
}

START_TEST(FastaEvidenceFilter, "$Id$")

START_SECTION(FastaEvidenceFilter())
{
  FastaEvidenceFilter filter;
  TEST_EQUAL(filter.getParameters().getValue("aggregation_method").toString(), "any")
  TEST_EQUAL(filter.getParameters().getValue("Stage2:mode").toString(), "qvalue")
  TEST_EQUAL(static_cast<Int>(filter.getParameters().getValue("Protein:min_confirmed_peptides")), 1)
}
END_SECTION

START_SECTION((std::string makeCanonicalPeptideKey(const std::string&, int)))
{
  TEST_EQUAL(FastaEvidenceFilter::makeCanonicalPeptideKey("PEPTIDE", 2), "PEPTIDE/2")
  TEST_EQUAL(FastaEvidenceFilter::makeCanonicalPeptideKey("AAC(Carbamidomethyl)DM(Oxidation)K", 3),
             "AAC(Carbamidomethyl)DM(Oxidation)K/3")
}
END_SECTION

START_SECTION((std::vector<PeptideEntry> generatePeptideEntries(const std::vector<FASTAFile::FASTAEntry>&) const))
{
  FastaEvidenceFilter filter = makeFilterForSingleChargePeptides();
  Param params = filter.getParameters();
  params.setValue("SearchSpace:min_size", 6);
  params.setValue("SearchSpace:max_size", 6);
  filter.setParameters(params);

  const vector<FASTAFile::FASTAEntry> fasta_entries{
    makeFastaEntry("protA", "AACDMK", "Protein A OS=Homo sapiens GN=GENEA PE=1 SV=1")
  };
  const auto peptides = filter.generatePeptideEntries(fasta_entries);

  TEST_EQUAL(peptides.size(), 2)

  set<string> modified_sequences;
  for (const auto& peptide : peptides)
  {
    modified_sequences.insert(peptide.modified_peptide_sequence);
    TEST_EQUAL(peptide.peptide_sequence, "AACDMK")
    TEST_EQUAL(peptide.precursor_charge, 2)
    TEST_EQUAL(peptide.protein_refs.size(), 1)
    TEST_EQUAL(peptide.protein_refs[0], "protA")
    TEST_EQUAL(peptide.protein_gene_names_by_accession.size(), 1)
    TEST_EQUAL(peptide.protein_gene_names_by_accession.at("protA"), "GENEA")
    TEST_EQUAL(peptide.fragments.size(), 4)
    TEST_TRUE(!peptide.canonical_key.empty())
    for (const auto& fragment : peptide.fragments)
    {
      TEST_TRUE(fragment.product_mz > 0.0)
      TEST_TRUE(fragment.product_charge == 1)
      TEST_TRUE(!fragment.product_type.empty())
      TEST_TRUE(fragment.product_ordinal > 2)
    }
  }

  TEST_TRUE(modified_sequences.count("AAC(Carbamidomethyl)DMK") == 1)
  TEST_TRUE(modified_sequences.count("AAC(Carbamidomethyl)DM(Oxidation)K") == 1)
}
END_SECTION

START_SECTION((std::unordered_set<std::string> selectSupportedProteins(const std::vector<PeptideEntry>&, Size, bool)))
{
  FastaEvidenceFilter::PeptideEntry shared;
  shared.protein_refs = {"protA", "protB"};
  FastaEvidenceFilter::PeptideEntry unique;
  unique.protein_refs = {"protA"};
  const vector<FastaEvidenceFilter::PeptideEntry> confirmed_peptides{shared, unique};

  const auto supported_all = FastaEvidenceFilter::selectSupportedProteins(confirmed_peptides, 1, false);
  TEST_EQUAL(supported_all.size(), 2)
  TEST_TRUE(supported_all.count("protA") == 1)
  TEST_TRUE(supported_all.count("protB") == 1)

  const auto supported_unique = FastaEvidenceFilter::selectSupportedProteins(confirmed_peptides, 1, true);
  TEST_EQUAL(supported_unique.size(), 1)
  TEST_TRUE(supported_unique.count("protA") == 1)

  const auto supported_min2 = FastaEvidenceFilter::selectSupportedProteins(confirmed_peptides, 2, false);
  TEST_EQUAL(supported_min2.size(), 1)
  TEST_TRUE(supported_min2.count("protA") == 1)
}
END_SECTION

START_SECTION((selectConfirmedPeptides() - raw_score and qvalue modes))
{
  FastaEvidenceFilter::PeptideEntry target_peptide;
  target_peptide.internal_key = "PEPTIDE/2";
  target_peptide.canonical_key = "PEPTIDE/2";
  target_peptide.modified_peptide_sequence = "PEPTIDE";
  target_peptide.precursor_charge = 2;

  const vector<FastaEvidenceFilter::PeptideEntry> target_peptides{target_peptide};
  const unordered_map<string, double> best_scores{{"PEPTIDE/2", 6.0}};

  FastaEvidenceFilter raw_filter;
  Param raw_params = raw_filter.getParameters();
  raw_params.setValue("Stage2:mode", "raw_score");
  raw_params.setValue("Stage2:min_matched_ions", 3);
  raw_filter.setParameters(raw_params);

  const auto raw_confirmed = raw_filter.selectConfirmedPeptides(target_peptides, best_scores, {});
  TEST_EQUAL(raw_confirmed.size(), 1)

  FastaEvidenceFilter qvalue_keep_filter;
  Param qvalue_keep_params = qvalue_keep_filter.getParameters();
  qvalue_keep_params.setValue("Stage2:mode", "qvalue");
  qvalue_keep_params.setValue("Stage2:max_qvalue", 0.05);
  qvalue_keep_params.setValue("Stage2:min_matched_ions", 3);
  qvalue_keep_filter.setParameters(qvalue_keep_params);

  const vector<FastaEvidenceFilter::PeptideScoreRecord> keep_records{
    {"PEPTIDE/2", 6.0, false},
    {"DECOY_PEPTIDE/2", 4.0, true}
  };
  const auto qvalue_confirmed = qvalue_keep_filter.selectConfirmedPeptides(target_peptides, best_scores, keep_records);
  TEST_EQUAL(qvalue_confirmed.size(), 1)

  FastaEvidenceFilter qvalue_drop_filter;
  Param qvalue_drop_params = qvalue_drop_filter.getParameters();
  qvalue_drop_params.setValue("Stage2:mode", "qvalue");
  qvalue_drop_params.setValue("Stage2:max_qvalue", 0.05);
  qvalue_drop_params.setValue("Stage2:min_matched_ions", 3);
  qvalue_drop_filter.setParameters(qvalue_drop_params);

  const vector<FastaEvidenceFilter::PeptideScoreRecord> drop_records{
    {"DECOY_PEPTIDE/2", 7.0, true},
    {"PEPTIDE/2", 6.0, false}
  };
  const auto qvalue_dropped = qvalue_drop_filter.selectConfirmedPeptides(target_peptides, best_scores, drop_records);
  TEST_EQUAL(qvalue_dropped.size(), 0)
}
END_SECTION

START_SECTION((filter() - stage1 support can be removed by stage2))
{
  FastaEvidenceFilter filter = makeFilterForSingleChargePeptides();
  Param params = filter.getParameters();
  params.setValue("SearchSpace:min_size", 7);
  params.setValue("SearchSpace:max_size", 7);
  params.setValue("Stage1:evidence_sources", "ms1");
  params.setValue("Stage1:min_supported_precursors", 1);
  params.setValue("Stage2:mode", "raw_score");
  params.setValue("Stage2:decoys", "false");
  params.setValue("Stage2:min_matched_ions", 1);
  filter.setParameters(params);

  const vector<FASTAFile::FASTAEntry> fasta_entries{makeFastaEntry("protA", "AAAAAAK")};
  const auto peptides = filter.generatePeptideEntries(fasta_entries);
  TEST_EQUAL(peptides.size(), 1)

  const double precursor_mz = peptides[0].precursor_mz;
  vector<OpenSwath::SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(true, 0.0, 0.0, {makeSpectrum(10.0, {{precursor_mz, 1000.0}})}));
  swath_maps.push_back(makeSwathMap(false, precursor_mz - 10.0, precursor_mz + 10.0,
                                    {makeSpectrum(12.0, {{50.0, 1000.0}, {60.0, 900.0}, {70.0, 800.0}})}));

  FastaEvidenceFilter::RunData run;
  run.swath_maps = std::move(swath_maps);
  run.pasef = false;

  const auto result = filter.filter({run}, fasta_entries, makeExtractParams(0.01), makeExtractParams(0.01), 1);
  TEST_EQUAL(result.stage1_supported_precursors, 1)
  TEST_EQUAL(result.stage2_confirmed_precursors, 0)
  TEST_EQUAL(result.retained_proteins, 0)
  TEST_EQUAL(result.confirmed_peptides.size(), 0)
  TEST_EQUAL(result.filtered_fasta.size(), 0)
  TEST_EQUAL(result.summary.hasSubstring("0 of 1"), true)
}
END_SECTION

START_SECTION((filter() - stage2 can confirm a supported precursor))
{
  FastaEvidenceFilter filter = makeFilterForSingleChargePeptides();
  Param params = filter.getParameters();
  params.setValue("SearchSpace:min_size", 7);
  params.setValue("SearchSpace:max_size", 7);
  params.setValue("Stage1:evidence_sources", "ms1");
  params.setValue("Stage1:min_supported_precursors", 1);
  params.setValue("Stage2:mode", "raw_score");
  params.setValue("Stage2:decoys", "false");
  params.setValue("Stage2:min_matched_ions", 1);
  filter.setParameters(params);

  const vector<FASTAFile::FASTAEntry> fasta_entries{makeFastaEntry("protA", "AAAAAAK")};
  const auto peptides = filter.generatePeptideEntries(fasta_entries);
  TEST_EQUAL(peptides.size(), 1)
  TEST_TRUE(!peptides[0].fragments.empty())

  vector<pair<double, double>> stage2_peaks;
  for (const auto& fragment : peptides[0].fragments)
  {
    stage2_peaks.emplace_back(fragment.product_mz, 1000.0);
  }

  const double precursor_mz = peptides[0].precursor_mz;
  vector<OpenSwath::SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(true, 0.0, 0.0, {makeSpectrum(10.0, {{precursor_mz, 1000.0}})}));
  swath_maps.push_back(makeSwathMap(false, precursor_mz - 10.0, precursor_mz + 10.0,
                                    {makeSpectrum(12.0, stage2_peaks)}));

  FastaEvidenceFilter::RunData run;
  run.swath_maps = std::move(swath_maps);
  run.pasef = false;

  const auto result = filter.filter({run}, fasta_entries, makeExtractParams(0.01), makeExtractParams(0.01), 1);
  TEST_EQUAL(result.stage1_supported_precursors, 1)
  TEST_EQUAL(result.stage2_confirmed_precursors, 1)
  TEST_EQUAL(result.retained_proteins, 1)
  TEST_EQUAL(result.confirmed_peptides.size(), 1)
  TEST_EQUAL(result.filtered_fasta.size(), 1)
  TEST_EQUAL(result.filtered_fasta[0].identifier, "protA")
  TEST_EQUAL(result.summary.hasSubstring("1 of 1"), true)
}
END_SECTION

START_SECTION((filter() - stage2 multithreaded scoring matches single-threaded scoring))
{
  FastaEvidenceFilter filter = makeFilterForSingleChargePeptides();
  Param params = filter.getParameters();
  params.setValue("SearchSpace:min_size", 7);
  params.setValue("SearchSpace:max_size", 7);
  params.setValue("Stage1:evidence_sources", "ms1");
  params.setValue("Stage1:min_supported_precursors", 1);
  params.setValue("Stage2:mode", "raw_score");
  params.setValue("Stage2:decoys", "false");
  params.setValue("Stage2:min_matched_ions", 1);
  filter.setParameters(params);

  const vector<FASTAFile::FASTAEntry> fasta_entries{
    makeFastaEntry("protA", "AAAAAAK"),
    makeFastaEntry("protB", "CCCCCCK")
  };
  const auto peptides = filter.generatePeptideEntries(fasta_entries);
  TEST_EQUAL(peptides.size(), 2)

  vector<pair<double, double>> ms1_peaks;
  vector<pair<double, double>> stage2_peaks_a;
  vector<pair<double, double>> stage2_peaks_b;
  for (const auto& peptide : peptides)
  {
    ms1_peaks.emplace_back(peptide.precursor_mz, 1000.0);
    vector<pair<double, double>>& stage2_peaks = peptide.protein_refs[0] == "protA" ? stage2_peaks_a : stage2_peaks_b;
    for (const auto& fragment : peptide.fragments)
    {
      stage2_peaks.emplace_back(fragment.product_mz, 1000.0);
    }
  }

  double precursor_mz_a = 0.0;
  double precursor_mz_b = 0.0;
  for (const auto& peptide : peptides)
  {
    if (peptide.protein_refs[0] == "protA")
    {
      precursor_mz_a = peptide.precursor_mz;
    }
    else
    {
      precursor_mz_b = peptide.precursor_mz;
    }
  }

  vector<OpenSwath::SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(true, 0.0, 0.0, {makeSpectrum(10.0, ms1_peaks)}));
  swath_maps.push_back(makeSwathMap(false, precursor_mz_a - 10.0, precursor_mz_a + 10.0,
                                    {makeSpectrum(12.0, stage2_peaks_a), makeSpectrum(12.5, stage2_peaks_a)}));
  swath_maps.push_back(makeSwathMap(false, precursor_mz_b - 10.0, precursor_mz_b + 10.0,
                                    {makeSpectrum(13.0, stage2_peaks_b), makeSpectrum(13.5, stage2_peaks_b)}));
  swath_maps.push_back(makeSwathMap(false, precursor_mz_a - 10.0, precursor_mz_a + 10.0,
                                    {makeSpectrum(14.0, stage2_peaks_a), makeSpectrum(14.5, stage2_peaks_a)}));
  swath_maps.push_back(makeSwathMap(false, precursor_mz_b - 10.0, precursor_mz_b + 10.0,
                                    {makeSpectrum(15.0, stage2_peaks_b), makeSpectrum(15.5, stage2_peaks_b)}));

  FastaEvidenceFilter::RunData run;
  run.swath_maps = std::move(swath_maps);
  run.pasef = false;

  const auto single_thread = filter.filter({run}, fasta_entries, makeExtractParams(0.01), makeExtractParams(0.01), 1);
  const auto multi_thread = filter.filter({run}, fasta_entries, makeExtractParams(0.01), makeExtractParams(0.01), 4);

  TEST_EQUAL(single_thread.stage1_supported_precursors, 2)
  TEST_EQUAL(single_thread.stage2_confirmed_precursors, 2)
  TEST_EQUAL(single_thread.retained_proteins, 2)
  TEST_EQUAL(single_thread.stage1_supported_precursors, multi_thread.stage1_supported_precursors)
  TEST_EQUAL(single_thread.stage2_confirmed_precursors, multi_thread.stage2_confirmed_precursors)
  TEST_EQUAL(single_thread.retained_proteins, multi_thread.retained_proteins)

  set<string> single_keys;
  set<string> multi_keys;
  for (const auto& peptide : single_thread.confirmed_peptides)
  {
    single_keys.insert(peptide.internal_key);
  }
  for (const auto& peptide : multi_thread.confirmed_peptides)
  {
    multi_keys.insert(peptide.internal_key);
  }
  TEST_EQUAL(single_keys.size(), multi_keys.size())
  TEST_TRUE(single_keys == multi_keys)

  set<string> single_proteins;
  set<string> multi_proteins;
  for (const auto& entry : single_thread.filtered_fasta)
  {
    single_proteins.insert(entry.identifier.c_str());
  }
  for (const auto& entry : multi_thread.filtered_fasta)
  {
    multi_proteins.insert(entry.identifier.c_str());
  }
  TEST_EQUAL(single_proteins.size(), multi_proteins.size())
  TEST_TRUE(single_proteins == multi_proteins)
}
END_SECTION

START_SECTION((filter() - stage1 multi-run aggregation preserves support with concurrent runs))
{
  FastaEvidenceFilter filter = makeFilterForSingleChargePeptides();
  Param params = filter.getParameters();
  params.setValue("aggregation_method", "all");
  params.setValue("SearchSpace:min_size", 7);
  params.setValue("SearchSpace:max_size", 7);
  params.setValue("Stage1:evidence_sources", "ms1");
  params.setValue("Stage1:min_supported_precursors", 1);
  params.setValue("Stage2:mode", "raw_score");
  params.setValue("Stage2:decoys", "false");
  params.setValue("Stage2:min_matched_ions", 1);
  filter.setParameters(params);

  const vector<FASTAFile::FASTAEntry> fasta_entries{makeFastaEntry("protA", "AAAAAAK")};
  const auto peptides = filter.generatePeptideEntries(fasta_entries);
  TEST_EQUAL(peptides.size(), 1)
  TEST_TRUE(!peptides[0].fragments.empty())

  vector<pair<double, double>> stage2_peaks;
  for (const auto& fragment : peptides[0].fragments)
  {
    stage2_peaks.emplace_back(fragment.product_mz, 1000.0);
  }

  const double precursor_mz = peptides[0].precursor_mz;

  vector<OpenSwath::SwathMap> swath_maps_a;
  swath_maps_a.push_back(makeSwathMap(true, 0.0, 0.0, {makeSpectrum(10.0, {{precursor_mz, 1000.0}})}));
  swath_maps_a.push_back(makeSwathMap(false, precursor_mz - 10.0, precursor_mz + 10.0,
                                      {makeSpectrum(12.0, stage2_peaks)}));

  vector<OpenSwath::SwathMap> swath_maps_b;
  swath_maps_b.push_back(makeSwathMap(true, 0.0, 0.0, {makeSpectrum(20.0, {{precursor_mz, 900.0}})}));
  swath_maps_b.push_back(makeSwathMap(false, precursor_mz - 10.0, precursor_mz + 10.0,
                                      {makeSpectrum(22.0, stage2_peaks)}));

  FastaEvidenceFilter::RunData run_a;
  run_a.swath_maps = std::move(swath_maps_a);
  run_a.pasef = false;

  FastaEvidenceFilter::RunData run_b;
  run_b.swath_maps = std::move(swath_maps_b);
  run_b.pasef = false;

  const auto result = filter.filter({run_a, run_b}, fasta_entries, makeExtractParams(0.01), makeExtractParams(0.01), 4);
  TEST_EQUAL(result.stage1_supported_precursors, 1)
  TEST_EQUAL(result.stage2_confirmed_precursors, 1)
  TEST_EQUAL(result.retained_proteins, 1)
  TEST_EQUAL(result.filtered_fasta.size(), 1)
  TEST_EQUAL(result.filtered_fasta[0].identifier, "protA")
}
END_SECTION

START_SECTION((filter() - stage1 precursor batching preserves support across chunks))
{
  FastaEvidenceFilter filter = makeFilterForSingleChargePeptides();
  Param params = filter.getParameters();
  params.setValue("SearchSpace:min_size", 7);
  params.setValue("SearchSpace:max_size", 7);
  params.setValue("Stage1:evidence_sources", "ms1");
  params.setValue("Stage1:min_supported_precursors", 1);
  params.setValue("Stage1:precursor_batch_size", 1);
  params.setValue("Stage2:mode", "raw_score");
  params.setValue("Stage2:decoys", "false");
  params.setValue("Stage2:min_matched_ions", 1);
  filter.setParameters(params);

  const vector<FASTAFile::FASTAEntry> fasta_entries{
    makeFastaEntry("protA", "AAAAAAK"),
    makeFastaEntry("protB", "CCCCCCK")
  };
  const auto peptides = filter.generatePeptideEntries(fasta_entries);
  TEST_EQUAL(peptides.size(), 2)

  vector<pair<double, double>> ms1_peaks;
  vector<pair<double, double>> stage2_peaks;
  double min_precursor_mz = numeric_limits<double>::max();
  double max_precursor_mz = numeric_limits<double>::lowest();
  for (const auto& peptide : peptides)
  {
    ms1_peaks.emplace_back(peptide.precursor_mz, 1000.0);
    min_precursor_mz = min(min_precursor_mz, peptide.precursor_mz);
    max_precursor_mz = max(max_precursor_mz, peptide.precursor_mz);
    for (const auto& fragment : peptide.fragments)
    {
      stage2_peaks.emplace_back(fragment.product_mz, 1000.0);
    }
  }

  vector<OpenSwath::SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(true, 0.0, 0.0, {makeSpectrum(10.0, ms1_peaks)}));
  swath_maps.push_back(makeSwathMap(false, min_precursor_mz - 10.0, max_precursor_mz + 10.0,
                                    {makeSpectrum(12.0, stage2_peaks)}));

  FastaEvidenceFilter::RunData run;
  run.swath_maps = std::move(swath_maps);
  run.pasef = false;

  const auto result = filter.filter({run}, fasta_entries, makeExtractParams(0.01), makeExtractParams(0.01), 1);
  TEST_EQUAL(result.stage1_supported_precursors, 2)
  TEST_EQUAL(result.stage2_confirmed_precursors, 2)
  TEST_EQUAL(result.retained_proteins, 2)
}
END_SECTION

START_SECTION((filter() - stage1 ms2 batching skips precursors outside SWATH coverage))
{
  FastaEvidenceFilter filter = makeFilterForSingleChargePeptides();
  Param params = filter.getParameters();
  params.setValue("SearchSpace:min_size", 7);
  params.setValue("SearchSpace:max_size", 7);
  params.setValue("Stage1:evidence_sources", "ms2");
  params.setValue("Stage1:min_supported_precursors", 1);
  params.setValue("Stage1:precursor_batch_size", 1);
  params.setValue("Stage2:mode", "raw_score");
  params.setValue("Stage2:decoys", "false");
  params.setValue("Stage2:min_matched_ions", 1);
  filter.setParameters(params);

  const vector<FASTAFile::FASTAEntry> fasta_entries{
    makeFastaEntry("protA", "AAAAAAK"),
    makeFastaEntry("protB", "CCCCCCK")
  };
  const auto peptides = filter.generatePeptideEntries(fasta_entries);
  TEST_EQUAL(peptides.size(), 2)

  double precursor_mz_in_range = 0.0;
  vector<pair<double, double>> in_range_fragments;
  for (const auto& peptide : peptides)
  {
    if (peptide.protein_refs[0] == "protA")
    {
      precursor_mz_in_range = peptide.precursor_mz;
      for (const auto& fragment : peptide.fragments)
      {
        in_range_fragments.emplace_back(fragment.product_mz, 1000.0);
      }
    }
  }
  TEST_TRUE(precursor_mz_in_range > 0.0)
  TEST_TRUE(!in_range_fragments.empty())

  vector<OpenSwath::SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(true, 0.0, 0.0, {makeSpectrum(10.0, {{50.0, 10.0}})}));
  swath_maps.push_back(makeSwathMap(false, precursor_mz_in_range - 5.0, precursor_mz_in_range + 5.0,
                                    {makeSpectrum(12.0, in_range_fragments)}));

  FastaEvidenceFilter::RunData run;
  run.swath_maps = std::move(swath_maps);
  run.pasef = false;

  const auto result = filter.filter({run}, fasta_entries, makeExtractParams(0.01), makeExtractParams(0.01), 4);
  TEST_EQUAL(result.stage1_supported_precursors, 1)
  TEST_EQUAL(result.stage2_confirmed_precursors, 1)
  TEST_EQUAL(result.retained_proteins, 1)
  TEST_EQUAL(result.filtered_fasta.size(), 1)
  TEST_EQUAL(result.filtered_fasta[0].identifier, "protA")
}
END_SECTION

END_TEST
