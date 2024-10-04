// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------
#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/KERNEL/MobilityPeak1D.h>
#include <OpenMS/ANALYSIS/OPENSWATH/IonMobilityScoring.h>

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/LogStream.h>

// scoring
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SONARScoring.h>

// auxiliary
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerMobilogram.h>

// #define DEBUG_IMSCORING

namespace OpenMS
{
  std::vector<double> IonMobilityScoring::computeGrid_(const std::vector< Mobilogram >& mobilograms, double eps)
  {
    // Extract all ion mobility values across all transitions and produce a
    // grid of all permitted ion mobility values
    std::vector<double> im_grid;
    std::vector< double > mobilityValues;
    for (const auto & im_profile : mobilograms)
    {
      mobilityValues.reserve(mobilityValues.size() + im_profile.size());
      for (const auto & k : im_profile) mobilityValues.push_back(k.getMobility());
    }

    // sort all extracted values
    std::sort(mobilityValues.begin(), mobilityValues.end());

    // Reduce mobility values to grid (consider equal if closer than eps)
    //
    // In some cases there are not enough datapoints available (one of the
    // transitions has no datapoints)
    if (!mobilityValues.empty())
    {
      im_grid.push_back( mobilityValues[0] );
      for (Size k = 1; k < mobilityValues.size(); k++)
      {
        double diff = fabs(mobilityValues[k] - mobilityValues[k-1]);
        if (diff > eps)
        {
          im_grid.push_back( mobilityValues[k] );
        }
      }
    }
    return im_grid;
  }

  void IonMobilityScoring::alignToGrid_(const Mobilogram& profile,
               const std::vector<double>& im_grid,
               Mobilogram & aligned_profile,
               double eps,
               Size & max_peak_idx)
  {
    auto pr_it = profile.begin();
    max_peak_idx = 0;
    double max_int = 0;
    for (Size k = 0; k < im_grid.size(); k++)
    {
      MobilityPeak1D mobi_peak;
      // In each iteration, the IM value of pr_it should be equal to or
      // larger than the master container. If it is equal, we add the current
      // data point, if it is larger we add zero and advance the counter k.
      if (pr_it != profile.end() && fabs(pr_it->getMobility() - im_grid[k] ) < eps*10)
      {
        mobi_peak.setIntensity(pr_it->getIntensity());
        mobi_peak.setMobility(pr_it->getMobility());

        ++pr_it;
      }
      else
      {
        mobi_peak.setIntensity(0.0);
        mobi_peak.setMobility(im_grid[k]);
      }
      // OPENMS_LOG_DEBUG << "grid position " << im_grid[k] << " profile position " << pr_it->first << std::endl;

      // check that we did not advance past
      if (pr_it != profile.end() && (im_grid[k] - pr_it->getMobility()) > eps*10)
      {
        std::cout << " This should never happen, pr_it has advanced past the master container: " << im_grid[k]  << "  / " <<  pr_it->getMobility()  << std::endl;
        throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }

      // collect maxima
      if (pr_it != profile.end() && pr_it->getIntensity() > max_int)
      {
        max_int = pr_it->getIntensity();
        max_peak_idx = k;
      }

      aligned_profile.push_back(mobi_peak);
    }
  }

  std::vector<double> IonMobilityScoring::extractIntensities(const Mobilogram& mobilogram) {
    std::vector<double> mobility_int;
    mobility_int.reserve(mobilogram.size());
    for (const auto& k : mobilogram) {
      mobility_int.push_back(k.getIntensity());
    }
    return mobility_int;
  }

  void IonMobilityScoring::extractIntensities(const std::vector<Mobilogram>& mobilograms,
                                 std::vector<std::vector<double>>& int_values) {
    int_values.clear();
    int_values.reserve(mobilograms.size());
    for (const auto& mobilogram : mobilograms) {
      int_values.emplace_back(extractIntensities(mobilogram));
    }
  }

  // compute ion mobilogram as well as im weighted average. This is based off of integrateWindows() in DIAHelper.cpp
  void IonMobilityScoring::computeIonMobilogram(const SpectrumSequence& spectra,
                              const RangeMZ& mz_range,
                              const RangeMobility& im_range,
                              double & im,
                              double & intensity,
                              Mobilogram & res,
                              double eps)
  {

    // rounding multiplier for the ion mobility value
    // TODO: how to improve this -- will work up to 42949.67296
    double IM_IDX_MULT = 1/eps;

    // We need to store all values that map to the same ion mobility in the
    // same spot in the ion mobilogram (they are not sorted by ion mobility in
    // the input data), therefore create a map to map to bins.
    std::map< int, double> im_chrom;

    for (const auto& spectrum : spectra)
    {
      OPENMS_PRECONDITION(spectrum->getDriftTimeArray() != nullptr, "Cannot filter by drift time if no drift time is available.");
      OPENMS_PRECONDITION(spectrum->getMZArray()->data.size() == spectrum->getIntensityArray()->data.size(), "MZ and Intensity array need to have the same length.");
      OPENMS_PRECONDITION(spectrum->getMZArray()->data.size() == spectrum->getDriftTimeArray()->data.size(), "MZ and Drift Time array need to have the same length.");

      auto mz_arr_end = spectrum->getMZArray()->data.end();
      auto int_it = spectrum->getIntensityArray()->data.begin();
      auto im_it = spectrum->getDriftTimeArray()->data.begin();

      // this assumes that the spectra are sorted!
      auto mz_it = std::lower_bound(spectrum->getMZArray()->data.begin(), mz_arr_end, mz_range.getMin());
      // auto mz_it_end = std::lower_bound(mz_it, mz_arr_end, mz_end);

      // also advance intensity and ion mobility iterator now
      auto iterator_pos = std::distance(spectrum->getMZArray()->data.begin(), mz_it);
      std::advance(int_it, iterator_pos);
      std::advance(im_it, iterator_pos);

      // Start iteration from mz start, end iteration when mz value is larger than mz_end, only store only storing ion mobility values that are in the range
      double mz_end = mz_range.getMax();
      while ((mz_it < mz_arr_end) && (*mz_it < mz_end))
      {
        if (im_range.contains(*im_it))
        {
          intensity += (*int_it);
          im += (*int_it) * (*im_it);
          im_chrom[ int((*im_it)*IM_IDX_MULT) ] += *int_it;
        }
        ++mz_it;
        ++int_it;
        ++im_it;
      }
    }

    // compute the weighted average ion mobility
    if (intensity > 0.)
    {
      im /= intensity;
    }
    else
    {
      im = -1;
      intensity = 0;
    }

    res.reserve(res.size() + im_chrom.size());
    for (const auto& k : im_chrom)
    {
      MobilityPeak1D mobi_peak;
      mobi_peak.setIntensity(k.second);
      mobi_peak.setMobility(k.first / IM_IDX_MULT);
      res.push_back(mobi_peak);
//      res.emplace_back(k.first / IM_IDX_MULT, k.second );
    }
  }

    std::vector<double> sumAlignedIntensities(const std::vector<std::vector<double>>& result) {
        if (result.empty()) {
            return {};
        }

        // Initialize the output vector with the same size as the inner vectors
        std::vector<double> summedIntensities(result[0].size(), 0.0);

        // Sum the intensities across all vectors
        for (const auto& innerVector : result) {
            std::transform(summedIntensities.begin(), summedIntensities.end(),
                           innerVector.begin(), summedIntensities.begin(),
                           std::plus<double>());
        }

        return summedIntensities;
    }

    OpenMS::Mobilogram sumAlignedMobilograms(const std::vector<OpenMS::Mobilogram>& aligned_mobilograms) {
      if (aligned_mobilograms.empty()) {
        return OpenMS::Mobilogram();
      }

      OpenMS::Mobilogram summed_mobilogram;

      // Use the first mobilogram to set the structure
      const auto& first_mobilogram = aligned_mobilograms[0];

      for (size_t j = 0; j < first_mobilogram.size(); ++j) {
        OpenMS::MobilityPeak1D summed_peak;
        summed_peak.setMobility(first_mobilogram[j].getMobility());
        summed_peak.setIntensity(0.0); // Initialize intensity to 0

        // Sum intensities from all mobilograms
        for (const auto& mobilogram : aligned_mobilograms) {
          if (j < mobilogram.size()) {
            summed_peak.setIntensity(summed_peak.getIntensity() + mobilogram[j].getIntensity());
          }
        }

        summed_mobilogram.push_back(summed_peak);
      }

      // Concatenate names of all mobilograms
      std::string concatenated_name;
      for (size_t i = 0; i < aligned_mobilograms.size(); ++i) {
        if (i > 0) {
          concatenated_name += ";";
        }
        concatenated_name += aligned_mobilograms[i].getName();
      }

      // Set the concatenated name for the summed mobilogram
      summed_mobilogram.setName(concatenated_name);

      return summed_mobilogram;
    }



//    std::tuple<size_t, size_t, size_t> findHighestPeak(const std::vector<double>& data, double threshold = 0.1) {
//        if (data.empty()) {
//            return {0, 0, 0};
//        }
//
//        // Find the highest point
//        auto maxIt = std::max_element(data.begin(), data.end());
//        size_t maxIndex = std::distance(data.begin(), maxIt);
//        double maxValue = *maxIt;
//
//        // Find left boundary
//        size_t leftBoundary = maxIndex;
//        for (size_t i = maxIndex; i > 0; --i) {
//            if (data[i] < maxValue * threshold) {
//                leftBoundary = i + 1;
//                break;
//            }
//        }
//
//        // Find right boundary
//        size_t rightBoundary = maxIndex;
//        for (size_t i = maxIndex; i < data.size(); ++i) {
//            if (data[i] < maxValue * threshold) {
//                rightBoundary = i - 1;
//                break;
//            }
//        }
//
//        return {leftBoundary, maxIndex, rightBoundary};
//    }

    std::tuple<size_t, size_t, size_t> findHighestPeak(const PeakPickerMobilogram& picker) {
      const auto& intensities = picker.integrated_intensities_;

      // Check if vectors are empty or of different sizes
      if (intensities.empty() ||
          intensities.size() != picker.left_width_.size() ||
          intensities.size() != picker.right_width_.size()) {
        // Return an "invalid" tuple if there's an issue
        return std::make_tuple(std::numeric_limits<size_t>::max(), 0, 0);
      }

      // Find the iterator pointing to the maximum element
      auto max_it = std::max_element(intensities.begin(), intensities.end());

      // Get the index of the maximum element
      size_t max_index = std::distance(intensities.begin(), max_it);

      // Return the tuple
      return std::make_tuple(max_index,
                             picker.left_width_[max_index],
                             picker.right_width_[max_index]);
    }

    void filterPeakIntensities(OpenMS::Mobilogram& mobilogram,
                               size_t left_index,
                               size_t right_index) {
      // Check if indices are valid (indicating peaks were found)
      bool peaksFound = (left_index != std::numeric_limits<size_t>::max() &&
                         right_index != 0 &&
                         left_index <= right_index);

      // If no peaks were found, return without filtering
      if (!peaksFound) {
        return; // No changes made to the mobilogram
      }

//      std::cout << "Filtering single mobilogram." << std::endl;

      // Create a temporary vector to hold the filtered peaks
      std::vector<OpenMS::MobilityPeak1D> filtered_peaks;

      // Ensure the indices are within bounds
      size_t start = std::max(left_index, static_cast<size_t>(0));
      size_t end = std::min(right_index, mobilogram.size() - 1);
//      std::cout << "Filtering mobilogram. start: " << start << ", end: " << end << std::endl;

      for (size_t i = start; i <= end; ++i) {
        const auto& peak = mobilogram[i];
        filtered_peaks.push_back(peak); // Collect the peaks within the range
      }

      // Clear existing data and replace with filtered peaks
      mobilogram.clear(); // Clear existing data in the Mobilogram
      for (const auto& peak : filtered_peaks) {
        mobilogram.push_back(peak); // Add filtered peaks back to the Mobilogram
      }
    }

    void filterPeakIntensities(std::vector<OpenMS::Mobilogram>& mobilograms,
                               size_t left_index,
                               size_t right_index) {
      // Check if indices are valid (indicating peaks were found)
      bool peaksFound = (left_index != std::numeric_limits<size_t>::max() &&
                         right_index != 0 &&
                         left_index <= right_index);

      // If no peaks were found, return without filtering
      if (!peaksFound) {
        return; // No changes made to mobilograms
      }

//      std::cout << "Number of mobilograms: " << mobilograms.size() << std::endl;

      for (auto& mobilogram : mobilograms) {
        // Create a temporary vector to hold the filtered peaks
        std::vector<OpenMS::MobilityPeak1D> filtered_peaks;

        // Ensure the indices are within bounds
        size_t start = std::max(left_index, static_cast<size_t>(0));
        size_t end = std::min(right_index, mobilogram.size() - 1);
//        std::cout << "Filtering mobilograms. start: " << start << ", end: " << end << std::endl;

        for (size_t i = start; i <= end; ++i) {
          const auto& peak = mobilogram[i];
          filtered_peaks.push_back(peak); // Collect the peaks within the range
        }

        // Clear existing data and replace with filtered peaks
        mobilogram.clear(); // Clear existing data in the Mobilogram
        for (const auto& peak : filtered_peaks) {
          mobilogram.push_back(peak); // Add filtered peaks back to the Mobilogram
        }
      }
    }



    void plotVectorWithPeak(const std::vector<double>& data, size_t left, size_t max, size_t right, int height = 20, int width = 80) {
        if (data.empty()) return;

        double min = *std::min_element(data.begin(), data.end());
        double max_val = *std::max_element(data.begin(), data.end());
        double range = max_val - min;

        std::vector<std::vector<char>> plot(height, std::vector<char>(width, ' '));

        for (size_t i = 0; i < data.size() && i < static_cast<size_t>(width); ++i) {
            int y = static_cast<int>((height - 1) * (data[i] - min) / range);
            char symbol = '*';
            if (i == left || i == right) symbol = '|';
            if (i == max) symbol = '^';
            plot[y][i] = symbol;
        }

        for (int i = height - 1; i >= 0; --i) {
            for (int j = 0; j < width; ++j) {
                std::cout << plot[i][j];
            }
            std::cout << '\n';
        }
    }

//    std::vector<std::vector<double>> filterPeakIntensities(const std::vector<std::vector<double>>& result,
//                                                           size_t leftBoundary, size_t rightBoundary) {
//        std::vector<std::vector<double>> filteredResult;
//        filteredResult.reserve(result.size());
//
//        for (const auto& innerVector : result) {
//            std::vector<double> filteredInnerVector;
//            filteredInnerVector.reserve(rightBoundary - leftBoundary + 1);
//
//            std::copy(innerVector.begin() + leftBoundary,
//                      innerVector.begin() + rightBoundary + 1,
//                      std::back_inserter(filteredInnerVector));
//
//            filteredResult.push_back(std::move(filteredInnerVector));
//        }
//
//        return filteredResult;
//    }
//
//    void filterPeakIntensitiesInPlace(std::vector<std::vector<double>>& result,
//                                      size_t leftBoundary, size_t rightBoundary) {
//        for (auto& innerVector : result) {
//            innerVector.erase(innerVector.begin() + rightBoundary + 1, innerVector.end());
//            innerVector.erase(innerVector.begin(), innerVector.begin() + leftBoundary);
//        }
//    }
//

    /// Constructor
  IonMobilityScoring::IonMobilityScoring() = default;

  /// Destructor
  IonMobilityScoring::~IonMobilityScoring() = default;

  void IonMobilityScoring::driftScoringMS1Contrast(const SpectrumSequence& spectra, const SpectrumSequence& ms1spectrum,
                                                   const std::vector<TransitionType> & transitions,
                                                   OpenSwath_Scores & scores,
                                                   RangeMobility im_range,
                                                   const double dia_extract_window_,
                                                   const bool dia_extraction_ppm_,
                                                   const double drift_extra)
  {
    OPENMS_PRECONDITION(!spectra.empty(), "Spectra cannot be empty")
    OPENMS_PRECONDITION(!ms1spectrum.empty(), "MS1 spectrum cannot be empty")
    OPENMS_PRECONDITION(!transitions.empty(), "Need at least one transition");

    //TODO not sure what error format is best
    for (const auto& s:spectra)
    {
      if (s->getDriftTimeArray() == nullptr)
      {
        OPENMS_LOG_DEBUG << " ERROR: Drift time is missing in ion mobility spectrum!" << std::endl;
        return;
      }
    }

    for (const auto& s:ms1spectrum)
    {
      if (s->getDriftTimeArray() == nullptr)
      {
        OPENMS_LOG_DEBUG << " ERROR: Drift time is missing in MS1 ion mobility spectrum!" << std::endl;
        return;
      }
    }

    double eps = 1e-5; // eps for two grid cells to be considered equal

    // extend IM range by drift_extra
    im_range.scaleBy(drift_extra * 2. + 1); // multiple by 2 because want drift extra to be extended by that amount on either side

    // Step 1: MS2 extraction
    std::vector< OpenMS::Mobilogram > ms2_mobilograms;
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      double im(0), intensity(0);
      Mobilogram res;
      const TransitionType transition = transitions[k];
      // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
      RangeMZ mz_range = DIAHelpers::createMZRangePPM(transition.getProductMZ(), dia_extract_window_, dia_extraction_ppm_);

      computeIonMobilogram(spectra, mz_range, im_range, im, intensity, res, eps);
      ms2_mobilograms.push_back(std::move(res));

    }

    // Step 2: MS1 extraction
    double im(0), intensity(0);
    Mobilogram ms1_profile;
    std::vector< OpenMS::Mobilogram > ms1_mobilograms;
    RangeMZ mz_range = DIAHelpers::createMZRangePPM(transitions[0].getPrecursorMZ(), dia_extract_window_, dia_extraction_ppm_);

    computeIonMobilogram(ms1spectrum, mz_range, im_range, im, intensity, ms1_profile, eps); // TODO: aggregate over isotopes
    ms2_mobilograms.push_back(ms1_profile);

    std::vector<double> im_grid = computeGrid_(ms2_mobilograms, eps); // ensure grid is based on all profiles!
    ms2_mobilograms.pop_back();

    // Step 3: Align the IonMobilogram vectors to the grid
    std::vector< OpenMS::Mobilogram > aligned_ms2_mobilograms;
    for (const auto & mobilogram : ms2_mobilograms)
    {
      Mobilogram aligned_mobilogram;
      Size max_peak_idx = 0;
      alignToGrid_(mobilogram, im_grid, aligned_mobilogram, eps, max_peak_idx);
      aligned_ms2_mobilograms.push_back(aligned_mobilogram);
    }

    Mobilogram aligned_ms1_mobilograms;
    Size max_peak_idx = 0;
    alignToGrid_(ms1_profile, im_grid, aligned_ms1_mobilograms, eps, max_peak_idx);
    std::vector<double> ms1_int_values;
    ms1_int_values.reserve(aligned_ms1_mobilograms.size());
    for (const auto & k : aligned_ms1_mobilograms) ms1_int_values.push_back(k.getIntensity());

    // Step 4: MS1 contrast scores
    std::vector< std::vector< double > > aligned_int_vec;
    extractIntensities(aligned_ms2_mobilograms, aligned_int_vec);
    {
      OpenSwath::MRMScoring mrmscore_;
      mrmscore_.initializeXCorrPrecursorContrastMatrix({ms1_int_values}, aligned_int_vec);
      OPENMS_LOG_DEBUG << "all-all: Contrast Scores : coelution precursor : " << mrmscore_.calcXcorrPrecursorContrastCoelutionScore() << " / shape  precursor " <<
        mrmscore_.calcXcorrPrecursorContrastShapeScore() << std::endl;
      scores.im_ms1_contrast_coelution = mrmscore_.calcXcorrPrecursorContrastCoelutionScore();
      scores.im_ms1_contrast_shape = mrmscore_.calcXcorrPrecursorContrastShapeScore();
    }

    // Step 5: contrast precursor vs summed fragment ions
    std::vector<double> fragment_values;
    fragment_values.resize(ms1_int_values.size(), 0);
    for (Size k = 0; k < fragment_values.size(); k++)
    {
      for (Size i = 0; i < aligned_int_vec.size(); i++)
      {
        fragment_values[k] += aligned_int_vec[i][k];
      }
    }

    OpenSwath::MRMScoring mrmscore_;
    // horribly broken: provides vector of length 1, but expects at least length 2 in calcXcorrPrecursorContrastCoelutionScore()
    mrmscore_.initializeXCorrPrecursorContrastMatrix({ms1_int_values}, {fragment_values});
    OPENMS_LOG_DEBUG << "Contrast Scores : coelution precursor : " << mrmscore_.calcXcorrPrecursorContrastSumFragCoelutionScore() << " / shape  precursor " <<
       mrmscore_.calcXcorrPrecursorContrastSumFragShapeScore() << std::endl;

    // in order to prevent assertion error call calcXcorrPrecursorContrastSumFragCoelutionScore, same as calcXcorrPrecursorContrastCoelutionScore() however different assertion
    scores.im_ms1_sum_contrast_coelution = mrmscore_.calcXcorrPrecursorContrastSumFragCoelutionScore();

    // in order to prevent assertion error call calcXcorrPrecursorContrastSumFragShapeScore(), same as calcXcorrPrecursorContrastShapeScore() however different assertion.
    scores.im_ms1_sum_contrast_shape = mrmscore_.calcXcorrPrecursorContrastSumFragShapeScore();
  }

  void IonMobilityScoring::driftScoringMS1(const SpectrumSequence & spectra,
                                           const std::vector<TransitionType> & transitions,
                                           OpenSwath_Scores & scores,
                                           const double drift_target,
                                           RangeMobility im_range,
                                           const double dia_extract_window_,
                                           const bool dia_extraction_ppm_,
                                           const bool /* use_spline */,
                                           const double drift_extra)
  {
    OPENMS_PRECONDITION(!spectra.empty(), "Spectra cannot be empty")
    OPENMS_PRECONDITION(!transitions.empty(), "Need at least one transition");

    for (auto s:spectra){
      if (s->getDriftTimeArray() == nullptr)
      {
        OPENMS_LOG_DEBUG << " ERROR: Drift time is missing in ion mobility spectrum!" << std::endl;
        return;
      }
    }

    im_range.scaleBy(drift_extra * 2. + 1); // multiple by 2 because want drift extra to be extended by that amount on either side

    double im(0), intensity(0), mz(0);
    RangeMZ mz_range = DIAHelpers::createMZRangePPM(transitions[0].getPrecursorMZ(), dia_extract_window_, dia_extraction_ppm_);

    DIAHelpers::integrateWindow(spectra, mz, im, intensity, mz_range, im_range);

    // Record the measured ion mobility
    scores.im_ms1_drift = im;

    // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
    scores.im_ms1_delta_score = fabs(drift_target - im);
    scores.im_ms1_delta = drift_target - im;
  }

  void IonMobilityScoring::driftScoring(const SpectrumSequence& spectra,
                                        const std::vector<TransitionType> & transitions,
                                        OpenSwath_Scores & scores,
                                        const double drift_target,
                                        RangeMobility im_range,
                                        const double dia_extract_window_,
                                        const bool dia_extraction_ppm_,
                                        const bool /* use_spline */,
                                        const double drift_extra,
                                        const bool apply_im_peak_picking)
  {
    OPENMS_PRECONDITION(!spectra.empty(), "Spectra cannot be empty");
    for (auto s:spectra)
    {
      if (s->getDriftTimeArray() == nullptr)
      {
        OPENMS_LOG_DEBUG << " ERROR: Drift time is missing in ion mobility spectrum!" << std::endl;
        return;
      }
    }

    double eps = 1e-5; // eps for two grid cells to be considered equal

    im_range.scaleBy(drift_extra * 2. + 1); // multiple by 2 because want drift extra to be extended by that amount on either side

    double delta_drift = 0;
    double delta_drift_abs = 0;
    double computed_im = 0;
    double computed_im_weighted = 0;
    double sum_intensity = 0;
    int tr_used = 0;

    // Step 1: MS2 extraction
    std::vector< OpenMS::Mobilogram > ms2_mobilograms;
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      const TransitionType transition = transitions[k];
      Mobilogram res;
      double im(0), intensity(0);

      // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
      RangeMZ mz_range = DIAHelpers::createMZRangePPM(transition.getProductMZ(), dia_extract_window_, dia_extraction_ppm_);

      //double left(transition.getProductMZ()), right(transition.getProductMZ());
      //DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
      computeIonMobilogram(spectra, mz_range, im_range, im, intensity, res, eps);
      res.setName(transition.getNativeID());
      ms2_mobilograms.push_back(std::move(res));

      // TODO what do to about those that have no signal ?
      if (intensity <= 0.0) {continue;} // note: im is -1 then

      tr_used++;

      delta_drift_abs += fabs(drift_target - im);
      delta_drift += drift_target - im;
      OPENMS_LOG_DEBUG << "  -- have delta drift time " << fabs(drift_target -im ) << " with im " << im << std::endl;
      computed_im += im;
      computed_im_weighted += im * intensity;
      sum_intensity += intensity;
      // delta_drift_weighted += delta_drift * normalized_library_intensity[k];
      // weights += normalized_library_intensity[k];
    }

    if (tr_used != 0)
    {
      delta_drift /= tr_used;
      delta_drift_abs /= tr_used;
      computed_im /= tr_used;
      computed_im_weighted /= sum_intensity;
    }
    else
    {
      delta_drift = -1;
      delta_drift_abs = -1;
      computed_im = -1;
      computed_im_weighted = -1;
    }

    OPENMS_LOG_DEBUG << " Scoring delta drift time " << delta_drift << std::endl;
    OPENMS_LOG_DEBUG << " Scoring weighted delta drift time " << computed_im_weighted << " -> get difference " << std::fabs(computed_im_weighted - drift_target)<< std::endl;
    scores.im_delta_score = delta_drift_abs;
    scores.im_delta = delta_drift;
    scores.im_drift = computed_im;
    scores.im_drift_weighted = computed_im_weighted;

    if (delta_drift == -1)
    {
      OPENMS_LOG_DEBUG << "No drift time found for any of the transitions" << std::endl;
    }

    // Step 2: Align the IonMobilogram vectors to the grid
    std::vector<double> im_grid = computeGrid_(ms2_mobilograms, eps);
    std::vector< OpenMS::Mobilogram > aligned_ms2_mobilograms;
    for (const auto & mobilogram : ms2_mobilograms)
    {
      Mobilogram aligned_mobilogram;
      Size max_peak_idx = 0;
      alignToGrid_(mobilogram, im_grid, aligned_mobilogram, eps, max_peak_idx);
        aligned_mobilogram.setName(mobilogram.getName());
      if (!aligned_mobilogram.empty()) aligned_ms2_mobilograms.push_back(std::move(aligned_mobilogram));
    }
    size_t left = 0, max = 0, right = 0;
    if ( apply_im_peak_picking ) {
//        std::vector<double> summedIntensities = sumAlignedIntensities(aligned_ms2_mobilograms);
//        auto [left, max, right] = findHighestPeak(summedIntensities);
////        plotVectorWithPeak(summedIntensities, left, max, right);
//
//        scores.im_drift_left = im_grid[left];
//        scores.im_drift_right = im_grid[right];
//
//        // Filter the original data and overwrite
//        aligned_ms2_mobilograms = filterPeakIntensities(aligned_ms2_mobilograms, left, right);
        if ( !aligned_ms2_mobilograms.empty())
        {
//          std::cout << "In progress..." << std::endl;
          OpenMS::Mobilogram summed_mobilogram = sumAlignedMobilograms(aligned_ms2_mobilograms);
          PeakPickerMobilogram picker_;
          Param picker_params = picker_.getParameters();
          picker_params.setValue("method", "corrected");
          picker_.setParameters(picker_params);
          Mobilogram picked_mobilogram, smoothed_mobilogram;
          picker_.pickMobilogram(summed_mobilogram, picked_mobilogram, smoothed_mobilogram);
          //              size_t max, left, right;
          std::tie(max, left, right) = findHighestPeak(picker_);
          scores.im_drift_left = im_grid[left];
          scores.im_drift_right = im_grid[right];
          filterPeakIntensities(aligned_ms2_mobilograms, left, right);
//          std::cout << "Done picking..." << std::endl;
        }
        else
        {
          scores.im_drift_left = -1;
          scores.im_drift_right = -1;
        }
    }

    // Step 3: Compute cross-correlation scores based on ion mobilograms
    if (aligned_ms2_mobilograms.size() < 2)
    {
      scores.im_xcorr_coelution_score = 0;
      scores.im_xcorr_shape_score = std::numeric_limits<double>::quiet_NaN();
      return;
    }

    std::vector< std::vector< double > > aligned_int_vec;
    extractIntensities(aligned_ms2_mobilograms, aligned_int_vec);
    OpenSwath::MRMScoring mrmscore_;
    mrmscore_.initializeXCorrMatrix(aligned_int_vec);

    double xcorr_coelution_score = mrmscore_.calcXcorrCoelutionScore();
    double xcorr_shape_score = mrmscore_.calcXcorrShapeScore(); // can be nan!

    scores.im_xcorr_coelution_score = xcorr_coelution_score;
    scores.im_xcorr_shape_score = xcorr_shape_score;
  }

  void IonMobilityScoring::driftIdScoring(const SpectrumSequence& spectra,
                                          const std::vector<TransitionType> & transition,
                                          MRMTransitionGroupType& trgr_detect,
                                          OpenSwath_Scores &scores,
                                          const double drift_target,
                                          RangeMobility im_range,
                                          const double dia_extract_window_,
                                          const bool dia_extraction_ppm_,
                                          const bool /* use_spline */,
                                          const double drift_extra,
                                          const bool apply_im_peak_picking)
  {
      // OPENMS_PRECONDITION(spectrum != nullptr, "Spectrum cannot be null");
      // OPENMS_PRECONDITION(!transition.empty(), "Need at least one transition");

      // if (spectrum->getDriftTimeArray() == nullptr)
      // {
      //   OPENMS_LOG_DEBUG << " ERROR: Drift time is missing in ion mobility spectrum!" << std::endl;
      //   return;
      // }

      OPENMS_PRECONDITION(!spectra.empty(), "Spectra cannot be empty");
      for (auto s:spectra)
      {
        if (s->getDriftTimeArray() == nullptr)
        {
          OPENMS_LOG_DEBUG << " ERROR: Drift time is missing in ion mobility spectrum!" << std::endl;
          return;
        }
      }

      double eps = 1e-5; // eps for two grid cells to be considered equal

      im_range.scaleBy(drift_extra * 2. + 1); // multiple by 2 because want drift extra to be extended by that amount on either side

      Mobilogram res;
      double im(0), intensity(0);
      // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
      // double left(transition[0].getProductMZ()), right(transition[0].getProductMZ());
      // DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
      // DIAHelpers::integrateDriftSpectrum(spectrum, left, right, im, intensity, drift_lower_used, drift_upper_used);

      RangeMZ mz_range = DIAHelpers::createMZRangePPM(transition[0].getProductMZ(), dia_extract_window_, dia_extraction_ppm_);
      computeIonMobilogram(spectra, mz_range, im_range, im, intensity, res, eps);

      // Record the measured ion mobility
      scores.im_drift = im;

      // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
      scores.im_delta_score = fabs(drift_target - im);
      scores.im_delta = drift_target - im;

      scores.im_log_intensity = std::log(intensity+1);
      OPENMS_LOG_DEBUG << "Identification Transition IM Scoring for " << transition[0].transition_name << " range (" << im_range.getMin() << " - " << im_range.getMax() << ") IM = " << im << " im_delta = " << drift_target - im << " int = " << intensity << " log int = " << std::log(intensity+1) << std::endl;


      //////////////////////////////////////////////////////////////////////////////////////////
      // Cross-Correlation of Identification against Detection Mobilogram Features

      // double eps = 1e-5; // eps for two grid cells to be considered equal

      // IonMobilogram: a data structure that holds points <im_value, intensity>
      std::vector< OpenMS::Mobilogram > mobilograms;

      // Step 1: MS2 detection transitions extraction
      for (std::size_t k = 0; k < trgr_detect.getTransitions().size(); k++)
      {
        double detection_im(0), detection_intensity(0);
        Mobilogram detection_mobilograms;
        const TransitionType detection_transition = trgr_detect.getTransitions()[k];
        // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
        // double detection_left(detection_transition.getProductMZ()), detection_right(detection_transition.getProductMZ());
        // DIAHelpers::adjustExtractionWindow(detection_right, detection_left, dia_extract_window_, dia_extraction_ppm_);

        // integrateDriftSpectrum(spectrum, detection_left, detection_right, detection_im, detection_intensity, detection_mobilograms, eps, drift_lower_used, drift_upper_used);

        RangeMZ detection_mz_range = DIAHelpers::createMZRangePPM(detection_transition.getProductMZ(), dia_extract_window_, dia_extraction_ppm_);
        computeIonMobilogram(spectra, detection_mz_range, im_range, detection_im, detection_intensity, detection_mobilograms, eps);
        detection_mobilograms.setName(detection_transition.getNativeID());
        mobilograms.push_back( std::move(detection_mobilograms) );
      }

      // Step 2: MS2 single identification transition extraction
      double identification_im(0), identification_intensity(0);
      Mobilogram identification_mobilogram;
      // double identification_left(transition[0].getProductMZ()), identification_right(transition[0].getProductMZ());
      // DIAHelpers::adjustExtractionWindow(identification_right, identification_left, dia_extract_window_, dia_extraction_ppm_);
      // integrateDriftSpectrum(spectrum, identification_left, identification_right, identification_im, identification_intensity, identification_mobilogram, eps, drift_lower_used, drift_upper_used); // TODO: aggregate over isotopes

      RangeMZ identification_mz_range = DIAHelpers::createMZRangePPM(transition[0].getProductMZ(), dia_extract_window_, dia_extraction_ppm_);
      computeIonMobilogram(spectra, identification_mz_range, im_range, identification_im, identification_intensity, identification_mobilogram, eps);
      identification_mobilogram.setName(transition[0].getNativeID());
      mobilograms.push_back(identification_mobilogram);

      // Check to make sure IM of identification is not -1, otherwise assign 0 for scores
      if ( identification_im!=-1 )
      {

        std::vector<double> im_grid = computeGrid_(mobilograms, eps); // ensure grid is based on all profiles!
        mobilograms.pop_back();

        // Step 3: Align the IonMobilogram vectors to the grid
        std::vector< OpenMS::Mobilogram > aligned_mobilograms;
        for (const auto &mobilogram : mobilograms)
        {
          Mobilogram aligned_mobilogram;
          Size max_peak_idx = 0;
          alignToGrid_(mobilogram, im_grid, aligned_mobilogram, eps, max_peak_idx);
          aligned_mobilogram.setName(mobilogram.getName());
          aligned_mobilograms.push_back(std::move(aligned_mobilogram));
        }

        size_t left = 0, max = 0, right = 0;
        if ( apply_im_peak_picking ) {
//            std::vector<double> summedIntensities = sumAlignedIntensities(aligned_mobilograms);
//            std::tie(left, max, right) = findHighestPeak(summedIntensities);
////            plotVectorWithPeak(summedIntensities, left, max, right);
//
//            scores.im_drift_left = im_grid[left];
//            scores.im_drift_right = im_grid[right];
//
//            // Filter the original data and overwrite
//            aligned_mobilograms = filterPeakIntensities(aligned_mobilograms, left, right);
              if ( !aligned_mobilograms.empty())
              {
                OpenMS::Mobilogram summed_mobilogram = sumAlignedMobilograms(aligned_mobilograms);
  //              std::cout << "In progress..." << std::endl;
  //              std::vector< OpenMS::Mobilogram > picked_mobilograms;
  //              std::vector< OpenMS::Mobilogram > smoothed_mobilograms;
  //              for (const auto &mobilogram : aligned_mobilograms)
  //              {
  //                PeakPickerMobilogram picker_;
  //                Param picker_params = picker_.getParameters();
  //                picker_params.setValue("method", "corrected");
  //                picker_.setParameters(picker_params);
  //                Mobilogram picked_mobilogram, smoothed_mobilogram;
  //                picker_.pickMobilogram(mobilogram, picked_mobilogram, smoothed_mobilogram);
  //                //              picked_mobilograms.sortByIntensity();
  //                picked_mobilograms.push_back(std::move(picked_mobilogram));
  //                smoothed_mobilograms.push_back(std::move(smoothed_mobilogram));
  //              }
  //
                PeakPickerMobilogram picker_;
                Param picker_params = picker_.getParameters();
                picker_params.setValue("method", "corrected");
                picker_.setParameters(picker_params);
                Mobilogram picked_mobilogram, smoothed_mobilogram;
                picker_.pickMobilogram(summed_mobilogram, picked_mobilogram, smoothed_mobilogram);
  //              size_t max, left, right;
                std::tie(max, left, right) = findHighestPeak(picker_);
                scores.im_drift_left = im_grid[left];
                scores.im_drift_right = im_grid[right];
  //              std::cout << "Highest peak found:" << std::endl;
  //              std::cout << "Index: " << max << std::endl;
  //              std::cout << "Left boundary: " << left << std::endl;
  //              std::cout << "Right boundary: " << right << std::endl;
  //              std::cout << "Intensity: " << picker_.integrated_intensities_[max] << std::endl;
                // Filter the mobilograms
                filterPeakIntensities(aligned_mobilograms, left, right);
  //              std::cout << "Done picking..." << std::endl;
              }
              else
              {
                scores.im_drift_left = -1;
                scores.im_drift_right = -1;
              }
        }

        Mobilogram aligned_identification_mobilogram;
        Size max_peak_idx = 0;
        alignToGrid_(identification_mobilogram,
                    im_grid,
                    aligned_identification_mobilogram,
                    eps,
                    max_peak_idx);
        aligned_identification_mobilogram.setName(transition[0].getNativeID());
        if ( apply_im_peak_picking )
        {
////            auto [left, max, right] = findHighestPeak(identification_int_values);
//            // based filtering on left and right width form detecting ion mobilograms
//            // TODO: Would inidividual boundaries help?
//            // TODO: pass identficaition int as nested double vector.
//            // Filter the original data and overwrite
//            std::vector <std::vector<double>> identification_int_values_filtered = filterPeakIntensities({identification_int_values}, left, right);
//            identification_int_values = identification_int_values_filtered[0];
//            std::cout << "In progress..." << std::endl;

//            PeakPickerMobilogram picker_;
//            Mobilogram picked_id_mobilogram, smoothed_id_mobilogram;
//            picker_.pickMobilogram(aligned_identification_mobilogram, picked_id_mobilogram, smoothed_id_mobilogram);
//
//            size_t max, left, right;
//            std::tie(max, left, right) = findHighestPeak(picker_);
//            scores.im_drift_left = im_grid[left];
//            scores.im_drift_right = im_grid[right];
//            std::cout << "Highest peak found:" << std::endl;
//            std::cout << "Index: " << max << std::endl;
//            std::cout << "Left boundary: " << left << std::endl;
//            std::cout << "Right boundary: " << right << std::endl;
//            std::cout << "Intensity: " << picker_.integrated_intensities_[max] << std::endl;
            // Filter the mobilograms
            if ( !aligned_identification_mobilogram.empty())
            {
              filterPeakIntensities(aligned_identification_mobilogram, left, right);
            }

//            std::cout << "Done picking..." << std::endl;
        }

        std::vector< std::vector< double > > aligned_int_vec;
        extractIntensities(aligned_mobilograms, aligned_int_vec);
        std::vector< double > identification_int_values = extractIntensities(aligned_identification_mobilogram);
        // Step 4: MS1 contrast scores
        {
          OpenSwath::MRMScoring mrmscore_;
          mrmscore_.initializeXCorrPrecursorContrastMatrix({identification_int_values}, aligned_int_vec);
          OPENMS_LOG_DEBUG << "all-all: Contrast Scores : coelution identification transition : "
                           << mrmscore_.calcXcorrPrecursorContrastCoelutionScore()
                           << " / shape  identification transition " <<
                           mrmscore_.calcXcorrPrecursorContrastShapeScore() << std::endl;
          scores.im_ms1_contrast_coelution = mrmscore_.calcXcorrPrecursorContrastCoelutionScore();
          scores.im_ms1_contrast_shape = mrmscore_.calcXcorrPrecursorContrastShapeScore();
        }

        // Step 5: contrast precursor vs summed fragment ions
        std::vector<double> fragment_values;
        fragment_values.resize(identification_int_values.size(), 0);
        for (Size k = 0; k < fragment_values.size(); k++)
        {
          for (Size i = 0; i < aligned_int_vec.size(); i++)
          {
            fragment_values[k] += aligned_int_vec[i][k];
          }
        }

        OpenSwath::MRMScoring mrmscore_;
        // horribly broken: provides vector of length 1, but expects at least length 2 in calcXcorrPrecursorContrastCoelutionScore()
        mrmscore_.initializeXCorrPrecursorContrastMatrix({identification_int_values}, {fragment_values});
        OPENMS_LOG_DEBUG << "Contrast Scores : coelution identification transition : "
                         << mrmscore_.calcXcorrPrecursorContrastSumFragCoelutionScore()
                         << " / shape  identification transition " <<
                         mrmscore_.calcXcorrPrecursorContrastSumFragShapeScore() << std::endl;

        // in order to prevent assertion error call calcXcorrPrecursorContrastSumFragCoelutionScore, same as calcXcorrPrecursorContrastCoelutionScore() however different assertion
        scores.im_ms1_sum_contrast_coelution = mrmscore_.calcXcorrPrecursorContrastSumFragCoelutionScore();

        // in order to prevent assertion error call calcXcorrPrecursorContrastSumFragShapeScore(), same as calcXcorrPrecursorContrastShapeScore() however different assertion.
        scores.im_ms1_sum_contrast_shape = mrmscore_.calcXcorrPrecursorContrastSumFragShapeScore();
      } else {
        OPENMS_LOG_DEBUG << "Identification Transition IM Scoring for " << transition[0].transition_name << " was -1. There was most likely no drift spectrum for the transition, setting cross-correlation scores to 0!" << std::endl;
        scores.im_ms1_contrast_coelution = 0;
        scores.im_ms1_contrast_shape = 0;
        scores.im_ms1_sum_contrast_coelution = 0;
        scores.im_ms1_sum_contrast_shape = 0;
      }
  }

}
