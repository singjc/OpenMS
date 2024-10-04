// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerMobilogram.h>
#include <iostream>
#include <iomanip>  // Add this for std::setw
#include <vector>
#include <algorithm>  // For std::min_element and std::max_element
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

namespace OpenMS
{
    PeakPickerMobilogram::PeakPickerMobilogram() :
            DefaultParamHandler("PeakPickerMobilogram")
//            ProgressLogger()
    {
        defaults_.setValue("sgolay_frame_length", 9, "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.");
        defaults_.setValue("sgolay_polynomial_order", 3, "Order of the polynomial that is fitted.");
        defaults_.setValue("gauss_width", 0.002, "Gaussian width in seconds, estimated peak size.");
        defaults_.setValue("use_gauss", "true", "Use Gaussian filter for smoothing (alternative is Savitzky-Golay filter)");
        defaults_.setValidStrings("use_gauss", {"false","true"});

        defaults_.setValue("peak_width", -1.0, "Force a certain minimal peak_width on the data (e.g. extend the peak at least by this amount on both sides) in seconds. -1 turns this feature off.");
        defaults_.setValue("signal_to_noise", 1.0, "Signal-to-noise threshold at which a peak will not be extended any more. Note that setting this too high (e.g. 1.0) can lead to peaks whose flanks are not fully captured.");
        defaults_.setMinFloat("signal_to_noise", 0.0);

        defaults_.setValue("sn_win_len", 1, "Signal to noise window length.");
        defaults_.setValue("sn_bin_count", 4, "Signal to noise bin count.");
        defaults_.setValue("write_sn_log_messages", "false", "Write out log messages of the signal-to-noise estimator in case of sparse windows or median in rightmost histogram bin");
        defaults_.setValidStrings("write_sn_log_messages", {"true","false"});

        defaults_.setValue("remove_overlapping_peaks", "false", "Try to remove overlapping peaks during peak picking");
        defaults_.setValidStrings("remove_overlapping_peaks", {"false","true"});

        defaults_.setValue("method", "corrected", "Which method to choose for mobilogram peak-picking (OpenSWATH legacy on raw data, corrected picking on smoothed mobilogram).");
        defaults_.setValidStrings("method", {"legacy","corrected","crawdad"});

        // write defaults into Param object param_
        defaultsToParam_();
        updateMembers_();

        // PeakPickerHiRes pp_;
        Param pepi_param = pp_.getDefaults();
        pepi_param.setValue("signal_to_noise", signal_to_noise_);
        // disable spacing constraints, since we're dealing with chromatograms
        pepi_param.setValue("spacing_difference", 0.0);
        pepi_param.setValue("spacing_difference_gap", 0.0);
        pepi_param.setValue("report_FWHM", "true");
        pepi_param.setValue("report_FWHM_unit", "absolute");
        pp_.setParameters(pepi_param);
    }

#include <iostream>
#include <iomanip>  // For std::setw
#include <vector>
#include <algorithm>  // For std::min_element and std::max_element

    void plotMobilogram(const OpenMS::Mobilogram& mobilogram, int height = 10, int width = 50) {
      if (mobilogram.empty()) {
        std::cout << "Mobilogram is empty.\n";
        return;
      }

      std::vector<double> intensities;
      std::vector<double> positions;

      for (const auto& peak : mobilogram) {
        intensities.push_back(peak.getIntensity());
        positions.push_back(peak.getPosition()[0]);
      }

      if (intensities.empty() || positions.empty()) {
        std::cout << "No valid data in mobilogram.\n";
        return;
      }

      double min_intensity = *std::min_element(intensities.begin(), intensities.end());
      double max_intensity = *std::max_element(intensities.begin(), intensities.end());
      double min_position = *std::min_element(positions.begin(), positions.end());
      double max_position = *std::max_element(positions.begin(), positions.end());

      if (min_intensity == max_intensity) {
        std::cout << "All intensities are the same: " << min_intensity << "\n";
        return;
      }

      if (min_position == max_position) {
        std::cout << "All positions are the same: " << min_position << "\n";
        return;
      }

      std::vector<std::vector<char>> plot(height, std::vector<char>(width, ' '));
      std::vector<int> index_map(width, -1);

      for (size_t i = 0; i < intensities.size(); ++i) {
        double normalized_intensity = (intensities[i] - min_intensity) / (max_intensity - min_intensity);
        double normalized_position = (positions[i] - min_position) / (max_position - min_position);

        int y = static_cast<int>(std::floor(normalized_intensity * (height - 1)));
        int x = static_cast<int>(std::floor(normalized_position * (width - 1)));

        if (y >= 0 && y < height && x >= 0 && x < width) {
          plot[height - 1 - y][x] = '*';
          index_map[x] = i;
        }
      }

      // Plot mobilogram with fixed-width formatting for alignment
      std::cout << "Mobilogram (max intensity: " << max_intensity << ")\n";
      for (const auto& row : plot) {
        for (char c : row) {
          std::cout << std::setw(2) << c;  // Ensure consistent width for each character
        }
        std::cout << '\n';
      }

      // Print indices below the intensity plot with a separator
      for (int idx : index_map) {
        if (idx != -1) {
          std::cout << std::setw(2) << idx << "-";  // Use "-" as separator between numbers
        } else {
          std::cout << "-";  // Ensure spacing and separator when no index is available
        }
      }
      std::cout << '\n';

      // Add a line to separate the plot and the index
      std::cout << std::string(width * 2, '-') << '\n';  // Adjust width of separator line
      std::cout << "Position range: " << min_position << " - " << max_position << "\n";
    }

    void writeMobilogramToCSV(const OpenMS::Mobilogram& mobilogram, const std::string& filename) {
      std::ofstream outFile(filename, std::ios::app);  // Open file in append mode

      if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
      }

      // Write header if the file is empty
      outFile.seekp(0, std::ios::end);
      if (outFile.tellp() == 0) {
        outFile << "Name,Mobility,Intensity" << std::endl;
      }

      // Write data
      for (size_t i = 0; i < mobilogram.size(); ++i) {
        const auto& peak = mobilogram[i];
        outFile << mobilogram.getName() << ","
                << peak.getMobility() << ","
                << peak.getIntensity() << std::endl;
      }

      outFile.close();

//      std::cout << "Mobilogram data appended to " << filename << std::endl;
    }

    void writeIntegratedDataToCSV(const OpenMS::Mobilogram& mobilogram,
                                  const std::vector<double>& integrated_intensities,
                                  const std::vector<Size>& left_width,
                                  const std::vector<Size>& right_width,
                                  const std::string& filename)
    {
      std::ofstream outFile(filename, std::ios::app);  // Open file in append mode

      if (!outFile.is_open())
      {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
      }

      // Write header if the file is empty
      outFile.seekp(0, std::ios::end);
      if (outFile.tellp() == 0) {
        outFile << "Name,Index,Integrated_Intensity,Left_Width,Right_Width,Left_Mobility,Right_Mobility" << std::endl;
      }

      // Write data
      for (Size i = 0; i < integrated_intensities.size(); i++)
      {
        outFile << mobilogram.getName() << ","
                << i << ","
                << std::fixed << std::setprecision(4) << integrated_intensities[i] << ","
                << left_width[i] << ","
                << right_width[i] << ","
                << std::fixed << std::setprecision(6) << mobilogram[left_width[i]].getMobility() << ","
                << std::fixed << std::setprecision(6) << mobilogram[right_width[i]].getMobility()
                << std::endl;
      }

      outFile.close();

//      std::cout << "Integrated data appended to " << filename << std::endl;
    }

    std::string generateRandomID(size_t length) {
      const char charset[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
      std::random_device rd;  // Obtain a random number from hardware
      std::mt19937 generator(rd());  // Seed the generator
      std::uniform_int_distribution<> distribution(0, sizeof(charset) - 2); // Exclude the null terminator

      std::string randomID;
      for (size_t i = 0; i < length; ++i) {
        randomID += charset[distribution(generator)];
      }
      return randomID;
    }


    void PeakPickerMobilogram::pickMobilogram(const Mobilogram& mobilogram, Mobilogram& picked_mobilogram)
    {
      Mobilogram s;
      pickMobilogram(mobilogram, picked_mobilogram, s);
    }

    void PeakPickerMobilogram::pickMobilogram(Mobilogram mobilogram, Mobilogram& picked_mobilogram, Mobilogram& smoothed_mobilogram)
    {
      if (!mobilogram.isSorted())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Mobilogram needs to be sorted by position.");
      }

      if (mobilogram.empty())
      {
          OPENMS_LOG_DEBUG << " ====  Mobilogram empty. Skip picking.";
          return;
      }
      else
      {
        OPENMS_LOG_DEBUG << " ====  Picking mobilogram for " << mobilogram.getName() << " at RT " << mobilogram.getRT() <<
        " with " << mobilogram.size() << " peaks (start at IM " << mobilogram[0].getMobility() << " to IM " << mobilogram.back().getMobility() << ") "
        "using method \'" << method_ << "\'" << std::endl;
      }
      picked_mobilogram.clear();

      // Generate a random ID and concatenate it with the names
      std::string concatenated_name;
      concatenated_name += mobilogram.getName();
      // Append a random unique ID
      std::string randomID = generateRandomID(8); // Generate an 8-character ID
      concatenated_name += ";" + randomID;

      // Set the concatenated name for the summed mobilogram
      mobilogram.setName(concatenated_name);

      // Smooth the mobilogram
      smoothed_mobilogram = mobilogram;
//      Mobilogram smoothed_mobilogram_gauss;
//      smoothed_mobilogram_gauss = mobilogram;
      if (!use_gauss_)
      {
          sgolay_.filter(smoothed_mobilogram);
//          gauss_.filter(smoothed_mobilogram);
      }
      else
      {
          gauss_.filter(smoothed_mobilogram);
      }

//      plotMobilogram(mobilogram);
//      writeMobilogramToCSV(mobilogram, "mobilogram_raw_data.csv");
//      plotMobilogram(smoothed_mobilogram);
//      writeMobilogramToCSV(smoothed_mobilogram_gauss, "mobilogram_sgolay_smooth_data.csv");
//      writeMobilogramToCSV(smoothed_mobilogram, "mobilogram_gauss_smooth_data.csv");
      // Find initial seeds (peak picking)
      pp_.pick(smoothed_mobilogram, picked_mobilogram);
//      std::cout << "Picked " << picked_mobilogram.size() << " mobilogram peaks." << std::endl;

      if (method_ == "legacy")
      {
        // Legacy is to use the original chromatogram for peak-detection
        pickMobilogram_(mobilogram, picked_mobilogram);
        if (remove_overlapping_)
          removeOverlappingPeaks_(mobilogram, picked_mobilogram);

        // for peak integration, we want to use the raw data
        integratePeaks_(mobilogram);
      }
      else if (method_ == "corrected")
      {
        // use the smoothed chromatogram to derive the peak boundaries
        pickMobilogram_(smoothed_mobilogram, picked_mobilogram);
        if (remove_overlapping_)
          removeOverlappingPeaks_(smoothed_mobilogram, picked_mobilogram);

        // for peak integration, we want to use the raw data
        integratePeaks_(mobilogram);
      }

//      // Store the result in the picked_chromatogram
//      OPENMS_POSTCONDITION(picked_mobilogram.getFloatDataArrays().size() == 1 &&
//                             picked_mobilogram.getFloatDataArrays()[IDX_FWHM].getName() == "FWHM", "Swath: PeakPicking did not deliver FWHM attributes.")

      picked_mobilogram.getFloatDataArrays().resize(SIZE_OF_FLOATINDICES);
      picked_mobilogram.getFloatDataArrays()[IDX_ABUNDANCE].setName("IntegratedIntensity");
      picked_mobilogram.getFloatDataArrays()[IDX_LEFTBORDER].setName("leftWidth");
      picked_mobilogram.getFloatDataArrays()[IDX_RIGHTBORDER].setName("rightWidth");
      // just copy FWHM from initial peak picking
//      std::cout << "Size of picked_mobilogram: " << picked_mobilogram.size() << std::endl;
      picked_mobilogram.getFloatDataArrays()[IDX_ABUNDANCE].reserve(picked_mobilogram.size());
      picked_mobilogram.getFloatDataArrays()[IDX_LEFTBORDER].reserve(picked_mobilogram.size());
      picked_mobilogram.getFloatDataArrays()[IDX_RIGHTBORDER].reserve(picked_mobilogram.size());
      for (Size i = 0; i < picked_mobilogram.size(); i++)
      {
//        std::cout << "Size of picked_mobilogram: " << picked_mobilogram.size() << " integrated_intensitie: " << integrated_intensities_[i] << " left_width: " << left_width_[i] << " right_width: " << right_width_[i] << std::endl;
        picked_mobilogram.getFloatDataArrays()[IDX_ABUNDANCE].push_back(integrated_intensities_[i]);
        picked_mobilogram.getFloatDataArrays()[IDX_LEFTBORDER].push_back((float)mobilogram[left_width_[i]].getMobility());
        picked_mobilogram.getFloatDataArrays()[IDX_RIGHTBORDER].push_back((float)mobilogram[right_width_[i]].getMobility());
      }
//      writeIntegratedDataToCSV(mobilogram, integrated_intensities_, left_width_, right_width_, "peak_picking_data.csv");

    }

    void PeakPickerMobilogram::pickMobilogram_(const Mobilogram& mobilogram, Mobilogram& picked_mobilogram)
    {

        integrated_intensities_.clear();
        left_width_.clear();
        right_width_.clear();
        integrated_intensities_.reserve(picked_mobilogram.size());
        left_width_.reserve(picked_mobilogram.size());
        right_width_.reserve(picked_mobilogram.size());

//        if (signal_to_noise_ > 0.0)
//        {
//          snt_.init(mobilogram);
//        }
        Size current_peak = 0;
        for (Size i = 0; i < picked_mobilogram.size(); i++)
        {
          const double central_peak_rt = picked_mobilogram[i].getMobility();
          current_peak = findClosestPeak_(mobilogram, central_peak_rt, current_peak);
          const Size min_i = current_peak;

          // peak core found, now extend it to the left
          Size k = 2;
          while ((min_i - k + 1) > 0
                 //&& std::fabs(chromatogram[min_i-k].getMZ() - peak_raw_data.begin()->first) < spacing_difference*min_spacing
                 && (mobilogram[min_i - k].getIntensity() < mobilogram[min_i - k + 1].getIntensity()
                     || (peak_width_ > 0.0 && std::fabs(mobilogram[min_i - k].getMobility() - central_peak_rt) < peak_width_)))
//                 && (signal_to_noise_ <= 0.0 || snt_.getSignalToNoise(min_i - k) >= signal_to_noise_))
          {
            ++k;
          }
          int left_idx = min_i - k + 1;

          // to the right
          k = 2;
          while ((min_i + k) < mobilogram.size()
                 //&& std::fabs(chromatogram[min_i+k].getMZ() - peak_raw_data.rbegin()->first) < spacing_difference*min_spacing
                 && (mobilogram[min_i + k].getIntensity() < mobilogram[min_i + k - 1].getIntensity()
                     || (peak_width_ > 0.0 && std::fabs(mobilogram[min_i + k].getMobility() - central_peak_rt) < peak_width_)))
//                 && (signal_to_noise_ <= 0.0 || snt_.getSignalToNoise(min_i + k) >= signal_to_noise_) )
          {
            ++k;
          }
          int right_idx = min_i + k - 1;

          left_width_.push_back(left_idx);
          right_width_.push_back(right_idx);
          integrated_intensities_.push_back(0);

          OPENMS_LOG_DEBUG << "Found peak at " << central_peak_rt << " with intensity "  << picked_mobilogram[i].getIntensity()
                           << " and borders " << mobilogram[left_width_[i]].getMobility() << " " << mobilogram[right_width_[i]].getMobility() <<
            " (" << mobilogram[right_width_[i]].getMobility() - mobilogram[left_width_[i]].getMobility() << ") "
                           << 0 << " weighted IM " << /* weighted_mz << */ std::endl;
        }

    }

    void PeakPickerMobilogram::integratePeaks_(const Mobilogram& mobilogram)
    {
      for (Size i = 0; i < left_width_.size(); i++)
      {
        const int current_left_idx = left_width_[i];
        const int current_right_idx = right_width_[i];

        // Also integrate the intensities
        integrated_intensities_[i] = 0;
        for (int k = current_left_idx; k <= current_right_idx; k++)
        {
          integrated_intensities_[i] += mobilogram[k].getIntensity();
        }
      }
    }

    Size PeakPickerMobilogram::findClosestPeak_(const Mobilogram& mobilogram, double target_im, Size current_peak)
    {
      while (current_peak < mobilogram.size())
      {
        // check if we have walked past the RT of the peak
        if (target_im < mobilogram[current_peak].getMobility())
        {
          // see which one is closer, the current one or the one before
          if (current_peak > 0 &&
              std::fabs(target_im - mobilogram[current_peak - 1].getMobility()) <
                std::fabs(target_im - mobilogram[current_peak].getMobility()))
          {
            current_peak--;
          }

          return current_peak;
        }
        current_peak++;
      }
      return current_peak;
    }

    void PeakPickerMobilogram::removeOverlappingPeaks_(const Mobilogram& mobilogram, Mobilogram& picked_mobilogram)
    {
      if (picked_mobilogram.empty()) {return; }
      OPENMS_LOG_DEBUG << "Remove overlapping peaks now (size " << picked_mobilogram.size() << ")" << std::endl;
      Size current_peak = 0;
      // Find overlapping peaks
      for (Size i = 0; i < picked_mobilogram.size() - 1; i++)
      {
        // Check whether the current right overlaps with the next left
        // See whether we can correct this and find some border between the two
        // features ...
        if (right_width_[i] > left_width_[i + 1])
        {
          const int current_left_idx = left_width_[i];
          const int current_right_idx = right_width_[i];
          const int next_left_idx = left_width_[i + 1];
          const int next_right_idx = right_width_[i + 1];
          OPENMS_LOG_DEBUG << " Found overlapping " << i << " : " << current_left_idx << " " << current_right_idx << std::endl;
          OPENMS_LOG_DEBUG << "                   -- with  " << i + 1 << " : " << next_left_idx << " " << next_right_idx << std::endl;

          // Find the peak width and best RT
          double central_peak_mz = picked_mobilogram[i].getMobility();
          double next_peak_mz = picked_mobilogram[i + 1].getMobility();
          current_peak = findClosestPeak_(mobilogram, central_peak_mz, current_peak);
          Size next_peak = findClosestPeak_(mobilogram, next_peak_mz, current_peak);

          // adjust the right border of the current and left border of next
          Size k = 1;
          while ((current_peak + k) < mobilogram.size()
                 && (mobilogram[current_peak + k].getIntensity() < mobilogram[current_peak + k - 1].getIntensity()))
          {
            ++k;
          }
          Size new_right_border = current_peak + k - 1;
          k = 1;
          while ((next_peak - k + 1) > 0
                 && (mobilogram[next_peak - k].getIntensity() < mobilogram[next_peak - k + 1].getIntensity()))
          {
            ++k;
          }
          Size new_left_border = next_peak - k + 1;

          // assert that the peaks are now not overlapping any more ...
          if (new_left_border < new_right_border)
          {
            std::cerr << "Something went wrong, peaks are still overlapping!" << " - new left border " << new_left_border << " vs " << new_right_border << " -- will take the mean" << std::endl;
            new_left_border = (new_left_border + new_right_border) / 2;
            new_right_border = (new_left_border + new_right_border) / 2;

          }

          OPENMS_LOG_DEBUG << "New peak l: " << mobilogram[current_left_idx].getMobility() << " " << mobilogram[new_right_border].getMobility() << " int " << integrated_intensities_[i] << std::endl;
          OPENMS_LOG_DEBUG << "New peak r: " << mobilogram[new_left_border].getMobility() << " " << mobilogram[next_right_idx].getMobility() << " int " << integrated_intensities_[i + 1] << std::endl;


          right_width_[i] = new_right_border;
          left_width_[i + 1] = new_left_border;

        }
      }
    }

    void PeakPickerMobilogram::updateMembers_()
    {
        sgolay_frame_length_ = (UInt)param_.getValue("sgolay_frame_length");
        sgolay_polynomial_order_ = (UInt)param_.getValue("sgolay_polynomial_order");
        gauss_width_ = (double)param_.getValue("gauss_width");
        peak_width_ = (double)param_.getValue("peak_width");
        signal_to_noise_ = (double)param_.getValue("signal_to_noise");
        sn_win_len_ = (double)param_.getValue("sn_win_len");
        sn_bin_count_ = (UInt)param_.getValue("sn_bin_count");
        // TODO make list, not boolean
        use_gauss_ = (bool)param_.getValue("use_gauss").toBool();
        write_sn_log_messages_ = (bool)param_.getValue("write_sn_log_messages").toBool();
        method_ = (String)param_.getValue("method").toString();

        Param sg_filter_parameters = sgolay_.getParameters();
        sg_filter_parameters.setValue("frame_length", sgolay_frame_length_);
        sg_filter_parameters.setValue("polynomial_order", sgolay_polynomial_order_);
        sgolay_.setParameters(sg_filter_parameters);

        Param gfilter_parameters = gauss_.getParameters();
        gfilter_parameters.setValue("gaussian_width", gauss_width_);
        gauss_.setParameters(gfilter_parameters);

        Param snt_parameters = snt_.getParameters();
        snt_parameters.setValue("win_len", sn_win_len_);
        snt_parameters.setValue("bin_count", sn_bin_count_);
        snt_parameters.setValue("write_log_messages", param_.getValue("write_sn_log_messages"));
        snt_.setParameters(snt_parameters);
    }
}


