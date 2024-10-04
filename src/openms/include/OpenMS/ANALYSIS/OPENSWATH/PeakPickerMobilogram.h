// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/KERNEL/MobilityPeak1D.h>

#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>

#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <vector>

namespace OpenMS
{
    /**
    @brief The PeakPickerChromatogram finds peaks a single mobilogram.

    @htmlinclude OpenMS_PeakPickerChromatogram.parameters

    It uses the PeakPickerHiRes internally to find interesting seed candidates.
    These candidates are then expanded and a right/left border of the peak is
    searched.
    Additionally, overlapping peaks can be removed.

    */

    class OPENMS_DLLAPI PeakPickerMobilogram :
            public DefaultParamHandler
    {
    public:
        //@{
        /// Constructor
        PeakPickerMobilogram();

        /// Destructor
        ~PeakPickerMobilogram() override {}
        //@}

        /// indices into FloatDataArrays of resulting picked mobilogram
        enum FLOATINDICES { IDX_FWHM = 0, IDX_ABUNDANCE = 1, IDX_LEFTBORDER = 2, IDX_RIGHTBORDER = 3, SIZE_OF_FLOATINDICES };

        /**
            @brief Finds peaks in a single mobilogram and annotates left/right borders

            It uses a modified algorithm of the PeakPickerHiRes

            This function will return a picked mobilogram
        */
        void pickMobilogram(const Mobilogram& mobilogram, Mobilogram& picked_mobilogram);

        /**
            @brief Finds peaks in a single mobilogram and annotates left/right borders

            It uses a modified algorithm of the PeakPickerHiRes

            This function will return a picked mobilogram and a smoothed mobilogram
        */
        void pickMobilogram(Mobilogram mobilogram, Mobilogram& picked_mobilogram, Mobilogram& smoothed_mobilogram);

        /// Temporary vector to hold the integrated intensities
        std::vector<double> integrated_intensities_;

        /// Temporary vector to hold the peak left widths
        std::vector<Size> left_width_;

        /// Temporary vector to hold the peak right widths
        std::vector<Size> right_width_;

      protected:
        void pickMobilogram_(const Mobilogram& mobilogram, Mobilogram& picked_mobilogram);

        /**
          @brief Compute peak area (peak integration)
        */
        void integratePeaks_(const Mobilogram& mobilogram);

        /**
          @brief Helper function to find the closest peak in a mobilogram to "target_im"

          The search will start from the index current_peak, so the function is
          assuming the closest peak is to the right of current_peak.

          It will return the index of the closest peak in the mobilogram.
        */
        Size findClosestPeak_(const Mobilogram& mobilogram, double target_im, Size current_peak = 0);

        /**
          @brief Helper function to remove overlapping peaks in a single Chromatogram

        */
        void removeOverlappingPeaks_(const Mobilogram& mobilogram, Mobilogram& picked_mobilogram);

        /// Synchronize members with param class
        void updateMembers_() override;

        // Members
        /// Frame length for the SGolay smoothing
        UInt sgolay_frame_length_;
        /// Polynomial order for the SGolay smoothing
        UInt sgolay_polynomial_order_;
        /// Width of the Gaussian smoothing
        double gauss_width_;
        /// Whether to use Gaussian smoothing
        bool use_gauss_;
        /// Whether to resolve overlapping peaks
        bool remove_overlapping_;

        /// Forced peak with
        double peak_width_;
        /// Signal to noise threshold
        double signal_to_noise_;

        /// Signal to noise window length
        double sn_win_len_;
        /// Signal to noise bin count
        UInt sn_bin_count_;
        /// Whether to write out log messages of the SN estimator
        bool write_sn_log_messages_;
        /// Peak picker method
        String method_;

        PeakPickerHiRes pp_;
        SavitzkyGolayFilter sgolay_;
        GaussFilter gauss_;
        SignalToNoiseEstimatorMedian<Mobilogram > snt_;
    };
}


