// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg $
// --------------------------------------------------------------------------

#pragma once
#include <OpenMS/QC/QCBase.h>

namespace OpenMS
{
  class FeatureMap;
  class MSExperiment;
  class TransformationDescription;

  /**
    @brief  QC metric to determine the number of MS2 scans per MS1 scan over RT

    TopNoverRT adds all unidentified MS2 scans to the unassignedPeptideIdentifications specifying the retention time value and
    setting the metaValues "ScanEventNumber" and "identified" for all PeptideIdentifications.

    "ScanEventNumber": number of the MS2 scans after the MS1 scan
    "identified": All PeptideIdentifications of the FeatureMap are marked with '+' and all unidentified MS2-Spectra with '-'.

    Once all this data is included, you can determine the number of MS2 scans after one MS1 scan and the part of identified MS2-Scans per "ScanEventNumber".
    **/
  class OPENMS_DLLAPI TopNoverRT : public QCBase
  {
  public:

    struct ScanEvent
    {
      ScanEvent(const UInt32 sen, const bool ms2_present)
        : scan_event_number(sen), ms2_presence(ms2_present) {}
      UInt32 scan_event_number;
      bool ms2_presence;
    };

    /// Constructor
    TopNoverRT() = default;

    /// Destructor
    virtual ~TopNoverRT() = default;

    /**
      @brief Calculate the ScanEventNumber, find all unidentified MS2-Spectra and add them to unassigned PeptideIdentifications,
             write meta values "ScanEventNumber" and "identified" in PeptideIdentification.
      @param exp Imported calibrated MzML file as MSExperiment
      @param features Imported featureXML file after FDR as FeatureMap
      @return unassigned peptide identifications newly generated from unidentified MS2-Spectra
      @throws MissingInformation If exp is empty
      @throws IllegalArgument If retention time of the MzML and featureXML file does not match
      @throws IllegalArgument If a peptide identification does not have a corresponding MS2 scan
    **/
    std::vector<PeptideIdentification> compute(const MSExperiment& exp, FeatureMap& features);

    /// returns the name of the metric
    const String& getName() const override;
    /// define the required input file: featureXML after FDR (=POSTFDRFEAT), MzML-file (MSExperiment) with all MS2-Spectra (=RAWMZML)
    Status requires() const override;

  private:

    /// name of the metric
    const String name_ = "TopNoverRT";
    /// ms2_included_ contains for every spectrum the information "ScanEventNumber" and presence MS2-scan in PeptideIDs
    std::vector<ScanEvent> ms2_included_{};

    /// compute "ScanEventNumber" for every spectrum: MS1=0, MS2=1-n, write into ms2_included_
    void setScanEventNumber_(const MSExperiment& exp);

    /// set ms2_included_ bool to true, if PeptideID exist and set "ScanEventNumber" for every PeptideID
    void setPresenceAndScanEventNumber_(PeptideIdentification& peptide_ID, const MSExperiment& exp);

    /// return all unidentified MS2-Scans as unassignedPeptideIDs, these contain only Information about RT and "ScanEventNumber"
    std::vector<PeptideIdentification> getUnassignedPeptideIdentifications_(const MSExperiment& exp);
  };
}
