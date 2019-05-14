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
// $Authors: Swenja Wagner, Patricia Scheil $
// --------------------------------------------------------------------------

#include <OpenMS/QC/QCBase.h>

namespace OpenMS
{
  const std::string QCBase::names_of_requires[] = {"fail", "raw.mzML", "postFDR.featureXML", "preFDR.featureXML", "contaminants.fasta", "trafoAlign.trafoXML"};

  QCBase::SpectraMap::SpectraMap() = default;

  QCBase::SpectraMap::SpectraMap(const MSExperiment& exp)
  {
    calculateMap(exp);
  }

  QCBase::SpectraMap::~SpectraMap() = default;

  void QCBase::SpectraMap::calculateMap(const MSExperiment& exp)
  {
    map_to_index_.clear();
    for (Size i = 0; i < exp.size(); ++i)
    {
      map_to_index_[exp[i].getNativeID()] = i;
    }
  }

  const UInt64& QCBase::SpectraMap::at(const String& identifier) const
  {
    if (map_to_index_.find(identifier) == map_to_index_.end())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No spectrum with given identifier in MSExperiment!");
    }
    return map_to_index_.at(identifier);
  }

  void QCBase::SpectraMap::clear()
  {
    map_to_index_.clear();
  }

  bool QCBase::SpectraMap::empty() const
  {
    return map_to_index_.empty();
  }
  
  Size QCBase::SpectraMap::size() const
  {
    return map_to_index_.size();
  }
} //namespace OpenMS
