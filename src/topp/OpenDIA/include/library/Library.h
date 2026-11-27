// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace OpenMS::OpenDIA
{
  /// One library fragment (a PQP TRANSITION row with its precursor's m/z and 1/K0).
  struct LibraryTransition
  {
    std::string transition_name;   ///< PQP TRANSITION.ID.
    std::string peptide_ref;       ///< PQP PRECURSOR.ID.
    double library_intensity = 0.0;
    double product_mz = 0.0;
    double precursor_mz = 0.0;
    double precursor_im = -1.0;    ///< Library 1/K0, -1 if absent.
    std::int8_t fragment_charge = 0;
    std::int16_t fragment_nr = -1;
    std::string fragment_type;     ///< Normalised ion type ("b", "y", "b-H2O", "prec", "" none, "unknown").
    bool decoy = false;

    bool getDecoy() const { return decoy; }
    void setDecoy(bool d) { decoy = d; }
    /// "y7^2", "b3", "prec", or "" without an ion type.
    std::string getAnnotation() const;
  };

  /// One library precursor.
  struct LibraryCompound
  {
    std::string id;                ///< PQP PRECURSOR.ID.
    double rt = std::numeric_limits<double>::quiet_NaN();  ///< PRECURSOR.LIBRARY_RT (0 when NULL in the file).
    double drift_time = -1.0;      ///< PRECURSOR.LIBRARY_DRIFT_TIME, -1 if absent.
    int charge = 0;
    std::string sequence;          ///< Modified sequence (unmodified if the modified one is empty); empty for a non-peptide compound.
    std::vector<std::string> protein_refs;
  };

  /**
    @brief An OpenDIA spectral library: precursors and fragments as read from a PQP or OpenDIA Parquet file.

    Only what OpenDIA uses is kept. Fragments are in the order of the file's transition query (explicitly sorted by
    precursor m/z, product m/z, RT, transition and precursor ID for libraries of at most 100,000 transitions).
    @p traml_ids holds (PRECURSOR.ID, PRECURSOR.TRAML_ID) for every precursor of the file, including precursors
    without fragments, in the file's PRECURSOR order (the insertion order readPrecursorTraMLIDs uses).
  */
  struct Library
  {
    std::vector<LibraryCompound> compounds;
    std::vector<LibraryTransition> transitions;
    std::vector<std::pair<std::string, std::string>> traml_ids;

    const std::vector<LibraryCompound>& getCompounds() const { return compounds; }
    const std::vector<LibraryTransition>& getTransitions() const { return transitions; }
  };

  /// Normalise a PQP ion type ("b-H20" is read as "b-H2O"; an unrecognised type becomes "unknown").
  std::string normaliseFragmentType(const std::string& type);

  /**
    @brief Read a PQP library (read-only SQLite access, the same query and row order as the OpenMS PQP reader).
    @throws Exception::FileNotReadable If the file cannot be opened.
    @throws Exception::InvalidValue If required tables are missing.
  */
  Library readLibrary(const std::string& pqp);

  /// PRECURSOR.ID to PRECURSOR.TRAML_ID of a PQP file.
  std::unordered_map<std::string, std::string> readPrecursorTraMLIDs(const std::string& pqp);

  /// The map readPrecursorTraMLIDs returns, built from Library::traml_ids (same insertion order, so the same iteration order).
  std::unordered_map<std::string, std::string> traMLIDMap(const std::vector<std::pair<std::string, std::string>>& traml_ids);

  /**
    @brief Write @p library as an OpenDIA Parquet library (schema version 1, zstd, row groups of at most 2^20 rows).

    One self-contained row per fragment in library order (precursor fields repeated, which zstd compresses away), so
    readParquetLibrary rebuilds the same Library and a later reader can stream row groups. Each row carries its
    precursor's position in Library::traml_ids. Entries of Library::traml_ids without fragments (and repeated
    precursor IDs) are one row each without a fragment (null transition_id), so they survive the round trip.
    @p metadata goes into the file footer (key, value). Writes to a temporary name next to @p path and renames, so an
    interrupted or concurrent write never leaves a partial file under @p path.
    @throws Exception::FileNotWritable If the file cannot be written.
    @throws Exception::InvalidValue If a fragment's precursor is not in Library::traml_ids.
  */
  void writeParquetLibrary(const Library& library, const std::string& path,
                           const std::vector<std::pair<std::string, std::string>>& metadata = {});

  /**
    @brief Read an OpenDIA Parquet library written by writeParquetLibrary.
    @throws Exception::FileNotReadable If the file cannot be opened.
    @throws Exception::InvalidValue If it is not an OpenDIA library or its schema version is not supported.
  */
  Library readParquetLibrary(const std::string& path);

  /// Footer metadata of an OpenDIA Parquet library; empty if @p path is missing or not a Parquet file.
  std::unordered_map<std::string, std::string> readParquetLibraryMetadata(const std::string& path);

  /// True if @p path names an OpenDIA Parquet library by its extension (.parquet).
  bool isParquetLibrary(const std::string& path);

  /**
    @brief Read the library at @p path, PQP or Parquet.

    A PQP is read and, without any option, converted once: its Parquet copy parquetCachePath(@p path) is written next
    to it and read instead on later calls while the PQP keeps its size and modification time. Libraries with fewer
    than @p cache_min_transitions fragments are not converted (they load in well under a second). If the copy cannot
    be written (read-only directory, full disk), the PQP is used as is. @p log receives one line per decision.
  */
  Library readAnyLibrary(const std::string& path, const std::function<void(const std::string&)>& log = {},
                         std::size_t cache_min_transitions = 100000);

  /// The Parquet copy readAnyLibrary uses for the PQP @p pqp: <directory>/<stem>.opendia.parquet.
  std::string parquetCachePath(const std::string& pqp);

  /// One library fragment as streamed by visitParquetLibrary / visitAnyLibrary (the fields of LibraryTransition).
  struct FragmentView
  {
    std::string_view precursor_id;
    std::string_view transition_id;
    std::string_view fragment_type;
    double precursor_mz = 0.0;
    double product_mz = 0.0;
    double library_intensity = 0.0;
    double precursor_im = -1.0;
    std::int8_t fragment_charge = 0;
    std::int16_t fragment_nr = -1;
    bool decoy = false;
  };

  /// The precursor of a streamed fragment (the fields its LibraryCompound gets when this is its first fragment).
  struct PrecursorView
  {
    std::string_view id;
    std::string_view sequence;
    double rt = 0.0;
    double drift_time = -1.0;
    int charge = 0;
    const std::vector<std::string_view>* proteins = nullptr;  ///< Never null.
  };

  /// getAnnotation() of a fragment given by its fields.
  std::string fragmentAnnotation(std::string_view type, std::int16_t ordinal, std::int8_t charge);

  /**
    @brief Receives the fragments of a library one at a time, in library order, without the library being held in
    memory (visitParquetLibrary, visitAnyLibrary).
  */
  class LibraryVisitor
  {
  public:
    virtual ~LibraryVisitor() = default;
    /// Number of fragments that follow, if known (a capacity hint).
    virtual void expect(std::size_t /*fragments*/) {}
    /// One fragment with its precursor. The views are only valid during the call.
    virtual void fragment(const FragmentView& fragment, const PrecursorView& precursor) = 0;
    /// Forget every fragment received so far (a Parquet copy that fails while being read is replaced by its PQP).
    virtual void reset() = 0;
  };

  /// Builds a Library from a stream of fragments (readParquetLibrary, readAnyLibrary).
  class LibraryBuilder final : public LibraryVisitor
  {
  public:
    Library library;  ///< traml_ids are set by the caller.
    void expect(std::size_t fragments) override;
    void fragment(const FragmentView& fragment, const PrecursorView& precursor) override;
    void reset() override;

  private:
    std::unordered_map<std::string, std::uint32_t> seen_;
  };

  /// Stream @p library into @p visitor in its fragment order (each fragment with the compound of its precursor).
  void visitLibrary(const Library& library, LibraryVisitor& visitor);

  /**
    @brief Stream an OpenDIA Parquet library into @p visitor (one row group in memory at a time).
    @return Library::traml_ids.
    @throws Exception::FileNotReadable, Exception::InvalidValue As readParquetLibrary.
  */
  std::vector<std::pair<std::string, std::string>> visitParquetLibrary(const std::string& path, LibraryVisitor& visitor);

  /**
    @brief readAnyLibrary, streamed: the same fragments in the same order, and the same Parquet copy decisions. A PQP
    without a usable copy is read (and converted) in full first; a Parquet library or copy is streamed.
    @return Library::traml_ids.
  */
  std::vector<std::pair<std::string, std::string>> visitAnyLibrary(const std::string& path, LibraryVisitor& visitor,
                                                                   const std::function<void(const std::string&)>& log = {},
                                                                   std::size_t cache_min_transitions = 100000);
}
