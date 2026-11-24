// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIALibrary.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <sqlite3.h>

#include <memory>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS::OpenDIA
{
  namespace
  {
    /// Libraries up to this many transitions are read in an explicit order (as the OpenMS PQP reader).
    constexpr std::size_t STABLE_ORDER_LIMIT = 100000;

    struct Database
    {
      sqlite3* db = nullptr;
      std::string path;

      explicit Database(const std::string& file) : path(file)
      {
        if (sqlite3_open_v2(file.c_str(), &db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK)
        {
          const std::string message = db ? sqlite3_errmsg(db) : "cannot open";
          sqlite3_close(db);
          throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file + ": " + message);
        }
      }
      ~Database() { sqlite3_close(db); }
      Database(const Database&) = delete;
      Database& operator=(const Database&) = delete;

      void exec(const std::string& sql)
      {
        char* error = nullptr;
        if (sqlite3_exec(db, sql.c_str(), nullptr, nullptr, &error) != SQLITE_OK)
        {
          const std::string message = error ? error : "error";
          sqlite3_free(error);
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path + ": " + message, sql);
        }
      }
    };

    struct Statement
    {
      sqlite3_stmt* stmt = nullptr;
      Statement(Database& d, const std::string& sql)
      {
        if (sqlite3_prepare_v2(d.db, sql.c_str(), -1, &stmt, nullptr) != SQLITE_OK)
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, d.path + ": " + sqlite3_errmsg(d.db), sql);
      }
      ~Statement() { sqlite3_finalize(stmt); }
      Statement(const Statement&) = delete;
      Statement& operator=(const Statement&) = delete;
      bool step() { return sqlite3_step(stmt) == SQLITE_ROW; }
      bool null(int c) const { return sqlite3_column_type(stmt, c) == SQLITE_NULL; }
      void get(int c, double& v) const { if (!null(c)) v = sqlite3_column_double(stmt, c); }
      void get(int c, int& v) const { if (!null(c)) v = sqlite3_column_int(stmt, c); }
      void get(int c, std::string& v) const { if (!null(c)) v = reinterpret_cast<const char*>(sqlite3_column_text(stmt, c)); }
    };

    bool tableExists(Database& d, const std::string& table)
    {
      Statement s(d, "SELECT 1 FROM sqlite_master WHERE type='table' AND name='" + table + "';");
      return s.step();
    }

    bool tableHasRows(Database& d, const std::string& table)
    {
      Statement s(d, "SELECT 1 FROM " + table + " LIMIT 1;");
      return s.step();
    }

    bool columnExists(Database& d, const std::string& table, const std::string& column)
    {
      Statement s(d, "PRAGMA table_info(" + table + ")");
      while (s.step())
      {
        std::string name;
        s.get(1, name);
        if (name == column) return true;
      }
      return false;
    }

    /// Same pragmas as the OpenMS PQP reader: a bounded page cache and mmap window for one read-only stream.
    void configure(Database& d)
    {
      d.exec("PRAGMA cache_size = -524288;");
      d.exec("PRAGMA mmap_size = 8589934592;");
      d.exec("PRAGMA query_only = ON;");
#ifdef _OPENMP
      d.exec("PRAGMA threads = " + std::to_string(std::max(1, omp_get_max_threads())) + ";");
#endif
    }
  }

  std::string normaliseFragmentType(const std::string& s)
  {
    if (s.empty()) return "";
    if (s == "a" || s == "b" || s == "c" || s == "x" || s == "y" || s == "z" || s == "z'" || s == "z." || s == "b-H2O" ||
        s == "y-H2O" || s == "b-NH3" || s == "y-NH3")
      return s;
    if (s == "prec" || s == "precursor") return "prec";
    if (s == "b-H20") return "b-H2O";
    if (s == "y-H20") return "y-H2O";
    return "unknown";
  }

  std::string fragmentAnnotation(std::string_view type, std::int16_t ordinal, std::int8_t charge)
  {
    if (type.empty()) return "";
    if (type == "prec") return "prec";
    std::string result(type);
    if (ordinal >= 0) result += std::to_string(ordinal);
    if (charge > 0) result += "^" + std::to_string(static_cast<int>(charge));
    return result;
  }

  std::string LibraryTransition::getAnnotation() const { return fragmentAnnotation(fragment_type, fragment_nr, fragment_charge); }

  Library readLibrary(const std::string& pqp)
  {
    Database d(pqp);
    configure(d);
    for (const char* table : {"PRECURSOR", "TRANSITION", "TRANSITION_PRECURSOR_MAPPING", "PRECURSOR_PEPTIDE_MAPPING", "PEPTIDE", "PROTEIN",
                              "PEPTIDE_PROTEIN_MAPPING"})
    {
      if (!tableExists(d, table)) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pqp + ": missing PQP table", table);
    }
    const bool compounds = tableExists(d, "PRECURSOR_COMPOUND_MAPPING") && tableHasRows(d, "PRECURSOR_COMPOUND_MAPPING");

    std::size_t transitions = 0;
    {
      Statement count(d, "SELECT COUNT(*) FROM TRANSITION;");
      if (count.step()) transitions = static_cast<std::size_t>(sqlite3_column_int64(count.stmt, 0));
    }
    // The query of the OpenMS PQP reader (same joins and, for small libraries, the same explicit order), so the row
    // order and with it every downstream ordering is unchanged. Unused columns are kept to keep the plan identical.
    const bool drift_time = columnExists(d, "PRECURSOR", "LIBRARY_DRIFT_TIME");
    const bool gene = tableExists(d, "GENE") && tableHasRows(d, "GENE");
    const bool peptidoforms = tableExists(d, "TRANSITION_PEPTIDE_MAPPING") && tableHasRows(d, "TRANSITION_PEPTIDE_MAPPING");
    const std::string annotation = columnExists(d, "TRANSITION", "ANNOTATION") ? "TRANSITION.ANNOTATION AS Annotation, " : "'' AS Annotation, ";
    std::string sql = "SELECT PRECURSOR.PRECURSOR_MZ AS precursor, TRANSITION.PRODUCT_MZ AS product, PRECURSOR.LIBRARY_RT AS rt_calibrated, "
                      "TRANSITION.ID AS transition_name, -1 AS CE, TRANSITION.LIBRARY_INTENSITY AS library_intensity, "
                      "PRECURSOR.ID AS group_id, TRANSITION.DECOY AS decoy, PEPTIDE.UNMODIFIED_SEQUENCE AS PeptideSequence, "
                      "PROTEIN_AGGREGATED.PROTEIN_ACCESSION AS ProteinName, " + annotation +
                      "PEPTIDE.MODIFIED_SEQUENCE AS FullPeptideName, NULL AS CompoundName, NULL AS SMILES, NULL AS SumFormula, "
                      "NULL AS Adducts, PRECURSOR.CHARGE AS precursor_charge, PRECURSOR.GROUP_LABEL AS peptide_group_label, "
                      "NULL AS label_type, TRANSITION.CHARGE AS fragment_charge, TRANSITION.ORDINAL AS fragment_nr, "
                      "NULL AS fragment_mzdelta, NULL AS fragment_modification, TRANSITION.TYPE AS fragment_type, NULL AS uniprot_id, "
                      "TRANSITION.DETECTING AS detecting_transition, TRANSITION.IDENTIFYING AS identifying_transition, "
                      "TRANSITION.QUANTIFYING AS quantifying_transition, ";
    sql += peptidoforms ? "PEPTIDE_AGGREGATED.PEPTIDOFORMS AS peptidoforms " : "NULL AS peptidoforms ";
    if (drift_time) sql += ", PRECURSOR.LIBRARY_DRIFT_TIME AS drift_time ";
    if (gene) sql += ", GENE_AGGREGATED.GENE_NAME AS gene_name ";
    sql += "FROM PRECURSOR ";
    if (gene)
    {
      sql += "INNER JOIN (SELECT PEPTIDE.ID AS PEPTIDE_ID, COALESCE(GROUP_CONCAT(GENE.GENE_NAME,';'),'NA') AS GENE_NAME FROM PEPTIDE "
             "LEFT OUTER JOIN PEPTIDE_GENE_MAPPING ON PEPTIDE.ID = PEPTIDE_GENE_MAPPING.PEPTIDE_ID "
             "LEFT OUTER JOIN GENE ON GENE.ID = PEPTIDE_GENE_MAPPING.GENE_ID GROUP BY PEPTIDE.ID) "
             "AS GENE_AGGREGATED ON PEPTIDE.ID = GENE_AGGREGATED.PEPTIDE_ID ";
    }
    sql += "INNER JOIN TRANSITION_PRECURSOR_MAPPING ON PRECURSOR.ID = TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID "
           "INNER JOIN TRANSITION ON TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID = TRANSITION.ID "
           "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID "
           "INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID "
           "INNER JOIN (SELECT PEPTIDE_ID, GROUP_CONCAT(PROTEIN_ACCESSION,';') AS PROTEIN_ACCESSION FROM PROTEIN "
           "INNER JOIN PEPTIDE_PROTEIN_MAPPING ON PROTEIN.ID = PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID GROUP BY PEPTIDE_ID) "
           "AS PROTEIN_AGGREGATED ON PEPTIDE.ID = PROTEIN_AGGREGATED.PEPTIDE_ID ";
    if (peptidoforms)
    {
      sql += "LEFT OUTER JOIN (SELECT TRANSITION_ID, GROUP_CONCAT(MODIFIED_SEQUENCE,'|') AS PEPTIDOFORMS FROM TRANSITION_PEPTIDE_MAPPING "
             "INNER JOIN PEPTIDE ON TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID GROUP BY TRANSITION_ID) "
             "AS PEPTIDE_AGGREGATED ON TRANSITION.ID = PEPTIDE_AGGREGATED.TRANSITION_ID ";
    }
    if (compounds)
    {
      // Non-peptide compounds (as the OpenMS reader's second branch): no sequence or proteins.
      const std::string adducts = columnExists(d, "COMPOUND", "ADDUCTS") ? "COMPOUND.ADDUCTS AS Adducts, " : "'' AS Adducts, ";
      sql += "UNION ALL SELECT PRECURSOR.PRECURSOR_MZ AS precursor, TRANSITION.PRODUCT_MZ AS product, PRECURSOR.LIBRARY_RT AS rt_calibrated, "
             "TRANSITION.ID AS transition_name, -1 AS CE, TRANSITION.LIBRARY_INTENSITY AS library_intensity, "
             "PRECURSOR.ID AS group_id, TRANSITION.DECOY AS decoy, NULL AS PeptideSequence, NULL AS ProteinName, " + annotation +
             "NULL AS FullPeptideName, COMPOUND.COMPOUND_NAME AS CompoundName, COMPOUND.SMILES AS SMILES, "
             "COMPOUND.SUM_FORMULA AS SumFormula, " + adducts +
             "PRECURSOR.CHARGE AS precursor_charge, PRECURSOR.GROUP_LABEL AS peptide_group_label, NULL AS label_type, "
             "TRANSITION.CHARGE AS fragment_charge, TRANSITION.ORDINAL AS fragment_nr, NULL AS fragment_mzdelta, "
             "NULL AS fragment_modification, TRANSITION.TYPE AS fragment_type, NULL AS uniprot_id, "
             "TRANSITION.DETECTING AS detecting_transition, TRANSITION.IDENTIFYING AS identifying_transition, "
             "TRANSITION.QUANTIFYING AS quantifying_transition, NULL AS peptidoforms ";
      if (drift_time) sql += ", PRECURSOR.LIBRARY_DRIFT_TIME AS drift_time ";
      if (gene) sql += ", 'NA' AS gene_name ";
      sql += "FROM PRECURSOR "
             "INNER JOIN TRANSITION_PRECURSOR_MAPPING ON PRECURSOR.ID = TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID "
             "INNER JOIN TRANSITION ON TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID = TRANSITION.ID "
             "INNER JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID "
             "INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID ";
    }
    if (transitions <= STABLE_ORDER_LIMIT) sql += "ORDER BY precursor, product, rt_calibrated, transition_name, group_id ";
    sql += "; ";

    Library lib;
    lib.transitions.reserve(transitions);
    lib.compounds.reserve(transitions / 6 + 1);
    std::unordered_set<std::string> seen;
    seen.reserve(transitions / 6 + 1);
    Statement s(d, sql);
    while (s.step())
    {
      if (s.null(0)) break;
      LibraryTransition t;
      double rt = 0.0, im = -1.0;
      int decoy = 0, charge = 0, fragment_charge = 0, fragment_nr = -1;
      std::string peptide, proteins, modified, type, compound_name;
      s.get(0, t.precursor_mz);
      s.get(1, t.product_mz);
      s.get(2, rt);
      s.get(3, t.transition_name);
      s.get(5, t.library_intensity);
      s.get(6, t.peptide_ref);
      s.get(7, decoy);
      s.get(8, peptide);
      s.get(9, proteins);
      s.get(11, modified);
      s.get(12, compound_name);
      s.get(16, charge);
      s.get(19, fragment_charge);
      s.get(20, fragment_nr);
      s.get(23, type);
      if (drift_time) s.get(29, im);
      t.precursor_im = im;
      t.fragment_charge = static_cast<std::int8_t>(fragment_charge);
      t.fragment_nr = static_cast<std::int16_t>(fragment_nr);
      t.fragment_type = normaliseFragmentType(type);
      t.decoy = decoy != 0;
      if (seen.insert(t.peptide_ref).second)
      {
        LibraryCompound c;
        c.id = t.peptide_ref;
        c.rt = rt;
        c.drift_time = im;
        c.charge = charge;
        if (compound_name.empty())
        {
          c.sequence = modified.empty() ? peptide : modified;
          if (!proteins.empty())
          {
            std::vector<std::string> parts;
            StringUtils::split(proteins, ';', parts);
            c.protein_refs.assign(parts.begin(), parts.end());
          }
        }
        lib.compounds.push_back(std::move(c));
      }
      lib.transitions.push_back(std::move(t));
    }
    Statement ids(d, "SELECT ID, TRAML_ID FROM PRECURSOR;");
    while (ids.step())
    {
      std::string id, traml;
      ids.get(0, id);
      ids.get(1, traml);
      lib.traml_ids.emplace_back(std::move(id), std::move(traml));
    }
    return lib;
  }

  void LibraryBuilder::expect(std::size_t fragments)
  {
    library.transitions.reserve(fragments);
    library.compounds.reserve(fragments / 6 + 1);
  }

  void LibraryBuilder::fragment(const FragmentView& f, const PrecursorView& p)
  {
    LibraryTransition t;
    t.transition_name = f.transition_id;
    t.peptide_ref = f.precursor_id;
    t.library_intensity = f.library_intensity;
    t.product_mz = f.product_mz;
    t.precursor_mz = f.precursor_mz;
    t.precursor_im = f.precursor_im;
    t.fragment_charge = f.fragment_charge;
    t.fragment_nr = f.fragment_nr;
    t.fragment_type = f.fragment_type;
    t.decoy = f.decoy;
    if (seen_.try_emplace(t.peptide_ref, static_cast<std::uint32_t>(library.compounds.size())).second)
    {
      LibraryCompound c;
      c.id = p.id;
      c.rt = p.rt;
      c.drift_time = p.drift_time;
      c.charge = p.charge;
      c.sequence = p.sequence;
      c.protein_refs.assign(p.proteins->begin(), p.proteins->end());
      library.compounds.push_back(std::move(c));
    }
    library.transitions.push_back(std::move(t));
  }

  void LibraryBuilder::reset()
  {
    library = Library();
    seen_.clear();
  }

  void visitLibrary(const Library& library, LibraryVisitor& visitor)
  {
    std::unordered_map<std::string_view, std::uint32_t> compound;
    compound.reserve(library.compounds.size());
    for (std::size_t i = 0; i < library.compounds.size(); ++i) compound.try_emplace(library.compounds[i].id, static_cast<std::uint32_t>(i));
    visitor.expect(library.transitions.size());
    std::vector<std::string_view> proteins;
    std::uint32_t proteins_of = UINT32_MAX;
    for (const auto& t : library.transitions)
    {
      const auto it = compound.find(t.peptide_ref);
      if (it == compound.end())
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Fragment of a precursor without compound", t.peptide_ref);
      const LibraryCompound& c = library.compounds[it->second];
      if (it->second != proteins_of)
      {
        proteins.assign(c.protein_refs.begin(), c.protein_refs.end());
        proteins_of = it->second;
      }
      FragmentView f;
      f.precursor_id = t.peptide_ref;
      f.transition_id = t.transition_name;
      f.fragment_type = t.fragment_type;
      f.precursor_mz = t.precursor_mz;
      f.product_mz = t.product_mz;
      f.library_intensity = t.library_intensity;
      f.precursor_im = t.precursor_im;
      f.fragment_charge = t.fragment_charge;
      f.fragment_nr = t.fragment_nr;
      f.decoy = t.decoy;
      PrecursorView p;
      p.id = c.id;
      p.sequence = c.sequence;
      p.rt = c.rt;
      p.drift_time = c.drift_time;
      p.charge = c.charge;
      p.proteins = &proteins;
      visitor.fragment(f, p);
    }
  }

  std::unordered_map<std::string, std::string> traMLIDMap(const std::vector<std::pair<std::string, std::string>>& traml_ids)
  {
    // No reserve: the same inserts as readPrecursorTraMLIDs give the same bucket layout and iteration order.
    std::unordered_map<std::string, std::string> out;
    for (const auto& [id, traml] : traml_ids) out.emplace(id, traml);
    return out;
  }

  std::unordered_map<std::string, std::string> readPrecursorTraMLIDs(const std::string& pqp)
  {
    Database d(pqp);
    Statement s(d, "SELECT ID, TRAML_ID FROM PRECURSOR;");
    std::unordered_map<std::string, std::string> out;
    while (s.step())
    {
      std::string id, traml;
      s.get(0, id);
      s.get(1, traml);
      out.emplace(std::move(id), std::move(traml));
    }
    return out;
  }
}
