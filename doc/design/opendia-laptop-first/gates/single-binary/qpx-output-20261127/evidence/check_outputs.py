#!/usr/bin/env python3
"""Independent OpenDIA output acceptance checks. Requires pyarrow, duckdb, numpy.

Usage: check_outputs.py OUTPUT_DIR [--refdir BASELINE_DIR]
       [--peakgroup-fdr .01 --peptide-fdr .01 --protein-fdr .01]

This validates stored evidence and algebra, not empirical FDR calibration or raw-point
provenance. Matrix gates are independent by entity level. Per-feature pass flags require
peak-group and the corresponding entity-level cuts. Missing evidence causes a failure
or an explicit unavailable message,
never an inferred success. DuckDB spills large joins to scratch instead of building
Python dictionaries containing the complete candidate population.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


REPORTS = ("precursors", "precursor_matrix", "peptides", "peptide_matrix",
           "protein_groups", "protein_group_matrix")
PRODUCER_COMPOSITE = ("run_file_name", "opendia_precursor_id", "opendia_candidate_row")
EXTENSION_TYPES = {
    "opendia_precursor_id": "VARCHAR", "opendia_candidate_row": "BIGINT",
    "opendia_winner": "BOOLEAN", "opendia_score_kept": "BOOLEAN",
    "opendia_network_score_kept": "BOOLEAN", "opendia_decoy_class": "VARCHAR",
}
for _extension in ("rt_apex", "rt_start", "rt_end", "area", "mz", "library_rt", "score", "qvalue",
                   "network_score", "network_qvalue", "peptide_global_qvalue", "pg_qvalue"):
    EXTENSION_TYPES["opendia_" + _extension] = "DOUBLE"


def producer_id(values, unordered=()):
    normalized = list(values)
    for position in unordered:
        encoded = {json.dumps(value, separators=(",", ":"), sort_keys=True, ensure_ascii=True): value
                   for value in normalized[position]}
        normalized[position] = [encoded[key] for key in sorted(encoded)]
    data = json.dumps(normalized, separators=(",", ":"), sort_keys=True, ensure_ascii=True).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(data, digest_size=8).digest(), "big", signed=True)


def check_producer_ids(path):
    source = pq.ParquetFile(path)
    metadata = source.schema_arrow.metadata or {}
    names = source.schema_arrow.names
    composite = metadata.get(b"opendia_identity_composite", metadata.get(b"identity_composite", b"")).decode().split(",")
    require(tuple(composite) == PRODUCER_COMPOSITE, f"{path}: unknown producer identity composite {composite}")
    if "x_opendia_precursor_id" in names:
        require(metadata.get(b"identity_source") == b"supplied", f"{path}: canonical opaque IDs must declare identity_source=supplied")
        require(b"opendia_identity_composite" in metadata, f"{path}: missing producer composite footer")
    physical = [name if name in names else "x_" + name for name in PRODUCER_COMPOSITE]
    checked = 0
    for batch in source.iter_batches(batch_size=8192, columns=["feature_id", *physical]):
        for row in batch.to_pylist():
            values = [row[name] for name in physical]
            require(all(value is not None for value in values), f"{path}: missing producer identity input")
            values[2] = int(values[2])
            require(row["feature_id"] == producer_id(values), f"{path}: independent producer feature_id hash mismatch at row {checked}")
            checked += 1
    print(json.dumps({"independent_producer_ids": checked, "file": str(path)}, sort_keys=True))


def check_pg_ids(path):
    source = pq.ParquetFile(path)
    checked = 0
    for batch in source.iter_batches(batch_size=8192, columns=["pg_id", "pg_accessions", "grouped_runs", "label"]):
        for row in batch.to_pylist():
            values = [row["pg_accessions"], row["grouped_runs"], row["label"]]
            require(row["pg_id"] == producer_id(values, (0, 1)), f"{path}: independent protein-group pg_id hash mismatch at row {checked}")
            checked += 1
    print(json.dumps({"independent_pg_ids": checked, "file": str(path)}, sort_keys=True))



def require(condition, message):
    if not condition:
        raise AssertionError(message)


def literal(value):
    return "'" + str(value).replace("'", "''") + "'"


def identifier(value):
    return '"' + str(value).replace('"', '""') + '"'


def gate(run_q, global_q, threshold):
    return (f"coalesce(isfinite({run_q}) AND isfinite({global_q}) "
            f"AND {run_q} BETWEEN 0 AND {threshold!r} "
            f"AND {global_q} BETWEEN 0 AND {threshold!r}, false)")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def tsv_parquet_roundtrip(tsv, parquet):
    """Compare each scalar exactly using the authoritative Parquet type, in row order."""
    source = pq.ParquetFile(parquet)
    nulls = {name: 0 for name in source.schema_arrow.names}
    count = 0
    with tsv.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t", quotechar='"', strict=True)
        require(next(reader, None) == source.schema_arrow.names, f"{tsv}: header/schema mismatch")
        for batch in source.iter_batches(batch_size=4096):
            columns = [column.to_pylist() for column in batch.columns]
            for index in range(batch.num_rows):
                text = next(reader, None)
                require(text is not None and len(text) == len(columns), f"{tsv}: incomplete row {count}")
                for field, values, value in zip(batch.schema, columns, text):
                    expected = values[index]
                    where = f"{tsv.name} row={count} column={field.name}"
                    if expected is None:
                        nulls[field.name] += 1
                        require(value == "NA", f"{where}: null is not NA")
                    elif pa.types.is_floating(field.type):
                        require(math.isfinite(expected), f"{where}: nonfinite report value instead of null")
                        require(float(value) == expected, f"{where}: numeric roundtrip {value!r} != {expected!r}")
                    elif pa.types.is_integer(field.type):
                        require(int(value) == expected, f"{where}: integer mismatch")
                    elif pa.types.is_boolean(field.type):
                        require(value == ("true" if expected else "false"), f"{where}: boolean mismatch")
                    elif pa.types.is_string(field.type):
                        require(value == expected, f"{where}: string/escaping mismatch")
                    else:
                        raise AssertionError(f"{where}: unsupported report type {field.type}")
                count += 1
        require(next(reader, None) is None, f"{tsv}: extra rows")
    require(count == source.metadata.num_rows, f"{tsv}: row count mismatch")
    for group in range(source.metadata.num_row_groups):
        require(source.metadata.row_group(group).num_rows <= 4096, f"{parquet}: row group exceeds 4096")
    print(json.dumps({"table": parquet.stem, "rows": count, "nulls": nulls,
                      "tsv_bytes": tsv.stat().st_size, "parquet_bytes": parquet.stat().st_size}, sort_keys=True))


class Checks:
    def __init__(self, connection):
        self.db = connection

    def count(self, query):
        return self.db.execute(query).fetchone()[0]

    def zero(self, query, message):
        sample = self.db.execute(query + " LIMIT 5").fetchall()
        require(not sample, f"{message}: sample={sample!r}")
        print("PASS", message)

    def unique(self, table, keys):
        names = ",".join(identifier(key) for key in keys)
        null_key = " OR ".join(identifier(key) + " IS NULL" for key in keys)
        self.zero(f"SELECT {names} FROM {table} WHERE {null_key}", f"nonnull identity {table}({','.join(keys)})")
        self.zero(f"SELECT {names},count(*) FROM {table} GROUP BY {names} HAVING count(*) != 1",
                  f"unique {table}({','.join(keys)})")

    def parquet_view(self, name, paths):
        require(paths, f"No Parquet files for {name}")
        expression = "[" + ",".join(literal(path) for path in paths) + "]"
        self.db.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet({expression}, union_by_name=true)")

    def canonical_view(self, name, path):
        self.parquet_view(name + "_raw", [path])
        fields = pq.read_schema(path)
        projection, aliases = [], set()
        encoded = 0
        for field in fields:
            target = field.name
            expression = identifier(field.name)
            if field.name.startswith("x_opendia_"):
                require(pa.types.is_string(field.type) and field.nullable,
                        f"{path}: canonical extension {field.name} must be nullable UTF8")
                extension = field.name[2:]
                target = {"opendia_peptide_global_qvalue": "peptide_global_qvalue",
                          "opendia_pg_qvalue": "pg_qvalue"}.get(extension, extension)
                expression = f"cast({expression} AS {EXTENSION_TYPES.get(extension, 'VARCHAR')})"
                encoded += 1
            require(target not in aliases, f"{path}: ambiguous typed/canonical extension alias {target}")
            aliases.add(target)
            projection.append(expression + " AS " + identifier(target))
        self.db.execute(f"CREATE VIEW {name} AS SELECT {','.join(projection)} FROM {name}_raw")
        print(json.dumps({"canonical_view": name, "string_extensions_decoded": encoded}, sort_keys=True))

    def same(self, left, right, keys, pairs, label, rtol=0.0, atol=0.0):
        using = ",".join(identifier(key) for key in keys)
        mismatch = []
        for lcol, rcol in pairs:
            left_col, right_col = f"a.{identifier(lcol)}", f"b.{identifier(rcol)}"
            if rtol or atol:
                mismatch.append(f"(({left_col} IS NULL) != ({right_col} IS NULL) OR "
                                f"abs({left_col}-{right_col}) > {atol!r} + {rtol!r}*abs({right_col}))")
            else:
                mismatch.append(f"{left_col} IS DISTINCT FROM {right_col}")
        condition = " OR ".join(mismatch) or "false"
        first = identifier(keys[0])
        self.zero(f"SELECT coalesce(a.{first},b.{first}) FROM {left} a FULL OUTER JOIN {right} b USING ({using}) "
                  f"WHERE a.{first} IS NULL OR b.{first} IS NULL OR {condition}", label)

    def matrix(self, stem, id_column, expected, runs, metadata):
        self.unique(stem, [id_column])
        columns = [row[0] for row in self.db.execute(f"DESCRIBE {stem}").fetchall()]
        require(columns == metadata + runs, f"{stem}: metadata/run column order differs: {columns}")
        if runs:
            values = ",".join(f"({literal(run)},m.{identifier(run)})" for run in runs)
            self.db.execute(f"CREATE VIEW {stem}_cells AS SELECT m.{identifier(id_column)} AS id, v.run, v.area "
                            f"FROM {stem} m CROSS JOIN LATERAL (VALUES {values}) v(run,area) WHERE v.area IS NOT NULL")
        else:
            self.db.execute(f"CREATE VIEW {stem}_cells AS SELECT cast(NULL AS VARCHAR) id, cast(NULL AS VARCHAR) run, "
                            "cast(NULL AS DOUBLE) area WHERE false")
        self.unique(expected, ["id", "run"])
        self.same(f"{stem}_cells", expected, ["id", "run"], [("area", "area")], f"{stem}: exact eligible cells and real zeros")
        self.zero(f"SELECT m.{identifier(id_column)} FROM {stem} m ANTI JOIN "
                  f"(SELECT DISTINCT id FROM {expected}) e ON m.{identifier(id_column)}=e.id",
                  f"{stem}: no rows without eligible cells")
        print(json.dumps({"matrix": stem, "eligible_cells": self.count(f"SELECT count(*) FROM {expected}"),
                          "real_zero_cells": self.count(f"SELECT count(*) FROM {expected} WHERE area=0")}, sort_keys=True))


def check_legacy(check, output, refdir, rtol, atol):
    legacy = output / "legacy"
    files = sorted(legacy.glob("*.precursors.tsv"))
    require(files, f"Missing mandatory per-run legacy evidence in {legacy}")
    # std::to_chars without precision prints shortest round-trip numbers, so exact equality is
    # the default. Nonzero tolerances must be explicitly requested and are printed in the summary.
    reads = [f"SELECT *, row_number() OVER () AS input_order FROM read_csv({literal(path)}, "
             "delim='\t',header=true,all_varchar=true)" for path in files]
    check.db.execute("CREATE VIEW legacy_candidates AS " + " UNION ALL ".join(reads))
    check.db.execute("CREATE VIEW legacy_winners AS SELECT run,precursor_id,sequence,cast(charge AS BIGINT) charge, "
                     "cast(candidate_row AS BIGINT) candidate_row,cast(is_decoy AS BOOLEAN) decoy," +
                     ",".join(f"cast({name} AS DOUBLE) {name}" for name in ("rt_apex", "rt_start", "rt_end", "area", "score", "q_value")) +
                     " FROM legacy_candidates QUALIFY row_number() OVER "
                     "(PARTITION BY run,precursor_id ORDER BY cast(score AS DOUBLE) DESC,input_order)=1")
    check.unique("legacy_winners", ["run", "precursor_id"])
    check.db.execute("CREATE VIEW retained_candidate_ids AS SELECT run_file_name run,opendia_precursor_id precursor_id, "
                     "opendia_candidate_row candidate_row FROM feature_winners")
    check.same("retained_candidate_ids", "legacy_winners", ["run", "precursor_id"], [("candidate_row", "candidate_row")],
               "retained candidate identity versus legacy candidate rows")
    check.same("precursors", "legacy_winners", ["run", "precursor_id"],
               [("sequence", "sequence"), ("charge", "charge"), ("decoy", "decoy")], "retained identity and metadata versus legacy")
    check.same("precursors", "legacy_winners", ["run", "precursor_id"],
               [("area", "area"), ("rt_apex", "rt_apex"), ("rt_start", "rt_start"), ("rt_end", "rt_end"),
                ("discriminant_score", "score"), ("run_q_value", "q_value")],
               "retained full-precision area/RT/score/q versus legacy", rtol, atol)
    combined = legacy / "precursors.tsv"
    require(combined.is_file(), "Missing legacy combined precursors.tsv")
    check.db.execute(f"CREATE VIEW legacy_combined AS SELECT precursor_id id,run,cast(area AS DOUBLE) area "
                     f"FROM read_csv({literal(combined)},delim='\t',header=true,all_varchar=true)")
    check.unique("legacy_combined", ["id", "run"])
    check.same("precursor_matrix_cells", "legacy_combined", ["id", "run"], [("area", "area")],
               "precursor matrix exactly equals legacy combined target cells", rtol, atol)
    if refdir:
        references = [path for path in refdir.iterdir() if path.name.endswith((".precursors.tsv", ".calibration.json")) or path.name == "precursors.tsv"]
        require(references, f"No reference output files in {refdir}")
        for reference in sorted(references):
            actual = output / reference.name if reference.name.endswith(".calibration.json") else legacy / reference.name
            require(actual.is_file(), f"Missing reference counterpart {actual}")
            require(reference.stat().st_size == actual.stat().st_size and sha256(reference) == sha256(actual),
                    f"Baseline byte identity failed: {reference.name}")
            print("BYTE_IDENTICAL", reference.name)


def check_xipm(check, output, required):
    files = sorted(output.rglob("*.xipm"))
    require(files or not required, "--require-xipm requested but no XIPM files exist")
    if not files:
        print("UNAVAILABLE optional XIPM checks: no .xipm files")
        return
    check.parquet_view("xipm", files)
    check.zero("SELECT x.feature_id FROM xipm x LEFT JOIN feature f USING(feature_id) "
               "WHERE f.feature_id IS NULL OR NOT f.opendia_winner "
               "OR x.run_file_name IS DISTINCT FROM f.run_file_name "
               "OR x.opendia_precursor_id IS DISTINCT FROM f.opendia_precursor_id "
               "OR x.opendia_candidate_row IS DISTINCT FROM f.opendia_candidate_row",
               "every XIPM map links to the actual exported QPX winner candidate")
    for path in files:
        counts = {"maps": 0, "points": 0, "missing_mobility_points": 0, "ms1_maps": 0, "ms2_maps": 0}
        for batch in pq.ParquetFile(path).iter_batches(batch_size=1):
            for row in batch.to_pylist():
                arrays = {}
                for prefix in ("MZ", "RT", "MOBILITY", "INTENSITY"):
                    require(row[prefix + "_COMPRESSION"] == 0, f"{path}: unsupported compressed XIPM data")
                    data = row[prefix + "_DATA"]
                    require(data is not None and len(data) % 8 == 0, f"{path}: invalid {prefix} binary length")
                    arrays[prefix] = np.frombuffer(data, dtype="=f8")
                length = len(arrays["MZ"])
                require(all(len(value) == length for value in arrays.values()), f"{path}: unequal XIPM point array lengths")
                for name in ("MZ", "RT", "INTENSITY"):
                    require(np.isfinite(arrays[name]).all(), f"{path}: nonfinite {name} point")
                # Native OpenMS uses NaN for unobserved mobility, not a fabricated numeric value.
                require(not np.isinf(arrays["MOBILITY"]).any(), f"{path}: infinite mobility")
                require((arrays["MZ"] > 0).all() and (arrays["INTENSITY"] >= 0).all(), f"{path}: invalid point coordinates/intensities")
                require(((arrays["RT"] >= row["RT_START"]) & (arrays["RT"] <= row["RT_END"])).all(), f"{path}: point outside recorded RT interval")
                require(row["MS_LEVEL"] in (1, 2), f"{path}: unexpected MS level")
                require(row["PEAKMAP_TYPE"] == ("precursor" if row["MS_LEVEL"] == 1 else "transition"), f"{path}: map type/MS level mismatch")
                counts["maps"] += 1
                counts["points"] += length
                counts["missing_mobility_points"] += int(np.isnan(arrays["MOBILITY"]).sum())
                counts[f"ms{row['MS_LEVEL']}_maps"] += 1
        print(json.dumps({"xipm": str(path), "bytes": path.stat().st_size, **counts}, sort_keys=True))


def run(args):
    output = args.output.resolve()
    require(output.is_dir(), f"Output directory does not exist: {output}")
    qpx = output / "qpx"
    feature_path = qpx / "opendia.feature.parquet"
    require(feature_path.is_file(), f"Missing canonical QPX feature view: {feature_path}")
    metadata = pq.read_schema(feature_path).metadata or {}
    cuts = {}
    for level in ("peakgroup", "peptide", "protein"):
        supplied = getattr(args, level + "_fdr")
        key = ("opendia_" + level + "_fdr").encode()
        require(supplied is not None or key in metadata, f"No {level} threshold argument or stored metadata")
        value = supplied if supplied is not None else float(metadata[key])
        require(math.isfinite(value) and 0 <= value <= 1, f"Invalid {level} threshold")
        if key in metadata:
            require(value == float(metadata[key]), f"Supplied {level} threshold contradicts stored metadata")
        cuts[level] = value
    print(json.dumps({"output": str(output), "thresholds": cuts,
                      "legacy_rtol": args.legacy_rtol, "legacy_atol": args.legacy_atol}, sort_keys=True))
    for path in output.rglob("*"):
        require(".partial" not in path.name, f"Orphan partial output: {path}")
        require(not path.name.startswith("opendia_cache_"), f"Cache directory not removed: {path}")
    for stem in REPORTS:
        for extension in ("tsv", "parquet"):
            require((output / f"{stem}.{extension}").is_file(), f"Missing {stem}.{extension}")
        tsv_parquet_roundtrip(output / f"{stem}.tsv", output / f"{stem}.parquet")

    with tempfile.TemporaryDirectory(prefix="check-duckdb-", dir=Path(__file__).resolve().parent) as temporary:
        connection = duckdb.connect()
        connection.execute("SET memory_limit='4GB'")
        connection.execute(f"SET temp_directory={literal(temporary)}")
        check = Checks(connection)
        for stem in REPORTS:
            check.parquet_view(stem, [output / f"{stem}.parquet"])
        check.canonical_view("feature", feature_path)
        check_producer_ids(feature_path)
        check.parquet_view("candidates", sorted((qpx / "opendia").glob("*.candidates.parquet")))
        pg_path = qpx / "opendia.pg.parquet"
        check.canonical_view("qpx_pg", pg_path)
        check_pg_ids(pg_path)
        check.parquet_view("qpx_run", [qpx / "opendia.run.parquet"])
        require((qpx / "opendia.provenance.parquet").is_file(), "Missing QPX provenance view")
        runs = [row[0] for row in connection.execute("SELECT run_file_name FROM qpx_run").fetchall()]
        require(len(runs) == len(set(runs)), "Duplicate QPX run names")
        check.unique("feature", ["feature_id"])
        check.unique("candidates", ["feature_id"])
        check.unique("precursors", ["feature_id"])
        check.unique("precursors", ["run", "precursor_id"])
        check.unique("feature", ["run_file_name", "opendia_precursor_id", "opendia_candidate_row"])
        feature_columns = {row[0] for row in connection.execute("DESCRIBE feature").fetchall()}
        candidate_columns = [row[0] for row in connection.execute("DESCRIBE candidates").fetchall()]
        shared = ["run_file_name", "is_decoy"]
        for name in candidate_columns:
            if name.startswith("opendia_"):
                require(name in feature_columns, f"Canonical feature view is missing diagnostic extension {name}")
                shared.append(name)
        check.same("feature", "candidates", ["feature_id"], [(name, name) for name in shared], "canonical QPX versus full diagnostic candidate IDs and exact values")
        connection.execute("CREATE VIEW feature_winners AS SELECT * FROM feature WHERE opendia_winner")
        check.same("precursors", "feature_winners", ["feature_id"],
                   [("run", "run_file_name"), ("precursor_id", "opendia_precursor_id"), ("decoy", "is_decoy"),
                    ("area", "opendia_area"), ("rt_apex", "opendia_rt_apex"), ("rt_start", "opendia_rt_start"),
                    ("rt_end", "opendia_rt_end"), ("discriminant_score", "opendia_score"), ("run_q_value", "opendia_qvalue"),
                    ("global_q_value", "global_qvalue"), ("peptide_run_q_value", "peptide_qvalue"),
                    ("peptide_global_q_value", "peptide_global_qvalue"), ("protein_group_run_q_value", "pg_qvalue"),
                    ("protein_group_global_q_value", "pg_global_qvalue")], "unfiltered retained report equals canonical QPX winners")
        peak = gate("run_q_value", "global_q_value", cuts["peakgroup"])
        peptide = gate("peptide_run_q_value", "peptide_global_q_value", cuts["peptide"])
        protein = gate("protein_group_run_q_value", "protein_group_global_q_value", cuts["protein"])
        check.zero(f"SELECT feature_id FROM precursors WHERE pass_peakgroup IS DISTINCT FROM (NOT decoy AND {peak}) "
                   f"OR pass_peptide IS DISTINCT FROM (NOT decoy AND {peak} AND {peptide}) "
                   f"OR pass_protein_group IS DISTINCT FROM (NOT decoy AND {peak} AND {protein})",
                   "per-feature peak-group plus corresponding-level pass flags")
        check.zero("SELECT feature_id FROM precursors WHERE quantity_type != 'raw_ms2_area'", "raw MS2 area semantics")
        connection.execute(f"CREATE VIEW expected_precursor AS SELECT precursor_id id,run,area FROM precursors "
                           f"WHERE NOT decoy AND {peak} AND area IS NOT NULL")
        check.matrix("precursor_matrix", "precursor_id", "expected_precursor", runs,
                     ["precursor_id", "sequence", "unmodified_sequence", "charge", "proteins", "genes"])
        for stem, table, id_column, level in (("peptide", "peptides", "unmodified_sequence", "peptide"),
                                              ("protein_group", "protein_groups", "protein_group", "protein")):
            eligible = gate("run_q_value", "global_q_value", cuts[level])
            key = identifier(id_column) if stem == "peptide" else f"coalesce(nullif(proteins,''),{identifier(id_column)})"
            connection.execute(f"CREATE VIEW expected_{stem} AS SELECT {key} id,run,area FROM {table} "
                               f"WHERE NOT decoy AND {eligible} AND area IS NOT NULL")
            check.unique(table, [id_column, "run", "decoy"])
            check.matrix(stem + "_matrix", id_column, "expected_" + stem, runs, [id_column, "proteins"])
            check.zero(f"SELECT {identifier(id_column)} FROM {table} WHERE pass IS DISTINCT FROM (NOT decoy AND {eligible})",
                       f"{table}: independent level pass flags")
        check.zero(f"SELECT feature_id FROM precursors WHERE NOT decoy AND {peak} AND unmodified_sequence IS NULL",
                   "eligible target peptide metadata available for independent sum validation")
        connection.execute(f"CREATE VIEW peptide_sums AS SELECT unmodified_sequence id,run,sum(area) area "
                           f"FROM precursors WHERE NOT decoy AND {peak} AND area IS NOT NULL GROUP BY ALL")
        connection.execute("CREATE VIEW target_peptide_areas AS SELECT unmodified_sequence id,run,area FROM peptides WHERE NOT decoy AND area IS NOT NULL")
        check.same("peptide_sums", "target_peptide_areas", ["id", "run"], [("area", "area")],
                   "peptide raw areas sum eligible precursor winners exactly once", rtol=2e-14)
        connection.execute(f"CREATE VIEW protein_sums AS SELECT array_to_string(list_transform(f.pg_accessions,lambda x: x.accession),';') id, "
                           f"p.run,sum(p.area) area FROM precursors p JOIN feature f USING(feature_id) "
                           f"WHERE NOT p.decoy AND {gate('p.run_q_value','p.global_q_value',cuts['peakgroup'])} "
                           "AND p.area IS NOT NULL AND len(f.pg_accessions)>0 GROUP BY ALL")
        connection.execute("CREATE VIEW target_group_areas AS SELECT coalesce(nullif(proteins,''),protein_group) id,run,area "
                           "FROM protein_groups WHERE NOT decoy AND area IS NOT NULL")
        check.same("protein_sums", "target_group_areas", ["id", "run"], [("area", "area")],
                   "resolved protein areas sum eligible unique-group precursor evidence independently of peptide q", rtol=2e-14)
        pg_metadata = pq.read_schema(pg_path).metadata or {}
        decoy_prefix = pg_metadata.get(b"opendia_decoy_accession_prefix", b"").decode()
        accessions = "pg_accessions"
        if decoy_prefix:
            check.zero(f"SELECT pg_id FROM qpx_pg WHERE is_decoy AND "
                       f"len(list_filter(pg_accessions,lambda x: NOT starts_with(x,{literal(decoy_prefix)})))>0",
                       "canonical decoy group accessions carry the declared prefix")
            accessions = (f"CASE WHEN is_decoy THEN list_transform(pg_accessions,lambda x: "
                          f"substring(x,{len(decoy_prefix) + 1})) ELSE pg_accessions END")
        connection.execute(f"CREATE VIEW pg_long AS SELECT array_to_string({accessions},';') id,grouped_runs[1] run,is_decoy decoy, "
                           "opendia_area area,pg_qvalue run_q,global_qvalue global_q,opendia_score score FROM qpx_pg")
        connection.execute("CREATE VIEW group_long AS SELECT coalesce(nullif(proteins,''),protein_group) id,run,decoy,area, "
                           "run_q_value run_q,global_q_value global_q,discriminant_score score FROM protein_groups")
        check.unique("pg_long", ["id", "run", "decoy"])
        check.same("pg_long", "group_long", ["id", "run", "decoy"], [(name,name) for name in ("area","run_q","global_q","score")],
                   "canonical QPX protein-group evidence matches unfiltered entity long table")
        check_legacy(check, output, args.refdir.resolve() if args.refdir else None, args.legacy_rtol, args.legacy_atol)
        check_xipm(check, output, args.require_xipm)
        print(json.dumps({"qpx_files": {str(path.relative_to(output)): path.stat().st_size for path in sorted(qpx.rglob("*.parquet"))},
                          "candidate_rows": check.count("SELECT count(*) FROM candidates"),
                          "retained_winners": check.count("SELECT count(*) FROM precursors"),
                          "retained_decoys": check.count("SELECT count(*) FROM precursors WHERE decoy"),
                          "peptide_entity_runs": check.count("SELECT count(*) FROM peptides"),
                          "protein_group_entity_runs": check.count("SELECT count(*) FROM protein_groups")}, sort_keys=True))
        connection.close()
    print("ACCEPTANCE_PASS: numeric evidence, independent matrices, nulls and file integrity")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output", type=Path)
    parser.add_argument("--refdir", type=Path)
    for level in ("peakgroup", "peptide", "protein"):
        parser.add_argument("--" + level + "-fdr", type=float)
    parser.add_argument("--legacy-rtol", type=float, default=0.0)
    parser.add_argument("--legacy-atol", type=float, default=0.0)
    parser.add_argument("--require-xipm", action="store_true")
    args = parser.parse_args()
    require(all(math.isfinite(value) and value >= 0 for value in (args.legacy_rtol, args.legacy_atol)),
            "Legacy tolerances must be finite and nonnegative")
    try:
        run(args)
    except Exception as error:
        print("ACCEPTANCE_FAIL:", error)
        raise


if __name__ == "__main__":
    main()
