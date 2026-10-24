#!/usr/bin/env python3
"""Read completed outputs without modifying them.

Schema 1: provenance identifies dataset/job, artifacts contain final apparent bytes,
SQLite sidecars and row counts. Each count has status and rows (null when absent
or unavailable, never a fabricated zero). Fanout includes FEATUREs with no
FEATURE_TRANSITION rows. Native Parquet counts use footer metadata, fanout reads
only feature_id columns. Resources/rescore are supplied measurements, not runs.
All ZIP members are streamed to anonymous controlled temporary files, never
extracted using archive paths. Measure quiescent outputs for consistent results.
The canonical absolute job identity is not automatically remapped after package
relocation. Retain captured reports as evidence, or explicitly remeasure/rebind
the relocated job before rerendering. Foreign resource summaries are rejected.
"""

import argparse
from collections import Counter
from contextlib import closing
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sqlite3
import tempfile
import zipfile


EXPECTED_TABLES = ("RUN", "PRECURSOR", "PEPTIDE", "PROTEIN", "TRANSITION",
                   "FEATURE", "FEATURE_MS1", "FEATURE_MS2", "FEATURE_TRANSITION")


def count_result(rows=None, status="absent"):
    return {"status": status, "rows": rows}


def fanout_summary(values):
    values = list(values)
    return {"status": "ok", "feature_count": len(values),
            "zero_transition_features": values.count(0),
            "min": min(values) if values else None,
            "max": max(values) if values else None,
            "mean": sum(values) / len(values) if values else None}


def quote_identifier(name):
    return '"' + name.replace('"', '""') + '"'


def measure_sqlite(path, temp_dir=None):
    # A read-only SQLite connection can still create/update WAL shared-memory
    # sidecars. Read a controlled snapshot so even those input bytes stay intact.
    # Outputs must be quiescent while their database and WAL are copied.
    with tempfile.TemporaryDirectory(prefix="opendia-measure-", dir=temp_dir) as directory:
        snapshot = Path(directory) / "snapshot.sqlite"
        shutil.copyfile(path, snapshot)
        for suffix in ("-wal", "-shm", "-journal"):
            sidecar = Path(str(path) + suffix)
            if sidecar.is_file():
                shutil.copyfile(sidecar, Path(str(snapshot) + suffix))
        return read_sqlite_snapshot(snapshot)


def read_sqlite_snapshot(path):
    # Do not use immutable=1, which would hide uncheckpointed WAL rows.
    with closing(sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)) as db:
        db.execute("PRAGMA query_only=ON")
        db.execute("BEGIN")
        names = [row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table'")]
        tables = {name: count_result() for name in EXPECTED_TABLES}
        for name in names:
            rows = db.execute(f"SELECT COUNT(*) FROM {quote_identifier(name)}").fetchone()[0]
            tables[name.upper()] = count_result(rows, "ok")
        fanout = {"status": "absent", "reason": "FEATURE or FEATURE_TRANSITION absent"}
        if {"FEATURE", "FEATURE_TRANSITION"}.issubset(name.upper() for name in names):
            rows = db.execute("SELECT COUNT(ft.FEATURE_ID) FROM FEATURE f LEFT JOIN "
                              "FEATURE_TRANSITION ft ON ft.FEATURE_ID=f.ID GROUP BY f.ID")
            fanout = fanout_summary(row[0] for row in rows)
    return {"tables": tables, "feature_transition_fanout": fanout}


def measure_archive(path, temp_dir=None):
    tables = {name: count_result() for name in ("FEATURE", "FEATURE_TRANSITION")}
    members = []
    try:
        import pyarrow.parquet as pq
    except ImportError:
        pq = None
    features, transitions = set(), Counter()
    with zipfile.ZipFile(path) as archive:
        for member in archive.infolist():
            if member.is_dir() or not member.filename.endswith(".parquet"):
                continue
            result = {"member": member.filename, "uncompressed_bytes": member.file_size,
                      **count_result(status="unavailable" if pq is None else "ok")}
            name = Path(member.filename).name
            logical = {"features.parquet": "FEATURE", "feature_transition.parquet": "FEATURE_TRANSITION"}.get(name)
            if pq is not None:
                # No extract/extractall: even hostile ../ and absolute names stay
                # inside an anonymous tempfile. Copying streams bounds RAM use.
                with tempfile.TemporaryFile(dir=temp_dir) as temporary:
                    with archive.open(member) as source:
                        shutil.copyfileobj(source, temporary)
                    temporary.seek(0)
                    parquet = pq.ParquetFile(temporary)
                    result["rows"] = parquet.metadata.num_rows
                    if logical:
                        counts = tables[logical]
                        counts["rows"] = (counts["rows"] or 0) + result["rows"]
                        counts["status"] = "ok"
                        partition = str(Path(member.filename).parent)
                        for batch in parquet.iter_batches(columns=["feature_id"]):
                            ids = batch.column(0).to_pylist()
                            if logical == "FEATURE":
                                features.update((partition, value) for value in ids)
                            else:
                                transitions.update((partition, value) for value in ids)
            elif logical:
                tables[logical] = count_result(status="unavailable")
            members.append(result)
    fanout = {"status": "unavailable" if pq is None else "absent"}
    if all(value["status"] == "ok" for value in tables.values()):
        fanout = fanout_summary(transitions[feature] for feature in features)
    return {"parquet_members": members, "tables": tables,
            "feature_transition_fanout": fanout,
            "row_count_method": "pyarrow.parquet.ParquetFile.metadata.num_rows" if pq else None,
            "note": "FEATURE rows are native features.parquet rows. Other denormalized columns are not separate SQLite table counts."}


def measure_artifact(path, temp_dir=None):
    path = Path(path).resolve()
    if not path.is_file():
        return {"path": str(path), "status": "absent", "final_bytes": None, "sidecars": []}
    result = {"path": str(path), "status": "ok", "final_bytes": path.stat().st_size,
              "sidecars": []}
    if path.suffix.lower() == ".oswpq":
        result["format"] = "oswpq"
        result.update(measure_archive(path, temp_dir))
    else:
        result["format"] = "sqlite"
        result.update(measure_sqlite(path, temp_dir))
        for suffix in ("-wal", "-shm", "-journal"):
            sidecar = Path(str(path) + suffix)
            present = sidecar.is_file()
            result["sidecars"].append({"path": str(sidecar), "status": "present" if present else "absent",
                                       "final_bytes": sidecar.stat().st_size if present else None})
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--job", required=True, type=Path, help="Existing job directory, stored as canonical absolute provenance")
    parser.add_argument("--temp-dir", required=True, type=Path,
                        help="Existing controlled scratch directory with space for database/sidecar snapshots and uncompressed archive members")
    parser.add_argument("--resources-json", type=Path, help="Existing monitor summary, copied verbatim")
    parser.add_argument("--rescore-json", type=Path, help="Existing standalone rescore measurement, copied verbatim")
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args(argv)
    if not args.job.is_dir() or not args.temp_dir.is_dir():
        parser.error("--job and --temp-dir must be existing directories")
    # Refuse to overwrite any measured input or its SQLite sidecars.
    protected = {Path(str(path.resolve()) + suffix) for path in args.inputs
                 for suffix in ("", "-wal", "-shm", "-journal")}
    protected.update(path.resolve() for path in (args.resources_json, args.rescore_json) if path)
    if args.output_json.resolve() in protected or (args.output_json.exists() and any(
            path.exists() and args.output_json.samefile(path) for path in protected)):
        parser.error("output JSON must not overwrite a measurement input")
    result = {"schema_version": 1, "measured_at_utc": datetime.now(timezone.utc).isoformat(),
              "provenance": {"dataset": args.dataset, "job": str(args.job.resolve())},
              "temp_directory": str(args.temp_dir.resolve()),
              "artifacts": [measure_artifact(path, args.temp_dir) for path in args.inputs],
              "resources": json.loads(args.resources_json.read_text()) if args.resources_json else None,
              "standalone_rescore": json.loads(args.rescore_json.read_text()) if args.rescore_json else None,
              "stages": {"status": "unavailable", "reason": "Clean merged engine has no stage profiler"},
              "interpretation": "Storage row counts, not biological performance. Peptide pseudo-protein groups are not biological proteins."}
    args.output_json.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
