#!/usr/bin/env python3
"""Positive and deliberate-corruption tests for scratch check_outputs.py."""
import csv
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent


def write_table(directory, name, table):
    pq.write_table(table, directory / (name + ".parquet"), row_group_size=4096)
    if (directory / (name + ".tsv")).exists():
        with (directory / (name + ".tsv")).open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(table.schema.names)
            for row in table.to_pylist():
                writer.writerow(["NA" if row[key] is None else "true" if row[key] is True else "false" if row[key] is False else row[key] for key in table.schema.names])


def change(directory, name, column, index, value):
    path = directory / (name + ".parquet")
    table = pq.read_table(path)
    values = table[column].to_pylist()
    values[index] = value
    table = table.set_column(table.schema.get_field_index(column), column, pa.array(values, type=table[column].type))
    write_table(directory, name, table)


def invoke(directory, expected, message, extra=()):
    result = subprocess.run([sys.executable, str(HERE / "check_outputs.py"), str(directory), *extra], capture_output=True, text=True)
    output = result.stdout + result.stderr
    if expected:
        assert result.returncode == 0 and "ACCEPTANCE_PASS" in output, output
    else:
        assert result.returncode != 0 and "ACCEPTANCE_FAIL" in output and message in output, output
    print("SELFTEST_PASS", directory.name, message or "all gates")


with tempfile.TemporaryDirectory(prefix="check-selftest-", dir=HERE) as workspace:
    workspace = Path(workspace)
    baseline = workspace / "baseline"
    subprocess.run([sys.executable, str(HERE / "create_check_fixture.py"), str(baseline)], check=True, capture_output=True)
    invoke(baseline, True, "", ("--require-xipm", "--refdir", str(baseline / "legacy")))
    old_typed = workspace / "old_typed"
    subprocess.run([sys.executable, str(HERE / "create_check_fixture.py"), str(old_typed), "--old-typed"], check=True, capture_output=True)
    invoke(old_typed, True, "", ("--require-xipm",))

    def variant(name):
        target = workspace / name
        shutil.copytree(baseline, target)
        return target

    case = variant("producer_hash_changed")
    for relative in ("qpx/opendia.feature.parquet", "qpx/opendia/runA.candidates.parquet"):
        path = case / relative
        table = pq.read_table(path)
        values = table["feature_id"].to_pylist()
        values[0] = 999
        table = table.set_column(table.schema.get_field_index("feature_id"), "feature_id", pa.array(values, type=pa.int64()))
        pq.write_table(table, path)
    invoke(case, False, "independent producer feature_id hash mismatch")

    case = variant("canonical_decimal_changed")
    change(case / "qpx", "opendia.feature", "x_opendia_area", 0, "0.0000000000000001")
    invoke(case, False, "canonical QPX versus full diagnostic candidate IDs and exact values")

    case = variant("canonical_bad_decimal")
    change(case / "qpx", "opendia.feature", "x_opendia_area", 0, "not-a-number")
    invoke(case, False, "Could not convert string")

    case = variant("pg_hash_changed")
    change(case / "qpx", "opendia.pg", "pg_id", 0, 999)
    invoke(case, False, "independent protein-group pg_id hash mismatch")

    case = variant("null_becomes_zero")
    change(case, "precursor_matrix", "runB", 0, 0.0)
    invoke(case, False, "precursor_matrix: exact eligible cells")

    case = variant("candidate_missing")
    path = case / "qpx/opendia/runA.candidates.parquet"
    pq.write_table(pq.read_table(path).slice(1), path)
    invoke(case, False, "canonical QPX versus full diagnostic candidate IDs")

    case = variant("peptide_double_count")
    change(case, "peptides", "area", 2, 6.123456789012345)
    change(case, "peptide_matrix", "runA", 1, 6.123456789012345)
    invoke(case, False, "peptide raw areas sum eligible precursor winners")

    case = variant("protein_crossgated_by_peptide")
    change(case, "protein_groups", "area", 0, 0.0)
    change(case, "protein_group_matrix", "runA", 0, 0.0)
    invoke(case, False, "resolved protein areas sum eligible unique-group precursor evidence independently of peptide q")

    case = variant("legacy_candidate_changed")
    path = case / "legacy/runA.precursors.tsv"
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    rows[1][1] = "999"
    with path.open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
    invoke(case, False, "retained candidate identity versus legacy candidate rows")

    case = variant("bad_xipm_lengths")
    path = case / "runA.xipm"
    table = pq.read_table(path)
    table = table.set_column(table.schema.get_field_index("RT_DATA"), "RT_DATA", pa.array([b"12345678"], type=pa.binary()))
    pq.write_table(table, path)
    invoke(case, False, "unequal XIPM point array lengths")

    case = variant("bad_xipm_link")
    path = case / "runA.xipm"
    table = pq.read_table(path)
    table = table.set_column(table.schema.get_field_index("feature_id"), "feature_id", pa.array([999], type=pa.int64()))
    pq.write_table(table, path)
    invoke(case, False, "every XIPM map links to the actual exported QPX winner")

    case = variant("optional_xipm_unavailable")
    (case / "runA.xipm").rename(case / "runA.xipm.hidden")
    invoke(case, True, "")
    invoke(case, False, "--require-xipm requested", ("--require-xipm",))

    invoke(baseline, False, "threshold contradicts stored metadata", ("--peptide-fdr", "0.02"))

print("HARNESS_SELFTEST_PASS: new canonical and old typed positive fixtures plus thirteen corruption/unavailable gates")
