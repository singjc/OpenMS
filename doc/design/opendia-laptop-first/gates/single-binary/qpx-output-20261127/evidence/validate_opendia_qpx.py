#!/usr/bin/env python3
"""Run the unmodified, pinned upstream QPX dataset conformance gate.

Usage: qpx-conformance-env/bin/python validate_opendia_qpx.py OUTPUT/qpx
A nonzero exit means conformance failed, including warnings. No issues are filtered.
The separately reported producer-ID check covers the requested custom candidate
composite even when this upstream pin rejects non-YAML identity fields.
"""
import argparse
import json
import subprocess
from pathlib import Path

import pyarrow.parquet as pq
import qpx
from qpx.core.data.identity import derive_composite_id

PIN = "2dadc0e120cf8b868aab0c4b6489f047efc19561"
SOURCE = Path("/home/singjc/.jcode/scratch/qpx-output-20261127/qpx-spec")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("--json-output", type=Path)
    args = parser.parse_args()
    revision = subprocess.check_output(["git", "-C", str(SOURCE), "rev-parse", "HEAD"], text=True).strip()
    if revision != PIN:
        raise RuntimeError(f"Upstream checkout moved: expected {PIN}, got {revision}")
    report = {"upstream_revision": revision, "dataset": str(args.dataset.resolve()), "views": {}}
    with qpx.open_dataset(args.dataset) as dataset:
        results = dataset.validate(strict=True, semantic=True, check_identity=True)
    passed = bool(results)
    for name, result in results.items():
        print(result.summary)
        for issue in result.issues:
            print(f"  {issue.severity}: {issue.check} [{issue.column}] {issue.message}")
        report["views"][name] = {"is_valid": result.is_valid, "issues": [vars(issue) for issue in result.issues]}
        passed = passed and result.is_valid and not result.warnings

    # This supplements, never replaces or suppresses, the upstream check.
    checked = stale = 0
    for file in sorted(args.dataset.glob("*.feature.parquet")):
        parquet = pq.ParquetFile(file)
        metadata = parquet.schema_arrow.metadata or {}
        composite = tuple(metadata.get(b"opendia_identity_composite", metadata[b"identity_composite"]).decode().split(","))
        physical = {name: name if name in parquet.schema_arrow.names else "x_" + name for name in composite}
        for batch in parquet.iter_batches(batch_size=8192, columns=["feature_id", *physical.values()]):
            for raw in batch.to_pylist():
                row = {name: raw[column] for name, column in physical.items()}
                if isinstance(row.get("opendia_candidate_row"), str):
                    row["opendia_candidate_row"] = int(row["opendia_candidate_row"])
                checked += 1
                stale += raw["feature_id"] != derive_composite_id(row, composite)
    report["producer_identity"] = {"rows": checked, "stale": stale}
    report["conformance_passed"] = passed and stale == 0
    print(f"Producer candidate identity: {checked} rows checked, {stale} stale")
    print("Overall: " + ("VALID" if report["conformance_passed"] else "INVALID"))
    if args.json_output:
        args.json_output.write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["conformance_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
