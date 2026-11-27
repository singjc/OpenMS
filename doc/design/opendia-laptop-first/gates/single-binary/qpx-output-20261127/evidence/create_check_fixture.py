#!/usr/bin/env python3
"""Synthetic self-test data for check_outputs.py, not empirical acceptance evidence."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(sys.argv[1])
CANONICAL = "--old-typed" not in sys.argv[2:]


def derive_id(values, unordered=()):
    values = list(values)
    for index in unordered:
        by_encoding = {json.dumps(item, separators=(",", ":"), sort_keys=True, ensure_ascii=True): item for item in values[index]}
        values[index] = [by_encoding[key] for key in sorted(by_encoding)]
    encoded = json.dumps(values, separators=(",", ":"), sort_keys=True, ensure_ascii=True).encode()
    return int.from_bytes(hashlib.blake2b(encoded, digest_size=8).digest(), "big", signed=True)


def canonical(table, feature=False):
    if not CANONICAL:
        return table
    arrays, fields = [], []
    for field in table.schema:
        name = field.name
        values = table[name]
        if name.startswith("opendia_") or (feature and name in ("peptide_global_qvalue", "pg_qvalue")):
            name = "x_" + name if name.startswith("opendia_") else "x_opendia_" + name
            values = pa.array([None if value is None else "true" if value is True else "false" if value is False
                               else format(value, ".17g") if isinstance(value, float) else str(value)
                               for value in values.to_pylist()], type=pa.string())
            fields.append(pa.field(name, pa.string(), nullable=True))
        else:
            fields.append(field)
        arrays.append(values)
    return pa.Table.from_arrays(arrays, schema=pa.schema(fields, metadata=table.schema.metadata))

ROOT.mkdir(parents=True, exist_ok=True)
(ROOT / "qpx/opendia").mkdir(parents=True, exist_ok=True)
(ROOT / "legacy").mkdir(exist_ok=True)
metadata = {b"opendia_peakgroup_fdr": b"0.01", b"opendia_peptide_fdr": b"0.01", b"opendia_protein_fdr": b"0.01"}


def write_report(stem, fields, rows):
    schema = pa.schema(fields)
    table = pa.Table.from_pylist(rows, schema=schema)
    pq.write_table(table, ROOT / (stem + ".parquet"), row_group_size=4096)
    with (ROOT / (stem + ".tsv")).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(schema.names)
        for row in table.to_pylist():
            writer.writerow(["NA" if row[name] is None else "true" if row[name] is True else "false" if row[name] is False else row[name] for name in schema.names])


text = lambda name: pa.field(name, pa.string())
real = lambda name: pa.field(name, pa.float64())
boolean = lambda name: pa.field(name, pa.bool_())
integer = lambda name: pa.field(name, pa.int64())
meta_fields = [text("precursor_id"), text("sequence"), text("unmodified_sequence"), integer("charge"), text("proteins"), text("genes")]
fields = [integer("feature_id"), text("run")] + meta_fields + [real(name) for name in ("rt_apex", "rt_start", "rt_end", "area")] + [text("quantity_type")] + [real(name) for name in ("discriminant_score", "run_q_value", "global_q_value")] + [boolean("decoy")] + [real(name) for name in ("peptide_run_q_value", "peptide_global_q_value", "protein_group_run_q_value", "protein_group_global_q_value")] + [boolean(name) for name in ("pass_peakgroup", "pass_peptide", "pass_protein_group")]
rows = []
for feature, run, precursor, seq, area, q, pepq, group, decoy in [
    (101, "runA", "p0", "PEPTIDE", 0.0, .001, .001, True, False),
    (102, "runB", "p0", "PEPTIDE", 7.0, .02, .001, True, False),
    (103, "runA", "p1", "OTHER", 5.123456789012345, .001, .001, False, False),
    (104, "runA", "p2", "FAILPEP", 2.0, .001, .02, True, False),
    (105, "runA", "d0", "DECOYPEP", 0.0, 0.0, 0.0, True, True),
]:
    peak = not decoy and q <= .01
    rows.append(dict(feature_id=derive_id([run, precursor, 10 + len(rows)]), run=run, precursor_id=precursor, sequence=seq, unmodified_sequence=seq,
                     charge=2, proteins="P1;P2" if group else "P1;P3", genes='GENE\t"quoted"\nline' if precursor == "p0" else None,
                     rt_apex=11.123456789012345, rt_start=10.0, rt_end=13.0, area=area, quantity_type="raw_ms2_area",
                     discriminant_score=float(feature), run_q_value=q, global_q_value=.001, decoy=decoy,
                     peptide_run_q_value=pepq, peptide_global_q_value=pepq, protein_group_run_q_value=.001 if group else None,
                     protein_group_global_q_value=.001 if group else None, pass_peakgroup=peak,
                     pass_peptide=peak and pepq <= .01, pass_protein_group=peak and group))
write_report("precursors", fields, rows)
wide = []
for row in rows:
    if row["run"] == "runA" and not row["decoy"]:
        wide.append({**{field.name: row[field.name] for field in meta_fields}, "runA": row["area"], "runB": None})
write_report("precursor_matrix", meta_fields + [real("runA"), real("runB")], wide)
entity_fields = lambda key: [text(key), text("proteins"), text("run"), real("discriminant_score"), real("run_q_value"), real("global_q_value"), boolean("decoy"), real("area"), text("quantity_type"), boolean("pass")]
peptides = []
for row in rows:
    peptides.append(dict(unmodified_sequence=("DECOY_" if row["decoy"] else "") + row["unmodified_sequence"],
                         proteins=row["proteins"], run=row["run"], discriminant_score=row["discriminant_score"],
                         run_q_value=row["peptide_run_q_value"], global_q_value=row["peptide_global_q_value"],
                         decoy=row["decoy"], area=row["area"] if row["run"] == "runA" else None,
                         quantity_type="raw_ms2_area", **{"pass": not row["decoy"] and row["peptide_run_q_value"] <= .01}))
write_report("peptides", entity_fields("unmodified_sequence"), peptides)
write_report("peptide_matrix", [text("unmodified_sequence"), text("proteins"), real("runA"), real("runB")],
             [dict(unmodified_sequence=row["unmodified_sequence"], proteins=row["proteins"], runA=row["area"], runB=None)
              for row in peptides if row["pass"] and row["area"] is not None])
groups = [dict(protein_group=("DECOY_P1;P2" if decoy else "P1;P2"), proteins="P1;P2", run=run,
               discriminant_score=score, run_q_value=.001, global_q_value=.001, decoy=decoy, area=area,
               quantity_type="raw_ms2_area", **{"pass": not decoy})
          for run, decoy, score, area in [("runA", False, 104.0, 2.0), ("runB", False, 102.0, None), ("runA", True, 105.0, 0.0)]]
write_report("protein_groups", entity_fields("protein_group"), groups)
write_report("protein_group_matrix", [text("protein_group"), text("proteins"), real("runA"), real("runB")],
             [dict(protein_group="P1;P2", proteins="P1;P2", runA=2.0, runB=None)])
features = []
for index, row in enumerate(rows):
    features.append(dict(feature_id=row["feature_id"], run_file_name=row["run"], opendia_precursor_id=row["precursor_id"],
                         opendia_candidate_row=10 + index, opendia_rt_apex=row["rt_apex"], opendia_rt_start=row["rt_start"],
                         opendia_rt_end=row["rt_end"], opendia_area=row["area"], opendia_score=row["discriminant_score"],
                         opendia_qvalue=row["run_q_value"], opendia_winner=True, is_decoy=row["decoy"],
                         global_qvalue=row["global_q_value"], peptide_qvalue=row["peptide_run_q_value"],
                         peptide_global_qvalue=row["peptide_global_q_value"], pg_qvalue=row["protein_group_run_q_value"],
                         pg_global_qvalue=row["protein_group_global_q_value"],
                         pg_accessions=[dict(accession=p) for p in ("P1", "P2")] if row["protein_group_run_q_value"] is not None else None))
features.append({**features[0], "feature_id": derive_id(["runA", "p0", 99]), "opendia_candidate_row": 99, "opendia_winner": False, "opendia_score": -1.0})
feature_metadata = {**metadata, b"opendia_identity_composite": b"run_file_name,opendia_precursor_id,opendia_candidate_row",
                    b"identity_composite": b"run_file_name,cv_params" if CANONICAL else b"run_file_name,opendia_precursor_id,opendia_candidate_row",
                    b"identity_source": b"supplied" if CANONICAL else b"derived"}
feature_table = pa.Table.from_pylist(features).replace_schema_metadata(feature_metadata)
pq.write_table(canonical(feature_table, feature=True), ROOT / "qpx/opendia.feature.parquet")
pq.write_table(feature_table, ROOT / "qpx/opendia/runA.candidates.parquet")
pg_rows = []
for row in groups:
    accessions = [("DECOY_" if row["decoy"] and CANONICAL else "") + name for name in ("P1", "P2")]
    pg_rows.append(dict(pg_id=derive_id([accessions, [row["run"]], "LFQ"], (0, 1)), label="LFQ",
                        pg_accessions=accessions, grouped_runs=[row["run"]], is_decoy=row["decoy"],
                        opendia_decoy_class="decoy" if row["decoy"] else "target", opendia_area=row["area"],
                        pg_qvalue=row["run_q_value"], global_qvalue=row["global_q_value"], opendia_score=row["discriminant_score"]))
pg_metadata = {b"opendia_decoy_accession_prefix": b"DECOY_"} if CANONICAL else {}
pq.write_table(canonical(pa.Table.from_pylist(pg_rows).replace_schema_metadata(pg_metadata)), ROOT / "qpx/opendia.pg.parquet")
pq.write_table(pa.table({"run_file_name": ["runA", "runB"]}), ROOT / "qpx/opendia.run.parquet")
pq.write_table(pa.table({"step": ["synthetic"]}), ROOT / "qpx/opendia.provenance.parquet")
legacy_fields = ["run", "candidate_row", "precursor_id", "is_decoy", "sequence", "charge", "rt_apex", "rt_start", "rt_end", "area", "transition_support", "score", "q_value"]
for run in ("runA", "runB"):
    with (ROOT / f"legacy/{run}.precursors.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=legacy_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, row in enumerate(rows):
            if row["run"] == run:
                record = {key: row[key] for key in ("run", "precursor_id", "sequence", "charge", "rt_apex", "rt_start", "rt_end", "area")}
                record.update(candidate_row=index + 10, is_decoy=int(row["decoy"]), transition_support=1,
                              score=row["discriminant_score"], q_value=row["run_q_value"])
                writer.writerow(record)
with (ROOT / "legacy/precursors.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["run", "precursor_id", "area"])
    for row in rows:
        if row["pass_peakgroup"]:
            writer.writerow([row["run"], row["precursor_id"], row["area"]])
xipm = dict(feature_id=rows[0]["feature_id"], run_file_name="runA", opendia_precursor_id="p0", opendia_candidate_row=10,
            MS_LEVEL=1, PEAKMAP_TYPE="precursor", RT_START=10.0, RT_END=13.0)
for prefix, values in {"MZ": [500., 500.01], "RT": [11., 12.], "MOBILITY": [np.nan, np.nan], "INTENSITY": [0., 3.]}.items():
    xipm[prefix + "_COMPRESSION"] = 0
    xipm[prefix + "_DATA"] = np.asarray(values, dtype="=f8").tobytes()
pq.write_table(pa.Table.from_pylist([xipm]), ROOT / "runA.xipm")
print(ROOT)
