#!/usr/bin/env python3
"""Synthetic test for nc_remedy.py: tiny precursors.tsv + library parquet with one planted neighbour."""
import os, subprocess, sys, json, tempfile
import pyarrow as pa, pyarrow.parquet as pq

d = tempfile.mkdtemp()
H = "precursor_id\tunmodified_sequence\tproteins\trt_apex\trt_start\trt_end\tdiscriminant_score\trun_q_value\tdecoy\n"
rows = []
lib = []
# 300 strong human target pairs (ids 2k target, 2k+1 decoy), well separated in m/z so they do not interact
for k in range(300):
    rows.append((2 * k, f"HUM{k}", "P1", 100 + k, 99 + k, 101 + k, 10.0 + k * 0.01, 0.0, "false"))
    rows.append((2 * k + 1, f"DEC{k}", "DECOY_P1", 100 + k, 99 + k, 101 + k, -5.0, 1.0, "true"))
    for r in range(1, 13):
        lib.append((str(2 * k), 400.0 + k * 100, 1.0, 1000.0 + k * 100 + r, r))
        lib.append((str(2 * k + 1), 400.0 + k * 100, 1.0, 5000.0 + k * 100 + r, r))
# planted entrapment target 1000 sharing 4 of 6 fragments with HUM0, co-eluting, lower score but accepted
rows.append((1000, "ENT", "X_p_target", 100.2, 100.0, 100.4, 9.6, 0.0, "false"))
rows.append((1001, "TNE", "DECOY_X_p_target", 100.2, 100.0, 100.4, -6.0, 1.0, "true"))
for r in range(1, 13):
    lib.append(("1000", 401.0, 1.02, 1000.0 + r if r <= 4 else 9000.0 + r, r))
    lib.append(("1001", 401.0, 1.02, 9500.0 + r, r))
# a few decoy winners so q is not trivially 0
for k in range(3):
    rows.append((2000 + 2 * k, f"W{k}", "P9", 500, 499, 501, -9.0, 1.0, "false"))
    rows.append((2001 + 2 * k, f"WD{k}", "DECOY_P9", 500, 499, 501, 9.5, 0.0, "true"))
open(f"{d}/p.tsv", "w").write(H + "".join("\t".join(map(str, r)) + "\n" for r in rows))
t = pa.table({"precursor_id": [x[0] for x in lib], "precursor_mz": [x[1] for x in lib], "precursor_im": [x[2] for x in lib],
              "product_mz": [x[3] for x in lib], "fragment_rank": pa.array([x[4] for x in lib], pa.int32())})
pq.write_table(t, f"{d}/l.parquet")
here = os.path.dirname(os.path.abspath(__file__))
subprocess.run([sys.executable, f"{here}/nc_remedy.py", f"{d}/o.json", f"fdrbench:t:{d}/p.tsv:{d}/l.parquet:1.0"], check=True)
R = json.load(open(f"{d}/o.json"))["t"]
b, a5, a67 = R["before"], R["share_ge_0.500"], R["share_ge_0.667"]
assert b["N_e"] == 1 and b["N"] == 301, b
assert a5["N_e"] == 0 and a5["N"] == 300 and a5["removed_accepted_targets_null"] == 1, a5   # 4/6 >= 0.5: removed
assert a67["N_e"] == 0, a67                                                                 # 4/6 >= 0.667: removed
assert R["share_ge_0.333"]["removed_q05_decoys"] == 0
print("nc_remedy synthetic test OK")
