# Diagnostic: plain TDC q = (D+1)/T (pi0 = 1) from the discriminant score vs the tool's run_q_value, for Pe' and Fe.
import csv, gzip, re, sys
import numpy as np
SPLIT = re.compile(r"[;,]")
ent = lambda p: bool([x for x in SPLIT.split(p) if x]) and all(x.startswith("ENTRAP_") for x in SPLIT.split(p) if x)
for name, path, r0, rlib in (("Pe'", "v4/slim_pe/runs_v4/Pep/precursors.slim.tsv.gz", 0.5902, 0.1592),
                             ("Fe", "v2/slim/runs_v2/Fe/precursors.slim.tsv.gz", 0.5894, 0.6056)):
    s, d, e, q = [], [], [], []
    for r in csv.DictReader(gzip.open(path, "rt"), delimiter="\t"):
        s.append(float(r["discriminant_score"])); d.append(r["decoy"] == "true"); e.append(ent(r["proteins"])); q.append(float(r["run_q_value"]))
    s, d, e, q = map(np.array, (s, d, e, q))
    o = np.argsort(-s, kind="stable"); d, e, q = d[o], e[o], q[o]
    T = np.cumsum(~d); D = np.cumsum(d)
    fdr = (D + 1) / np.maximum(T, 1)
    qt = np.minimum.accumulate(fdr[::-1])[::-1]
    for lab, acc in (("tool q", q <= 0.01), ("TDC (D+1)/T", qt <= 0.01)):
        t = acc & ~d; N = int(t.sum()); ne = int((t & e).sum())
        print(f"{name:4s} {lab:12s} N={N:6d} N_e={ne:3d} decoys={int((acc & d).sum()):4d} "
              f"FDP_c(r0)={ne * (1 + 1 / r0) / N:.4f} FDP_c(r_lib)={ne * (1 + 1 / rlib) / N:.4f}")
