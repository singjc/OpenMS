#!/usr/bin/env python3
"""ADDENDUM_H4 fixed-depth check (descriptive): target-win share among the top 100/200/400 null pairs by pair score.
Usage: d4_depth.py arm slim.tsv.gz ...  (FDRBench label)"""
import csv, gzip, math, re, sys

SPLIT = re.compile(r"[;,]")
a = sys.argv[1:]
for i in range(0, len(a), 2):
    arm, path = a[i], a[i + 1]
    t, d, ent = {}, {}, {}
    with gzip.open(path, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            k = int(r["precursor_id"]) // 2
            s = float(r["discriminant_score"])
            if r["decoy"] == "true":
                d[k] = s
            else:
                t[k] = s
                ps = [p for p in SPLIT.split(r["proteins"]) if p]
                ent[k] = bool(ps) and all(p.endswith("_p_target") for p in ps)
    ks = sorted((k for k in t if k in d and ent[k]), key=lambda k: -max(t[k], d[k]))
    out = []
    for n in (100, 200, 400):
        w = sum(t[k] > d[k] for k in ks[:n]); p = w / n
        out.append(f"top{n}: {p:.3f}±{1.96 * math.sqrt(p * (1 - p) / n):.3f}")
    print(f"{arm:8s} null pairs={len(ks):8d} " + "; ".join(out))
