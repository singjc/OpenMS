#!/usr/bin/env python3
"""S1 (diagnostic, existing outputs only): is the target advantage among null pairs a property of the SELECTED pairs
(present even under K = 0 scoring) or created by SCORING in the prefiltered library?

Pair indices are the same across FDRBench arms (same FASTA and enumeration). Using the K = 0 arm's scores and q (FbF),
count accepted (q <= 0.01) entrapment pair winners, split by whether the pair is in a prefilter arm's library.
Under H2 (selection), T_e/D_e inside the kept subset is high under FbF scoring too.
Usage: s1_selection.py FbF_slim.tsv.gz kept_FbV1.txt [kept_FbP.txt ...]
kept_*.txt: one pair index per line (pairs in the arm's final library), written on codon by s1_kept.py.
"""
import csv, gzip, math, re, sys

SPLIT = re.compile(r"[;,]")
is_e = lambda p: p.endswith("_p_target")
best, ent = {}, {}
with gzip.open(sys.argv[1], "rt") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        k = int(row["precursor_id"]) // 2
        s, dec, q = float(row["discriminant_score"]), row["decoy"] == "true", float(row["run_q_value"])
        ps = [p[6:] if p.startswith("DECOY_") else p for p in SPLIT.split(row["proteins"]) if p]
        if not dec:
            ent[k] = bool(ps) and all(is_e(p) for p in ps)
        if k not in best or (s, dec) > best[k][:2]:
            best[k] = (s, dec, q)


def ci(t, d):
    n = t + d
    if not n: return None
    p = t / n; h = 1.96 * math.sqrt(p * (1 - p) / n)
    lo, hi = max(p - h, 1e-9), min(p + h, 1 - 1e-9)
    return [round(lo / (1 - lo), 2), round(hi / (1 - hi), 2)]


for path in sys.argv[2:]:
    kept = set(int(x) for x in open(path))
    c = {True: [0, 0, 0, 0], False: [0, 0, 0, 0]}  # in kept: [T_e, D_e, T_h, D_h]
    for k, (s, dec, q) in best.items():
        if q > 0.01 or k not in ent:
            continue
        c[k in kept][(0 if ent[k] else 2) + (1 if dec else 0)] += 1
    for side in (True, False):
        te, de, th, dh = c[side]
        print(f"{path}: {'in kept' if side else 'not kept'}: T_e={te} D_e={de} T_e/D_e={te / de if de else None} "
              f"ci={ci(te, de)} T_h={th} D_h={dh}")
