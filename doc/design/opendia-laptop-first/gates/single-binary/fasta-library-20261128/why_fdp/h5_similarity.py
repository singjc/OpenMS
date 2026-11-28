#!/usr/bin/env python3
"""Post hoc (descriptive, after ADDENDUM_H4): fragment similarity to the human parent, entrapment target vs its
OpenDIA decoy. For every null pair in an arm's output (or only accepted ones), the share of the target's / decoy's
singly charged b2..b(L-1) and y1..y(L-1) fragment masses (0.01 Da) that also occur in the parent's fragment set.
Usage: h5_similarity.py human_fdrbench.fasta arm slim.tsv.gz ...
"""
import csv, gzip, re, sys

MASS = {"G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841, "T": 101.04768, "C": 160.03065,
        "L": 113.08406, "I": 113.08406, "N": 114.04293, "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259,
        "M": 131.04049, "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931}
SPLIT = re.compile(r"[;,]")


def frags(s):
    b, out = 0.0, set()
    for i, a in enumerate(s[:-1]):
        b += MASS[a]
        if i >= 1: out.add(("b", round(b + 1.00728, 2)))
    y = 18.01056
    for a in reversed(s[1:]):
        y += MASS[a]
        out.add(("y", round(y + 1.00728, 2)))
    return out


prot, acc = {}, None
for line in open(sys.argv[1]):
    if line.startswith(">"):
        parts = line[1:].split()[0].split("|"); acc = parts[1] if len(parts) >= 3 else parts[0]; prot[acc] = []
    else:
        prot[acc].append(line.strip())
prot = {k: "".join(v) for k, v in prot.items()}


def parent(pep, ps):
    for p in ps:
        es, hs = prot.get(p), prot.get(p[:-9] if p.endswith("_p_target") else None)
        if es and hs and len(es) == len(hs):
            i = es.find(pep)
            if i >= 0: return hs[i:i + len(pep)]
    return None


a = sys.argv[2:]
for i in range(0, len(a), 2):
    arm, path = a[i], a[i + 1]
    m = {}
    with gzip.open(path, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            k = int(r["precursor_id"]) // 2
            m.setdefault(k, {})[r["decoy"] == "true"] = (r["unmodified_sequence"], [p for p in SPLIT.split(r["proteins"]) if p],
                                                         float(r["discriminant_score"]), float(r["run_q_value"]))
    for scope in ("all null pairs", "accepted, target won", "accepted, decoy won"):
        st = sd = n = 0
        hiT = hiD = 0
        for k, v in m.items():
            if False not in v or True not in v: continue
            t, d = v[False], v[True]
            if not (t[1] and all(p.endswith("_p_target") for p in t[1])): continue
            win_dec = d[2] >= t[2]; q = d[3] if win_dec else t[3]
            if scope != "all null pairs":
                if q > 0.01 or win_dec != (scope == "accepted, decoy won"): continue
            par = parent(t[0], t[1])
            if not par or set(MASS).isdisjoint(par): continue
            fp = frags(par); ft, fd = frags(t[0]), frags(d[0])
            a_, b_ = len(ft & fp) / len(ft), len(fd & fp) / len(fd)
            st += a_; sd += b_; n += 1; hiT += a_ >= 0.5; hiD += b_ >= 0.5
        if n:
            print(f"{arm:8s} {scope:22s} n={n:7d} mean shared-with-parent: target {st / n:.3f} decoy {sd / n:.3f}; "
                  f"share >= 50 %: target {hiT / n:.3f} decoy {hiD / n:.3f}")
