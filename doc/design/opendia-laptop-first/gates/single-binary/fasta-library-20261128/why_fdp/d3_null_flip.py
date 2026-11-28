#!/usr/bin/env python3
"""D3 (diagnostic, existing outputs only): target-decoy exchangeability across ALL entrapment (null) pairs.

For every entrapment pair with both members scored, does the target beat its own decoy? Under exchangeable decoys the
share is 0.5 at every score level. Binned by the pair's max score (deciles of entrapment pairs, plus the top 1 % and
top 0.1 %). Human pairs shown for contrast (they contain true targets, so > 0.5 is expected there).
Usage: d3_null_flip.py label arm slim.tsv.gz ...
"""
import csv, gzip, math, re, sys

SPLIT = re.compile(r"[;,]")
args = sys.argv[1:]
for i in range(0, len(args), 3):
    mode, arm, path = args[i:i + 3]
    is_e = (lambda p: p.startswith("ENTRAP_")) if mode == "entrap" else (lambda p: p.endswith("_p_target"))
    t, d, ent = {}, {}, {}
    with gzip.open(path, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            k = int(row["precursor_id"]) // 2
            s = float(row["discriminant_score"])
            if row["decoy"] == "true":
                d[k] = s
            else:
                t[k] = s
                ps = [p for p in SPLIT.split(row["proteins"]) if p]
                ent[k] = bool(ps) and all(is_e(p) for p in ps)
    both = [k for k in t if k in d]
    for name, sel in (("entrapment", True), ("human", False)):
        ks = sorted((k for k in both if ent[k] == sel), key=lambda k: -max(t[k], d[k]))
        n = len(ks)
        out = []
        for lab, lo, hi in (("top0.1%", 0, n // 1000), ("top1%", 0, n // 100), ("top10%", 0, n // 10), ("all", 0, n),
                            ("10-50%", n // 10, n // 2), ("50-100%", n // 2, n)):
            seg = ks[lo:hi]
            w = sum(t[k] > d[k] for k in seg); m = len(seg)
            p = w / m if m else float("nan"); h = 1.96 * math.sqrt(p * (1 - p) / m) if m else 0
            out.append(f"{lab}: {p:.3f}±{h:.3f} (n={m})")
        print(f"{arm:5s} {name:10s} pairs={n:8d} target-wins share: " + "; ".join(out))
