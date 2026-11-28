#!/usr/bin/env python3
"""H3 read-out (PLAN_PREFILTER_FDP): from <run>.features.tsv, per precursor the best window (by the given score),
per pair the winner (ties to the decoy), q by TDC (D+1)/T with pi0 = 1 (D1 showed the tool's q equals this), then at
q <= 0.01: T, T_e, D_e, FDP_c(r0), and the top-0.1 % null target-win share. Done for score (final) and network_score
(pre-GBM). Entrapment status from the slim precursor table (protein labels per precursor_id).
Usage: h3_network.py <features.tsv> <slim.tsv.gz> <r0>
Run on codon (features tables are large).
"""
import csv, gzip, re, sys
import numpy as np

feat, slim, r0 = sys.argv[1], sys.argv[2], float(sys.argv[3])
SPLIT = re.compile(r"[;,]")
ent = {}
with gzip.open(slim, "rt") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["decoy"] == "false":
            ps = [p for p in SPLIT.split(row["proteins"]) if p]
            ent[int(row["precursor_id"]) // 2] = bool(ps) and all(p.endswith("_p_target") for p in ps)
cols = {}
best = {"score": {}, "network_score": {}}
with open(feat) as f:
    head = f.readline().rstrip("\n").split("\t")
    ip, idc, isc, ins = head.index("precursor_id"), head.index("is_decoy"), head.index("score"), head.index("network_score")
    for line in f:
        r = line.rstrip("\n").split("\t")
        p = int(r[ip])
        for name, i in (("score", isc), ("network_score", ins)):
            v = float(r[i])
            if v <= -1e30:
                continue
            b = best[name]
            if p not in b or v > b[p]:
                b[p] = v
for name, b in best.items():
    pairs = {}
    for p, v in b.items():
        pairs.setdefault(p // 2, {})[p % 2] = v
    win = []
    for k, m in pairs.items():
        t, d = m.get(0, -np.inf), m.get(1, -np.inf)
        dec = d >= t
        win.append((max(t, d), dec, k, t, d))
    win.sort(key=lambda x: (-x[0], -x[1]))
    T = D = 0
    fdr = []
    for s, dec, k, t, d in win:
        D += dec; T += (not dec)
        fdr.append((D + 1) / max(T, 1))
    q = np.minimum.accumulate(np.array(fdr)[::-1])[::-1]
    acc = [w for w, qq in zip(win, q) if qq <= 0.01]
    Tn = sum(1 for w in acc if not w[1]); Te = sum(1 for w in acc if not w[1] and ent.get(w[2], False))
    De = sum(1 for w in acc if w[1] and ent.get(w[2], False))
    nulls = [w for w in win if ent.get(w[2], False) and np.isfinite(w[3]) and np.isfinite(w[4])]
    top = nulls[: max(1, len(nulls) // 1000)]
    share = sum(1 for w in top if w[3] > w[4]) / len(top)
    print(f"{name:14s} T={Tn} T_e={Te} D_e={De} FDP_c={Te * (1 + 1 / r0) / max(Tn, 1):.4f} "
          f"top0.1%_null_target_share={share:.3f} (n={len(top)})")
