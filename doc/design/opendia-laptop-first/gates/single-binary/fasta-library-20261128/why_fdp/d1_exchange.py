#!/usr/bin/env python3
"""D1 (diagnostic, existing outputs only): within-null target-decoy exchangeability per arm.

For each pair (precursor_id // 2) the pair winner is the member with the higher discriminant_score (both members of a
pair share a fold, so the per-fold standardisation keeps their order). Accepted = pair winners with run_q_value <= q.
Pair type: entrapment if every protein of the pair (DECOY_ stripped) carries the entrapment label, else human.

Reported at q = 0.01 (and 0.005, 0.02):
  T_e, D_e  entrapment-pair target / decoy winners accepted. Under valid TDC (null pairs exchangeable) T_e ~ D_e.
  T_h, D_h  human-pair target / decoy winners accepted.
  D = D_e + D_h, the decoys TDC counts. FDR_tdc = (D + 1) / T.
  FDP_c = T_e (1 + 1/r0) / T as before; entrapment-implied false human targets T_e / r0, compared with D_h.
  losers_q: decoy rows with q <= 0.01 that lost their pair (sanity: should be 0 if q is assigned to pair winners only).
Usage: d1_exchange.py OUT.json label:arm:path:r0 ...   (label = entrap | fdrbench)
"""
import csv, gzip, json, math, re, sys

SPLIT = re.compile(r"[;,]")
out = sys.argv[1]
R = {}
for spec in sys.argv[2:]:
    mode, arm, path, r0 = spec.split(":", 3)
    r0 = float(r0)
    is_e = (lambda p: p.startswith("ENTRAP_")) if mode == "entrap" else (lambda p: p.endswith("_p_target"))
    pairs = {}
    with gzip.open(path, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            k = int(row["precursor_id"]) // 2
            ps = [p[6:] if p.startswith("DECOY_") else p for p in SPLIT.split(row["proteins"]) if p]
            ent = bool(ps) and all(is_e(p) for p in ps)
            pairs.setdefault(k, []).append((float(row["discriminant_score"]), row["decoy"] == "true", float(row["run_q_value"]), ent))
    res = {}
    losers = 0
    for q in (0.005, 0.01, 0.02):
        c = {"T_e": 0, "D_e": 0, "T_h": 0, "D_h": 0}
        for k, m in pairs.items():
            w = max(m, key=lambda x: (x[0], x[1]))  # ties to the decoy, as pooledPairedTDC
            if q == 0.01:
                losers += sum(1 for x in m if x is not w and x[1] and x[2] <= 0.01)
            if w[2] > q:
                continue
            ent = any(x[3] for x in m) if not w[1] else any(x[3] for x in m)
            c[("D_" if w[1] else "T_") + ("e" if ent else "h")] += 1
        T, D = c["T_e"] + c["T_h"], c["D_e"] + c["D_h"]
        ratio = c["T_e"] / c["D_e"] if c["D_e"] else None
        # 95 % CI of T_e / D_e from a binomial on T_e / (T_e + D_e)
        n = c["T_e"] + c["D_e"]
        if n:
            p = c["T_e"] / n; h = 1.96 * math.sqrt(p * (1 - p) / n)
            lo, hi = max(p - h, 1e-9), min(p + h, 1 - 1e-9)
            ci = [lo / (1 - lo), hi / (1 - hi)]
        else:
            ci = None
        c.update(T=T, D=D, FDR_tdc=(D + 1) / T if T else None, FDP_c=c["T_e"] * (1 + 1 / r0) / T if T else None,
                 Te_over_De=ratio, Te_over_De_ci95=ci, implied_false_human=c["T_e"] / r0,
                 Dh_over_implied=c["D_h"] / (c["T_e"] / r0) if c["T_e"] else None)
        res[str(q)] = c
    res["decoy_losers_with_q_le_0.01"] = losers
    R[arm] = res
    c = res["0.01"]
    print(f"{arm:5s} T={c['T']:6d} T_e={c['T_e']:4d} D_e={c['D_e']:4d} T_e/D_e={c['Te_over_De'] and round(c['Te_over_De'],2)} "
          f"ci={c['Te_over_De_ci95'] and [round(x,2) for x in c['Te_over_De_ci95']]} D_h={c['D_h']:4d} "
          f"impl_false_h={c['implied_false_human']:.0f} FDRtdc={c['FDR_tdc']:.4f} FDP_c={c['FDP_c']:.4f} losers={losers}")
json.dump(R, open(out, "w"), indent=1)
