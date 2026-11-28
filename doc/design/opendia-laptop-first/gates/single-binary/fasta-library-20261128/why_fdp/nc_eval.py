#!/usr/bin/env python3
"""Apply PLAN_NEIGHBOUR_COMPETITION criteria E1-E4 (frozen, aac76d23cc) to nc_main.json. Usage: nc_eval.py nc_main.json"""
import json, math, sys

R = json.load(open(sys.argv[1]))
PRE = ["FbV1", "FbV1_s1", "FbV1_s2", "FbV1_g1", "FbP", "FbP_g1", "FbPp"]
CUSTOM = ["Fe", "Pe", "Pep"]
LEVELS = ["share_ge_0.333", "share_ge_0.500", "share_ge_0.667"]
PRIMARY = "share_ge_0.500"


def ratio_ci(te, de):
    n = te + de
    if not n or not de: return None
    p = te / n; h = 1.96 * math.sqrt(p * (1 - p) / n)
    lo, hi = max(p - h, 1e-9), min(p + h, 1 - 1e-9)
    return [lo / (1 - lo), hi / (1 - hi)]


def diff(a, b):
    d = a["FDP_c"] - b["FDP_c"]; h = 1.96 * math.sqrt(a["FDP_c_var"] + b["FDP_c_var"])
    return [d, d - h, d + h]


out = {}
for lev in ["before"] + LEVELS:
    g = lambda arm: R[arm][lev]
    ref = g("FbF")
    mean_pre = sum(g(a)["FDP_c"] for a in PRE) / len(PRE)
    te = sum(g(a)["N_e"] for a in PRE); de = sum(g(a)["D_e"] for a in PRE)
    out[lev] = {"FbF_FDP_c": ref["FDP_c"], "mean_prefilter_FDP_c": mean_pre, "excess": mean_pre - ref["FDP_c"],
                "FbV1_minus_FbF": diff(g("FbV1"), ref),
                "pooled_Te_De": te / de if de else None, "pooled_Te_De_ci95": ratio_ci(te, de),
                "FbF_Te_De": ref["T_e_over_D_e"],
                "arms": {a: {k: g(a)[k] for k in ("N", "N_e", "N_human", "D_e", "FDP_c", "FDP_c_ci95")}
                         for a in ["FbF"] + PRE + CUSTOM}}

b = out["before"]
verdicts = {}
for lev in LEVELS:
    a = out[lev]
    e1 = a["excess"] <= 0.5 * b["excess"] and a["FbV1_minus_FbF"][1] <= 0 <= a["FbV1_minus_FbF"][2]
    ci = a["pooled_Te_De_ci95"]
    e2 = bool(ci and a["FbF_Te_De"] is not None and ci[0] <= a["FbF_Te_De"] <= ci[1])
    e3_ratio = {x: R[x][lev]["N_human"] / R[x]["before"]["N_human"] for x in PRE}
    e3 = all(v >= 0.95 for v in e3_ratio.values())
    e4_arms = {x: {"dFDP_c": R[x][lev]["FDP_c"] - R[x]["before"]["FDP_c"],
                   "human_ratio": R[x][lev]["N_human"] / R[x]["before"]["N_human"]} for x in CUSTOM}
    e4 = all(v["dFDP_c"] <= 0.0025 and v["human_ratio"] >= 0.95 for v in e4_arms.values())
    if e1 and e3 and e4: v = "NC works"
    elif e1 and not e3: v = "removes the excess at an unacceptable ID cost"
    elif not e1: v = "rejected in this form (E1 fails)"
    else: v = "E1 and E3 pass, E4 fails (harm on custom entrapment)"
    verdicts[lev] = {"E1": e1, "E2": e2, "E3": e3, "E4": e4, "E3_human_ratio": e3_ratio, "E4": e4, "E4_arms": e4_arms,
                     "verdict": v, "primary": lev == PRIMARY}
out["criteria"] = verdicts
json.dump(out, open(sys.argv[1].replace(".json", "_eval.json"), "w"), indent=1)

f = lambda x: f"{100 * x:.2f} %"
print(f"{'level':16s} {'FbF':>8s} {'mean pre':>9s} {'excess':>8s} {'FbV1-FbF (CI)':>26s} {'pooled Te/De (CI)':>24s}")
for lev in ["before"] + LEVELS:
    a = out[lev]; d = a["FbV1_minus_FbF"]; ci = a["pooled_Te_De_ci95"]
    print(f"{lev:16s} {f(a['FbF_FDP_c']):>8s} {f(a['mean_prefilter_FDP_c']):>9s} {100 * a['excess']:+7.2f}  "
          f"{100 * d[0]:+.2f} ({100 * d[1]:+.2f}, {100 * d[2]:+.2f})      {a['pooled_Te_De']:.2f} ({ci[0]:.2f}-{ci[1]:.2f})"
          f"  FbF {a['FbF_Te_De']:.2f}")
print()
for lev in ["before"] + LEVELS:
    print(lev, " ".join(f"{x}:N={out[lev]['arms'][x]['N']},Ne={out[lev]['arms'][x]['N_e']},{f(out[lev]['arms'][x]['FDP_c'])}"
                        for x in ["FbF"] + PRE + CUSTOM))
print()
for lev, v in verdicts.items():
    print(lev, "(PRIMARY)" if v["primary"] else "", {k: v[k] for k in ("E1", "E2", "E3", "E4")}, "->", v["verdict"])
    print("   E3 human ratio:", {k: round(x, 4) for k, x in v["E3_human_ratio"].items()})
    print("   E4:", {k: {kk: round(x, 4) for kk, x in vv.items()} for k, vv in v["E4_arms"].items()})
