#!/usr/bin/env python3
"""PLAN_CLAIM_COMPETITION evaluation (frozen a5e180dac9; written before any Cc result). Applies K1-K3 mechanically.

Usage: eval_cc.py <slim_dir> <ratios_cc.json> <out.json>
  slim_dir/runs_cc/<ARM>/precursors.slim.tsv.gz (+ <ARM>.time); ratios_cc.json from ratio_cc.sh: candidates_<ARM>.
Per arm: N (target precursor winners at run q <= 0.01), N_e, N_human, FDP_c(r0) = N_e (1 + 1/r0) / N with each arm's
own r0, 95 % Poisson CI on N_e, wall, peak RSS. Same estimator as eval_fb.py / eval_fasta.py.
Baselines (before) are the recorded values of the matching existing arms (eval_fb.json, d1_pf.json, eval_v4.json).
"""
import csv, glob, gzip, json, math, re, sys

slim, ratios, out = sys.argv[1:4]
SPLIT = re.compile(r"[;,]")
LABEL = {"CcV1": "fdrbench", "CcV1_s1": "fdrbench", "CcP": "fdrbench", "CcPp": "fdrbench", "CcF": "fdrbench",
         "Ce": "entrap", "CcPep": "entrap"}
# Before (recorded; N, N_e, FDP_c(r0), N_human). FbV1_s1 human = T - T_e from d1_pf.json.
BEFORE = {"FbV1": (14994, 129, 0.017212, 14865), "FbV1_s1": (14774, 118, 0.015979, 14656), "FbP": (11120, 85, 0.015291, 11035),
          "FbPp": (20258, 149, 0.014715, 20109), "FbF": (8074, 49, 0.012140, 8025), "Fe": (9270, 56, 0.016291, 9214),
          "Pep": (23546, 82, 0.009383, 23464)}
MATCH = {"CcV1": "FbV1", "CcV1_s1": "FbV1_s1", "CcP": "FbP", "CcPp": "FbPp", "CcF": "FbF", "Ce": "Fe", "CcPep": "Pep"}


def member(mode):
    return (lambda p: p.endswith("_p_target")) if mode == "fdrbench" else (lambda p: p.startswith("ENTRAP_"))


def timefield(arm, key):
    for line in open(f"{slim}/runs_cc/{arm}.time"):
        if key in line: return line.rsplit(" ", 1)[1].strip()


def evaluate(arm, r0):
    is_e = member(LABEL[arm])
    rows = [r for r in csv.DictReader(gzip.open(f"{slim}/runs_cc/{arm}/precursors.slim.tsv.gz", "rt"), delimiter="\t")
            if r["decoy"] == "false" and float(r["run_q_value"]) <= 0.01]
    n = len(rows)
    ne = sum(1 for r in rows if (lambda ps: bool(ps) and all(is_e(p) for p in ps))([p for p in SPLIT.split(r["proteins"]) if p]))
    k = 1 + 1 / r0
    t = [float(x) for x in timefield(arm, "Elapsed (wall clock)").split(":")]
    return {"N": n, "N_e": ne, "N_human": n - ne, "r0": r0, "FDP_c": ne * k / n,
            "ci95": [max(ne - 1.96 * math.sqrt(ne), 0) * k / n, (ne + 1.96 * math.sqrt(ne)) * k / n],
            "var": ne * (k / n) ** 2, "wall_s": sum(v * 60 ** i for i, v in enumerate(reversed(t))),
            "maxrss_gb": int(timefield(arm, "Maximum resident")) / 1024 ** 2}


RT = json.load(open(ratios))
R = {a: evaluate(a, RT[f"candidates_{a}"]["ratio"]) for a in LABEL}
PRE_AFTER, PRE_BEFORE = ["CcV1", "CcV1_s1", "CcP", "CcPp"], ["FbV1", "FbV1_s1", "FbP", "FbPp"]
ex_before = sum(BEFORE[a][2] for a in PRE_BEFORE) / 4 - BEFORE["FbF"][2]
ex_after = sum(R[a]["FDP_c"] for a in PRE_AFTER) / 4 - R["CcF"]["FDP_c"]
d = R["CcV1"]["FDP_c"] - R["CcF"]["FDP_c"]; h = 1.96 * math.sqrt(R["CcV1"]["var"] + R["CcF"]["var"])
K1 = {"excess_before": ex_before, "excess_after": ex_after, "limit": 0.5 * ex_before,
      "CcV1_minus_CcF": [d, d - h, d + h], "pass": ex_after <= 0.5 * ex_before and d - h <= 0 <= d + h}
K2 = {a: {"after": R[a]["FDP_c"], "before": BEFORE[MATCH[a]][2], "pass": R[a]["FDP_c"] <= BEFORE[MATCH[a]][2] + 0.0025}
      for a in ("CcF", "Ce", "CcPep")}
K2["pass"] = all(v["pass"] for v in K2.values())
K3 = {a: {"human_after": R[a]["N_human"], "human_before": BEFORE[MATCH[a]][3],
          "ratio": R[a]["N_human"] / BEFORE[MATCH[a]][3]} for a in LABEL}
for v in K3.values(): v["pass"] = v["ratio"] >= 0.97
K3["pass"] = all(v["pass"] for k, v in K3.items() if k != "pass")
if K1["pass"] and K2["pass"] and K3["pass"]: verdict = "claim competition works"
elif not K1["pass"]: verdict = "fails (K1)"
elif not K3["pass"]: verdict = "costs IDs (K1 passes, K3 fails)"
else: verdict = "K1 and K3 pass, K2 fails (FDP harm)"
R["criteria"] = {"K1": K1, "K2": K2, "K3": K3, "verdict": verdict}
json.dump(R, open(out, "w"), indent=1)
for a in LABEL:
    r = R[a]; b = BEFORE[MATCH[a]]
    print(f"{a:8s} N {r['N']:6d} (before {b[0]:6d})  N_e {r['N_e']:4d} ({b[1]:4d})  FDP_c {100 * r['FDP_c']:.2f} % "
          f"({100 * r['ci95'][0]:.2f}-{100 * r['ci95'][1]:.2f}; before {100 * b[2]:.2f} %)  human x{r['N_human'] / b[3]:.3f}  "
          f"wall {r['wall_s']:.0f} s  rss {r['maxrss_gb']:.1f} GB")
print(json.dumps(R["criteria"], indent=1))
