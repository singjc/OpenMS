#!/usr/bin/env python3
"""v4 evaluation (P', Pe', H'; V1+pool). Runs the frozen eval_fasta.py on an arm mapping, then adds the
AMEND_SURVEY_POOL estimator FDP_c(r_lib) with r_lib from Pe''s final library.

Usage: eval_v4.py <v4_slim_dir> <v2_slim_dir> <ratios_v4.json> <ratios_v2.json> <workdir> <out.json>
  P <- Pp, Pe <- Pep, H <- Hp (v4); F, Fe <- v2 (references; the survey options do not apply at K = 0).
"""
import json, os, subprocess, sys

v4, v2, rat4, rat2, work, out = sys.argv[1:7]
here = os.path.dirname(os.path.abspath(__file__))
runs = os.path.join(work, "runs"); cache = os.path.join(work, "cache")
os.makedirs(runs); os.makedirs(cache)
for arm, (src, name) in {"P": (v4, "Pp"), "Pe": (v4, "Pep"), "H": (v4, "Hp"), "F": (v2, "F"), "Fe": (v2, "Fe")}.items():
    if os.path.isdir(f"{src}/runs_v4/{name}") or os.path.isdir(f"{src}/runs_v2/{name}"):
        d = "runs_v4" if src == v4 else "runs_v2"
        os.symlink(os.path.abspath(f"{src}/{d}/{name}"), f"{runs}/{arm}")
        os.symlink(os.path.abspath(f"{src}/{d}/{name}.time"), f"{runs}/{arm}.time")
r4, r2 = json.load(open(rat4)), json.load(open(rat2))
# Wiring fix (written after the ratios, before any FDP was computed): under V1 the survey library drops pairs without
# predicted survey fragments, so Pe''s candidate set differs slightly from Fe's (r0 0.5902 vs 0.5894) and the frozen
# assert candidates_Fe == candidates_Pe would fire. Pe' uses its own candidates for r0 (eval_fasta's rule); Fe is
# re-evaluated below with its own r0 (v2), which reproduces FDP_c(Fe) = 0.0163.
json.dump({"candidates_Pe": r4["candidates_Pep"], "library_Pe": r4["library_Pep"]}, open(f"{cache}/ratios.json", "w"))
subprocess.run([sys.executable, f"{here}/eval_fasta.py", runs, cache, f"{work}/eval_frozen.json"], check=True,
               stdout=subprocess.DEVNULL)
R = json.load(open(f"{work}/eval_frozen.json"))
R["candidates_Fe"], R["final_library_ratio_Fe_sensitivity"] = r2["candidates_Fe"], r2["library_Fe"]

# AMEND_SURVEY_POOL: FDP_c(r_lib) for Pe', same procedure as eval_fasta.evaluate with k = 1 + 1/r_lib.
sys.argv = ["eval_fasta.py", runs, cache, os.devnull]
ns = {"__name__": "eval_fasta_lib"}
src = open(f"{here}/eval_fasta.py").read().split("\narms = {")[0]  # function definitions only
exec(compile(src, "eval_fasta.py", "exec"), ns)
r_lib = R["final_library_ratio_Pe_sensitivity"]["ratio"]
pe_lib = ns["evaluate"]("Pe", 1 + 1 / r_lib)
R["Pe_r_lib"] = {"r_lib": r_lib, **pe_lib}
R["Fe"] = ns["evaluate"]("Fe", 1 + 1 / r2["candidates_Fe"]["ratio"])
R["Fe_r_lib"] = {"r_lib": r2["library_Fe"]["ratio"], **ns["evaluate"]("Fe", 1 + 1 / r2["library_Fe"]["ratio"])}
fe = R["Fe"]["FDP_c"]
c3 = R["criteria"]["3_ids"]
c3["N_at_FDPc_Fe"] = R["Fe"]["N_at_FDPc_1pct"]
c3["ratio_Pe_Fe"] = c3["N_at_FDPc_Pe"] / c3["N_at_FDPc_Fe"]
c3["pass_N_at_FDPc"] = c3["N_at_FDPc_Pe"] >= 0.97 * c3["N_at_FDPc_Fe"]
c3["pass"] = c3["pass_N"] and c3["pass_N_at_FDPc"]
c3["N_at_FDPc_rlib_Pe"] = pe_lib["N_at_FDPc_1pct"]
c3["pass_N_at_FDPc_rlib"] = pe_lib["N_at_FDPc_1pct"] >= 0.97 * c3["N_at_FDPc_Fe"]
R["criteria"]["2_fdr"] = {"FDP_c_Pe": R["Pe"]["FDP_c"], "FDP_c_Fe": fe, "pass": R["Pe"]["FDP_c"] <= fe + 0.003}
R["criteria"]["2_fdr_amended"] = {
    "FDP_c_r0_Pe": R["Pe"]["FDP_c"], "FDP_c_rlib_Pe": pe_lib["FDP_c"], "FDP_c_Fe": fe,
    "pass_r0": R["Pe"]["FDP_c"] <= fe + 0.003, "pass_rlib": pe_lib["FDP_c"] <= fe + 0.003,
    "pass": R["Pe"]["FDP_c"] <= fe + 0.003 and pe_lib["FDP_c"] <= fe + 0.003,
    "N_at_FDPc_rlib_Pe": pe_lib["N_at_FDPc_1pct"]}
if "H" in R:
    R["criteria"]["H_prime"] = {"N_Hp": R["H"]["N"], "N_H_v2": 37830, "N_FH": 52847, "ratio_Hp_FH": R["H"]["N"] / 52847}
R["mapping"] = "P=Pp, Pe=Pep, H=Hp (v4, V1+pool); F, Fe from v2"
json.dump(R, open(out, "w"), indent=1)
print(json.dumps(R["criteria"], indent=1))
