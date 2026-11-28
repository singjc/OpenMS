#!/usr/bin/env python3
"""PLAN_FASTA evaluation (frozen criteria 2-4 plus the reported goal, H, seed spread and library ratio).

Usage: eval_fasta.py <runs_dir> <cache_dir> <out.json>
  runs_dir/<ARM>/precursors.tsv   top-level precursor-winner table (targets and decoys, run_q_value)
  runs_dir/<ARM>.time             /usr/bin/time -v
  cache_dir/<ARM>/*.candidates.parquet, *.library.parquet   (Pe, Fe) for r0 and the final-library ratio

N = target precursor winners with run_q_value <= 0.01. Entrapment = every protein starts with ENTRAP_.
FDP_c = N_e (1 + 1/r0) / N, where r0 = entrapment/human target precursors in the candidate set (before the
prefilter). N@FDP_c 1 % = the largest prefix of target winners in the tool's q order (ties keep the
discriminant-score order) with FDP_c <= 0.01. The procedure is the same as pgent/eval_pgent.py.
"""
import csv, glob, json, re, sys
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

runs, cache, out = sys.argv[1:4]
SPLIT = re.compile(r"[;,]")


def is_entrap(proteins):
    ps = [p for p in SPLIT.split(proteins) if p]
    return bool(ps) and all(p.startswith("ENTRAP_") for p in ps)


def wall(arm):
    for line in open(f"{runs}/{arm}.time"):
        if "Elapsed (wall clock)" in line:
            t = [float(x) for x in line.rsplit(" ", 1)[1].strip().split(":")]
            return sum(v * 60 ** i for i, v in enumerate(reversed(t)))


def maxrss_gb(arm):
    for line in open(f"{runs}/{arm}.time"):
        if "Maximum resident" in line:
            return int(line.split()[-1]) / 1024 ** 2


def library_ratio(path):
    """Target precursors that are entrapment vs human. Fragment rows of a precursor are contiguous."""
    t = pq.read_table(path, columns=["precursor_index", "decoy", "proteins"]).combine_chunks()
    idx = t.column("precursor_index").to_numpy()
    first = np.ones(len(idx), bool)
    first[1:] = idx[1:] != idx[:-1]
    assert len(np.unique(idx)) == first.sum(), "precursor rows are not contiguous"
    keep = first & ~t.column("decoy").to_numpy(zero_copy_only=False)
    prot = t.column("proteins").filter(pa.array(keep)).combine_chunks()
    flat = pc.starts_with(pc.list_flatten(prot), "ENTRAP_").to_numpy(zero_copy_only=False)
    parent = pc.list_parent_indices(prot).to_numpy()
    n = len(prot)
    human_any = np.bincount(parent, weights=~flat, minlength=n) > 0
    has = np.bincount(parent, minlength=n) > 0
    e = int((has & ~human_any).sum())
    return {"entrapment_targets": e, "human_targets": n - e, "ratio": e / (n - e)}


def winners(arm):
    rows = [r for r in csv.DictReader(open(f"{runs}/{arm}/precursors.tsv"), delimiter="\t")
            if r["decoy"] == "false"]
    rows.sort(key=lambda r: -float(r["discriminant_score"]))
    rows.sort(key=lambda r: float(r["run_q_value"]))
    return rows


def evaluate(arm, k=None):
    rows = winners(arm)
    acc = [r for r in rows if float(r["run_q_value"]) <= 0.01]
    res = {"N": len(acc), "wall_s": wall(arm), "maxrss_gb": maxrss_gb(arm)}
    if k is not None:
        e = [is_entrap(r["proteins"]) for r in rows]
        ne = sum(e[: len(acc)])
        res.update(N_e=ne, FDP_c=ne * k / len(acc) if acc else None)
        best, run_e = 0, 0
        for i, x in enumerate(e, 1):
            run_e += x
            if run_e * k / i <= 0.01:
                best = i
        res["N_at_FDPc_1pct"] = best
        res["N_human"] = len(acc) - ne
    return res


arms = {a for a in ("P", "F", "Pe", "Fe", "P1", "P2", "H") if glob.glob(f"{runs}/{a}/precursors.tsv")}
R = {"arms_present": sorted(arms)}
for a in ("Pe", "Fe"):
    c = glob.glob(f"{cache}/{a}/*.candidates.parquet")
    if c:
        R[f"candidates_{a}"] = library_ratio(c[0])
if "candidates_Pe" in R:
    r0 = R["candidates_Pe"]["ratio"]
    if "candidates_Fe" in R:
        assert R["candidates_Fe"] == R["candidates_Pe"], "Pe and Fe candidate sets differ"
    R["r0"] = r0
    k = 1 + 1 / r0
    lib = glob.glob(f"{cache}/Pe/*.library.parquet")
    if lib:
        R["final_library_ratio_Pe_sensitivity"] = library_ratio(lib[0])
for a in sorted(arms):
    R[a] = evaluate(a, k if a in ("Pe", "Fe") and "r0" in R else None)

C = {}
if {"Pe", "Fe"} <= arms and "r0" in R:
    C["2_fdr"] = {"FDP_c_Pe": R["Pe"]["FDP_c"], "FDP_c_Fe": R["Fe"]["FDP_c"],
                  "pass": R["Pe"]["FDP_c"] <= R["Fe"]["FDP_c"] + 0.003}
if {"P", "F"} <= arms:
    c3 = {"N_P": R["P"]["N"], "N_F": R["F"]["N"], "ratio_P_F": R["P"]["N"] / R["F"]["N"]}
    c3["pass_N"] = R["P"]["N"] >= 0.97 * R["F"]["N"]
    if "2_fdr" in C:
        c3["N_at_FDPc_Pe"], c3["N_at_FDPc_Fe"] = R["Pe"]["N_at_FDPc_1pct"], R["Fe"]["N_at_FDPc_1pct"]
        c3["ratio_Pe_Fe"] = c3["N_at_FDPc_Pe"] / c3["N_at_FDPc_Fe"]
        c3["pass_N_at_FDPc"] = c3["N_at_FDPc_Pe"] >= 0.97 * c3["N_at_FDPc_Fe"]
        c3["pass"] = c3["pass_N"] and c3["pass_N_at_FDPc"]
    C["3_ids"] = c3
    C["4_cost"] = {"wall_P_s": R["P"]["wall_s"], "wall_F_s": R["F"]["wall_s"],
                   "ratio": R["P"]["wall_s"] / R["F"]["wall_s"], "pass": R["P"]["wall_s"] <= 0.5 * R["F"]["wall_s"]}
    C["goal_N_P_ge_24669"] = {"N_P": R["P"]["N"], "met": R["P"]["N"] >= 24669}
if {"P", "P1", "P2"} <= arms:
    ns = [R[a]["N"] for a in ("P", "P1", "P2")]
    C["seed_spread_P"] = {"N": ns, "range": max(ns) - min(ns)}
R["criteria"] = C
json.dump(R, open(out, "w"), indent=1)
print(json.dumps(R, indent=1))
