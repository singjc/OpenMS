#!/usr/bin/env python3
"""PLAN_ENTRAP_FDRBENCH evaluation: the eval_fasta.py estimator with a configurable entrapment label, plus Poisson CIs.

Usage: eval_fb.py <label_mode> <runs_dir> <ratios.json> <out.json> ARM...
  label_mode: "entrap" (custom: every protein starts with ENTRAP_) or "fdrbench" (every protein ends with _p_target)
  runs_dir/<ARM>/precursors.slim.tsv.gz (or precursors.tsv), runs_dir/<ARM>.time
  ratios.json: {"candidates_<ARM>": {...ratio}, "library_<ARM>": {...ratio}} from ratio_fb.sh / ratio_v2.sh

Per arm: N (target precursor winners at run q <= 0.01), N_e, FDP_c(r0) and FDP_c(r_lib) with
FDP_c = N_e (1 + 1/r) / N, 95 % Poisson CI on N_e, N@FDP_c 1 % (same procedure as eval_fasta.py), wall, peak RSS.
"""
import csv, glob, gzip, json, math, re, sys

mode, runs, ratios, out = sys.argv[1:5]
arms = sys.argv[5:]
SPLIT = re.compile(r"[;,]")
if mode == "entrap":
    member = lambda p: p.startswith("ENTRAP_")
elif mode == "fdrbench":
    member = lambda p: p.endswith("_p_target")
else:
    sys.exit("label_mode must be entrap or fdrbench")


def is_entrap(proteins):
    ps = [p for p in SPLIT.split(proteins) if p]
    return bool(ps) and all(member(p) for p in ps)


def timefield(arm, key):
    for line in open(f"{runs}/{arm}.time"):
        if key in line:
            return line.rsplit(" ", 1)[1].strip()


def wall(arm):
    t = [float(x) for x in timefield(arm, "Elapsed (wall clock)").split(":")]
    return sum(v * 60 ** i for i, v in enumerate(reversed(t)))


def winners(arm):
    slim = f"{runs}/{arm}/precursors.slim.tsv.gz"
    f = gzip.open(slim, "rt") if glob.glob(slim) else open(f"{runs}/{arm}/precursors.tsv")
    rows = [r for r in csv.DictReader(f, delimiter="\t") if r["decoy"] == "false"]
    rows.sort(key=lambda r: -float(r["discriminant_score"]))
    rows.sort(key=lambda r: float(r["run_q_value"]))
    return rows


def fdp(ne, n, r):
    k = 1 + 1 / r
    lo, hi = max(ne - 1.96 * math.sqrt(ne), 0), ne + 1.96 * math.sqrt(ne)
    return {"FDP_c": ne * k / n, "ci95": [lo * k / n, hi * k / n]}


def n_at(e, r):
    k, best, run_e = 1 + 1 / r, 0, 0
    for i, x in enumerate(e, 1):
        run_e += x
        if run_e * k / i <= 0.01:
            best = i
    return best


RT = json.load(open(ratios))
R = {"label_mode": mode}
for a in arms:
    rows = winners(a)
    acc = sum(float(r["run_q_value"]) <= 0.01 for r in rows)
    e = [is_entrap(r["proteins"]) for r in rows]
    ne = sum(e[:acc])
    res = {"N": acc, "N_e": ne, "N_human": acc - ne, "wall_s": wall(a),
           "maxrss_gb": int(timefield(a, "Maximum resident")) / 1024 ** 2}
    for kind in ("candidates", "library"):
        key = f"{kind}_{a}"
        if key in RT:
            r = RT[key]["ratio"]
            tag = "r0" if kind == "candidates" else "r_lib"
            res[tag] = r
            res[f"{tag}_fdp"] = fdp(ne, acc, r)
            res[f"{tag}_N_at_FDPc_1pct"] = n_at(e, r)
    R[a] = res
json.dump(R, open(out, "w"), indent=1)
print(json.dumps(R, indent=1))
