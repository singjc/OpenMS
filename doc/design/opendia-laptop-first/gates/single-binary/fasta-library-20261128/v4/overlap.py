# Diagnostic: P' (V1+pool) IDs vs F (K = 0) IDs on D, by own V1 survey rank; and whether P'-only IDs are
# "related" to another P' ID (nested sequence: one unmodified sequence contains the other, e.g. missed cleavage).
import csv, gzip, re, json
import numpy as np
def ids(path):
    s = {}
    for r in csv.DictReader(gzip.open(path, "rt"), delimiter="\t"):
        if r["decoy"] == "false" and float(r["run_q_value"]) <= 0.01: s[int(r["precursor_id"]) // 2] = float(r["run_q_value"])
    return s
P = ids("v4/slim/runs_v4/Pp/precursors.slim.tsv.gz"); F = ids("v2/slim/runs_v2/F/precursors.slim.tsv.gz")
traml, S = [], []
with gzip.open("survey_v3/V1-D/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.survey.tsv.gz", "rt") as f:
    next(f)
    for line in f:
        r = line.rstrip("\n").split("\t"); traml.append(r[1]); S.append(float(r[4]))
S = np.array(S); n = len(S)
rank = np.empty(n, np.int64); rank[np.lexsort((np.arange(n), -S))] = np.arange(n)
plain = lambda k: re.sub(r"\(UniMod:\d+\)", "", traml[k].rsplit("/", 1)[0])
pseq = {plain(k) for k in P}
def related(k):
    s = plain(k)
    return any((s != t) and (s in t or t in s) for t in pseq if abs(len(s) - len(t)) <= 12)
R = {"N_P": len(P), "N_F": len(F), "both": len(P.keys() & F.keys())}
for name, lo, hi in (("own<=800k", 0, 800_000), ("own>800k", 800_000, n + 1)):
    ks = [k for k in P if lo <= rank[k] < hi]
    only = [k for k in ks if k not in F]
    rel_only = sum(related(k) for k in only)
    both = [k for k in ks if k in F]
    import random; random.seed(0); smp = random.sample(both, min(len(both), len(only)))
    R[name] = {"N": len(ks), "in_F": len(both), "P_only": len(only), "P_only_frac": round(len(only) / max(len(ks), 1), 3),
               "P_only_nested_related": rel_only, "in_F_nested_related_sample": f"{sum(related(k) for k in smp)}/{len(smp)}"}
F_only = [k for k in F if k not in P]
R["F_only"] = len(F_only)
print(json.dumps(R, indent=1))
