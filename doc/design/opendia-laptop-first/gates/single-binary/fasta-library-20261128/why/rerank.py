#!/usr/bin/env python3
"""Offline exploration (no rule, after results): pair-level, target-decoy-symmetric re-rankings of the P survey.
Each statistic is a function of pair statistics S_k = max(T_k, D_k) and label-free pair attributes (charge,
sequence family, proteins), so for a null pair swapping T_k and D_k leaves the selection unchanged.
Metric: recall at K of F's target IDs (q <= 0.01, full library) and of the empirical IDs E_c. F's IDs are a proxy
truth (about 1-1.6 % false; biased toward what the full search finds).
"""
import csv, gzip, json, re, collections
import numpy as np
import pyarrow as pa, pyarrow.compute as pc, pyarrow.parquet as pq

S = "/home/singjc/.jcode/scratch/fasta-library-20261128"
Ks = [300_000, 500_000, 800_000, 1_200_000]

traml, T, D = [], [], []
with gzip.open(f"{S}/why/survey_P.tsv.gz", "rt") as f:
    r = csv.reader(f, delimiter="\t"); next(r)
    for row in r:
        traml.append(row[1]); T.append(float(row[2])); D.append(float(row[3]))
T = np.array(T); D = np.array(D); n = len(traml)
S0 = np.maximum(T, D)
charge = np.array([int(t.rsplit("/", 1)[1]) for t in traml])
modseq = [t.rsplit("/", 1)[0] for t in traml]
plain = [re.sub(r"\(UniMod:\d+\)", "", s) for s in modseq]

# proteins of each pair's target (first row of each target precursor; pair k = id // 2)
t = pq.read_table(f"{S}/design_D/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.candidates.parquet",
                  columns=["precursor_index", "precursor_id", "decoy", "proteins"]).combine_chunks()
idx = t.column("precursor_index").to_numpy()
first = np.ones(len(idx), bool); first[1:] = idx[1:] != idx[:-1]
keep = first & ~t.column("decoy").to_numpy(zero_copy_only=False)
sub = t.filter(pa.array(keep))
pids = np.array([int(x) for x in sub.column("precursor_id").to_pylist()])
prots = sub.column("proteins").to_pylist()
del t, sub
pair_prot = [None] * n
for p, ps in zip(pids, prots):
    pair_prot[p // 2] = ps


def accepted(arm):
    out = set()
    for r in csv.DictReader(gzip.open(f"{S}/v2/slim/runs_v2/{arm}/precursors.slim.tsv.gz", "rt"), delimiter="\t"):
        if r["decoy"] == "false" and float(r["run_q_value"]) <= 0.01:
            out.add(int(r["precursor_id"]) // 2)
    return out


truthF = np.zeros(n, bool); truthF[list(accepted("F"))] = True
emp = set(r["sequence"] + "/" + r["charge"] for r in csv.DictReader(open(
    "/home/singjc/.jcode/scratch/entity-inference-20261127/out_pasef_owned/legacy/precursors.tsv"), delimiter="\t"))
truthE = np.array([t in emp for t in traml])


def group_max(keys, x):
    best = collections.defaultdict(float)
    for k, v in zip(keys, x):
        if v > best[k]: best[k] = v
    return np.array([best[k] for k in keys])


def within_charge_quantile(x):
    q = np.zeros(n)
    for z in np.unique(charge):
        m = charge == z
        q[m] = (np.argsort(np.argsort(x[m], kind="stable"), kind="stable") + 1) / m.sum()
    return q


def protein_evidence(x, top=3):
    """Mean of the top-3 pair statistics of each protein, excluding the pair itself is skipped (cheap proxy)."""
    per = collections.defaultdict(list)
    for k, ps in enumerate(pair_prot):
        for p in ps or ():
            per[p].append(x[k])
    score = {p: float(np.mean(sorted(v, reverse=True)[:top])) for p, v in per.items()}
    return np.array([max((score[p] for p in ps), default=0.0) for ps in pair_prot])


def recall(stat):
    order = np.lexsort((np.arange(n), -stat))
    out = {}
    for K in Ks:
        kept = np.zeros(n, bool); kept[order[:K]] = True
        out[str(K)] = {"F": round(float(truthF[kept].sum() / truthF.sum()), 4),
                       "E": round(float(truthE[kept].sum() / truthE.sum()), 4),
                       "F_charge34": round(float((truthF & kept & (charge >= 3)).sum() / (truthF & (charge >= 3)).sum()), 4)}
    return out


q0 = within_charge_quantile(S0)
sib = group_max(modseq, S0)            # same modified sequence, all charges
fam = group_max(plain, S0)             # same unmodified sequence (mods + charges)
famq = group_max(plain, q0)
prot = protein_evidence(S0)
rules = {
    "R0 current S": S0,
    "R1 within-charge quantile": q0,
    "R2 charge siblings max": sib,
    "R3 peptide family max (mods+charges)": fam,
    "R4 family max of within-charge quantile": famq,
    "R5 S + 0.5 x protein top-3 mean": S0 + 0.5 * prot,
    "R6 family quantile + 0.5 x protein quantile": famq + 0.5 * within_charge_quantile(prot),
}
res = {name: recall(x) for name, x in rules.items()}
for name, x in rules.items():
    order = np.lexsort((np.arange(n), -x))
    res[name]["excess_target_wins"] = {}
    for K in Ks:
        k = order[:K]
        res[name]["excess_target_wins"][str(K)] = int((T[k] > D[k]).sum() - (D[k] > T[k]).sum())
res["_excess_all_pairs"] = int((T > D).sum() - (D > T).sum())
res["_truth"] = {"F_ids": int(truthF.sum()), "E_c": int(truthE.sum()), "pairs": n}
json.dump(res, open(f"{S}/why/rerank.json", "w"), indent=1)
for k, v in res.items():
    print(k, v)
