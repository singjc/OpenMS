#!/usr/bin/env python3
"""PLAN_SURVEY evaluation. Usage: eval_survey.py <survey_v3 dir> <truth dir> <candidates.parquet> <out.json>

<survey_v3 dir>/<V>-<RUN>/*.survey.tsv.gz and *.survey.features.tsv.gz (V0-V3, D/H).
<truth dir>/{F,FH}/precursors.slim.tsv.gz: the full-library IDs on D (F) and H (FH), target, q <= 0.01.
Pair k is precursors 2k (target) and 2k+1 (decoy); survey pair statistics are max(T, D).

Primary: recall at K = 800k of the truth IDs. Secondary: recall of charge 3-4, the empirical E_c (D only), the
excess target wins among kept pairs. +pool: max over charge siblings of S, plus 0.5 x the protein's top-3 mean.
+learn: target (1) vs decoy (0) logistic score on the survey features, 3-fold cross-fitted by family (FNV-1a of
the unmodified sequence mod 3); the pair statistic is max(f(T), f(D)), a symmetric function of the pair. It uses
no truth labels.
"""
import csv, gzip, json, re, sys, collections, glob
import numpy as np
import pyarrow as pa, pyarrow.parquet as pq

root, truth_dir, cand, out = sys.argv[1:5]
K = 800_000


def fnv1a(s):
    h = 0xcbf29ce484222325
    for b in s.encode():
        h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h


def load_survey(path):
    traml, T, D = [], [], []
    with gzip.open(path, "rt") as f:
        r = csv.reader(f, delimiter="\t"); next(r)
        for row in r:
            traml.append(row[1]); T.append(float(row[2])); D.append(float(row[3]))
    return traml, np.array(T), np.array(D)


def truth(arm, n):
    t = np.zeros(n, bool)
    for r in csv.DictReader(gzip.open(f"{truth_dir}/{arm}/precursors.slim.tsv.gz", "rt"), delimiter="\t"):
        if r["decoy"] == "false" and float(r["run_q_value"]) <= 0.01:
            t[int(r["precursor_id"]) // 2] = True
    return t


# proteins per pair (from the candidate library; IDs are deterministic across variants: same enumeration key part)
tb = pq.read_table(cand, columns=["precursor_index", "precursor_id", "decoy", "proteins"]).combine_chunks()
idx = tb.column("precursor_index").to_numpy()
first = np.ones(len(idx), bool); first[1:] = idx[1:] != idx[:-1]
sub = tb.filter(pa.array(first & ~tb.column("decoy").to_numpy(zero_copy_only=False)))
pair_prot = {}
for p, ps in zip(sub.column("precursor_id").to_pylist(), sub.column("proteins").to_pylist()):
    pair_prot[int(p) // 2] = ps
del tb, sub

emp = set(r["sequence"] + "/" + r["charge"] for r in csv.DictReader(open(
    "/home/singjc/.jcode/scratch/entity-inference-20261127/out_pasef_owned/legacy/precursors.tsv"), delimiter="\t"))


def group_max(keys, x):
    best = collections.defaultdict(float)
    for k, v in zip(keys, x):
        if v > best[k]: best[k] = v
    return np.array([best[k] for k in keys])


def protein_top3(x, n):
    per = collections.defaultdict(list)
    for k in range(n):
        for p in pair_prot.get(k, ()):
            per[p].append(x[k])
    score = {p: float(np.mean(sorted(v, reverse=True)[:3])) for p, v in per.items()}
    return np.array([max((score[p] for p in pair_prot.get(k, ())), default=0.0) for k in range(n)])


def metrics(stat, T, D, tr, charge, empm):
    n = len(stat)
    order = np.lexsort((np.arange(n), -stat))
    kept = np.zeros(n, bool); kept[order[:K]] = True
    kept &= stat > 0
    k = np.flatnonzero(kept)
    res = {"recall": round(float((tr & kept).sum() / tr.sum()), 4),
           "recall_charge34": round(float((tr & kept & (charge >= 3)).sum() / (tr & (charge >= 3)).sum()), 4),
           "excess_target_wins": int((T[k] > D[k]).sum() - (D[k] > T[k]).sum()),
           "kept": int(kept.sum())}
    if empm is not None:
        res["recall_E"] = round(float((empm & kept).sum() / empm.sum()), 4)
    return res


def learned(feat_path, n, charge, plain):
    """3-fold cross-fitted logistic target-vs-decoy score per precursor; pair stat = max(f(T), f(D))."""
    cols = ["coelution", "support", "matched_scans", "apex_rt_error_s", "fwhm_s", "apex_intensity", "peaks_in_window"]
    X = np.zeros((2 * n, len(cols) + 4)); have = np.zeros(2 * n, bool)
    with gzip.open(feat_path, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            i = int(r["precursor_id"])
            v = [float(r[c]) for c in cols]
            z = int(r["charge"])
            X[i] = [v[0], v[1], np.log1p(v[2]), abs(v[3]), v[4], np.log1p(v[5]), np.log1p(v[6]), z == 2, z == 3, z >= 4, 1.0]
            have[i] = True
    y = np.zeros(2 * n); y[0::2] = 1.0
    fold = np.array([fnv1a(plain[k]) % 3 for k in range(n)]).repeat(2)
    mu = X[have].mean(0); sd = X[have].std(0); sd[sd == 0] = 1; mu[-1] = 0; sd[-1] = 1
    Z = (X - mu) / sd
    score = np.full(2 * n, -np.inf)
    for f in range(3):
        tr_m = have & (fold != f); te = have & (fold == f)
        w = np.zeros(Z.shape[1])
        A, b = Z[tr_m], y[tr_m]
        for _ in range(200):  # Newton-free gradient steps on the logistic loss, L2 1e-3
            p = 1 / (1 + np.exp(-A @ w))
            g = A.T @ (p - b) / len(b) + 1e-3 * w
            w -= 2.0 * g
        score[te] = Z[te] @ w
    st = np.maximum(score[0::2], score[1::2])
    st[~(have[0::2] | have[1::2])] = -np.inf
    # shift so that every surveyed pair is > 0 (the metric keeps stat > 0 only)
    finite = np.isfinite(st)
    st[finite] = st[finite] - st[finite].min() + 1e-6
    st[~finite] = 0.0
    return st, score[0::2], score[1::2]


R = {"K": K}
ref_traml = {}  # pair index -> TRAML_ID from V0 (pairs dropped by a variant have an empty TRAML_ID and S = 0)
for run, tarm in (("D", "F"), ("H", "FH")):
    for V in ("V0", "V1", "V2", "V3"):
        d = f"{root}/{V}-{run}"
        s = glob.glob(f"{d}/*.survey.tsv.gz")
        if not s or not glob.glob(f"{truth_dir}/{tarm}/precursors.slim.tsv.gz"):
            continue
        traml, T, D = load_survey(s[0])
        if V == "V0":
            ref_traml[run] = traml
        traml = [t if t else ref_traml[run][k] for k, t in enumerate(traml)]
        n = len(traml)
        charge = np.array([int(t.rsplit("/", 1)[1]) for t in traml])
        modseq = [t.rsplit("/", 1)[0] for t in traml]
        plain = [re.sub(r"\(UniMod:\d+\)", "", x) for x in modseq]
        tr = truth(tarm, n)
        empm = np.array([t in emp for t in traml]) if run == "D" else None
        S0 = np.maximum(T, D)
        r = {"base": metrics(S0, T, D, tr, charge, empm), "truth_ids": int(tr.sum())}
        pool = group_max(modseq, S0) + 0.5 * protein_top3(S0, n)
        r["pool"] = metrics(pool, T, D, tr, charge, empm)
        fp = glob.glob(f"{d}/*.survey.features.tsv.gz")
        if fp:
            L, fT, fD = learned(fp[0], n, charge, plain)
            r["learn"] = metrics(L, T, D, tr, charge, empm)
            # excess target wins under the learned score itself (label-free check of selection symmetry)
            order = np.lexsort((np.arange(n), -L)); k = order[:K]
            r["learn"]["excess_target_wins_learned_score"] = int((fT[k] > fD[k]).sum() - (fD[k] > fT[k]).sum())
            lp = group_max(modseq, L) + 0.5 * protein_top3(L, n)
            r["learn+pool"] = metrics(lp, T, D, tr, charge, empm)
        tf = glob.glob(f"{d}.time")
        if tf:
            for line in open(tf[0]):
                if "Elapsed" in line:
                    r["wall"] = line.rsplit(" ", 1)[1].strip()
        R[f"{V}-{run}"] = r
        print(V, run, json.dumps(r))
json.dump(R, open(out, "w"), indent=1)
