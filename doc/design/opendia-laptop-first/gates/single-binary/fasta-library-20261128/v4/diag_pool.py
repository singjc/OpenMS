#!/usr/bin/env python3
"""Diagnostic (not a rule): where do Pe' (V1+pool) IDs and entrapment hits come from, and what do pooling variants
cost in entrapment asymmetry? Run on codon: uv run --with pyarrow --with numpy python diag_pool.py out.json

1. Strata by the pair's OWN V1 survey rank (no pooling): own-admitted (rank <= 800k) vs boost-admitted. Per stratum:
   library entrapment/human targets (r_b), target IDs and entrapment IDs at q <= 0.01, decoy winners at q <= 0.01
   (a per-stratum TDC view), and a stratified entrapment FDP: sum_b N_e,b (1 + 1/r_b) / N.
2. Trade-off at K = 800k on the Pe' survey: base, sibling-only, and sibling + w x protein for several w. Per variant:
   r_sel = entrapment/human targets among kept pairs, and recall of F's human IDs (matched by TRAML_ID).
"""
import csv, gzip, json, re, sys
import numpy as np, pyarrow as pa, pyarrow.compute as pc, pyarrow.parquet as pq

W = "/nfs/research/juan/singj/opendia/fasta-20261128"
K = 800_000
out = sys.argv[1]
SPLIT = re.compile(r"[;,]")


def is_entrap(proteins):
    ps = [p for p in SPLIT.split(proteins) if p]
    return bool(ps) and all(p.startswith("ENTRAP_") for p in ps)


# Pe' survey: own statistic per pair, TRAML_ID of the target
traml, S = [], []
with open(f"{W}/cache_v4/Pep/human_entrapment.survey.tsv") as f:
    next(f)
    for line in f:
        r = line.rstrip("\n").split("\t"); traml.append(r[1]); S.append(float(r[4]))
S = np.array(S); n = len(S)
present = np.array([bool(t) for t in traml])

# proteins per pair (targets, first row per precursor) from the candidates
tb = pq.read_table(f"{W}/cache_v4/Pep/human_entrapment.candidates.parquet",
                   columns=["precursor_index", "precursor_id", "decoy", "proteins"]).combine_chunks()
idx = tb.column("precursor_index").to_numpy()
first = np.ones(len(idx), bool); first[1:] = idx[1:] != idx[:-1]
sub = tb.filter(pa.array(first & ~tb.column("decoy").to_numpy(zero_copy_only=False)))
del tb
pair_of_row = np.array([int(x) // 2 for x in sub.column("precursor_id").to_pylist()])
prot = sub.column("proteins").combine_chunks()
flat = pc.list_flatten(prot)
parent = pair_of_row[pc.list_parent_indices(prot).to_numpy()]
entrap_flat = pc.starts_with(flat, "ENTRAP_").to_numpy(zero_copy_only=False)
human_any = np.bincount(parent, weights=~entrap_flat, minlength=n) > 0
has = np.bincount(parent, minlength=n) > 0
entrap_pair = has & ~human_any
pid = pc.dictionary_encode(flat).indices.to_numpy()
del sub, prot, flat

# sibling max (same target modified sequence), present pairs only
seq = np.array([t.rsplit("/", 1)[0] if t else "" for t in traml])
_, seq_id = np.unique(seq, return_inverse=True)
sib = np.zeros(seq_id.max() + 1); np.maximum.at(sib, seq_id[present], S[present])
sib_pair = np.where(present, sib[seq_id], 0.0)
# protein top-3 mean of S over present pairs
m = present[parent]
pp, pr, ps = parent[m], pid[m], S[parent[m]]
o = np.lexsort((-ps, pr)); pr_s, ps_s = pr[o], ps[o]
start = np.r_[0, np.flatnonzero(pr_s[1:] != pr_s[:-1]) + 1]
rank = np.arange(len(pr_s)) - np.repeat(start, np.diff(np.r_[start, len(pr_s)]))
top = rank < 3
P = pid.max() + 1
ev = np.bincount(pr_s[top], weights=ps_s[top], minlength=P) / np.maximum(np.bincount(pr_s[top], minlength=P), 1)
prot_pair = np.zeros(n); np.maximum.at(prot_pair, pp, ev[pr])


def kept(stat):
    k = np.zeros(n, bool); k[np.lexsort((np.arange(n), -stat))[:K]] = True
    return k & (stat > 0)


# F's human IDs as TRAML_IDs (human candidate pair -> TRAML_ID from the V0-D survey)
htraml = []
with gzip.open(f"{W}/survey_v3/V0-D/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.survey.tsv.gz", "rt") as f:
    next(f)
    for line in f: htraml.append(line.split("\t", 2)[1])
truth = set()
for r in csv.DictReader(open(f"{W}/runs_v2/F/precursors.tsv"), delimiter="\t"):
    if r["decoy"] == "false" and float(r["run_q_value"]) <= 0.01: truth.add(htraml[int(r["precursor_id"]) // 2])
tr = np.array([t in truth for t in traml])

R = {"pairs": int(n), "truth_F": len(truth), "truth_in_Pe_candidates": int(tr.sum())}
variants = {"base": S, "sibling_only": np.where(present, sib_pair, 0.0)}
for w in (0.1, 0.25, 0.5):
    variants[f"sibling+{w}xprotein"] = np.where(present, sib_pair + w * prot_pair, 0.0)
variants["protein_only_0.5"] = np.where(present, S + 0.5 * prot_pair, 0.0)
cand_r = entrap_pair.sum() / (has & ~entrap_pair).sum()
R["r0_candidates"] = float(cand_r)
T = {}
for name, st in variants.items():
    k = kept(st)
    e, h = int((k & entrap_pair).sum()), int((k & has & ~entrap_pair).sum())
    T[name] = {"recall_F": round(float((k & tr).sum() / tr.sum()), 4), "r_sel": round(e / h, 4),
               "r_sel_over_r0": round(e / h / cand_r, 3)}
R["tradeoff"] = T

# strata for the actual Pe' run
own_rank = np.empty(n, np.int64); own_rank[np.lexsort((np.arange(n), -S))] = np.arange(n)
lib = pq.read_table(f"{W}/cache_v4/Pep/human_entrapment.library.parquet", columns=["precursor_id"]).column("precursor_id")
in_lib = np.zeros(n, bool); in_lib[np.unique(np.array([int(x) // 2 for x in pc.unique(lib).to_pylist()]))] = True
edges = [0, 800_000, 1_600_000, 3_200_000, n + 1]
names = ["own<=800k", "800k-1.6M", "1.6M-3.2M", ">3.2M"]
bin_of = np.searchsorted(edges, own_rank, side="right") - 1
ids = {b: [0, 0, 0] for b in range(4)}  # targets, entrapment targets, decoys at q <= 0.01
for r in csv.DictReader(open(f"{W}/runs_v4/Pep/precursors.tsv"), delimiter="\t"):
    if float(r["run_q_value"]) > 0.01: continue
    b = bin_of[int(r["precursor_id"]) // 2]
    if r["decoy"] == "true": ids[b][2] += 1
    else:
        ids[b][0] += 1; ids[b][1] += is_entrap(r["proteins"])
N = sum(v[0] for v in ids.values())
st, strat = {}, 0.0
for b in range(4):
    e = int((in_lib & (bin_of == b) & entrap_pair).sum()); h = int((in_lib & (bin_of == b) & has & ~entrap_pair).sum())
    rb = e / h if h else None
    contrib = ids[b][1] * (1 + 1 / rb) if rb else None
    if contrib: strat += contrib
    st[names[b]] = {"lib_entrap": e, "lib_human": h, "r_b": rb, "N": ids[b][0], "N_e": ids[b][1], "decoys_q01": ids[b][2],
                    "est_false": contrib}
R["strata"] = st
R["N"] = N
R["FDP_stratified"] = strat / N
# Fine strata: 20 equal-count bins of the library pairs by own survey rank (finer control of selection on own noise).
lib_rank = own_rank[in_lib]
qs = np.quantile(lib_rank, np.linspace(0, 1, 21)); qs[-1] = n + 1
fine = np.searchsorted(qs, own_rank, side="right") - 1
fine = np.clip(fine, 0, 19)
fe_ids = np.zeros((20, 3), int)
for r in csv.DictReader(open(f"{W}/runs_v4/Pep/precursors.tsv"), delimiter="\t"):
    if float(r["run_q_value"]) > 0.01: continue
    b = fine[int(r["precursor_id"]) // 2]
    if r["decoy"] == "true": fe_ids[b, 2] += 1
    else:
        fe_ids[b, 0] += 1; fe_ids[b, 1] += is_entrap(r["proteins"])
fs_, rows = 0.0, []
for b in range(20):
    sel = in_lib & (fine == b)
    e = int((sel & entrap_pair).sum()); h = int((sel & has & ~entrap_pair).sum())
    c = fe_ids[b, 1] * (1 + h / e) if e else 0.0
    fs_ += c
    rows.append({"own_rank_from": int(qs[b]), "r_b": round(e / h, 4) if h else None, "N": int(fe_ids[b, 0]),
                 "N_e": int(fe_ids[b, 1]), "decoys_q01": int(fe_ids[b, 2]), "est_false": round(c, 1)})
R["fine_strata"] = rows
R["FDP_stratified_fine"] = fs_ / N
R["decoys_q01_total"] = int(fe_ids[:, 2].sum())
R["note"] = "strata with r_b ~ 0 and N_e = 0 are not covered by entrapment"
json.dump(R, open(out, "w"), indent=1)
print(json.dumps(R, indent=1))
