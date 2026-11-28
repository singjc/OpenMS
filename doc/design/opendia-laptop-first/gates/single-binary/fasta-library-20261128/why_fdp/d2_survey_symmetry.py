#!/usr/bin/env python3
"""D2 (diagnostic, existing outputs only): does the prefilter's selection make entrapment targets beat their decoys?

The prefilter keeps a pair if max(T, D) of its survey co-elution is large. Within a NULL pair (entrapment, no true
signal), T and D are exchangeable a priori. Selecting on max(T, D) keeps them exchangeable, so P(survey target > survey
decoy | kept) should be 0.5. If that holds but the final search still prefers entrapment targets (T_e/D_e > 1, see D1),
the extra target wins are produced by the final search, not by which pairs were selected.

Run on codon (survey and library caches live there). For each FDRBench prefilter arm:
  - kept entrapment pairs (in the final library): fraction with survey target_coelution > decoy_coelution
  - human pairs, for contrast
  - accepted (q <= 0.01) entrapment-pair winners: survey target_coelution > decoy_coelution fraction,
    split by whether the target or the decoy won the final competition.
Usage: d2_survey_symmetry.py OUT.json
"""
import csv, json, re
import numpy as np, pyarrow as pa, pyarrow.compute as pc, pyarrow.parquet as pq

W = "/nfs/research/juan/singj/opendia/fasta-20261128"
SPLIT = re.compile(r"[;,]")
R = {}
for arm in ("FbP", "FbV1", "FbPp"):
    st, sd = {}, {}
    with open(f"{W}/cache_fb/{arm}/human_fdrbench.survey.tsv") as f:
        next(f)
        for line in f:
            k, _, t, d, _ = line.rstrip("\n").split("\t")
            st[int(k)] = float(t); sd[int(k)] = float(d)
    lib = pq.read_table(f"{W}/cache_fb/{arm}/human_fdrbench.library.parquet", columns=["precursor_index", "precursor_id", "decoy", "proteins"]).combine_chunks()
    idx = lib.column("precursor_index").to_numpy()
    first = np.ones(len(idx), bool); first[1:] = idx[1:] != idx[:-1]
    sub = lib.filter(pa.array(first & ~lib.column("decoy").to_numpy(zero_copy_only=False)))
    del lib
    pid = np.array([int(x) // 2 for x in sub.column("precursor_id").to_pylist()])
    prot = sub.column("proteins").combine_chunks()
    flat = pc.ends_with(pc.list_flatten(prot), "_p_target").to_numpy(zero_copy_only=False)
    par = pc.list_parent_indices(prot).to_numpy()
    human_any = np.bincount(par, weights=~flat, minlength=len(prot)) > 0
    ent = dict(zip(pid.tolist(), (~human_any).tolist()))
    res = {}
    for name, sel in (("entrapment", True), ("human", False)):
        ks = [k for k, e in ent.items() if e == sel]
        tw = sum(st[k] > sd[k] for k in ks); dw = sum(sd[k] > st[k] for k in ks)
        res[f"kept_{name}"] = {"pairs": len(ks), "survey_target_gt_decoy": tw, "survey_decoy_gt_target": dw,
                               "frac_target": tw / max(tw + dw, 1)}
    # accepted entrapment winners
    best = {}
    with open(f"{W}/runs_fb/{arm}/precursors.tsv") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            k = int(row["precursor_id"]) // 2
            s = float(row["discriminant_score"]); dec = row["decoy"] == "true"
            if k not in best or (s, dec) > best[k][:2]:
                best[k] = (s, dec, float(row["run_q_value"]))
    acc = {"target_won": [0, 0], "decoy_won": [0, 0]}
    for k, (s, dec, q) in best.items():
        if q > 0.01 or not ent.get(k, False):
            continue
        key = "decoy_won" if dec else "target_won"
        acc[key][0 if st[k] > sd[k] else 1] += 1
    res["accepted_entrapment"] = {k: {"survey_target_gt_decoy": v[0], "survey_decoy_ge_target": v[1]} for k, v in acc.items()}
    R[arm] = res
    print(arm, json.dumps(res))
json.dump(R, open("why_fdp/d2.json", "w"), indent=1)
