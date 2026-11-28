#!/usr/bin/env python3
"""PLAN_RESIDUAL_DIAG: where in the pipeline do accepted null (entrapment) target winners get their edge over decoys?

Per arm, every accepted null pair winner (run q <= 0.01; FDRBench entrapment pair) is put in exactly one category,
using the winning member's own apex and fragments (same rule for target and decoy winners):
  A      the FDRBench parent (any charge) is an accepted target in this arm and its RT span contains the winner apex
  B_lib  not A; the parent is accepted in the K = 0 reference FbF with a span containing the apex (the parent is
         present in the sample and co-elutes), and the parent is in this arm's searched library (scored, not accepted)
  B_cut  as B_lib, but the parent is not in this arm's searched library (the survey prefilter cut it)
  C      not A/B; another accepted target of this arm co-elutes (span contains apex), |d m/z| <= 25, |d 1/K0| <= 0.05,
         and shares >= 1 of the winner's 6 search fragments within 10 ppm
  D      none of the above
Reported per arm and pooled: T_c, D_c, T_c - D_c and the share (T_c - D_c) / (T - D).
Q2 (arms with -features_out): interference features of the winner window (i_explained, i_nclaim, ...) and the
standardised mean difference of all features, target vs decoy winners, by category group (A+B, C+D).
Q3 (V1 arms with survey.tsv): share of accepted null pairs whose survey target co-elution exceeds the decoy's.
Q4 (category A): fragment share with the parent and whether the parent outscores the winner (why post hoc NC missed).
DRY=1 prints join and mapping rates only and never splits by target/decoy.
Usage: diag_residual.py OUT.json ARM:precursors.tsv:library.search.parquet:survey.tsv|-:features.tsv|- ...
"""
import bisect, csv, json, math, os, re, sys
import numpy as np
import pyarrow as pa, pyarrow.compute as pc, pyarrow.csv as pcsv, pyarrow.parquet as pq

W = "/nfs/research/juan/singj/opendia/fasta-20261128"
REF = f"{W}/runs_fb/FbF/precursors.tsv"
MZ, IM, PPM = 25.0, 0.05, 10e-6
DRY = os.environ.get("DRY", "0") == "1"
SPLIT = re.compile(r"[;,]")
MOD = re.compile(r"\(UniMod:\d+\)|\[[^\]]*\]")
COLS = ["precursor_id", "sequence", "unmodified_sequence", "charge", "proteins", "rt_apex", "rt_start", "rt_end",
        "discriminant_score", "run_q_value", "decoy"]

prot, acc = {}, None
for line in open(f"{W}/fasta/human_fdrbench.fasta"):
    if line.startswith(">"):
        parts = line[1:].split()[0].split("|")
        acc = parts[1] if len(parts) >= 3 and parts[0] in ("sp", "tr") else parts[0]
        prot[acc] = []
    else:
        prot[acc].append(line.strip())
prot = {k: "".join(v) for k, v in prot.items()}


def parents(pep, proteins):
    res = set()
    for p in proteins:
        if not p.endswith("_p_target"): continue
        es, hs = prot.get(p), prot.get(p[:-len("_p_target")])
        if not es or not hs or len(es) != len(hs): continue
        i = es.find(pep)
        while i >= 0:
            res.add(hs[i:i + len(pep)]); i = es.find(pep, i + 1)
    return res


def load(path):
    t = pcsv.read_csv(path, parse_options=pcsv.ParseOptions(delimiter="\t"),
                      convert_options=pcsv.ConvertOptions(include_columns=COLS, strings_can_be_null=False,
                                                          column_types={"precursor_id": pa.int64(), "decoy": pa.string(),
                                                                        "proteins": pa.string(), "sequence": pa.string()}))
    return {c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS}


def accepted_spans(x):
    """unmodified sequence -> [(start, end, score, precursor_id)] of accepted (q <= 0.01) target rows."""
    out = {}
    m = (x["run_q_value"] <= 0.01) & (x["decoy"] == "false")
    for i in np.where(m)[0]:
        out.setdefault(x["unmodified_sequence"][i], []).append(
            (x["rt_start"][i], x["rt_end"][i], x["discriminant_score"][i], int(x["precursor_id"][i])))
    return out


def fragments(lib, ids):
    want = pa.array([str(i) for i in ids])
    f = pq.ParquetFile(lib); out = {}
    for g in range(f.num_row_groups):
        t = f.read_row_group(g, columns=["precursor_id", "precursor_mz", "precursor_im", "product_mz", "fragment_rank"])
        t = t.filter(pc.is_in(t.column("precursor_id"), value_set=want))
        for pid, mz, im, fm, rk in zip(*(t.column(c).to_pylist() for c in t.column_names)):
            out.setdefault(int(pid), [mz, im, []])[2].append((rk, fm))
    for e in out.values():
        e[2] = np.sort(np.array([m for r, m in sorted(e[2]) if r <= 6]))
    return out


def nshared(cf, af):
    if not len(cf) or not len(af): return 0
    i = np.searchsorted(af, cf)
    lo = np.abs(af[np.clip(i - 1, 0, len(af) - 1)] - cf); hi = np.abs(af[np.clip(i, 0, len(af) - 1)] - cf)
    return int((np.minimum(lo, hi) <= PPM * cf).sum())


print("loading reference FbF", flush=True)
REFSPANS = accepted_spans(load(REF))
CATS = ["A", "B_lib", "B_cut", "C", "D"]
R = {"params": {"MZ": MZ, "IM": IM, "PPM": PPM, "dry": DRY}}


def run(arm, path, lib, survey, feats):
    x = load(path)
    n = len(x["precursor_id"])
    dec = x["decoy"] == "true"
    s = x["discriminant_score"].astype(float)
    pair = x["precursor_id"] // 2
    best = {}
    for i in range(n):
        k = (int(pair[i]), bool(dec[i]))
        if k not in best or s[i] > s[best[k]]: best[k] = i
    pairs = {}
    for (p, d), i in best.items(): pairs.setdefault(p, [None, None])[int(d)] = i
    spans = accepted_spans(x)
    lib_targets = set(x["unmodified_sequence"][~dec].tolist())
    acc_t = np.where((x["run_q_value"] <= 0.01) & ~dec)[0]
    rows = []  # accepted null winners
    for p, (t, d) in pairs.items():
        if t is None: continue
        ps = [q for q in SPLIT.split(x["proteins"][t]) if q]
        if not ps or not all(q.endswith("_p_target") for q in ps): continue
        w = d if (d is not None and s[d] >= s[t]) else t
        if x["run_q_value"][w] > 0.01: continue
        rows.append((p, t, d, w))
    ids = {int(x["precursor_id"][w]) for _, _, _, w in rows} | {int(x["precursor_id"][i]) for i in acc_t}
    fr = fragments(lib, ids)
    anchors = sorted((fr[int(x["precursor_id"][i])][0], i) for i in acc_t if int(x["precursor_id"][i]) in fr)
    amz = [m for m, _ in anchors]
    res = {"accepted_null_pairs": len(rows), "missing_fragments": sum(int(x["precursor_id"][w]) not in fr for *_, w in rows)}
    per = []
    mapped = 0
    for p, t, d, w in rows:
        apex = x["rt_apex"][w]; wid = int(x["precursor_id"][w]); isdec = bool(dec[w])
        par = parents(x["unmodified_sequence"][t], [q for q in SPLIT.split(x["proteins"][t]) if q])
        mapped += bool(par)
        arm_sp = [sp for pp in par for sp in spans.get(pp, [])]
        ref_sp = [sp for pp in par for sp in REFSPANS.get(pp, [])]
        co_arm = [sp for sp in arm_sp if sp[0] <= apex <= sp[1]]
        co_ref = [sp for sp in ref_sp if sp[0] <= apex <= sp[1]]
        in_lib = any(pp in lib_targets for pp in par)
        q4 = None
        if co_arm:
            cat = "A"
            sp = max(co_arm, key=lambda z: z[2])
            share = nshared(fr[wid][2], fr[sp[3]][2]) if wid in fr and sp[3] in fr else None
            q4 = {"parent_outscores": bool(sp[2] > s[w]), "shared_of_6": share}
        elif co_ref:
            cat = "B_lib" if in_lib else "B_cut"
        else:
            cat = "D"
            if wid in fr:
                mz, im, cf = fr[wid]
                lo, hi = bisect.bisect_left(amz, mz - MZ), bisect.bisect_right(amz, mz + MZ)
                tseq = x["unmodified_sequence"][t]
                for _, a in anchors[lo:hi]:
                    if a == w or x["unmodified_sequence"][a] == tseq: continue
                    if not (x["rt_start"][a] <= apex <= x["rt_end"][a]): continue
                    fa = fr[int(x["precursor_id"][a])]
                    if abs(fa[1] - im) > IM: continue
                    if nshared(cf, fa[2]) >= 1: cat = "C"; break
        per.append({"pid": wid, "target_pid": int(x["precursor_id"][t]), "decoy_won": isdec, "cat": cat, "q4": q4,
                    "tseq": x["sequence"][t], "charge": int(x["charge"][t]), "apex": float(apex)})
    res["parent_mapped"] = mapped
    # Q3 survey join
    if survey != "-":
        sv = {}
        for r in csv.DictReader(open(survey), delimiter="\t"):
            sv[r["traml_id"]] = (float(r["target_coelution"]), float(r["decoy_coelution"]))
        for e in per:
            v = sv.get(f"{e['tseq']}/{e['charge']}")
            e["survey_t_gt_d"] = None if v is None else bool(v[0] > v[1])
        res["survey_joined"] = sum(e["survey_t_gt_d"] is not None for e in per)
    # Q2 features join (winner window = same precursor, nearest apex)
    if feats != "-":
        want = {e["pid"]: e["apex"] for e in per}
        rd = pcsv.open_csv(feats, parse_options=pcsv.ParseOptions(delimiter="\t"),
                           read_options=pcsv.ReadOptions(block_size=1 << 26))
        names, got = None, {}
        for b in rd:
            if names is None:
                names = [c for c in b.schema.names if c not in ("candidate_row", "precursor_id", "is_decoy", "fold",
                                                                "score", "q_value", "network_score", "network_q_value")]
            m = pc.is_in(b.column("precursor_id"), value_set=pa.array(list(want), pa.int64()))
            b = b.filter(m)
            if not b.num_rows: continue
            pid = b.column("precursor_id").to_pylist(); ap = b.column("rt_apex").to_pylist()
            vals = np.column_stack([np.asarray(b.column(c).to_numpy(zero_copy_only=False), float) for c in names])
            for k in range(len(pid)):
                dd = abs(ap[k] - want[pid[k]])
                if pid[k] not in got or dd < got[pid[k]][0]: got[pid[k]] = (dd, vals[k])
        ok = 0
        for e in per:
            g = got.get(e["pid"])
            if g is not None and g[0] < 1e-3: e["features"] = g[1]; ok += 1
        res["features_joined"] = ok; res["feature_names"] = names
    if DRY:
        out = {k: v for k, v in res.items() if k != "feature_names"}
        print(arm, json.dumps(out), flush=True)
        return out
    # counts
    T = sum(not e["decoy_won"] for e in per); D = len(per) - T
    res["T"], res["D"] = T, D
    res["by_cat"] = {}
    for c in CATS:
        tc = sum(e["cat"] == c and not e["decoy_won"] for e in per); dc = sum(e["cat"] == c and e["decoy_won"] for e in per)
        res["by_cat"][c] = {"T": tc, "D": dc, "excess": tc - dc, "share_of_excess": (tc - dc) / (T - D) if T != D else None}
    a = [e["q4"] for e in per if e["cat"] == "A" and not e["decoy_won"] and e["q4"]]
    res["Q4_A_targets"] = {"n": len(a), "parent_outscores": sum(z["parent_outscores"] for z in a),
                           "shared_hist": {str(k): sum(z["shared_of_6"] == k for z in a) for k in range(7)}}
    if survey != "-":
        res["Q3"] = {}
        for grp, cs in (("all", CATS), ("A+B", ["A", "B_lib", "B_cut"]), ("C+D", ["C", "D"]), ("C", ["C"]), ("D", ["D"])):
            for side in ("T", "D"):
                v = [e["survey_t_gt_d"] for e in per if e["cat"] in cs and e["decoy_won"] == (side == "D")
                     and e["survey_t_gt_d"] is not None]
                res["Q3"][f"{grp}_{side}"] = {"n": len(v), "survey_target_higher": sum(v) / len(v) if v else None}
    if feats != "-":
        res["Q2"] = {}
        for grp, cs in (("A+B", ["A", "B_lib", "B_cut"]), ("C+D", ["C", "D"]), ("C", ["C"]), ("D", ["D"]), ("all", CATS)):
            ft = np.array([e["features"] for e in per if "features" in e and e["cat"] in cs and not e["decoy_won"]])
            fd = np.array([e["features"] for e in per if "features" in e and e["cat"] in cs and e["decoy_won"]])
            if len(ft) < 3 or len(fd) < 3:
                res["Q2"][grp] = {"n_T": len(ft), "n_D": len(fd)}; continue
            sd = np.sqrt((ft.var(0, ddof=1) + fd.var(0, ddof=1)) / 2)
            smd = np.where(sd > 0, (ft.mean(0) - fd.mean(0)) / np.where(sd > 0, sd, 1), 0.0)
            order = np.argsort(-np.abs(smd))
            res["Q2"][grp] = {"n_T": len(ft), "n_D": len(fd),
                              "top_smd": [[names[k], round(float(smd[k]), 3), round(float(ft[:, k].mean()), 4),
                                           round(float(fd[:, k].mean()), 4)] for k in order[:12]],
                              "smd_all": {names[k]: round(float(smd[k]), 3) for k in range(len(names))},
                              "i_nclaim_pos_T": float((ft[:, names.index("i_nclaim")] > 0).mean()),
                              "i_nclaim_pos_D": float((fd[:, names.index("i_nclaim")] > 0).mean()),
                              "i_explained_T": float(ft[:, names.index("i_explained")].mean()),
                              "i_explained_D": float(fd[:, names.index("i_explained")].mean())}
    out = {k: v for k, v in res.items() if k != "feature_names"}
    out["pairs"] = [{k: v for k, v in e.items() if k != "features"} for e in per]
    print(arm, json.dumps({k: out[k] for k in ("T", "D", "by_cat")}), flush=True)
    return out


out = sys.argv[1]
for spec in sys.argv[2:]:
    arm, path, lib, survey, feats = spec.split(":")
    R[arm] = run(arm, path, lib, survey, feats)
    json.dump(R, open(out, "w"), indent=1, default=str)
if not DRY:
    groups = {"V1 (FbV1 s0/s1/s2/g1)": ["FbV1", "FbV1_s1", "FbV1_s2", "FbV1_g1"], "V0 (FbP, FbP_g1)": ["FbP", "FbP_g1"],
              "V1+pool (FbPp)": ["FbPp"], "K = 0 (FbF)": ["FbF"]}
    R["pooled"] = {}
    for g, arms in groups.items():
        arms = [a for a in arms if a in R]
        if not arms: continue
        T = sum(R[a]["T"] for a in arms); D = sum(R[a]["D"] for a in arms)
        R["pooled"][g] = {"T": T, "D": D, "by_cat": {c: {"T": sum(R[a]["by_cat"][c]["T"] for a in arms),
                                                         "D": sum(R[a]["by_cat"][c]["D"] for a in arms)} for c in CATS}}
        for c in CATS:
            v = R["pooled"][g]["by_cat"][c]; v["excess"] = v["T"] - v["D"]
            v["share_of_excess"] = v["excess"] / (T - D) if T != D else None
    json.dump(R, open(out, "w"), indent=1, default=str)
    for g, v in R["pooled"].items():
        print(g, "T", v["T"], "D", v["D"], {c: (z["T"], z["D"], z["excess"], z["share_of_excess"] and round(z["share_of_excess"], 2))
                                            for c, z in v["by_cat"].items()})
