#!/usr/bin/env python3
"""PLAN_NEIGHBOUR_COMPETITION: label-blind neighbour competition applied after scoring, then plain TDC again.

Rule (applied identically to target and decoy rows; no entrapment label is read until evaluation):
  anchors   = accepted target pair winners (run q <= 0.01 in the arm's own output).
  candidate = any pair winner with run q <= 0.05.
  A candidate c is *explained* by an anchor a when all hold:
    score(a) > score(c);
    unmodified_sequence(a) != the target sequence of c's pair (a charge sibling of the pair's target is not a neighbour);
    rt_start(a) <= rt_apex(c) <= rt_end(a);
    |precursor_mz(a) - precursor_mz(c)| <= MZ (25 Th, the diaPASEF window width);
    |precursor_im(a) - precursor_im(c)| <= IM (0.05, library 1/K0);
    share(c, a) >= SHARE, share = fraction of c's top-6 search-library fragments within PPM (10 ppm) of any of a's
    search-library fragments.
  An explained candidate gets score -inf: its pair partner, if any, becomes the pair winner with its own score.
  q is then recomputed by plain TDC, q = cummin((D + 1) / T), on all pair winners ranked by discriminant_score.
Evaluation: N, N_e, FDP_c = N_e (1 + 1/r0) / N at q <= 0.01, T_e / D_e, before and after.
Usage: nc_remedy.py OUT.json label:arm:precursors.tsv:library.search.parquet:r0 ...
"""
import bisect, json, math, re, sys
import numpy as np
import pyarrow as pa, pyarrow.compute as pc, pyarrow.csv as pcsv, pyarrow.parquet as pq

MZ, IM, PPM, SHARE_LEVELS, PRIMARY = 25.0, 0.05, 10e-6, (2 / 6, 0.5, 4 / 6), 0.5
SPLIT = re.compile(r"[;,]")
DRY = bool(int(__import__("os").environ.get("NC_DRY", "0")))


def tdc(score, decoy):
    o = np.lexsort((~decoy, -score))  # descending score, ties: decoy first
    d = decoy[o]
    T, D = np.cumsum(~d), np.cumsum(d)
    f = np.minimum(1.0, (D + 1) / np.maximum(T, 1))
    q = np.empty_like(f)
    q[o] = np.minimum.accumulate(f[::-1])[::-1]
    return q


def share(cf, af):
    """fraction of candidate fragments cf within PPM of any anchor fragment (af sorted)."""
    if not len(cf): return 0.0
    i = np.searchsorted(af, cf)
    lo = np.abs(af[np.clip(i - 1, 0, len(af) - 1)] - cf)
    hi = np.abs(af[np.clip(i, 0, len(af) - 1)] - cf)
    return float(np.mean(np.minimum(lo, hi) <= PPM * cf))


def load(path):
    cols = ["precursor_id", "unmodified_sequence", "proteins", "rt_apex", "rt_start", "rt_end", "discriminant_score",
            "run_q_value", "decoy"]
    t = pcsv.read_csv(path, parse_options=pcsv.ParseOptions(delimiter="\t"),
                      convert_options=pcsv.ConvertOptions(include_columns=cols, strings_can_be_null=False,
                                                          column_types={"precursor_id": pa.int64(), "decoy": pa.string(),
                                                                        "proteins": pa.string()}))
    return {c: t.column(c).to_numpy(zero_copy_only=False) for c in cols}


def fragments(lib, ids):
    want = pa.array([str(i) for i in ids])
    f = pq.ParquetFile(lib)
    out = {}
    for g in range(f.num_row_groups):
        t = f.read_row_group(g, columns=["precursor_id", "precursor_mz", "precursor_im", "product_mz", "fragment_rank"])
        t = t.filter(pc.is_in(t.column("precursor_id"), value_set=want))
        for pid, mz, im, fm, rk in zip(*(t.column(c).to_pylist() for c in t.column_names)):
            e = out.setdefault(int(pid), [mz, im, []])
            e[2].append((rk, fm))
    for e in out.values():
        e[2].sort()
    return out


def evaluate(win_score, win_decoy, win_ent, r0):
    q = tdc(win_score, win_decoy)
    acc = q <= 0.01
    T = int((acc & ~win_decoy).sum()); Te = int((acc & ~win_decoy & win_ent).sum())
    De = int((acc & win_decoy & win_ent).sum()); D = int((acc & win_decoy).sum())
    k = 1 + 1 / r0
    return {"N": T, "N_e": Te, "N_human": T - Te, "D": D, "D_e": De, "T_e_over_D_e": Te / De if De else None,
            "FDP_c": Te * k / T if T else None,
            "FDP_c_ci95": [max(Te - 1.96 * math.sqrt(Te), 0) * k / T, (Te + 1.96 * math.sqrt(Te)) * k / T] if T else None,
            "FDP_c_var": Te * (k / T) ** 2 if T else None}


def run(mode, arm, path, lib, r0):
    is_e = (lambda p: p.startswith("ENTRAP_")) if mode == "entrap" else (lambda p: p.endswith("_p_target"))
    x = load(path)
    n = len(x["precursor_id"])
    dec = x["decoy"] == "true"
    ent = np.array([bool(ps) and all(is_e(p[6:] if p.startswith("DECOY_") else p) for p in ps)
                    for ps in ([p for p in SPLIT.split(s) if p] for s in x["proteins"])])
    s = x["discriminant_score"].astype(float)
    pair = x["precursor_id"] // 2
    # pair members: best row per (pair, side)
    rows = {}
    for i in range(n):
        k = (int(pair[i]), bool(dec[i]))
        if k not in rows or s[i] > s[rows[k]]: rows[k] = i
    pairs = {}
    for (p, d), i in rows.items(): pairs.setdefault(p, [None, None])[int(d)] = i
    keys = np.fromiter(pairs.keys(), dtype=np.int64)
    tgt = np.array([pairs[p][0] if pairs[p][0] is not None else -1 for p in keys])
    dcy = np.array([pairs[p][1] if pairs[p][1] is not None else -1 for p in keys])
    st = np.where(tgt >= 0, s[np.maximum(tgt, 0)], -np.inf); sd = np.where(dcy >= 0, s[np.maximum(dcy, 0)], -np.inf)
    win = np.where(sd >= st, dcy, tgt)  # ties to the decoy
    tseq = np.where(tgt >= 0, x["unmodified_sequence"][np.maximum(tgt, 0)], "")
    pair_ent = np.array([(ent[t] if t >= 0 else False) or (ent[d] if d >= 0 else False) for t, d in zip(tgt, dcy)])
    q0 = tdc(s[win], dec[win])
    if DRY: pair_ent = np.zeros(len(keys), bool)  # never look at labels in the dry run
    res = {"label_mode": mode, "r0": r0, "before": evaluate(s[win], dec[win], pair_ent, r0),
           "tool_N": int(((x["run_q_value"] <= 0.01) & ~dec).sum())}
    anc = np.where((q0 <= 0.01) & ~dec[win])[0]  # pair indices of anchors
    cand = np.where(q0 <= 0.05)[0]
    ids = set(x["precursor_id"][win[cand]].tolist())
    fr = fragments(lib, ids)
    am = [(fr[x["precursor_id"][win[j]]][0], j) for j in anc if x["precursor_id"][win[j]] in fr]
    am.sort()
    amz = [m for m, _ in am]
    afr = {j: np.sort(np.array([m for _, m in fr[x["precursor_id"][win[j]]][2]])) for _, j in am}
    best = np.zeros(len(keys))
    missing = 0
    for j in cand:
        pid = x["precursor_id"][win[j]]
        if pid not in fr: missing += 1; continue
        mz, im, fl = fr[pid]
        cf = np.array([m for r, m in fl if r <= 6])
        r_apex, sc = x["rt_apex"][win[j]], s[win[j]]
        lo, hi = bisect.bisect_left(amz, mz - MZ), bisect.bisect_right(amz, mz + MZ)
        b = 0.0
        for m_, a in am[lo:hi]:
            ra = win[a]
            if a == j or s[ra] <= sc or x["unmodified_sequence"][ra] == tseq[j]: continue
            if not (x["rt_start"][ra] <= r_apex <= x["rt_end"][ra]): continue
            if abs(fr[x["precursor_id"][ra]][1] - im) > IM: continue
            b = max(b, share(cf, afr[a]))
        best[j] = b
    res["candidates"] = int(len(cand)); res["anchors"] = int(len(anc)); res["missing_fragments"] = missing
    if DRY:  # plumbing check only: no removal results
        res = {k: res[k] for k in ("tool_N", "candidates", "anchors", "missing_fragments")} | {"recomputed_N": res["before"]["N"]}
        print(arm, json.dumps(res), flush=True)
        return res
    for lev in SHARE_LEVELS:
        rm = best >= lev - 1e-9
        ns = np.where(rm, np.where(dec[win], st, sd), np.maximum(st, sd))
        nd = np.where(rm, ~dec[win], dec[win])  # partner promoted (score -inf for the removed member)
        ok = np.isfinite(ns)
        r = evaluate(ns[ok], nd[ok], pair_ent[ok], r0)
        acc0 = q0 <= 0.01
        r["removed_accepted_targets"] = int((rm & acc0 & ~dec[win]).sum())
        r["removed_accepted_targets_null"] = int((rm & acc0 & ~dec[win] & pair_ent).sum())
        r["removed_accepted_decoys"] = int((rm & acc0 & dec[win]).sum())
        r["removed_accepted_decoys_null"] = int((rm & acc0 & dec[win] & pair_ent).sum())
        r["removed_q05_targets"] = int((rm & ~dec[win]).sum()); r["removed_q05_decoys"] = int((rm & dec[win]).sum())
        res[f"share_ge_{lev:.3f}"] = r
    print(arm, json.dumps({k: (v if not isinstance(v, dict) else {kk: v[kk] for kk in ("N", "N_e", "D_e", "FDP_c")})
                           for k, v in res.items()}), flush=True)
    return res


if __name__ == "__main__":
    out = sys.argv[1]
    R = {"params": {"MZ": MZ, "IM": IM, "PPM": PPM, "primary_share": PRIMARY}}
    for spec in sys.argv[2:]:
        mode, arm, path, lib, r0 = spec.split(":")
        R[arm] = run(mode, arm, path, lib, float(r0))
        json.dump(R, open(out, "w"), indent=1)
