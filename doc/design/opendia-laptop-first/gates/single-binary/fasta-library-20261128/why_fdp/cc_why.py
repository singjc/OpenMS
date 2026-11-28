#!/usr/bin/env python3
"""Design evidence for PLAN_CLAIM_COMPETITION (descriptive, D only): for accepted null pair winners in one arm, the
best co-eluting accepted target anchor (same constraints as nc_remedy.py), how many of the winner's top-6 search
fragments it shares, which fragments (type+ordinal), the apex difference and the anchor's library intensity rank on
the shared fragment. Target winners vs decoy winners. Usage: cc_why.py arm precursors.tsv library.search.parquet"""
import bisect, json, re, sys
from collections import Counter
import numpy as np, pyarrow as pa, pyarrow.compute as pc, pyarrow.csv as pcsv, pyarrow.parquet as pq
MZ, IM, PPM = 25.0, 0.05, 10e-6
SPLIT = re.compile(r"[;,]")
arm, path, lib = sys.argv[1:4]
cols = ["precursor_id", "unmodified_sequence", "proteins", "rt_apex", "rt_start", "rt_end", "discriminant_score", "run_q_value", "decoy"]
t = pcsv.read_csv(path, parse_options=pcsv.ParseOptions(delimiter="\t"), convert_options=pcsv.ConvertOptions(include_columns=cols,
    column_types={"precursor_id": pa.int64(), "decoy": pa.string(), "proteins": pa.string()}))
x = {c: t.column(c).to_numpy(zero_copy_only=False) for c in cols}
dec = x["decoy"] == "true"; s = x["discriminant_score"]
best = {}
for i in range(len(s)):
    k = (int(x["precursor_id"][i]) // 2, bool(dec[i]))
    if k not in best or s[i] > s[best[k]]: best[k] = i
pairs = {}
for (p, d), i in best.items(): pairs.setdefault(p, [None, None])[int(d)] = i
acc_t = np.where((x["run_q_value"] <= 0.01) & ~dec)[0]
rows = []
for p, (ti, di) in pairs.items():
    if ti is None: continue
    ps = [q for q in SPLIT.split(x["proteins"][ti]) if q]
    if not ps or not all(q.endswith("_p_target") for q in ps): continue
    w = di if (di is not None and s[di] >= s[ti]) else ti
    if x["run_q_value"][w] <= 0.01: rows.append((ti, w))
ids = {int(x["precursor_id"][w]) for _, w in rows} | {int(x["precursor_id"][i]) for i in acc_t}
want = pa.array([str(i) for i in ids]); f = pq.ParquetFile(lib); fr = {}
for g in range(f.num_row_groups):
    b = f.read_row_group(g, columns=["precursor_id", "precursor_mz", "precursor_im", "product_mz", "fragment_rank", "fragment_type", "fragment_ordinal", "library_intensity"])
    b = b.filter(pc.is_in(b.column("precursor_id"), value_set=want))
    for pid, mz, im, fm, rk, ft, fo, li in zip(*(b.column(c).to_pylist() for c in b.column_names)):
        fr.setdefault(int(pid), [mz, im, []])[2].append((rk, fm, f"{ft}{fo}", li))
anc = sorted((fr[int(x["precursor_id"][i])][0], i) for i in acc_t if int(x["precursor_id"][i]) in fr)
amz = [m for m, _ in anc]
out = {"T": Counter(), "D": Counter()}; types = {"T": Counter(), "D": Counter()}; dap = {"T": [], "D": []}; arank = {"T": Counter(), "D": Counter()}
nanch = {"T": [], "D": []}
for ti, w in rows:
    side = "D" if dec[w] else "T"
    wid = int(x["precursor_id"][w]); mz, im, fl = fr[wid]
    cf = [(m, ty) for r, m, ty, _ in fl if r <= 6]
    lo, hi = bisect.bisect_left(amz, mz - MZ), bisect.bisect_right(amz, mz + MZ)
    bestn, bestinfo, n_co = 0, None, 0
    for _, a in anc[lo:hi]:
        if a == w or x["unmodified_sequence"][a] == x["unmodified_sequence"][ti]: continue
        if not (x["rt_start"][a] <= x["rt_apex"][w] <= x["rt_end"][a]): continue
        fa = fr[int(x["precursor_id"][a])]
        if abs(fa[1] - im) > IM: continue
        n_co += 1
        af = [(m, r) for r, m, _, _ in fa[2] if r <= 6]
        sh = [(ty, min((r for am, r in af if abs(am - m) <= PPM * m), default=None)) for m, ty in cf]
        sh = [(ty, r) for ty, r in sh if r is not None]
        if len(sh) > bestn or (len(sh) == bestn and bestinfo and s[a] > bestinfo[2]):
            bestn, bestinfo = len(sh), (sh, abs(x["rt_apex"][a] - x["rt_apex"][w]), s[a])
    out[side][bestn] += 1; nanch[side].append(n_co)
    if bestn >= 3:
        for ty, r in bestinfo[0]: types[side][ty] += 1; arank[side][r] += 1
        dap[side].append(bestinfo[1])
for side in "TD":
    n = sum(out[side].values())
    print(f"{arm} {side} n={n} co-eluting anchors median={np.median(nanch[side]):.0f} best shared-of-6 hist={dict(sorted(out[side].items()))}")
    if dap[side]:
        print(f"   share>=3: |d apex| median {np.median(dap[side]):.2f} s, p75 {np.percentile(dap[side], 75):.2f}; shared fragment types {types[side].most_common(10)}; anchor rank of shared {dict(sorted(arank[side].items()))}")
