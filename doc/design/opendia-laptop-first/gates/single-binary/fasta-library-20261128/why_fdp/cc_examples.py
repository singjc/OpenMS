import bisect, re, sys
import numpy as np, pyarrow as pa, pyarrow.compute as pc, pyarrow.csv as pcsv, pyarrow.parquet as pq
exec(open("cc_why.py").read().split("out = {")[0])  # reuse loading (arm, pairs, rows, fr, anc)
seq = pcsv.read_csv(path, parse_options=pcsv.ParseOptions(delimiter="\t"), convert_options=pcsv.ConvertOptions(include_columns=["sequence"])).column("sequence").to_pylist()
shown = {"T": 0, "D": 0}
for ti, w in rows:
    side = "D" if dec[w] else "T"
    if shown[side] >= 8: continue
    wid = int(x["precursor_id"][w]); mz, im, fl = fr[wid]
    cf = [(m, ty) for r, m, ty, _ in fl if r <= 6]
    lo, hi = bisect.bisect_left(amz, mz - MZ), bisect.bisect_right(amz, mz + MZ)
    for _, a in anc[lo:hi]:
        if a == w or x["unmodified_sequence"][a] == x["unmodified_sequence"][ti]: continue
        if not (x["rt_start"][a] <= x["rt_apex"][w] <= x["rt_end"][a]): continue
        fa = fr[int(x["precursor_id"][a])]
        if abs(fa[1] - im) > IM: continue
        sh = [ty for m, ty in cf if any(abs(am - m) <= PPM * m for r, am, _, _ in fa[2] if r <= 6)]
        if len(sh) >= 3:
            print(f"{side} winner {seq[w]:28s} mz {mz:.3f} | anchor {seq[a]:28s} mz {fa[0]:.3f} apex {x['rt_apex'][w]:.2f}/{x['rt_apex'][a]:.2f} shared {sh}")
            shown[side] += 1; break
