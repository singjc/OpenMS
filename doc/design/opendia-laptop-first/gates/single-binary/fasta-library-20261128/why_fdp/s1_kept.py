#!/usr/bin/env python3
"""Write the pair indices in each FDRBench prefilter arm's final library (run on codon)."""
import pyarrow.compute as pc, pyarrow.parquet as pq

W = "/nfs/research/juan/singj/opendia/fasta-20261128"
for arm in ("FbP", "FbV1", "FbPp"):
    ids = pc.unique(pq.read_table(f"{W}/cache_fb/{arm}/human_fdrbench.library.parquet", columns=["precursor_id"]).column("precursor_id"))
    pairs = sorted({int(x) // 2 for x in ids.to_pylist()})
    with open(f"{W}/why_fdp/kept_{arm}.txt", "w") as f:
        f.write("\n".join(map(str, pairs)) + "\n")
    print(arm, len(pairs))
