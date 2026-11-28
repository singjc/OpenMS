#!/bin/bash
# Run on codon (pst_prd): entrapment/human target ratio of the Pe and Fe candidate sets (r0) and of the final libraries
# (sensitivity), with the same library_ratio as eval_fasta.py.
set -euo pipefail
cd /nfs/research/juan/singj/opendia/fasta-20261128
cat > ${TMPDIR:-/tmp}/ratio_$$.py <<'EOF'
import json, sys, glob
import numpy as np, pyarrow as pa, pyarrow.compute as pc, pyarrow.parquet as pq
def library_ratio(path):
    t = pq.read_table(path, columns=["precursor_index", "decoy", "proteins"]).combine_chunks()
    idx = t.column("precursor_index").to_numpy()
    first = np.ones(len(idx), bool); first[1:] = idx[1:] != idx[:-1]
    assert len(np.unique(idx)) == first.sum()
    keep = first & ~t.column("decoy").to_numpy(zero_copy_only=False)
    prot = t.column("proteins").filter(pa.array(keep)).combine_chunks()
    flat = pc.starts_with(pc.list_flatten(prot), "ENTRAP_").to_numpy(zero_copy_only=False)
    parent = pc.list_parent_indices(prot).to_numpy(); n = len(prot)
    human_any = np.bincount(parent, weights=~flat, minlength=n) > 0
    has = np.bincount(parent, minlength=n) > 0
    e = int((has & ~human_any).sum())
    return {"entrapment_targets": e, "human_targets": n - e, "ratio": e / (n - e)}
R = {}
for a in ("Pe", "Fe"):
    for kind in ("candidates", "library"):
        f = f"cache_v2/{a}/human_entrapment.{kind}.parquet"
        if glob.glob(f): R[f"{kind}_{a}"] = library_ratio(f)
print(json.dumps(R, indent=1))
EOF
uv run -q --with pyarrow --with numpy python ${TMPDIR:-/tmp}/ratio_$$.py
rm -f ${TMPDIR:-/tmp}/ratio_$$.py
