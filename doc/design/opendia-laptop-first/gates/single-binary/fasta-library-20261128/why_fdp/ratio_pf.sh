#!/bin/bash
# Run on codon as a Slurm job: entrapment/human target ratio of the candidates (r0) and final libraries (r_lib) of the
# FDRBench arms. Same function as ratio_v2.sh, with entrapment = every protein ends with _p_target.
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -t 1:00:00
set -euo pipefail
cd /nfs/research/juan/singj/opendia/fasta-20261128
cat > ${TMPDIR:-/tmp}/ratio_$$.py <<'PY'
import json, glob
import numpy as np, pyarrow as pa, pyarrow.compute as pc, pyarrow.parquet as pq
def library_ratio(path):
    t = pq.read_table(path, columns=["precursor_index", "decoy", "proteins"]).combine_chunks()
    idx = t.column("precursor_index").to_numpy()
    first = np.ones(len(idx), bool); first[1:] = idx[1:] != idx[:-1]
    assert len(np.unique(idx)) == first.sum()
    keep = first & ~t.column("decoy").to_numpy(zero_copy_only=False)
    prot = t.column("proteins").filter(pa.array(keep)).combine_chunks()
    flat = pc.ends_with(pc.list_flatten(prot), "_p_target").to_numpy(zero_copy_only=False)
    parent = pc.list_parent_indices(prot).to_numpy(); n = len(prot)
    human_any = np.bincount(parent, weights=~flat, minlength=n) > 0
    has = np.bincount(parent, minlength=n) > 0
    e = int((has & ~human_any).sum())
    return {"entrapment_targets": e, "human_targets": n - e, "ratio": e / (n - e)}
R = {}
for a in ("FbV1_g1", "FbP_g1", "FbV1_s1", "FbV1_s2"):
    for kind in ("candidates", "library"):
        f = f"cache_pf/{a}/human_fdrbench.{kind}.parquet"
        if glob.glob(f): R[f"{kind}_{a}"] = library_ratio(f)
print(json.dumps(R, indent=1))
PY
uv run -q --with pyarrow --with numpy python ${TMPDIR:-/tmp}/ratio_$$.py
rm -f ${TMPDIR:-/tmp}/ratio_$$.py
