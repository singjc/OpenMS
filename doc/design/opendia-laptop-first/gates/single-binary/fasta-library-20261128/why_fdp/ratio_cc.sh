#!/bin/bash
# Run on codon as a Slurm job: entrapment/human target ratio of the candidates (r0) and final libraries (r_lib) of the
# PLAN_CLAIM_COMPETITION arms (FDRBench: _p_target suffix; custom Ce/CcPep: ENTRAP_ prefix).
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -t 1:00:00
set -euo pipefail
cd /nfs/research/juan/singj/opendia/fasta-20261128
cat > ${TMPDIR:-/tmp}/ratio_$$.py <<'PY'
import json, glob
import numpy as np, pyarrow as pa, pyarrow.compute as pc, pyarrow.parquet as pq
def library_ratio(path, custom):
    t = pq.read_table(path, columns=["precursor_index", "decoy", "proteins"]).combine_chunks()
    idx = t.column("precursor_index").to_numpy()
    first = np.ones(len(idx), bool); first[1:] = idx[1:] != idx[:-1]
    assert len(np.unique(idx)) == first.sum()
    keep = first & ~t.column("decoy").to_numpy(zero_copy_only=False)
    prot = t.column("proteins").filter(pa.array(keep)).combine_chunks()
    flat = (pc.starts_with(pc.list_flatten(prot), "ENTRAP_") if custom else pc.ends_with(pc.list_flatten(prot), "_p_target")).to_numpy(zero_copy_only=False)
    parent = pc.list_parent_indices(prot).to_numpy(); n = len(prot)
    human_any = np.bincount(parent, weights=~flat, minlength=n) > 0
    has = np.bincount(parent, minlength=n) > 0
    e = int((has & ~human_any).sum())
    return {"entrapment_targets": e, "human_targets": n - e, "ratio": e / (n - e)}
R = {}
import sys
if len(sys.argv) > 1 and sys.argv[1] == "check":
    # Self-check of both label branches on the existing reference caches (recorded: Fe 0.58938, Pe' 0.59024, FbV1 0.99940).
    for p, custom in (("cache_v2/Fe/human_entrapment.candidates.parquet", True), ("cache_v4/Pep/human_entrapment.candidates.parquet", True),
                      ("cache_fb/FbV1/human_fdrbench.candidates.parquet", False)):
        print(p, library_ratio(p, custom)["ratio"])
    sys.exit(0)
for a in ("CcV1", "CcV1_s1", "CcP", "CcPp", "CcF", "Ce", "CcPep"):
    for kind in ("candidates", "library"):
        stem = "human_entrapment" if a in ("Ce", "CcPep") else "human_fdrbench"
        f = f"cache_cc/{a}/{stem}.{kind}.parquet"
        if glob.glob(f): R[f"{kind}_{a}"] = library_ratio(f, a in ("Ce", "CcPep"))
print(json.dumps(R, indent=1))
PY
$HOME/.local/bin/uv run -q --with pyarrow --with numpy python ${TMPDIR:-/tmp}/ratio_$$.py ${1:-}
rm -f ${TMPDIR:-/tmp}/ratio_$$.py
