# Recall at K=800k under the exact C++ pooledStatistic rule (absent pairs excluded), vs the evaluator's rule.
import gzip, csv, collections, sys
import numpy as np, pyarrow as pa, pyarrow.parquet as pq
run, tarm = sys.argv[1], sys.argv[2]
K = 800_000
traml, T, D = [], [], []
with gzip.open(f"survey_v3/V1-{run}/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.survey.tsv.gz", "rt") as f:
    next(f)
    for line in f:
        r = line.rstrip("\n").split("\t"); traml.append(r[1]); T.append(float(r[2])); D.append(float(r[3]))
S = np.maximum(np.array(T), np.array(D)); n = len(S)
tb = pq.read_table("design_D/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.candidates.parquet",
                   columns=["precursor_index", "precursor_id", "decoy", "proteins"]).combine_chunks()
idx = tb.column("precursor_index").to_numpy(); first = np.ones(len(idx), bool); first[1:] = idx[1:] != idx[:-1]
sub = tb.filter(pa.array(first & ~tb.column("decoy").to_numpy(zero_copy_only=False)))
prot = [()] * n
for p, ps in zip(sub.column("precursor_id").to_pylist(), sub.column("proteins").to_pylist()): prot[int(p) // 2] = ps
present = np.array([bool(t) for t in traml])
seq = [t.rsplit("/", 1)[0] if t else "" for t in traml]
sib = collections.defaultdict(float); per = collections.defaultdict(list)
for k in range(n):
    if not present[k]: continue
    sib[seq[k]] = max(sib[seq[k]], S[k])
    for p in prot[k]: per[p].append(S[k])
ev = {p: float(np.mean(sorted(v, reverse=True)[:3])) for p, v in per.items()}
pool = np.array([(sib[seq[k]] + 0.5 * max((ev[p] for p in prot[k]), default=0.0)) if present[k] else 0.0 for k in range(n)])
tie = None
tr = np.zeros(n, bool)
for r in csv.DictReader(gzip.open(f"v2/slim/runs_v2/{tarm}/precursors.slim.tsv.gz", "rt"), delimiter="\t"):
    if r["decoy"] == "false" and float(r["run_q_value"]) <= 0.01: tr[int(r["precursor_id"]) // 2] = True
order = np.lexsort((np.arange(n), -pool)); kept = np.zeros(n, bool); kept[order[:K]] = True; kept &= pool > 0
print(run, "cpp-rule recall", round(tr[kept].sum() / tr.sum(), 4), "kept", int(kept.sum()), "truth", int(tr.sum()))
