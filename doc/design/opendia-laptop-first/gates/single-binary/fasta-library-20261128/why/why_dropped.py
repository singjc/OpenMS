#!/usr/bin/env python3
"""Why does the survey prefilter drop peptides that are in the data?

Uses, read only:
  why/survey_P.tsv.gz                   codon P survey (pair, traml_id, target/decoy coelution, S = max)
  v2/slim/runs_v2/{P,F}/precursors.slim.tsv.gz
  design_D/*.candidates.parquet         local candidate library: precursor_id -> traml_id, sequence, charge, mz, rt, im
  entity-inference out_pasef_owned legacy/precursors.tsv  (empirical IDs, read only)

Groups of true positives (target IDs at q <= 0.01 by F, the full library):
  kept    : in P's top-K by S
  dropped : not kept
For each group: the S distribution, the share with S = 0 (no survey candidate), target-vs-decoy (T > D) wins,
charge, length, missed cleavages, oxidation, predicted RT position, and whether the empirical library has it.
"""
import csv, gzip, json, sys, collections
import numpy as np
import pyarrow.parquet as pq

S = "/home/singjc/.jcode/scratch/fasta-library-20261128"
K = 800_000


def fnv1a(s):
    h = 0xcbf29ce484222325
    for b in s.encode():
        h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h


pairs = []
with gzip.open(f"{S}/why/survey_P.tsv.gz", "rt") as f:
    r = csv.reader(f, delimiter="\t"); next(r)
    for row in r:
        pairs.append((int(row[0]), row[1], float(row[2]), float(row[3]), float(row[4])))
stat = {p[1]: p for p in pairs}
ranked = sorted((p for p in pairs if p[4] > 0), key=lambda p: (-p[4], fnv1a(p[1]), p[0]))
rank = {p[1]: i for i, p in enumerate(ranked)}
kept = {p[1] for p in ranked[:K]}
s_cut = ranked[K - 1][4]

# precursor_id -> traml (target only); IDs match because candidate enumeration is deterministic (same key)
t = pq.read_table(f"{S}/design_D/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.candidates.parquet",
                  columns=["precursor_id", "traml_id", "decoy", "sequence", "precursor_charge", "library_rt", "precursor_im"])
pid = t.column("precursor_id").to_numpy(zero_copy_only=False)
first = np.ones(len(pid), bool); first[1:] = pid[1:] != pid[:-1]
info = {}
for p, tr, d, seq, z, rt in zip(pid[first], t.column("traml_id").to_numpy(zero_copy_only=False)[first],
                               t.column("decoy").to_numpy(zero_copy_only=False)[first],
                               t.column("sequence").to_numpy(zero_copy_only=False)[first],
                               t.column("precursor_charge").to_numpy()[first], t.column("library_rt").to_numpy()[first]):
    if not d:
        info[p] = (tr, seq, int(z), float(rt))
del t


def accepted(arm):
    out = {}
    for r in csv.DictReader(gzip.open(f"{S}/v2/slim/runs_v2/{arm}/precursors.slim.tsv.gz", "rt"), delimiter="\t"):
        if r["decoy"] == "false" and float(r["run_q_value"]) <= 0.01:
            out[r["precursor_id"]] = float(r["discriminant_score"])
    return out


F = accepted("F"); P = accepted("P")
emp = set()
for r in csv.DictReader(open("/home/singjc/.jcode/scratch/entity-inference-20261127/out_pasef_owned/legacy/precursors.tsv"), delimiter="\t"):
    emp.add(r["sequence"] + "/" + r["charge"])

rts = np.array([v[3] for v in info.values()])
groups = collections.defaultdict(list)
for p, score in F.items():
    tr, seq, z, rt = info[p]
    groups["kept" if tr in kept else "dropped"].append((p, tr, seq, z, rt, score))


def describe(rows):
    s = np.array([stat[r[1]][4] for r in rows]); tc = np.array([stat[r[1]][2] for r in rows]); dc = np.array([stat[r[1]][3] for r in rows])
    rk = np.array([rank.get(r[1], len(ranked)) for r in rows])
    import re
    plain = [re.sub(r"\(UniMod:\d+\)", "", r[2]) for r in rows]
    mc = [sum(1 for i, a in enumerate(x[:-1]) if a in "KR" and x[i + 1] != "P") for x in plain]
    rtq = np.searchsorted(np.sort(rts), [r[4] for r in rows]) / len(rts)
    fscore = np.array([r[5] for r in rows])
    return {
        "n": len(rows),
        "S_median": float(np.median(s)), "S_quartiles": [float(np.quantile(s, .25)), float(np.quantile(s, .75))],
        "S_zero_share": float((s == 0).mean()),
        "target_beats_decoy_share": float((tc > dc).mean()),
        "decoy_sets_S_share": float((dc > tc).mean()),
        "target_coelution_median": float(np.median(tc)), "decoy_coelution_median": float(np.median(dc)),
        "rank_median": float(np.median(rk)), "rank_quartiles": [float(np.quantile(rk, .25)), float(np.quantile(rk, .75))],
        "charge": dict(collections.Counter(r[3] for r in rows)),
        "length_median": float(np.median([len(x) for x in plain])),
        "missed_cleavage": dict(collections.Counter(mc)),
        "oxidation_share": float(np.mean(["UniMod:35" in r[2] for r in rows])),
        "predicted_rt_quantile_median": float(np.median(rtq)),
        "predicted_rt_edge_share_(<5%_or_>95%)": float(np.mean((rtq < .05) | (rtq > .95))),
        "in_empirical_library_ids": float(np.mean([r[1] in emp for r in rows])),
        "F_discriminant_median": float(np.median(fscore)),
    }


background = [stat[p[1]] for p in pairs]
bs = np.array([p[4] for p in background]); btc = np.array([p[2] for p in background]); bdc = np.array([p[3] for p in background])
res = {
    "K": K, "S_at_cut": s_cut, "pairs": len(pairs), "pairs_S_gt_0": len(ranked),
    "all_pairs": {"S_median": float(np.median(bs)), "S_zero_share": float((bs == 0).mean()),
                  "target_beats_decoy_share": float((btc > bdc).mean()), "decoy_sets_S_share": float((bdc > btc).mean())},
    "F_true_ids_kept_by_prefilter": describe(groups["kept"]),
    "F_true_ids_dropped_by_prefilter": describe(groups["dropped"]),
    "S_histogram_edges": [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6],
}
for name, rows in (("kept", groups["kept"]), ("dropped", groups["dropped"])):
    res[f"S_histogram_{name}"] = np.histogram([stat[r[1]][4] for r in rows], bins=res["S_histogram_edges"])[0].tolist()
res["S_histogram_all_pairs"] = np.histogram(bs, bins=res["S_histogram_edges"])[0].tolist()
# How large would K need to be to keep 95 % of F's IDs?
fr = np.sort([rank.get(info[p][0], len(ranked)) for p in F])
res["K_needed_for_F_recall"] = {str(q): int(fr[int(q * len(fr)) - 1]) + 1 for q in (0.5, 0.8, 0.9, 0.95)}
res["F_ids_with_S_0"] = int(sum(1 for p in F if stat[info[p][0]][4] == 0))
json.dump(res, open(f"{S}/why/why.json", "w"), indent=1)
print(json.dumps(res, indent=1))
