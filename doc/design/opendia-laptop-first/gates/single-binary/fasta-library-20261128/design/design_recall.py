#!/usr/bin/env python3
"""PLAN_FASTA design rule: recall(K) of the empirical D IDs among the kept survey pairs.

E_c = empirical target IDs at q <= 0.01 (legacy/precursors.tsv, read only) whose (UniMod sequence, charge) is a
candidate target. recall(K) = |E_c among the top-K pairs| / |E_c|. Pairs are ranked as selectPairs ranks them:
S descending, FNV-1a 64 of the target TRAML_ID ascending, pair index ascending; S = 0 is never kept.
K* = the smallest K in {100k, 200k, 300k, 500k, 800k} with recall >= 0.95, else 800k.
"""
import csv
import json
import sys

survey, empirical, out = sys.argv[1], sys.argv[2], sys.argv[3]


def fnv1a(s):
    h = 0xcbf29ce484222325
    for b in s.encode():
        h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h


pairs = []
with open(survey) as f:
    r = csv.reader(f, delimiter="\t")
    next(r)
    for row in r:
        pairs.append((int(row[0]), row[1], float(row[2]), float(row[3]), float(row[4])))
rank_of = {}
ranked = sorted((p for p in pairs if p[4] > 0), key=lambda p: (-p[4], fnv1a(p[1]), p[0]))
for i, p in enumerate(ranked):
    rank_of[p[1]] = i
candidates = {p[1] for p in pairs}

ids = set()
with open(empirical) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        ids.add(row["sequence"] + "/" + row["charge"])
e_c = sorted(i for i in ids if i in candidates)
grid = [100_000, 200_000, 300_000, 500_000, 800_000]
recall = {k: sum(1 for i in e_c if rank_of.get(i, 10**12) < k) / len(e_c) for k in grid}
k_star = next((k for k in grid if recall[k] >= 0.95), 800_000)
target_wins = sum(1 for p in pairs if p[2] > p[3])
decoy_wins = sum(1 for p in pairs if p[3] > p[2])
result = {
    "pairs": len(pairs),
    "pairs_with_candidate": len(ranked),
    "empirical_ids": len(ids),
    "empirical_in_candidates": len(e_c),
    "empirical_with_candidate": sum(1 for i in e_c if i in rank_of),
    "recall": {str(k): recall[k] for k in grid},
    "k_star": k_star,
    "target_stat_greater": target_wins,
    "decoy_stat_greater": decoy_wins,
}
json.dump(result, open(out, "w"), indent=1)
print(json.dumps(result, indent=1))
