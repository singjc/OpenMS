#!/usr/bin/env python3
"""ADDENDUM_H4 test: do accepted null (entrapment) target winners co-elute with their identified human parent more
often than accepted null decoy winners do? Run on codon.
Usage: h4_parent.py OUT.json DIR:ARM ...   (DIR/ARM/precursors.tsv; FASTA = fasta/human_fdrbench.fasta)
"""
import csv, json, math, sys

W = "/nfs/research/juan/singj/opendia/fasta-20261128"
out = sys.argv[1]
prot = {}
acc = None
for line in open(f"{W}/fasta/human_fdrbench.fasta"):
    if line.startswith(">"):
        ident = line[1:].split()[0]
        parts = ident.split("|")
        acc = parts[1] if len(parts) >= 3 and parts[0] in ("sp", "tr") else ident
        prot[acc] = []
    else:
        prot[acc].append(line.strip())
prot = {k: "".join(v) for k, v in prot.items()}


def parents(pep, proteins):
    res = set()
    for p in proteins:
        if not p.endswith("_p_target"):
            continue
        es, hs = prot.get(p), prot.get(p[: -len("_p_target")])
        if not es or not hs or len(es) != len(hs):
            continue
        i = es.find(pep)
        while i >= 0:
            res.add(hs[i:i + len(pep)])
            i = es.find(pep, i + 1)
    return res


def ci(k, n):
    if not n: return [None, None]
    p = k / n; h = 1.96 * math.sqrt(p * (1 - p) / n)
    return [round(max(p - h, 0), 4), round(min(p + h, 1), 4)]


R = {}
for spec in sys.argv[2:]:
    d, arm = spec.split(":")
    rows = {}
    human = {}
    with open(f"{W}/{d}/{arm}/precursors.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            pid = int(r["precursor_id"])
            dec = r["decoy"] == "true"
            ps = [p for p in r["proteins"].replace(",", ";").split(";") if p]
            if dec:
                ps = [p[6:] if p.startswith("DECOY_") else p for p in ps]
            ent = bool(ps) and all(p.endswith("_p_target") for p in ps)
            q = float(r["run_q_value"])
            rec = (float(r["discriminant_score"]), dec, q, ent, r["unmodified_sequence"], int(r["charge"]), ps,
                   float(r["rt_apex"]), float(r["rt_start"]), float(r["rt_end"]))
            rows[pid] = rec
            if not dec and not ent and q <= 0.01:
                human.setdefault((r["unmodified_sequence"], int(r["charge"])), []).append((rec[8], rec[9]))
    pairs = {}
    for pid, rec in rows.items():
        pairs.setdefault(pid // 2, {})[1 if rec[1] else 0] = rec
    c = {"T": [0, 0, 0], "D": [0, 0, 0]}  # [n, parent accepted (any RT), parent co-eluting]
    for k, m in pairs.items():
        t = m.get(0)
        if t is None or not t[3]:
            continue
        w = max(m.values(), key=lambda x: (x[0], x[1]))
        if w[2] > 0.01:
            continue
        side = "D" if w[1] else "T"
        c[side][0] += 1
        par = parents(t[4], t[6])
        spans = [s for pp in par for s in human.get((pp, t[5]), [])]
        if spans:
            c[side][1] += 1
            if any(a <= w[7] <= b for a, b in spans):
                c[side][2] += 1
    nT, aT, cT = c["T"]; nD, aD, cD = c["D"]
    pT, pD = (cT / nT if nT else 0), (cD / nD if nD else 0)
    h = 1.96 * math.sqrt((pT * (1 - pT) / nT if nT else 0) + (pD * (1 - pD) / nD if nD else 0))
    R[arm] = {"T_e": nT, "D_e": nD, "T_parent_accepted": aT, "D_parent_accepted": aD, "T_coelute": cT, "D_coelute": cD,
              "p_T": pT, "p_T_ci": ci(cT, nT), "p_D": pD, "p_D_ci": ci(cD, nD), "diff": pT - pD, "diff_ci": [pT - pD - h, pT - pD + h],
              "Te_De_all": nT / nD if nD else None, "Te_De_without_coelute": (nT - cT) / (nD - cD) if nD - cD else None}
    print(arm, json.dumps(R[arm]))
json.dump(R, open(out, "w"), indent=1)
