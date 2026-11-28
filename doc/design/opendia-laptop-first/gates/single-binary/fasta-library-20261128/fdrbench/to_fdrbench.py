#!/usr/bin/env python3
"""Convert an OpenDIA slim precursor table to FDRBench's precursor-level FDP input (README v1.1.1 format), targets only
(decoys removed as FDRBench requires). The entrapment label in `protein` follows FDRBench's convention: a precursor is
an entrapment if every protein carries the label; such precursors get `<first protein>` (already ending _p_target);
others get the first non-entrapment protein.

Usage: to_fdrbench.py <label_mode entrap|fdrbench> <slim.tsv.gz> <out.tsv>
For label_mode entrap (custom), ENTRAP_ proteins are relabelled to end with _p_target so FDRBench can recognise them.
score = discriminant_score (higher is better; run with -score 'score:1').
"""
import csv, gzip, re, sys

mode, slim, out = sys.argv[1:4]
QMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0  # rows above QMAX cannot change FDP at q <= QMAX
SPLIT = re.compile(r"[;,]")
is_e = (lambda p: p.startswith("ENTRAP_")) if mode == "entrap" else (lambda p: p.endswith("_p_target"))
with gzip.open(slim, "rt") as f, open(out, "w", newline="") as g:
    r = csv.DictReader(f, delimiter="\t")
    w = csv.writer(g, delimiter="\t")
    w.writerow(["run", "peptide", "mod_peptide", "charge", "q_value", "PEP", "protein", "score"])
    for row in r:
        if row["decoy"] != "false" or float(row["run_q_value"]) > QMAX:
            continue
        ps = [p for p in SPLIT.split(row["proteins"]) if p]
        ent = bool(ps) and all(is_e(p) for p in ps)
        prot = (ps[0] if mode == "fdrbench" else ps[0] + "_p_target") if ent else next(p for p in ps if not is_e(p))
        pid = row["precursor_id"]
        pep = row.get("unmodified_sequence") or f"P{pid}"
        mod = row.get("sequence") or f"P{pid}"
        w.writerow([row.get("run") or "D", pep, mod, row.get("charge") or "0", row["run_q_value"], "0", prot,
                    row["discriminant_score"]])
