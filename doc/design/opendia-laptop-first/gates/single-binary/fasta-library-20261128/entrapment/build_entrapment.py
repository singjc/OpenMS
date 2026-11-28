#!/usr/bin/env python3
"""PLAN_FASTA entrapment FASTA (built and committed before any search).

Source: foreign4.fasta (Arabidopsis, rice, maize, yeast; reviewed). Every protein is digested with the OpenDIA
defaults (Trypsin/P, 1 missed cleavage, length 7-30, standard residues only). A foreign protein is dropped if any of
its peptides is within Hamming distance <= 2 (I = L) of a same-length human candidate peptide (pigeonhole index: 3
segments, a match within distance 2 shares one segment exactly). Surviving proteins are drawn in a random order (seed
20261128) until their distinct target peptides (not in the human set, not already drawn) reach 0.5 x the human
target peptides. Accessions are renamed ENTRAP_<accession>. Output: the entrapment FASTA and the concatenated
human + entrapment FASTA, with a JSON summary.

--peptide (AMEND_ENTRAPMENT, before any search): the protein rule leaves 973 of 27,299 proteins (ratio 0.0067), so the
filter is applied per peptide instead. A digestion product of a foreign protein is admissible if it and every product
of its own digestion (its pieces at 0 missed cleavages) are at Hamming distance > 2 (I = L) from every same-length
human peptide. Each admissible peptide is one FASTA entry ENTRAP_<accession>_<n>, so the tool's digestion of the
entrapment FASTA gives admissible peptides only. Proteins are drawn as above and contribute all their admissible
peptides until the distinct entrapment peptides reach 0.5 x the human target peptides.

Usage: build_entrapment.py HUMAN.fasta FOREIGN.fasta OUT_DIR [--peptide]
"""
import json
import os
import random
import re
import sys

SEED = 20261128
RATIO = 0.5
MIN_LEN, MAX_LEN, MISSED = 7, 30, 1
STANDARD = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def read_fasta(path):
    entries, header, seq = [], None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq)))
                header, seq = line[1:], []
            elif line:
                seq.append(line.strip())
    if header is not None:
        entries.append((header, "".join(seq)))
    return entries


def accession(header):
    ident = header.split()[0]
    parts = ident.split("|")
    if len(parts) >= 3 and parts[0] in ("sp", "tr"):
        return parts[1]
    return ident


def digest(seq):
    """Trypsin/P (cleave after K or R, also before P), up to MISSED missed cleavages, as OpenMS digests."""
    sites = [i + 1 for i, c in enumerate(seq[:-1]) if c in "KR"]
    bounds = [0] + sites + [len(seq)]
    out = set()
    for i in range(len(bounds) - 1):
        for j in range(i + 1, min(i + 2 + MISSED, len(bounds))):
            p = seq[bounds[i]:bounds[j]]
            if MIN_LEN <= len(p) <= MAX_LEN and STANDARD.match(p):
                out.add(p)
    return out


def segments(n):
    a, b = n // 3, 2 * n // 3
    return [(0, a), (a, b), (b, n)]


def main():
    human_path, foreign_path, out_dir = sys.argv[1:4]
    os.makedirs(out_dir, exist_ok=True)
    human = read_fasta(human_path)
    human_peptides = set()
    for _, s in human:
        human_peptides |= digest(s)
    human_il = {p.replace("I", "L") for p in human_peptides}
    # Pigeonhole index over I=L human peptides: (length, segment, text) -> peptides.
    index = {}
    for p in human_il:
        for k, (a, b) in enumerate(segments(len(p))):
            index.setdefault((len(p), k, p[a:b]), []).append(p)

    def near_human(p):
        q = p.replace("I", "L")
        if q in human_il:
            return True
        for k, (a, b) in enumerate(segments(len(q))):
            for h in index.get((len(q), k, q[a:b]), ()):
                d = 0
                for x, y in zip(q, h):
                    if x != y:
                        d += 1
                        if d > 2:
                            break
                if d <= 2:
                    return True
        return False

    foreign = read_fasta(foreign_path)
    if len(sys.argv) > 4 and sys.argv[4] == "--peptide":
        return peptide_level(human, human_peptides, foreign, near_human, out_dir, human_path)
    survivors, dropped = [], 0
    for header, s in foreign:
        peptides = digest(s)
        if any(near_human(p) for p in peptides):
            dropped += 1
            continue
        survivors.append((header, s, peptides))
    rng = random.Random(SEED)
    order = list(range(len(survivors)))
    rng.shuffle(order)
    goal = RATIO * len(human_peptides)
    chosen, entrap_peptides = [], set()
    for i in order:
        if len(entrap_peptides) >= goal:
            break
        header, s, peptides = survivors[i]
        chosen.append((header, s))
        entrap_peptides |= peptides
    entrap_fasta = os.path.join(out_dir, "entrapment.fasta")
    with open(entrap_fasta, "w") as f:
        for header, s in chosen:
            acc = accession(header)
            f.write(">ENTRAP_" + acc + " " + header + "\n")
            for k in range(0, len(s), 60):
                f.write(s[k:k + 60] + "\n")
    combined = os.path.join(out_dir, "human_entrapment.fasta")
    with open(combined, "w") as f, open(human_path) as h, open(entrap_fasta) as e:
        f.write(h.read())
        f.write(e.read())
    summary = {
        "seed": SEED,
        "human_proteins": len(human),
        "human_peptides": len(human_peptides),
        "foreign_proteins": len(foreign),
        "foreign_dropped_hamming_le2": dropped,
        "foreign_survivors": len(survivors),
        "entrapment_proteins": len(chosen),
        "entrapment_peptides": len(entrap_peptides),
        "peptide_ratio": len(entrap_peptides) / len(human_peptides),
    }
    json.dump(summary, open(os.path.join(out_dir, "entrapment_summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


def peptide_level(human, human_peptides, foreign, near_human, out_dir, human_path):
    clean = {}

    def is_clean(p):
        if p not in clean:
            clean[p] = not near_human(p)
        return clean[p]

    rng = random.Random(SEED)
    order = list(range(len(foreign)))
    rng.shuffle(order)
    goal = RATIO * len(human_peptides)
    entries, seen, proteins_used = [], set(), 0
    rejected = 0
    for i in order:
        if len(seen) >= goal:
            break
        header, s = foreign[i]
        acc = accession(header)
        added = 0
        for p in sorted(digest(s)):
            if p in seen:
                continue
            if not (is_clean(p) and all(is_clean(q) for q in digest(p))):
                rejected += 1
                continue
            seen.add(p)
            added += 1
            entries.append(("ENTRAP_%s_%d" % (acc, added), p))
        proteins_used += added > 0
    entrap_fasta = os.path.join(out_dir, "entrapment.fasta")
    with open(entrap_fasta, "w") as f:
        for name, p in entries:
            f.write(">" + name + "\n" + p + "\n")
    combined = os.path.join(out_dir, "human_entrapment.fasta")
    with open(combined, "w") as f, open(human_path) as h, open(entrap_fasta) as e:
        f.write(h.read())
        f.write(e.read())
    lengths = {}
    for _, p in entries:
        lengths[len(p)] = lengths.get(len(p), 0) + 1
    summary = {
        "rule": "peptide (AMEND_ENTRAPMENT)",
        "seed": SEED,
        "human_proteins": len(human),
        "human_peptides": len(human_peptides),
        "foreign_proteins": len(foreign),
        "foreign_proteins_used": proteins_used,
        "rejected_peptides_seen": rejected,
        "entrapment_peptides": len(entries),
        "peptide_ratio": len(entries) / len(human_peptides),
        "entrapment_length_histogram": dict(sorted(lengths.items())),
    }
    json.dump(summary, open(os.path.join(out_dir, "entrapment_summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
