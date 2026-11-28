# AMEND_ENTRAPMENT: peptide-level entrapment filter (PLAN_FASTA)

Written 28 November 2026. This is before any entrapment search, and before any FASTA-library ID count, FDP or survey
statistic has been seen. It amends only the "Entrapment FASTA" section of PLAN_FASTA. The criteria, arms, r0
definition, tolerances and outcome rules are unchanged.

## Why
The frozen protein-level rule is infeasible on `foreign4.fasta` against the human reviewed FASTA. A foreign protein
is dropped if **any** of its peptides is within Hamming ≤ 2 (I = L) of a same-length human candidate peptide.

| | value |
|---|---|
| foreign proteins | 27,299 |
| dropped by the protein rule | 26,326 |
| surviving proteins | 973 |
| their distinct peptides | 8,875 |
| ratio to the 1,317,899 human peptides | 0.0067 (target 0.5) |

This happens because short peptides almost always have a near neighbour. The admissible share of foreign peptides by
length is 1.5 % at length 7, 31.7 % at 8, 80.9 % at 9, and above 95 % from length 10. A protein almost always has at
least one short peptide, so it is almost always dropped. At a ratio of 0.0067, the expected number of entrapment IDs
at 1 % FDR is a handful, so criterion 2 could not be measured. Output of the frozen rule: `entrapment/v1_infeasible/`.

## Amended rule (peptide level)
1. Every protein of `foreign4.fasta` is digested with the tool's defaults: Trypsin/P, 1 missed cleavage,
   length 7-30, standard residues only.
2. A product is **admissible** only if both of these are at Hamming distance > 2 (I = L) from every same-length human
   candidate peptide:
   - the product itself;
   - every product of re-digesting it with the same settings. This covers the pieces a 1-missed-cleavage peptide
     contains, so OpenDIA's own digestion of the entrapment FASTA cannot create a non-admissible peptide.
3. Each admissible peptide is written as its own FASTA entry, `ENTRAP_<accession>_<n>`.
4. Proteins are drawn in a random order (seed 20261128). Each drawn protein adds all its admissible peptides that
   have not been seen yet. Drawing stops once the distinct entrapment peptides reach 0.5 × the human target peptides.
5. The combined FASTA is the human FASTA followed by the entrapment FASTA. A precursor is an entrapment if all its
   proteins start with `ENTRAP_`. This is unchanged.

## Consequences, stated before results
- **Composition differs from the human set.** Short entrapment peptides (7-9 aa) are under-represented compared with
  human peptides. Missed-cleavage forms appear only when their pieces are admissible too. This is why r0 is taken
  from the **candidate set**, as frozen: it is the entrapment/human ratio of target precursors that OpenDIA actually
  enumerates. The same r0 applies to both arms Pe and Fe, so the comparison is not biased by composition.
- **Proteins.** Entrapment proteins are single-peptide entries. Protein-level counts are not used for entrapment.
  Precursor FDP_c is the only entrapment estimate, as frozen.
- The script and the SHA-256 checksums of the source, entrapment and combined FASTAs are committed with this
  amendment, before any Pe or Fe search.
