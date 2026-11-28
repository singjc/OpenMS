# ADDENDUM_H4: near-neighbour signal borrowing (written after the PLAN_PREFILTER_FDP causal arms, before this analysis)

## Where PLAN_PREFILTER_FDP left things (its rules, applied as frozen; full table in RESULT_PREFILTER_FDP)
- **H1 (extra GBM rounds) rejected.** FbV1_g1 (one GBM round) FDP_c is 1.68 %, at or above 1.60 %. FbP_g1 is 1.47 %.
- **H3 partly supported.** The pre-GBM network score is already elevated: FbV1_s1 1.68 %, FbV1_s2 1.55 %, FbV1_g1
  1.62 %, FbP_g1 1.50 %. One of the two named seeds meets ≥ 1.60 %, the mean of the two is 1.615 %.
- **Seed noise:** FbV1 is 1.72 / 1.60 / 1.61 % (seeds 0/1/2), all above 1.51 %, so the elevation is robust to the
  seed.

## H4 (new, analysis of existing outputs only)
An entrapment target wins its pair because it borrows real signal from a co-eluting true human peptide: same
isolation window, overlapping fragment ions. Its decoy borrows less. TDC decoys model random matches, not
near-neighbour matches, so these extra false targets are invisible to q. The prefilter raises FDP because it keeps
near-neighbour null pairs (strong survey signal, from the real neighbour) and removes most random null pairs. That
raises the near-neighbour share of the null pairs that compete near the cut. D2 cannot see this: near-neighbour
pairs are a tiny fraction of all kept null pairs.

FDRBench lets this be tested directly. Each entrapment protein shuffles its human parent's tryptic peptides with the
C-terminus fixed (no internal K/R, so peptide boundaries are preserved). So every entrapment peptide has a known
**parent** human peptide: same protein offset and length, identical composition and so identical precursor mass,
same y1. Its OpenDIA decoy (pseudo-reverse, C-terminus fixed) also has the same mass and y1. So mass and window
alone are symmetric between target and decoy. Only the fragment-level resemblance differs.

## Test (`h4_parent.py`, run on codon on runs_fb and runs_pf precursors.tsv)
- For every accepted null pair (q ≤ 0.01, entrapment pair winner): find the parent of the pair's entrapment target
  peptide (unmodified, same charge). Call the pair **parent-co-eluting** if a parent precursor is an accepted human
  target (q ≤ 0.01) and the winner's rt_apex lies within the parent's [rt_start, rt_end].
- Report the parent-co-eluting rate among target winners (p_T) and among decoy winners (p_D), with 95 % CIs, per arm.
  Also T_e/D_e after removing parent-co-eluting pairs.

## Predictions and rules, fixed now
- **H4 supported** if, pooled over the FDRBench arms with a prefilter (FbP, FbV1, FbPp, FbV1_s1, FbV1_s2, FbV1_g1,
  FbP_g1; pairs counted per arm), p_T − p_D ≥ 0.10 and the 95 % CI of the difference excludes 0. **Rejected** if the
  CI of p_T − p_D includes 0.
- **H4 explains the prefilter excess** if additionally:
  - without parent-co-eluting pairs, the pooled prefilter-arm T_e/D_e has a 95 % CI including 1.3 (FbF's value),
    or lower;
  - the parent-co-eluting share of accepted null pairs is higher in the pooled prefilter arms than in FbF.
- Fixed-depth check (descriptive, `d4_depth.py`): the target-win share among the top 100, 200 and 400 null pairs by
  pair score, per arm. "Top 0.1 %" in D3 compared different absolute depths (2,330 null pairs in FbF against 369 in
  FbV1).
