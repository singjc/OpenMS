# RESULT_PREFILTER_FDP: why prefiltered searches show higher entrapment FDP (28 November 2026)

Plans: PLAN_PREFILTER_FDP (67cea0d45b; causal arms, jobs 45366000-03) and ADDENDUM_H4 (5368f2a9ce). Each was written
before the analysis it governs. FDRBench database, D, SIF v4.

## Causal arms (FDP_c(r0), `eval_fb.py` procedure via `d1_exchange.py`)
| arm | N | T_e / D_e | FDP_c | pre-GBM network FDP_c | top-0.1 % null target share |
|---|---|---|---|---|---|
| FbF (K = 0, ref) | 8,074 | 49 / 38 | 1.21 % | | 0.509 |
| FbV1 s0 / s1 / s2 | 14,994 / 14,774 / 14,760 | 129/62, 118/64, 119/59 | 1.72 / 1.60 / 1.61 % | - / 1.68 / 1.55 % | 0.615 / 0.602 / 0.604 |
| FbV1_g1 (1 GBM round) | 14,757 | 124 / 68 | 1.68 % | 1.62 % | 0.591 |
| FbP / FbP_g1 | 11,120 / 11,307 | 85/51, 83/51 | 1.53 / 1.47 % | - / 1.50 % | 0.590 / 0.558 |

- **H1 (extra GBM rounds): rejected** (FbV1_g1 1.68 % ≥ 1.60 %). The confound is real (K = 0 runs one round), but it
  is not the cause.
- **H3 (already tilted before the GBM): partly supported.** Network-score FDP_c is 1.55-1.68 % in V1 arms.
- **Seed:** the FbV1 excess is robust (1.60-1.72 %).
- **Original guess, survey selection tilts null pairs: rejected** by the diagnostics. Kept entrapment pairs have a
  survey target share of 0.502-0.506, and TDC itself is exact (FDR_tdc = 1.00 %).

## H4: near-neighbour signal borrowing (ADDENDUM_H4 rules)
FDRBench entrapment peptides are composition-preserving shuffles of a human parent: same mass, same C-terminus.
- Accepted null target winners co-elute with their identified human parent in **15.6 %** of pairs (126 / 807,
  pooled prefilter arms). Accepted null decoy winners: **0.2 %** (1 / 450). Difference 0.154 (95 % CI 0.128-0.179),
  in every arm separately too. **H4 supported.**
- Parent-co-eluting share of accepted null pairs: **10.1 %** with a prefilter against **2.3 %** at K = 0 (met).
- Without parent-co-eluting pairs, pooled T_e/D_e falls from 1.79 to 1.52 (CI 1.35-1.71). It does not reach FbF's
  1.32, so the "explains the excess" rule is **not fully met**. Against FbF's ratio, these pairs account for about
  **60 %** of the excess target wins (125 of 211). The rest is unexplained. Neighbours whose parent was not itself
  identified, or other similar peptides, are plausible but untested.
- Fragment similarity to the parent (post hoc, `h5_similarity.py`): an entrapment target shares 22-23 % of its b/y
  fragment masses with its parent, its OpenDIA decoy 5 %. Among accepted target winners: 32-34 %, with 23-26 % of
  them sharing ≥ 50 %. Decoys never reach 50 %.

## Mechanism (what the data support)
1. Decoys model random matches. A false target that shares fragments with a co-eluting present peptide beats its
   decoy systematically. Near-neighbour nulls are not exchangeable with their decoys.
2. At K = 0, millions of random nulls set a high score threshold, so near-neighbour nulls are a small share of what
   crosses it (2.3 %).
3. The prefilter removes random nulls, targets and decoys alike (symmetric). It keeps near-neighbour nulls, because
   their target has real borrowed survey signal. The threshold drops, and near-neighbour nulls become about 10 % of
   accepted null pairs. TDC cannot see them, so FDP rises above q.

So the prefilter does not create false targets. It changes the mix of nulls toward the kind decoys do not model.

## How much this applies to real searches (open)
- **FDRBench makes near-neighbour nulls by construction.** Every entrapment peptide has an isobaric parent sharing
  about 23 % of fragments, and many parents are present in the sample.
- **The custom entrapment removed them by construction** (Hamming > 2 from any human peptide). It showed no
  prefilter increase (Pe 1.34 % against Fe 1.63 %). This reconciles the C2 sign disagreement in
  RESULT_ENTRAP_FDRBENCH.
- **Real false human targets sit between the two.** Paralogs, nested missed-cleavage forms and isobaric variants
  are near neighbours, but most absent peptides are not.
- **A crude pure-human check** (overlap.py, earlier, own-rank ≤ 800k stratum): IDs gained by V1+pool and absent
  from F are nested-related to another ID at a rate no higher than a same-size sample of shared IDs (8.8 % against
  11.6 %).
- The real-data FDP of the prefilter is therefore bounded by these two entrapments, not measured.

## Remedy candidates (not implemented; each needs its own frozen plan)
- **Neighbour competition:** a target whose matched fragments are explained by a co-eluting, better-scoring
  precursor in the same window loses them (or is removed). This is standard in DIA tools, and it is target-decoy
  symmetric in rule. It applies mostly to targets because decoys rarely have such neighbours, which is exactly the
  missing correction.
- **A shared-fragment feature:** a scoring feature for the fraction of a candidate's fragments shared with
  co-eluting accepted IDs.
- Either can be tested on the existing FDRBench arms (FDP against N), then on the custom entrapment.
