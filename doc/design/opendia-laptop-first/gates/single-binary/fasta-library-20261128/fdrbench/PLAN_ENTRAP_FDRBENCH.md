# PLAN_ENTRAP_FDRBENCH: cross-check the custom entrapment against FDRBench (28 November 2026)

Written before any search on the FDRBench database. Requested by the user: "create an FDRBench entrapment proteome
(with the same fasta), run the OpenDIA prefilter + prediction, and see whether we get the same numbers and the same
FDP problem". Agreement would raise confidence in the custom entrapment (AMEND_ENTRAPMENT); disagreement would show
which one is off.

## Database (built, hashed, committed with this plan)
- FDRBench v1.1.1, the official GitHub release (`fdrbench-1.1.1.zip` sha256 d1740d90...dceed, jar 45817391...18f5).
- The README's documented protein-level shuffled command, unchanged:
  `java -jar fdrbench-1.1.1.jar -level protein -db human.fasta -o human_I2L_entrapment.fasta -I2L -diann -uniprot -fix_nc c -check`
  Input: the same human FASTA (b205ed99...3e63). Defaults: Trypsin (no P rule), seed 2000, fold 1.
- Output `human_fdrbench.fasta` (4d6272d4...f20f):
  - 20,431 human proteins with I converted to L;
  - 20,431 paired entrapment proteins (`<accession>_p_target`), each with the same composition as its human protein,
    tryptic peptides shuffled with the C-terminus fixed (3 identical to their target).
  - Two builds gave identical bytes.
- OpenDIA reads it unchanged. An entrapment precursor is one whose proteins all end with `_p_target`. OpenDIA's
  own target-decoy machinery is untouched (decoys `DECOY_`).

## Arms (D, seed 0, SIF v4 e1915832...047b, code 143d3ac7c6, cold caches; `arms_fb.sbatch`)
| arm | options | custom-entrapment counterpart |
|---|---|---|
| FbF | `-prefilter_pairs 0` | Fe |
| FbP | default (K* = 800k, V0 survey) | Pe |
| FbV1 | `-survey_predicted` | none (candidate default) |
| FbPp | `-survey_predicted -survey_pool` | Pe' |

## Estimators (`eval_fb.py`, label mode `fdrbench`; smoke-tested in mode `entrap` on Pe, Fe and Pe', reproducing
13,501 / 67 / 0.0134, 9,270 / 56 / 0.0163 and 23,546 / 82 / 0.0094 / 0.0254 exactly)
- Primary, as frozen in PLAN_FASTA: FDP_c(r0) = N_e (1 + 1/r0) / N, r0 = entrapment/human target precursors in each
  arm's candidates. For FDRBench, r0 should be about 1 by construction.
- Also FDP_c(r_lib) (AMEND_SURVEY_POOL), 95 % Poisson CIs, and N@FDP_c 1 %.
- Secondary: FDRBench's own combined estimator at protein-pair level (`-level precursor -fold 1`, its README input
  format, decoys removed). Paired estimation needs a peptide-level pair file, which protein-level databases do not
  have, so it is not used.

## Comparisons fixed now (reported; no default decision rests on them alone)
- **C1 base agreement:** FbF and Fe FDP_c(r0) agree if their 95 % CIs overlap.
- **C2 same prefilter effect:** the change FDP_c(prefilter arm) - FDP_c(K = 0 arm) has the same sign and overlapping
  CIs between FbP and Pe (V0) and between FbPp and Pe' (V1+pool).
- **C3 pooling asymmetry:** r_lib / r0 for FbPp compared with Pe' (0.27). If FDRBench entrapment proteins also lose
  out under the protein term, r_lib / r0 will be well below 1 there too. That would confirm the mechanism and that no
  entrapment design escapes it.
- **C4 ID counts:** N of FbP, FbV1 and FbPp against FbF, as a ratio; compared with P/F (0.634) and P'/F (0.992) in the
  pure-human arms.

## Predictions, written before results
- FDRBench entrapment proteins have realistic compositions, multiple peptides and missed-cleavage forms. So, unlike
  the single-peptide custom entries, they can collect protein evidence, but only from entrapment-level (noise)
  signal.
- Prediction: r_lib / r0 for FbPp is between 0.27 and 1, closer to 1 than for the custom entrapment, but still
  below 1.
- If C1 fails, the 1.6× base miscalibration of Fe is an artefact of the custom database and must be re-examined.

## Caveats
- I to L in the human half changes the target sequences, so N is not directly comparable with the pure-human arms.
  Only ratios within the FDRBench arms are compared.
- r0 = 1 doubles the search space relative to pure human, against 1.59 for the custom entrapment. The ID counts and
  the target-decoy competition differ for that reason too.
