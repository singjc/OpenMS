# DESIGN_K: the prefilter size K* (PLAN_FASTA design stage)

Written 28 November 2026, after the survey-only design run on HeLa D and before any final-library search. The rule
was frozen in PLAN_FASTA and is applied unchanged.

## Run
- OpenDIA at `ea50326e31`, plus the `-fixed_mods`/`-variable_mods` list-copy fix (`d3dd1e8f1a`; no logic change).
  Run locally with 16 threads and `-survey_only`, at the defaults.
- FASTA: human reviewed, sha256 `b205ed99…63`. Run D: `…S3-A1_1_2737.d`.
- Empirical reference: read-only `entity-inference-20261127/out_pasef_owned/legacy/precursors.tsv`, sha256
  `d01f027c…762a`, 27,410 IDs.
- Scripts: `design_survey.sh`, `design_recall.py`, `design_diag.py`. Outputs: `design_recall.json`,
  `design_diag.json`.

| stage | value |
|---|---|
| target peptidoforms | 1,725,322 |
| pairs | 3,433,234 (6,866,468 precursors) |
| dropped: no decoy / decoy collision | 1,235 / 772 |
| PeptDeep RT, 3.36M peptidoforms | 588 s |
| PeptDeep CCS, 6.87M precursors | 1,335 s |
| survey library write | 102 s, 588 MB |
| calibration panel pairs | 38,372 |
| used anchors | 608 |
| RT tolerance | 35.9 s |
| search | 23.6M candidates, 121 s |
| wall | 37:33 |
| peak RSS | 16.8 GB |

## Recall (frozen rule)
E_c = 26,185 of the 27,410 empirical IDs are candidates (modified sequence and charge). Of these, 26,122 have at
least one survey candidate.

| K | 100k | 200k | 300k | 500k | 800k |
|---|---|---|---|---|---|
| recall(K) | 0.245 | 0.322 | 0.380 | 0.473 | 0.581 |

No K reaches 0.95, so **K\* = 800,000** (the rule's fallback). It is now the default of `-prefilter_pairs`.

## What this says, stated before any final-library result
- **The survey statistic separates weakly.**
  - The median rank of an E_c pair is 570k of 2.82M pairs with a candidate. The quartiles are 106k and 1.40M.
  - For E_c pairs, the target beats its own decoy in 70 % of pairs: median target co-elution 1.91, median decoy 1.32.
  - For all other pairs, the medians are 0.99 and 0.98.
  - The max(target, decoy) rule is symmetric by construction. It ranks by the larger member, so a real pair often
    sits below null pairs whose decoy happened to co-elute well. That is the price of symmetry.
  - Theoretical fragments without intensities, and a 36 s RT window, give every precursor many chances.
- **Calibration found 608 anchors, against 7,188 with the empirical library.** The predicted-library run of 18
  November (`entlib-odp5 D_cur_s0`) found 140 anchors with a 41 s window, so a wide window from predicted RT is
  typical, not a bug of this feature.
- **Expectation for criterion 3.** With recall 0.58 at K\*, N(P) ≥ 0.97 N(F) is unlikely unless F itself loses
  most of these IDs to the larger library. Under PLAN_FASTA, a failure of criterion 3 is reported and the user
  decides. Nothing here is retuned (K, statistic, fragment rule, panel).
- **Speed note.** The RT and CCS rates in the pipeline were 5.7K/s and 5.1K/s. The isolated benchmark, rerun
  today, gave 24.8K/s and 12.6K/s for unmodified peptides. Both arms P and F pay the same cost, so criterion 4 is
  unaffected. The gap is a follow-up and is not investigated in this plan.
