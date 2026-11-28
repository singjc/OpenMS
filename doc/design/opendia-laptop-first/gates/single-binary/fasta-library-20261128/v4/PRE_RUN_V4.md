# PRE_RUN_V4: P', Pe', H' (V1+pool), written before submission (28 November 2026)

- Code 143d3ac7c6 adds the opt-in `-survey_pool`: `Preparation::pooledStatistic` applied to the survey statistic before
  `selectPairs`. The rule is the one evaluated in RESULT_SURVEY_SELECT: per pair, the max over charge siblings (same
  target modified sequence) of S, plus 0.5 x the max over the pair's proteins of that protein's top-3 mean S. Pairs
  absent from the survey library (no survey fragments; 26,291 under V1 on D, none of them in F's IDs) get 0.
  The final-library cache key gains `survey_pool=1`. Unit test added; ctest -R OpenDIA 29/29.
- Offline check that the C++ rule equals the evaluated one: `cpp_rule_recall.py` on V1-D gives recall 0.8786 at
  K = 800k (evaluator: 0.878).
- Criterion 1 at defaults: `identity_v4_compare.txt`, 3 BYTE_IDENTICAL against out_pasef_owned and 15 files the same
  as the v1 identity run.
- SIF v4: `sif.sha256` (e1915832...047b).
- Arms (`arms_v4.sbatch`): Pp = P', Pep = Pe', Hp = H', each cold, `-survey_predicted -survey_pool`, default K.
- Evaluation (no retuning; rule from PLAN_FASTA as amended):
  - criterion 2: FDP_c of Pe' with r0 (0.5894) and with r_lib (Pe' final library, `ratio_v4.sh`), both against Fe's
    FDP_c = 0.0163 (AMEND_SURVEY_POOL);
  - criterion 3: P' against F (25,121) and Pe' N@FDP_c 1 % against Fe (8,543);
  - criterion 4: P' wall <= 0.5 x F wall (3,981 s), i.e. <= 1,990 s;
  - H': N against H (37,830); FH (52,847) is the K = 0 reference on H.
