# RESULT_FASTA: OpenDIA FASTA library with a survey prefilter (28 November 2026)

Plan: [PLAN_FASTA](PLAN_FASTA.md), frozen at docs 63cbf06a49. Amendments, both written before the runs they affect:
- [AMEND_ENTRAPMENT](AMEND_ENTRAPMENT.md): per-peptide entrapment, because the frozen protein rule kept only 0.0067 × human.
- [AMEND_CALIBRATION](AMEND_CALIBRATION.md): the final search calibrates on the panel pairs only, variant A1. The user did not choose this; I adopted it myself. Code is at OpenMS 7fb7a6fb8d, SIF `8ce13018…3e5e`.

The design rule fixed K* = 800,000, the fallback, since no K reached recall 0.95 ([DESIGN_K](DESIGN_K.md)).

## Outcome
**Criterion 3 fails (IDs).** By the plan, a failure of 3 is reported and the user decides. Nothing was retuned. `-fasta` works, and the prefilter is the frozen default K* = 800,000 at this commit. **The user must decide whether to ship the prefilter on by default.** The evidence is mixed. On the human FASTA, K = 0 gives 58 % more IDs (25,121 against 15,917). On the 1.6× larger human + entrapment space, K* gives more (13,132 against 8,543 at FDP_c 1 %). See "Reading" below.

| criterion | rule | result | |
|---|---|---|---|
| 1 no regression | 29/29 OpenDIA CTests; empirical path on D byte-identical | 29/29. Legacy outputs and calibration byte-identical to `out_pasef_owned`, at a445c1892e and at 7fb7a6fb8d. v1 and v2 give identical bytes in all 17 output files. | **pass** |
| 2 FDR | FDP_c(Pe) ≤ FDP_c(Fe) + 0.003 | 0.0134 vs 0.0163 | **pass** |
| 3 IDs | N(P) ≥ 0.97 N(F), and N@FDP_c 1 %(Pe) ≥ 0.97 × the same for Fe | N: 15,917 vs 25,121, **0.634**. N@FDP_c: 13,132 vs 8,543, 1.54 | **fail** (first part) |
| 4 cost | wall(P) ≤ 0.5 wall(F), cold cache | 1,311 s vs 3,981 s, 0.33 | **pass** |
| goal | N(P) ≥ 24,669 | 15,917 | **not met** |

## Arms (codon, 16 cores, D = S3-A1, H = S1-C8, cold library cache except P1/P2)
| arm | FASTA | K | N (q ≤ 0.01) | N_e | FDP_c | N@FDP_c 1 % | wall | peak RSS | final-search calibration pairs / anchors |
|---|---|---|---|---|---|---|---|---|---|
| P | human | 800k | **15,917** | - | - | - | 21:51 | 16.4 GB | 13,186 / 540 |
| P1 (seed 1) | human | 800k (P's library) | 16,030 | - | - | - | 15:01 | 12.4 GB | |
| P2 (seed 2) | human | 800k (P's library) | 16,025 | - | - | - | 14:11 | 12.4 GB | |
| F | human | 0 | **25,121** | - | - | - | 1:06:21 | 36.0 GB | 38,372 / 608 |
| Pe | human + entrapment | 800k | 13,501 | 67 | 0.0134 | **13,132** | 25:50 | 26.0 GB | 9,157 / 314 |
| Fe | human + entrapment | 0 | 9,270 | 56 | 0.0163 | **8,543** | 2:04:13 | 57.4 GB | 38,856 / 342 |
| H (held out) | human | 800k | 37,830 | - | - | - | 29:58 | 23.8 GB | 16,065 / 1,171 |

- **r0** comes from the candidate set: 2,023,491 entrapment / 3,433,231 human targets = 0.5894, the same for Pe and Fe. That gives FDP_c = 2.697 N_e / N.
- **Library ratio** (sensitivity only): Pe final library 0.585, Fe final library 0.606. So the prefilter removes entrapments and absent human peptides at nearly the same rate.
- **Seed spread** for P: 15,917 / 16,030 / 16,025, a range of 113 (0.7 %).
- **v1** (frozen code a445c1892e, [INTERIM_V1](v1_frozen_code/INTERIM_V1.md)): P = 15,942 and F failed calibration. The amendment moved P by -25 (-0.16 %).
- **Reference points** on D: the empirical library gives 27,410 at q ≤ 0.01. DIA-NN gives 26,006 at FDP_c 1 % (BOARD).

## Reading (stated after the results)
- **The prefilter costs IDs on the human FASTA.** F finds 25,121 and P 15,917. Of the 11,116 F IDs that P misses, 10,141 (91 %) are not in P's 800k-pair library at all. So the survey statistic ranks them below K*. This matches the design-stage recall of 0.58. The prefilter, not scoring, is the loss.
- **With entrapment, the order reverses.** Fe (5.46M pairs, of which 1.84M entrapment targets) drops to 9,270 at q ≤ 0.01. Its N@FDP_c 1 % is 8,543, below Pe's 13,132. So on the 1.6× larger search space, the full library loses more power than the prefilter does. This is the classic search-space cost. Both arms keep FDP_c in the known diaPASEF excess band (1.3-1.6 %), and the prefilter does not worsen FDR (criterion 2).
- **No FDR damage from the double use of the runs.** FDP_c(Pe) 0.0134 ≤ FDP_c(Fe) 0.0163, and the library ratio stays at r0. Together these support the symmetric (max) selection argument on these data.
- **Cost.** F runs 3× longer than P, mostly MS2 prediction plus the larger search and reports. Fe peaks at 57 GB.
- **The held-out run H** reaches 37,830 with P. That is 90 % of DIA-NN's 41,272 on H at FDP_c 1 % (not the same metric). The survey on H calibrates with 1,218 anchors, against 608 on D.

## What the user decides (per PLAN_FASTA "Outcome": 3 fails → reported, user decides)
1. **Default K.** The options are K = 0 (F: more IDs on a human-only FASTA, 3× slower, about 36 GB) or K* = 800k (faster, but -37 % IDs on D). Today's code default is K* = 800k, as frozen. A follow-up plan could raise K, or improve the survey statistic, since 0.58 recall is the limit. That would be a new plan, not a retune.
2. Whether AMEND_CALIBRATION, adopted without an explicit choice, is accepted. If it is rejected, v1 stands: P = 15,942, and F is infeasible.

## Evidence
- `v2/eval_v2.json`: `eval_fasta.py` output on slim per-precursor tables, with the columns the evaluator reads. Full-table sha256 are in `v2/full.sha256`; the full outputs stay on codon in `fasta-20261128/runs_v2/`.
- `v2/ratios.json`: r0 and the library ratios, computed on codon by `codon_v2/ratio_v2.sh` with the same function.
- `v2/missed_composition.json`, `v2/logs/`: per-arm logs and `/usr/bin/time`.
- `codon_v2/`: sbatch, submit, SIF hash, container build, the criterion 1 recheck.

## Diagnostic: why the prefilter drops real peptides (after the results, no rule; `why/why.json`)
"True" here means F's targets at q ≤ 0.01 (25,121; FDP_c about 1-1.6 %). P's survey statistic S is the largest co-elution over 6 theoretical survey fragments, the sum of positive Pearson r, range 0-6.
- **40 % of F's IDs are cut** (10,140), and nearly all of them have S > 0 (only 0.9 % have no survey candidate). They are seen, but score low. Their median S is 1.43, against 1.95 at the K* cut and 2.97 for the kept IDs. Most absent candidates sit at S 1.0-2.0 too (median of all pairs 1.42).
- **Their survey signal is at noise level.** In the cut group, the target's co-elution beats its own decoy's only 58 % of the time, against 85 % for the kept group and 42 % over all pairs.
- **Precursor charge is the strongest factor.** 70 % of the cut IDs are charge 3 or 4, against 30 % of the kept IDs. The survey fragments are y ions at charge 1 only, which fit high-charge precursors badly. Missed cleavages, oxidation, length and predicted-RT position do not differ.
- The cut IDs are harder ones: 29 % are in the empirical library's IDs, against 61 % of the kept. F still scores them confidently, with median discriminant 15.8 against 16.5.
- Keeping 95 % of F's IDs would need K ≈ 2.49M of 3.43M pairs (90 %: 2.20M, 80 %: 1.70M). The current statistic cannot reach 0.97 at any useful K.
