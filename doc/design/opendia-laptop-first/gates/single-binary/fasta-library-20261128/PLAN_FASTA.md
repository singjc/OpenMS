# PLAN_FASTA: FASTA in silico library with a data-driven, target-decoy-symmetric prefilter

Written 28 November 2026, before any code for the feature exists. No FASTA-library search has been run, so no ID count,
FDP or survey statistic for it has been seen. The numbers that were known when this plan was written:
- the feasibility benchmarks (`~/.jcode/scratch/fasta-feasibility-20261128/`);
- the empirical-library reference on HeLa D (27,410 precursors at q ≤ 0.01);
- the composition of those 27,410 IDs (see Defaults).

## Goal
`OpenDIA -fasta human.fasta -in runs... -out_dir out` builds one fixed library and analyses every run with it, with no
spectral library. The steps:
1. enumerate candidates from the FASTA;
2. predict the cheap RT and CCS for all of them;
3. survey every run with theoretical fragments;
4. keep the best pairs;
5. predict the expensive MS2 spectra only for the pairs that are kept;
6. write and cache `library.parquet`.

## Pipeline (all of it under `library/preparation/` and `workflow/`, no OpenSWATH code)
1. **Candidates.** OpenMS `ProteaseDigestion` then `ModifiedPeptideGenerator`, giving targets (peptidoform, charge),
   each with every accession it maps to.
   - Charges from `-min_charge` to `-max_charge`. Precursor m/z must be in [`-min_mz`, `-max_mz`].
   - **Decoys:** one pseudo-reversed decoy per target peptidoform, the existing OpenDIA rule factored out of
     `decoySequences`, so the result is byte-identical.
   - A target whose decoy is impossible is dropped, together with all its charges. That covers terminal modifications,
     and decoys that collide with a target or an earlier decoy. The library therefore holds only complete pairs.
   - **Modifications:** only UniMod modifications with a known formula are accepted. PeptDeep encodes them by element
     composition. Terminal modifications are rejected, because they cannot get a decoy.
2. **Coordinates.** PeptDeep predicts RT for every target and decoy peptidoform. CCS is predicted for every precursor
   when ion mobility is on, then converted to 1/K0.
   - Library RT = 100 × the PeptDeep RT. This is a unit only, since the calibration maps any monotone scale.
3. **Survey library.** Every pair gets up to 6 theoretical *survey fragments* with no MS2 model. The rule:
   - y ions at charge 1, ordinal 3 to n-1, with m/z in [`-fragment_min_mz`, `-fragment_max_mz`];
   - the 6 closest to ordinal 7, lower ordinal first on ties;
   - topped up with b ions at charge 1 (ordinal ≥ 3, closest to 5) when fewer than 6 y ions qualify.

   On the empirical HeLa library this rule hits 3.76 of each precursor's 6 empirical fragments on average. The "ymid"
   row of the pre-plan check covered 5 rules, measured on fragment annotations only, with no search.
   - **Calibration panel:** up to 40,000 pairs whose target is unmodified (fixed mods allowed), charge 2, 0 missed
     cleavages and length 7-20.
     - They are taken at an even stride over predicted RT, and their decoys come along.
     - They get predicted MS2 (top 6) so the identification-free calibration sees intensities. Only these pairs are
       calibration pairs.
     - The selection is label-free and uses no data.
4. **Survey (data-driven prefilter).** For each run, `RunAnalysis::survey` runs:
   - the unchanged calibration;
   - the unchanged `StreamSearch`, with the calibrated RT, m/z and 1/K0 tolerances.

   There is no extraction and no scoring.
   - Statistic per precursor and run: the largest `coelution` of its search candidates (0 with no candidate).
   - Pair statistic: S = max over runs and over {target, decoy}.
   - The `-prefilter_pairs K` pairs with the largest S > 0 are kept. Ties are broken by FNV-1a of the target's
     TRAML_ID.
   - `K = 0` turns the prefilter off, and every pair is predicted.
5. **Final library.**
   - PeptDeep MS2 for every member of every kept pair, in batches.
   - Fragments: b and y ions of ordinal ≥ 3, fragment charge ≤ max(1, min(2, z-1)), m/z within the fragment range,
     intensity > 0, top `-fragments N` by intensity (default 12).
   - A pair is dropped if either member has fewer than 3 fragments. The rule is symmetric and needs no data.
   - The search uses the top 6 fragments, scoring uses the top min(N, 12). There are two Parquet files with the same
     IDs: `<stem>.library.parquet` (N fragments) and `<stem>.library.search.parquet`.
6. **Cache.**
   - `<stem>.candidates.parquet` is keyed by: the FASTA SHA-1, every enumeration and prediction setting, the SHA-1 of
     the 3 model files, and the OpenMS version and revision.
   - The final library has the candidates key plus: the fragment settings, K, the survey settings, and every run's
     canonical path and size.
   - A matching `.key` file means the file is reused. Writing goes to a temporary name and is then renamed.
   - The analysis then runs unchanged on the final library.

## Defaults (chosen before results)
- Trypsin/P, 1 missed cleavage, length 7-30, charge 2-4, m/z 400-1200, fixed Carbamidomethyl (C), variable
  Oxidation (M), at most 1 variable modification.
- Fragment m/z 200-1800, N = 12.
- Instrument and NCE:
  - `auto` means timsTOF for `.d` input, SciexTOF for `.wiff`, and QE otherwise;
  - NCE is 30;
  - `-library_im auto` predicts CCS only for `.d` input.
- From the 27,410 empirical IDs on D:
  - charge 2/3/4/1 = 20,547 / 6,233 / 601 / 29;
  - 27,382 have length 7-30;
  - missed cleavages 0/1/2 = 22,123 / 4,764 / 523;
  - 414 carry N-terminal acetyl. These are lost by design, because terminal modifications get no decoy.

  About 1.5 % (acetyl) plus 1.9 % (2 missed cleavages) of those IDs are outside the defaults.
- `-prefilter_pairs K`: K = K* is fixed by the design rule below.

## FDR reasoning and caveats (stated before results)
- **Why the selection is safe in principle:** a pair is selected by a *symmetric* function of its target and decoy
  statistics (max), and with the same rule for both. For a null pair, the target and decoy are exchangeable before
  selection, and they stay exchangeable after it. Selecting on the target alone would push false targets above their
  decoys, so that is not done.
- **What this rests on:**
  - The decoy's survey fragments, RT and 1/K0 come from its own sequence by the same models and rules as the target's.
  - The same runs are used for the survey and for the search (double use). Under exchangeability this is harmless,
    but it is the main risk, and the entrapment arms test it.
- **The known diaPASEF excess remains:** FDP_c of about 1.4-1.6 % at q ≤ 0.01 (PGENT2). Absolute FDP below 1 % is
  therefore *not* a criterion. The criterion is relative: prefilter against no prefilter.
- **Entrapment estimator:** FDP_c = N_e (1 + 1/r0) / N. Here r0 is the entrapment/human ratio of target precursors in
  the *candidate* set, before the prefilter.
  - The prefilter removes absent human peptides and entrapments alike, so r0 is the ratio among false candidates.
  - The ratio in the final library is reported as a sensitivity check only.

## Design stage (fixes K*, before any final-library search)
- On HeLa D, run the survey alone with the human FASTA.
- E_c is the 27,410 empirical IDs that exist in the candidate set, matched by modified sequence (UniMod) and charge.
  The empirical result files are read-only.
- recall(K) = |E_c among the kept targets| / |E_c|.
- **Rule:** K* is the smallest K in {100k, 200k, 300k, 500k, 800k} with recall ≥ 0.95. If none reaches it,
  K* = 800k.
- K* is written into `DESIGN_K.md`, committed, then made the code default. It is not changed afterwards.

## Entrapment FASTA (built before any search)
- Source: `foreign4.fasta` on codon (Arabidopsis, rice, maize and yeast, reviewed).
- Each foreign protein is digested with the tool's defaults. A protein is dropped if any of its peptides is within
  Hamming distance ≤ 2 (I = L) of a same-length human candidate peptide.
- Proteins are drawn at random (seed 20261128) until entrapment target peptides are about 0.5 × the human target
  peptides.
- Accessions are renamed `ENTRAP_<accession>`. A precursor is an entrapment if every one of its proteins starts
  with `ENTRAP_`.
- The script and the SHA-256 checksums are committed with this plan.

## Arms (codon, OpenDIA at the feature commit, seed 0 unless noted)
| arm | FASTA | prefilter | runs |
|---|---|---|---|
| P | human | K* | D (seeds 0, 1, 2), H |
| F | human | off (K = 0) | D |
| Pe | human + entrapment | K* | D |
| Fe | human + entrapment | off | D |

D = `20200505_..._S3-A1_1_2737.d` (development), H = `20201207_..._S1-C8_1_22474.d` (held out).
The library is built once per arm and FASTA. The seeds only change scoring.

## Criteria (all must hold)
1. **No regression.**
   - All OpenDIA CTests pass.
   - The empirical-library path on D, run locally, gives legacy outputs and calibration byte-identical to
     `entity-inference-20261127/out_pasef_owned`.
2. **FDR, prefilter against full.** At q ≤ 0.01 on D, FDP_c(Pe) ≤ FDP_c(Fe) + 0.003.
3. **IDs, prefilter against full.**
   - N(P, D, seed 0) ≥ 0.97 × N(F, D) at q ≤ 0.01.
   - Also N@FDP_c 1 %(Pe) ≥ 0.97 × N@FDP_c 1 %(Fe).
4. **Cost.** On D, P's wall time from FASTA to results (cold cache) ≤ 0.5 × F's.

**Goal, reported but not a pass condition:** N(P, D, seed 0) ≥ 24,669 (0.9 × 27,410, the empirical-library result).

## Outcome
- **All pass:** `-fasta` ships with prefilter K*.
- **2 fails:** the prefilter becomes opt-in (default K = 0) with a stated warning, and the user is asked.
- **3 or 4 fails:** reported, and the user decides.
- Nothing is retuned after results. H, the seed spread and the library ratio are reported either way.
- A missed goal is reported honestly, together with the composition of the missed IDs. Any follow-up gets a new plan.
