# OpenDIA WP-00 benchmarks

Build, run and check the two local datasets. Keep each run in a new folder. The scripts do not overwrite a previous baseline.

The `plan/opendia_rewrite` branch is for planning and exploration, not for merging. Commit only reusable scripts and documentation. Keep datasets, run results, plots and PDF reports local. This folder ignores generated files, and the repository ignores the `opendia-benchmark/` output folder used below.

## Build and run

```bash
tools/scripts/opendia_benchmark/build.sh ../openms-build

export BUILD_DIR="$PWD/../openms-build" THREADS=4
# Optional Linux job limits. These limits affect the measured run.
export BENCH_MEMORY_MAX=20G BENCH_MEMORY_SWAP_MAX=6G

tools/scripts/opendia_benchmark/run_timsim.sh \
  /home/sing/Documents/datasets/OpenMS-TimSim-Fixture/packages/opendia_timsim_fixture_validation_3x3 \
  "$PWD/opendia-benchmark/timsim"

tools/scripts/opendia_benchmark/run_astral.sh \
  /home/sing/Documents/datasets/proteobench-astral-entrapment \
  "$PWD/opendia-benchmark/astral"
```

Run the datasets one at a time. TimSim processes all six diaPASEF files in two separate searches: target-only, then target plus absent sequences (entrapment). Its optional third argument is `target` or `entrapment`. Astral processes all three RAW files. Set `FASTA=/path/to/file.fasta` if the Astral folder has more than one FASTA file. Extra Astral arguments can follow the output folder.

An existing configured build is reconfigured. A new build uses `PRESET=linux-x64-release`. The build script enables OpenTIMS and Thermo RAW support, builds `OpenDIA` and `OpenSwathExport`, checks both executables and saves a build manifest. Use `OPEN_DIA` to select another executable.

## Shared settings and filters

Both runners use `cacheWorkingInMemory`, `batchSize=0`, `innerBatchSize=-1`, `maxConcurrentSwaths=-1` and `outer_loop_threads=-1`. The scheduler selects waves automatically. Instrument and library settings differ because the inputs differ. These settings are not proof of the fastest possible run.

The analysis threshold is fixed at 1%. Both runners compute global peptide and protein scores. Run-specific and experiment-wide peptide and protein inference are off. Gene inference is off. The standard analysis uses all three filters:

- Precursor peak-group q-value: `m_score <= 0.01`.
- Global peptide q-value: `m_score_peptide_global <= 0.01`.
- Global protein q-value: `m_score_protein_global <= 0.01`.

The raw results retain peak groups through 5% precursor q-value and retain decoys. Global peptide and protein score annotation is enabled with raw export limits of `1.01`. The native OSW exporter rejects scores greater than or equal to its limit, so `1.01` retains every valid global q-value, including `1`, while appending the scores. It is not an analyst acceptance threshold. Missing globals are still not passing values. Matrix exports and analyst standard views keep all three 1% filters unchanged. This permits checks at lower precursor limits. The metrics select the best peak group for each run, modified peptide and charge. They do not count multiple peak groups as separate precursor IDs. A missing global score is not a passing score.

## Metrics and report

The runners save `timsim/validation.json` and `astral/metrics.json`. You can repeat the analysis without rerunning OpenDIA:

```bash
python3 tools/scripts/opendia_benchmark/validate_timsim.py DATASET_DIR \
  TIMSIM_DIR/target_only/results/OpenDIA.results.tsv \
  TIMSIM_DIR/entrapment/results/OpenDIA.results.tsv --output TIMSIM_DIR/validation.json
python3 tools/scripts/opendia_benchmark/metrics_astral.py \
  ASTRAL_DIR/results/OpenDIA.results.tsv --output ASTRAL_DIR/metrics.json
```

Generate the report from completed jobs. A previous baseline is optional:

```bash
python3 tools/scripts/opendia_benchmark/make_report.py TIMSIM_DIR \
  --astral ASTRAL_DIR --previous-timsim PREVIOUS_TIMSIM_DIR --output-dir REPORT_DIR
python3 tools/scripts/opendia_benchmark/render_report.py REPORT_DIR/FINAL_REPORT.md
```

ID plots show individual runs and mean with sample standard deviation. TimSim has three control and three treatment runs. Astral has three repeat runs. Quantification checks show repeat variation and completeness. TimSim also has known abundance changes for an accuracy check. Astral has no abundance truth. Its FASTA entries are peptide sequences, not biological protein accessions. Thus Astral protein-inference groups do not establish biological protein performance.

Entrapment plots use `m_score` at each precursor q-value limit. They show absent-sequence fractions, not a search-space-adjusted FDR estimate. They do not prove calibration. No accepted IDs means an unavailable fraction, not zero error. Do not substitute a global peptide q-value for `m_score`.

### ProteoBench Astral analysis alignment

The Astral report now replaces its raw-share distribution with per-run FDP bounds
and target ID counts at each cutoff. The formulas were checked against ProteoBench
revision `d2f03663f0f606ef300211907ef2f350517edbc1`:

- [Astral entrapment methodology](https://github.com/proteobench/proteobench/blob/d2f03663f0f606ef300211907ef2f350517edbc1/docs/modules/dia/entrapment-dia-astral.md)
- [Estimator implementation](https://github.com/proteobench/proteobench/blob/d2f03663f0f606ef300211907ef2f350517edbc1/proteobench/score/entrapmentscores.py)

For accepted distinct precursors, `E` is the entrapment count and `T` the target
count. Lower FDP is `E / (T + E)`. Combined FDP is
`(1 + 1/r) * E / (T + E)`, where `r` is the effective entrapment:target
search-space ratio. ProteoBench uses `r=1` for its full matched database. We do
not assume that ratio for a randomly sampled FASTA or infer it from accepted
identifications. Combined FDP and its ratio to the cutoff are unavailable unless
the analyst supplies an independently verified ratio:

```bash
# Only for an independently verified 1:1 effective search space:
python3 tools/scripts/opendia_benchmark/metrics_astral.py \
  ASTRAL_DIR/results/OpenDIA.results.tsv --entrapment-fold 1 \
  --output ASTRAL_DIR/metrics.json
```

The default runner deliberately leaves the ratio unknown. A sampled FASTA is not
automatically the official benchmark. A count ratio alone does not verify that
sampling, modifications, charge and m/z filters preserve the matched search space.

**Scope matters:** ProteoBench tests **global precursor-ion FDR** and deduplicates
precursors across runs, keeping the best global precursor q-value. Our exported
`m_score` is run-level, so this adaptation computes within-run diagnostics only.
It does not pool run-level scores into a purported global calibration test, apply
extra peptide/protein filters to FDP counts, or issue ProteoBench valid/invalid
verdicts. Global peptide/protein scores are not substitutes for global precursor
scores. Paired FDP is explicitly unavailable: it requires the exact modified-peptide
pair mapping and appropriate score ordering, not guessed pairs from FASTA names.
No IDs gives unavailable FDP. Bounds are not clipped to 100%.

The plot compares lower/combined FDP with the nominal 1% cutoff for each run,
while the table retains all cutoffs and target IDs to show identification yield alongside
error. CV, completeness and correlations remain secondary repeatability checks,
not evidence of FDR calibration or biological protein quantification. Historical
metrics without the new `proteobench_entrapment` block still render as raw-share
diagnostics. Regenerate metrics and then the report to use the new analysis.

## Time, memory and host load

Each search saves its command, source revision/status, binary hash and input hashes under `meta/`. `logs/resources.txt` records GNU time elapsed time and peak resident process memory. This does not add the memory of simultaneous child processes. Time covers OpenDIA, including prediction when used. Input checksum checks, archive extraction and metric generation are outside that time. A complete multi-file search is one timing observation, not three or six independent measurements. The two TimSim searches are not timing replicates.

The completed Astral search took 35m18s. Its original raw TSV had blank global scores because the earlier raw export configuration disabled score annotation, although inference and native standard matrices completed successfully. Native `OpenSwathExport` re-export from the retained `workflow.osw` supplies those scores without repeating the search. That recovery is outside the original search timing. Preserve the original TSV, capture the re-export command and verify the exported row identities before using the score-annotated replacement.

The Linux monitor saves `logs/host_workload.jsonl` and a summary in `logs/host_workload.json`. Samples contain CPU use, estimated CPU use by other processes, available RAM, swap, load, pressure and free disk space. Host CPU includes the benchmark. Other-process CPU is an estimate. These observations show contention but do not correct elapsed time. The monitor stops only the measured job if free disk falls below `MIN_FREE_DISK_GIB` (default 3 GiB).

The report includes the machine and CPU models, total RAM, host logical CPU threads,
CPU affinity count and configured benchmark threads for each search. Available RAM
is Linux `MemAvailable`, recorded at monitor start and as mean/min/max during the
search. Free CPU is idle time excluding I/O wait, shown as a percentage and as
idle logical-thread equivalents (idle fraction times host thread count). These
are observed host resources, not reserved capacity. Missing historical hardware
or idle samples remain `n/a`, never substituted with the report-rendering host.
Hardware metadata does not record serial numbers, user names or host names.

Optional `BENCH_MEMORY_MAX` and `BENCH_MEMORY_SWAP_MAX` use a user systemd scope. Limits are recorded in `meta/resource_limits.txt`. A limited or busy run must not be presented as an isolated maximum-speed measurement. Previous runs without the same inference, search space, limits and workload are context only, not a controlled speed comparison.

## Checks

The runners use Bash, Python's standard library, CMake, GNU time and the dependencies of OpenDIA. PDF rendering also requires Python `markdown` and Chrome or Chromium.

```bash
python3 tools/scripts/opendia_benchmark/test_benchmark.py
python3 tools/scripts/opendia_benchmark/test_metrics.py
python3 tools/scripts/opendia_benchmark/test_make_report.py
python3 tools/scripts/opendia_benchmark/test_monitor_host.py
```

Pass captured `meta/command.sh` files to `test_benchmark.py` to check the shared speed, inference and filter settings. The test always sources the current `lib.sh` arrays to check score annotation with `1.01` raw limits. Captured older commands with both raw score flags disabled are explicitly reported as historical raw-export mode, not claimed to use the current settings. Set `OPEN_DIA` to also check the current sourced settings and captured command modes with the real executable. This parameter check deliberately stops before data processing. It does not replace a full dataset run. Do not use `-write_ini` as a parameter-resolution test: it writes defaults.

Generate `FINAL_REPORT.md`, its plots and `FINAL_REPORT.pdf` in your local `REPORT_DIR`. These outputs are not tracked in Git. Keep the original result folders and JSON files when making a new comparison.
