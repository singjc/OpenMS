# WP-00: Acceptance baseline and datasets

- **Track:** baseline.
- **Status:** proposed, updated 21 October 2026. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10662](https://github.com/OpenMS/OpenMS/issues/10662).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Establish reproducible scientific, runtime and storage comparisons without blocking independent implementation.

## Scope

- Pin an upstream/CI `develop` revision and the applicable #10642 revision/state, commands, dependencies and input checksums.
- Pin the established TimSim and exploratory Astral inputs below. Select separate known-site data for localization validation, which neither current dataset supplies.
- Record cold preparation and warm reuse, stage/total runtime, peak memory, actual peak disk and output row cardinalities.
- Inventory OSW `FEATURE_TRANSITION` fanout and rescoring cost for WP-08, alongside identification, localization and quantitative error.
- The local 600pt fixture is a regression aid, not the baseline. Repairing local `libnethost` is not a prerequisite.

## Established datasets and scripts

| Dataset | Current search and validation scope |
|---|---|
| TimSim fixture `validation_3x3` | Six diaPASEF inputs: three controls and three treatments. Run target-only and target-plus-entrapment transition libraries against packaged realized identification/abundance and entrapment truth. |
| Astral HeLa 50 ng, 15 min | Three technical replicates, currently searched with a randomly sampled 100k-peptide FASTA. Exploratory identification, within-run entrapment diagnostics and repeatability only, not official full ProteoBench global precursor validation or biological protein performance. No abundance truth is supplied. |

Use [tools/scripts/opendia_benchmark](../../../../tools/scripts/opendia_benchmark/README.md): `build.sh`, `run_timsim.sh`, `run_astral.sh`, `validate_timsim.py`, `metrics_astral.py`, `make_report.py` and `render_report.py`. The runners capture commands, source/binary/input provenance, elapsed time, process peak RSS and sampled output-directory disk use. `monitor_host.py` and optional systemd memory/swap limits add system-resource context, not proof of whole-process/cgroup peaks or stage attribution. Machine, total RAM/logical threads, available RAM and idle CPU are now captured per search. Fresh target-only/mixed six-file TimSim and three-file Astral searches completed successfully on 21 October 2026 at source `d3b00d382efc`, with fresh logs, metrics and a regenerated report/PDF. TimSim target/mixed wall times were 6 min 55 s / 11 min 11 s with 0.84 / 0.85 GiB process peak RSS. Astral took 3 h 8 min 31 s with 10.41 GiB process peak RSS. All used four configured threads on a 30.80 GiB, 14-thread laptop, with other host activity, so these are not isolated speed comparisons. Results and reporting artifacts are tracked in [#10662](https://github.com/OpenMS/OpenMS/issues/10662). Repeat-tolerance, stage/cgroup/peak-disk and scientific acceptance gates below remain open.

The standard view requires 1% run-level precursor, global peptide and global protein q-values. Astral FDP diagnostics instead use run-level precursor `m_score` without the extra peptide/protein filters. Global precursor q-values and the exact modified-peptide pair mapping are unavailable, so global precursor validation and paired FDP remain blocked. Combined FDP also requires an independently verified effective entrapment:target search-space ratio, which the sampled FASTA does not establish. Neither current dataset closes the known-site localization gate.

## Dependencies

- **Start / local merge:** none.
- Supply pinned inputs and thresholds to WP-10. WP-01 and the WP-07 early experiment can proceed independently.

## Local merge criteria

- [ ] A manifest records exact source/CI and #10642 provenance, checksums, parameters and executable commands for both representative datasets.
- [ ] Fresh successful producer logs and reference outputs are attached for Astral and diaPASEF, with no stale-output reuse.
- [ ] A repeat baseline run reproduces identities and scores/quantities within recorded numerical tolerances.
- [ ] The measurement table includes stage runtime and memory, process/cgroup peak memory, rescoring share, table row counts, final bytes and observed peak disk including sidecars/staging.
- [ ] Numerical quality, runtime and results/rescore size gates are recorded before candidate-engine evaluation, with unavailable truth explicitly blocked.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
