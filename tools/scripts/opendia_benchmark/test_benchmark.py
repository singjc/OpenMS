#!/usr/bin/env python3
"""Small regression check. Optional arguments are captured meta/command.sh files."""

import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from make_report import parse_seconds, read_resources, write_chart
from validate_timsim import best_results, fdr_curve


def main():
    script = Path(__file__).with_name('lib.sh').resolve()
    with tempfile.TemporaryDirectory(dir=os.environ.get('JCODE_SCRATCH_DIR')) as temporary:
        root = Path(temporary)
        out = root / 'job'
        environment = dict(os.environ, OPEN_DIA='/usr/bin/true', MIN_FREE_DISK_GIB='0',
                           BENCH_CGROUP='0', ENGINE_SOURCE_DIR=str(script.parents[3]),
                           DOTNET_ROOT='/test/dotnet', BENCH_TEST_SECRET='must-not-be-captured')
        environment.pop('BENCH_MEMORY_MAX', None)
        environment.pop('BENCH_MEMORY_SWAP_MAX', None)
        environment.pop('THREADS', None)
        environment.pop('BENCH_FROZEN_HASHES', None)
        ran = subprocess.run(['bash', '-c',
            'source "$1"; prepare_output "$2"; printf "%s %s" "$THREADS" "$PERCOLATOR_THREADS" > "$2/meta/threads.txt"; '
            'run_measured "$2" env LD_LIBRARY_PATH=/test/runtime sh -c \'printf "%s" "$TMPDIR"; sleep 0.2\'',
            'bash', str(script), str(out)], env=environment, capture_output=True, text=True, timeout=10)
        assert ran.returncode == 0, ran.stderr
        assert (out / 'logs/stdout.log').read_text() == str(out / 'tmp')
        runtime = json.loads((out / 'meta/runtime_environment.json').read_text())
        assert set(runtime) == {'DOTNET_ROOT', 'DOTNET_ROOT_X64', 'OPENMS_DATA_PATH',
                                'PATH', 'LD_LIBRARY_PATH', 'TMPDIR'}
        assert runtime['DOTNET_ROOT'] == '/test/dotnet'
        assert runtime['LD_LIBRARY_PATH'] == '/test/runtime'
        assert runtime['TMPDIR'] == str(out / 'tmp')
        threads = max(1, len(os.sched_getaffinity(0)) - 3)
        assert (out / 'meta/threads.txt').read_text() == f'{threads} {min(3, threads)}'
        for name in ('source_revision.txt', 'source_status.txt', 'harness_revision.txt',
                     'harness_status.txt', 'source_locations.txt'):
            assert (out / 'meta' / name).is_file()
        assert (out / 'meta/exit_code.txt').read_text().strip() == '0'
        copied = root / 'relocated/tools/scripts/opendia_benchmark'
        copied.mkdir(parents=True)
        for name in ('lib.sh', 'monitor_host.py'):
            shutil.copyfile(script.with_name(name), copied / name)
        relocated_out = root / 'relocated_job'
        relocated_env = dict(environment, HARNESS_SOURCE_DIR=str(script.parents[3]))
        relocated = subprocess.run(['bash', '-c',
            'source "$1"; prepare_output "$2"; run_measured "$2" true',
            'bash', str(copied / 'lib.sh'), str(relocated_out)], env=relocated_env,
            capture_output=True, text=True, timeout=10)
        assert relocated.returncode == 0, relocated.stderr
        assert (relocated_out / 'meta/harness_revision.txt').read_text() == (out / 'meta/harness_revision.txt').read_text()
        locations = (relocated_out / 'meta/source_locations.txt').read_text()
        assert f'harness_source_dir={script.parents[3]}\n' in locations
        assert f'engine_source_dir={script.parents[3]}\n' in locations
        fixture = root / 'frozen_input'
        fixture.write_text('original')
        frozen_hashes = root / 'frozen.sha256'
        frozen_hashes.write_text(f'{hashlib.sha256(fixture.read_bytes()).hexdigest()}  {fixture}\n')
        frozen_env = dict(environment, BENCH_FROZEN_HASHES=str(frozen_hashes))
        fixture.write_text('changed')
        for mode in ('before', 'during'):
            checked_out = root / f'frozen_{mode}'
            if mode == 'during':
                fixture.write_text('original')
            checked = subprocess.run(['bash', '-c',
                'source "$1"; prepare_output "$2"; run_measured "$2" sh -c '
                "'printf changed > \"$1\"; touch \"$2/meta/executed\"' sh \"$3\" \"$2\"",
                'bash', str(script), str(checked_out), str(fixture)], env=frozen_env,
                capture_output=True, text=True, timeout=10)
            assert checked.returncode == 2, checked.stderr
            assert (checked_out / 'meta/executed').exists() == (mode == 'during')
            failed_log = checked_out / 'logs' / f'frozen-{"before" if mode == "before" else "after"}.log'
            assert 'FAILED' in failed_log.read_text()
        # Explicit opt-in keeps the public regression portable to non-systemd hosts.
        if os.environ.get('BENCH_TEST_CGROUP') == '1':
            for mode in ('cgroup', 'swap_only', 'nonzero', 'low_disk', 'orphan'):
                scoped = root / mode
                scoped_env = dict(environment)
                scoped_env.update({'BENCH_MEMORY_SWAP_MAX': '0'} if mode == 'swap_only'
                                  else {'BENCH_CGROUP': '1'})
                if mode == 'low_disk':
                    scoped_env['MIN_FREE_DISK_GIB'] = '1000000000'
                program = ("import time; time.sleep(10)" if mode == 'low_disk' else
                           "import time; time.sleep(0.15); a=bytearray(8*1024*1024)")
                if mode == 'nonzero':
                    program += '; raise SystemExit(3)'
                if mode == 'orphan':
                    program = ("import subprocess; from pathlib import Path; "
                               "child=subprocess.Popen(['sleep','10']); "
                               f"Path({str(scoped / 'meta/orphan_pid.txt')!r}).write_text(str(child.pid))")
                checked = subprocess.run(['bash', '-c',
                    'source "$1"; prepare_output "$2"; run_measured "$2" python3 -c '
                    + shlex.quote(program),
                    'bash', str(script), str(scoped)], env=scoped_env,
                    capture_output=True, text=True, timeout=10)
                assert checked.returncode == (2 if mode in ('low_disk', 'orphan') else 3 if mode == 'nonzero' else 0), checked.stderr
                summary = json.loads((scoped / 'logs/host_workload.json').read_text())
                if mode == 'low_disk':
                    assert summary['abort_reason']
                    continue
                if mode == 'orphan':
                    assert 'survived' in summary['abort_reason']
                    orphan_pid = int((scoped / 'meta/orphan_pid.txt').read_text())
                    try:
                        state = Path(f'/proc/{orphan_pid}/stat').read_text().rsplit(')', 1)[1].split()[0]
                        assert state == 'Z'
                    except FileNotFoundError:
                        pass
                    continue
                assert summary['cgroup_final']['cgroup_memory_peak_bytes'] >= 8 * 1024**2
                assert summary['abort_reason'] is None
                assert (scoped / 'meta/cgroup_release.txt').read_text() == 'captured\n'
                if mode == 'nonzero':
                    assert (scoped / 'meta/exit_code.txt').read_text().strip() == '3'
                else:
                    assert read_resources(scoped / 'logs/resources.txt')['seconds'] < 0.8
                assert (scoped / 'logs/scope_resources.txt').is_file()
    assert parse_seconds("1:27:54") == 5274
    assert parse_seconds("1:14.09") == 74.09
    points = fdr_curve([(0.0005, False), (0.001, True), (0.02, True)], 0.01)
    assert points[0]["entrapment_share"] is None
    assert points[1]["entrapment_share"] == 0
    assert points[2]["entrapment_share"] == 0.5
    assert points[-1]["rows"] == 2
    assert fdr_curve([], 0.01)[-1]["entrapment_share"] is None

    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        resources = root / "resources.txt"
        resources.write_text("Elapsed (wall clock) time (h:mm:ss or m:ss): 1:14.09\nMaximum resident set size (kbytes): 2048\nExit status: 0\n")
        assert read_resources(resources) == {"seconds": 74.09, "memory_mib": 2, "status": 0}
        resources.write_text(resources.read_text().replace("status: 0", "status: 1"))
        try:
            read_resources(resources)
        except SystemExit:
            pass
        else:
            raise AssertionError("A failed run was accepted")

        result = root / "results.tsv"
        columns = ["run_name", "Sequence", "Charge", "Intensity", "RT", "EXP_IM", "m_score", "FullPeptideName", "ProteinName"]
        with result.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(columns)
            for index, q in enumerate(("nan", "", "inf", "-0.1", "0.02", "0.001")):
                writer.writerow(["run", f"PEPTIDE{index}", 2, 10, 1, 1, q, f"PEPTIDE{index}", "TARGET"])
        rows, q_column = best_results(result, 0.01, ("m_score",))
        assert len(rows) == 1 and q_column == "m_score"
        metrics_script = Path(__file__).with_name("metrics_astral.py")
        completed = subprocess.run([sys.executable, metrics_script, result], capture_output=True, text=True, check=True)
        assert json.loads(completed.stdout)["identification"]["target_rows"] == 1
        rejected = subprocess.run([sys.executable, metrics_script, result, "--q-threshold", "nan"], capture_output=True, text=True)
        assert rejected.returncode == 2 and "between 0 and 1" in rejected.stderr

        # Chrome can exit zero without writing a PDF. Keep the old report intact.
        fake_browser = root / "google-chrome"
        fake_browser.write_text("#!/bin/sh\nexit 0\n")
        fake_browser.chmod(0o755)
        report = root / "report.md"
        report.write_text("# Report\n")
        pdf = report.with_suffix(".pdf")
        pdf.write_bytes(b"%PDF-old report")
        environment = dict(os.environ, PATH=f"{root}:{os.environ.get('PATH', '')}")
        rendered = subprocess.run([
            sys.executable, Path(__file__).with_name("render_report.py"), report,
        ], env=environment, capture_output=True, text=True)
        assert rendered.returncode != 0
        assert "did not create a PDF" in rendered.stderr
        assert pdf.read_bytes() == b"%PDF-old report"

        for values in ([0, 1, None], [-1, 0, 1], [None, None, None]):
            chart = root / "chart.svg"
            write_chart(chart, "A & B", "Value", ["a", "b", "c"], [("Data", "blue", values)])
            tree = ET.parse(chart)
            for rectangle in tree.findall("{http://www.w3.org/2000/svg}rect"):
                if not rectangle.get("height", "").endswith("%"):
                    assert float(rectangle.get("height", "0")) >= 0
            if None in values:
                assert "n/a" in chart.read_text()

    # Verify actual captured commands, not just script text.
    expected = {
        "readOptions": "cacheWorkingInMemory",
        "batchSize": "0",
        "innerBatchSize": "-1",
        "maxConcurrentSwaths": "-1",
        "outer_loop_threads": "-1",
    }
    analysis_expected = {
        "-Inference:peptide:global": "true",
        "-Inference:protein:global": "true",
        "-Export:results:max_rs_peakgroup_qvalue": "0.05",
        "-Export:results:use_peptide_scores": "true",
        "-Export:results:use_protein_scores": "true",
        "-Export:results:max_global_peptide_qvalue": "1.01",
        "-Export:results:max_global_protein_qvalue": "1.01",
        "-Export:results:exclude_decoys": "false",
        "-Export:matrix:max_rs_peakgroup_qvalue": "0.01",
        "-Export:matrix:use_peptide_scores": "true",
        "-Export:matrix:use_protein_scores": "true",
        "-Export:matrix:max_global_peptide_qvalue": "0.01",
        "-Export:matrix:max_global_protein_qvalue": "0.01",
    }
    for level in ("peptide", "protein", "gene"):
        for context in ("run_specific", "experiment_wide"):
            analysis_expected[f"-Inference:{level}:{context}"] = "false"
    analysis_expected["-Inference:gene:global"] = "false"
    # Source the live shared arrays, independent of immutable historical captures.
    sourced = subprocess.run(['bash', '-c', 'source "$1"; printf "%s\\0" "${OPEN_DIA_SPEED_ARGS[@]}" "${OPEN_DIA_ANALYSIS_ARGS[@]}"',
                              'bash', str(Path(__file__).with_name('lib.sh').resolve())],
                             capture_output=True, check=True)
    current_args = sourced.stdout.decode().rstrip('\0').split('\0')
    current_command = ['env', f"LD_LIBRARY_PATH={os.environ.get('LD_LIBRARY_PATH', '')}",
                       os.environ.get('OPEN_DIA', 'OpenDIA'), '-in', 'missing.mzML', '-tr', 'missing.tsv', *current_args]
    commands = [('current sourced lib.sh', current_command)]
    commands.extend((argument, shlex.split(Path(argument).read_text())) for argument in sys.argv[1:])
    historical = 0
    for argument, command in commands:
        command_expected = dict(analysis_expected)
        score_flags = [command[command.index(f'-Export:results:use_{level}_scores') + 1] for level in ('peptide', 'protein')]
        if score_flags == ['false', 'false'] and argument != 'current sourced lib.sh':
            # Historical raw exports disabled annotation, not inference or 1% matrices.
            historical += 1
            for level in ('peptide', 'protein'):
                command_expected[f'-Export:results:use_{level}_scores'] = 'false'
                del command_expected[f'-Export:results:max_global_{level}_qvalue']
            print(f'Historical raw export mode (global annotation disabled): {argument}')
        for option, value in command_expected.items():
            assert command.count(option) == 1, (argument, option)
            assert command[command.index(option) + 1] == value, (argument, option)
        for name, value in expected.items():
            option = f"-TargetedDataExtraction:{name}"
            assert command.count(option) == 1, (argument, option)
            assert command[command.index(option) + 1] == value, (argument, option)
        if os.environ.get("OPEN_DIA"):
            # Run the real parser, then stop at a missing library before costly I/O.
            # -write_ini is not a dry run: it writes defaults, not resolved values.
            with tempfile.TemporaryDirectory() as temporary:
                # Captures may start with env TMPDIR=... before the original env.
                executable_index = command.index('-in') - 1
                command[executable_index] = str(Path(os.environ["OPEN_DIA"]).resolve())
                library = Path(command[command.index("-tr") + 1])
                command[command.index("-tr") + 1] = str(Path(temporary) / ("missing" + library.suffix))
                command.extend(["-debug", "2"])
                accepted = subprocess.run(command, capture_output=True, text=True, timeout=30,
                                          env=dict(os.environ, OPENMS_DISABLE_UPDATE_CHECK="ON"))
                log = accepted.stdout + accepted.stderr
                assert accepted.returncode != 0 and "Error: File not found" in log
                assert "Parameters passed to" not in log and "Unknown option(s)" not in log
                for name, value in expected.items():
                    assert f'"TargetedDataExtraction|{name}" -> "{value}"' in log, (argument, name)
                for option, value in command_expected.items():
                    prefix, name = option[1:].rsplit(":", 1)
                    key = f"{prefix}|{name}"
                    assert f'"{key}" -> "{value}"' in log, (argument, option)
    print(f"PASS: current sourced export settings, {historical} historical raw-export captures, metrics, invalid q-values, resource failures, SVG edge cases, PDF failures, {len(sys.argv) - 1} captured commands")
    if os.environ.get("OPEN_DIA"):
        print("PASS: real OpenDIA resolved current shared settings and captured command modes; this parser check does not process the data")


if __name__ == "__main__":
    main()
