#!/usr/bin/env python3
"""Exercise real Linux workload sampling and benchmark-only disk termination."""

import json
import os
import resource
import signal
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

from monitor_host import cgroup_memory, disk_bytes, job_memory, process_ticks, process_tree


def running(pid):
    tree = process_tree(pid)
    return pid in tree and tree[pid][0] != "Z"


def test_cpu_accounting():
    before = process_ticks(os.getpid())[0]
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    subprocess.run([
        sys.executable, "-c",
        "import time; end = time.process_time() + 0.25\nwhile time.process_time() < end: pass",
    ], check=True)
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    actual_cpu = usage_after.ru_utime + usage_after.ru_stime - usage_before.ru_utime - usage_before.ru_stime
    observed_cpu = (process_ticks(os.getpid())[0] - before) / os.sysconf("SC_CLK_TCK")
    assert observed_cpu >= actual_cpu - 0.04, (observed_cpu, actual_cpu)

    # A subprocess belongs to the creating thread's children list, not the main thread's.
    ready, done = threading.Event(), threading.Event()
    children = []

    def worker():
        child = subprocess.Popen(["sleep", "30"])
        children.append(child)
        ready.set()
        try:
            done.wait(timeout=5)
        finally:
            child.terminate()
            child.wait(timeout=5)

    thread = threading.Thread(target=worker)
    thread.start()
    try:
        assert ready.wait(timeout=5)
        _, descendants = process_ticks(os.getpid())
        assert descendants.count(children[0].pid) == 1
    finally:
        done.set()
        thread.join(timeout=5)
        assert not thread.is_alive()


def test_sampling(root, mode):
    monitor = Path(__file__).with_name("monitor_host.py")
    output = root / mode
    output.mkdir()
    guard = mode != "normal"
    child = None
    if mode == "orphan_thread":
        # TERM kills this parent but not its worker-thread child. The guard must retain it.
        program = """import subprocess, sys, threading, time

def worker():
    child = subprocess.Popen([sys.executable, '-c',
        'import os,signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); print(os.getpid(),flush=True); time.sleep(30)'])
    child.wait()
threading.Thread(target=worker).start()
time.sleep(30)
"""
    elif mode == "resistant":
        program = "import os,signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); print(os.getpid(),flush=True); time.sleep(30)"
    else:
        program = f"import os,time; print(os.getpid(),flush=True); time.sleep({30 if guard else 0.5})"
    job = subprocess.Popen([
        "/usr/bin/time", "-v", "-o", str(output / "resources.txt"), sys.executable, "-c", program,
    ], stdout=subprocess.PIPE, text=True)
    sampler = None
    try:
        child = int(job.stdout.readline())
        sampler = subprocess.Popen([
            sys.executable, monitor, str(job.pid), output, "--interval", "0.1",
            "--min-free-disk-gib", "1000000000" if guard else "0",
        ])
        status = job.wait(timeout=8)
        assert sampler.wait(timeout=8) == (2 if guard else 0)
        assert not running(child), (mode, child)
        data = json.loads((output / "host_workload.json").read_text())
        assert data["sample_count"] > 0
        assert data['resource_schema_version'] == 1
        assert data['job_tree_rss_bytes_sampled']['max'] > 0
        assert data['job_largest_process_hwm_bytes_observed']['max'] > 0
        assert data['output_apparent_bytes_sampled']['max'] >= 0
        assert data['temp_apparent_bytes_sampled']['max'] >= 0
        assert data['cgroup_path_requested'] is None
        assert 'cgroup_memory_peak_bytes' not in data
        assert data["logical_cpus"] >= 1
        assert data['system']['logical_cpus'] == data['logical_cpus']
        assert data['system']['affinity_logical_cpus'] >= 1
        assert data['system']['mem_total_mib'] > 0
        assert 0 < data['system']['mem_available_mib_at_start'] <= data['system']['mem_total_mib']
        assert data['started_at'] <= data['finished_at']
        assert not {'hostname', 'serial_number', 'username'} & data['system'].keys()
        assert data["mem_available_mib"]["min"] > 0
        assert bool(data["abort_reason"]) == guard
        assert (status != 0) == guard
        raw = [json.loads(line) for line in (output / "host_workload.jsonl").read_text().splitlines()]
        assert len(raw) == data["sample_count"]
        # GNU time itself must survive to write its resource/failure record.
        assert "Maximum resident set size" in (output / "resources.txt").read_text()
        if not guard:
            assert 0 <= data["host_cpu_busy_percent"]["mean"] <= 100
            assert 0 <= data['host_cpu_idle_percent']['mean'] <= 100
            assert abs(data['host_cpu_idle_thread_equivalents']['mean'] -
                       data['logical_cpus'] * data['host_cpu_idle_percent']['mean'] / 100) < 1e-9
    finally:
        if sampler is not None and sampler.poll() is None:
            sampler.kill()
            sampler.wait(timeout=5)
        # Only these synthetic processes are candidates for test cleanup.
        descendants = process_tree(job.pid)
        if child is not None and running(child):
            descendants.update(process_tree(child))
        for pid in reversed(descendants):
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        job.wait(timeout=5)
        job.stdout.close()


def main():
    test_cpu_accounting()
    memory = job_memory(os.getpid())
    assert memory['job_tree_rss_bytes_sampled'] > 0
    assert memory['job_largest_process_hwm_bytes_observed'] > 0
    assert cgroup_memory(None, os.getpid()) == {}
    assert cgroup_memory(Path('/sys/fs/cgroup'), os.getpid()) == {}
    with tempfile.TemporaryDirectory(dir=os.environ.get("JCODE_SCRATCH_DIR")) as temporary:
        files = Path(temporary) / 'files'
        files.mkdir()
        original = files / 'original'
        original.write_bytes(b'x' * 123)
        os.link(original, files / 'hardlink')
        (files / 'symlink').symlink_to(original)
        assert disk_bytes(files) == 123
        assert disk_bytes(files / 'missing') == 0
        output = Path(temporary) / 'failed_cgroup'
        output.mkdir()
        done = output / 'done'
        release = output / 'release'
        done.write_text('0\n')
        with subprocess.Popen(['sleep', '0.3']) as job:
            checked = subprocess.run([sys.executable, Path(__file__).with_name('monitor_host.py'),
                str(job.pid), output, '--interval', '0.05', '--min-free-disk-gib', '0',
                '--cgroup', '/sys/fs/cgroup/opendia-bench-nonexistent.scope',
                '--cgroup-done', done, '--cgroup-release', release], timeout=5)
            # Reap the child independently, since a zombie also has a /proc entry.
            job.wait(timeout=5)
        assert checked.returncode == 2
        assert release.read_text() == 'failed\n'
        assert json.loads((output / 'host_workload.json').read_text())['abort_reason']
        # An unrelated process must remain alive through every guard test.
        unrelated = subprocess.Popen(["sleep", "30"])
        try:
            for mode in ("normal", "cooperative", "resistant", "orphan_thread"):
                test_sampling(Path(temporary), mode)
                assert unrelated.poll() is None
        finally:
            unrelated.terminate()
            unrelated.wait(timeout=5)
    print("PASS: host sampling, waited-child CPU, threaded descendants, TERM/KILL and orphan guards, GNU time records, unrelated process safety")


if __name__ == "__main__":
    main()
